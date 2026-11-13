"""Loading and validating setup.yaml files.

`load_setup` returns a typed `Setup`. Every error is a `ConfigError` that names
the offending key. schema/setup-1.schema.json describes the same structure.

Two analysis types: `chsh` (detectors A+, A-, B+, B-; four setting pairs) and
`hom` (detectors A and B on the beamsplitter outputs; acquisitions "dip" at
zero path delay and "reference" at `reference_delay_ps`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from dadl import _checks
from dadl.analysis.chsh import CHANNELS, Window
from dadl.analysis.hom import HOM_CHANNELS, HOM_SETTINGS
from dadl.errors import ConfigError
from dadl.sim import CHSH_SETTINGS, DetectorParams, PairSource

SUPPORTED_SCHEMA_VERSIONS = (1,)
REQUIRED_KEYS = ("schema_version", "name", "source", "detectors", "analysis")
MAX_FILE_BYTES = 1_000_000

TOP_LEVEL_KEYS = frozenset(REQUIRED_KEYS) | {"backend", "seed", "channels", "policy", "hardware"}
#: Source models DADL implements. Pulsed sources need slots locked to the laser (planned).
SOURCE_TYPES = ("spdc-cw",)
HARDWARE_KEYS = frozenset({"tagger", "detectors"})
DETECTOR_KEYS = frozenset(
    {"family", "efficiency", "dark_rate_cps", "dead_time_ps", "jitter_sigma_ps", "afterpulse_prob",
     "afterpulse_tau_ps", "delay_ps"}
)
SOURCE_KEYS = frozenset({"type", "pair_rate_per_s", "visibility", "coherence_time_ps"})
ANALYSIS_KEYS = frozenset({"type", "settings_deg", "duration_s_per_setting", "window", "reference_delay_ps"})
ANALYSIS_TYPES = ("chsh", "hom")
WINDOW_KEYS = frozenset({"width_ps", "source", "delays_ps"})


class _UniqueKeyLoader(yaml.SafeLoader):
    """SafeLoader that refuses duplicate mapping keys instead of keeping the last one."""


def _construct_mapping(loader: _UniqueKeyLoader, node: yaml.MappingNode, deep: bool = False) -> dict[Any, Any]:
    loader.flatten_mapping(node)
    seen: set[Any] = set()
    for key_node, _ in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in seen:
            raise ConfigError(f"duplicate key {key!r} at line {key_node.start_mark.line + 1}")
        seen.add(key)
    return yaml.SafeLoader.construct_mapping(loader, node, deep=deep)


_UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


@dataclass(frozen=True)
class Setup:
    """A validated setup file."""

    name: str
    backend: str
    seed: int | None
    source: PairSource
    detectors: dict[str, DetectorParams]
    settings_deg: dict[str, float]
    duration_s_per_setting: float
    window: Window
    #: Maps channel numbers in a tag file to channel names (CHSH: A+, A-, B+, B-; HOM: A, B).
    channels: dict[int, str]
    path: Path | None = None
    data: dict[str, Any] = field(default_factory=dict, repr=False)
    analysis_type: str = "chsh"
    #: Tagger serial and detector serials per channel name (hardware setups; recorded with calibrations).
    tagger: str | None = None
    detector_serials: dict[str, str] = field(default_factory=dict)
    #: HOM only: the path delay of the reference acquisition, far from the dip.
    reference_delay_ps: float | None = None

    def channel_id(self, name: str):
        """The ChannelId of a channel name: tagger serial (or "sim") and input number, with role and detector."""
        from dadl.data import ChannelId

        number = {n: k for k, n in self.channels.items()}
        if self.backend == "sim":
            device = "sim"
        elif self.tagger:
            device = self.tagger
        else:
            raise ConfigError(
                f"{self.name}: calibrations of hardware data need hardware.tagger (the tagger's serial) in "
                "setup.yaml, so the record identifies the physical device, not just the channel name"
            )
        detector = self.detector_serials.get(name)
        if self.backend != "sim" and not detector:
            # Calibration-lab review, cycle 3: the tagger port was identified but not the detector
            # being calibrated.
            raise ConfigError(
                f"{self.name}: calibrations of hardware data need the detector's serial for channel {name!r} "
                f"under hardware.detectors in setup.yaml, so the record identifies the calibrated item"
            )
        return ChannelId(device, number[name], role=name, detector=detector)

    @property
    def channel_names(self) -> tuple[str, ...]:
        return HOM_CHANNELS if self.analysis_type == "hom" else CHANNELS

    @property
    def setting_keys(self) -> tuple:
        """Acquisition keys in tag-file setting order."""
        if self.analysis_type == "hom":
            return HOM_SETTINGS
        from dadl.tagfile import SETTING_PAIRS

        return SETTING_PAIRS


def _mapping(value: Any, where: str) -> dict[Any, Any]:
    if not isinstance(value, dict):
        raise ConfigError(f"{where} must be a mapping, got {type(value).__name__}")
    return value


def _no_unknown(mapping: dict[Any, Any], allowed: frozenset[str], where: str) -> None:
    unknown = sorted(str(k) for k in mapping if k not in allowed)
    if unknown:
        raise ConfigError(f"unknown key(s) in {where}: {', '.join(unknown)}; allowed: {', '.join(sorted(allowed))}")


def _wrap(where: str, fn, *args, **kwargs):
    try:
        return fn(*args, **kwargs)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{where}: {exc}") from exc


def read_yaml(path: str | Path) -> dict[str, Any]:
    """Read a YAML mapping, refusing duplicate keys, oversize files and bad encodings."""
    path = Path(path)
    try:
        size = path.stat().st_size
        if size > MAX_FILE_BYTES:
            raise ConfigError(f"setup file is {size} bytes; the limit is {MAX_FILE_BYTES}")
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError as exc:
        raise ConfigError(f"setup file not found: {path}") from exc
    except IsADirectoryError as exc:
        raise ConfigError(f"setup path is a directory, not a file: {path}") from exc
    except UnicodeDecodeError as exc:
        raise ConfigError(f"setup file is not UTF-8 text: {path}") from exc
    except OSError as exc:
        raise ConfigError(f"cannot read setup file {path}: {exc}") from exc
    try:
        data = yaml.load(text, Loader=_UniqueKeyLoader)  # noqa: S506 - SafeLoader subclass
    except yaml.YAMLError as exc:
        raise ConfigError(f"setup file is not valid YAML: {exc}") from exc
    return _mapping(data, "the setup file's top level")


def load_setup(path: str | Path) -> Setup:
    """Load and validate a setup file."""
    path = Path(path)
    data = read_yaml(path)
    setup = parse_setup(data)
    return Setup(**{**setup.__dict__, "path": path})


def parse_setup(data: dict[str, Any]) -> Setup:
    """Validate an already-parsed setup mapping."""
    missing = [k for k in REQUIRED_KEYS if k not in data]
    if missing:
        raise ConfigError(f"setup file is missing keys: {', '.join(missing)}")
    _no_unknown(data, TOP_LEVEL_KEYS, "the setup file")

    version = data["schema_version"]
    if type(version) is not int or version not in SUPPORTED_SCHEMA_VERSIONS:
        raise ConfigError(
            f"schema_version {version!r} is not supported; supported: {list(SUPPORTED_SCHEMA_VERSIONS)}"
        )
    name = data["name"]
    if not isinstance(name, str) or not name.strip():
        raise ConfigError("name must be a non-empty string")
    backend = data.get("backend", "sim")
    if backend not in ("sim", "file"):
        raise ConfigError(f"backend must be 'sim' or 'file' in phase 1, got {backend!r}")
    seed = data.get("seed")
    if seed is not None:
        seed = _wrap("seed", _checks.integer, seed, "seed", minimum=0)
    if data.get("policy", "default") != "default":
        raise ConfigError(
            "policy must be 'default' in phase 1; loosen rules in code with Policy.allow(analysis, correction, reason)"
        )

    # source
    src = _mapping(data["source"], "source")
    _no_unknown(src, SOURCE_KEYS, "source")
    if src.get("type", "spdc-cw") not in SOURCE_TYPES:
        raise ConfigError(
            f"source.type {src.get('type')!r} is not supported; DADL models {list(SOURCE_TYPES)} only "
            "(pulsed sources need coincidence slots locked to the laser, which is planned)"
        )
    source = _wrap(
        "source",
        PairSource,
        pair_rate_per_s=src.get("pair_rate_per_s", PairSource.pair_rate_per_s),
        visibility=src.get("visibility", PairSource.visibility),
        coherence_time_ps=src.get("coherence_time_ps", PairSource.coherence_time_ps),
    )

    # detectors
    dets = _mapping(data["detectors"], "detectors")
    if not dets:
        raise ConfigError("detectors must list at least one channel")
    detectors: dict[str, DetectorParams] = {}
    for ch, params in dets.items():
        params = _mapping(params, f"detectors.{ch}")
        _no_unknown(params, DETECTOR_KEYS, f"detectors.{ch}")
        kwargs = {k: v for k, v in params.items() if k != "family"}
        detectors[str(ch)] = _wrap(f"detectors.{ch}", DetectorParams, **kwargs)

    # analysis
    analysis = _mapping(data["analysis"], "analysis")
    _no_unknown(analysis, ANALYSIS_KEYS, "analysis")
    kind = analysis.get("type")
    if kind not in ANALYSIS_TYPES:
        raise ConfigError(f"analysis.type must be one of {list(ANALYSIS_TYPES)}, got {kind!r}")
    names = HOM_CHANNELS if kind == "hom" else CHANNELS
    missing_ch = [ch for ch in names if ch not in detectors]
    if missing_ch:
        raise ConfigError(f"a {kind} setup needs detectors {list(names)}; missing {missing_ch}")
    settings_deg: dict[str, float] = {}
    reference_delay = None
    if kind == "chsh":
        if "reference_delay_ps" in analysis:
            raise ConfigError("analysis.reference_delay_ps belongs to hom setups")
        settings = _mapping(analysis.get("settings_deg", dict(CHSH_SETTINGS)), "analysis.settings_deg")
        missing_s = [s for s in CHSH_SETTINGS if s not in settings]
        if missing_s:
            raise ConfigError(f"analysis.settings_deg is missing {missing_s}")
        settings_deg = {str(k): _wrap(f"analysis.settings_deg.{k}", _checks.finite, v, k) for k, v in settings.items()}
    else:
        if "settings_deg" in analysis:
            raise ConfigError("analysis.settings_deg belongs to chsh setups")
        if "reference_delay_ps" not in analysis:
            raise ConfigError("a hom setup needs analysis.reference_delay_ps (a path delay far from the dip)")
        reference_delay = _wrap(
            "analysis.reference_delay_ps", _checks.finite, analysis["reference_delay_ps"], "reference_delay_ps"
        )
        if reference_delay == 0:
            raise ConfigError("analysis.reference_delay_ps must be away from the dip (non-zero)")
    duration = _wrap(
        "analysis.duration_s_per_setting",
        _checks.finite,
        analysis.get("duration_s_per_setting", 1.0),
        "duration_s_per_setting",
        minimum=0.0,
    )
    if duration <= 0:
        raise ConfigError("analysis.duration_s_per_setting must be positive")

    win = _mapping(analysis.get("window"), "analysis.window")
    _no_unknown(win, WINDOW_KEYS, "analysis.window")
    if win.get("source") != "plan":
        raise ConfigError("analysis.window.source must be 'plan': windows are declared before acquisition")
    if "width_ps" not in win:
        raise ConfigError("analysis.window.width_ps is required")
    delays = _mapping(win.get("delays_ps", {}), "analysis.window.delays_ps")
    window = _wrap("analysis.window", Window, width_ps=win["width_ps"], source="plan", delays_ps=delays)

    # channel map for tag files: number -> name; default 0, 1, ... in channel-name order
    raw_map = _mapping(data.get("channels", dict(enumerate(names))), "channels")
    channels: dict[int, str] = {}
    for num, ch in raw_map.items():
        num = _wrap("channels", _checks.integer, num, "channel number", minimum=0)
        if ch not in names:
            raise ConfigError(f"channels: {num} maps to {ch!r}; names must be one of {list(names)}")
        channels[num] = ch
    if sorted(channels.values()) != sorted(names):
        raise ConfigError(f"channels must map exactly one number to each of {list(names)}")

    hw = _mapping(data.get("hardware", {}), "hardware")
    _no_unknown(hw, HARDWARE_KEYS, "hardware")
    tagger = hw.get("tagger")
    if tagger is not None and (not isinstance(tagger, str) or not tagger.strip()):
        raise ConfigError("hardware.tagger must be the tagger's serial, e.g. HydraHarp-1234567")
    serials = {str(k): str(v) for k, v in _mapping(hw.get("detectors", {}), "hardware.detectors").items()}
    unknown_roles = sorted(set(serials) - set(names))
    if unknown_roles:
        raise ConfigError(f"hardware.detectors names unknown channels {unknown_roles}")

    return Setup(
        name=name,
        backend=backend,
        seed=seed,
        source=source,
        detectors=detectors,
        settings_deg=settings_deg,
        duration_s_per_setting=duration,
        window=window,
        channels=channels,
        data=data,
        analysis_type=kind,
        tagger=tagger,
        detector_serials=serials,
        reference_delay_ps=reference_delay,
    )
