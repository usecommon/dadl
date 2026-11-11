"""Loading and validating setup.yaml files.

`load_setup` returns a typed `Setup`. Every error is a `ConfigError` that names
the offending key. A JSON Schema for setup files (dadl check-config) and
migrations (dadl migrate) arrive in phase 2.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from dadl import _checks
from dadl.analysis.chsh import CHANNELS, Window
from dadl.errors import ConfigError
from dadl.sim import CHSH_SETTINGS, DetectorParams, PairSource

SUPPORTED_SCHEMA_VERSIONS = (1,)
REQUIRED_KEYS = ("schema_version", "name", "source", "detectors", "analysis")
MAX_FILE_BYTES = 1_000_000

TOP_LEVEL_KEYS = frozenset(REQUIRED_KEYS) | {"backend", "seed", "channels", "policy"}
DETECTOR_KEYS = frozenset(
    {"family", "efficiency", "dark_rate_cps", "dead_time_ps", "jitter_sigma_ps", "afterpulse_prob",
     "afterpulse_tau_ps", "delay_ps"}
)
SOURCE_KEYS = frozenset({"type", "pair_rate_per_s", "visibility"})
ANALYSIS_KEYS = frozenset({"type", "settings_deg", "duration_s_per_setting", "window"})
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
    #: Maps channel numbers in a tag file to channel names (A+, A-, B+, B-).
    channels: dict[int, str]
    path: Path | None = None
    data: dict[str, Any] = field(default_factory=dict, repr=False)


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
    source = _wrap(
        "source",
        PairSource,
        pair_rate_per_s=src.get("pair_rate_per_s", PairSource.pair_rate_per_s),
        visibility=src.get("visibility", PairSource.visibility),
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
    if analysis.get("type") != "chsh":
        raise ConfigError(f"analysis.type must be 'chsh' in phase 1, got {analysis.get('type')!r}")
    missing_ch = [ch for ch in CHANNELS if ch not in detectors]
    if missing_ch:
        raise ConfigError(f"a chsh setup needs detectors {list(CHANNELS)}; missing {missing_ch}")
    settings = _mapping(analysis.get("settings_deg", dict(CHSH_SETTINGS)), "analysis.settings_deg")
    missing_s = [s for s in CHSH_SETTINGS if s not in settings]
    if missing_s:
        raise ConfigError(f"analysis.settings_deg is missing {missing_s}")
    settings_deg = {str(k): _wrap(f"analysis.settings_deg.{k}", _checks.finite, v, k) for k, v in settings.items()}
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

    # channel map for tag files: number -> name; default 0..3 -> A+, A-, B+, B-
    raw_map = _mapping(data.get("channels", dict(enumerate(CHANNELS))), "channels")
    channels: dict[int, str] = {}
    for num, ch in raw_map.items():
        num = _wrap("channels", _checks.integer, num, "channel number", minimum=0)
        if ch not in CHANNELS:
            raise ConfigError(f"channels: {num} maps to {ch!r}; names must be one of {list(CHANNELS)}")
        channels[num] = ch
    if sorted(channels.values()) != sorted(CHANNELS):
        raise ConfigError(f"channels must map exactly one number to each of {list(CHANNELS)}")

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
    )
