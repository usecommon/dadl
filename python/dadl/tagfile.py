"""Reading and writing CHSH time-tag files (phase 1 generic format).

A tag file is a table with three integer columns:

    setting        CHSH: 0..3, the setting pair, in the order
                   (a, b), (a, b'), (a', b), (a', b'); HOM: 0 dip, 1 reference
    channel        detector input number, mapped to A+, A-, B+, B- by the
                   `channels:` map in setup.yaml (default 0, 1, 2, 3)
    timestamp_ps   int64 time tag in picoseconds

Two encodings are accepted:

- `.npz`: a NumPy archive with arrays named `setting`, `channel`, `timestamp_ps`;
- `.csv`: a header line `setting,channel,timestamp_ps`, then one row per tag.

A third form is a folder of PicoQuant PTU files (T2 or T3), one per
acquisition, named setting0.ptu, setting1.ptu, ... in setting order; PTU input
numbers are the channel numbers. In T3 files, input 0 (the syncs that started
a recorded measurement, derived by the adapter) is read only when mapped. Other vendor formats (Swabian .ttbin, Photon-HDF5) need
exporting to one of these.
"""

from __future__ import annotations

import math
import zipfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np

from dadl import _checks
from dadl.analysis.chsh import CHANNELS, SIGNS
from dadl.errors import ConfigError

COLUMNS = ("setting", "channel", "timestamp_ps")
SETTING_PAIRS: tuple[tuple[str, str], ...] = tuple(SIGNS)

Runs = dict[Any, dict[str, np.ndarray]]


def ptu_files(path: str | Path) -> list[Path]:
    """The PTU files of a tag path: a folder's setting*.ptu in order, a single .ptu file, else []."""
    path = Path(path)
    if path.is_dir():
        return sorted(path.glob("setting*.ptu"))
    if path.suffix.lower() == ".ptu" and path.is_file():
        return [path]
    return []


def read_tags(
    path: str | Path,
    channels: dict[int, str] | None = None,
    settings: Sequence = SETTING_PAIRS,
    only: Sequence[int] | None = None,
) -> Runs:
    """Read a tag file into the `runs` structure that `dadl.chsh` takes.

    `settings` names the setting indices in order (CHSH: the four setting
    pairs; HOM: "dip", "reference"). For PTU input, `only` reads just those
    setting indices (the others are empty), and a single .ptu file is read as
    setting 0 (council round 9, Tw6: one dark file needed copies as setting1-3).
    """
    path = Path(path)
    channels = channels or dict(enumerate(CHANNELS))
    settings = tuple(settings)
    if path.is_dir():
        return _read_ptu_dir(path, channels, settings, only)
    if path.suffix.lower() == ".ptu" and path.is_file():
        return _read_ptu_dir(path.parent, channels, settings, (0,), single=path)
    if path.exists() and not path.is_file():
        raise ConfigError(f"{path}: not a regular file")
    try:
        if path.suffix == ".npz":
            _check_npz_size(path)
            with np.load(path, allow_pickle=False) as npz:
                missing = [c for c in COLUMNS if c not in npz.files]
                if missing:
                    raise ConfigError(f"{path}: missing array(s) {missing}; need {list(COLUMNS)}")
                cols = {c: npz[c] for c in COLUMNS}
        elif path.suffix == ".csv":
            # utf-8-sig: Excel's "CSV UTF-8" starts with a byte-order mark.
            with path.open(encoding="utf-8-sig") as fh:
                first = fh.readline()
                skip = 1
                if first.strip() == SIMULATED_CSV_MARK:
                    first, skip = fh.readline(), 2
                header = [h.strip() for h in first.split(",")]
            if header != list(COLUMNS):
                raise ConfigError(f"{path}: header must be {','.join(COLUMNS)}, got {','.join(header)}")
            # Read as text and check every cell is an integer: some NumPy versions
            # silently truncate "1.5" when asked for int64 directly.
            text = np.loadtxt(path, delimiter=",", skiprows=skip, dtype=str, ndmin=2, encoding="utf-8-sig")
            if text.size:
                text = np.char.strip(text)
                # ASCII only: str.isdigit and int() accept full-width digits (archivist review).
                try:
                    text.astype(np.bytes_)
                except UnicodeEncodeError as exc:
                    raise ConfigError(f"{path}: non-ASCII characters in the tag columns") from exc
                digits = np.char.isdigit(np.char.lstrip(text, "+-"))
                if not digits.all():
                    row, col = map(int, np.argwhere(~digits)[0])
                    raise ConfigError(
                        f"{path}: row {row + 2}, column {COLUMNS[col]}: {text[row, col]!r} is not an integer"
                    )
            cols = {
                c: text[:, i].astype(np.int64) if text.size else np.array([], dtype=np.int64)
                for i, c in enumerate(COLUMNS)
            }
        else:
            raise ConfigError(f"{path}: tag files must be .npz, .csv or .ptu (or a folder of settingN.ptu files)")
    except FileNotFoundError as exc:
        raise ConfigError(f"tag file not found: {path}") from exc
    except PermissionError as exc:
        raise ConfigError(f"cannot read {path}: permission denied") from exc
    except (zipfile.BadZipFile, EOFError) as exc:
        raise ConfigError(f"{path}: not a readable .npz file (truncated or corrupt): {exc}") from exc
    except (ValueError, OverflowError) as exc:
        raise ConfigError(f"{path}: cannot parse as integer columns: {exc}") from exc

    try:
        setting = _checks.timestamps(cols["setting"], "setting")
        channel = _checks.timestamps(cols["channel"], "channel")
        stamps = _checks.timestamps(cols["timestamp_ps"], "timestamp_ps")
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{path}: {exc}") from exc
    if not (setting.size == channel.size == stamps.size):
        raise ConfigError(f"{path}: columns have different lengths")
    bad_s = np.setdiff1d(np.unique(setting), np.arange(len(settings)))
    if bad_s.size:
        raise ConfigError(f"{path}: setting values must be 0..{len(settings) - 1}, found {bad_s.tolist()}")
    bad_c = np.setdiff1d(np.unique(channel), np.array(sorted(channels), dtype=np.int64))
    if bad_c.size:
        raise ConfigError(f"{path}: channel numbers {bad_c.tolist()} are not in the setup's channel map")

    runs: Runs = {}
    for i, pair in enumerate(settings):
        in_setting = setting == i
        runs[pair] = {}
        for num, name in channels.items():
            arr = stamps[in_setting & (channel == num)]
            # Each channel's tags must already be in time order in the file: sorting them hid a
            # tagger clock reset, which read as one overlapping run at twice the rate (operator
            # review, cycle 8).
            back = np.flatnonzero(np.diff(arr) < 0)
            if back.size:
                raise ConfigError(
                    f"{path}: channel {name} of acquisition {pair}: timestamps go back in time "
                    f"({int(arr[back[0]])} then {int(arr[back[0] + 1])} ps). A tagger clock reset, or files "
                    "joined out of order? Split the acquisition at the reset; DADL does not sort tags."
                )
            runs[pair][name] = arr
    return runs


#: Largest decompressed size, and compression ratio, accepted for a .npz tag file: a 584 KB file
#: of zeros decompressed to 1.2 GB (security review, cycle 9). Tag columns compress about 2-4x.
MAX_NPZ_BYTES = 16 * 2**30
MAX_NPZ_RATIO = 200


def _check_npz_size(path: Path) -> None:
    try:
        with zipfile.ZipFile(path) as zf:
            infos = zf.infolist()
    except zipfile.BadZipFile as exc:
        raise ConfigError(f"{path}: not a readable .npz file (truncated or corrupt): {exc}") from exc
    total = sum(i.file_size for i in infos)
    packed = sum(i.compress_size for i in infos) or 1
    if total > MAX_NPZ_BYTES or (total > 2**26 and total / packed > MAX_NPZ_RATIO):
        raise ConfigError(
            f"{path}: decompresses to {total} bytes ({total / packed:.0f}x its stored size); refused as a likely "
            f"decompression bomb (limits {MAX_NPZ_BYTES} bytes, {MAX_NPZ_RATIO}x above 64 MiB)"
        )


#: First line of a CSV tag file written from a simulation.
SIMULATED_CSV_MARK = "# dadl-simulated"


def write_tags(
    path: str | Path,
    runs: Runs,
    channels: dict[int, str] | None = None,
    settings: Sequence = SETTING_PAIRS,
    *,
    simulated: bool,
) -> Path:
    """Write `runs` to a .npz or .csv tag file (the inverse of `read_tags`).

    `simulated` is required: a read-trim-write round trip must not drop the marking (red-team
    review, cycle 5). Simulated CSV files start with the line `# dadl-simulated`.
    """
    path = Path(path)
    channels = channels or dict(enumerate(CHANNELS))
    number = {name: num for num, name in channels.items()}
    parts = []
    for i, pair in enumerate(settings):
        for name, tags in runs[pair].items():
            tags = np.asarray(tags, dtype=np.int64)
            parts.append(np.stack([np.full(tags.size, i), np.full(tags.size, number[name]), tags], axis=1))
    table = np.concatenate(parts) if parts else np.empty((0, 3), dtype=np.int64)
    table = table[np.lexsort((table[:, 2], table[:, 1], table[:, 0]))]
    if path.suffix == ".npz":
        extra = {"dadl_simulated": np.array(True)} if simulated else {}
        np.savez_compressed(path, **{c: table[:, i] for i, c in enumerate(COLUMNS)}, **extra)
    elif path.suffix == ".csv":
        header = (SIMULATED_CSV_MARK + "\n" if simulated else "") + ",".join(COLUMNS)
        np.savetxt(path, table, fmt="%d", delimiter=",", header=header, comments="")
    else:
        raise ConfigError(f"{path}: tag files must be .npz or .csv")
    return path


def ptu_file(directory: str | Path, index: int) -> Path:
    """Where acquisition `index` lives in a folder of PTU files."""
    return Path(directory) / f"setting{index}.ptu"


def _read_ptu_dir(
    directory: Path, channels: dict[int, str], settings: tuple, only: Sequence[int] | None = None,
    single: Path | None = None,
) -> Runs:
    from dadl.adapters.ptu import PtuFile
    from dadl.data import RunTags
    from dadl.errors import AdapterError

    runs: Runs = {}
    for i, key in enumerate(settings):
        if only is not None and i not in only:
            runs[key] = RunTags({name: np.empty(0, dtype=np.int64) for name in channels.values()})
            continue
        path = single if single is not None else ptu_file(directory, i)
        if not path.is_file():
            raise ConfigError(f"{directory}: missing {path.name} (one PTU file per acquisition, setting0.ptu, ...)")
        try:
            reader = PtuFile(path)
            tagset = reader.read()
        except AdapterError as exc:
            raise ConfigError(str(exc)) from exc
        short_range = reader.dtime_range_short()
        if short_range:
            # Council round 7 (G4): rates and dark calibrations need every photon of the period.
            raise AdapterError(f"{path}: {short_range}. DADL does not analyse such files; record with a coarser "
                               "dtime resolution, so the range covers the sync period")
        unknown_records = reader.record_counts().get("unknown", 0)
        if unknown_records:
            # Tagger engineer review, cycle 2: 1000 unrecognised special records were dropped
            # silently and the analysis printed S as normal.
            raise ConfigError(
                f"{path}: {unknown_records} special records of an unrecognised kind; DADL does not analyse files "
                "it cannot fully account for. Inspect them with dadl.adapters.PtuFile(path).events()."
            )
        divided = [e.value for e in reader.events() if e.kind == "sync-divider"]
        if divided and 0 in channels:
            raise ConfigError(
                f"{path}: input 0 (sync) is divided by {divided[0]} (HWSync_Divider) but is mapped to "
                f"channel {channels[0]!r}; a divided input records one event in {divided[0]} and cannot be a detector"
            )
        if tagset.metadata.get("incomplete"):
            raise ConfigError(
                f"{path}: the acquisition lost data ({tagset.metadata.get('data_loss')}, from TTResult_StopReason); "
                "DADL does not analyse incomplete files. Read it with dadl.adapters.PtuFile to inspect it."
            )
        by_input = {ch.input: tagset[ch] for ch in tagset.channels}
        if reader.mode == "T3" and 0 not in channels:
            # Derived from the tag records, not recorded: no record goes unread when it is left out.
            by_input.pop(0, None)
        unknown = sorted(set(by_input) - set(channels))
        if unknown:
            raise ConfigError(
                f"{path}: inputs {unknown} are not in the setup's channel map. On HydraHarp, MultiHarp and "
                "TimeHarp 260 files input 0 is the sync input and detector inputs are numbered from 1 (as in "
                "PicoQuant's software); map them under `channels:` in setup.yaml"
            )
        meta: dict[str, Any] = {"mode": reader.mode, "tick_ps": reader.tick_ps, "lattice_ps": reader.lattice_ps}
        if reader.mode == "T3":
            # Council round 7 (G3): analyses label T3 input and refuse a corrected value unless the
            # source is declared CW.
            meta["sync_period_ps"] = tagset.metadata["sync_period_ps"]
            if reader.sync_divider is not None:
                meta["sync_divider"] = reader.sync_divider
        runs[key] = RunTags(
            {name: by_input.get(num, np.empty(0, dtype=np.int64)) for num, name in channels.items()}, metadata=meta
        )
    return runs


def write_ptu_dir(
    directory: str | Path,
    runs: Runs,
    channels: dict[int, str] | None = None,
    settings: Sequence = SETTING_PAIRS,
    *,
    device: str = "hydraharp",
) -> list[Path]:
    """Write each acquisition of `runs` as a simulated PTU file (see `dadl.adapters.ptu_write`).

    A run whose earliest tag is negative (jitter near t = 0) is shifted as a
    whole so it starts at 0, keeping every time difference.
    """
    from dadl.adapters.ptu_write import write_ptu

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    channels = channels or dict(enumerate(CHANNELS))
    number = {name: num for num, name in channels.items()}
    out = []
    for i, key in enumerate(settings):
        tags = {number[name]: np.asarray(arr, dtype=np.int64) for name, arr in runs[key].items()}
        first = min((int(a[0]) for a in tags.values() if a.size), default=0)
        shift = -first if first < 0 else 0
        out.append(
            write_ptu(
                ptu_file(directory, i),
                {k: a + shift for k, a in tags.items()},
                device=device,
                comment=f"DADL simulation, not a recording; acquisition {key!r}; shifted by {shift} ps",
            )
        )
    return out


def tag_tick_ps(path: str | Path) -> int | None:
    """The tagger tick in whole ps for a folder of PTU files (from the header), else None.

    T2: the whole-ps tick, or None when it is not whole. T3: the dtime resolution (plus 1 ps when
    the tags are rounded), whatever the sync period: a recorded gap is within one dtime resolution
    of the true gap (council round 7, G1; it was gcd(resolution, period), or None, which a
    measured sync period reduced to 1 ps).
    """
    path = Path(path)
    files = ptu_files(path)
    if not files:
        return None
    from dadl.adapters.ptu import PtuFile

    return PtuFile(files[0]).tick_ps


def tag_lattice_ps(path: str | Path) -> int | None:
    """The lattice of pair delays in whole ps for a folder of PTU files (`PtuFile.lattice_ps`),
    else None. Histograms of pair delays bin on it; `tag_tick_ps` is the error bound of a gap and
    serves dead-time intervals only (council round 8)."""
    path = Path(path)
    files = ptu_files(path)
    if not files:
        return None
    from dadl.adapters.ptu import PtuFile

    return PtuFile(files[0]).lattice_ps


def t3_folder(path: str | Path) -> dict[str, Any] | None:
    """For a folder of PTU T3 files: the sync period (ps) and the sync divider (None when not
    divided) of its first file; None for T2 or a tag file."""
    path = Path(path)
    files = ptu_files(path)
    if not files:
        return None
    from dadl.adapters.ptu import PtuFile

    reader = PtuFile(files[0])
    if reader.mode != "T3":
        return None
    assert reader.sync_period_s is not None
    return {"sync_period_ps": reader.sync_period_s * 1e12, "sync_divider": reader.sync_divider}


#: PTU header entries that define where a tag lands in time: input offsets, trigger/CFD levels and
#: edges, dead-time settings, sync offset and divider, and the tick.
TAGGER_SETTING_PREFIXES = ("HWInpChan_", "HWSync_", "HW_InpChannels", "MeasDesc_GlobalResolution", "HW_ExternalRefClock")


def tagger_settings(path: str | Path) -> dict[str, float]:
    """Numeric tagger settings from a folder of PTU files, for an artifact's `conditions`.

    A calibrated delay holds only while these stay as they were (tagger engineer review, cycle 2:
    input offsets and trigger levels were never recorded). Files of one folder must agree, each
    setting to its tolerance (`dadl.propagate.setting_rtol`: in T3 the sync period may drift by
    1e-4; council round 7, G2).
    """
    from dadl.propagate import setting_rtol

    path = Path(path)
    files = ptu_files(path)
    if not files:
        return {}
    from dadl.adapters.ptu import PtuFile

    found: dict[str, float] | None = None
    first = files[0].name if files else ""
    for f in files:
        reader = PtuFile(f)
        header = reader.header
        # In T3 the dtime resolution places tags too (MeasDesc_GlobalResolution is the sync period).
        prefixes = TAGGER_SETTING_PREFIXES + (("MeasDesc_Resolution",) if reader.mode == "T3" else ())
        settings = {
            f"tagger_{k}": float(v)
            for k, v in sorted(header.items())
            if k.startswith(prefixes) and isinstance(v, (int, float)) and not isinstance(v, bool)
        }
        # Which meaning MeasDesc_GlobalResolution has: the T2 tick or the T3 sync period (council
        # round 9, Tw2). Compared only between artifacts and datasets of the same mode.
        settings["tagger_mode"] = 3.0 if reader.mode == "T3" else 2.0
        if found is None:
            found = settings
            continue
        changed = sorted(
            k for k in set(found) | set(settings)
            if k not in found or k not in settings
            or not math.isclose(found[k], settings[k], rel_tol=setting_rtol(k), abs_tol=0.0)
        )
        if changed:
            shown = ", ".join(f"{k}: {_fmt(found.get(k))} in {first}, {_fmt(settings.get(k))} in {f.name}"
                              for k in changed[:5])
            raise ConfigError(f"{path}: the tagger settings differ between files ({shown})")
    return found or {}


def _fmt(value: float | None) -> str:
    return "absent" if value is None else f"{value:.10g}"


def tags_are_simulated(path: str | Path) -> bool:
    """True when the tags were written by DADL's simulator: a folder of PTU files whose headers say
    so, or an .npz written by `write_tags` with simulated=True.

    A folder mixing simulated and recorded files is refused: its provenance cannot be one label.
    """
    path = Path(path)
    if path.suffix == ".npz" and path.is_file():
        with np.load(path, allow_pickle=False) as npz:
            return "dadl_simulated" in npz.files and bool(npz["dadl_simulated"])
    if path.suffix == ".csv" and path.is_file():
        with path.open(encoding="utf-8-sig") as fh:
            return fh.readline().strip() == SIMULATED_CSV_MARK
    files = ptu_files(path)
    if not files:
        return False
    from dadl.adapters.ptu import PtuFile

    headers = [PtuFile(f).header for f in files]
    flags = {h.get("DADL_Simulated") is True or h.get("CreatorSW_Name") == "DADL" for h in headers}
    if len(flags) > 1:
        raise ConfigError(f"{path} mixes simulated and recorded PTU files; keep them in separate folders")
    return bool(files) and flags == {True}
