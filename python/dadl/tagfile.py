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

A third form is a folder of PicoQuant PTU files (T2), one per acquisition,
named setting0.ptu, setting1.ptu, ... in setting order; PTU input numbers are
the channel numbers. Other vendor formats (Swabian .ttbin, Photon-HDF5) need
exporting to one of these.
"""

from __future__ import annotations

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


def read_tags(path: str | Path, channels: dict[int, str] | None = None, settings: Sequence = SETTING_PAIRS) -> Runs:
    """Read a tag file into the `runs` structure that `dadl.chsh` takes.

    `settings` names the setting indices in order (CHSH: the four setting
    pairs; HOM: "dip", "reference").
    """
    path = Path(path)
    channels = channels or dict(enumerate(CHANNELS))
    settings = tuple(settings)
    if path.is_dir():
        return _read_ptu_dir(path, channels, settings)
    try:
        if path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as npz:
                missing = [c for c in COLUMNS if c not in npz.files]
                if missing:
                    raise ConfigError(f"{path}: missing array(s) {missing}; need {list(COLUMNS)}")
                cols = {c: npz[c] for c in COLUMNS}
        elif path.suffix == ".csv":
            with path.open(encoding="utf-8") as fh:
                header = [h.strip() for h in fh.readline().split(",")]
            if header != list(COLUMNS):
                raise ConfigError(f"{path}: header must be {','.join(COLUMNS)}, got {','.join(header)}")
            # Read as text and check every cell is an integer: some NumPy versions
            # silently truncate "1.5" when asked for int64 directly.
            text = np.loadtxt(path, delimiter=",", skiprows=1, dtype=str, ndmin=2)
            if text.size:
                text = np.char.strip(text)
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
            raise ConfigError(f"{path}: tag files must be .npz or .csv")
    except FileNotFoundError as exc:
        raise ConfigError(f"tag file not found: {path}") from exc
    except ValueError as exc:
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
        runs[pair] = {
            name: np.sort(stamps[in_setting & (channel == num)]) for num, name in channels.items()
        }
    return runs


def write_tags(
    path: str | Path,
    runs: Runs,
    channels: dict[int, str] | None = None,
    settings: Sequence = SETTING_PAIRS,
    *,
    simulated: bool = False,
) -> Path:
    """Write `runs` to a .npz or .csv tag file (the inverse of `read_tags`)."""
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
        np.savetxt(path, table, fmt="%d", delimiter=",", header=",".join(COLUMNS), comments="")
    else:
        raise ConfigError(f"{path}: tag files must be .npz or .csv")
    return path


def ptu_file(directory: str | Path, index: int) -> Path:
    """Where acquisition `index` lives in a folder of PTU files."""
    return Path(directory) / f"setting{index}.ptu"


def _read_ptu_dir(directory: Path, channels: dict[int, str], settings: tuple) -> Runs:
    from dadl.adapters.ptu import PtuFile
    from dadl.errors import AdapterError

    runs: Runs = {}
    for i, key in enumerate(settings):
        path = ptu_file(directory, i)
        if not path.exists():
            raise ConfigError(f"{directory}: missing {path.name} (one PTU file per acquisition, setting0.ptu, ...)")
        try:
            tagset = PtuFile(path).read()
        except AdapterError as exc:
            raise ConfigError(str(exc)) from exc
        if tagset.metadata.get("incomplete"):
            raise ConfigError(
                f"{path}: the acquisition lost data ({tagset.metadata.get('data_loss')}, from TTResult_StopReason); "
                "DADL does not analyse incomplete files. Read it with dadl.adapters.PtuFile to inspect it."
            )
        by_input = {ch.input: tagset[ch] for ch in tagset.channels}
        unknown = sorted(set(by_input) - set(channels))
        if unknown:
            raise ConfigError(
                f"{path}: inputs {unknown} are not in the setup's channel map. On HydraHarp, MultiHarp and "
                "TimeHarp 260 files input 0 is the sync input and detector inputs are numbered from 1 (as in "
                "PicoQuant's software); map them under `channels:` in setup.yaml"
            )
        runs[key] = {name: by_input.get(num, np.empty(0, dtype=np.int64)) for num, name in channels.items()}
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
    """The tagger tick in whole ps for a folder of PTU files (from the header), else None."""
    path = Path(path)
    if not path.is_dir():
        return None
    from dadl.adapters.ptu import PtuFile

    files = sorted(path.glob("setting*.ptu"))
    if not files:
        return None
    res = PtuFile(files[0]).resolution_s * 1e12
    whole = round(res)
    return whole if whole >= 1 and abs(res - whole) <= 1e-9 * res else None


def tags_are_simulated(path: str | Path) -> bool:
    """True when the tags were written by DADL's simulator: a folder of PTU files whose headers say
    so, or an .npz written by `write_tags` with simulated=True.

    A folder mixing simulated and recorded files is refused: its provenance cannot be one label.
    """
    path = Path(path)
    if path.suffix == ".npz" and path.is_file():
        with np.load(path, allow_pickle=False) as npz:
            return "dadl_simulated" in npz.files and bool(npz["dadl_simulated"])
    if not path.is_dir():
        return False
    from dadl.adapters.ptu import PtuFile

    files = sorted(path.glob("setting*.ptu"))
    flags = {PtuFile(f).header.get("DADL_Simulated") is True or PtuFile(f).header.get("CreatorSW_Name") == "DADL" for f in files}
    if len(flags) > 1:
        raise ConfigError(f"{path} mixes simulated and recorded PTU files; keep them in separate folders")
    return bool(files) and flags == {True}
