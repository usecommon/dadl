"""Reading and writing CHSH time-tag files (phase 1 generic format).

A tag file is a table with three integer columns:

    setting        0..3, the CHSH setting pair, in the order
                   (a, b), (a, b'), (a', b), (a', b')
    channel        detector input number, mapped to A+, A-, B+, B- by the
                   `channels:` map in setup.yaml (default 0, 1, 2, 3)
    timestamp_ps   int64 time tag in picoseconds

Two encodings are accepted:

- `.npz`: a NumPy archive with arrays named `setting`, `channel`, `timestamp_ps`;
- `.csv`: a header line `setting,channel,timestamp_ps`, then one row per tag.

Vendor formats (PTU, Swabian .ttbin, Photon-HDF5) arrive as adapters in
phase 2; until then, export to one of these.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from dadl import _checks
from dadl.analysis.chsh import CHANNELS, SIGNS
from dadl.errors import ConfigError

COLUMNS = ("setting", "channel", "timestamp_ps")
SETTING_PAIRS: tuple[tuple[str, str], ...] = tuple(SIGNS)

Runs = dict[tuple[str, str], dict[str, np.ndarray]]


def read_tags(path: str | Path, channels: dict[int, str] | None = None) -> Runs:
    """Read a tag file into the `runs` structure that `dadl.chsh` takes."""
    path = Path(path)
    channels = channels or dict(enumerate(CHANNELS))
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
            table = np.loadtxt(path, delimiter=",", skiprows=1, dtype=np.int64, ndmin=2)
            cols = {c: table[:, i] if table.size else np.array([], dtype=np.int64) for i, c in enumerate(COLUMNS)}
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
    bad_s = np.setdiff1d(np.unique(setting), np.arange(len(SETTING_PAIRS)))
    if bad_s.size:
        raise ConfigError(f"{path}: setting values must be 0..3, found {bad_s.tolist()}")
    bad_c = np.setdiff1d(np.unique(channel), np.array(sorted(channels), dtype=np.int64))
    if bad_c.size:
        raise ConfigError(f"{path}: channel numbers {bad_c.tolist()} are not in the setup's channel map")

    runs: Runs = {}
    for i, pair in enumerate(SETTING_PAIRS):
        in_setting = setting == i
        runs[pair] = {
            name: np.sort(stamps[in_setting & (channel == num)]) for num, name in channels.items()
        }
    return runs


def write_tags(path: str | Path, runs: Runs, channels: dict[int, str] | None = None) -> Path:
    """Write `runs` to a .npz or .csv tag file (the inverse of `read_tags`)."""
    path = Path(path)
    channels = channels or dict(enumerate(CHANNELS))
    number = {name: num for num, name in channels.items()}
    parts = []
    for i, pair in enumerate(SETTING_PAIRS):
        for name, tags in runs[pair].items():
            tags = np.asarray(tags, dtype=np.int64)
            parts.append(np.stack([np.full(tags.size, i), np.full(tags.size, number[name]), tags], axis=1))
    table = np.concatenate(parts) if parts else np.empty((0, 3), dtype=np.int64)
    table = table[np.lexsort((table[:, 2], table[:, 1], table[:, 0]))]
    if path.suffix == ".npz":
        np.savez_compressed(path, **{c: table[:, i] for i, c in enumerate(COLUMNS)})
    elif path.suffix == ".csv":
        np.savetxt(path, table, fmt="%d", delimiter=",", header=",".join(COLUMNS), comments="")
    else:
        raise ConfigError(f"{path}: tag files must be .npz or .csv")
    return path
