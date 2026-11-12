"""Dataset types: channels, tag sets and datasets.

These are thin, validated containers. A `TagSet` holds the recorded tags of
one acquisition, per channel, as sorted int64 picoseconds on one declared
clock. A `Dataset` groups acquisitions with the setup metadata they need.
Per-tag work stays in the Rust core and NumPy; these types never hold
per-tag Python objects.
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any

import numpy as np

from dadl import _checks


@dataclass(frozen=True, order=True)
class ChannelId:
    """One detector input: the device (serial or simulator ID) and the tagger input number."""

    device: str
    input: int

    def __post_init__(self) -> None:
        if not isinstance(self.device, str) or not self.device.strip():
            raise ValueError("device must be a non-empty string")
        object.__setattr__(self, "input", _checks.integer(self.input, "input", minimum=0))

    def as_record(self) -> dict[str, object]:
        return {"device": self.device, "input": self.input}

    def __str__(self) -> str:
        return f"{self.device}:{self.input}"


def tags_sha256(tags: np.ndarray) -> str:
    """SHA-256 of a tag array as little-endian int64 bytes (the hash recorded in artifacts)."""
    arr = np.ascontiguousarray(_checks.timestamps(tags), dtype="<i8")
    return hashlib.sha256(arr.tobytes()).hexdigest()


@dataclass(frozen=True)
class TagSet:
    """The recorded tags of one acquisition, per channel, sorted int64 ps on one clock."""

    tags: Mapping[ChannelId, np.ndarray]
    #: Acquisition length in ps, if known; otherwise the span of the tags is used.
    duration_ps: int | None = None
    #: Which clock the tags are on (e.g. "tagger", "sim").
    clock: str = "tagger"
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        checked: dict[ChannelId, np.ndarray] = {}
        for ch, arr in self.tags.items():
            if not isinstance(ch, ChannelId):
                raise TypeError(f"TagSet keys must be ChannelId, got {ch!r}")
            arr = _checks.timestamps(arr, str(ch))
            _checks.check_sorted(arr, str(ch))
            _checks.check_span(arr, str(ch))
            arr = arr.copy()
            arr.setflags(write=False)
            checked[ch] = arr
        object.__setattr__(self, "tags", MappingProxyType(checked))
        if self.duration_ps is not None:
            object.__setattr__(self, "duration_ps", _checks.integer(self.duration_ps, "duration_ps", minimum=1))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def channels(self) -> tuple[ChannelId, ...]:
        return tuple(sorted(self.tags))

    def __getitem__(self, ch: ChannelId) -> np.ndarray:
        return self.tags[ch]

    def duration(self, ch: ChannelId | None = None) -> int:
        """Declared duration, or the span of channel `ch` (or of all channels)."""
        if self.duration_ps is not None:
            return self.duration_ps
        arrays = [self.tags[ch]] if ch is not None else [a for a in self.tags.values() if a.size]
        if not arrays or all(a.size < 2 for a in arrays):
            raise ValueError("duration is unknown: declare duration_ps or provide at least two tags")
        return int(max(int(a[-1]) for a in arrays) - min(int(a[0]) for a in arrays))


@dataclass(frozen=True)
class Dataset:
    """Acquisitions plus the setup they came from."""

    acquisitions: Mapping[str, TagSet]
    setup_name: str | None = None
    backend: str = "hardware"

    def __post_init__(self) -> None:
        if self.backend not in ("sim", "hardware"):
            raise ValueError(f"backend must be 'sim' or 'hardware', got {self.backend!r}")
        object.__setattr__(self, "acquisitions", MappingProxyType(dict(self.acquisitions)))
