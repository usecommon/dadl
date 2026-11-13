"""Tagger adapters: vendor files and devices to `TagSet`s.

This is the Python subset of the Rust `Adapter` trait in docs/design.md
("Adapter contract"). Every adapter:

- converts timestamps to int64 ps on a declared clock and handles overflow
  and rollover itself;
- reports buffer overflows, markers and anything else it does not turn into
  tags as typed `Event`s, and accounts for every record it reads (a silent
  drop is a test failure);
- declares its capabilities rather than letting callers assume them;
- records the device configuration it saw, as provenance.

Phase 2 ships one file adapter, `PtuFile` (PicoQuant PTU, T2 mode), and a simulated live device, `SimulatedTagger`. Writing a new
one: see adapters/TEMPLATE.md in the repository.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any, Protocol, runtime_checkable

from dadl.data import TagSet


@dataclass(frozen=True)
class DeviceInfo:
    """What the adapter is reading from."""

    model: str
    serial: str | None
    firmware: str | None
    #: Input numbers that can carry tags.
    channels: tuple[int, ...]
    #: Timestamp resolution in ps (the size of one tick).
    resolution_ps: float


@dataclass(frozen=True)
class Capabilities:
    """Declared, not assumed: what the core may ask of this adapter."""

    live: bool = False
    replay: bool = True
    vendor_histograms: bool = False
    settable_delays: bool = False
    gating: bool = False


@dataclass(frozen=True)
class Event:
    """Something in the stream that is not a tag: an overflow, a marker, a dropped batch."""

    kind: str
    #: Time on the adapter's clock in ps, when the event has one.
    time_ps: int | None = None
    #: Kind-specific value (e.g. marker bits, number of overflows).
    value: int = 0
    detail: Mapping[str, Any] = field(default_factory=dict)


@runtime_checkable
class Adapter(Protocol):
    """Python Protocol matching the Rust trait (read-only subset)."""

    def describe(self) -> DeviceInfo: ...

    def capabilities(self) -> Capabilities: ...

    def config(self) -> Mapping[str, Any]:
        """Device configuration as recorded (trigger levels, input offsets, ...), for provenance."""
        ...

    def read(self) -> TagSet:
        """All tags, as int64 ps per channel on the declared clock."""
        ...

    def events(self) -> tuple[Event, ...]:
        """Typed events met while reading (empty before `read`)."""
        ...


from dadl.adapters.ptu import PtuFile, read_ptu  # noqa: E402


def __getattr__(name: str):  # SimulatedTagger loads the simulator, so import it on first use
    if name in ("SimulatedTagger", "Faults", "check_accounting"):
        from dadl.adapters import simdevice

        return getattr(simdevice, name)
    raise AttributeError(name)


__all__ = [
    "Adapter",
    "Capabilities",
    "DeviceInfo",
    "Event",
    "Faults",
    "PtuFile",
    "SimulatedTagger",
    "check_accounting",
    "read_ptu",
]
