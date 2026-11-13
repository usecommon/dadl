"""A simulated live time tagger: the device half of DADL's hardware twin.

`SimulatedTagger` behaves like a live tagger behind the `Adapter` protocol.
It streams batches of tags from the physics simulator (`dadl.simsetup`),
quantized to a PicoQuant device's tick, and injects the faults a real stream
has, each reported as a typed `Event`:

- `buffer-overflow`: the tagger's buffer filled, and every tag in a time
  window was lost (`value` = tags lost, `detail` = window and per-input counts);
- `dropped-batch`: a whole batch never reached the host;
- `config-change`: a setting changed mid-run (recorded as provenance).

Nothing is ever lost silently: for every input, tags delivered plus tags
reported lost equal the tags generated (`truth()`), which is the accounting
rule every adapter must meet (adapters/TEMPLATE.md, rule 3) and what
`check_accounting` verifies. Everything it returns is marked simulated.
"""

from __future__ import annotations

import time
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from dadl import __version__
from dadl.adapters import Capabilities, DeviceInfo, Event
from dadl.adapters.ptu_write import DEVICES, MAX_INPUT
from dadl.config import Setup
from dadl.data import ChannelId, TagSet
from dadl.simsetup import simulate_setup

PS_PER_S = 10**12


@dataclass(frozen=True)
class Faults:
    """Faults to inject. Times are seconds from the start of the acquisition."""

    #: (start_s, length_s) windows in which every tag is lost.
    buffer_overflows: tuple[tuple[float, float], ...] = ()
    #: Indices of batches that never arrive.
    dropped_batches: tuple[int, ...] = ()
    #: (time_s, key, new value) configuration changes.
    config_changes: tuple[tuple[float, str, Any], ...] = ()


class SimulatedTagger:
    """A live tagger emulated from a setup file (one acquisition), with injectable faults."""

    def __init__(
        self,
        setup: Setup,
        acquisition: Any | None = None,
        *,
        device: str = "hydraharp",
        batch_s: float = 0.1,
        faults: Faults = Faults(),
        realtime: bool = False,
    ) -> None:
        if device not in DEVICES:
            raise ValueError(f"device must be one of {sorted(DEVICES)}, got {device!r}")
        if not batch_s > 0:
            raise ValueError("batch_s must be positive")
        spans = sorted(faults.buffer_overflows)
        if any(d <= 0 for _, d in spans) or any(s0 + d0 > s1 for (s0, d0), (s1, _) in zip(spans, spans[1:])):
            raise ValueError("buffer_overflows must have positive lengths and must not overlap")
        key = setup.setting_keys[0] if acquisition is None else acquisition
        if key not in setup.setting_keys:
            raise ValueError(f"acquisition must be one of {list(setup.setting_keys)}, got {key!r}")
        _, self._tick_ps, hw_type = DEVICES[device]
        run = simulate_setup(setup)[key]
        number = {name: num for num, name in setup.channels.items()}
        self._truth = {
            number[name]: (np.asarray(arr, dtype=np.int64) // self._tick_ps) * self._tick_ps for name, arr in run.items()
        }
        first = min((int(a[0]) for a in self._truth.values() if a.size), default=0)
        if first < 0:  # jitter near t = 0: start the clock at the earliest tag
            self._truth = {k: a - first for k, a in self._truth.items()}
        self._duration_ps = int(round(setup.duration_s_per_setting * PS_PER_S))
        self._batch_ps = int(round(batch_s * PS_PER_S))
        # Twin review, cycle 9: the PTU writer enforced each device's inputs, the live twin did not;
        # faults outside the acquisition removed tags without an event, or were ignored.
        too_high = sorted(k for k in self._truth if k > MAX_INPUT[device])
        if too_high:
            raise ValueError(f"{hw_type} has inputs 0-{MAX_INPUT[device]}; the setup maps inputs {too_high}")
        duration_s = self._duration_ps / PS_PER_S
        if any(s0 < 0 or s0 + d0 > duration_s for s0, d0 in spans):
            raise ValueError(f"buffer_overflows must lie within the {duration_s:g} s acquisition")
        n_batches = max(1, -(-self._duration_ps // self._batch_ps))
        if any(not 0 <= b < n_batches for b in faults.dropped_batches):
            raise ValueError(f"dropped_batches must be batch indices 0-{n_batches - 1}")
        if any(not 0 <= t < duration_s for t, _, _ in faults.config_changes):
            raise ValueError(f"config_changes must fall within the {duration_s:g} s acquisition")
        self._faults = faults
        self._realtime = realtime
        self._setup_name = setup.name
        self._acquisition = key
        self._info = DeviceInfo(
            model=f"{hw_type} (simulated)",
            serial="SIM-0001",
            firmware=f"dadl {__version__}",
            channels=tuple(sorted(self._truth)),
            resolution_ps=float(self._tick_ps),
        )
        self._config: dict[str, Any] = {
            "simulated": True,
            "setup": setup.name,
            "acquisition": str(key),
            **{f"HWInpChan_TrgLevel({k})": -50 for k in sorted(self._truth)},
        }
        self._events: list[Event] = []

    # --- Adapter protocol --------------------------------------------------------------------

    def describe(self) -> DeviceInfo:
        return self._info

    def capabilities(self) -> Capabilities:
        return Capabilities(live=True, replay=True)

    def config(self) -> Mapping[str, Any]:
        """The configuration at the start of the run; changes arrive as `config-change` events."""
        return dict(self._config)

    def events(self) -> tuple[Event, ...]:
        return tuple(self._events)

    def read(self) -> TagSet:
        """Run the whole acquisition and return every delivered tag."""
        parts: dict[int, list[np.ndarray]] = {k: [] for k in self._truth}
        for batch in self.stream():
            for ch in batch.channels:
                parts[ch.input].append(batch[ch])
        tags = {
            self._channel(k): (np.concatenate(v) if v else np.empty(0, dtype=np.int64)) for k, v in parts.items()
        }
        lost = [e for e in self._events if e.kind in ("buffer-overflow", "dropped-batch")]
        metadata = self._metadata()
        if lost:
            # As a PTU file stopped by data loss is (twin review, cycle 9).
            metadata.update({"incomplete": True, "data_loss": ", ".join(sorted({e.kind for e in lost}))})
        return TagSet(tags, duration_ps=self._duration_ps, clock=self._clock(), metadata=metadata)

    # --- live streaming ----------------------------------------------------------------------

    def stream(self) -> Iterator[TagSet]:
        """Batches of `batch_s` in time order; faults remove tags and add events."""
        self._events = []
        overflow_lost = {k: np.zeros(a.size, dtype=bool) for k, a in self._truth.items()}
        windows = [(int(s * PS_PER_S), int((s + d) * PS_PER_S)) for s, d in self._faults.buffer_overflows]
        changes = sorted(self._faults.config_changes, key=lambda c: c[0])
        n_batches = max(1, -(-self._duration_ps // self._batch_ps))
        dropped = set(self._faults.dropped_batches)
        for b in range(n_batches):
            lo, hi = b * self._batch_ps, min((b + 1) * self._batch_ps, self._duration_ps)
            if b == n_batches - 1:
                hi = max(hi, max((int(a[-1]) + 1 for a in self._truth.values() if a.size), default=hi))
            for t_s, key, value in changes:
                if lo <= int(t_s * PS_PER_S) < hi:
                    self._events.append(
                        Event("config-change", time_ps=int(t_s * PS_PER_S), detail={"key": key, "value": value})
                    )
            in_batch = {k: (a >= lo) & (a < hi) for k, a in self._truth.items()}
            for w_lo, w_hi in windows:
                if lo <= w_lo < hi:
                    # Report each overflow once, in the batch where it starts, counting
                    # every tag it removed outside dropped batches.
                    per_input = {}
                    for k, a in self._truth.items():
                        in_w = (a >= w_lo) & (a < w_hi)
                        batch_of = np.minimum(a // self._batch_ps, n_batches - 1)
                        not_dropped = ~np.isin(batch_of, list(dropped)) if dropped else True
                        per_input[k] = int(np.sum(in_w & not_dropped))
                    self._events.append(
                        Event("buffer-overflow", time_ps=w_lo, value=sum(per_input.values()),
                              detail={"until_ps": w_hi, "per_input": per_input})
                    )
            if b in dropped:
                per_input = {k: int(m.sum()) for k, m in in_batch.items()}
                self._events.append(
                    Event("dropped-batch", time_ps=lo, value=sum(per_input.values()),
                          detail={"until_ps": hi, "per_input": per_input, "batch": b})
                )
                continue
            batch_tags = {}
            for k, a in self._truth.items():
                keep = in_batch[k].copy()
                for w_lo, w_hi in windows:
                    lost = keep & (a >= w_lo) & (a < w_hi)
                    overflow_lost[k] |= lost
                    keep &= ~lost
                batch_tags[self._channel(k)] = a[keep]
            if self._realtime:
                time.sleep((hi - lo) / PS_PER_S)
            yield TagSet(batch_tags, duration_ps=max(1, hi - lo), clock=self._clock(),
                         metadata={**self._metadata(), "batch": b, "start_ps": lo})

    # --- for tests: what really happened -----------------------------------------------------

    def truth(self) -> dict[int, np.ndarray]:
        """Every tag generated, per input, before faults (not available from a real device)."""
        return {k: a.copy() for k, a in self._truth.items()}

    # --- helpers -----------------------------------------------------------------------------

    def _channel(self, inp: int) -> ChannelId:
        return ChannelId(f"{self._info.model.split(' (')[0]}-SIM-0001", inp)

    def _clock(self) -> str:
        return f"{self._info.model} clock"

    def _metadata(self) -> dict[str, Any]:
        return {"source": "simulated-tagger", "simulated": True, "setup": self._setup_name,
                "acquisition": str(self._acquisition), "resolution_ps": float(self._tick_ps)}


def check_accounting(delivered: TagSet, events: tuple[Event, ...], generated: Mapping[int, np.ndarray]) -> dict[int, int]:
    """For each input, the number of tags not accounted for; all zeros means nothing was lost silently.

    Compares times, not just counts: every delivered tag must be a generated tag (no duplicates,
    shifts or extra inputs), and every generated tag that was not delivered must fall inside a
    reported loss (a buffer overflow window or a dropped batch). Needs the generated truth, so it
    applies to the simulated tagger; for a real device the adapter's rule is record accounting
    plus the vendor's own loss flags (adapters/TEMPLATE.md).
    """
    devices = {ch.device for ch in delivered.channels}
    if len(devices) > 1:
        raise ValueError(f"delivered tags come from several devices {sorted(devices)}; check each device separately")
    windows = [
        (int(e.time_ps), int(e.detail["until_ps"]))
        for e in events
        if e.kind in ("buffer-overflow", "dropped-batch") and e.time_ps is not None and "until_ps" in e.detail
    ]
    got = {ch.input: np.asarray(delivered[ch]) for ch in delivered.channels}
    out: dict[int, int] = {}
    for k in sorted(set(generated) | set(got)):
        truth = np.sort(np.asarray(generated.get(k, np.empty(0, dtype=np.int64))))
        seen = np.sort(got.get(k, np.empty(0, dtype=np.int64)))
        values, t_count = np.unique(truth, return_counts=True)
        s_values, s_count = np.unique(seen, return_counts=True)
        idx = np.searchsorted(values, s_values)
        idx_ok = idx < values.size
        match = np.zeros(s_values.size, dtype=bool)
        match[idx_ok] = values[idx[idx_ok]] == s_values[idx_ok]
        # Delivered tags that were never generated, or delivered more often than generated.
        extra = int(s_count[~match].sum())
        over = np.zeros(values.size, dtype=np.int64)
        over[idx[match]] = s_count[match]
        extra += int(np.clip(over - t_count, 0, None).sum())
        missing_counts = np.clip(t_count - over, 0, None)
        missing = np.repeat(values, missing_counts)
        covered = np.zeros(missing.size, dtype=bool)
        for lo, hi in windows:
            covered |= (missing >= lo) & (missing < hi)
        out[k] = extra + int((~covered).sum())
    return out
