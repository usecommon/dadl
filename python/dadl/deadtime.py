"""Dead-time filters (pure-Python reference; the Rust core has the fast versions).

Timestamps are int64 picoseconds, sorted ascending. Unsorted input is refused:
a dead-time filter applied to unsorted tags gives a plausible-looking but
wrong answer.
"""

from __future__ import annotations

import numpy as np

from dadl import _checks


def apply_nonparalyzable(times: np.ndarray, tau_ps: int) -> np.ndarray:
    """Drop arrivals within `tau_ps` of the last recorded event."""
    tau_ps = _checks.integer(tau_ps, "dead time tau_ps", minimum=0)
    times = _checks.timestamps(times)
    _checks.check_sorted(times)
    if tau_ps == 0 or times.size == 0:
        return times.copy()
    out = np.empty_like(times)
    n = 0
    next_allowed: int | None = None
    for t in times.tolist():  # Python ints: t + tau never wraps
        if next_allowed is None or t >= next_allowed:
            out[n] = t
            n += 1
            next_allowed = t + tau_ps
    return out[:n]


def apply_paralyzable(times: np.ndarray, tau_ps: int) -> np.ndarray:
    """Record an arrival only if the previous arrival (recorded or not) is at least `tau_ps` earlier."""
    tau_ps = _checks.integer(tau_ps, "dead time tau_ps", minimum=0)
    times = _checks.timestamps(times)
    _checks.check_sorted(times)
    _checks.check_span(times)
    if times.size == 0:
        return times.copy()
    gaps = np.diff(times)
    keep = np.concatenate(([True], gaps >= tau_ps))
    return times[keep]
