"""Fixed-grid (clock-locked) coincidence counting.

Time is cut into slots of a width declared before the data is seen. A
coincidence is a slot occupied by both channels after each channel's declared
delay is removed. This is the only coincidence method allowed for Bell-type
analyses: windows centred on detection events open the coincidence-time
loophole.

Inputs need not be sorted. Widths and delays are integer picoseconds; a
delay that would move a tag outside the int64 range is refused.
"""

from __future__ import annotations

import numpy as np

from dadl import _backend, _checks


def occupied_slots(times: np.ndarray, width_ps: int, delay_ps: int = 0) -> np.ndarray:
    """Sorted, de-duplicated slot indices occupied by `times` after removing `delay_ps`."""
    core = _backend.core()
    if core is not None:
        return core.occupied_slots(
            _checks.timestamps(times), _checks.integer(width_ps, "width_ps"), _checks.integer(delay_ps, "delay_ps")
        )
    width_ps = _checks.integer(width_ps, "width_ps", minimum=1)
    delay_ps = _checks.integer(delay_ps, "delay_ps")
    times = _checks.timestamps(times)
    if times.size:
        lo, hi = int(times.min()) - delay_ps, int(times.max()) - delay_ps
        if lo < _checks.INT64_MIN or hi > _checks.INT64_MAX:
            raise OverflowError("removing the declared delay moves a tag outside the int64 range")
    return np.unique(np.floor_divide(times - np.int64(delay_ps), np.int64(width_ps)))


def grid_coincidences(
    a: np.ndarray,
    b: np.ndarray,
    width_ps: int,
    delay_a_ps: int = 0,
    delay_b_ps: int = 0,
) -> int:
    """Number of slots occupied by both channels."""
    core = _backend.core()
    if core is not None:
        return int(
            core.grid_coincidences(
                _checks.timestamps(a, "a"),
                _checks.timestamps(b, "b"),
                _checks.integer(width_ps, "width_ps"),
                _checks.integer(delay_a_ps, "delay_a_ps"),
                _checks.integer(delay_b_ps, "delay_b_ps"),
            )
        )
    slots_a = occupied_slots(a, width_ps, delay_a_ps)
    slots_b = occupied_slots(b, width_ps, delay_b_ps)
    return int(np.intersect1d(slots_a, slots_b, assume_unique=True).size)
