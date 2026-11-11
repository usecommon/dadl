"""Fixed-grid (clock-locked) coincidence counting.

Time is cut into slots of a width declared before the data is seen. A
coincidence is a slot occupied by both channels after each channel's declared
delay is removed. This is the only coincidence method allowed for Bell-type
analyses: windows centred on detection events open the coincidence-time
loophole.
"""

from __future__ import annotations

import numpy as np


def occupied_slots(times: np.ndarray, width_ps: int, delay_ps: int = 0) -> np.ndarray:
    if width_ps <= 0:
        raise ValueError("slot width must be positive")
    times = np.asarray(times, dtype=np.int64)
    return np.unique(np.floor_divide(times - np.int64(delay_ps), np.int64(width_ps)))


def grid_coincidences(
    a: np.ndarray,
    b: np.ndarray,
    width_ps: int,
    delay_a_ps: int = 0,
    delay_b_ps: int = 0,
) -> int:
    """Number of slots occupied by both channels."""
    slots_a = occupied_slots(a, width_ps, delay_a_ps)
    slots_b = occupied_slots(b, width_ps, delay_b_ps)
    return int(np.intersect1d(slots_a, slots_b, assume_unique=True).size)
