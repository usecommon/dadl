"""Checks that one channel's acquisition was uninterrupted.

A latched SNSPD, an overload shutdown or a dropped stretch of data leaves one gap between
consecutive tags far longer than the rest. The SNSPD specialist's review found a latch in 20% of a
HOM dip run raising the visibility by 13 sigma, and a 10% latch on a heralding signal channel
reading eta 0.657 against 0.730, both reported as normal.

The test compares the longest gap with the gaps above the median: past the median gap, dead time
and afterpulsing are over (for any detector this package models) and the gaps of a stationary
Poisson stream are exponential, so the excess over the median has a rate estimated from the data.
This makes the test robust to dead time and afterpulsing, which the plain rate n/T is not.
"""

from __future__ import annotations

import math

import numpy as np

#: Below this many gaps the test has no power and returns p = 1.
MIN_GAPS = 20


def interruption_pvalue(tags: np.ndarray) -> tuple[float, float]:
    """(p, longest gap in ps): p that the longest gap is this long in an uninterrupted stream."""
    tags = np.asarray(tags)
    if tags.size < MIN_GAPS + 1:
        return 1.0, 0.0
    gaps = np.diff(tags).astype(np.float64)
    med = float(np.median(gaps))
    excess = gaps[gaps > med] - med
    longest = float(gaps.max())
    if excess.size < MIN_GAPS // 2 or float(excess.mean()) <= 0:
        return 1.0, longest
    # Rate from the excesses without the largest one, so a latch does not dilute its own test.
    rest = np.sort(excess)[:-1]
    rate = 1.0 / float(rest.mean())
    m = int(excess.size)
    g = longest - med
    p = float(-np.expm1(m * np.log1p(-math.exp(-rate * g)))) if g > 0 else 1.0
    return min(1.0, max(p, 0.0)), longest


def channels_pvalue(*streams: np.ndarray) -> tuple[float, float]:
    """Bonferroni-combined p over several channels, and the longest gap seen (ps)."""
    results = [interruption_pvalue(s) for s in streams]
    if not results:
        return 1.0, 0.0
    p = min(1.0, len(results) * min(r[0] for r in results))
    return p, max(r[1] for r in results)


def interrupted(*streams: np.ndarray, alpha: float = 0.001) -> bool:
    """True when any stream fails the interruption test (Bonferroni over streams at `alpha`)."""
    return channels_pvalue(*streams)[0] < alpha
