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


def interruption_pvalue(
    tags: np.ndarray, start: int | None = None, end: int | None = None
) -> tuple[float, float]:
    """(p, longest gap in ps): p that the longest gap is this long in an uninterrupted stream.

    `start` and `end` bound the acquisition (e.g. the earliest and latest tag over all its
    channels): the gaps from `start` to the first tag and from the last tag to `end` count too, so
    a channel that comes up late or goes dark early is caught (QKD review, cycle 4: a trailing 30%
    latch read heralding 0.457 for 0.653 with status ok).
    """
    tags = np.asarray(tags)
    if tags.size and start is not None and start < int(tags[0]):
        tags = np.concatenate([[start], tags])
    if tags.size and end is not None and end > int(tags[-1]):
        tags = np.concatenate([tags, [end]])
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
    """Bonferroni-combined p over several channels, and the longest gap seen (ps).

    The streams are taken as one acquisition: its bounds are the earliest and latest tag of all.
    """
    nonempty = [np.asarray(s) for s in streams if np.asarray(s).size]
    start = min((int(s[0]) for s in nonempty), default=None)
    end = max((int(s[-1]) for s in nonempty), default=None)
    results = [interruption_pvalue(s, start, end) for s in streams]
    if not results:
        return 1.0, 0.0
    p = min(1.0, len(results) * min(r[0] for r in results))
    return p, max(r[1] for r in results)


def interrupted(*acquisitions: list[np.ndarray], alpha: float = 0.001) -> bool:
    """True when any channel of any acquisition fails the interruption test.

    Each acquisition is a list of its channels' tag arrays, bounded by its own earliest and latest
    tag; Bonferroni over all channels at `alpha`.
    """
    n = sum(len(acq) for acq in acquisitions)
    if n == 0:
        return False
    worst = min(channels_pvalue(*acq)[0] / max(len(acq), 1) for acq in acquisitions if len(acq))
    return min(1.0, n * worst) < alpha
