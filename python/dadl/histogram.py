"""Histograms used by calibration fits (reference for core/src/histogram.rs)."""

from __future__ import annotations

import numpy as np

from dadl import _checks


def interarrival(times: np.ndarray, bin_width_ps: int, n_bins: int) -> np.ndarray:
    """Histogram of gaps between consecutive events on one channel.

    Bin k counts gaps in [k * bin_width_ps, (k + 1) * bin_width_ps). Gaps
    beyond the last bin are dropped. `times` must be sorted ascending.
    """
    bin_width_ps = _checks.integer(bin_width_ps, "bin_width_ps", minimum=1)
    n_bins = _checks.integer(n_bins, "n_bins", minimum=0)
    times = _checks.timestamps(times)
    _checks.check_sorted(times)
    _checks.check_span(times)
    hist = np.zeros(n_bins, dtype=np.uint64)
    if times.size < 2 or n_bins == 0:
        return hist
    k = np.diff(times) // bin_width_ps
    k = k[k < n_bins]
    np.add.at(hist, k, 1)
    return hist


#: Pairs processed per chunk in `cross_correlation`, to bound memory.
_CHUNK_PAIRS = 1 << 22


def cross_correlation(a: np.ndarray, b: np.ndarray, bin_width_ps: int, range_ps: int) -> np.ndarray:
    """Histogram of delays b[j] - a[i] over all pairs of tags within `range_ps`.

    Delays d with -range_ps <= d < range_ps are counted. Bin k covers
    [-range_ps + k * bin_width_ps, -range_ps + (k + 1) * bin_width_ps); the
    last of the ceil(2 * range_ps / bin_width_ps) bins is cut at range_ps.
    Both inputs must be sorted ascending. Reference for
    core/src/histogram.rs, held to the same vectors.
    """
    bin_width_ps = _checks.integer(bin_width_ps, "bin_width_ps")
    range_ps = _checks.integer(range_ps, "range_ps")
    if bin_width_ps <= 0 or range_ps <= 0:
        raise ValueError("bin_width_ps and range_ps must be positive")
    a = _checks.timestamps(a, "a")
    b = _checks.timestamps(b, "b")
    _checks.check_sorted(a, "a")
    _checks.check_sorted(b, "b")
    span = 2 * range_ps
    if span > _checks.INT64_MAX:
        raise OverflowError("2 * range_ps exceeds the int64 range")
    if b.size and (int(b.min()) - range_ps < _checks.INT64_MIN or int(b.max()) + range_ps > _checks.INT64_MAX):
        raise OverflowError("a tag plus or minus range_ps leaves the int64 range")
    n_bins = (span - 1) // bin_width_ps + 1
    hist = np.zeros(n_bins, dtype=np.uint64)
    if a.size == 0 or b.size == 0:
        return hist
    r = np.int64(range_ps)
    lo = np.searchsorted(a, b - r, side="right")  # first a with b - a < range
    hi = np.searchsorted(a, b + r, side="right")  # past the last a with b - a >= -range
    counts = hi - lo
    ends = np.cumsum(counts)
    j0 = 0
    while j0 < b.size:
        # Take as many b tags as fit in one chunk of pairs (at least one).
        base = ends[j0 - 1] if j0 else 0
        j1 = max(j0 + 1, int(np.searchsorted(ends, base + _CHUNK_PAIRS, side="right")))
        c = counts[j0:j1]
        total = int(c.sum())
        if total:
            jb = np.repeat(np.arange(j0, j1), c)
            offsets = np.arange(total) - np.repeat(np.cumsum(c) - c, c)
            ia = lo[jb] + offsets
            k = (b[jb] - a[ia] + r) // np.int64(bin_width_ps)
            hist += np.bincount(k, minlength=n_bins).astype(np.uint64)
        j0 = j1
    return hist
