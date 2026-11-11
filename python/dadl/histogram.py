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
