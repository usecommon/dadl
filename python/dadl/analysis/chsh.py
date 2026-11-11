"""CHSH from fixed-grid coincidences.

Headline S uses raw counts only. Every coincidence-normalized S is labelled
"fair sampling assumed", because normalizing by coincidences post-selects on
detection.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from dadl.coincidences import grid_coincidences
from dadl.errors import InsufficientCounts
from dadl.guardrails import Policy
from dadl.result import Estimate, LabelledValue

FAIR_SAMPLING = "fair sampling assumed (post-selected on coincidences)"

#: Below this many coincidences per setting the Gaussian error is not trusted.
MIN_COUNTS_PER_SETTING = 100

SIGNS: dict[tuple[str, str], int] = {("a", "b"): 1, ("a", "b'"): -1, ("a'", "b"): 1, ("a'", "b'"): 1}


@dataclass(frozen=True)
class Window:
    """A coincidence window. `source` must be "plan" for Bell-type analyses."""

    width_ps: int
    source: str = "plan"
    delays_ps: dict[str, int] = field(default_factory=dict)


def correlator(tags: dict[str, np.ndarray], window: Window) -> tuple[float, float, int]:
    """E = (N++ + N-- - N+- - N-+) / N and its variance (1 - E^2) / N."""
    d = window.delays_ps

    def n(a: str, b: str) -> int:
        return grid_coincidences(tags[a], tags[b], window.width_ps, d.get(a, 0), d.get(b, 0))

    npp, nmm, npm, nmp = n("A+", "B+"), n("A-", "B-"), n("A+", "B-"), n("A-", "B+")
    total = npp + nmm + npm + nmp
    if total < MIN_COUNTS_PER_SETTING:
        raise InsufficientCounts(f"{total} coincidences in a setting; need at least {MIN_COUNTS_PER_SETTING}.")
    e = (npp + nmm - npm - nmp) / total
    return e, (1.0 - e * e) / total, total


def chsh(
    runs: dict[tuple[str, str], dict[str, np.ndarray]],
    window: Window,
    policy: Policy | None = None,
) -> LabelledValue:
    """Raw CHSH S from one run per setting pair, with fixed windows from the analysis plan."""
    policy = policy or Policy.default()
    policy.check_window("chsh", window.source)
    s = 0.0
    var = 0.0
    for key, sign in SIGNS.items():
        if key not in runs:
            raise ValueError(f"missing run for setting pair {key}")
        e, v, _ = correlator(runs[key], window)
        s += sign * e
        var += v
    labels = (FAIR_SAMPLING,) + policy.labels_for("chsh")
    return LabelledValue(analysis="chsh", raw=Estimate(s, math.sqrt(var)), labels=labels)
