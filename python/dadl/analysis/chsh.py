"""CHSH from fixed-grid coincidences.

Headline S uses raw counts only. Every coincidence-normalized S is labelled
"fair sampling assumed", because normalizing by coincidences post-selects on
detection.

Phase 1 limits, stated so nobody over-reads the result:
- the window's `source` is declared by the caller, not verified; signed
  analysis plans (phase 3) will bind it to a plan fixed before acquisition;
- the tags are trusted as given (no dataset ID yet);
- the error is statistical only and Gaussian.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from dadl import _checks
from dadl.coincidences import grid_coincidences
from dadl.errors import InsufficientCounts
from dadl.guardrails import INTERRUPTED, Policy, threshold_labels
from dadl.stationarity import interrupted
from dadl.result import Estimate, LabelledValue

FAIR_SAMPLING = "fair sampling and fair coincidence assumed (post-selected on coincidences)"
#: A fixed grid closes the event-centred window loophole but not the fair-coincidence one: a
#: local model that records every photon and shifts timing with the setting reached S = 2.8
#: (expert review). Static settings leave locality and freedom of choice open.
STATIC_SETTINGS = "static settings: locality, freedom-of-choice and memory loopholes open; not a loophole-free test"

#: Below this many coincidences per setting the Gaussian error is not trusted.
MIN_COUNTS_PER_SETTING = 100

#: Channel names a CHSH run must contain.
CHANNELS: tuple[str, ...] = ("A+", "A-", "B+", "B-")

#: Setting pairs and their signs in S = E(a,b) - E(a,b') + E(a',b) + E(a',b').
SIGNS: dict[tuple[str, str], int] = {("a", "b"): 1, ("a", "b'"): -1, ("a'", "b"): 1, ("a'", "b'"): 1}


@dataclass(frozen=True)
class Window:
    """A fixed coincidence window. `source` must be "plan" for Bell-type analyses.

    `delays_ps` maps channel name to its declared delay in integer ps.
    """

    width_ps: int
    source: str = "plan"
    delays_ps: dict[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "width_ps", _checks.integer(self.width_ps, "window width_ps", minimum=1))
        if not isinstance(self.source, str):
            raise TypeError("window source must be a string")
        delays = {str(ch): _checks.integer(d, f"delay for {ch}") for ch, d in dict(self.delays_ps).items()}
        object.__setattr__(self, "delays_ps", delays)

    def as_record(self) -> dict[str, object]:
        return {"width_ps": self.width_ps, "delays_ps": dict(sorted(self.delays_ps.items())), "source": self.source}


def correlator(tags: dict[str, np.ndarray], window: Window) -> tuple[float, float, int]:
    """E = (N++ + N-- - N+- - N-+) / N and its variance (1 - E^2) / N."""
    missing = [ch for ch in CHANNELS if ch not in tags]
    if missing:
        raise ValueError(f"run is missing channel(s) {missing}; a CHSH run needs {list(CHANNELS)}")
    d = window.delays_ps

    def n(a: str, b: str) -> int:
        return grid_coincidences(tags[a], tags[b], window.width_ps, d.get(a, 0), d.get(b, 0))

    npp, nmm, npm, nmp = n("A+", "B+"), n("A-", "B-"), n("A+", "B-"), n("A-", "B+")
    total = npp + nmm + npm + nmp
    if total < MIN_COUNTS_PER_SETTING:
        raise InsufficientCounts(
            f"{total} coincidences in a setting; need at least {MIN_COUNTS_PER_SETTING}. "
            "Check that timestamps are int64 picoseconds and that the declared delays align the channels."
        )
    e = (npp + nmm - npm - nmp) / total
    return e, (1.0 - e * e) / total, total


def setting_counts(tags: dict[str, np.ndarray], window: Window) -> dict[str, float]:
    """Singles, coincidences and the coincidence-to-singles ratio of one setting."""
    _, _, total = correlator(tags, window)
    singles_a = int(tags["A+"].size + tags["A-"].size)
    singles_b = int(tags["B+"].size + tags["B-"].size)
    return {
        "coincidences": float(total),
        "singles_a": float(singles_a),
        "singles_b": float(singles_b),
        "coincidences_per_singles": total / max(1.0, min(singles_a, singles_b)),
    }


def chsh(
    runs: dict[tuple[str, str], dict[str, np.ndarray]],
    window: Window,
    policy: Policy | None = None,
) -> LabelledValue:
    """Raw CHSH S from one run per setting pair, with a fixed window from the analysis plan.

    `runs` maps each setting pair ("a", "b"), ("a", "b'"), ("a'", "b"),
    ("a'", "b'") to a dict of channel name -> int64 ps tags for A+, A-, B+, B-.
    """
    policy = policy or Policy.default()
    policy.check_window("chsh", window.source)
    s = 0.0
    var = 0.0
    for key, sign in SIGNS.items():
        if key not in runs:
            raise ValueError(f"missing run for setting pair {key}; need {list(SIGNS)}")
        e, v, _ = correlator(runs[key], window)
        s += sign * e
        var += v
    raw = Estimate(s, math.sqrt(var))
    labels = (FAIR_SAMPLING, STATIC_SETTINGS) + policy.labels_for("chsh") + threshold_labels("chsh", raw)
    if interrupted(*(runs[k][ch] for k in SIGNS for ch in CHANNELS)):
        labels += (INTERRUPTED,)
    counts = {f"{k[0]},{k[1]}": setting_counts(runs[k], window) for k in SIGNS}
    return LabelledValue(
        analysis="chsh",
        raw=raw,
        labels=labels,
        window=window.as_record(),
        counts=counts,
        **policy.provenance(),
    )
