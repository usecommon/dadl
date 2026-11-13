"""HOM dip visibility from time tags: two acquisitions, at the dip and far from it.

The two acquisitions share one fixed-grid window (declared width and delays).
C_min is the coincidence count at the dip (path delay 0), C_max at the
reference delay, many coherence times away. The raw visibility is
(C_max - C_min) / C_max. Accidentals in the window are estimated off-peak
from the dip acquisition (`dadl.calibration.accidentals`) and subtracted only
in the labelled secondary value, which the default policy allows for HOM.

The dip acquisition's accidentals are used for both counts: singles fall at
the dip (both photons on one detector count once), so these are the
accidentals that matter for C_min; the error this makes in C_max is of order
A / C_max.

Exposure: each acquisition's duration is its span of tags, and V compares
rates (`hom_visibility(..., exposure_ratio=T_dip / T_ref)`). The recorded
singles cannot normalise drift: with detectors that do not resolve photon
number, singles fall at the dip by design ((eta - eta^2/2) / (eta - eta^2/4)),
so singles normalisation reads V low (0.857 for 0.919).

Limits: two-fold HOM only (not four-fold heralded HOM between separate
sources); a 50:50 beamsplitter (R:T imbalance caps V at 2RT / (R^2 + T^2));
the reference delay must be small against the window minus the jitter, since
the fixed grid loses part of the split peak at +-delay (V reads low: 0.86 at
a 100 ps reference delay in a 200 ps window); CW sources.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np

from dadl.analysis.chsh import Window, check_delay_channels, check_no_duplicates, check_reference_channel
from dadl.analysis.visibility import hom_visibility
from dadl.coincidences import grid_coincidences
from dadl.guardrails import INTERRUPTED, Policy
from dadl.stationarity import interrupted
from dadl.result import LabelledValue

#: The visibility comes from two acquisitions, not a fitted scan: an off-centre dip (0.3 ps gave
#: 0.882 for 0.919) or a pump drift between them (+-5% gave +-3 sigma) moves it, unlabelled
#: before (expert review).
TWO_POINT = (
    "two-point visibility at the declared delays: a lower bound if the dip is off-centre; "
    "not normalised for rate drift between the acquisitions"
)

#: Channel names of a HOM setup: the two beamsplitter outputs.
HOM_CHANNELS: tuple[str, ...] = ("A", "B")
#: Acquisitions of a HOM setup, in tag-file setting order.
HOM_SETTINGS: tuple[str, ...] = ("dip", "reference")


def _span_ps(run: dict[str, np.ndarray]) -> float:
    """An acquisition's span: first to last tag over both channels."""
    arrays = [np.asarray(run[ch]) for ch in HOM_CHANNELS if np.asarray(run[ch]).size]
    if not arrays:
        return 0.0
    return float(max(int(a.max()) for a in arrays) - min(int(a.min()) for a in arrays))


def hom_from_tags(
    runs: dict[str, dict[str, np.ndarray]],
    window: Window,
    *,
    subtract_accidentals: bool = True,
    policy: Policy | None = None,
) -> tuple[LabelledValue, dict[str, float]]:
    """HOM visibility from `runs["dip"]` and `runs["reference"]` (each {"A": tags, "B": tags}).

    Returns the labelled result and the counts it used.
    """
    from dadl.calibration.accidentals import shifted_coincidences, side_offsets

    for name in HOM_SETTINGS:
        if name not in runs:
            raise ValueError(f"missing the {name!r} acquisition")
        missing = [ch for ch in HOM_CHANNELS if ch not in runs[name]]
        if missing:
            raise ValueError(f"the {name!r} acquisition is missing channel(s) {missing}")
    check_delay_channels(window, HOM_CHANNELS)
    check_reference_channel(window, HOM_CHANNELS)
    check_no_duplicates(runs, HOM_CHANNELS)
    d = window.delays_ps
    da, db = d.get("A", 0), d.get("B", 0)

    def count(run: dict[str, np.ndarray]) -> int:
        return grid_coincidences(run["A"], run["B"], window.width_ps, da, db)

    c_min, c_max = count(runs["dip"]), count(runs["reference"])
    offsets = side_offsets(window.width_ps)
    side = shifted_coincidences(np.sort(runs["dip"]["A"]), np.sort(runs["dip"]["B"]), window.width_ps, da, db, offsets)
    acc = float(side.sum()) / side.size
    acc_err = math.sqrt(float(side.sum())) / side.size
    counts = {
        "c_min": float(c_min),
        "c_max": float(c_max),
        "accidentals": acc,
        "accidentals_stderr": acc_err,
        "singles_dip": float(runs["dip"]["A"].size + runs["dip"]["B"].size),
        "singles_reference": float(runs["reference"]["A"].size + runs["reference"]["B"].size),
    }
    # Exposure from each acquisition's own span: a dip run half as long as the reference read
    # V = 0.960 for 0.92 (HOM review, cycle 3).
    t_dip, t_ref = _span_ps(runs["dip"]), _span_ps(runs["reference"])
    ratio = t_dip / t_ref if t_dip > 0 and t_ref > 0 else 1.0
    counts["duration_dip_s"] = t_dip / 1e12
    counts["duration_reference_s"] = t_ref / 1e12
    if subtract_accidentals:
        # hom_visibility decides whether the subtraction can be done, and labels it if not.
        result = hom_visibility(
            c_max, c_min, accidentals=acc, accidentals_stderr=acc_err, policy=policy, exposure_ratio=ratio
        )
    else:
        result = hom_visibility(c_max, c_min, policy=policy, exposure_ratio=ratio)
    counts["subtracted"] = float(result.corrected is not None)
    labels = result.labels + (TWO_POINT,)
    if interrupted(*([np.sort(runs[k][ch]) for ch in HOM_CHANNELS] for k in HOM_SETTINGS)):
        # SNSPD specialist review: a latch in 20% of the dip run raised V by 13 sigma, unlabelled.
        labels += (INTERRUPTED,)
    return dataclasses.replace(result, window=window.as_record(), labels=labels, counts=counts), counts
