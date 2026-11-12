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
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np

from dadl.analysis.chsh import Window
from dadl.analysis.visibility import hom_visibility
from dadl.coincidences import grid_coincidences
from dadl.guardrails import Policy
from dadl.result import LabelledValue

#: Channel names of a HOM setup: the two beamsplitter outputs.
HOM_CHANNELS: tuple[str, ...] = ("A", "B")
#: Acquisitions of a HOM setup, in tag-file setting order.
HOM_SETTINGS: tuple[str, ...] = ("dip", "reference")


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
    d = window.delays_ps
    da, db = d.get("A", 0), d.get("B", 0)

    def count(run: dict[str, np.ndarray]) -> int:
        return grid_coincidences(run["A"], run["B"], window.width_ps, da, db)

    c_min, c_max = count(runs["dip"]), count(runs["reference"])
    offsets = side_offsets(window.width_ps)
    side = shifted_coincidences(np.sort(runs["dip"]["A"]), np.sort(runs["dip"]["B"]), window.width_ps, da, db, offsets)
    acc = float(side.sum()) / side.size
    acc_err = math.sqrt(float(side.sum())) / side.size
    counts = {"c_min": float(c_min), "c_max": float(c_max), "accidentals": acc, "accidentals_stderr": acc_err}
    if subtract_accidentals and acc <= c_min:
        result = hom_visibility(c_max, c_min, accidentals=acc, accidentals_stderr=acc_err, policy=policy)
    else:
        result = hom_visibility(c_max, c_min, policy=policy)
    return dataclasses.replace(result, window=window.as_record()), counts
