"""Bell-type fringe visibility from tags (phase 3, M2).

Two acquisitions on one pair of channels: `runs["max"]` at the analyser settings of the fringe
maximum and `runs["min"]` at the minimum (for example 0/0 and 0/90 degrees), each
{"A": tags, "B": tags}. The raw value is `bell_visibility` on the fixed-grid coincidence counts,
which requires equal exposures (within 1%).

Corrections follow the Bell rule: under the default policy none is allowed (accidental subtraction
can move a non-violating dataset past 1/sqrt(2)). With `calibrations` and a policy that allows the
corrections, the corrected value comes from the forward model (`dadl.forward`), with calibration
uncertainty, never from the old subtraction.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from dadl.analysis._twopoint import check_calibration_args, forward_check, forward_from_counts, run_counts, span_ps
from dadl.analysis.chsh import Window, check_delay_channels, check_no_duplicates, check_reference_channel
from dadl.analysis.visibility import bell_visibility
from dadl.coincidences import grid_coincidences
from dadl.data import ChannelId
from dadl.guardrails import INTERRUPTED, Policy, threshold_labels
from dadl.propagate import CalibrationSet
from dadl.result import LabelledValue
from dadl.stationarity import interrupted

FRINGE_CHANNELS: tuple[str, ...] = ("A", "B")
FRINGE_SETTINGS: tuple[str, ...] = ("max", "min")
#: Largest relative difference of the two acquisitions' spans for the raw count comparison.
MAX_EXPOSURE_MISMATCH = 0.01


def fringe_from_tags(
    runs: dict[str, dict[str, np.ndarray]],
    window: Window,
    *,
    subtract_accidentals: bool = False,
    policy: Policy | None = None,
    calibrations: CalibrationSet | None = None,
    channels: tuple[ChannelId, ChannelId] | None = None,
    seed: int = 0,
) -> tuple[LabelledValue, dict[str, float]]:
    """Fringe visibility from `runs["max"]` and `runs["min"]`. Returns the result and its counts."""
    if calibrations is not None:
        check_calibration_args(subtract_accidentals, channels)
    from dadl.calibration.accidentals import shifted_coincidences, side_offsets

    for name in FRINGE_SETTINGS:
        if name not in runs:
            raise ValueError(f"missing the {name!r} acquisition")
        missing = [ch for ch in FRINGE_CHANNELS if ch not in runs[name]]
        if missing:
            raise ValueError(f"the {name!r} acquisition is missing channel(s) {missing}")
    check_delay_channels(window, FRINGE_CHANNELS)
    check_reference_channel(window, FRINGE_CHANNELS)
    check_no_duplicates(runs, FRINGE_CHANNELS)
    d = window.delays_ps
    da, db = d.get("A", 0), d.get("B", 0)
    offsets = side_offsets(window.width_ps)

    def reduce(name: str):
        run = runs[name]
        c = grid_coincidences(run["A"], run["B"], window.width_ps, da, db)
        side = shifted_coincidences(np.sort(run["A"]), np.sort(run["B"]), window.width_ps, da, db, offsets)
        return c, side, span_ps(run, FRINGE_CHANNELS)

    c_max, side_max, t_max = reduce("max")
    c_min, side_min, t_min = reduce("min")
    if t_max <= 0 or t_min <= 0:
        raise ValueError("each acquisition needs tags spanning a nonzero time")
    # bell_visibility takes whole counts at equal exposure; the forward model uses each span.
    if abs(t_min / t_max - 1.0) > MAX_EXPOSURE_MISMATCH:
        raise ValueError(
            f"the max and min acquisitions span {t_max / 1e12:.4g} s and {t_min / 1e12:.4g} s: fringe visibility "
            "compares counts, so acquire both for the same time (within 1%)"
        )
    acc = float(side_min.sum()) / side_min.size
    counts = {
        "c_max": float(c_max),
        "c_min": float(c_min),
        "accidentals": acc,
        "accidentals_max": float(side_max.sum()) / side_max.size,
        "duration_max_s": t_max / 1e12,
        "duration_min_s": t_min / 1e12,
    }
    if subtract_accidentals:
        result = bell_visibility(
            c_max, c_min, accidentals=acc, accidentals_stderr=float(np.sqrt(side_min.sum())) / side_min.size,
            policy=policy, window=window,
        )
    else:
        result = bell_visibility(c_max, c_min, policy=policy, window=window)
    labels = result.labels
    if calibrations is not None:
        check_p, check_labels = forward_check(
            "bell_visibility", {"max": (runs["max"], (0.0,)), "min": (runs["min"], (0.0,))},
            FRINGE_CHANNELS, tuple(channels), window_ps=window.width_ps, delays_ps=dict(window.delays_ps),
            calibrations=calibrations, policy=policy,
        )
        if check_p is not None:
            counts["forward_check_p"] = check_p
        high = run_counts(runs["max"], c_max, side_max, t_max)
        low = run_counts(runs["min"], c_min, side_min, t_min)
        fwd, extra = forward_from_counts(
            "bell_visibility", high, low, tuple(channels), calibrations=calibrations, policy=policy, seed=seed,
            window_ps=window.width_ps, declared_delay_ps=float(db - da),
        )
        counts.update(extra)
        labels = tuple(dict.fromkeys(
            labels + fwd.labels + check_labels + threshold_labels("bell_visibility", result.raw, fwd.estimate)
        ))
        result = dataclasses.replace(
            result, corrected=fwd.estimate, corrections=fwd.corrections, artifact_ids=fwd.artifact_ids
        )
    if interrupted(*([np.sort(runs[k][ch]) for ch in FRINGE_CHANNELS] for k in FRINGE_SETTINGS)):
        labels += (INTERRUPTED,)
    return dataclasses.replace(result, labels=labels, counts=counts), counts


__all__ = ["FRINGE_CHANNELS", "FRINGE_SETTINGS", "fringe_from_tags"]
