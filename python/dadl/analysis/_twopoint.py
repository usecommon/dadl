"""Shared tag-level step of the two-point visibilities (HOM and fringe): reduce two acquisitions
to `dadl.forward.RunCounts` and run the forward model with propagated calibrations."""

from __future__ import annotations

import numpy as np

from dadl.data import ChannelId
from dadl.guardrails import Policy
from dadl.propagate import CalibrationSet


def span_ps(run: dict[str, np.ndarray], channels: tuple[str, str]) -> float:
    """An acquisition's span: first to last tag over both channels."""
    arrays = [np.asarray(run[ch]) for ch in channels if np.asarray(run[ch]).size]
    if not arrays:
        return 0.0
    return float(max(int(a.max()) for a in arrays) - min(int(a.min()) for a in arrays))


def run_counts(
    run: dict[str, np.ndarray], coincidences: int, side: np.ndarray, span: float, channels=("A", "B"), offsets=(0.0,)
):
    from dadl.forward import RunCounts

    a, b = channels
    return RunCounts(
        coincidences=int(coincidences),
        offpeak_sum=int(side.sum()),
        offpeak_shifts=int(side.size),
        duration_s=span / 1e12,
        singles_a=int(np.asarray(run[a]).size),
        singles_b=int(np.asarray(run[b]).size),
        pair_offsets_ps=tuple(float(o) for o in offsets),
    )


def forward_from_counts(
    analysis: str,
    high,
    low,
    channels: tuple[ChannelId, ChannelId],
    *,
    calibrations: CalibrationSet,
    policy: Policy | None,
    seed: int,
    window_ps: int,
    declared_delay_ps: float,
) -> tuple[object, dict[str, float]]:
    """The forward-model correction and the counts it adds (the D3 range when a dead time is a bound,
    labelled when it is wider than the corrected value's standard error)."""
    import dataclasses

    from dadl.forward import DEAD_TIME_RANGE_WIDE, bound_range, corrected_visibility

    grid = {"window_ps": window_ps, "declared_delay_ps": declared_delay_ps}
    fwd = corrected_visibility(analysis, high, low, tuple(channels), calibrations=calibrations, policy=policy, seed=seed, **grid)
    extra: dict[str, float] = {}
    rng = bound_range(analysis, high, low, tuple(channels), calibrations=calibrations, policy=policy, **grid)
    if rng is not None:
        extra["corrected_at_zero_dead_time"] = rng["at_zero"]
        extra["corrected_at_bound_dead_time"] = rng["at_bound"]
        if abs(rng["at_bound"] - rng["at_zero"]) > fwd.estimate.stderr:
            fwd = dataclasses.replace(fwd, labels=fwd.labels + (DEAD_TIME_RANGE_WIDE,))
    return fwd, extra


def use_checks(
    analysis: str,
    runs: list[dict[str, np.ndarray]],
    names: tuple[str, str],
    channels: tuple[ChannelId, ChannelId],
    *,
    calibrations: CalibrationSet,
    policy: Policy | None,
) -> tuple[tuple[str, ...], dict[str, object]]:
    """The artifacts' detectors against the analysis channels (`propagate.check_detectors`), the
    analysed gaps against the dark dead times (`forward.check_dead_time_gaps`) and the dataset
    context the calibrations were checked against: (labels, counts)."""
    from dadl.forward import check_dead_time_gaps
    from dadl.propagate import check_detectors

    labels = check_detectors(calibrations, channels, analysis=analysis, policy=policy)
    gap_labels, counts = check_dead_time_gaps(calibrations, dict(zip(names, channels)), runs, analysis=analysis,
                                              policy=policy)
    labels = tuple(dict.fromkeys(labels + gap_labels))
    record = calibrations.dataset_record()
    if record is not None:
        counts["dataset"] = record
    return labels, counts


def range_ends(extra: dict[str, float]) -> tuple[float, ...]:
    """The corrected value at both ends of the D3 range, if any (for `threshold_labels`)."""
    return tuple(extra[k] for k in ("corrected_at_zero_dead_time", "corrected_at_bound_dead_time") if k in extra)


def check_calibration_args(subtract_accidentals: bool | None, channels) -> None:
    from dadl.errors import GuardrailRefusal

    if subtract_accidentals:
        raise GuardrailRefusal(
            "subtract_accidentals and calibrations both correct the accidentals: use the forward model "
            "(calibrations) alone, so no correction is counted twice"
        )
    if channels is None or len(channels) != 2:
        raise ValueError("with calibrations, pass channels=(channel of A, channel of B)")


def forward_check(
    analysis: str,
    runs: dict[str, tuple[dict[str, np.ndarray], tuple[float, ...]]],
    names: tuple[str, str],
    channels: tuple[ChannelId, ChannelId],
    *,
    window_ps: int,
    delays_ps: dict[str, int],
    calibrations: CalibrationSet,
    policy: Policy | None,
) -> tuple[float | None, tuple[str, ...]]:
    """The forward model's predictive check on each acquisition's histogram of the channel pair
    (`dadl.forward_check`): the adjusted p-value and labels; refuses a misfit the policy does not allow.

    Runs only when a delay-jitter calibration gives the peak shape the model uses (window capture);
    without one the model makes no claim about the shape (p is None)."""
    from dadl.forward_check import PairData, delay_start, enforce, run_check
    from dadl.guardrails import Policy as _Policy

    a, b = names
    declared = float(delays_ps.get(b, 0) - delays_ps.get(a, 0))
    d_cal, jitter = delay_start(calibrations, *channels)
    if jitter is None:
        return None, ()
    pairs = [
        PairData(name, np.sort(run[a]), np.sort(run[b]), declared, window_ps, d_cal - declared, jitter, offsets)
        for name, (run, offsets) in runs.items()
    ]
    result = run_check(pairs)
    return result.pvalue, enforce(result, analysis, policy or _Policy.default())
