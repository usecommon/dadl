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
    """The forward-model correction and the counts it adds (the D3 range when a dead time is a bound)."""
    from dadl.forward import bound_range, corrected_visibility

    grid = {"window_ps": window_ps, "declared_delay_ps": declared_delay_ps}
    fwd = corrected_visibility(analysis, high, low, tuple(channels), calibrations=calibrations, policy=policy, seed=seed, **grid)
    extra: dict[str, float] = {}
    rng = bound_range(analysis, high, low, tuple(channels), calibrations=calibrations, policy=policy, **grid)
    if rng is not None:
        extra["corrected_at_zero_dead_time"] = rng["at_zero"]
        extra["corrected_at_bound_dead_time"] = rng["at_bound"]
    return fwd, extra


def check_calibration_args(subtract_accidentals: bool | None, channels) -> None:
    from dadl.errors import GuardrailRefusal

    if subtract_accidentals:
        raise GuardrailRefusal(
            "subtract_accidentals and calibrations both correct the accidentals: use the forward model "
            "(calibrations) alone, so no correction is counted twice"
        )
    if channels is None or len(channels) != 2:
        raise ValueError("with calibrations, pass channels=(channel of A, channel of B)")
