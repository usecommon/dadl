"""Shared tag-level step of the two-point visibilities (HOM and fringe): reduce two acquisitions
to `dadl.forward.RunCounts` and run the forward model with propagated calibrations."""

from __future__ import annotations

import numpy as np

from dadl.coincidences import as_sorted as _as_sorted
from dadl.data import ChannelId, run_t3_period
from dadl.guardrails import Policy
from dadl.propagate import CalibrationSet


#: Off-peak flatness check (council round 8, Ph3): its share of the 1e-3 per-analysis refusal
#: budget, Bonferroni over the acquisitions (or acquisitions x outcome pairs) checked.
OFFPEAK_ALPHA = 2e-4
#: Sequential Monte Carlo of the exact dispersion test (council round 9, Pf1): stop at this many
#: exceedances (p = h / draws). A fixed 20,000 draws gave a smallest doubled p of 1e-4, which could
#: never reach the threshold with 16 sides (CHSH), only with 2 (HOM and fringe).
OFFPEAK_H = 20
#: n_max = OFFPEAK_N_MAX_FACTOR x 2 x sides / alpha draws: the smallest p is a fifth of a side's
#: dispersion threshold.
OFFPEAK_N_MAX_FACTOR = 4
OFFPEAK_REFUSAL = (
    "the off-peak windows are not flat (pulsed light or a periodic background): the accidentals cannot be read "
    "off-peak"
)
OFFPEAK_ACCEPTED = (
    "off-peak windows not flat, accepted: the accidentals read off-peak may be wrong (pulsed light or a periodic "
    "background; policy modified)"
)


def sequential_dispersion_pvalue(
    counts: np.ndarray,
    n_sides: int = 1,
    alpha: float = OFFPEAK_ALPHA,
    h: int = OFFPEAK_H,
    n_max_factor: int = OFFPEAK_N_MAX_FACTOR,
    stop_at: int | None = None,
) -> tuple[float, int]:
    """Exact dispersion p-value of shifted counts by sequential Monte Carlo (Besag and Clifford
    1991): (p, draws used).

    Multinomial null tables (the total spread evenly over the windows) are drawn in growing chunks.
    At the h-th draw whose statistic reaches the observed one, p = h / draws so far. Otherwise the
    draws stop at n_max = n_max_factor x 2 x n_sides / alpha, with p = (g + 1) / (n_max + 1) for g
    exceedances. Either way P(p <= u) <= u. The generator is seeded from the counts, so the same data
    give the same p, and different sides use unrelated streams (a fixed seed made the 5% level run at
    6%: council round 9). At a mean of `EXACT_BELOW` per window or more the chi-square law is used
    (0 draws). `stop_at` ends the draws at that many exceedances (p = stop_at / draws, an early
    stop used once the analysis is already refused)."""
    from scipy import stats as st

    from dadl.calibration.accidentals import EXACT_BELOW

    counts = np.asarray(counts)
    m = float(np.mean(counts)) if counts.size else 0.0
    if counts.size < 2 or m <= 0:
        return 1.0, 0
    stat = float(np.sum((counts - m) ** 2) / m)
    if m >= EXACT_BELOW:
        return float(st.chi2.sf(stat, counts.size - 1)), 0
    n_max = int(np.ceil(n_max_factor * 2 * n_sides / alpha))
    stop = h if stop_at is None else min(h, int(stop_at))
    total, k = int(counts.sum()), counts.size
    rng = np.random.default_rng(np.random.SeedSequence(counts.astype(np.uint32)))
    exceeded = drawn = 0
    chunk = 64
    while drawn < n_max:
        size = min(chunk, n_max - drawn)
        sims = rng.multinomial(total, np.full(k, 1.0 / k), size=size)
        hits = np.cumsum(np.sum((sims - m) ** 2, axis=1) / m >= stat - 1e-9)
        if exceeded + hits[-1] >= stop:
            at = int(np.searchsorted(hits, stop - exceeded)) + 1  # draws in this chunk up to the stop
            return stop / (drawn + at), drawn + at
        exceeded += int(hits[-1])
        drawn += size
        chunk = min(chunk * 2, 16384)
    return (exceeded + 1) / (n_max + 1), drawn


def _offpeak(side: np.ndarray, offsets: np.ndarray, n_sides: int, stop_at: int | None = None) -> tuple[float, int]:
    from dadl.calibration.accidentals import trend_pvalue

    p_disp, draws = sequential_dispersion_pvalue(side, n_sides, stop_at=stop_at)
    return min(1.0, 2.0 * min(p_disp, trend_pvalue(side, offsets))), draws


def offpeak_p(side: np.ndarray, offsets: np.ndarray, n_sides: int = 1) -> float:
    """The off-peak flatness p-value of one acquisition's shifted counts among `n_sides` checked:
    min(dispersion, trend), doubled, with the sequential dispersion test
    (`sequential_dispersion_pvalue`)."""
    return _offpeak(side, offsets, n_sides)[0]


def _refusing_exceedances(n_sides: int) -> int:
    """The fewest dispersion exceedances whose p at n_max, (g + 1) / (n_max + 1), no longer flags a
    side: once a side has this many, its dispersion half cannot flag it (4 for 2 or 16 sides)."""
    n_max = int(np.ceil(OFFPEAK_N_MAX_FACTOR * 2 * n_sides / OFFPEAK_ALPHA))
    g = 0
    while 2.0 * (g + 1) / (n_max + 1) < OFFPEAK_ALPHA / n_sides:
        g += 1
    return g


def offpeak_check(
    analysis: str, sides: dict[str, np.ndarray], offsets: np.ndarray, policy: Policy | None
) -> tuple[float, tuple[str, ...], int]:
    """Check that each acquisition's off-peak windows are flat, as the accidentals assume (council
    round 8, Ph3: pulsed light on T2 or plain tags read V 10-40 sd off, unlabelled). Returns the
    smallest p, Bonferroni-adjusted over `sides`, the labels and the Monte Carlo draws used. Flagged
    when the smallest p is below `OFFPEAK_ALPHA` / len(sides): refused (`Misspecified`) unless the
    policy allows `misfit_forward_model`, then labelled `OFFPEAK_ACCEPTED`. Call it only where an
    accidental-subtracted value is computed.

    Once one side is flagged, the remaining sides are drawn only until they cannot be flagged
    (a flagged side runs to n_max, about 3 s with 16 sides).

    Power (council round 9, L4): with about one off-peak coincidence per window, as in a typical
    CHSH run, the check reliably finds pulsed accidentals only when they are about half or more of
    one acquisition's off-peak counts. In council round 9's test it refused 30 of 30 runs at one
    half, 8-17 of 30 at 30%, and 0-2 of 30 at 20%. A smaller pulsed share passes unlabelled. The
    off-peak windows then miss part of the same-pulse accidentals, so V and |E| read low
    (conservative for the threshold). Longer acquisitions raise the counts per window and the power."""
    from dadl.errors import Misspecified
    from dadl.forward_check import OVERRIDE
    from dadl.guardrails import Policy as _Policy, override_hint

    if not sides:
        return 1.0, (), 0
    n = len(sides)
    ps: dict[str, float] = {}
    draws = 0
    stop_at = None
    for name, side in sides.items():
        ps[name], used = _offpeak(side, offsets, n, stop_at)
        draws += used
        if stop_at is None and ps[name] < OFFPEAK_ALPHA / n:
            stop_at = _refusing_exceedances(n)
    worst = min(ps, key=ps.get)
    adjusted = min(1.0, ps[worst] * n)
    if ps[worst] >= OFFPEAK_ALPHA / n:
        return adjusted, (), draws
    policy = policy or _Policy.default()
    if OVERRIDE not in policy.rule_for(analysis).allowed_corrections:
        raise Misspecified(
            f"{OFFPEAK_REFUSAL} ({worst}: p = {adjusted:.3g} after Bonferroni over {n}). "
            f"{override_hint(analysis, OVERRIDE)}"
        )
    return adjusted, (OFFPEAK_ACCEPTED,), draws


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
) -> tuple[object, dict[str, object]]:
    """The forward-model correction and the counts it adds (the Monte Carlo seed and draws; the D3
    range when a dead time is a bound, labelled when it is wider than the corrected value's standard
    error)."""
    import dataclasses

    from dadl.forward import DEAD_TIME_RANGE_WIDE, bound_range, corrected_visibility

    grid = {"window_ps": window_ps, "declared_delay_ps": declared_delay_ps}
    fwd = corrected_visibility(analysis, high, low, tuple(channels), calibrations=calibrations, policy=policy, seed=seed, **grid)
    extra: dict[str, object] = dict(fwd.counts)
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
    lattice_ps: int | None = None,
) -> tuple[float | None, tuple[str, ...]]:
    """The forward model's predictive check on each acquisition's histogram of the channel pair
    (`dadl.forward_check`): the adjusted p-value and labels; refuses a misfit the policy does not allow.

    Runs only when a delay-jitter calibration gives the peak shape the model uses (window capture);
    without one the model makes no claim about the shape (p is None, labelled `NO_PREDICTIVE_CHECK`;
    council round 10, Ux7). `lattice_ps` is the runs' lattice (`dadl.data.runs_lattice`); at 0 (no
    integer lattice) the check is skipped (p None)."""
    from dadl.forward_check import NO_PREDICTIVE_CHECK, PairData, delay_start, enforce, run_check
    from dadl.guardrails import Policy as _Policy

    a, b = names
    declared = float(delays_ps.get(b, 0) - delays_ps.get(a, 0))
    d_cal, jitter = delay_start(calibrations, *channels)
    if lattice_ps == 0:
        return None, ()
    if jitter is None:
        return None, (NO_PREDICTIVE_CHECK,)
    pairs = [
        PairData(name, _as_sorted(run[a]), _as_sorted(run[b]), declared, window_ps, d_cal - declared, jitter, offsets,
                 lattice_ps, run_t3_period(run))
        for name, (run, offsets) in runs.items()
    ]
    result = run_check(pairs)
    return result.pvalue, enforce(result, analysis, policy or _Policy.default())
