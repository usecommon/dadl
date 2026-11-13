"""Heralded detection probability of a signal channel from a pair source.

Method `herald-ratio` (v3). Coincidences are counted per fixed-grid slot,
present or absent, so a herald slot holds a coincidence when its partner is
detected (probability eta) or, failing that, an accidental falls in it
(probability p):

    q = P(coincidence | herald slot) = 1 - (1 - eta)(1 - p),
    eta = (q - p) / (1 - p).

p is estimated off-peak (`dadl.calibration.accidentals`): lambda accidental
coincidences expected over the S herald slots gives p = lambda / S. With C
coincidences, q has the Jeffreys Beta(C + 1/2, S - C + 1/2) posterior and
lambda the Jeffreys Gamma posterior from the shifted windows; the interval is
the posterior mixture over p, computed exactly, so it is deterministic, and draws sample p, then q.
(Methods v1-v2 used (C - lambda) / S, which tends to eta (1 - p): biased low
by a fraction p, with intervals that missed the truth from p of about 1%.)

What it measures: the probability that a herald slot holds the signal
photon's detection in the declared window, in situ ("as detected"): detector
efficiency times coupling, filters and connectors, times the share of pairs
the fixed window keeps at slot edges, times the signal channel's live
fraction given a recorded herald, diluted by dark counts on the herald
channel. That live fraction is higher than the channel's average by about
eta_a eta_b R min(dead times) (the coincidence rate times the dead time),
because a recorded herald means no herald was recorded just before it; do
not correct the result with the average dead fraction. It is not the
detector-only efficiency, which needs an external reference, and it holds at
the calibration's count rates, window and delays only (recorded in
`conditions`). The acquisition must have no analysers or beamsplitter between
source and detectors: projections enter the ratio directly.

Heralds are counted per slot, the unit C is counted in. The predictive check
is the off-peak check (`dadl.calibration.accidentals.offpeak_pvalue`).

Limits: with a CW source, independent (Poisson) multi-pair coincidences are in
the simulator and absorbed by the off-peak estimate, leaving a bias of about
-mu eta_a (1 - eta_b) / 2 for mu pairs per slot (-1.3% at mu = 0.1, eta 0.5;
source-physics review, cycle 7); pulsed and thermal pair statistics are not
modelled. Herald dark counts and afterpulses dilute the ratio (500 cps: 0.185
for 0.2; 20 kcps: 0.095; 8% herald afterpulsing: -7.6%). The source must be CW (pulsed sources need slots
locked to the laser, planned).
"""

from __future__ import annotations

import numpy as np
from scipy import optimize
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.accidentals import (
    GUARD_PS,
    N_SIDE,
    guard_for,
    offpeak_pvalue,
    shifted_coincidences,
    side_exposure,
    side_offsets,
)
from dadl.calibration.base import N_DRAWS, Calibration, with_interruption_check, draws_rng
from dadl.coincidences import grid_coincidences, occupied_slots
from dadl.data import ChannelId, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.stats import LEVEL, Interval

KIND = "heralding-efficiency"
METHOD = "herald-ratio"
METHOD_VERSION = "6"
MIN_HERALDS = 100

UNITS = {"heralding_efficiency": "1"}
DEFINITIONS = {
    "heralding_efficiency": "heralded detection probability per herald slot, as detected in the declared window "
    "(detector x coupling x window capture x live fraction given a recorded herald), diluted by the herald "
    "channel's dark counts and afterpulses"
}


#: Quadrature points over the accidental probability p.
P_POINTS = 200


#: The pair peak is searched for over +-max(PEAK_SEARCH_SLOTS windows, PEAK_SEARCH_PS) around the
#: declared slot (sized in time: 10 slots of a 200 ps window missed a 3 ns cable change; second
#: red-team review, cycle 9), at most PEAK_SEARCH_MAX_SLOTS slots each side.
PEAK_SEARCH_SLOTS = 10
PEAK_SEARCH_PS = 20_000
PEAK_SEARCH_MAX_SLOTS = 1000


def _peak_placement_pvalue(
    herald: np.ndarray, signal: np.ndarray, width_ps: int, d_h: int, d_s: int, c: int, notes: list[str]
) -> float:
    """0 when a neighbouring slot holds more coincidences than the declared one: the declared delay
    misses the pair peak (operator review, cycle 8: a 3 ns cable change read eta 0.0005, status ok).
    """
    reach = min(PEAK_SEARCH_MAX_SLOTS, max(PEAK_SEARCH_SLOTS, -(-PEAK_SEARCH_PS // width_ps)))
    shifts = np.array([k * width_ps for k in range(-reach, reach + 1) if k], dtype=np.int64)
    near = shifted_coincidences(herald, signal, width_ps, d_h, d_s, shifts)
    best = int(np.argmax(near))
    # A neighbour holding more than half the declared slot's coincidences (significantly above the
    # background) means the peak straddles or misses the slot: a half-slot misplacement read eta
    # 0.466 for 0.729.
    background = float(np.median(near))
    excess = near[best] - background
    significant = excess > 10 and excess > 5.0 * np.sqrt(near[best] + 1.0)
    signal = max(c - background, 0.0)
    if significant and abs(int(shifts[best])) == width_ps and excess <= signal:
        # Next to the declared slot and smaller than it: a misplaced peak leans one way, a window
        # narrow against the jitter spills into both neighbours alike (final statistician review,
        # cycle 10: W < 1.35 sigma was flagged although correctly placed).
        other = float(near[int(np.flatnonzero(shifts == -shifts[best])[0])]) - background
        if other >= 0.6 * excess:
            notes.append(
                "the peak spills into both neighbouring slots: the window is narrow against the timing jitter, "
                "so the heralding figure includes a small window capture"
            )
            return 1.0
    if significant and excess > 0.5 * signal:
        notes.append(
            f"the pair peak is {int(shifts[best])} ps from the declared slot ({int(near[best])} coincidences there "
            f"against {c}): the declared delays are stale; recalibrate the delay and declare it"
        )
        return 0.0
    return 1.0


def _p_nodes(side_total: int, n_windows: float, slots: int) -> np.ndarray:
    """Equal-weight quantiles of p = lambda / S, lambda ~ Gamma(side_total + 1/2, rate n_windows)."""
    q = (np.arange(P_POINTS) + 0.5) / P_POINTS
    return _st.gamma.ppf(q, side_total + 0.5, scale=1.0 / n_windows) / slots


def _interval(c: int, s: int, p_nodes: np.ndarray, level: float = LEVEL) -> tuple[float, float]:
    """Exact mixture interval for eta = (q - p) / (1 - p), clipped to [0, 1]."""
    a, b = c + 0.5, s - c + 0.5

    def cdf(x: float) -> float:
        # P(eta <= x) = E_p[ P(q <= p + x (1 - p)) ]
        with np.errstate(divide="ignore"):
            return float(np.mean(_st.beta.cdf(np.clip(p_nodes + x * (1 - p_nodes), 0, 1), a, b)))

    tail = (1 - level) / 2

    def solve(target: float) -> float:
        if cdf(0.0) >= target:
            return 0.0
        if cdf(1.0) <= target:
            return 1.0
        return optimize.brentq(lambda x: cdf(x) - target, 0.0, 1.0, xtol=1e-12)

    return solve(tail), solve(1 - tail)


def calibrate_heralding(
    herald: np.ndarray,
    signal: np.ndarray,
    channels: tuple[ChannelId, ChannelId],
    *,
    width_ps: int,
    delay_herald_ps: int = 0,
    delay_signal_ps: int = 0,
    rep_period_ps: int | None = None,
    n_side: int = N_SIDE,
    guard_ps: int = GUARD_PS,
    afterpulse_tau_ps: float | None = None,
    dead_time_ps: float | None = None,
    backend: str = "hardware",
    seed: int | None = None,
    n_draws: int = N_DRAWS,
    conditions: dict[str, float] | None = None,
) -> Calibration:
    """Heralding efficiency of `signal` given `herald`, with the declared fixed-grid window and delays."""
    if len(channels) != 2:
        raise ValueError("channels must be (herald channel, signal channel)")
    herald = _checks.timestamps(herald, "herald")
    signal = _checks.timestamps(signal, "signal")
    _checks.check_sorted(herald, "herald")
    _checks.check_sorted(signal, "signal")
    width_ps = _checks.integer(width_ps, "width_ps", minimum=1)
    s = int(occupied_slots(herald, width_ps, delay_herald_ps).size)
    shared = int(herald.size) - s
    if s < MIN_HERALDS:
        raise InsufficientCounts(f"need at least {MIN_HERALDS} herald slots, got {s}")
    c = grid_coincidences(herald, signal, width_ps, delay_herald_ps, delay_signal_ps)
    offsets = side_offsets(
        width_ps, rep_period_ps=rep_period_ps, n_side=n_side, guard_ps=guard_for(afterpulse_tau_ps, dead_time_ps, guard_ps)
    )
    side = shifted_coincidences(herald, signal, width_ps, delay_herald_ps, delay_signal_ps, offsets)
    side_total, m = int(side.sum()), side.size
    exposure = side_exposure(offsets, herald, signal)
    lam = side_total / exposure
    p_hat = min(lam / s, 0.999)

    lo, hi = _interval(c, s, _p_nodes(side_total, exposure, s))
    value = min(max((c / s - p_hat) / (1 - p_hat), 0.0), 1.0)
    rng = draws_rng(seed, herald, signal)
    p_draws = np.minimum(rng.gamma(side_total + 0.5, 1.0 / exposure, n_draws) / s, 0.999)
    q_draws = rng.beta(c + 0.5, s - c + 0.5, n_draws)
    draws = np.clip((q_draws - p_draws) / (1 - p_draws), 0.0, 1.0)

    notes = [
        f"{c} coincidences for {s} herald slots; {lam:.4g} expected accidentals from {m} shifted windows "
        f"(accidental probability {p_hat:.3g} per herald slot)",
        "heralded detection probability as detected: detector x coupling x window capture x live fraction, "
        "diluted by the herald channel's dark counts and afterpulses, "
        "at this calibration's rates, window and delays; not the detector-only efficiency",
    ]
    if dead_time_ps is not None and width_ps > 0.2 * float(dead_time_ps):
        notes.append(
            f"the {width_ps} ps window is a large fraction of the {dead_time_ps:g} ps dead time: misses caused by dead "
            "time then correlate with accidentals in the same slot, and the estimate reads high"
        )
    if dead_time_ps is not None and signal.size > 1:
        span = (int(signal[-1]) - int(signal[0])) or 1
        dead_fraction = signal.size * float(dead_time_ps) / span
        if dead_fraction > 0.05:
            notes.append(
                f"the signal channel is dead {dead_fraction:.0%} of the time: misses caused by dead time then "
                "correlate with accidentals in the same slot, and the estimate reads high"
            )
    if dead_time_ps is not None and herald.size > 1:
        span_h = (int(herald[-1]) - int(herald[0])) or 1
        conditional = (c - lam) / span_h * float(dead_time_ps)
        if conditional > 0.005:
            notes.append(
                f"the signal channel's live fraction given a recorded herald exceeds its average by about "
                f"{conditional:.1%} (coincidence rate x dead time): the result includes that, so do not correct it "
                "with the channel's average dead fraction"
            )
    if p_hat > 0.01:
        notes.append(f"accidentals are {p_hat:.1%} per herald slot; the per-slot estimator accounts for them")
    if shared:
        notes.append(f"{shared} heralds shared a slot with another herald; heralds are counted per slot")
    return Calibration(
        kind=KIND,
        method=METHOD,
        method_version=METHOD_VERSION,
        channels=tuple(channels),
        estimates={"heralding_efficiency": Interval(value, lo, hi)},
        units=dict(UNITS),
        definitions=dict(DEFINITIONS),
        # p is drawn jointly with eta (Bayesian review, cycle 3): propagating both from separate
        # artifacts would lose their correlation, which dominates eta's spread at high p.
        samples=np.column_stack([draws, p_draws]),
        sample_params=("heralding_efficiency", "accidental_probability"),
        converged=True,
        ppc_pvalue=with_interruption_check(
            min(offpeak_pvalue(side, offsets), _peak_placement_pvalue(herald, signal, width_ps, delay_herald_ps,
                                                                     delay_signal_ps, c, notes)),
            notes, herald, signal,
        ),
        inputs=((tags_sha256(herald), "raw-tags", str(channels[0])), (tags_sha256(signal), "raw-tags", str(channels[1]))),
        backend=backend,
        conditions={"width_ps": float(width_ps), **(conditions or {})},
        notes=tuple(notes),
    )


__all__ = ["KIND", "calibrate_heralding"]
