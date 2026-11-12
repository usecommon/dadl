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
the exact mixture over p, so it is deterministic, and draws sample p, then q.
(Methods v1-v2 used (C - lambda) / S, which tends to eta (1 - p): biased low
by a fraction p, with intervals that missed the truth from p of about 1%.)

What it measures: the probability that a herald slot holds the signal
photon's detection in the declared window, in situ ("as detected"): detector
efficiency times coupling, filters and connectors, times the share of pairs
the fixed window keeps at slot edges, times the signal channel's live
fraction, diluted by dark counts on the herald channel. It is not the
detector-only efficiency, which needs an external reference, and it holds at
the calibration's count rates, window and delays only (recorded in
`conditions`). The acquisition must have no analysers or beamsplitter between
source and detectors: projections enter the ratio directly.

Heralds are counted per slot, the unit C is counted in. The predictive check
is the off-peak check (`dadl.calibration.accidentals.offpeak_pvalue`).

Limits: multi-pair emission is not modelled (it lowers the ratio by about the
mean pair number per slot); the source must be CW (pulsed sources need slots
locked to the laser, planned).
"""

from __future__ import annotations

import numpy as np
from scipy import optimize
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.accidentals import GUARD_PS, N_SIDE, guard_for, offpeak_pvalue, shifted_coincidences, side_offsets
from dadl.calibration.base import N_DRAWS, Calibration
from dadl.coincidences import grid_coincidences, occupied_slots
from dadl.data import ChannelId, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.stats import LEVEL, Interval

KIND = "heralding-efficiency"
METHOD = "herald-ratio"
METHOD_VERSION = "3"
MIN_HERALDS = 100

UNITS = {"heralding_efficiency": "1 (per herald slot)"}


#: Quadrature points over the accidental probability p.
P_POINTS = 200


def _p_nodes(side_total: int, n_windows: int, slots: int) -> np.ndarray:
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
    lam = side_total / m
    p_hat = min(lam / s, 0.999)

    lo, hi = _interval(c, s, _p_nodes(side_total, m, s))
    value = min(max((c / s - p_hat) / (1 - p_hat), 0.0), 1.0)
    rng = np.random.default_rng(seed)
    p_draws = np.minimum(rng.gamma(side_total + 0.5, 1.0 / m, n_draws) / s, 0.999)
    q_draws = rng.beta(c + 0.5, s - c + 0.5, n_draws)
    draws = np.clip((q_draws - p_draws) / (1 - p_draws), 0.0, 1.0)

    notes = [
        f"{c} coincidences for {s} herald slots; {lam:.4g} expected accidentals from {m} shifted windows "
        f"(accidental probability {p_hat:.3g} per herald slot)",
        "heralded detection probability as detected: detector x coupling x window capture x live fraction, "
        "at this calibration's rates, window and delays; not the detector-only efficiency",
    ]
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
        samples=draws[:, None],
        sample_params=("heralding_efficiency",),
        converged=True,
        ppc_pvalue=offpeak_pvalue(side, offsets),
        inputs=((tags_sha256(herald), "raw-tags"), (tags_sha256(signal), "raw-tags")),
        backend=backend,
        conditions={"width_ps": float(width_ps), **(conditions or {})},
        notes=tuple(notes),
    )


__all__ = ["KIND", "calibrate_heralding"]
