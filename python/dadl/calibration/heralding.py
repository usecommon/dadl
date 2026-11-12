"""Heralding efficiency of a signal channel from a pair source.

Method `herald-ratio`: with S herald detections, C fixed-grid coincidences
in the declared window and lambda the expected accidentals in that window
(estimated off-peak, as in `dadl.calibration.accidentals`),

    eta = (C - lambda) / S,

the in-situ system efficiency of the signal arm "as detected" (detector times
coupling, filters and connectors; NIST IR 8486r1), including the share of
pairs the fixed window cuts at slot edges. It is not the detector-only
efficiency, which needs an external reference.

Posterior: the accidental part of C, A, has the Gamma-Poisson
(negative binomial) predictive from the off-peak counts; given A, the pair
coincidences C - A out of S heralds have the Jeffreys Beta posterior. The
interval is computed from the exact mixture CDF, so it is deterministic;
draws sample A, then eta.

The predictive check is the off-peak dispersion test.

Phase 2 limits: dark counts on the herald channel count as heralds (they
lower eta slightly, as in every heralding-ratio measurement); multi-pair
emission is not modelled.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.accidentals import GUARD_PS, N_SIDE, dispersion_pvalue, shifted_coincidences, side_offsets
from dadl.calibration.base import N_DRAWS, Calibration
from dadl.coincidences import grid_coincidences
from dadl.data import ChannelId, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.stats import LEVEL, Interval

KIND = "heralding-efficiency"
METHOD = "herald-ratio"
METHOD_VERSION = "1"
MIN_HERALDS = 100

UNITS = {"heralding_efficiency": "probability per herald"}


def _accidental_pmf(side_total: int, n_windows: int) -> tuple[np.ndarray, np.ndarray]:
    """Support and weights of the on-peak accidental count (Gamma-Poisson predictive)."""
    nb = _st.nbinom(side_total + 0.5, n_windows / (n_windows + 1.0))
    lo, hi = int(nb.ppf(1e-9)), int(nb.ppf(1 - 1e-9))
    k = np.arange(lo, hi + 1)
    w = nb.pmf(k)
    return k, w / w.sum()


def _interval(c: int, s: int, k: np.ndarray, w: np.ndarray, level: float = LEVEL) -> tuple[float, float]:
    pairs = np.clip(c - k, 0, s)
    alpha, beta = pairs + 0.5, s - pairs + 0.5

    def cdf(x: float) -> float:
        return float(np.sum(w * _st.beta.cdf(x, alpha, beta)))

    q = (1 - level) / 2
    lo = optimize.brentq(lambda x: cdf(x) - q, 0.0, 1.0, xtol=1e-12)
    hi = optimize.brentq(lambda x: cdf(x) - (1 - q), 0.0, 1.0, xtol=1e-12)
    return lo, hi


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
    s = int(herald.size)
    if s < MIN_HERALDS:
        raise InsufficientCounts(f"need at least {MIN_HERALDS} herald detections, got {s}")
    c = grid_coincidences(herald, signal, width_ps, delay_herald_ps, delay_signal_ps)
    offsets = side_offsets(width_ps, rep_period_ps=rep_period_ps, n_side=n_side, guard_ps=guard_ps)
    side = shifted_coincidences(herald, signal, width_ps, delay_herald_ps, delay_signal_ps, offsets)
    side_total, m = int(side.sum()), side.size
    lam = side_total / m

    k, w = _accidental_pmf(side_total, m)
    lo, hi = _interval(c, s, k, w)
    value = min(max((c - lam) / s, 0.0), 1.0)
    rng = np.random.default_rng(seed)
    acc = rng.choice(k, size=n_draws, p=w)
    pairs = np.clip(c - acc, 0, s)
    draws = rng.beta(pairs + 0.5, s - pairs + 0.5)

    notes = [
        f"{c} coincidences for {s} heralds; {lam:.4g} expected accidentals from {m} shifted windows",
        "system efficiency as detected (detector x coupling), including window losses at slot edges",
    ]
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
        ppc_pvalue=dispersion_pvalue(side),
        inputs=((tags_sha256(herald), "raw-tags"), (tags_sha256(signal), "raw-tags")),
        backend=backend,
        conditions={"width_ps": float(width_ps), **(conditions or {})},
        notes=tuple(notes),
    )


__all__ = ["KIND", "calibrate_heralding"]
