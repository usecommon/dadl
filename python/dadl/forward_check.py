"""Predictive check of the forward model on the analysed acquisitions (phase 3, M3).

The forward model (`dadl.forward`) claims that, within a coincidence window, the coincidences of
a channel pair are a peak with the calibrated IRF shape (a Gaussian, or the equal mixture at the
declared pair offsets, such as +-d for a HOM path delay) on an accidental background that is
flat apart from a step at the peak (dead time removes background on one side). Anything else
that lands in the window, crosstalk at its own delay above all, is counted as signal and biases
the corrected value with an interval that does not know it.

The check fits that shape to each analysed histogram of delays b - a, over |b - a - declared
delay| < W (every delay the window grid can count), by Poisson maximum likelihood with the delay,
jitter, peak area and the two background levels free, and takes a grouped Pearson p-value on
(groups - 1) degrees of freedom (`_grouped_pearson`, which says why). The analysis is flagged when
the smallest p falls below `MISFIT_ALPHA` / (number of histograms) (Bonferroni, so the false-alarm
rate per analysis is at most about `MISFIT_ALPHA` however many pairs it holds).

The check runs only where the model uses the peak shape: when a delay-jitter calibration of the
pair gives the window capture. A flagged analysis is refused (`Misspecified`) unless the policy allows `misfit_forward_model`;
then the corrected value carries `FORWARD_MISFIT_ACCEPTED`.

What it does not check: the calibrated delay and jitter against the analysed run's own (they are
refitted here, so a drift between calibration and run is not seen), and anything outside the
window.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy import optimize
from scipy import stats as _st

from dadl.calibration.delay import _g
from dadl.data import tag_tick
from dadl.errors import Misspecified
from dadl.guardrails import Policy
from dadl.histogram import cross_correlation

#: Family-wise false-alarm rate of the check per analysis.
MISFIT_ALPHA = 0.001
OVERRIDE = "misfit_forward_model"
FORWARD_MISFIT_ACCEPTED = "forward-model misfit accepted: the predictive check failed (policy modified)"
#: Fewest coincidences in the fit region for a histogram to be checked.
MIN_COUNTS = 50
#: Adjacent bins are merged until each expects at least this many counts.
MIN_EXPECTED = 5.0


def _grouped_pearson(n: np.ndarray, mu: np.ndarray) -> float:
    """Pearson chi-square p-value over bins merged until each expects >= MIN_EXPECTED, on
    (groups - 1) degrees of freedom.

    The parameters are fitted on the fine bins, not the merged ones, so the statistic lies between
    chi2(groups - params) and chi2(groups - 1) (Chernoff and Lehmann, 1954). With groups - params,
    clean simulated runs fell below p = 0.001 at about three times the nominal rate; groups - 1
    bounds it (about half the nominal rate at p = 0.01 and 0.05, over 1,500 HOM and 4,800 CHSH
    histograms).
    """
    obs, exp = [], []
    acc_n = acc_mu = 0.0
    for ni, mi in zip(n, mu):
        acc_n += ni
        acc_mu += mi
        if acc_mu >= MIN_EXPECTED:
            obs.append(acc_n)
            exp.append(acc_mu)
            acc_n = acc_mu = 0.0
    if acc_mu > 0 and exp:
        obs[-1] += acc_n
        exp[-1] += acc_mu
    dof = len(exp) - 1
    if dof < 1:
        return 1.0
    o, e = np.array(obs), np.array(exp)
    return float(_st.chi2.sf(float(np.sum((o - e) ** 2 / e)), dof))


@dataclass(frozen=True)
class PairData:
    """One analysed coincidence histogram's inputs."""

    label: str
    a: np.ndarray
    b: np.ndarray
    #: The window's declared delay of b relative to a (ps); the histogram is centred on it.
    declared_delay_ps: float
    window_ps: int
    #: Starting values: calibrated delay relative to the declared one, and jitter (ps).
    delay_ps: float = 0.0
    jitter_ps: float | None = None
    offsets_ps: tuple[float, ...] = (0.0,)


@dataclass(frozen=True)
class CheckResult:
    pvalues: dict[str, float]
    alpha: float

    @property
    def pvalue(self) -> float:
        """The smallest per-histogram p-value, Bonferroni-adjusted (capped at 1)."""
        if not self.pvalues:
            return 1.0
        return min(1.0, min(self.pvalues.values()) * len(self.pvalues))

    @property
    def flagged(self) -> bool:
        return self.pvalue < self.alpha

    def worst(self) -> str:
        return min(self.pvalues, key=self.pvalues.get) if self.pvalues else ""


class _ShapeModel:
    """Poisson NLL of the forward model's peak shape: mixture of Gaussians at d + offsets, on a
    background with one level each side of d. theta = (log bL, log A, d, log sigma, log bR)."""

    def __init__(self, counts: np.ndarray, edges: np.ndarray, tick: float, offsets: Sequence[float]) -> None:
        self.n = counts.astype(float)
        self.lo, self.hi = edges[:-1], edges[1:]
        self.w = self.hi - self.lo
        self.tick = float(tick)
        self.offsets = tuple(float(o) for o in offsets)

    def mu(self, theta: np.ndarray) -> np.ndarray:
        log_bl, log_a, d, log_s, log_br = (float(v) for v in theta)
        s = math.exp(min(max(log_s, -50.0), 700.0))
        t = self.tick / s
        g = np.zeros_like(self.n)
        for o in self.offsets:
            z_hi, z_lo = (self.hi - d - o) / s, (self.lo - d - o) / s
            g += (s / self.tick) * (_g(z_hi) - _g(z_hi - t) - _g(z_lo) + _g(z_lo - t))
        g /= len(self.offsets)
        right = np.clip(self.hi - d, 0.0, self.w)
        left = self.w - right
        mu = math.exp(min(log_bl, 700.0)) * left + math.exp(min(log_br, 700.0)) * right + math.exp(min(log_a, 700.0)) * g
        return np.maximum(mu, 1e-300)

    def nll(self, theta: np.ndarray) -> float:
        mu = self.mu(theta)
        return float(np.sum(mu - self.n * np.log(mu)))


def pair_pvalue(pair: PairData) -> float | None:
    """The grouped Pearson p-value of the forward model's shape on one histogram, or None when
    the region holds too few coincidences to check."""
    a = np.asarray(pair.a, dtype=np.int64)
    b = np.asarray(pair.b, dtype=np.int64)
    if a.size == 0 or b.size == 0:
        return None
    tick = tag_tick(a, b)
    shift = int(round(pair.declared_delay_ps / tick)) * tick
    half = int(pair.window_ps)
    half = -(-half // tick) * tick
    sigma0 = pair.jitter_ps if pair.jitter_ps and pair.jitter_ps > 0 else max(half / 10, tick)
    width = max(tick, int(sigma0 / 3) // tick * tick)
    # Keep the bin count bounded for wide windows against narrow peaks.
    width = max(width, -(-(2 * half) // 4000) // tick * tick or tick)
    hist = cross_correlation(a, b - shift, width, half).astype(float)
    if hist.sum() < MIN_COUNTS:
        return None
    edges = -half + width * np.arange(hist.size + 1, dtype=float)
    edges[-1] = min(edges[-1], float(half))
    model = _ShapeModel(hist, edges, tick, pair.offsets_ps)
    d0 = float(pair.delay_ps) + (pair.declared_delay_ps - shift)
    far = np.abs(0.5 * (edges[:-1] + edges[1:]) - d0) > 5 * sigma0 + max((abs(o) for o in pair.offsets_ps), default=0.0)
    bg = float(np.median(hist[far])) / width if far.any() else 0.5 / (2 * half)
    bg = max(bg, 0.5 / (2 * half))
    peak = max(float(hist.sum()) - bg * 2 * half, 1.0)
    total = float(hist.sum()) + 1.0
    lb = math.log(1e-6 / (2 * half))
    bounds = [
        (lb, math.log(10 * total / width)),
        (0.0, math.log(10 * total)),
        (-float(half), float(half)),
        (math.log(max(tick / 100, 0.01)), math.log(2 * half)),
        (lb, math.log(10 * total / width)),
    ]
    best = None
    for s_start in dict.fromkeys((sigma0, sigma0 * 2, sigma0 / 2)):
        x0 = np.array([math.log(bg), math.log(peak), d0, math.log(max(s_start, tick / 50)), math.log(bg)])
        x0 = np.clip(x0, [lo for lo, _ in bounds], [hi for _, hi in bounds])
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            r = optimize.minimize(model.nll, x0, method="L-BFGS-B", bounds=bounds, options={"ftol": 1e-14, "maxiter": 2000})
        if best is None or (np.isfinite(r.fun) and r.fun < best.fun):
            best = r
    return _grouped_pearson(model.n, model.mu(best.x))


def run_check(pairs: Sequence[PairData], alpha: float = MISFIT_ALPHA) -> CheckResult:
    pvalues = {}
    for pair in pairs:
        p = pair_pvalue(pair)
        if p is not None:
            pvalues[pair.label] = p
    return CheckResult(pvalues, alpha)


def enforce(result: CheckResult, analysis: str, policy: Policy) -> tuple[str, ...]:
    """Refuse a flagged analysis unless the policy allows the override; return its labels."""
    if not result.flagged:
        return ()
    if OVERRIDE not in policy.rule_for(analysis).allowed_corrections:
        raise Misspecified(
            f"the forward model does not fit the analysed coincidences ({result.worst()}: predictive check "
            f"p = {result.pvalue:.3g} after Bonferroni over {len(result.pvalues)} histograms): something in the "
            "window that the model lacks (crosstalk, an IRF tail, drift) would bias the corrected value. "
            f"Override with Policy.allow({analysis!r}, {OVERRIDE!r}, reason)."
        )
    return (FORWARD_MISFIT_ACCEPTED,)


def delay_start(calibrations, a_id, b_id) -> tuple[float, float | None]:
    """(delay of b after a, jitter) from the pair's delay-jitter artifact, or (0, None)."""
    from dadl.calibration.delay import KIND

    if calibrations is None:
        return 0.0, None
    for m in calibrations.members:
        if m.kind == KIND and set(m.channels) == {a_id, b_id}:
            cal = m.calibration
            sign = 1.0 if tuple(cal.channels) == (a_id, b_id) else -1.0
            return sign * float(cal.estimates["delay_ps"].value), float(cal.estimates["jitter_sigma_ps"].value)
    return 0.0, None
