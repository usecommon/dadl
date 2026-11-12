"""Dark count rate, afterpulsing and dead time from one dark acquisition.

Fitted jointly from the inter-arrival gaps of a single channel with no light
(shutter closed or source off), as one artifact of kind
`dark-afterpulse-deadtime` (docs/design.md, "Joint posteriors and data
reuse").

Method `interarrival-mle`:

- **Dead time.** No gap is shorter than the (effective, detector plus TDC)
  dead time, so its estimate is the shortest gap m. Its sampling
  distribution is m - tau ~ Exponential(n f(0)), with f(0) the fitted gap
  density just past the dead time; the point estimate is the median-unbiased
  m - ln(2)/(n f(0)) and the interval is [m - ln(40)/(n f(0)),
  m - ln(40/39)/(n f(0))] and the draws follow the same law.
- **Dark rate, afterpulse probability and time constant.** Maximum
  likelihood on the gaps past the dead time under `DarkAfterpulseModel`, with
  a Laplace posterior in (log r, logit p, log tau_a).
- **No afterpulsing detected.** The time constant is searched from 1 ns to
  a tenth of the mean gap. If the likelihood ratio against a pure Poisson
  model is below the 99.9% chi-square(2) point (a stringent level, because of
  that search), or the time constant fits to the edge of the range, the dark
  rate comes
  from its exact Gamma posterior, the afterpulse probability is reported as a
  one-sided 95% upper limit from the profile likelihood, and the time
  constant is marked unidentified.
- **Predictive check.** The fitted survival function maps the gaps to
  uniforms; a chi-square test on 20 equal bins (19 degrees of freedom,
  conservative because the fit is unbinned) gives the p-value recorded in the
  artifact.

Phase 2 limits: one exponential afterpulse component; tags must come from
one uninterrupted acquisition.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import optimize
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.base import N_DRAWS, Calibration
from dadl.data import ChannelId, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.models import DarkAfterpulseModel
from dadl.stats import Interval, Laplace, expit, gamma_posterior, numerical_hessian

KIND = "dark-afterpulse-deadtime"
METHOD = "interarrival-mle"
METHOD_VERSION = "1"
MIN_GAPS = 200
PPC_BINS = 20
PS_PER_S = 1e12
#: Afterpulse time constants the fit considers, in ps. Shorter afterpulses
#: cannot be told apart from dead time; the upper end is also capped at a
#: tenth of the mean gap, where afterpulses blur into the dark counts.
TAU_SEARCH_PS = (1e3, 1e7)
#: Likelihood-ratio threshold for claiming afterpulsing: chi-square(2) at 99.9%.
#: A stringent level, because the time constant is searched over a range.
DETECTION_LR = 13.82

UNITS = {
    "dark_rate_cps": "1/s",
    "afterpulse_prob": "probability per avalanche",
    "afterpulse_tau_ps": "ps",
    "dead_time_ps": "ps",
}

_MODEL = DarkAfterpulseModel()


def _tau_bounds(s: np.ndarray) -> tuple[float, float]:
    """Afterpulse time constants the data can identify, in ps: from 1 ns to a tenth of the mean gap."""
    return TAU_SEARCH_PS[0], min(TAU_SEARCH_PS[1], float(s.mean()) / 10)


def _fit(s: np.ndarray, r0: float) -> optimize.OptimizeResult | None:
    """Multi-start bounded fit of the negative log-likelihood; None if no start is usable."""
    lo_tau, hi_tau = _tau_bounds(s)
    if hi_tau <= lo_tau:
        return None

    def nll(theta: np.ndarray) -> float:
        return -_MODEL.log_likelihood(s, theta)

    def grad(theta: np.ndarray) -> np.ndarray:
        return _MODEL.gradient(s, theta)

    bounds = [(math.log(r0) - 5, math.log(r0) + 5), (-30.0, 10.0), (math.log(lo_tau), math.log(hi_tau))]
    best = None
    for log_tau0 in np.linspace(bounds[2][0], bounds[2][1], 4):
        for p0 in (0.01, 0.1):
            x0 = np.array([math.log(r0), math.log(p0 / (1 - p0)), log_tau0])
            with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
                res = optimize.minimize(nll, x0, jac=grad, method="L-BFGS-B", bounds=bounds,
                                        options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 1000})
            if np.isfinite(res.fun) and (best is None or res.fun < best.fun):
                best = res
    if best is not None:
        best.at_bound = bool(min(best.x[2] - bounds[2][0], bounds[2][1] - best.x[2]) < 1e-3)
    return best


def _ppc_pvalue(s: np.ndarray, theta: np.ndarray) -> float:
    """Chi-square test of the probability-integral transform on PPC_BINS equal bins.

    The parameters are fitted to the unbinned gaps, so the statistic is not
    chi-square(bins - 1 - fitted) but lies between that and chi-square(bins - 1)
    (Chernoff and Lehmann, 1954). Using bins - 1 keeps the test conservative.
    """
    u = 1.0 - _MODEL.survival(s, theta)
    counts, _ = np.histogram(u, bins=PPC_BINS, range=(0.0, 1.0))
    expected = s.size / PPC_BINS
    chi2 = float(np.sum((counts - expected) ** 2 / expected))
    return float(_st.chi2.sf(chi2, PPC_BINS - 1))


def _dead_time(m: int, n: int, f0: float, rng: np.random.Generator, n_draws: int) -> tuple[Interval, np.ndarray]:
    rate = n * f0  # per ps
    low = max(0.0, m - math.log(40.0) / rate)
    high = m - math.log(40.0 / 39.0) / rate
    draws = np.maximum(0.0, m - rng.exponential(1.0 / rate, n_draws))
    return Interval(max(0.0, m - math.log(2.0) / rate), low, high), draws


def _p_upper_limit(s: np.ndarray, nll_null: float, log_r0: float) -> tuple[float, float, float]:
    """One-sided 95% upper limit on p, the largest over a grid of tau_a (r profiled).

    Returns (limit, tau_min, tau_max): the limit holds for afterpulse time
    constants in [tau_min, tau_max] ps.
    """
    target = 2.706  # chi-square(1) 90% two-sided point = one-sided 95%
    lo_tau, hi_tau = _tau_bounds(s)
    if hi_tau <= lo_tau:
        return 1.0, lo_tau, lo_tau

    def nll_profile_r(logit_p: float, log_ta: float) -> float:
        res = optimize.minimize_scalar(
            lambda x: -_MODEL.log_likelihood(s, np.array([x, logit_p, log_ta])),
            bounds=(log_r0 - 1.0, log_r0 + 1.0),
            method="bounded",
            options={"xatol": 1e-7},
        )
        return float(res.fun)

    limit = 0.0
    for log_ta in np.linspace(math.log(lo_tau), math.log(hi_tau), 12):
        def excess(y: float, log_ta: float = log_ta) -> float:
            return 2.0 * (nll_profile_r(y, log_ta) - nll_null) - target

        lo, hi = -15.0, 6.0
        if excess(hi) < 0:
            return 1.0, lo_tau, hi_tau
        y = lo if excess(lo) > 0 else optimize.brentq(excess, lo, hi, xtol=1e-3)
        limit = max(limit, float(expit(y)))
    return limit, lo_tau, hi_tau


def calibrate_dark(
    tags: np.ndarray,
    channel: ChannelId,
    *,
    backend: str = "hardware",
    seed: int | None = None,
    n_draws: int = N_DRAWS,
    conditions: dict[str, float] | None = None,
) -> Calibration:
    """Fit dark rate, afterpulsing and dead time from one channel's dark tags (int64 ps, sorted)."""
    tags = _checks.timestamps(tags, "tags")
    _checks.check_sorted(tags, "tags")
    _checks.check_span(tags, "tags")
    gaps = np.diff(tags)
    if gaps.size < MIN_GAPS:
        raise InsufficientCounts(f"need at least {MIN_GAPS} gaps between dark counts, got {gaps.size}")
    rng = np.random.default_rng(seed)

    m = int(gaps.min())
    # The shortest gap defines m; it sits at s = 0 by construction, so it is left
    # out of the shape fit (otherwise it reads as a spurious, very fast afterpulse).
    s = np.delete(gaps, int(np.argmin(gaps))).astype(float) - m
    n = s.size
    exposure = float(s.sum())
    if exposure <= 0:
        raise InsufficientCounts("every gap equals the dead time; the dark rate is not identified")
    r0 = n / exposure
    nll_null = -(n * math.log(r0) - r0 * exposure)

    res = _fit(s, r0)
    notes: list[str] = []
    unidentified: tuple[str, ...] = ()
    lr = 0.0 if res is None else 2.0 * (nll_null - float(res.fun))
    detected = res is not None and lr > DETECTION_LR
    if detected and res.at_bound:
        notes.append("the afterpulse time constant fitted to the edge of its search range; treated as not identified")
        detected = False

    laplace = None
    if detected:
        try:
            hess = numerical_hessian(lambda th: _MODEL.gradient(s, th), res.x)
            laplace = Laplace.from_hessian(res.x, hess)
        except ValueError as exc:
            notes.append(f"afterpulse fit not identified ({exc}); reporting the no-afterpulse model")
            detected = False

    if detected and laplace is not None:
        theta = res.x
        est_r = laplace.interval(0, lambda x: math.exp(x) * PS_PER_S)
        est_p = laplace.interval(1, lambda y: float(expit(y)))
        est_ta = laplace.interval(2, math.exp)
        d = laplace.draws(rng, n_draws)
        r_draws = np.exp(d[:, 0]) * PS_PER_S
        p_draws = expit(d[:, 1])
        ta_draws = np.exp(d[:, 2])
        # BFGS often stops with "precision loss" at the optimum; judge by the gradient instead.
        converged = bool(res.success) or float(np.max(np.abs(_MODEL.gradient(s, res.x)))) < 1e-3 * math.sqrt(n)
    else:
        theta = np.array([math.log(r0), -40.0, 0.0])
        est_r, r_draws = gamma_posterior(n, exposure, rng, n_draws)
        est_r = Interval(est_r.value * PS_PER_S, est_r.low * PS_PER_S, est_r.high * PS_PER_S)
        r_draws = r_draws * PS_PER_S
        p_hi, t_lo, t_hi = _p_upper_limit(s, nll_null, math.log(r0))
        est_p = Interval(0.0, 0.0, p_hi)
        est_ta = None
        unidentified = ("afterpulse_tau_ps",)
        notes.append(
            f"no afterpulsing detected (likelihood ratio {lr:.2f}); afterpulse_prob is a one-sided 95% "
            f"upper limit for afterpulse time constants from {t_lo:.3g} to {t_hi:.3g} ps"
        )
        notes.append(
            "the dark rate counts every avalanche past the dead time, so undetected afterpulses (at most the "
            "upper limit on afterpulse_prob, as a fraction) are included in it"
        )
        converged = True

    f0 = _MODEL.density_at_zero(theta)
    est_tau, tau_draws = _dead_time(m, n, f0, rng, n_draws)

    estimates = {"dark_rate_cps": est_r, "afterpulse_prob": est_p}
    columns = [r_draws]
    names = ["dark_rate_cps"]
    if est_ta is not None:
        estimates["afterpulse_tau_ps"] = est_ta
        columns += [p_draws, ta_draws]
        names += ["afterpulse_prob", "afterpulse_tau_ps"]
    estimates["dead_time_ps"] = est_tau
    columns.append(tau_draws)
    names.append("dead_time_ps")

    return Calibration(
        kind=KIND,
        method=METHOD,
        method_version=METHOD_VERSION,
        channels=(channel,),
        estimates=estimates,
        units={k: UNITS[k] for k in estimates},
        samples=np.column_stack(columns),
        sample_params=tuple(names),
        converged=converged,
        ppc_pvalue=_ppc_pvalue(s, theta),
        inputs=((tags_sha256(tags), "raw-tags"),),
        backend=backend,
        conditions=conditions or {},
        unidentified=unidentified,
        notes=tuple(notes),
    )
