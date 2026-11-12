"""Dark count rate, afterpulsing and dead time from one dark acquisition.

Fitted jointly from the inter-arrival gaps of a single channel with no light
(shutter closed or source off), as one artifact of kind
`dark-afterpulse-deadtime` (docs/design.md, "Joint posteriors and data
reuse").

Method `interarrival-mle` (v4):

- **Dead time.** No gap is shorter than the (effective, detector plus TDC)
  dead time, so its estimate is the shortest gap m. m - tau is a whole number
  D of tags' ticks with P(D <= k) = 1 - exp(-lambda (k + 1/2)), lambda =
  n f(0) tick, f(0) the gap density just past the dead time. The interval
  uses these exact quantiles averaged over the posterior of f(0); the point
  estimate is the median-unbiased m - ln(2)/(n f(0)). When afterpulsing is
  not detected it is not ruled out either, and an undetected afterpulse
  usually supplies the shortest gap, so the interval then runs up to m. (v1
  ignored the ticks: 89% coverage with fast afterpulsing. v2-v3 treated f(0)
  as known and the dark-count law as certain without detected afterpulsing:
  as low as 0% coverage with a few hundred gaps.)
- **Dark rate, afterpulse probability and time constant.** Maximum
  likelihood on the gaps past the dead time under `DarkAfterpulseModel`, with
  a Laplace posterior in (log r, logit p, log tau_a) for the draws and the
  rate. Below PROFILE_BELOW gaps the intervals for p and tau_a come from the
  profile likelihood (chi-square(1) at 95%), which is accurate with few counts
  where the Gaussian shape is not (tau_a covered 82% at 50 gaps).
- **No afterpulsing detected.** The time constant is searched from 1 ns to
  a tenth of the mean gap. If the likelihood ratio against a pure Poisson
  model is below the 99.9% chi-square(2) point (a stringent level, because of
  that search), or the time constant fits to the edge of the range, the dark
  rate comes from its exact Gamma posterior, the time constant is marked
  unidentified, and p is reported from 0 to the larger of its profile
  interval's upper end and the largest upper limit over a grid of time
  constants. Undetected samples hold fewer afterpulses than expected, so a
  less conservative limit undercovers there.
- **Few counts.** At least MIN_GAPS (20) gaps are needed; below a few hundred
  the afterpulse parameters are rarely identified, but the rate and dead time
  keep their coverage (probed down to 12 gaps).
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
from dadl.data import ChannelId, tag_tick, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.models import DarkAfterpulseModel
from dadl.stats import Interval, Laplace, expit, gamma_posterior, numerical_hessian

KIND = "dark-afterpulse-deadtime"
METHOD = "interarrival-mle"
METHOD_VERSION = "4"
MIN_GAPS = 20
#: Posterior samples of f(0) averaged over in the dead-time interval.
F0_SAMPLES = 4000
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


#: Below this many gaps the afterpulse intervals come from the profile likelihood; above it the
#: Laplace (Gaussian) intervals are as accurate (93-95% coverage in the regime probes) and faster.
PROFILE_BELOW = 5000
#: Profile-likelihood threshold for a two-sided 95% interval: chi-square(1) at 95%.
PROFILE_TARGET = 3.841
#: logit(p) limits of the profile search: p from about 3e-7 (reported as 0) to 0.9975.
LOGIT_P_RANGE = (-15.0, 6.0)


def _profile(s: np.ndarray, r0: float, index: int, value: float, warm: np.ndarray | None) -> float:
    """Smallest negative log-likelihood with coordinate `index` of theta fixed at `value`."""
    lo_tau, hi_tau = _tau_bounds(s)
    bounds = [(math.log(r0) - 5, math.log(r0) + 5), LOGIT_P_RANGE, (math.log(lo_tau), math.log(hi_tau))]
    free = [j for j in range(3) if j != index]

    def full(x: np.ndarray) -> np.ndarray:
        th = np.empty(3)
        th[index] = value
        th[free] = x
        return th

    def nll(x: np.ndarray) -> float:
        return -_MODEL.log_likelihood(s, full(x))

    def grad(x: np.ndarray) -> np.ndarray:
        return _MODEL.gradient(s, full(x))[free]

    starts = []
    if warm is not None:
        starts.append(np.clip(np.asarray(warm, dtype=float)[free], [bounds[j][0] for j in free], [bounds[j][1] for j in free]))
    if index == 1:
        starts += [np.array([math.log(r0), z]) for z in np.linspace(bounds[2][0], bounds[2][1], 3)]
    else:
        starts += [np.array([math.log(r0), y]) for y in (math.log(0.01 / 0.99), math.log(0.1 / 0.9))]
    best = math.inf
    for x0 in starts:
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            res = optimize.minimize(nll, x0, jac=grad, method="L-BFGS-B", bounds=[bounds[j] for j in free],
                                    options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500})
        if np.isfinite(res.fun):
            best = min(best, float(res.fun))
    return best


def _profile_interval(
    s: np.ndarray, r0: float, index: int, centre: float, nll_min: float, limits: tuple[float, float],
    warm: np.ndarray | None,
) -> tuple[float, float]:
    """{v : 2 (profile(v) - nll_min) <= PROFILE_TARGET}, clipped to `limits` (in fit coordinates)."""

    def excess(v: float) -> float:
        return 2.0 * (_profile(s, r0, index, v, warm) - nll_min) - PROFILE_TARGET

    lo_lim, hi_lim = limits
    centre = min(max(centre, lo_lim), hi_lim)
    hi = hi_lim if excess(hi_lim) < 0 else optimize.brentq(excess, centre, hi_lim, xtol=1e-3)
    lo = lo_lim if centre <= lo_lim or excess(lo_lim) < 0 else optimize.brentq(excess, lo_lim, centre, xtol=1e-3)
    return lo, hi


def _dead_time(
    m: int,
    n: int,
    f0: float,
    f0_samples: np.ndarray,
    rng: np.random.Generator,
    n_draws: int,
    tick: int = 1,
    afterpulsing_ruled_out: bool = True,
) -> tuple[Interval, np.ndarray]:
    """Interval and draws for the dead time from the shortest gap m.

    m - tau is a whole number D of ticks with P(D <= k) = 1 - exp(-lambda (k + 1/2)),
    lambda = n f(0) tick. f(0) is itself estimated, so the law is averaged over
    `f0_samples` from its posterior: with few gaps f(0) is loose, and treating it
    as known made the interval far too narrow.

    When afterpulsing was not detected it is not ruled out either: an afterpulse
    below the detection threshold usually supplies the shortest gap, within
    nanoseconds of the dead time. The interval then runs up to m itself and puts
    all 5% in the lower tail, which is exact without afterpulsing and conservative
    with it. (Method v2 used the dark-count law alone there and covered 0% of
    runs with undetected afterpulsing.)
    """
    lam = n * np.asarray(f0_samples, dtype=float) * tick  # per tick, one per posterior sample

    def cdf(y: float) -> float:
        return float(np.mean(-np.expm1(-lam * y)))

    def quantile(q: float) -> int:
        """Smallest whole k with P(D <= k) >= q (y = k + 1/2)."""
        if cdf(0.5) >= q:
            return 0
        hi = 1.0
        while cdf(hi) < q:
            hi *= 2.0
        y = optimize.brentq(lambda v: cdf(v) - q, 0.5, hi, xtol=1e-9)
        return max(0, math.ceil(y - 0.5 - 1e-12))

    if afterpulsing_ruled_out:
        k_lo = quantile(0.025)
        k_hi = quantile(0.95 if k_lo == 0 and cdf(0.5) > 0.025 else 0.975)
    else:
        k_lo, k_hi = 0, quantile(0.95)
    low, high = max(0.0, m - k_hi * tick), float(m - k_lo * tick)
    # Median-unbiased when afterpulsing is detected; otherwise the maximum-likelihood value m,
    # since the dark-count law alone would put the estimate microseconds below m (often at 0).
    point = max(0.0, m - math.log(2.0) / (n * f0)) if afterpulsing_ruled_out else float(m)
    pick = rng.integers(0, lam.size, n_draws)
    draws = np.maximum(0.0, m - tick * np.floor(rng.exponential(1.0 / lam[pick]) + 0.5))
    return Interval(min(max(point, low), high), low, high), draws



def _p_upper_limit(s: np.ndarray, nll_null: float, log_r0: float) -> tuple[float, float, float]:
    """Upper limit on p, the largest over a grid of tau_a (r profiled), at the two-sided 95% threshold.

    Taking the largest limit over the grid keeps it conservative when afterpulsing went
    undetected: those are the samples with fewer afterpulses than expected, and a limit
    profiled over tau_a as well undercovered there (41% at 200 gaps).

    Returns (limit, tau_min, tau_max): the limit holds for afterpulse time
    constants in [tau_min, tau_max] ps.
    """
    target = PROFILE_TARGET
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

    # The afterpulse probability gets one profile-likelihood interval whether or not afterpulsing
    # is detected, so the choice between an interval and a limit never depends on the data
    # (switching on detection undercovered: 82-90% when afterpulsing went undetected).
    use_profile = n < PROFILE_BELOW
    nll_min = min(nll_null, float(res.fun)) if res is not None else nll_null
    p_centre = float(res.x[1]) if res is not None and res.fun <= nll_null else LOGIT_P_RANGE[0]
    if use_profile or not detected:
        p_lo, p_hi = _profile_interval(s, r0, 1, p_centre, nll_min, LOGIT_P_RANGE, res.x if res is not None else None)
        p_low = 0.0 if p_lo <= LOGIT_P_RANGE[0] else float(expit(p_lo))

    if detected and laplace is not None:
        theta = res.x
        est_r = laplace.interval(0, lambda x: math.exp(x) * PS_PER_S)
        if use_profile:
            est_p = Interval(float(expit(theta[1])), p_low, float(expit(p_hi)))
            lo_tau, hi_tau = _tau_bounds(s)
            z_lo, z_hi = _profile_interval(
                s, r0, 2, float(theta[2]), float(res.fun), (math.log(lo_tau), math.log(hi_tau)), theta
            )
            est_ta = Interval(math.exp(theta[2]), math.exp(z_lo), math.exp(z_hi))
            if z_lo <= math.log(lo_tau) + 1e-6 or z_hi >= math.log(hi_tau) - 1e-6:
                notes.append("the afterpulse time constant's interval reaches the edge of its search range")
        else:
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
        ul, t_lo, t_hi = _p_upper_limit(s, nll_null, math.log(r0))
        est_p = Interval(0.0, 0.0, max(ul, float(expit(p_hi))))
        est_ta = None
        unidentified = ("afterpulse_tau_ps",)
        notes.append(
            f"no afterpulsing detected (likelihood ratio {lr:.2f}); afterpulse_prob is reported with its 95% "
            f"profile-likelihood interval, an upper limit when it starts at 0, for afterpulse time constants "
            f"from {t_lo:.3g} to {t_hi:.3g} ps"
        )
        notes.append(
            "the dark rate counts every avalanche past the dead time, so undetected afterpulses (at most the "
            "upper limit on afterpulse_prob, as a fraction) are included in it"
        )
        converged = True

    f0 = _MODEL.density_at_zero(theta)
    # Posterior samples of f(0) for the dead-time law (fixed seed: the interval depends on the data only).
    if detected and laplace is not None:
        th = Laplace(laplace.mean, laplace.cov).draws(np.random.default_rng(0), F0_SAMPLES)
        r_s = np.exp(np.clip(th[:, 0], -80, 0))
        p_s = expit(np.clip(th[:, 1], -40, 40))
        ta_s = np.exp(np.clip(th[:, 2], -20, 45))
        f0_samples = r_s * (1 - p_s) + p_s * (r_s + 1 / ta_s)
    else:
        q = (np.arange(F0_SAMPLES) + 0.5) / F0_SAMPLES
        f0_samples = _st.gamma.ppf(q, n + 0.5, scale=1.0 / exposure)
    est_tau, tau_draws = _dead_time(
        m, n, f0, f0_samples, rng, n_draws, tag_tick(tags), afterpulsing_ruled_out=bool(detected)
    )
    if est_tau.low == 0.0:
        notes.append(
            f"the dead time is only bounded above (<= {est_tau.high:.6g} ps): too few gaps fall near it to "
            "bound it below; measure it with light on the detector"
        )

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
