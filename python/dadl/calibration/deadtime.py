"""Series dead time (paralyzable then non-paralyzable) from illuminated acquisitions (phase 3, M4).

Kind `deadtime-series` (phase 3 plan, decision D4), method `series-gap-mle`.

Dark data cannot identify a paralyzable stage: paralysis needs arrivals during the dead time.
This calibration takes two or more acquisitions of one channel under steady light at different
rates (`dadl simulate --deadtime-scan` writes them) and fits the gap distribution of the chain
in the simulator and the design ("Detection process"): every avalanche candidate restarts a
paralyzable stage tau_p, and what passes it meets a non-paralyzable tau_np >= tau_p.

The scan needs steady light: an attenuated CW laser with its intensity stable to about 1-2% over
each run, not pulsed light or a mains-powered lamp. Slow intensity noise looks like less
paralysis: a +-20% 100 Hz modulation biased tau_p 10% low with the predictive check passing;
+-40% is flagged.

**The gap density.** For a Poisson input of rate n, measure u = gap - tau_np. The next recorded
avalanche is the first input after tau_np whose preceding input is at least tau_p earlier
(inputs during the non-paralyzable dead time still restart the paralyzable stage). Its density
is g(u) = a S(u), a = n exp(-n tau_p), where S(w) = 1 for w < tau_p and S'(w) = -a S(w - tau_p):

    S(w) = sum_{j <= w / tau_p} (-a (w - j tau_p))^j / j!

summed exactly up to K tau_p and continued beyond as its two real modes exp(-lambda w),
lambda = a exp(lambda tau_p) (Lambert W branches 0 and -1), matched to S and S' there. The
complex modes have decayed by about (a tau_p / 2 pi)^K; the second real mode has not when
n tau_p is near 1, where the two real roots merge. The distribution function is
1 - S(u + tau_p), so the bin probabilities are exact. tau_p = 0 is the hard-edge
(non-paralyzable) model.

**Afterpulsing.** With probability p a recorded avalanche adds one afterpulse candidate at
u = E ~ Exp(tau_a) (the simulator's model and the dark fit's). Before it arrives the gap density
is the one above; it is recorded if no input came in the tau_p before it; otherwise it restarts
the paralyzable stage, and the next record follows the same law from there:

    g(u) = g0(u) + p [ -g0(u)(1 - e^{-u/tau_a}) + e^{-u/tau_a} S(u) e^{-n tau_p} / tau_a
                       + int_0^{u - tau_p} f_E(w) (S(w + tau_p) - S(w) e^{-n tau_p}) q(u - w) dw ]

with q(t) = a S(t - tau_p) for t >= tau_p (a discrete convolution on the fit's grid, q averaged
over each cell so the profile is smooth in tau_p). p and tau_a are held at a dark artifact's
estimates (`afterpulse=`); their uncertainty is carried by refitting at shifted values (a linear
sensitivity, of tau_p or of its upper limit) and added to the intervals and draws. An afterpulse
still pending when the next gap starts is not modelled (a fraction about p n tau_a / (1 + n tau_a)
of gaps, recorded as `afterpulse_pending_max`); above 0.5% a note says tau_p is biased low.

**Fit.** tau_np is the shortest gap over all runs, with the dark fit's exact tick-quantile
interval (`dadl.calibration.dark._dead_time`) on the summed gap density at the edge. The gaps
past it are binned (fine uniform bins, then geometric ones) and fitted by multinomial maximum
likelihood in (tau_p, n_1 ... n_K); given tau_p the runs are independent, so the profile
likelihood of tau_p is a sum of one-dimensional fits. The 95% interval for tau_p is the profile
likelihood's (chi-square(1)); when it reaches tau_p = 0 the hard-edge model is not rejected and
tau_p is an upper limit (unidentified; draws uniform below it, a stand-in). tau_p's point, upper
end and draws are capped at tau_np's. The likelihood ratio against the hard-edge model is
recorded (`lr_hard_edge`).

**Predictive check.** A Pearson chi-square over all runs' bins, merged until each expects at
least 5 gaps, on (groups - 1) degrees of freedom (conservative for a fit on the fine bins; see
`dadl.forward_check`), and the interruption check of the other calibrations, each run on its own
bounds (Bonferroni over runs). The check groups thousands of fine bins on (groups - 1) degrees of
freedom, so it sees gross shape misfits (soft recovery, reset-window paralysis, strong
modulation) but has little power against a discrepancy confined to the tail or the edge: a few-%
misfit there can pass while biasing tau_p.

Stays out (phase 3 plan, M4): soft (continuous) recovery of SNSPDs and passive-quench efficiency
recovery curves. With enough gaps such a detector shows as a misfit; with few, or when its
recovery is close to fully paralyzable, it can pass with an effective tau_p near tau_np (README).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from scipy import optimize, signal, special
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.base import N_DRAWS, Calibration, draws_rng, with_interruption_check_runs
from dadl.calibration.dark import _dead_time
from dadl.data import ChannelId, tag_tick, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.stats import Interval

KIND = "deadtime-series"
METHOD = "series-gap-mle"
METHOD_VERSION = "2"
MIN_RUNS = 2
MIN_GAPS = 1000
#: Terms of the exact series for S; beyond K tau_p its two real modes.
K_TERMS = 8
#: Fine bins: at most this fraction of a run's mean gap, and of tau_a.
BINS_PER_MEAN_GAP = 200
BINS_PER_TAU_A = 20
#: The fine grid covers this many mean gaps (and the afterpulse kernel), geometric bins beyond.
UNIFORM_MEAN_GAPS = 5.0
AFTERPULSE_TAUS = 40.0
GEOMETRIC_RATIO = 1.05
MAX_UNIFORM_BINS = 6000
MIN_EXPECTED = 5.0
PROFILE_TARGET = 3.841
F0_SAMPLES = 4000
#: Above this fraction of gaps starting with an afterpulse still pending (not modelled), a note.
PENDING_MAX = 0.005
PS_PER_S = 1e12

UNITS = {"dead_time_p_ps": "ps", "dead_time_ps": "ps"}
DEFINITIONS = {
    "dead_time_p_ps": "paralyzable stage: an avalanche candidate passes only if the previous candidate (passed or not) "
    "is at least this much earlier",
    "dead_time_ps": "non-paralyzable stage after the paralyzable one: blocks this long after each recorded avalanche; "
    "the shortest recorded gap, so it includes the tagger's timing jitter",
}


def _modes(a: float, tau: float) -> tuple[float, float]:
    """The two real roots l1 <= l2 of lambda = a exp(lambda tau) (Lambert W branches 0 and -1).

    They merge at a tau = 1/e (n tau = 1), where both are 1/tau; scipy gives nan there.
    """
    z = max(-a * tau, -1.0 / math.e)
    l1 = float(-special.lambertw(z, 0).real / tau)
    l2 = float(-special.lambertw(z, -1).real / tau)
    if not (math.isfinite(l1) and math.isfinite(l2)):
        return 1.0 / tau, 1.0 / tau
    return l1, l2


def survival(w: np.ndarray, n: float, tau: float) -> np.ndarray:
    """S(w): 1 below tau, S' = -a S(w - tau) above (see the module docstring)."""
    w = np.asarray(w, dtype=float)
    if tau <= 0.0:
        return np.exp(-n * np.maximum(w, 0.0))
    a = n * math.exp(-n * tau)
    u0 = K_TERMS * tau

    def series(x: np.ndarray) -> np.ndarray:
        out = np.zeros_like(x)
        for j in range(K_TERMS + 1):
            r = x - j * tau
            out += np.where(r >= 0.0, (-a * np.maximum(r, 0.0)) ** j / math.factorial(j), 0.0)
        return out

    # Past u0 the two real modes, matched to S(u0) and S'(u0) = -a S(u0 - tau). The divided
    # difference is written with expm1 so it stays exact as the modes merge (n tau near 1).
    s0 = float(series(np.array([u0]))[0])
    ds0 = -a * float(series(np.array([u0 - tau]))[0])
    l1, l2 = _modes(a, tau)
    t = np.maximum(w - u0, 0.0)
    dl = l2 - l1
    if abs(dl) * max(float(t.max(initial=0.0)), 1.0) < 1e-9:
        dd = t * np.exp(-l1 * t)
    else:
        dd = np.exp(-l1 * t) * (-np.expm1(-dl * t)) / dl
    out = s0 * np.exp(-l1 * t) + (ds0 + l1 * s0) * dd
    near = w <= u0
    if near.any():
        out[near] = series(w[near])
    return out


def gap_cdf(u: np.ndarray, n: float, tau: float) -> np.ndarray:
    """P(gap - tau_np <= u) without afterpulsing."""
    return 1.0 - survival(np.asarray(u, dtype=float) + tau, n, tau)


class _Run:
    """One acquisition's gaps (past the shortest), binned for the fit."""

    def __init__(self, s: np.ndarray, tick: int, tau_a: float | None) -> None:
        self.n_gaps = s.size
        self.mean = float(s.mean()) if s.size else 1.0
        h = self.mean / BINS_PER_MEAN_GAP
        if tau_a:
            h = min(h, tau_a / BINS_PER_TAU_A)
        h = max(float(tick), math.floor(h / tick) * tick)
        top = float(s.max()) + tick
        span = max(UNIFORM_MEAN_GAPS * self.mean, AFTERPULSE_TAUS * tau_a if tau_a else 0.0)
        n_uni = int(min(MAX_UNIFORM_BINS, max(1, math.ceil(min(span, top) / h))))
        edges = list(-0.5 * tick + h * np.arange(n_uni + 1))
        e = edges[-1]
        while e < top:
            e = max(e * GEOMETRIC_RATIO, e + h)
            edges.append(e)
        self.edges = np.array(edges)
        self.h = h
        self.n_uni = n_uni
        self.lo = np.maximum(self.edges[:-1], 0.0)
        self.hi = self.edges[1:]
        self.mid = 0.5 * (self.lo + self.hi)
        self.width = self.hi - self.lo
        self.counts = np.histogram(s, self.edges)[0].astype(float)
        self.nz = self.counts > 0

    def probs(self, n: float, tau: float, p: float = 0.0, tau_a: float | None = None) -> np.ndarray:
        cdf = gap_cdf(self.edges.clip(min=0.0), n, tau)
        prob = np.diff(cdf)
        if p > 0 and tau_a:
            prob = prob + p * self._afterpulse(n, tau, tau_a) * self.width
        return np.maximum(prob, 1e-300)

    def _afterpulse(self, n: float, tau: float, tau_a: float) -> np.ndarray:
        """Delta(u) at the bin midpoints (module docstring: the bracket multiplying p)."""
        u = self.mid
        e_t = math.exp(-n * tau)
        a = n * e_t
        s_u = survival(u, n, tau)
        decay = np.exp(-u / tau_a)
        out = -a * s_u * (1.0 - decay) + decay * s_u * e_t / tau_a
        h = self.h
        n_w = int(min(self.n_uni, math.ceil(AFTERPULSE_TAUS * tau_a / h)))
        w = (np.arange(n_w) + 0.5) * h
        f = np.exp(-w / tau_a) / tau_a * (survival(w + tau, n, tau) - survival(w, n, tau) * e_t) * h
        lags = np.arange(self.n_uni) * h

        def q_cell(t: np.ndarray) -> np.ndarray:
            # q averaged over the cell [t - h/2, t + h/2]: a S(t' - tau) for t' >= tau integrates to
            # S(lo) - S(hi) (S = 1 below tau), exact and smooth in tau. The point value snapped the
            # jump at tau to the grid and made the profile jagged (council round 1).
            lo, hi = np.maximum(t - 0.5 * h, 0.0), np.maximum(t + 0.5 * h, 0.0)
            return (survival(lo, n, tau) - survival(hi, n, tau)) / h

        conv = signal.fftconvolve(f, q_cell(lags))[: self.n_uni]
        out[: self.n_uni] += conv
        if u.size > self.n_uni:
            t = u[self.n_uni :, None] - w[None, :]
            out[self.n_uni :] += q_cell(t.ravel()).reshape(t.shape) @ f
        return out

    def nll(self, n: float, tau: float, p: float = 0.0, tau_a: float | None = None) -> float:
        prob = self.probs(n, tau, p, tau_a)
        return float(-np.sum(self.counts[self.nz] * np.log(prob[self.nz])))


class _Fit:
    def __init__(self, runs: list[_Run], p: float, tau_a: float | None) -> None:
        self.runs, self.p, self.tau_a = runs, p, tau_a
        self._n_cache: list[float] = [1.0 / r.mean for r in runs]

    def best_n(self, run_index: int, tau: float) -> tuple[float, float]:
        """The run's input rate maximizing the likelihood at this tau_p, by Newton steps in log n
        (the likelihood is close to quadratic there: its curvature is about the number of gaps)."""
        run = self.runs[run_index]

        def f(y: float) -> float:
            return run.nll(math.exp(y), tau, self.p, self.tau_a)

        y = math.log(self._n_cache[run_index])
        d = max(1e-4, 3.0 / math.sqrt(run.n_gaps))
        f0 = f(y)
        for _ in range(12):
            fm, fp = f(y - d), f(y + d)
            curv = (fp - 2.0 * f0 + fm) / (d * d)
            if not (np.isfinite(curv) and curv > 0):
                y = y - d if fm < fp else y + d
                f0 = min(fm, fp)
                continue
            step = float(np.clip(-(fp - fm) / (2.0 * d * curv), -0.5, 0.5))
            y_new = y + step
            f_new = f(y_new)
            while f_new > f0 and abs(step) > 1e-12:
                step *= 0.5
                y_new = y + step
                f_new = f(y_new)
            if f_new <= f0:
                y, f0 = y_new, f_new
            # A few parts in 1e5 of the rate: far inside its statistical error (about 1/sqrt(gaps)).
            if abs(step) < 2e-5:
                break
            d = max(min(d, 2.0 * abs(step)), 1e-6)
        if not np.isfinite(f0):
            raise InsufficientCounts("the gap fit failed")
        self._n_cache[run_index] = math.exp(y)
        return math.exp(y), f0

    def profile(self, tau: float) -> float:
        return sum(self.best_n(i, tau)[1] for i in range(len(self.runs)))

    def rates(self, tau: float) -> list[float]:
        return [self.best_n(i, tau)[0] for i in range(len(self.runs))]


def _fit_tau(fit: _Fit, tau_hi: float) -> tuple[float, float, float]:
    """(tau_p MLE, its NLL, the NLL at tau_p = 0)."""
    nll0 = fit.profile(0.0)
    grid = np.geomspace(max(tau_hi * 1e-3, 1.0), tau_hi, 12)
    vals = [fit.profile(float(t)) for t in grid]
    k = int(np.argmin(vals))
    if vals[k] >= nll0:
        return 0.0, nll0, nll0
    lo = math.log(grid[max(k - 1, 0)])
    hi = math.log(grid[min(k + 1, grid.size - 1)])
    res = optimize.minimize_scalar(lambda y: fit.profile(math.exp(y)), bounds=(lo, hi), method="bounded",
                                   options={"xatol": 1e-5})
    tau, val = (math.exp(float(res.x)), float(res.fun)) if res.fun < vals[k] else (float(grid[k]), vals[k])
    if val >= nll0:
        return 0.0, nll0, nll0
    return tau, val, nll0


def _refit_near(fit: _Fit, tau: float, tau_hi: float) -> float:
    """tau_p MLE near a previous one (the afterpulse sensitivity refits)."""
    lo, hi = math.log(max(tau * 0.7, 1.0)), math.log(min(tau * 1.3, tau_hi))
    res = optimize.minimize_scalar(lambda y: fit.profile(math.exp(y)), bounds=(lo, hi), method="bounded",
                                   options={"xatol": 2e-4})
    return math.exp(float(res.x))


def _profile_interval(fit: _Fit, tau: float, nll_min: float, nll0: float, tau_hi: float) -> tuple[float, float]:
    target = nll_min + PROFILE_TARGET / 2.0

    def excess(t: float) -> float:
        return fit.profile(t) - target

    if tau <= 0.0 or nll0 <= target:
        low = 0.0
    else:
        low = optimize.brentq(lambda y: excess(math.exp(y)), math.log(max(tau * 1e-6, 1e-3)), math.log(tau), xtol=1e-4) \
            if excess(max(tau * 1e-6, 1e-3)) > 0 else 0.0
        low = math.exp(low) if low else 0.0
    start = max(tau, 1.0)
    hi = start * 1.5
    while excess(hi) < 0 and hi < tau_hi:
        hi *= 1.5
    high = math.exp(optimize.brentq(lambda y: excess(math.exp(y)), math.log(start), math.log(min(hi, tau_hi)), xtol=1e-4)) \
        if excess(min(hi, tau_hi)) > 0 else tau_hi
    return low, high


def _pearson(runs: list[_Run], probs: list[np.ndarray]) -> float:
    obs, exp = [], []
    for run, pr in zip(runs, probs):
        mu = pr * run.n_gaps
        acc_n = acc_mu = 0.0
        for ni, mi in zip(run.counts, mu):
            acc_n += ni
            acc_mu += mi
            if acc_mu >= MIN_EXPECTED:
                obs.append(acc_n)
                exp.append(acc_mu)
                acc_n = acc_mu = 0.0
        if exp:
            obs[-1] += acc_n
            exp[-1] += acc_mu
    dof = len(exp) - 1
    if dof < 1:
        return 1.0
    o, e = np.array(obs), np.array(exp)
    return float(_st.chi2.sf(float(np.sum((o - e) ** 2 / e)), dof))


def _afterpulse_params(afterpulse: Calibration | None, notes: list[str]):
    """(p, tau_a, draws of (p, tau_a) or None) from a dark artifact."""
    if afterpulse is None:
        notes.append("no afterpulse artifact given: afterpulsing is taken as absent")
        return 0.0, None, None
    est = afterpulse.estimates
    if "afterpulse_tau_ps" not in est or "afterpulse_prob" in afterpulse.stand_in:
        limit = est["afterpulse_prob"].high if "afterpulse_prob" in est else float("nan")
        notes.append(
            f"the dark artifact detects no afterpulsing (afterpulse_prob <= {limit:.3g}): it is taken as absent"
        )
        return 0.0, None, None
    p, tau_a = float(est["afterpulse_prob"].value), float(est["afterpulse_tau_ps"].value)
    cols = list(afterpulse.sample_params)
    draws = np.column_stack([afterpulse.samples[:, cols.index("afterpulse_prob")],
                             afterpulse.samples[:, cols.index("afterpulse_tau_ps")]])
    notes.append(f"afterpulsing held at the dark artifact's p = {p:.4g}, tau_a = {tau_a:.4g} ps; its uncertainty is "
                 "propagated by a linear sensitivity")
    return p, tau_a, draws


def calibrate_deadtime_series(
    runs: Sequence[np.ndarray],
    channel: ChannelId,
    *,
    afterpulse: Calibration | None = None,
    backend: str = "hardware",
    seed: int | None = None,
    n_draws: int = N_DRAWS,
    conditions: dict[str, float] | None = None,
) -> Calibration:
    """Fit tau_p and tau_np from illuminated acquisitions of one channel at two or more rates.

    `runs` are sorted int64 ps tags, one array per acquisition, under steady light at
    different rates. `afterpulse` is the channel's dark artifact (`calibrate_dark`), whose
    afterpulse probability and time constant are held fixed with their uncertainty propagated.
    """
    if len(runs) < MIN_RUNS:
        raise ValueError(f"need at least {MIN_RUNS} acquisitions at different rates; paralysis shows in how the "
                         "gap distribution changes with rate")
    tags = []
    for i, t in enumerate(runs):
        t = _checks.timestamps(t, f"runs[{i}]")
        _checks.check_sorted(t, f"runs[{i}]")
        _checks.check_span(t, f"runs[{i}]")
        tags.append(t)
    gaps = [np.diff(t) for t in tags]
    short = [i for i, g in enumerate(gaps) if g.size < MIN_GAPS]
    if short:
        raise InsufficientCounts(f"need at least {MIN_GAPS} gaps per acquisition; run(s) {short} have fewer")
    tick = tag_tick(*tags)
    notes: list[str] = []
    rng = draws_rng(seed, *tags)
    p, tau_a, ap_draws = _afterpulse_params(afterpulse, notes)

    m = int(min(int(g.min()) for g in gaps))
    which = int(np.argmin([g.min() for g in gaps]))
    shifted = []
    for i, g in enumerate(gaps):
        s = g.astype(float) - m
        if i == which:
            s = np.delete(s, int(np.argmin(s)))
        shifted.append(s)
    binned = [_Run(s, tick, tau_a) for s in shifted]
    fit = _Fit(binned, p, tau_a)
    tau_hi = float(m)  # tau_p <= tau_np
    tau, nll_min, nll0 = _fit_tau(fit, tau_hi)
    lr = 2.0 * (nll0 - nll_min)
    low, high = _profile_interval(fit, tau, nll_min, nll0, tau_hi)
    rates = fit.rates(tau)
    identified = low > 0.0

    # Afterpulse uncertainty: d tau_p / d(p, tau_a) by refitting at shifted values; for an upper
    # limit, the sensitivity of the limit itself (whatever the MLE).
    shift_sd = 0.0
    sens = np.zeros(2)
    if ap_draws is not None:
        sd = ap_draws.std(axis=0)
        for k in range(2):
            if sd[k] <= 0:
                continue
            step = sd[k]
            p_k, ta_k = (p + step, tau_a) if k == 0 else (p, tau_a + step)
            other = _Fit(binned, p_k, ta_k)
            other._n_cache = list(fit._n_cache)
            if identified:
                sens[k] = (_refit_near(other, tau, tau_hi) - tau) / step
            else:
                t_k, nll_k, nll0_k = _fit_tau(other, tau_hi)
                sens[k] = (_profile_interval(other, t_k, nll_k, nll0_k, tau_hi)[1] - high) / step
        cov = np.cov(ap_draws.T)
        shift_sd = float(math.sqrt(max(sens @ cov @ sens, 0.0)))
        if shift_sd > 0:
            what = "tau_p" if identified else "tau_p's upper limit"
            notes.append(f"afterpulse uncertainty adds {shift_sd:.3g} ps (1 sd) to {what}")

    unidentified: tuple[str, ...] = ()
    stand_in: tuple[str, ...] = ()
    if identified:
        # Combine the profile interval (statistical) with the afterpulse shift, side by side.
        lo_hw, hi_hw = tau - low, high - tau
        lo_hw = math.hypot(lo_hw, 1.96 * shift_sd)
        hi_hw = math.hypot(hi_hw, 1.96 * shift_sd)
        est_p = Interval(tau, max(0.0, tau - lo_hw), tau + hi_hw)
        # Draws: a split normal matching the interval's two halves (in ps).
        z = rng.standard_normal(n_draws)
        tau_draws = np.maximum(0.0, tau + np.where(z < 0, z * lo_hw, z * hi_hw) / 1.96)
        notes.append(f"the hard-edge (non-paralyzable) model is rejected: likelihood ratio {lr:.1f}")
    else:
        upper = math.hypot(high, 1.96 * shift_sd) if shift_sd else high
        est_p = Interval(tau, 0.0, upper)
        tau_draws = rng.uniform(0.0, upper, n_draws)
        unidentified = ("dead_time_p_ps",)
        stand_in = ("dead_time_p_ps",)

    # tau_np from the shortest gap: the summed gap density at the edge, over the rates' uncertainty.
    def edge_density(n_vals: np.ndarray, t: np.ndarray) -> np.ndarray:
        e_t = np.exp(-n_vals * t)
        out = n_vals * e_t
        if p > 0 and tau_a:
            out = out + p * e_t / tau_a
        return out

    sizes = np.array([r.n_gaps for r in binned], dtype=float)
    n_arr = np.array(rates)
    f0 = float(np.sum(sizes * edge_density(n_arr, np.full_like(n_arr, tau))))
    srng = np.random.default_rng(0)
    rel = 1.0 / np.sqrt(sizes)
    n_s = n_arr[None, :] * np.exp(srng.standard_normal((F0_SAMPLES, n_arr.size)) * rel[None, :])
    t_s = np.resize(tau_draws, F0_SAMPLES)[:, None]
    f0_samples = np.sum(sizes[None, :] * edge_density(n_s, t_s), axis=1)
    est_np, np_draws = _dead_time(m, 1, f0, f0_samples, rng, n_draws, tick, afterpulsing_ruled_out=True)
    # tau_p <= tau_np: the point, the upper end and every draw (reviews S6, D8: the profile runs up
    # to the shortest gap, which can lie above tau_np's point and upper end).
    cap_high = min(est_p.high, est_np.high)
    est_p = Interval(min(est_p.value, est_np.value), min(est_p.low, cap_high), cap_high)
    tau_draws = np.minimum(tau_draws, np_draws)
    if not identified:
        notes.append(
            f"the hard-edge (non-paralyzable) model is not rejected (likelihood ratio {lr:.2f}): dead_time_p_ps is "
            f"an upper limit ({est_p.high:.4g} ps); draws are uniform below it, a stand-in"
        )
    if high >= tau_hi:
        notes.append("tau_p reaches tau_np: the data do not separate the stages (consistent with a fully "
                     "paralyzable dead time)")

    probs = [r.probs(n, tau, p, tau_a) for r, n in zip(binned, rates)]
    ppc = with_interruption_check_runs(_pearson(binned, probs), notes, tags)
    for i, n in enumerate(rates):
        notes.append(f"run {i}: input rate {n * PS_PER_S:.4g} cps (fitted), {binned[i].n_gaps} gaps")
    cond = {"runs": float(len(runs)), "lr_hard_edge": lr, "gaps": float(sizes.sum())}
    if p > 0:
        cond.update({"afterpulse_prob_held": p, "afterpulse_tau_ps_held": float(tau_a)})
        pending = max(p * n * tau_a / (1.0 + n * tau_a) for n in rates)
        cond["afterpulse_pending_max"] = pending
        if pending > PENDING_MAX:
            notes.append(f"about {pending:.1%} of gaps start with an afterpulse still pending, which the gap law "
                         "does not model: tau_p is biased low (by about 1 sd at 3%); scan at lower rates or with a "
                         "detector whose afterpulsing is faster")
    estimates = {"dead_time_p_ps": est_p, "dead_time_ps": est_np}
    return Calibration(
        kind=KIND,
        method=METHOD,
        method_version=METHOD_VERSION,
        channels=(channel,),
        estimates=estimates,
        units=dict(UNITS),
        definitions=dict(DEFINITIONS),
        samples=np.column_stack([tau_draws, np_draws]),
        sample_params=("dead_time_p_ps", "dead_time_ps"),
        converged=True,
        ppc_pvalue=ppc,
        inputs=tuple((tags_sha256(t), "raw-tags", str(channel)) for t in tags),
        backend=backend,
        conditions={**cond, **(conditions or {})},
        unidentified=unidentified,
        stand_in=stand_in,
        notes=tuple(notes),
    )
