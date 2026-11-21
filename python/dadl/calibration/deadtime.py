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
sensitivity, of tau_p or of its upper limit) and added to the intervals and draws.

An afterpulse still pending when the next gap starts (from the previous avalanche, Exp(tau_a)
after the record by memorylessness, a fraction pi of gaps recorded as `afterpulse_pending_max`) is
modelled to first order: one pending candidate, without its product with the record's own
afterpulse (council round 5). At v = r - tau_np in u, with r ~ Exp(tau_a) after the record, it has
no effect below -tau_p, restarts the paralyzable stage between -tau_p and 0, and acts as the
record's own afterpulse above 0:

    g(u) = g_ap(u) + pi Delta1(u),  K = e^{-tau_np/tau_a},  Delta_own = the bracket above,
    Delta1(u) = K Delta_own(u) - K (e^{tau_p/tau_a} - 1) g0(u) + (K/tau_a) int_0^tau_p e^{s/tau_a} q(u + s) ds
    pi_own = p int f_E(e) (1 - S(e + tau_p)) de,  pi = pi_own / (1 - K pi_own / p)

(Delta1 integrates to 0; pi counts a pending one that persists through the next gap.) The
remainder, about pi (p + pi) of gaps, biases tau_p low by an amount that does not shrink with data
(about 0.1 sd per unit of remainder x sqrt(gaps)); above 2 a note says so. Checked against the
simulator at 6e6 cps (p 2%, tau_a 30 and 90 ns) and 1e7 cps (p 5%, tau_a 200 ns; pi 6.8%, where
the remainder shows at 1.3e7 gaps, chi-square p 0.05).

**Two rates per gap law.** Without afterpulsing the law depends on n only through
a = n exp(-n tau_p), which takes each value twice (n tau_p below and above 1): past the count-rate
maximum the gaps do not tell which rate it is. With afterpulsing the branches differ a little.
Each run is taken below the count-rate maximum (n tau_p < 1, the designed scan) unless the fit
with each run on its better branch is better by more than `BRANCH_MARGIN` = 5 in NLL (a likelihood
ratio of about 150); then those runs are fitted past it. Within the margin the interval spans both
fits and a note names the runs; the draws are an equal-weight mixture of the two fits' split
normals (council round 6, B6: one split normal over the union put 18% of draws in the gap between). Below the maximum the high branch fits worse in expectation
(Asimov: 2.5-5.6 NLL for the CI and CHANGELOG scans) and wins by more than 3 in about 0.5% of scans
(9 of 2,000, round 5); with margin 5 the false-switch rate is at most about 8e-4 for any separation
(Cox's normal approximation for non-nested fits: Phi(-sqrt(2M))). Past it (1/3/60 Mcps) the margin
is passed in 28 of 48 and the union covers the rest. Refits for the afterpulse sensitivity hold
each run's branch. Past the maximum, coverage is checked by probes only (no CI regime; M10).

**Fit.** tau_np is the shortest gap over all runs, with the dark fit's exact tick-quantile
interval (`dadl.calibration.dark._dead_time`) on the summed gap density at the edge, widened by one
tick either side (council round 6, T1). Dead times are read on the tagger's tick: a recorded gap is
within one tick of the true gap, and a true tau_np off the tick grid (real dead times are analogue)
lies up to one tick above the shortest gap and, with the edge excess, below the quantile interval's
lower end. The interval therefore allows one tick either side, and the draws add the triangular
rounding of a difference of two floored times (tick x (u1 - u2)). The dark dead-time interval does
not; it spans many ticks in any dark run, so the difference is below its resolution. On the
simulator's grid (integer tau_np, 1 ps tags) the widening is pure conservatism. The tick is the
larger of `tick_ps`, when given, and the common divisor of the tags' offsets: on PTU T3 data the
tags are nsync x period + dtime x resolution and a measured (non-whole) period leaves no common
divisor, while each recorded gap is still within one dtime resolution of the true gap; the CLI
passes the file's tick (`dadl.tagfile.tag_tick_ps`; council round 7, G1: at 25 ps the 1 ps
widening covered 10 of 20). `tag_tick_ps` records the tick used, `tag_tick_declared` (1 or 0) and a note its source. The gaps
past it are binned (fine uniform bins, then geometric ones) and fitted by multinomial maximum
likelihood in (tau_p, n_1 ... n_K); given tau_p the runs are independent, so the profile
likelihood of tau_p is a sum of one-dimensional fits. The 95% interval for tau_p is the profile
likelihood's (chi-square(1)) on the decided branches (above); when it reaches tau_p = 0 the
hard-edge model is not rejected and tau_p is an upper limit (unidentified; draws uniform below
it, a stand-in). tau_p's point, upper end and draws are capped at tau_np's. When the profile
reaches the shortest gap, tau_p's upper end is raised to tau_np's (one tick above it; council round
10, Rg1). 2 ln(likelihood
ratio) against the hard-edge model is recorded (`two_ln_lr_hard_edge`).

**Intervals.** The 95% interval of an identified tau_p combines the profile half-widths
(statistical, recorded as `dead_time_p_stat_low_ps` and `dead_time_p_stat_high_ps`) with 1.96
times the afterpulse shift sd (`dead_time_p_u_afterpulse_ps`, from the sensitivities
`dtau_p_dp_ps` and `dtau_p_dtau_a` and the dark draws' covariance, correlation
`afterpulse_p_tau_corr_used`), side by side in quadrature. An upper limit on tau_p uses the
two-sided 95% threshold (2 ln LR = 3.84); on the tau_p >= 0 boundary that is a one-sided limit of
about 97.5%. No coverage regime is unidentified, so this level is not gated. The scanned range is
recorded (`scan_rate_low_cps`, `scan_rate_high_cps`: recorded count rates; `scan_duration_low_s`,
`scan_duration_high_s`); an effective tau_p describes those rates, not rates above them (checked at
use from M5).

**Conditions.** The dark artifact's recorded conditions are compared with the scan's
(`dadl.propagate.condition_mismatches`): afterpulsing depends on bias and temperature, so a
mismatch is refused (`StaleCalibration`) unless `allow_condition_mismatch` gives a reason (then
noted). The dark artifact's bias and temperature are copied as `dark_bias_V` and
`dark_temperature_K`.

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

import functools
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
METHOD_VERSION = "8"
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
#: Grid points of the tau_p profile before its local minima are refined.
PROFILE_GRID = 30
F0_SAMPLES = 4000
#: A note when the unmodelled pending remainder r (about pi (p + pi) of gaps) has r sqrt(gaps) above this.
UNMODELLED_Z = 2.0
#: NLL by which the fit with runs past the count-rate maximum must beat the one below it (2 ln LR = 10).
BRANCH_MARGIN = 5.0
PS_PER_S = 1e12

UNITS = {"dead_time_p_ps": "ps", "dead_time_ps": "ps"}
DEFINITIONS = {
    "dead_time_p_ps": "paralyzable stage: an avalanche candidate passes only if the previous candidate (passed or not) "
    "is at least this much earlier",
    "dead_time_ps": "non-paralyzable dead time as recorded: the lower edge of the gap distribution, at the tagger's "
    "tick, so it includes the tagger's timing jitter; estimated as the shortest gap minus ln 2 / (edge density), to "
    "within one tick (the interval allows one tick either side of the tick-quantile interval)",
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


def _series(x: np.ndarray, a: float, tau: float) -> np.ndarray:
    out = np.zeros_like(x)
    top = float(x.max(initial=-1.0))
    for j in range(K_TERMS + 1):
        if j * tau > top and j > 0:
            break  # every remaining term is masked to 0.0: adding it changes nothing
        r = x - j * tau
        out += np.where(r >= 0.0, (-a * np.maximum(r, 0.0)) ** j / math.factorial(j), 0.0)
    return out


@functools.lru_cache(maxsize=4096)
def _survival_consts(n: float, tau: float):
    a = n * math.exp(-n * tau)
    u0 = K_TERMS * tau
    # Past u0 the two real modes, matched to S(u0) and S'(u0) = -a S(u0 - tau). The divided
    # difference is written with expm1 so it stays exact as the modes merge (n tau near 1).
    s0 = float(_series(np.array([u0]), a, tau)[0])
    ds0 = -a * float(_series(np.array([u0 - tau]), a, tau)[0])
    l1, l2 = _modes(a, tau)
    return a, u0, s0, ds0, l1, l2


def survival(w: np.ndarray, n: float, tau: float) -> np.ndarray:
    """S(w): 1 below tau, S' = -a S(w - tau) above (see the module docstring)."""
    w = np.asarray(w, dtype=float)
    n, tau = float(n), float(tau)
    if tau <= 0.0:
        return np.exp(-n * np.maximum(w, 0.0))
    a, u0, s0, ds0, l1, l2 = _survival_consts(n, tau)

    def series(x: np.ndarray) -> np.ndarray:
        return _series(x, a, tau)
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


def _pending(p: float, x: float, k: float) -> float:
    """pi = pi_own / (1 - K pi_own / p), with x = pi_own / p and K = exp(-tau_np / tau_a)."""
    return p * x / (1.0 - k * x)


def pending_fraction(n: float, tau_p: float, tau_np: float, p: float, tau_a: float) -> float:
    """The fraction of gaps that start with an afterpulse candidate still pending (module
    docstring, "Afterpulsing"), on a fine grid of its own."""
    if p <= 0 or not tau_a:
        return 0.0
    cells = 8000
    h = AFTERPULSE_TAUS * tau_a / cells
    e = (np.arange(cells) + 0.5) * h
    x = float(np.sum(np.exp(-e / tau_a) / tau_a * h * (1.0 - survival(e + tau_p, n, tau_p))))
    x += math.exp(-AFTERPULSE_TAUS)
    return _pending(p, x, math.exp(-tau_np / tau_a))


def gap_cdf(u: np.ndarray, n: float, tau: float) -> np.ndarray:
    """P(gap - tau_np <= u) without afterpulsing."""
    n, tau = float(n), float(tau)
    return 1.0 - survival(np.asarray(u, dtype=float) + tau, n, tau)


class _Run:
    """One acquisition's gaps (past the shortest), binned for the fit."""

    def __init__(self, s: np.ndarray, tick: int, tau_a: float | None, tau_np: float | None = None) -> None:
        #: tau_np (the fit's shortest gap): with it, the pending-afterpulse term is added (p > 0).
        self.tau_np = tau_np
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
            own, delta1, x = self._afterpulse(n, tau, tau_a)
            prob = prob + p * own * self.width
            if delta1 is not None:
                prob = prob + _pending(p, x, math.exp(-self.tau_np / tau_a)) * delta1 * self.width
        return np.maximum(prob, 1e-300)

    def _afterpulse(self, n: float, tau: float, tau_a: float) -> tuple[np.ndarray, np.ndarray | None, float]:
        """At the bin midpoints: Delta_own(u) (module docstring: the bracket multiplying p), and with
        `tau_np` the pending term Delta1(u) and pi_own / p (else None and 0)."""
        u = self.mid
        e_t = math.exp(-n * tau)
        a = n * e_t
        s_u = survival(u, n, tau)
        decay = np.exp(-u / tau_a)
        out = -a * s_u * (1.0 - decay) + decay * s_u * e_t / tau_a
        h = self.h
        n_w = int(min(self.n_uni, math.ceil(AFTERPULSE_TAUS * tau_a / h)))
        w = (np.arange(n_w) + 0.5) * h
        s_wt = survival(w + tau, n, tau)
        f_e = np.exp(-w / tau_a) / tau_a * h
        f = f_e * (s_wt - survival(w, n, tau) * e_t)
        pending = self.tau_np is not None
        # The pending term reads q up to tau_p beyond the grid: ceil(tau_p / h) more lags.
        n_s = int(math.ceil(tau / h)) if pending and tau > 0 else 0
        lags = np.arange(self.n_uni + n_s + 1) * h

        def q_cell(t: np.ndarray) -> np.ndarray:
            # q averaged over the cell [t - h/2, t + h/2]: a S(t' - tau) for t' >= tau integrates to
            # S(lo) - S(hi) (S = 1 below tau), exact and smooth in tau. The point value snapped the
            # jump at tau to the grid and made the profile jagged (council round 1).
            lo, hi = np.maximum(t - 0.5 * h, 0.0), np.maximum(t + 0.5 * h, 0.0)
            return (survival(lo, n, tau) - survival(hi, n, tau)) / h

        q_lags = q_cell(lags)
        conv = signal.fftconvolve(f, q_lags[: self.n_uni])[: self.n_uni]
        out[: self.n_uni] += conv
        geo = u.size > self.n_uni
        if geo:
            t = u[self.n_uni :, None] - w[None, :]
            out[self.n_uni :] += q_cell(t.ravel()).reshape(t.shape) @ f
        if not pending:
            return out, None, 0.0
        # pi_own / p = int f_E(e) (1 - S(e + tau_p)) de; past the grid S is negligible (the grid
        # reaches 40 tau_a, the longest gap or 30 mean gaps), so the tail is exp(-W / tau_a).
        x = float(np.sum(f_e * (1.0 - s_wt))) + math.exp(-n_w * h / tau_a)
        k = math.exp(-self.tau_np / tau_a)
        delta1 = k * out - k * math.expm1(tau / tau_a) * a * s_u
        if n_s:
            # (K / tau_a) int_0^tau_p e^{s/tau_a} q(u + s) ds: cells of width h in s (the last one cut
            # at tau_p), q at u + s on the lag grid; a correlation with the lags already computed.
            s_edges = np.minimum(np.arange(n_s + 1) * h, tau)
            wgt = k * np.diff(np.exp(s_edges / tau_a))
            delta1[: self.n_uni] += signal.fftconvolve(q_lags[1:], wgt[::-1], mode="valid")[: self.n_uni]
            if geo:
                t = u[self.n_uni :, None] + ((np.arange(n_s) + 0.5) * h)[None, :]
                delta1[self.n_uni :] += q_cell(t.ravel()).reshape(t.shape) @ wgt
        return out, delta1, x

    def nll(self, n: float, tau: float, p: float = 0.0, tau_a: float | None = None) -> float:
        prob = self.probs(n, tau, p, tau_a)
        return float(-np.sum(self.counts[self.nz] * np.log(prob[self.nz])))


def _other_branch(n: float, tau: float) -> float:
    """The other input rate with the same a = n exp(-n tau) (Lambert W branch -1 below n tau = 1,
    branch 0 above)."""
    x = n * tau
    w = special.lambertw(-x * math.exp(-x), -1 if x < 1.0 else 0)
    return float(-w.real / tau)


class _Fit:
    def __init__(self, runs: list[_Run], p: float, tau_a: float | None) -> None:
        self.runs, self.p, self.tau_a = runs, p, tau_a
        self._n_cache: list[float] = [1.0 / r.mean for r in runs]
        #: Per run, the last call's other-branch (rate, NLL above the kept one), or None.
        self.other: list[tuple[float, float] | None] = [None] * len(runs)
        #: False: the low branch only (a warm start there).
        self.both_branches = True

    def best_n(self, run_index: int, tau: float) -> tuple[float, float]:
        """The run's input rate maximizing the likelihood at this tau_p: Newton from the cached
        low-branch rate and, when tau > 0, from the other branch's rate (`_other_branch`); the
        lower NLL is kept. The low-branch rate is cached (a warm start keeps the branch); a
        high-branch rate moves fast with tau, so it always starts from the low one. With
        `both_branches` False (the decided fit and the afterpulse refits) only the cached branch is
        searched; both are searched for the branch decision and the p = 0 rate notes."""
        n, f0 = self._newton(run_index, tau, self._n_cache[run_index])
        self.other[run_index] = None
        if tau <= 0.0 or not self.both_branches:
            self._n_cache[run_index] = n
            return n, f0
        n2, f2 = self._newton(run_index, tau, self._branch_start(run_index, n, tau))
        self._n_cache[run_index] = min(n, n2)
        if f2 < f0:
            n, f0, n2, f2 = n2, f2, n, f0
        if abs(math.log(n2 / n)) > 1e-3:
            self.other[run_index] = (n2, f2 - f0)
        return n, f0

    def _branch_start(self, run_index: int, n: float, tau: float) -> float:
        """The other branch's rate, or the recorded rate where n tau e^{-n tau} underflows."""
        other = _other_branch(n, tau)
        return other if math.isfinite(other) and other > 0 else 1.0 / self.runs[run_index].mean

    def _newton(self, run_index: int, tau: float, n0: float) -> tuple[float, float]:
        """Newton steps in log n from n0 (the likelihood is close to quadratic there: its curvature
        is about the number of gaps)."""
        run = self.runs[run_index]

        def f(y: float) -> float:
            return run.nll(math.exp(y), tau, self.p, self.tau_a)

        y = math.log(n0)
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
        return math.exp(y), f0

    def profile(self, tau: float) -> float:
        return sum(self.best_n(i, tau)[1] for i in range(len(self.runs)))

    def rates(self, tau: float) -> list[float]:
        return [self.best_n(i, tau)[0] for i in range(len(self.runs))]


def _minimize(fit: _Fit, grid: np.ndarray, vals: list[float], k: int) -> tuple[float, float]:
    """The profile's minimum near grid point k."""
    lo = math.log(grid[max(k - 1, 0)])
    hi = math.log(grid[min(k + 1, grid.size - 1)])
    res = optimize.minimize_scalar(lambda y: fit.profile(math.exp(y)), bounds=(lo, hi), method="bounded",
                                   options={"xatol": 1e-5})
    return (math.exp(float(res.x)), float(res.fun)) if res.fun < vals[k] else (float(grid[k]), vals[k])


def _fit_tau(fit: _Fit, tau_hi: float, held_min: tuple[float, float] | None = None,
             nll0: float | None = None) -> tuple[float, float, float, list[tuple[float, float]]]:
    """(tau_p MLE, its NLL, the NLL at tau_p = 0, the profile at grid points and minima).

    The profile can have two minima: a run whose rate may lie on either branch can fit about as
    well at a lower tau_p on the high branch (passive-quench, 6e6 cps: 30.8 ns there, 40 ns on the
    low branch, NLLs within 0.6). So the minima are sought twice, each refined from a grid: on the
    profile itself (every local minimum of a `PROFILE_GRID`-point grid) and with every run on the
    low branch (unimodal; its minimum is then refined on the profile). The lowest is the MLE; the
    interval spans every minimum inside it (`_profile_interval`). A fit held to the low branch
    (`both_branches` False) has the one minimum: that is the decided fit below the count-rate
    maximum. The both-branch search runs only to decide the branch (p > 0); `held_min` and `nll0`
    pass in the held fit's minimum and hard-edge NLL so they are not searched again."""
    nll0 = fit.profile(0.0) if nll0 is None else nll0
    both = fit.both_branches
    if held_min is not None and both:
        tau_low, low_best = held_min
    else:
        fit.both_branches = False
        coarse = np.geomspace(max(tau_hi * 1e-3, 1.0), tau_hi, 12)
        low_vals = [fit.profile(float(t)) for t in coarse]
        tau_low, low_best = _minimize(fit, coarse, low_vals, int(np.argmin(low_vals)))
        fit.both_branches = both
    if not both:
        points = [(float(t), v) for t, v in zip(coarse, low_vals)] + [(tau_low, low_best)]
        if low_best >= nll0:
            return 0.0, nll0, nll0, points
        fit.profile(tau_low)
        return tau_low, low_best, nll0, points
    grid = np.geomspace(max(tau_hi * 1e-3, 1.0), tau_hi, PROFILE_GRID)
    vals = [fit.profile(float(t)) for t in grid]
    points = [(float(t), v) for t, v in zip(grid, vals)]
    minima = [_minimize(fit, grid, vals, k) for k in range(grid.size)
              if vals[k] < nll0 and vals[k] <= min(vals[max(k - 1, 0)], vals[min(k + 1, grid.size - 1)])]
    near = np.array([tau_low * 0.9, tau_low, min(tau_low * 1.1, tau_hi)])
    minima.append(_minimize(fit, near, [fit.profile(float(t)) for t in near], 1))
    points += minima
    best_tau, best = min(minima, key=lambda m: m[1])
    if best >= nll0:
        return 0.0, nll0, nll0, points
    fit.profile(best_tau)  # leave the warm starts at the MLE
    return best_tau, best, nll0, points


def _refit_near(fit: _Fit, tau: float, tau_hi: float) -> float:
    """tau_p MLE near a previous one (the afterpulse sensitivity refits)."""
    lo, hi = math.log(max(tau * 0.7, 1.0)), math.log(min(tau * 1.3, tau_hi))
    res = optimize.minimize_scalar(lambda y: fit.profile(math.exp(y)), bounds=(lo, hi), method="bounded",
                                   options={"xatol": 2e-4})
    return math.exp(float(res.x))


def _profile_interval(
    fit: _Fit, tau: float, nll_min: float, nll0: float, tau_hi: float, points: list[tuple[float, float]]
) -> tuple[float, float]:
    """The 95% profile-likelihood interval for tau_p: from the lowest to the highest tau_p whose
    profile is within chi-square(1) of the minimum (`points`, the grid from `_fit_tau`, locate
    them; both minima are inside when the profile has two)."""
    target = nll_min + PROFILE_TARGET / 2.0

    def excess(t: float) -> float:
        return fit.profile(t) - target

    def root(t_out: float, t_in: float) -> float:
        a, b = sorted((t_out, t_in))
        if excess(t_out) <= 0 or excess(t_in) > 0:
            return t_in
        return math.exp(optimize.brentq(lambda y: excess(math.exp(y)), math.log(a), math.log(b), xtol=1e-4))

    pts = sorted(points + [(max(tau, 1.0), nll_min)])
    inside = [t for t, v in pts if v <= target]
    lo_in, hi_in = inside[0], inside[-1]
    if tau <= 0.0 or nll0 <= target:
        low = 0.0
    else:
        below = [t for t, v in pts if t < lo_in and v > target]
        t_out = below[-1] if below else max(lo_in * 1e-6, 1e-3)
        low = root(t_out, lo_in) if excess(t_out) > 0 else 0.0
    above = [t for t, v in pts if t > hi_in and v > target]
    if above:
        high = root(above[0], hi_in)
    else:
        high = root(tau_hi, hi_in) if excess(tau_hi) > 0 else tau_hi
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
    allow_condition_mismatch: str | None = None,
    tick_ps: int | None = None,
) -> Calibration:
    """Fit tau_p and tau_np from illuminated acquisitions of one channel at two or more rates.

    `runs` are sorted int64 ps tags, one array per acquisition, under steady light at
    different rates. `afterpulse` is the channel's dark artifact (`calibrate_dark`), whose
    afterpulse probability and time constant are held fixed with their uncertainty propagated.

    Refuses fewer than two runs, a run under 1,000 gaps, duplicate timestamps, the same run twice,
    an afterpulse artifact of another kind or channel, and one whose recorded conditions (bias,
    temperature, tagger settings) do not cover `conditions` (`StaleCalibration`), unless
    `allow_condition_mismatch` gives a reason, which is noted. Each run is fitted below the count-rate
    maximum unless the fit past it wins by `BRANCH_MARGIN` (module docstring, "Two rates per gap
    law"). `tick_ps` is the tagger's tick when the tags do not show it (PTU T3: the dtime resolution).

    On a coarse tag tick (PTU T3, about 25-80 ps), a dead time close to fully paralyzable has tau_p's
    upper end at tau_np's upper end, one tick above the shortest recorded gap, since either stage can
    lie up to a tick above it. The interval is then about one tick wider than the profile alone
    (council round 10, Rg1: coverage 60 / 60 against 7-8 / 60 before).
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
    for i, g in enumerate(gaps):
        if g.size and int(g.min()) <= 0:
            raise ValueError(f"duplicate timestamps in runs[{i}]: a detector channel records one event per tag")
    digests = [tags_sha256(t) for t in tags]
    for i, d in enumerate(digests):
        if d in digests[:i]:
            raise ValueError(f"runs[{i}] is the same acquisition as runs[{digests.index(d)}]: each run needs its own "
                             "rate")
    if afterpulse is not None:
        if afterpulse.kind != "dark-afterpulse-deadtime":
            raise ValueError(f"afterpulse must be a dark artifact (kind dark-afterpulse-deadtime), not "
                             f"{afterpulse.kind}")
        if channel not in afterpulse.channels:
            have = ", ".join(str(c) for c in afterpulse.channels)
            raise ValueError(f"the afterpulse artifact is for channel(s) {have}, not {channel}")
    notes: list[str] = []
    dark_conditions: dict[str, float] = {}
    if afterpulse is not None:
        # Council round 6 (T2): a 25 V / 233 K dark artifact was held for a 30 V / 300 K scan, silently.
        from dadl.errors import StaleCalibration
        from dadl.propagate import condition_mismatches

        dark_recorded = dict(getattr(afterpulse, "conditions", None) or {})
        mism = condition_mismatches(dark_recorded, dict(conditions or {}))
        if mism and not (allow_condition_mismatch and allow_condition_mismatch.strip()):
            raise StaleCalibration(
                f"the dark artifact's conditions do not match this scan's: {'; '.join(mism)}. Afterpulsing depends "
                "on bias and temperature; calibrate dark under the scan's conditions, or pass "
                "allow_condition_mismatch=REASON (CLI: --allow-condition-mismatch REASON)."
            )
        if mism:
            notes.append(f"the dark artifact's conditions differ from this scan's ({'; '.join(mism)}); used anyway: "
                         f"{allow_condition_mismatch.strip()}")
        # Names outside the use-time prefixes (bias*, temperature_K), so analyses do not demand them.
        for key, value in dark_recorded.items():
            if key.startswith("bias"):
                dark_conditions[f"dark_{key}"] = float(value)
            elif key == "temperature_K":
                dark_conditions["dark_temperature_K"] = float(value)
    if afterpulse is not None and afterpulse.status == "misfit":
        notes.append(f"the afterpulse artifact for {channel} is a misfit (predictive p = {afterpulse.ppc_pvalue:.2g}): "
                     "its p and tau_a are used as given")
    short = [i for i, g in enumerate(gaps) if g.size < MIN_GAPS]
    if short:
        raise InsufficientCounts(f"need at least {MIN_GAPS} gaps per acquisition; run(s) {short} have fewer")
    tick = tag_tick(*tags)
    tick_source = "tags"
    if tick_ps is not None:
        declared = _checks.integer(tick_ps, "tick_ps", minimum=1)
        if declared > tick:
            tick, tick_source = declared, "declared"
    if tick_source == "declared":
        notes.append(f"tag tick {tick} ps as declared (tick_ps), not the {tag_tick(*tags)} ps the tags show")
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
    binned = [_Run(s, tick, tau_a, float(m)) for s in shifted]
    # The branch of each run is decided once (council round 4): the low branch (n tau_p < 1, the
    # designed scan) unless the fit with each run on its better branch is better by more than
    # BRANCH_MARGIN in NLL; within the margin, the interval spans both fits.
    fit = _Fit(binned, p, tau_a)
    fit.both_branches = False
    tau_hi = float(m)  # tau_p <= tau_np
    tau, nll_min, nll0, points = _fit_tau(fit, tau_hi)
    lr = 2.0 * (nll0 - nll_min)
    low, high = _profile_interval(fit, tau, nll_min, nll0, tau_hi, points)
    rates = fit.rates(tau)
    branch_gain = 0.0
    #: Within the branch margin: the two fits (MLE, low, high), for the mixture draws.
    union: list[tuple[float, float, float]] = []
    if p > 0:
        glob = _Fit(binned, p, tau_a)
        t_g, nll_g, nll0_g, pts_g = _fit_tau(glob, tau_hi, held_min=(tau if tau > 0 else 1.0, nll_min), nll0=nll0)
        rates_g = glob.rates(t_g)
        past = [i for i, n in enumerate(rates_g) if n * t_g > 1.0]
        branch_gain = max(nll_min - nll_g, 0.0) if past else 0.0
        if past:
            low_g, high_g = _profile_interval(glob, t_g, nll_g, nll0_g, tau_hi, pts_g)
            if branch_gain > BRANCH_MARGIN:
                notes.append(f"run(s) {past} fitted past the count-rate maximum (n tau_p > 1): that fit is better "
                             f"by {branch_gain:.1f} in NLL than every run below it")
                tau, nll_min, nll0, low, high, rates = t_g, nll_g, nll0_g, low_g, high_g, rates_g
                lr = 2.0 * (nll0 - nll_min)
            else:
                notes.append(f"run(s) {past} fit about as well past the count-rate maximum ({branch_gain:.1f} in "
                             f"NLL, margin {BRANCH_MARGIN:g}): the interval spans that fit ({low_g:.4g} to "
                             f"{high_g:.4g} ps) and the one below it ({low:.4g} to {high:.4g} ps)")
                union = [(tau, low, high), (t_g, low_g, high_g)]
                low, high = min(low, low_g), max(high, high_g)
    identified = low > 0.0

    # Afterpulse uncertainty: d tau_p / d(p, tau_a) by refitting at shifted values; for an upper
    # limit, the sensitivity of the limit itself (whatever the MLE).
    shift_sd = 0.0
    sens = np.zeros(2)
    ap_corr = 0.0
    if ap_draws is not None:
        sd = ap_draws.std(axis=0)
        for k in range(2):
            if sd[k] <= 0:
                continue
            step = sd[k]
            p_k, ta_k = (p + step, tau_a) if k == 0 else (p, tau_a + step)
            other = _Fit(binned, p_k, ta_k)
            other.both_branches = False  # each run held on its decided branch (council round 4, DT1)
            other._n_cache = list(rates)
            if identified:
                sens[k] = (_refit_near(other, tau, tau_hi) - tau) / step
            else:
                t_k, nll_k, nll0_k, pts_k = _fit_tau(other, tau_hi)
                sens[k] = (_profile_interval(other, t_k, nll_k, nll0_k, tau_hi, pts_k)[1] - high) / step
        cov = np.cov(ap_draws.T)
        shift_sd = float(math.sqrt(max(sens @ cov @ sens, 0.0)))
        sd_p, sd_t = math.sqrt(max(cov[0, 0], 0.0)), math.sqrt(max(cov[1, 1], 0.0))
        ap_corr = float(cov[0, 1] / (sd_p * sd_t)) if sd_p > 0 and sd_t > 0 else 0.0
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
        if union:
            # Two fits within the branch margin: an equal-weight mixture of their split normals
            # (conservative), not one split normal over the union (council round 6, B6).
            second = rng.random(n_draws) < 0.5
            for k, (t_c, lo_c, hi_c) in enumerate(union):
                lo_w, hi_w = math.hypot(t_c - lo_c, 1.96 * shift_sd), math.hypot(hi_c - t_c, 1.96 * shift_sd)
                pick = second == bool(k)
                tau_draws[pick] = np.maximum(0.0, t_c + np.where(z[pick] < 0, z[pick] * lo_w, z[pick] * hi_w) / 1.96)
        notes.append(f"the hard-edge (non-paralyzable) model is rejected: 2 ln(likelihood ratio) = {lr:.1f}")
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
    # One tick either side (council round 6, T1): a recorded gap is within one tick of the true gap,
    # and an off-grid true tau_np was missed every time by the zero-width interval [m, m].
    np_low, np_high = max(0.0, est_np.low - tick), est_np.high + tick
    est_np = Interval(min(max(est_np.value, np_low), np_high), np_low, np_high)
    np_draws = np.maximum(0.0, np_draws + tick * (rng.random(n_draws) - rng.random(n_draws)))
    stat_low, stat_high = (low, high) if identified else (0.0, high)
    # tau_p <= tau_np: the point, the upper end and every draw (reviews S6, D8: the profile runs up
    # to the shortest gap, which can lie above tau_np's point and upper end).
    if high >= tau_hi:  # the profile reached the shortest gap: tau_p, like tau_np, may sit a tick above it (r10 Rg1)
        est_p = Interval(est_p.value, est_p.low, max(est_p.high, est_np.high))
        stat_high = max(stat_high, est_np.high)
    cap_high = min(est_p.high, est_np.high)
    est_p = Interval(min(est_p.value, est_np.value), min(est_p.low, cap_high), cap_high)
    stat_high = min(stat_high, est_np.high)
    stat_low = min(stat_low, stat_high)
    tau_draws = np.minimum(tau_draws, np_draws)
    if not identified:
        notes.append(
            f"the hard-edge (non-paralyzable) model is not rejected (2 ln(likelihood ratio) = {lr:.2f}): "
            f"dead_time_p_ps is an upper limit ({est_p.high:.4g} ps) at the two-sided 95% threshold (2 ln LR = "
            "3.84), on the tau_p >= 0 boundary a one-sided limit of about 97.5%; draws are uniform below it, a "
            "stand-in"
        )
    if high >= tau_hi:
        notes.append("tau_p reaches tau_np: the data do not separate the stages (consistent with a fully "
                     "paralyzable dead time)")

    probs = [r.probs(n, tau, p, tau_a) for r, n in zip(binned, rates)]
    ppc = with_interruption_check_runs(_pearson(binned, probs), notes, tags)
    if p == 0.0:
        fit.both_branches = True
        fit.rates(tau)  # each run's other branch, for the notes
    for i, n in enumerate(rates):
        other = fit.other[i]
        if p == 0.0 and other is not None and other[1] < 0.5:
            lo_n, hi_n = sorted((n, other[0]))
            notes.append(f"run {i}: avalanche-candidate rate (efficiency x photon rate + dark) {lo_n * PS_PER_S:.4g} "
                         f"or {hi_n * PS_PER_S:.4g} s^-1 (fitted; past the count-rate maximum the gaps do not tell "
                         f"which), {binned[i].n_gaps} gaps")
        else:
            notes.append(f"run {i}: avalanche-candidate rate (efficiency x photon rate + dark) {n * PS_PER_S:.4g} s^-1 "
                         f"(fitted), {binned[i].n_gaps} gaps")
    spans = np.array([float(t[-1] - t[0]) / PS_PER_S for t in tags])
    recorded = np.array([g.size for g in gaps]) / spans
    cond = {"runs": float(len(runs)), "two_ln_lr_hard_edge": lr, "gaps": float(sizes.sum()),
            "tag_tick_ps": float(tick), "tag_tick_declared": float(tick_source == "declared"),
            "dead_time_p_stat_low_ps": float(stat_low), "dead_time_p_stat_high_ps": float(stat_high),
            "dead_time_p_u_afterpulse_ps": shift_sd,
            # Recorded count rates and run durations of the scan (council round 6, T6); named outside
            # the _min_/_max_ bracket convention until M5 checks them at use.
            "scan_rate_low_cps": float(recorded.min()), "scan_rate_high_cps": float(recorded.max()),
            "scan_duration_low_s": float(spans.min()), "scan_duration_high_s": float(spans.max()),
            **dark_conditions}
    if p > 0:
        cond.update({"afterpulse_prob_held": p, "afterpulse_tau_ps_held": float(tau_a),
                     "branch_nll_gain": branch_gain})
        if ap_draws is not None:
            cond["afterpulse_prob_sd_held"] = float(ap_draws[:, 0].std())
            cond["afterpulse_tau_ps_sd_held"] = float(ap_draws[:, 1].std())
            cond["dtau_p_dp_ps"] = float(sens[0])
            cond["dtau_p_dtau_a"] = float(sens[1])
            cond["afterpulse_p_tau_corr_used"] = ap_corr
        pis = [pending_fraction(n, tau, m, p, tau_a) for n in rates]
        r = max(pi * (p + pi) for pi in pis)
        cond["afterpulse_pending_max"] = max(pis)
        cond["afterpulse_unmodelled"] = r
        z = r * math.sqrt(float(sizes.sum()))
        if z > UNMODELLED_Z:
            notes.append(f"about {r:.2%} of gaps carry a second pending afterpulse or one beside the record's own, "
                         "which the gap law does not model: tau_p is biased low by an amount that does not shrink "
                         f"with more data (about {0.1 * z:.1f} sd here, twice that for tau_a beyond twice tau_p); "
                         "scan at lower rates or with faster afterpulsing")
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
        # The dark artifact's own inputs too (as dcr does): two dark artifacts give different inputs.
        inputs=tuple(dict.fromkeys([*((d, "raw-tags", str(channel)) for d in digests),
                                    *(tuple(i) for i in (afterpulse.inputs if afterpulse is not None else ()))])),
        backend=backend,
        conditions={**cond, **(conditions or {})},
        unidentified=unidentified,
        stand_in=stand_in,
        notes=tuple(notes),
    )
