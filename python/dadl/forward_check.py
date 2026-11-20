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
rate per analysis is at most about `MISFIT_ALPHA` however many pairs it holds). `MISFIT_ALPHA` is
4e-4, this check's share of the 1e-3 per-analysis refusal budget (council round 8; README,
"Refusals and their rates").

Binning: the histogram bins are whole multiples of the lattice the pair delays sit on
(`PairData.lattice_ps` from the file's metadata, else the tags' common step), with their edges
half a lattice step off the modal residue of the delays, so a +-1 ps rounding never crosses an
edge (council round 8: on T3 tags with a measured sync period the tags' step is 1 ps and 46 ps
bins on a 25 ps lattice gave p about 1e-100). T3 input whose dtime resolution is not a whole
number of ps has no integer lattice: the check is skipped and the counts say so.

The check runs only where the model uses the peak shape: when a delay-jitter calibration of the
pair gives the window capture. A flagged analysis is refused (`Misspecified`) unless the policy allows `misfit_forward_model`;
then the corrected value carries `FORWARD_MISFIT_ACCEPTED`.

What it does not check: anything outside the window. For CHSH the calibrated delays are compared
with the analysed peaks through the window capture (`chsh`); for HOM and fringe a drift between
calibration and run is not seen (the delay and jitter are refitted here).
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
from dadl.guardrails import Policy, override_hint
from dadl.histogram import cross_correlation

#: Family-wise false-alarm rate of the check per analysis.
MISFIT_ALPHA = 4e-4
OVERRIDE = "misfit_forward_model"
FORWARD_MISFIT_ACCEPTED = "forward-model misfit accepted: the predictive check failed (policy modified)"
#: Recorded in the counts (`forward_check_note`) when the check is skipped for want of a lattice.
NO_LATTICE = "forward check skipped: no integer lattice (T3 dtime resolution not a whole number of ps)"
#: Fewest coincidences in the fit region for a histogram to be checked.
MIN_COUNTS = 50
#: Adjacent bins are merged until each expects at least this many counts.
MIN_EXPECTED = 5.0
#: Below this histogram p-value the fit is polished before the p-value is taken.
POLISH_BELOW = 0.01


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
    #: The lattice of the pair delays (ps), from the file's metadata; None takes the tags' step.
    lattice_ps: int | None = None
    #: The T3 sync period (ps), from the file's metadata: cross-sync delays are put back on the
    #: lattice (`dadl.calibration.delay.t3_group_histograms`; council round 9, Tw1).
    t3_period_ps: float | None = None


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

    def __init__(
        self, counts: np.ndarray, edges: np.ndarray, tick: float, offsets: Sequence[float],
        weights: Sequence[float] | None = None,
    ) -> None:
        self.n = counts.astype(float)
        self.lo, self.hi = edges[:-1], edges[1:]
        self.w = self.hi - self.lo
        self.tick = float(tick)
        self.offsets = tuple(float(o) for o in offsets)
        #: Mixture weights (T3 sync-count groups); None is the equal mixture.
        self.weights = None if weights is None else tuple(float(x) for x in weights)

    def mu(self, theta: np.ndarray) -> np.ndarray:
        log_bl, log_a, d, log_s, log_br = (float(v) for v in theta)
        s = math.exp(min(max(log_s, -50.0), 700.0))
        t = self.tick / s
        g = np.zeros_like(self.n)
        for i, o in enumerate(self.offsets):
            z_hi, z_lo = (self.hi - d - o) / s, (self.lo - d - o) / s
            part = (s / self.tick) * (_g(z_hi) - _g(z_hi - t) - _g(z_lo) + _g(z_lo - t))
            g += part if self.weights is None else self.weights[i] * part
        if self.weights is None:
            g /= len(self.offsets)
        right = np.clip(self.hi - d, 0.0, self.w)
        left = self.w - right
        mu = math.exp(min(log_bl, 700.0)) * left + math.exp(min(log_br, 700.0)) * right + math.exp(min(log_a, 700.0)) * g
        return np.maximum(mu, 1e-300)

    def nll(self, theta: np.ndarray) -> float:
        mu = self.mu(theta)
        return float(np.sum(mu - self.n * np.log(mu)))


def _tick(pairs: Sequence[PairData]) -> int:
    """The lattice the pairs' delays sit on: their declared `lattice_ps`, else the tags' step."""
    given = [int(p.lattice_ps) for p in pairs if p.lattice_ps]
    if given:
        return math.gcd(*given)
    return math.gcd(*(tag_tick(np.asarray(p.a, dtype=np.int64), np.asarray(p.b, dtype=np.int64)) for p in pairs))


def _histogram(pair: PairData, shift: int, width: int, half: int, tick: int) -> tuple[np.ndarray, dict[float, float]]:
    """A pair's histogram of b - a - shift, and the peak's T3 sync-count groups as {offset: count}:
    group k sits k r below the delay (council round 9, Tw1). Non-T3 input is one group at 0."""
    from dadl.calibration.delay import t3_group_histograms, t3_residue

    a = np.asarray(pair.a, dtype=np.int64)
    b = np.asarray(pair.b, dtype=np.int64)
    groups = t3_group_histograms(a, b, width, half, pair.t3_period_ps, tick, shift) if pair.t3_period_ps else None
    if not groups:
        hist = cross_correlation(a, b - shift, width, half).astype(float)
        return hist, {0.0: float(hist.sum())}
    r = t3_residue(float(pair.t3_period_ps), tick)
    return sum(groups.values()).astype(float), {-k * r: float(h.sum()) for k, h in groups.items()}


def _mixture(offsets: Sequence[float], groups: dict[float, float]) -> tuple[tuple[float, ...], tuple[float, ...] | None]:
    """The peak's offsets and weights: the declared offsets, each split over the T3 groups."""
    total = sum(groups.values())
    if set(groups) == {0.0} or total <= 0:
        return tuple(offsets), None
    pairs = [(o + g, n / total / len(offsets)) for o in offsets for g, n in groups.items()]
    return tuple(o for o, _ in pairs), tuple(w for _, w in pairs)


def _align(pairs: Sequence[PairData], shift: int, half: int, tick: int) -> tuple[int, int]:
    """(shift, tooth): the histogram shift moved so that the modal residue of the delays modulo
    the lattice sits `tooth` = tick // 2 ps above a bin edge (council round 8). A tooth at offset c
    in its bin reads as a delay c ps low in the model, so fitted delays get `tooth` added back.
    With a 1 ps lattice nothing moves."""
    if tick <= 1:
        return shift, 0
    fine = sum(_histogram(p, shift, 1, half, tick)[0] for p in pairs)
    x = -half + np.arange(fine.size)
    residue = np.bincount(np.mod(x, tick), weights=fine, minlength=tick)
    r0 = int(np.argmax(residue)) if residue.sum() > 0 else 0
    tooth = tick // 2
    return shift + (r0 - tooth) % tick, tooth


def _layout(pair: PairData, tick: int) -> tuple[int, int, int, float]:
    """(shift, half, bin width, starting jitter) of a pair's histogram at `tick`."""
    shift = int(round(pair.declared_delay_ps / tick)) * tick
    half = int(pair.window_ps)
    half = -(-half // tick) * tick
    sigma0 = pair.jitter_ps if pair.jitter_ps and pair.jitter_ps > 0 else max(half / 10, tick)
    width = max(tick, int(sigma0 / 3) // tick * tick)
    # Keep the bin count bounded for wide windows against narrow peaks.
    width = max(width, -(-(2 * half) // 4000) // tick * tick or tick)
    return shift, half, width, sigma0


def _edges(hist: np.ndarray, half: int, width: int) -> np.ndarray:
    edges = -half + width * np.arange(hist.size + 1, dtype=float)
    edges[-1] = min(edges[-1], float(half))
    return edges


def fit_peak(
    hist: np.ndarray,
    edges: np.ndarray,
    tick: float,
    d0: float,
    sigma0: float,
    offsets: Sequence[float],
    *,
    polish: bool = False,
    weights: Sequence[float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Poisson maximum-likelihood fit of the forward model's shape (`_ShapeModel`) to a histogram:
    (theta, expected counts per bin). theta = (log bL, log A, d, log sigma, log bR).

    Starts: the jitter at sigma0, twice and half it, with the background from the bins more than
    5 sigma0 from d0 (or a floor), and once more with the background from the three outer bins
    each side. The floor start alone stalled on a clean histogram whose peak fills the window (p
    about 1e-206 at one seed; council round 2).

    With `polish`, a Nelder-Mead search continues from the best start. L-BFGS-B often stops early
    on the background step at d, leaving the delay a few sd off (-386 ps for a -400 ps peak with
    6,000 pairs); `pair_pvalue` polishes when its p is low, and the window capture always does."""
    hist = np.asarray(hist, dtype=float)
    lo, hi = edges[:-1], edges[1:]
    width = float(np.max(hi - lo))
    half = float(-edges[0])
    model = _ShapeModel(hist, edges, tick, offsets, weights)
    floor = 0.5 / (2 * half)
    far = np.abs(0.5 * (lo + hi) - d0) > 5 * sigma0 + max((abs(o) for o in offsets), default=0.0)
    bg = float(np.median(hist[far])) / width if far.any() else floor
    bg = max(bg, floor)
    outer = max(float(np.median(np.r_[hist[:3], hist[-3:]])) / width, floor)
    total = float(hist.sum()) + 1.0
    lb = math.log(1e-6 / (2 * half))
    bounds = [
        (lb, math.log(10 * total / width)),
        (0.0, math.log(10 * total)),
        (-half, half),
        (math.log(max(tick / 100, 0.01)), math.log(2 * half)),
        (lb, math.log(10 * total / width)),
    ]
    starts = [(bg, s) for s in dict.fromkeys((sigma0, sigma0 * 2, sigma0 / 2))]
    if outer != bg:
        starts.append((outer, sigma0))
    best = None
    for b0, s_start in starts:
        peak = max(float(hist.sum()) - b0 * 2 * half, 1.0)
        x0 = np.array([math.log(b0), math.log(peak), d0, math.log(max(s_start, tick / 50)), math.log(b0)])
        x0 = np.clip(x0, [lo for lo, _ in bounds], [hi for _, hi in bounds])
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            r = optimize.minimize(model.nll, x0, method="L-BFGS-B", bounds=bounds, options={"ftol": 1e-14, "maxiter": 2000})
        if best is None or (np.isfinite(r.fun) and r.fun < best.fun):
            best = r
    if polish:
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            r = optimize.minimize(
                model.nll, best.x, method="Nelder-Mead", bounds=bounds,
                options={"maxiter": 5000, "xatol": 1e-4, "fatol": 1e-6},
            )
        if np.isfinite(r.fun) and r.fun < best.fun:
            best = r
    return best.x, model.mu(best.x)


def pair_pvalue(pair: PairData) -> float | None:
    """The grouped Pearson p-value of the forward model's shape on one histogram, or None when
    the region holds too few coincidences to check."""
    a = np.asarray(pair.a, dtype=np.int64)
    b = np.asarray(pair.b, dtype=np.int64)
    if a.size == 0 or b.size == 0:
        return None
    tick = _tick([pair])
    shift, half, width, sigma0 = _layout(pair, tick)
    shift, tooth = _align([pair], shift, half, tick)
    hist, groups = _histogram(pair, shift, width, half, tick)
    if hist.sum() < MIN_COUNTS:
        return None
    offsets, weights = _mixture(pair.offsets_ps, groups)
    edges = _edges(hist, half, width)
    d0 = float(pair.delay_ps) + (pair.declared_delay_ps - shift) - tooth
    theta, mu = fit_peak(hist, edges, tick, d0, sigma0, offsets, weights=weights)
    p = _grouped_pearson(hist, mu)
    if p >= POLISH_BELOW:
        return p
    # L-BFGS-B stalls when the start delay is tens of ps off (p 2e-6 on a clean histogram, 0.69
    # polished; council round 3). Polishing only lowers the NLL, so the test is the same with a
    # better optimiser; the highest bin is a second start.
    model = _ShapeModel(hist, edges, tick, offsets, weights)
    best, best_nll = (theta, mu), model.nll(theta)
    k = int(np.argmax(hist))
    for start in (d0, 0.5 * float(edges[k] + edges[k + 1])):
        th, m = fit_peak(hist, edges, tick, start, sigma0, offsets, polish=True, weights=weights)
        nll = model.nll(th)
        if nll < best_nll:
            best, best_nll = (th, m), nll
    return _grouped_pearson(hist, best[1])


@dataclass(frozen=True)
class PeakFit:
    """A pooled pair peak: offset from the declared delay and jitter (ps), with standard errors."""

    offset_ps: float
    sigma_ps: float
    offset_se_ps: float
    sigma_se_ps: float


#: Refinement passes of `pooled_peak` after the coarse fit.
POOL_PASSES = 3
#: Most bins in a refined pooled histogram.
POOL_MAX_BINS = 4000


def _summed(pairs: Sequence[PairData], shift: int, width: int, half: int, tick: int) -> tuple[np.ndarray, dict]:
    """The pairs' summed histogram and their pooled T3 groups (`_histogram`)."""
    total, groups = 0.0, {}
    for p in pairs:
        hist, found = _histogram(p, shift, width, half, tick)
        total = total + hist
        for g, n in found.items():
            groups[g] = groups.get(g, 0.0) + n
    return total, groups


def pooled_peak(pairs: Sequence[PairData]) -> PeakFit | None:
    """The pair peak of `pairs` (one channel pair and declared delay, such as one outcome pair over
    a Bell run's settings), fitted to their summed histograms, or None below `MIN_COUNTS` within
    the window.

    A coarse fit over |delay - declared| < W comes first (bins from the jitter, or W/30 without
    one; it starts at the pairs' given delay, or at the highest bin). Then up to `POOL_PASSES`
    refits, each from the previous one, until the bin width stops changing: the region extends
    5 sigma beyond the window, so the background is fitted outside the peak (a peak as wide as the
    window traded its offset against the background), and the bins are sigma/3 wide, so they
    resolve the peak (a 30 ps peak sat in one or two W/30 bins). Standard errors: sigma / sqrt(A)
    for the offset and sigma / sqrt(2 A) for the jitter, A the fitted peak area (council round 3).
    """
    pairs = [p for p in pairs if np.asarray(p.a).size and np.asarray(p.b).size]
    if not pairs:
        return None
    first = pairs[0]
    if any(p.declared_delay_ps != first.declared_delay_ps or p.window_ps != first.window_ps for p in pairs):
        raise ValueError("pooled pairs must share the declared delay and window")
    tick = _tick(pairs)
    shift, half, width, sigma0 = _layout(first, tick)
    shift, tooth = _align(pairs, shift, half, tick)
    hist, groups = _summed(pairs, shift, width, half, tick)
    if hist.sum() < MIN_COUNTS:
        return None
    offsets, weights = _mixture(first.offsets_ps, groups)
    edges = _edges(hist, half, width)
    if first.jitter_ps:
        d0 = float(first.delay_ps) + (first.declared_delay_ps - shift) - tooth
    else:
        k = int(np.argmax(hist))
        d0 = 0.5 * float(edges[k] + edges[k + 1])
    theta, _ = fit_peak(hist, edges, tick, d0, sigma0, offsets, polish=True, weights=weights)
    w_ps = int(first.window_ps)
    width = None
    for _ in range(POOL_PASSES):
        d, sigma = float(theta[2]), math.exp(float(theta[3]))
        sigma = min(max(sigma, tick / 2), 4.0 * w_ps)
        half = -(-int(w_ps + 5 * sigma) // tick) * tick
        new_width = max(tick, int(sigma / 3) // tick * tick, -(-(2 * half) // POOL_MAX_BINS) // tick * tick)
        if new_width == width:
            break
        width = new_width
        hist, groups = _summed(pairs, shift, width, half, tick)
        offsets, weights = _mixture(first.offsets_ps, groups)
        theta, _ = fit_peak(hist, _edges(hist, half, width), tick, d, sigma, offsets, polish=True, weights=weights)
    sigma = math.exp(float(theta[3]))
    area = math.exp(float(theta[1]))
    return PeakFit(
        float(theta[2]) + tooth - (first.declared_delay_ps - shift), sigma, sigma / math.sqrt(area),
        sigma / math.sqrt(2.0 * area),
    )


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
            f"window that the model lacks (crosstalk, an IRF tail, drift) would bias the corrected value. "
            f"{override_hint(analysis, OVERRIDE)}"
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
