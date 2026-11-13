"""Relative delay and combined timing jitter of two channels from a pair source.

Method `xcorr-gaussian-mle`: the histogram of delays b - a
(`dadl.histogram.cross_correlation`) is fitted with a Gaussian peak on a
background with one level on each side of the peak centre (dead time removes
background on one side only), with Poisson likelihood per bin:

    mu_k = beta_L * w_k^L + beta_R * w_k^R + A * [Phi((e_{k+1} - d) / sigma) - Phi((e_k - d) / sigma)]

(w_k^L and w_k^R the parts of bin k left and right of d).

where e_k are the bin edges, w_k the bin widths, d the relative delay and
sigma the combined jitter (detector, TDC and cable on both channels; sigma^2
is the sum of the two channels' variances). The fit runs on a window of the
histogram around the peak, so dead-time and afterpulse structure tens of
nanoseconds away does not enter. Intervals and draws come from a Laplace
posterior in (log beta, log A, d, log sigma); the background and peak area
are nuisance parameters and are not stored.

Tagger tick: delays come in whole ticks (1 ps for a HydraHarp, 5 ps for a
MultiHarp 160, 80 ps for a MultiHarp 150 N, 25 ps for a TimeHarp 260 P). The
tick is passed by the caller when the adapter knows it (`tick_ps`); otherwise
it is found from the data (the greatest common divisor of the tags, which a
non-tick offset between channels would spoil). The bin width is rounded up to a whole
number of ticks. Each channel's time is quantized to the tick with a phase
uniform over the tick, so the difference of two quantized times lands on tick
k with probability E[max(0, 1 - |D / tick - k|)]: the Gaussian peak
convolved with a triangle two ticks wide (method v5). The peak's share of a
bin [e, e + w) is then exact:

    (sigma / tick) [G(z_hi) - G(z_hi - tick / sigma) - G(z_lo) + G(z_lo - tick / sigma)],

with G(z) = z Phi(z) + phi(z) and z = (edge - d) / sigma, so sigma is the
physical jitter, without the tick's quantization noise. (v1-v4 used a box one
tick wide: with an 80 ps tick and 21 ps jitter the delay read up to 1.3 ps off
and covered 0-12%, with the check passing; tagger engineer review, cycle 2.)

The predictive check is a Pearson chi-square on the fit window, with
adjacent bins merged until each expects at least 5 counts.

Phase 2 limits: a Gaussian IRF only (an exponential tail biases the delay by
tens of ps and shows as a failing predictive check, which marks the
calibration a misfit); CW sources (for pulsed sources the side peaks must lie
outside the range); a delay calibrated at one count rate may shift at another
(record the rates; recalibrate after large rate changes).
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
from scipy import optimize
from scipy import special
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.base import N_DRAWS, Calibration, with_interruption_check, draws_rng
from dadl.data import ChannelId, tag_tick, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.histogram import cross_correlation
from dadl.stats import Interval, Laplace, numerical_hessian

KIND = "delay-jitter"
METHOD = "xcorr-gaussian-mle"
METHOD_VERSION = "6"
#: Coincidences above background the peak must hold before a fit is attempted.
MIN_PEAK = 100
#: Half-width of the fit window around the peak, in units of the starting sigma.
WINDOW_SIGMAS = 25
MIN_EXPECTED = 5.0

UNITS = {"delay_ps": "ps", "jitter_sigma_ps": "ps"}
DEFINITIONS = {
    "delay_ps": "relative delay of channel b after channel a, in the tagger's frame",
    "jitter_sigma_ps": "combined Gaussian timing jitter of both channels (sigma), without the tick's quantization",
}

_SQRT2PI = math.sqrt(2 * math.pi)


def _edges(n_bins: int, bin_width_ps: int, range_ps: int) -> np.ndarray:
    """Bin edges: bin [e, e + w) holds the measured delays e, e + tick, ..., e + w - tick."""
    e = -range_ps + bin_width_ps * np.arange(n_bins + 1, dtype=float)
    e[-1] = min(e[-1], float(range_ps))
    return e


def _g(z: np.ndarray) -> np.ndarray:
    """G(z) = z Phi(z) + phi(z), the antiderivative of the normal CDF."""
    return z * special.ndtr(z) + np.exp(-0.5 * z * z) / _SQRT2PI


class _PeakModel:
    """Poisson negative log-likelihood of the histogram: Gaussian peak on a two-level background.

    The background takes one level left of the peak centre and another right of it. Dead time
    makes this step: after a pair is detected the signal detector is dead for a dead time, which
    removes background on one side of the peak only. With a bright background a flat model
    shifted the delay by more than its interval (52% coverage at 8 Mcps; method v3).
    theta = (log background left, log peak area, delay, log sigma, log background right).
    """

    def __init__(self, counts: np.ndarray, edges: np.ndarray, tick: float = 1.0) -> None:
        self.n = counts.astype(float)
        self.lo = edges[:-1]
        self.hi = edges[1:]
        self.w = self.hi - self.lo
        self.tick = float(tick)

    def _parts(self, theta: np.ndarray):
        log_bl, log_a, d, log_s, log_br = (float(v) for v in theta)
        log_bl, log_br, log_a = min(log_bl, 700.0), min(log_br, 700.0), min(log_a, 700.0)
        log_s = min(max(log_s, -50.0), 700.0)
        s = math.exp(log_s)
        z_hi, z_lo = (self.hi - d) / s, (self.lo - d) / s
        t = self.tick / s
        g = (s / self.tick) * (_g(z_hi) - _g(z_hi - t) - _g(z_lo) + _g(z_lo - t))
        right = np.clip(self.hi - d, 0.0, self.w)
        left = self.w - right
        bl, br = math.exp(log_bl), math.exp(log_br)
        bg_l, bg_r = bl * left, br * right
        peak = math.exp(log_a) * g
        mu = np.maximum(bg_l + bg_r + peak, 1e-300)
        straddle = (self.lo < d) & (d < self.hi)
        return s, z_hi, z_lo, bg_l, bg_r, peak, mu, straddle, bl, br

    def mu(self, theta: np.ndarray) -> np.ndarray:
        return self._parts(theta)[6]

    def nll(self, theta: np.ndarray) -> float:
        mu = self.mu(theta)
        return float(np.sum(mu - self.n * np.log(mu)))

    def grad(self, theta: np.ndarray) -> np.ndarray:
        s, z_hi, z_lo, bg_l, bg_r, peak, mu, straddle, bl, br = self._parts(theta)
        r = 1.0 - self.n / mu
        a = math.exp(theta[1])
        t = self.tick / s

        def phi(z: np.ndarray) -> np.ndarray:
            return np.exp(-0.5 * z * z) / _SQRT2PI

        cdf = special.ndtr
        second = cdf(z_hi) - cdf(z_hi - t) - cdf(z_lo) + cdf(z_lo - t)
        dens = phi(z_hi) - phi(z_hi - t) - phi(z_lo) + phi(z_lo - t)
        # Moving d right moves bin width from the right level to the left one where d falls in a bin.
        d_d = -a * second / self.tick + np.where(straddle, bl - br, 0.0)
        d_logs = a * (s / self.tick) * dens
        return np.array([np.sum(r * bg_l), np.sum(r * peak), np.sum(r * d_d), np.sum(r * d_logs), np.sum(r * bg_r)])


def _grouped_pearson(n: np.ndarray, mu: np.ndarray, n_params: int) -> float:
    """Pearson chi-square p-value with adjacent bins merged until each expects >= MIN_EXPECTED."""
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
    dof = len(exp) - n_params
    if dof < 1:
        return _bootstrap_deviance(n, mu)
    o, e = np.array(obs), np.array(exp)
    return float(_st.chi2.sf(float(np.sum((o - e) ** 2 / e)), dof))


def _bootstrap_deviance(n: np.ndarray, mu: np.ndarray, n_sim: int = 1000) -> float:
    """Parametric-bootstrap p-value of the Poisson deviance, for peaks too sparse for a chi-square.

    Uses a fixed seed so the same data always give the same p-value. Parameters
    are not refitted per simulation, so the test is mildly conservative.
    """

    def deviance(k: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(k > 0, k * np.log(k / mu), 0.0)
        return 2.0 * np.sum(t - (k - mu), axis=-1)

    rng = np.random.default_rng(0)
    sims = deviance(rng.poisson(mu, size=(n_sim, mu.size)).astype(float))
    return float((1 + np.sum(sims >= deviance(n))) / (n_sim + 1))


def calibrate_delay(
    a: np.ndarray,
    b: np.ndarray,
    channels: tuple[ChannelId, ChannelId],
    *,
    range_ps: int = 20_000,
    bin_width_ps: int = 20,
    backend: str = "hardware",
    seed: int | None = None,
    n_draws: int = N_DRAWS,
    conditions: dict[str, float] | None = None,
    tick_ps: int | None = None,
    _refined: bool = False,
) -> Calibration:
    """Fit the delay of channel b relative to a (b - a, ps) and the combined jitter sigma.

    `a` and `b` are sorted int64 ps tags from one acquisition of a pair
    source. The peak must lie within `range_ps` of zero. When the fitted jitter
    is under two bin widths the fit is repeated once with bins of about a
    third of the jitter (whole tagger ticks): a peak one or two bins wide has a
    non-Gaussian likelihood and its intervals undercovered (92%).
    """
    if len(channels) != 2:
        raise ValueError("channels must be (channel a, channel b)")
    a = _checks.timestamps(a, "a")
    b = _checks.timestamps(b, "b")
    bin_width_ps = _checks.integer(bin_width_ps, "bin_width_ps", minimum=1)
    range_ps = _checks.integer(range_ps, "range_ps", minimum=1)
    found = tag_tick(a, b)
    notes: list[str] = []
    if tick_ps is None:
        tick = found
    else:
        tick = _checks.integer(tick_ps, "tick_ps", minimum=1)
        if found % tick:
            notes.append(
                f"the tags are not all on the {tick} ps tick grid (their common step is {found} ps): an offset "
                "between channels? The declared tick is used"
            )
    if bin_width_ps % tick:
        new = -(-bin_width_ps // tick) * tick
        notes.append(f"bin width raised from {bin_width_ps} to {new} ps, a whole number of {tick} ps tagger ticks")
        bin_width_ps = new
    if range_ps % tick:
        range_ps = -(-range_ps // tick) * tick
    hist = cross_correlation(a, b, bin_width_ps, range_ps)
    edges = _edges(hist.size, bin_width_ps, range_ps)
    counts = hist.astype(float)
    background = float(np.median(counts))
    excess = counts - background
    if float(excess.sum()) < MIN_PEAK or counts.max() - background < 3 * math.sqrt(max(background, 1.0)):
        raise InsufficientCounts(
            f"no coincidence peak: need at least {MIN_PEAK} coincidences above background within +-{range_ps} ps"
        )

    # Starting values from the histogram's moments around its highest bin, with the moment window
    # widened until it holds +-4 sigma: a fixed +-50 bins badly underestimated broad peaks.
    k = int(np.argmax(counts))
    centres = 0.5 * (edges[:-1] + edges[1:])
    d0, half_bins = float(centres[k]), 50
    for _ in range(8):
        near = slice(max(0, k - half_bins), min(hist.size, k + half_bins + 1))
        wts = np.clip(excess[near], 0, None)
        if wts.sum() <= 0:
            break
        d0 = float(np.sum(wts * centres[near]) / wts.sum())
        s0 = max(float(math.sqrt(np.sum(wts * (centres[near] - d0) ** 2) / wts.sum())), bin_width_ps / 2)
        need = int(math.ceil(4 * s0 / bin_width_ps))
        if need <= half_bins:
            break
        half_bins = min(need * 2, hist.size)
    half = max(WINDOW_SIGMAS * s0, 50 * bin_width_ps)
    sel = (edges[1:] > d0 - half) & (edges[:-1] < d0 + half)
    idx = np.flatnonzero(sel)
    model = _PeakModel(counts[idx], edges[idx[0] : idx[-1] + 2], tick)
    widths = model.w
    b0 = max(background / bin_width_ps, 0.5 / float(widths.sum()))
    a0 = max(float(np.sum(counts[idx] - b0 * widths)), 1.0)
    x0 = np.array([math.log(b0), math.log(a0), d0, math.log(s0), math.log(b0)])

    lo_e, hi_e = float(model.lo[0]), float(model.hi[-1])
    total = float(counts[idx].sum()) + 1.0
    bounds = [
        (math.log(1e-6 / float(widths.sum())), math.log(10 * total / bin_width_ps)),  # background per ps
        (0.0, math.log(10 * total)),  # peak area
        (lo_e, hi_e),  # delay
        (math.log(max(tick / 100, 0.01)), math.log(hi_e - lo_e)),  # sigma (Bayesian review, cycle 3: tick/4 pinned it)
        (math.log(1e-6 / float(widths.sum())), math.log(10 * total / bin_width_ps)),  # background right
    ]
    x0 = np.clip(x0, [b[0] for b in bounds], [b[1] for b in bounds])
    # Several starting jitters: against a coarse tick the likelihood has a second optimum at
    # sigma -> 0 (the bare tick triangle), which a single start can fall into.
    res = None
    for log_s in sorted({float(x0[3]), math.log(max(tick / 2, 1.0)), math.log(max(tick / 8, 0.5))}):
        start = x0.copy()
        start[3] = min(max(log_s, bounds[3][0]), bounds[3][1])
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            r = optimize.minimize(model.nll, start, jac=model.grad, method="L-BFGS-B", bounds=bounds,
                                  options={"ftol": 1e-15, "gtol": 1e-8, "maxiter": 5000})
        if res is None or (np.isfinite(r.fun) and r.fun < res.fun - 1e-9):
            res = r
    theta = res.x
    # Absolute steps: the delay's in ps scaled to the jitter, not to |d| (a 1e-4 relative step was
    # 100 ps at a 1 us delay, wider than the peak).
    steps = np.full(5, 1e-4)
    steps[2] = 1e-3 * max(math.exp(theta[3]), 0.01)
    hess = numerical_hessian(model.grad, theta, steps=steps)
    sigma_held = False
    try:
        laplace = Laplace.from_hessian(theta, hess)
    except ValueError:
        # Condition on what the data do not identify, in turn: the background (no background to
        # speak of; it barely affects the peak), then the jitter (far below the tick, the peak is
        # the tick's triangle and sigma has no curvature; Bayesian review, cycle 3).
        laplace = None
        for keep, held in (([1, 2, 3], "background"), ([0, 1, 2, 4], "jitter"), ([1, 2], "both")):
            try:
                inner = Laplace.from_hessian(theta[keep], hess[np.ix_(keep, keep)])
            except ValueError:
                continue
            cov = np.zeros((5, 5))
            cov[np.ix_(keep, keep)] = inner.cov
            laplace = Laplace(theta, cov)
            sigma_held = held in ("jitter", "both")
            if held in ("background", "both"):
                notes.append(
                    "the background is not identified (it fits to zero); the delay and jitter are conditional on it"
                )
            break
        if laplace is None:
            raise InsufficientCounts("the delay fit is not identified by these data (no curvature in the delay)")
    at_floor = theta[3] - bounds[3][0] < 1e-3
    # A log-jitter spread above 1 (a factor e) is a flat direction, not a measurement.
    sigma_held = sigma_held or at_floor or float(laplace.cov[3, 3]) > 1.0
    n_peak = math.exp(theta[1])
    converged = bool(res.success) or float(np.max(np.abs(model.grad(theta)))) < 1e-3 * math.sqrt(n_peak)
    sigma = math.exp(theta[3])
    if not (edges[0] + 5 * sigma < theta[2] < edges[-1] - 5 * sigma):
        notes.append("the peak lies within 5 sigma of the edge of the histogram; widen range_ps")
        converged = False
    finer = max(tick, int(sigma / 3) // tick * tick)
    if not _refined and sigma < 2 * bin_width_ps and finer < bin_width_ps:
        refined = calibrate_delay(
            a, b, channels, range_ps=range_ps, bin_width_ps=finer, backend=backend, seed=seed,
            n_draws=n_draws, conditions=conditions, tick_ps=tick_ps, _refined=True,
        )
        note = f"bins refined from {bin_width_ps} to {finer} ps because the jitter ({sigma:.3g} ps) is under two bin widths"
        return dataclasses.replace(refined, notes=(note, *refined.notes))
    if sigma < bin_width_ps:
        notes.append(f"the fitted jitter ({sigma:.3g} ps) is below the bin width; use narrower bins")
    if tick > sigma / 2:
        notes.append(
            f"the {tick} ps tick is coarse against the {sigma:.3g} ps jitter: the fit models the tick's "
            "quantization exactly, assuming each channel's phase within the tick is uniform"
        )

    rng = draws_rng(seed, a, b)
    unidentified: tuple[str, ...] = ()
    if sigma_held:
        # The jitter is far below the tick: report a profile-likelihood upper limit, condition the
        # delay on the fitted jitter, and draw the jitter uniformly below the limit.
        if not at_floor:
            cov = laplace.cov.copy()
            cov[3, :] = cov[:, 3] = 0.0
            laplace = Laplace(theta, cov)
        upper = _sigma_upper_limit(model, theta, bounds)
        draws = laplace.draws(rng, n_draws)
        draws[:, 3] = np.log(rng.uniform(math.exp(bounds[3][0]), upper, n_draws))
        estimates = {
            "delay_ps": laplace.interval(2, lambda x: x),
            "jitter_sigma_ps": Interval(upper, 0.0, upper),
        }
        unidentified = ("jitter_sigma_ps",)
        notes.append(
            f"the jitter is not resolved against the {tick} ps tick: 95% upper limit {upper:.3g} ps (profile "
            "likelihood); the delay is conditional on the fitted jitter, and jitter draws are uniform below the limit"
        )
    else:
        draws = laplace.draws(rng, n_draws)
        estimates = {
            "delay_ps": laplace.interval(2, lambda x: x),
            "jitter_sigma_ps": laplace.interval(3, math.exp),
        }
    notes.append(
        f"peak {n_peak:.4g} coincidences on a background of {math.exp(theta[0]) * bin_width_ps:.3g} (left) and "
        f"{math.exp(theta[4]) * bin_width_ps:.3g} (right) per "
        f"{bin_width_ps} ps bin; fit window {model.lo[0]:.0f} to {model.hi[-1]:.0f} ps"
    )
    return Calibration(
        kind=KIND,
        method=METHOD,
        method_version=METHOD_VERSION,
        channels=tuple(channels),
        estimates=estimates,
        units=dict(UNITS),
        definitions=dict(DEFINITIONS),
        samples=np.column_stack([draws[:, 2], np.exp(draws[:, 3])]),
        sample_params=("delay_ps", "jitter_sigma_ps"),
        converged=converged,
        ppc_pvalue=with_interruption_check(_grouped_pearson(model.n, model.mu(theta), 5), notes, a, b),
        inputs=((tags_sha256(a), "raw-tags", str(channels[0])), (tags_sha256(b), "raw-tags", str(channels[1]))),
        backend=backend,
        conditions={"tick_ps": float(tick), **(conditions or {})},
        unidentified=unidentified,
        stand_in=unidentified,
        notes=tuple(notes),
    )


def _sigma_upper_limit(model: _PeakModel, theta: np.ndarray, bounds: list) -> float:
    """One-sided 95% upper limit on sigma (ps): where the profile negative log-likelihood rises by
    1.35 (chi-square(1) at 90%, halved); 1.92 would be a 97.5% limit."""
    best = model.nll(theta)
    free = [0, 1, 2, 4]
    x = theta[free].copy()

    def profile(log_s: float) -> float:
        nonlocal x

        def f(y: np.ndarray) -> float:
            t = theta.copy()
            t[free] = y
            t[3] = log_s
            return model.nll(t)

        def g(y: np.ndarray) -> np.ndarray:
            t = theta.copy()
            t[free] = y
            t[3] = log_s
            return model.grad(t)[free]

        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            r = optimize.minimize(f, x, jac=g, method="L-BFGS-B", bounds=[bounds[i] for i in free])
        x = r.x
        return float(r.fun) - best

    lo = float(theta[3])
    hi = lo
    step = math.log(1.5)
    while profile(hi) < 1.35 and hi < bounds[3][1]:
        lo, hi = hi, hi + step
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if profile(mid) < 1.35:
            lo = mid
        else:
            hi = mid
    return math.exp(hi)


__all__ = ["KIND", "calibrate_delay"]
