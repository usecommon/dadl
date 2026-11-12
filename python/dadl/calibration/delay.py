"""Relative delay and combined timing jitter of two channels from a pair source.

Method `xcorr-gaussian-mle`: the histogram of delays b - a
(`dadl.histogram.cross_correlation`) is fitted with a Gaussian peak on a flat
background, with Poisson likelihood per bin:

    mu_k = beta * w_k + A * [Phi((e_{k+1} - d) / sigma) - Phi((e_k - d) / sigma)]

where e_k are the bin edges, w_k the bin widths, d the relative delay and
sigma the combined jitter (detector, TDC and cable on both channels; sigma^2
is the sum of the two channels' variances). The fit runs on a window of the
histogram around the peak, so dead-time and afterpulse structure tens of
nanoseconds away does not enter. Intervals and draws come from a Laplace
posterior in (log beta, log A, d, log sigma); the background and peak area
are nuisance parameters and are not stored.

The predictive check is a Pearson chi-square on the fit window, with
adjacent bins merged until each expects at least 5 counts.

Phase 2 limits: a Gaussian IRF only (no exponential tail, which shows as a
failing predictive check); a flat background, which suits CW sources (for
pulsed sources the side peaks must lie outside the range).
"""

from __future__ import annotations

import math

import numpy as np
from scipy import optimize
from scipy import special
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.base import N_DRAWS, Calibration
from dadl.data import ChannelId, tags_sha256
from dadl.errors import InsufficientCounts
from dadl.histogram import cross_correlation
from dadl.stats import Laplace, numerical_hessian

KIND = "delay-jitter"
METHOD = "xcorr-gaussian-mle"
METHOD_VERSION = "2"
#: Coincidences above background the peak must hold before a fit is attempted.
MIN_PEAK = 100
#: Half-width of the fit window around the peak, in units of the starting sigma.
WINDOW_SIGMAS = 25
MIN_EXPECTED = 5.0

UNITS = {"delay_ps": "ps", "jitter_sigma_ps": "ps"}

_SQRT2PI = math.sqrt(2 * math.pi)


def _edges(n_bins: int, bin_width_ps: int, range_ps: int) -> np.ndarray:
    """Model bin edges. Delays are whole ps, so the bin [e, e + w) holds the integers e .. e + w - 1,
    which a continuous delay fills from e - 1/2 to e + w - 1/2 (ignoring this biased v1 by +0.5 ps)."""
    e = -range_ps + bin_width_ps * np.arange(n_bins + 1, dtype=float)
    e[-1] = min(e[-1], float(range_ps))
    return e - 0.5


class _PeakModel:
    """Poisson negative log-likelihood of the histogram under peak + flat background."""

    def __init__(self, counts: np.ndarray, edges: np.ndarray) -> None:
        self.n = counts.astype(float)
        self.lo = edges[:-1]
        self.hi = edges[1:]
        self.w = self.hi - self.lo

    def _parts(self, theta: np.ndarray):
        log_b, log_a, d, log_s = theta
        s = math.exp(log_s)
        z_hi, z_lo = (self.hi - d) / s, (self.lo - d) / s
        g = special.ndtr(z_hi) - special.ndtr(z_lo)
        bg = math.exp(log_b) * self.w
        peak = math.exp(log_a) * g
        mu = np.maximum(bg + peak, 1e-300)
        return s, z_hi, z_lo, bg, peak, mu

    def mu(self, theta: np.ndarray) -> np.ndarray:
        return self._parts(theta)[-1]

    def nll(self, theta: np.ndarray) -> float:
        mu = self.mu(theta)
        return float(np.sum(mu - self.n * np.log(mu)))

    def grad(self, theta: np.ndarray) -> np.ndarray:
        s, z_hi, z_lo, bg, peak, mu = self._parts(theta)
        r = 1.0 - self.n / mu
        a = math.exp(theta[1])
        phi_hi = np.exp(-0.5 * z_hi**2) / _SQRT2PI
        phi_lo = np.exp(-0.5 * z_lo**2) / _SQRT2PI
        d_d = -a * (phi_hi - phi_lo) / s
        d_logs = -a * (z_hi * phi_hi - z_lo * phi_lo)
        return np.array([np.sum(r * bg), np.sum(r * peak), np.sum(r * d_d), np.sum(r * d_logs)])


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
        return float("nan")
    o, e = np.array(obs), np.array(exp)
    return float(_st.chi2.sf(float(np.sum((o - e) ** 2 / e)), dof))


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
) -> Calibration:
    """Fit the delay of channel b relative to a (b - a, ps) and the combined jitter sigma.

    `a` and `b` are sorted int64 ps tags from one acquisition of a pair
    source. The peak must lie within `range_ps` of zero.
    """
    if len(channels) != 2:
        raise ValueError("channels must be (channel a, channel b)")
    a = _checks.timestamps(a, "a")
    b = _checks.timestamps(b, "b")
    hist = cross_correlation(a, b, bin_width_ps, range_ps)
    edges = _edges(hist.size, bin_width_ps, range_ps)
    counts = hist.astype(float)
    background = float(np.median(counts))
    excess = counts - background
    if float(excess.sum()) < MIN_PEAK or counts.max() - background < 3 * math.sqrt(max(background, 1.0)):
        raise InsufficientCounts(
            f"no coincidence peak: need at least {MIN_PEAK} coincidences above background within +-{range_ps} ps"
        )

    # Starting values from the histogram's moments around its highest bin.
    k = int(np.argmax(counts))
    centres = 0.5 * (edges[:-1] + edges[1:])
    near = slice(max(0, k - 50), min(hist.size, k + 51))
    wts = np.clip(excess[near], 0, None)
    d0 = float(np.sum(wts * centres[near]) / wts.sum())
    s0 = max(float(math.sqrt(np.sum(wts * (centres[near] - d0) ** 2) / wts.sum())), bin_width_ps / 2)
    half = max(WINDOW_SIGMAS * s0, 50 * bin_width_ps)
    sel = (edges[1:] > d0 - half) & (edges[:-1] < d0 + half)
    idx = np.flatnonzero(sel)
    model = _PeakModel(counts[idx], edges[idx[0] : idx[-1] + 2])
    widths = model.w
    b0 = max(background / bin_width_ps, 0.5 / float(widths.sum()))
    a0 = max(float(np.sum(counts[idx] - b0 * widths)), 1.0)
    x0 = np.array([math.log(b0), math.log(a0), d0, math.log(s0)])

    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        res = optimize.minimize(model.nll, x0, jac=model.grad, method="BFGS", options={"gtol": 1e-6, "maxiter": 2000})
    notes: list[str] = []
    theta = res.x
    hess = numerical_hessian(model.grad, theta)
    laplace = Laplace.from_hessian(theta, hess)
    n_peak = math.exp(theta[1])
    converged = bool(res.success) or float(np.max(np.abs(model.grad(theta)))) < 1e-3 * math.sqrt(n_peak)
    sigma = math.exp(theta[3])
    if not (edges[0] + 5 * sigma < theta[2] < edges[-1] - 5 * sigma):
        notes.append("the peak lies within 5 sigma of the edge of the histogram; widen range_ps")
        converged = False
    if sigma < bin_width_ps:
        notes.append(f"the fitted jitter ({sigma:.3g} ps) is below the bin width; use narrower bins")

    rng = np.random.default_rng(seed)
    draws = laplace.draws(rng, n_draws)
    estimates = {
        "delay_ps": laplace.interval(2, lambda x: x),
        "jitter_sigma_ps": laplace.interval(3, math.exp),
    }
    notes.append(
        f"peak {n_peak:.4g} coincidences on a background of {math.exp(theta[0]) * bin_width_ps:.3g} per "
        f"{bin_width_ps} ps bin; fit window {model.lo[0]:.0f} to {model.hi[-1]:.0f} ps"
    )
    return Calibration(
        kind=KIND,
        method=METHOD,
        method_version=METHOD_VERSION,
        channels=tuple(channels),
        estimates=estimates,
        units=dict(UNITS),
        samples=np.column_stack([draws[:, 2], np.exp(draws[:, 3])]),
        sample_params=("delay_ps", "jitter_sigma_ps"),
        converged=converged,
        ppc_pvalue=_grouped_pearson(model.n, model.mu(theta), 4),
        inputs=((tags_sha256(a), "raw-tags"), (tags_sha256(b), "raw-tags")),
        backend=backend,
        conditions=conditions or {},
        notes=tuple(notes),
    )


__all__ = ["KIND", "calibrate_delay"]
