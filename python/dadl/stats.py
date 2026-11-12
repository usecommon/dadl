"""Interval estimates, posteriors and coverage checks used by the calibrations.

Phase 2 uses maximum likelihood with SciPy. Posterior samples come from the
exact conjugate posterior where one exists (Beta for a binomial ratio, Gamma
for a Poisson rate) and otherwise from a Laplace approximation: a Gaussian at
the maximum-likelihood point in transformed coordinates (log or logit), so
samples respect the parameter's range. Full MCMC (NUTS) arrives with the
joint model in phase 3.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy import stats as _st

#: Two-sided level of every interval DADL reports.
LEVEL = 0.95


@dataclass(frozen=True)
class Interval:
    """A point estimate with a two-sided interval at `level`."""

    value: float
    low: float
    high: float
    level: float = LEVEL

    def contains(self, x: float) -> bool:
        return self.low <= x <= self.high

    def __str__(self) -> str:
        return f"{self.value:.6g} [{self.low:.6g}, {self.high:.6g}]"


def _z(level: float) -> float:
    return float(_st.norm.ppf(0.5 + level / 2))


def wilson(k: int, n: int, level: float = LEVEL) -> Interval:
    """Wilson score interval for a binomial proportion k/n."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError(f"need 0 <= k <= n and n > 0, got k={k}, n={n}")
    z = _z(level)
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return Interval(p, max(0.0, centre - half), min(1.0, centre + half), level)


def clopper_pearson(k: int, n: int, level: float = LEVEL) -> Interval:
    """Exact (conservative) binomial interval for k/n."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError(f"need 0 <= k <= n and n > 0, got k={k}, n={n}")
    a = 1 - level
    low = 0.0 if k == 0 else float(_st.beta.ppf(a / 2, k, n - k + 1))
    high = 1.0 if k == n else float(_st.beta.ppf(1 - a / 2, k + 1, n - k))
    return Interval(k / n, low, high, level)


def beta_posterior(k: int, n: int, rng: np.random.Generator, n_draws: int, level: float = LEVEL) -> tuple[Interval, np.ndarray]:
    """Jeffreys Beta(k + 1/2, n - k + 1/2) posterior for a proportion: interval and draws."""
    dist = _st.beta(k + 0.5, n - k + 0.5)
    a = 1 - level
    return Interval(k / n, float(dist.ppf(a / 2)), float(dist.ppf(1 - a / 2)), level), dist.rvs(n_draws, random_state=rng)


def gamma_posterior(
    count: int, exposure: float, rng: np.random.Generator, n_draws: int, level: float = LEVEL
) -> tuple[Interval, np.ndarray]:
    """Jeffreys Gamma(count + 1/2, exposure) posterior for a Poisson rate: interval and draws."""
    if exposure <= 0:
        raise ValueError("exposure must be positive")
    dist = _st.gamma(count + 0.5, scale=1.0 / exposure)
    a = 1 - level
    return (
        Interval(count / exposure, float(dist.ppf(a / 2)), float(dist.ppf(1 - a / 2)), level),
        dist.rvs(n_draws, random_state=rng),
    )


def numerical_hessian(grad: Callable[[np.ndarray], np.ndarray], x: np.ndarray, rel_step: float = 1e-4) -> np.ndarray:
    """Symmetrized central-difference Hessian from an analytic gradient."""
    x = np.asarray(x, dtype=float)
    k = x.size
    h = rel_step * np.maximum(1.0, np.abs(x))
    hess = np.empty((k, k))
    for i in range(k):
        e = np.zeros(k)
        e[i] = h[i]
        hess[:, i] = (grad(x + e) - grad(x - e)) / (2 * h[i])
    return 0.5 * (hess + hess.T)


@dataclass(frozen=True)
class Laplace:
    """Gaussian approximation N(mean, cov) in the fit's transformed coordinates."""

    mean: np.ndarray
    cov: np.ndarray

    @classmethod
    def from_hessian(cls, mean: np.ndarray, hessian: np.ndarray) -> Laplace:
        """Raise ValueError when the Hessian is not positive definite (the fit is not identified)."""
        try:
            np.linalg.cholesky(hessian)
        except np.linalg.LinAlgError as exc:
            raise ValueError("the Hessian is not positive definite; a parameter is not identified by the data") from exc
        return cls(np.asarray(mean, dtype=float), np.linalg.inv(hessian))

    def interval(self, i: int, transform: Callable[[float], float], level: float = LEVEL) -> Interval:
        """Wald interval in coordinate i, mapped through a monotone increasing `transform`."""
        sd = float(np.sqrt(self.cov[i, i]))
        z = _z(level)
        m = float(self.mean[i])
        return Interval(transform(m), transform(m - z * sd), transform(m + z * sd), level)

    def draws(self, rng: np.random.Generator, n_draws: int) -> np.ndarray:
        return rng.multivariate_normal(self.mean, self.cov, size=n_draws, method="cholesky")


def coverage_band(n_replicates: int, level: float = LEVEL, confidence: float = 0.99) -> tuple[float, float]:
    """Range of observed coverage consistent with true coverage `level` (binomial, two-sided `confidence`).

    For 1,000 replicates at 95% and 99% confidence this is 93.1% to 96.7%;
    `dadl.calibration.coverage.gate` splits the confidence across parameters.
    """
    a = 1 - confidence
    lo = _st.binom.ppf(a / 2, n_replicates, level) / n_replicates
    hi = _st.binom.ppf(1 - a / 2, n_replicates, level) / n_replicates
    return float(lo), float(hi)


def logit(p: float | np.ndarray) -> float | np.ndarray:
    return np.log(p) - np.log1p(-p)


def expit(x: float | np.ndarray) -> float | np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))
