"""Detector models used by the calibrations.

A `DetectorModel` names its parameters, simulates data from them and scores
data against them. The calibration fits use the model's likelihood; the
simulator in `dadl.sim` is an independent, event-by-event implementation of
the same physics, so a model that disagrees with it shows up as a coverage
failure rather than passing by construction.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class DetectorModel(Protocol):
    """The interface every calibration model implements."""

    #: Parameter names, in the order used by `log_likelihood`.
    params: tuple[str, ...]

    def log_likelihood(self, data: np.ndarray, theta: np.ndarray) -> float:
        """Log-likelihood of reduced data (e.g. inter-arrival gaps) at parameter vector `theta`."""
        ...

    def simulate(self, theta: np.ndarray, size: int, rng: np.random.Generator) -> np.ndarray:
        """Draw reduced data from the model (used for predictive checks)."""
        ...


class DarkAfterpulseModel:
    """Inter-arrival gaps of a dark detector past its dead time.

    After each recorded avalanche and its dead time, the next avalanche is the
    first of a Poisson dark count (rate r) and, with probability p, that
    avalanche's own afterpulse after an exponential delay (mean tau_a). With
    s = gap - dead time >= 0:

        S(s) = exp(-r s) (1 - p + p exp(-s / tau_a))
        f(s) = exp(-r s) [r (1 - p) + p (r + 1 / tau_a) exp(-s / tau_a)]

    Parameters are fitted in transformed coordinates theta = (log r, logit p,
    log tau_a), with r in 1/ps and tau_a in ps. Afterpulses still pending
    from earlier avalanches are ignored; their effect is of order r * tau_a.
    """

    params = ("log_r", "logit_p", "log_tau_ap")

    @staticmethod
    def natural(theta: np.ndarray) -> tuple[float, float, float]:
        # Clip to a wide physical range so line searches cannot overflow.
        x, y, z = (float(v) for v in np.clip(theta, (-80.0, -40.0, -20.0), (0.0, 40.0, 45.0)))
        return float(np.exp(x)), float(1.0 / (1.0 + np.exp(-y))), float(np.exp(z))

    def _terms(self, s: np.ndarray, theta: np.ndarray):
        r, p, ta = self.natural(theta)
        e = np.exp(-s / ta)
        a = r * (1 - p)
        b = p * (r + 1 / ta)
        d = a + b * e
        return r, p, ta, e, d

    def log_likelihood(self, data: np.ndarray, theta: np.ndarray) -> float:
        s = np.asarray(data, dtype=float)
        r, _p, _ta, _e, d = self._terms(s, theta)
        return float(np.sum(np.log(d)) - r * s.sum())

    def gradient(self, data: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """Gradient of the negative log-likelihood in transformed coordinates."""
        s = np.asarray(data, dtype=float)
        r, p, ta, e, d = self._terms(s, theta)
        b = p * (r + 1 / ta)
        d_r = (1 - p) + p * e
        d_p = -r + (r + 1 / ta) * e
        d_ta = (e / ta**2) * (-p + b * s)
        return np.array(
            [
                r * (s.sum() - np.sum(d_r / d)),
                -p * (1 - p) * np.sum(d_p / d),
                -ta * np.sum(d_ta / d),
            ]
        )

    def survival(self, s: np.ndarray, theta: np.ndarray) -> np.ndarray:
        r, p, ta = self.natural(theta)
        s = np.asarray(s, dtype=float)
        return np.exp(-r * s) * (1 - p + p * np.exp(-s / ta))

    def density_at_zero(self, theta: np.ndarray) -> float:
        r, p, ta = self.natural(theta)
        return r * (1 - p) + p * (r + 1 / ta)

    def simulate(self, theta: np.ndarray, size: int, rng: np.random.Generator) -> np.ndarray:
        r, p, ta = self.natural(theta)
        dark = rng.exponential(1 / r, size)
        ap = np.where(rng.random(size) < p, rng.exponential(ta, size), np.inf)
        return np.minimum(dark, ap)
