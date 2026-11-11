"""Time-tag simulator with known, injected detector imperfections.

Every parameter that the calibration engine will later estimate is set
explicitly here, so tests can compare estimates against the truth.

Units: time tags are int64 picoseconds; rates are per second.

Phase 1 approximations, to be replaced by the self-exciting detection model
in docs/design.md:
- afterpulses are generated for one generation only (no cascades);
- crosstalk and breakdown flash are not simulated yet.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from dadl.deadtime import apply_nonparalyzable

PS_PER_S = 10**12

#: Analyzer angles in degrees for the standard CHSH settings on a Phi+ state.
CHSH_SETTINGS: dict[str, float] = {"a": 0.0, "a'": 45.0, "b": 22.5, "b'": 67.5}
CHSH_PAIRS: tuple[tuple[str, str], ...] = (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))
CHANNELS: tuple[str, ...] = ("A+", "A-", "B+", "B-")


@dataclass(frozen=True)
class DetectorParams:
    """Ground truth for one detector channel."""

    efficiency: float = 0.5  # system efficiency (detector x coupling)
    dark_rate: float = 100.0  # counts per second
    dead_time_ps: int = 22_000  # non-paralyzable
    jitter_ps: float = 350.0  # Gaussian sigma
    afterpulse_prob: float = 0.0  # per recorded avalanche
    afterpulse_tau_ps: float = 50_000.0  # exponential delay after the dead time
    delay_ps: int = 0  # fixed channel delay

    def __post_init__(self) -> None:
        if not 0.0 <= self.efficiency <= 1.0:
            raise ValueError("efficiency must be in [0, 1]")
        if self.dark_rate < 0 or self.dead_time_ps < 0 or self.jitter_ps < 0:
            raise ValueError("rates, dead time and jitter must be non-negative")
        if not 0.0 <= self.afterpulse_prob < 1.0:
            raise ValueError("afterpulse probability must be in [0, 1)")


@dataclass(frozen=True)
class PairSource:
    """A CW polarization-entangled SPDC source producing a Werner-like Phi+ state."""

    pair_rate: float = 1e5  # pairs per second at the crystal
    visibility: float = 0.95  # E(alpha, beta) = V cos 2(alpha - beta)

    def __post_init__(self) -> None:
        if self.pair_rate < 0:
            raise ValueError("pair rate must be non-negative")
        if not 0.0 <= self.visibility <= 1.0:
            raise ValueError("visibility must be in [0, 1]")


def poisson_times(rate: float, duration_s: float, rng: np.random.Generator) -> np.ndarray:
    """Sorted float64 arrival times (ps) of a homogeneous Poisson process."""
    n = rng.poisson(rate * duration_s)
    return np.sort(rng.uniform(0.0, duration_s * PS_PER_S, n))


def detect(
    photon_times_ps: np.ndarray,
    det: DetectorParams,
    duration_s: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Turn photon arrival times into recorded int64 time tags for one detector."""
    photons = np.asarray(photon_times_ps, dtype=np.float64)
    kept = photons[rng.random(photons.size) < det.efficiency]
    if det.jitter_ps > 0:
        kept = kept + rng.normal(0.0, det.jitter_ps, kept.size)
    kept = kept + det.delay_ps
    darks = poisson_times(det.dark_rate, duration_s, rng)
    arrivals = np.sort(np.rint(np.concatenate([kept, darks]))).astype(np.int64)
    recorded = apply_nonparalyzable(arrivals, det.dead_time_ps)

    if det.afterpulse_prob > 0 and recorded.size:
        mask = rng.random(recorded.size) < det.afterpulse_prob
        n_ap = int(mask.sum())
        if n_ap:
            delays = det.dead_time_ps + rng.exponential(det.afterpulse_tau_ps, n_ap)
            afterpulses = np.rint(recorded[mask] + delays).astype(np.int64)
            merged = np.sort(np.concatenate([recorded, afterpulses]))
            recorded = apply_nonparalyzable(merged, det.dead_time_ps)
    return recorded


def simulate_chsh_setting(
    source: PairSource,
    alpha_deg: float,
    beta_deg: float,
    detectors: dict[str, DetectorParams],
    duration_s: float,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """Recorded tags on channels A+, A-, B+, B- for one pair of analyzer angles."""
    missing = set(CHANNELS) - set(detectors)
    if missing:
        raise ValueError(f"missing detector parameters for {sorted(missing)}")
    t = poisson_times(source.pair_rate, duration_s, rng)
    n = t.size
    a_plus = rng.random(n) < 0.5
    p_same = 0.5 * (1.0 + source.visibility * math.cos(2.0 * math.radians(alpha_deg - beta_deg)))
    same = rng.random(n) < p_same
    b_plus = np.where(same, a_plus, ~a_plus)
    photons = {"A+": t[a_plus], "A-": t[~a_plus], "B+": t[b_plus], "B-": t[~b_plus]}
    return {ch: detect(photons[ch], detectors[ch], duration_s, rng) for ch in CHANNELS}


def simulate_chsh(
    source: PairSource,
    detectors: dict[str, DetectorParams],
    duration_s: float,
    rng: np.random.Generator,
    settings: dict[str, float] | None = None,
) -> dict[tuple[str, str], dict[str, np.ndarray]]:
    """One acquisition per CHSH setting pair, keyed like ("a", "b'")."""
    angles = settings or CHSH_SETTINGS
    return {
        (sa, sb): simulate_chsh_setting(source, angles[sa], angles[sb], detectors, duration_s, rng)
        for sa, sb in CHSH_PAIRS
    }
