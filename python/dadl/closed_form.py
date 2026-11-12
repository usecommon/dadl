"""Closed-form results used as model-independent test oracles."""

from __future__ import annotations

import math


def nonparalyzable_measured_rate(true_rate: float, tau_s: float) -> float:
    """m = n / (1 + n tau)."""
    return true_rate / (1.0 + true_rate * tau_s)


def paralyzable_measured_rate(true_rate: float, tau_s: float) -> float:
    """m = n exp(-n tau)."""
    return true_rate * math.exp(-true_rate * tau_s)


def series_measured_rate(true_rate: float, tau_p_s: float, tau_np_s: float) -> float:
    """Exact rate for a paralyzable stage (tau_p) followed by a non-paralyzable stage (tau_np).

    With m1 = n exp(-n tau_p) the output of the paralyzable stage:
    m = m1 if tau_np <= tau_p, else m = m1 / (1 + m1 (tau_np - tau_p)).
    """
    m1 = paralyzable_measured_rate(true_rate, tau_p_s)
    if tau_np_s <= tau_p_s:
        return m1
    return m1 / (1.0 + m1 * (tau_np_s - tau_p_s))


def patil_usman_hybrid_rate(true_rate: float, tau_p_s: float, tau_np_s: float) -> float:
    """Empirical hybrid m = n exp(-n tau_p) / (1 + n tau_np) (Patil and Usman).

    A phenomenological fit form, not the exact series model: it differs from a
    simulated paralyzable -> non-paralyzable chain by up to ~30 % at high rate.
    Do not use it as a test oracle; see `series_measured_rate`.
    """
    return true_rate * math.exp(-true_rate * tau_p_s) / (1.0 + true_rate * tau_np_s)


def dark_afterpulse_rate(dark_rate: float, afterpulse_prob: float, afterpulse_tau_s: float, dead_time_s: float) -> float:
    """Recorded rate of a dark detector with non-paralyzable dead time and cascading afterpulses.

    Each recorded avalanche blocks `dead_time_s`; after it, the next event is
    the first of a Poisson dark count (rate r) and, with probability p, its own
    afterpulse after an exponential delay (mean tau_a). The gaps are a renewal
    process with survival S(s) = exp(-r s) (1 - p + p exp(-s / tau_a)) past the
    dead time, so the mean gap is tau + (1 - p) / r + p / (r + 1 / tau_a).
    """
    r, p, ta = dark_rate, afterpulse_prob, afterpulse_tau_s
    mean_gap = dead_time_s + (1.0 - p) / r + p / (r + 1.0 / ta)
    return 1.0 / mean_gap


def chsh_werner(visibility: float) -> float:
    """S = 2 sqrt(2) V for a Werner-like Phi+ state at the optimal settings."""
    return 2.0 * math.sqrt(2.0) * visibility


def chsh_with_accidentals(visibility: float, true_coinc: float, accidentals: float) -> float:
    """Uniform accidentals scale every correlator, so S = 2 sqrt(2) V C / (C + A)."""
    return chsh_werner(visibility) * true_coinc / (true_coinc + accidentals)


def grid_capture_probability(jitter_sigma_ps: float, width_ps: float) -> float:
    """Probability that both photons of a pair land in the same fixed slot.

    Each arm has independent Gaussian jitter of sigma `jitter_sigma_ps`, the
    declared delays are exact, and the pair time is uniform within a slot.
    With D the arrival difference (sigma s = sigma * sqrt(2)) and w the width,
    P = E[(1 - |D| / w)+] = erf(u / sqrt 2) - (2 s / (w sqrt(2 pi))) (1 - exp(-u^2 / 2)), u = w / s.
    """
    if width_ps <= 0:
        raise ValueError("width must be positive")
    if not jitter_sigma_ps >= 0:
        raise ValueError("jitter sigma must be non-negative")
    s = jitter_sigma_ps * math.sqrt(2.0)
    if s == 0:
        return 1.0
    u = width_ps / s
    return math.erf(u / math.sqrt(2.0)) - (2.0 * s / (width_ps * math.sqrt(2.0 * math.pi))) * (
        1.0 - math.exp(-0.5 * u * u)
    )


def chsh_expected_grid(
    *,
    pair_rate_per_s: float,
    duration_s: float,
    efficiency: float,
    dark_rate_cps: float,
    jitter_sigma_ps: float,
    width_ps: float,
    visibility: float,
) -> float:
    """Expected raw S for the simulator with identical detectors, no dead time or afterpulsing.

    True coincidences per setting: C = N T eta^2 P_capture (summed over the four
    outcome combinations). Accidentals per combination: r^2 w T with singles
    rate r = N eta / 2 + d; uniform across combinations, so they add nothing to
    the correlator numerator.
    """
    c = pair_rate_per_s * duration_s * efficiency**2 * grid_capture_probability(jitter_sigma_ps, width_ps)
    r = pair_rate_per_s * efficiency / 2.0 + dark_rate_cps
    a = 4.0 * r * r * (width_ps * 1e-12) * duration_s
    return chsh_with_accidentals(visibility, c, a)
