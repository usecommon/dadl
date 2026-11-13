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
    settings_deg: dict[str, float] | None = None,
) -> float:
    """Expected raw S for the simulator with identical detectors, no dead time or afterpulsing.

    With `settings_deg` the ideal S is sum of sign x V cos 2(alpha - beta) over the four setting
    pairs (the simulator's Phi+ correlator); without it, the optimal 2 sqrt(2) V. (Bell review,
    cycle 2: other angles printed an expectation of 2.69 beside a measured S of 0.01.)

    True coincidences per setting: C = N T eta^2 P_capture (summed over the four
    outcome combinations). Accidentals per combination: r^2 w T with singles
    rate r = N eta / 2 + d; uniform across combinations, so they add nothing to
    the correlator numerator.
    """
    c = pair_rate_per_s * duration_s * efficiency**2 * grid_capture_probability(jitter_sigma_ps, width_ps)
    r = pair_rate_per_s * efficiency / 2.0 + dark_rate_cps
    a = 4.0 * r * r * (width_ps * 1e-12) * duration_s
    if settings_deg is None:
        return chsh_with_accidentals(visibility, c, a)
    return chsh_ideal(visibility, settings_deg) * c / (c + a)


def chsh_ideal(visibility: float, settings_deg: dict[str, float]) -> float:
    """S = sum of sign x V cos 2(alpha - beta) with the fixed signs of `dadl.analysis.chsh.SIGNS`."""
    from dadl.analysis.chsh import SIGNS

    return sum(
        sign * visibility * math.cos(2.0 * math.radians(settings_deg[a] - settings_deg[b])) for (a, b), sign in SIGNS.items()
    )


#: Boltzmann constant in eV/K.
K_B_EV = 8.617333262e-5


def arrhenius_dark_rate(dcr_ref_cps: float, activation_energy_ev: float, temperature_k: float, reference_k: float) -> float:
    """Dark count rate at `temperature_k` for Arrhenius generation: DCR_ref * exp(-Ea/k (1/T - 1/T_ref))."""
    if temperature_k <= 0 or reference_k <= 0:
        raise ValueError("temperatures must be positive (kelvin)")
    return dcr_ref_cps * math.exp(-activation_energy_ev / K_B_EV * (1.0 / temperature_k - 1.0 / reference_k))


def _recorded(
    rate_in: float, dead_time_s: float, afterpulse_prob: float = 0.0, afterpulse_tau_s: float = 0.0
) -> tuple[float, float]:
    """Recorded rate and dead fraction of a non-paralyzable detector fed a Poisson rate.

    With afterpulsing the recorded rate follows `dark_afterpulse_rate` (the input need only be
    Poisson), and the dead fraction is recorded rate x dead time: a photon arriving at a random
    time finds the detector dead with that probability.
    """
    if afterpulse_prob > 0:
        rate = dark_afterpulse_rate(rate_in, afterpulse_prob, afterpulse_tau_s, dead_time_s)
        return rate, rate * dead_time_s
    x = rate_in * dead_time_s
    return rate_in / (1.0 + x), x / (1.0 + x)


def accidentals_expected(
    rate_a_in: float,
    rate_b_in: float,
    dead_time_a_s: float,
    dead_time_b_s: float,
    width_ps: float,
    duration_s: float,
    afterpulse_a: tuple[float, float] = (0.0, 0.0),
    afterpulse_b: tuple[float, float] = (0.0, 0.0),
) -> float:
    """Expected fixed-grid accidental coincidences: R1 R2 width duration with the recorded rates.

    Exact for stationary streams with at most one tag per slot per channel
    (width below the dead time), counted far from the pair peak.
    """
    ra, _ = _recorded(rate_a_in, dead_time_a_s, *afterpulse_a)
    rb, _ = _recorded(rate_b_in, dead_time_b_s, *afterpulse_b)
    return ra * rb * width_ps * 1e-12 * duration_s


def heralding_ratio_expected(
    pair_rate_per_s: float,
    efficiency_a: float,
    efficiency_b: float,
    dark_a_cps: float,
    dark_b_cps: float,
    dead_time_a_s: float,
    dead_time_b_s: float,
    jitter_combined_ps: float,
    width_ps: float,
    afterpulse_a: tuple[float, float] = (0.0, 0.0),
    afterpulse_b: tuple[float, float] = (0.0, 0.0),
) -> float:
    """True heralded detection probability per herald slot: the measurand of `calibrate_heralding`.

    For herald a and signal b, a CW pair source and exact declared delays.

    Pair photons are recorded with probability eta (1 - dead fraction) on each
    arm; both land in one slot with grid_capture_probability. The signal arm's
    live fraction is conditional on a recorded herald: the herald arm was live,
    so no pair with a recorded herald arrived within min(dead times) before it,
    and the signal arm is less often dead. To first order this raises the ratio
    by eta_a eta_b R min(tau_a, tau_b) (the statistician's review measured +1.24%
    at 1 Mpairs/s, 20 ns, eta 0.8, against +1.28% from this term; the
    unconditional form read 0/30 inside the interval). Remaining terms are of
    order (R tau)^2: at 1 Mpairs/s and 20-100 ns the truth sits about 0.0005 above
    the simulated mean, half an interval at 1 s; use it for runs up to about 1 s at
    such rates.
    """
    ra, fa = _recorded(efficiency_a * pair_rate_per_s + dark_a_cps, dead_time_a_s, *afterpulse_a)
    _, fb = _recorded(efficiency_b * pair_rate_per_s + dark_b_cps, dead_time_b_s, *afterpulse_b)
    capture = grid_capture_probability(jitter_combined_ps / math.sqrt(2.0), width_ps)
    conditional = 1.0 + efficiency_a * efficiency_b * pair_rate_per_s * min(dead_time_a_s, dead_time_b_s)
    return pair_rate_per_s * efficiency_a * (1 - fa) * efficiency_b * (1 - fb) * conditional * capture / ra


def hom_split_probability(visibility: float, path_delay_ps: float, coherence_time_ps: float) -> float:
    """Probability that the two photons leave a 50:50 beamsplitter by different ports."""
    return 0.5 * (1.0 - visibility * math.exp(-((path_delay_ps / coherence_time_ps) ** 2)))


def hom_visibility_expected(visibility: float, reference_delay_ps: float, coherence_time_ps: float) -> float:
    """Dip visibility 1 - C(0) / C(reference) without accidentals, dead time or window losses."""
    return 1.0 - hom_split_probability(visibility, 0.0, coherence_time_ps) / hom_split_probability(
        visibility, reference_delay_ps, coherence_time_ps
    )
