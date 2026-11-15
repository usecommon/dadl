"""Crosstalk in the simulator: closed-form counts and the dead-time interaction.

Each test uses a setting where the expected count is exact or a binomial / geometric
closed form, with tolerances of about 4 standard deviations.
"""

import math

import numpy as np
import pytest

from dadl.sim import Crosstalk, DetectorParams, PairSource, detect, detect_channels, poisson_times, simulate_hom

IDEAL = DetectorParams(efficiency=1.0, dark_rate_cps=0.0, dead_time_ps=1_000, jitter_sigma_ps=0.0)
DURATION_S = 0.5
RATE = 1e5  # mean gap 10 us: dead times of ns block about 1e-4 of events


def _offsets(a: np.ndarray, b: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """All b - a in [lo, hi)."""
    out = []
    for shift in range(-3, 4):
        idx = np.searchsorted(a, b - hi) + shift
        ok = (idx >= 0) & (idx < a.size)
        d = b[ok] - a[idx[ok]]
        out.append(d[(d >= lo) & (d < hi)])
    return np.concatenate(out)


def _spaced(n: int, rng: np.random.Generator) -> np.ndarray:
    """n photon times at least 200 ns apart: no dead time reaches from one to the next."""
    return np.cumsum(rng.exponential(1e7, n) + 200_000.0)


def _within(observed: float, expected: float, sd: float, k: float = 4.0) -> None:
    assert abs(observed - expected) <= k * sd, (observed, expected, sd)


def test_without_crosstalk_results_are_unchanged():
    photons = {"A": poisson_times(RATE, 0.1, np.random.default_rng(0)), "B": poisson_times(RATE, 0.1, np.random.default_rng(1))}
    det = DetectorParams(afterpulse_prob=0.02)
    got = detect_channels(photons, {"A": det, "B": det}, 0.1, np.random.default_rng(7))
    rng = np.random.default_rng(7)
    want = {ch: detect(photons[ch], det, 0.1, rng) for ch in photons}
    for ch in photons:
        np.testing.assert_array_equal(got[ch], want[ch])
    src = PairSource(pair_rate_per_s=5e4, visibility=0.9, coherence_time_ps=2.0)
    a = simulate_hom(src, det, det, 0.0, 0.05, np.random.default_rng(3))
    b = simulate_hom(src, det, det, 0.0, 0.05, np.random.default_rng(3), crosstalk=())
    for ch in a:
        np.testing.assert_array_equal(a[ch], b[ch])


def test_crosstalk_on_a_dark_channel_is_a_thinned_shifted_copy():
    rng = np.random.default_rng(11)
    p, delay = 0.05, 2_500
    photons = {"A": poisson_times(RATE, DURATION_S, rng), "B": np.empty(0)}
    out = detect_channels(photons, {"A": IDEAL, "B": IDEAL}, DURATION_S, rng,
                          [Crosstalk("A", "B", p, delay_ps=delay)])
    a, b = out["A"], out["B"]
    # Every crosstalk tag sits exactly `delay` after a recorded A tag.
    assert np.isin(b - delay, a).all()
    _within(b.size, p * a.size, math.sqrt(a.size * p * (1 - p)))


def test_cross_correlation_excess_is_p_times_source_singles():
    """Two independent Poisson channels: the B - A histogram has a crosstalk peak of p N_A over a flat floor."""
    rng = np.random.default_rng(12)
    p, delay, sigma = 0.02, 3_000, 150.0
    photons = {"A": poisson_times(RATE, DURATION_S, rng), "B": poisson_times(RATE, DURATION_S, rng)}
    out = detect_channels(photons, {"A": IDEAL, "B": IDEAL}, DURATION_S, rng,
                          [Crosstalk("A", "B", p, delay_ps=delay, jitter_sigma_ps=sigma)])
    a, b = out["A"], out["B"]
    half = int(5 * sigma)
    peak = _offsets(a, b, delay - half, delay + half).size
    # Off-peak floor from the same histogram, 20 windows of the same width far from the peak.
    floor = np.mean([_offsets(a, b, lo, lo + 2 * half).size for lo in range(50_000, 50_000 + 40 * half, 2 * half)])
    # Crosstalk lost in B's dead time (another B tag within 1 ns before) is about RATE x dead.
    expected = p * a.size * (1 - RATE * IDEAL.dead_time_ps * 1e-12)
    _within(peak - floor, expected, math.sqrt(expected + 2 * floor + 1))
    # Nothing at negative offsets: crosstalk never precedes its cause.
    mirror = _offsets(a, b, -delay - half, -delay + half).size
    _within(mirror, floor, math.sqrt(floor + 1))


def test_crosstalk_cascades_both_ways():
    """A <-> B with probability p each way: geometric chains, A gets N p^2 / (1 - p^2) echoes at 2 x delay."""
    rng = np.random.default_rng(13)
    p, delay = 0.3, 2_000
    photons = {"A": poisson_times(RATE, DURATION_S, rng), "B": np.empty(0)}
    links = [Crosstalk("A", "B", p, delay_ps=delay), Crosstalk("B", "A", p, delay_ps=delay)]
    out = detect_channels(photons, {"A": IDEAL, "B": IDEAL}, DURATION_S, rng, links)
    n = photons["A"].size
    a, b = out["A"], out["B"]
    echoes_a = a.size - n  # A's own photons are all recorded at this rate (up to ~1e-4)
    q = p * p
    # Each primary starts an independent chain; per chain, A echoes are geometric with mean q / (1 - q).
    _within(echoes_a, n * q / (1 - q), math.sqrt(n * q / (1 - q) ** 2) + 0.002 * n)
    _within(b.size, n * p / (1 - q), math.sqrt(n * p * (1 + q) / (1 - q) ** 2) + 0.002 * n)
    # Echoes sit at exactly 2 x delay after an A tag, and B tags at delay after one.
    assert np.isin(b - delay, a).all()
    assert _offsets(a, a, 2 * delay, 2 * delay + 1).size >= echoes_a - 0.002 * n


def test_crosstalk_inside_the_target_dead_time_is_lost():
    rng = np.random.default_rng(14)
    t = _spaced(20_000, rng)
    det_b = DetectorParams(efficiency=1.0, dark_rate_cps=0.0, dead_time_ps=50_000, jitter_sigma_ps=0.0)
    photons = {"A": t, "B": t}
    # B records its own photon at the same time as A's: crosstalk 10 ns later falls in B's dead time.
    blocked = detect_channels(photons, {"A": IDEAL, "B": det_b}, DURATION_S, rng,
                              [Crosstalk("A", "B", 0.5, delay_ps=10_000)])
    np.testing.assert_array_equal(blocked["B"], blocked["A"])
    # 60 ns later it is recorded: B gains p N tags at the crosstalk offset.
    seen = detect_channels(photons, {"A": IDEAL, "B": det_b}, DURATION_S, rng,
                           [Crosstalk("A", "B", 0.5, delay_ps=60_000)])
    n = seen["A"].size
    extra = seen["B"].size - n
    assert np.isin(seen["B"], np.concatenate([seen["A"], seen["A"] + 60_000])).all()
    _within(extra, 0.5 * n, math.sqrt(0.25 * n))


def test_recorded_crosstalk_starts_a_dead_time():
    """B's own photon arrives inside the dead time of a crosstalk avalanche: exactly one of the two is recorded."""
    rng = np.random.default_rng(15)
    p, delay, dead = 0.3, 5_000, 40_000
    t = _spaced(20_000, rng)
    det_b = DetectorParams(efficiency=1.0, dark_rate_cps=0.0, dead_time_ps=dead, jitter_sigma_ps=0.0,
                           delay_ps=delay + dead // 2)
    out = detect_channels({"A": t, "B": t}, {"A": IDEAL, "B": det_b}, DURATION_S, rng,
                          [Crosstalk("A", "B", p, delay_ps=delay)])
    a, b = out["A"], out["B"]
    at_crosstalk = np.isin(b - delay, a).sum()
    at_photon = np.isin(b - delay - dead // 2, a).sum()
    assert at_crosstalk + at_photon == a.size == b.size
    _within(at_crosstalk, p * a.size, math.sqrt(a.size * p * (1 - p)))


def test_crosstalk_avalanches_afterpulse():
    """A crosstalk avalanche is an avalanche: B's afterpulses scale with its crosstalk count."""
    rng = np.random.default_rng(16)
    p_ct, p_ap = 0.2, 0.25
    det_b = DetectorParams(efficiency=1.0, dark_rate_cps=0.0, dead_time_ps=1_000, jitter_sigma_ps=0.0,
                           afterpulse_prob=p_ap, afterpulse_tau_ps=2_000.0)
    photons = {"A": poisson_times(RATE, DURATION_S, rng), "B": np.empty(0)}
    out = detect_channels(photons, {"A": IDEAL, "B": det_b}, DURATION_S, rng, [Crosstalk("A", "B", p_ct)])
    n_ct = photons["A"].size * p_ct
    # Each crosstalk avalanche starts a geometric afterpulse chain: B = N_ct / (1 - p_ap).
    expected = n_ct / (1 - p_ap)
    _within(out["B"].size, expected, math.sqrt(expected * (1 + p_ap) / (1 - p_ap)) + 0.01 * expected)


def test_bad_links_are_refused():
    with pytest.raises(ValueError, match="two different channels"):
        Crosstalk("A", "A", 0.1)
    with pytest.raises(ValueError, match=r"\[0, 1\)"):
        Crosstalk("A", "B", 1.0)
    with pytest.raises(ValueError):
        Crosstalk("A", "B", 0.1, delay_ps=-5)
    with pytest.raises(ValueError, match="crosstalk names"):
        detect_channels({"A": np.empty(0)}, {"A": IDEAL}, 0.1, np.random.default_rng(0), [Crosstalk("A", "C", 0.1)])
