"""Dark rate, afterpulsing and dead-time calibration (milestone 1 of phase 2)."""

import math

import numpy as np
import pytest

from dadl.calibration import calibrate_dark
from dadl.calibration.coverage import dark_coverage, gate
from dadl.closed_form import dark_afterpulse_rate
from dadl.data import ChannelId
from dadl.errors import InsufficientCounts
from dadl.sim import DetectorParams, simulate_dark

CH = ChannelId("sim", 0)
TRUTH = DetectorParams(dark_rate_cps=1000.0, dead_time_ps=22_000, afterpulse_prob=0.03, afterpulse_tau_ps=50_000.0)


def dark(det=TRUTH, seconds=20.0, seed=1):
    return simulate_dark(det, seconds, np.random.default_rng(seed))


# --- simulator against closed form (the afterpulse closed-form test) -------------------------


@pytest.mark.parametrize(
    "rate, p, tau_ps",
    [(1000.0, 0.0, 50_000.0), (1000.0, 0.05, 50_000.0), (2e4, 0.2, 200_000.0), (5e5, 0.1, 20_000.0)],
)
def test_simulated_dark_rate_matches_closed_form(rate, p, tau_ps):
    det = DetectorParams(dark_rate_cps=rate, dead_time_ps=22_000, afterpulse_prob=p, afterpulse_tau_ps=tau_ps)
    seconds = 4e5 / rate
    tags = dark(det, seconds, seed=11)
    expected = dark_afterpulse_rate(rate, p, tau_ps * 1e-12, 22e-9)
    assert tags.size / seconds == pytest.approx(expected, rel=0.02)


def test_afterpulses_cascade():
    # With p = 0.5 the share of gaps in the afterpulse region exceeds a single generation's p.
    det = DetectorParams(dark_rate_cps=100.0, dead_time_ps=22_000, afterpulse_prob=0.5, afterpulse_tau_ps=10_000.0)
    tags = dark(det, 200.0, seed=3)
    gaps = np.diff(tags)
    fast = np.mean(gaps < 22_000 + 200_000)
    assert fast == pytest.approx(0.5, abs=0.02)  # renewal: P(afterpulse wins) = p for r * tau << 1
    assert tags.size > 1.8 * 100 * 200  # each primary yields 1 / (1 - p) recorded avalanches


# --- the fit -----------------------------------------------------------------------------------


def test_recovers_injected_parameters():
    cal = calibrate_dark(dark(), CH, backend="sim", seed=2)
    est = cal.estimates
    assert est["dark_rate_cps"].contains(1000.0)
    assert est["afterpulse_prob"].contains(0.03)
    assert est["afterpulse_tau_ps"].contains(50_000.0)
    assert est["dead_time_ps"].contains(22_000.0)
    assert cal.converged
    assert cal.ppc_pvalue > 0.001
    assert cal.kind == "dark-afterpulse-deadtime"
    assert cal.samples.shape == (2000, 4)
    assert cal.sample_params == ("dark_rate_cps", "afterpulse_prob", "afterpulse_tau_ps", "dead_time_ps")
    assert cal.inputs[0][1] == "raw-tags" and len(cal.inputs[0][0]) == 64


def test_samples_match_intervals():
    cal = calibrate_dark(dark(), CH, backend="sim", seed=2)
    for i, name in enumerate(cal.sample_params):
        iv = cal.estimates[name]
        lo, hi = np.quantile(cal.samples[:, i], [0.025, 0.975])
        width = iv.high - iv.low
        assert lo == pytest.approx(iv.low, abs=0.1 * width)
        assert hi == pytest.approx(iv.high, abs=0.1 * width)


def test_seed_makes_draws_reproducible():
    tags = dark()
    a = calibrate_dark(tags, CH, seed=5)
    b = calibrate_dark(tags, CH, seed=5)
    np.testing.assert_array_equal(a.samples, b.samples)


def test_no_afterpulsing_gives_an_upper_limit():
    det = DetectorParams(dark_rate_cps=1000.0, dead_time_ps=22_000, afterpulse_prob=0.0)
    cal = calibrate_dark(dark(det), CH, seed=1)
    assert cal.unidentified == ("afterpulse_tau_ps",)
    assert "afterpulse_tau_ps" not in cal.estimates
    p = cal.estimates["afterpulse_prob"]
    assert p.low == 0.0 and 0.0 < p.high < 0.01
    assert cal.estimates["dark_rate_cps"].contains(1000.0)
    assert any("upper limit" in n for n in cal.notes)
    assert cal.sample_params == ("dark_rate_cps", "dead_time_ps")


def test_misfit_is_flagged_by_the_predictive_check():
    # Two afterpulse time constants (5 ns and 1 us): the one-exponential model must not pass quietly.
    rng = np.random.default_rng(4)
    n, r = 60_000, 1000e-12
    darks = rng.exponential(1 / r, n)
    comp = rng.random(n)
    ap = np.where(comp < 0.08, rng.exponential(5_000, n), np.where(comp < 0.16, rng.exponential(1e6, n), np.inf))
    gaps = 22_000 + np.rint(np.minimum(darks, ap)).astype(np.int64)
    tags = np.cumsum(gaps)
    cal = calibrate_dark(tags, CH, seed=1)
    assert cal.ppc_pvalue < 0.01


def test_too_few_counts_is_refused():
    with pytest.raises(InsufficientCounts):
        calibrate_dark(np.arange(0, 100 * 50_000, 50_000, dtype=np.int64), CH)


def test_float_and_unsorted_tags_are_refused():
    with pytest.raises(TypeError):
        calibrate_dark(np.linspace(0, 1, 500), CH)
    with pytest.raises(ValueError, match="sorted"):
        calibrate_dark(dark()[::-1], CH)


def test_backend_must_be_known():
    with pytest.raises(ValueError):
        calibrate_dark(dark(), CH, backend="lab")


# --- the phase 2 gate ----------------------------------------------------------------------------


def test_gate_is_bonferroni_corrected():
    lo1, hi1 = gate(1000, 1)
    lo4, hi4 = gate(1000, 4)
    assert lo4 < lo1 < 0.95 < hi1 < hi4
    assert math.isclose(lo1, 0.931, abs_tol=0.002) and math.isclose(hi1, 0.967, abs_tol=0.002)


@pytest.mark.coverage
def test_dark_calibration_coverage():
    result = dark_coverage(1000)
    print(result)
    assert result.passed(), str(result)
