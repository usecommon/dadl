"""CHSH from the simulator against closed-form expectations, and error-bar calibration.

Werner-like Phi+ state with Poisson accidentals: S = 2 sqrt(2) V C / (C + A).
"""

import math

import numpy as np
import pytest

from dadl.analysis.chsh import FAIR_SAMPLING, Window, chsh
from dadl.closed_form import chsh_expected_grid, chsh_werner, grid_capture_probability
from dadl.errors import InsufficientCounts
from dadl.sim import CHANNELS, DetectorParams, PairSource, simulate_chsh

PARAMS = dict(pair_rate_per_s=1e5, duration_s=0.5, efficiency=0.5, jitter_sigma_ps=50.0, width_ps=10_000, visibility=0.95)


def run(dark_rate_cps: float, seed: int, **overrides):
    p = {**PARAMS, **overrides}
    det = DetectorParams(
        efficiency=p["efficiency"], dark_rate_cps=dark_rate_cps, dead_time_ps=0, jitter_sigma_ps=p["jitter_sigma_ps"]
    )
    runs = simulate_chsh(
        PairSource(pair_rate_per_s=p["pair_rate_per_s"], visibility=p["visibility"]),
        {ch: det for ch in CHANNELS},
        p["duration_s"],
        np.random.default_rng(seed),
    )
    result = chsh(runs, Window(width_ps=p["width_ps"]))
    expected = chsh_expected_grid(dark_rate_cps=dark_rate_cps, **p)
    return result, expected


def test_capture_probability_limits():
    assert grid_capture_probability(0.0, 1000) == 1.0
    assert grid_capture_probability(50.0, 10_000) == pytest.approx(1 - 70.71 * math.sqrt(2 / math.pi) / 10_000, rel=1e-3)
    assert 0.0 < grid_capture_probability(350.0, 2_000) < 1.0
    with pytest.raises(ValueError):
        grid_capture_probability(-50.0, 1000)


def test_negligible_accidentals_gives_werner_value():
    result, expected = run(dark_rate_cps=0.0, seed=1)
    assert expected == pytest.approx(chsh_werner(0.95), rel=0.01)
    assert result.raw.value == pytest.approx(expected, abs=5 * result.raw.stderr + 0.01)
    assert FAIR_SAMPLING in result.labels
    assert result.window == {"width_ps": 10_000, "delays_ps": {}, "source": "plan"}


def test_accidentals_reduce_s_as_predicted():
    # Dark counts chosen so accidentals are ~30% of true coincidences: S falls from ~2.69 to ~2.08.
    result, expected = run(dark_rate_cps=4e5, seed=2)
    assert 1.9 < expected < 2.2
    assert result.raw.value == pytest.approx(expected, abs=5 * result.raw.stderr + 0.01)
    assert result.corrected is None, "Bell headline numbers are raw only"


def test_missing_channel_named():
    result_runs = {k: {"A+": np.array([0])} for k in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    with pytest.raises(ValueError, match="A-"):
        chsh(result_runs, Window(width_ps=1000))


def test_too_few_counts_hints_at_units():
    empty = {ch: np.array([1, 2, 3], dtype=np.int64) for ch in CHANNELS}
    runs = {k: empty for k in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    with pytest.raises(InsufficientCounts, match="picoseconds"):
        chsh(runs, Window(width_ps=1000))


def test_float_tags_refused():
    floats = {ch: np.array([0.1, 0.2]) for ch in CHANNELS}
    runs = {k: floats for k in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    with pytest.raises(TypeError, match="picoseconds"):
        chsh(runs, Window(width_ps=1000))


@pytest.mark.parametrize("width", [0, -5, 0.5])
def test_window_width_validated(width):
    with pytest.raises((TypeError, ValueError)):
        Window(width_ps=width)


@pytest.mark.slow
def test_stderr_is_calibrated_by_pull_study():
    """Over 200 seeds the reported stderr must match the scatter of S."""
    pulls = []
    for seed in range(200):
        result, expected = run(dark_rate_cps=2e4, seed=1000 + seed, duration_s=0.02, pair_rate_per_s=5e4)
        pulls.append((result.raw.value - expected) / result.raw.stderr)
    pulls = np.array(pulls)
    assert 0.85 <= pulls.std(ddof=1) <= 1.15
    assert abs(pulls.mean()) < 0.25
