"""CHSH from the simulator against closed-form expectations.

Werner-like Phi+ state with Poisson accidentals: S = 2 sqrt(2) V C / (C + A).
"""

import math

import numpy as np
import pytest

from dadl.analysis.chsh import FAIR_SAMPLING, Window, chsh
from dadl.closed_form import chsh_expected_grid, chsh_werner, grid_capture_probability
from dadl.sim import CHANNELS, DetectorParams, PairSource, simulate_chsh


def run(dark_rate: float, seed: int):
    params = dict(pair_rate=1e5, duration_s=0.5, efficiency=0.5, jitter_sigma_ps=50.0, width_ps=10_000, visibility=0.95)
    det = DetectorParams(
        efficiency=params["efficiency"],
        dark_rate=dark_rate,
        dead_time_ps=0,
        jitter_ps=params["jitter_sigma_ps"],
    )
    runs = simulate_chsh(
        PairSource(pair_rate=params["pair_rate"], visibility=params["visibility"]),
        {ch: det for ch in CHANNELS},
        params["duration_s"],
        np.random.default_rng(seed),
    )
    result = chsh(runs, Window(width_ps=params["width_ps"]))
    expected = chsh_expected_grid(dark_rate=dark_rate, **params)
    return result, expected


def test_capture_probability_limits():
    assert grid_capture_probability(0.0, 1000) == 1.0
    assert grid_capture_probability(50.0, 10_000) == pytest.approx(1 - 70.71 * math.sqrt(2 / math.pi) / 10_000, rel=1e-3)
    assert 0.0 < grid_capture_probability(350.0, 2_000) < 1.0


def test_negligible_accidentals_gives_werner_value():
    result, expected = run(dark_rate=0.0, seed=1)
    assert expected == pytest.approx(chsh_werner(0.95), rel=0.01)
    assert result.raw.value == pytest.approx(expected, abs=5 * result.raw.stderr + 0.01)
    assert FAIR_SAMPLING in result.labels


def test_accidentals_reduce_s_as_predicted():
    # Dark counts chosen so accidentals are ~30% of true coincidences: S falls from ~2.69 to ~2.08.
    result, expected = run(dark_rate=4e5, seed=2)
    assert 1.9 < expected < 2.2
    assert result.raw.value == pytest.approx(expected, abs=5 * result.raw.stderr + 0.01)
    assert result.corrected is None, "Bell headline numbers are raw only"
