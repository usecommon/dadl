"""Dead-time filters against closed-form Poisson results."""

import numpy as np
import pytest

from dadl.closed_form import (
    nonparalyzable_measured_rate,
    paralyzable_measured_rate,
    patil_usman_hybrid_rate,
    series_measured_rate,
)
from dadl.deadtime import apply_nonparalyzable, apply_paralyzable
from dadl.sim import poisson_times

TRUE_RATE = 1e6  # per second
DURATION_S = 0.1
TAU_S = 500e-9  # n * tau = 0.5, so the two models differ by ~10%


@pytest.fixture(scope="module")
def arrivals():
    rng = np.random.default_rng(7)
    return np.rint(poisson_times(TRUE_RATE, DURATION_S, rng)).astype(np.int64)


def test_nonparalyzable_rate_matches_closed_form(arrivals):
    measured = apply_nonparalyzable(arrivals, int(TAU_S * 1e12)).size / DURATION_S
    expected = nonparalyzable_measured_rate(TRUE_RATE, TAU_S)
    assert measured == pytest.approx(expected, rel=0.015)


def test_paralyzable_rate_matches_closed_form(arrivals):
    measured = apply_paralyzable(arrivals, int(TAU_S * 1e12)).size / DURATION_S
    expected = paralyzable_measured_rate(TRUE_RATE, TAU_S)
    assert measured == pytest.approx(expected, rel=0.015)


def test_models_are_distinguishable(arrivals):
    tau_ps = int(TAU_S * 1e12)
    assert apply_paralyzable(arrivals, tau_ps).size < 0.95 * apply_nonparalyzable(arrivals, tau_ps).size


@pytest.mark.parametrize(
    "rate, tau_p_s, tau_np_s",
    [(1e7, 10e-9, 40e-9), (2e7, 22e-9, 22e-9), (4e7, 5e-9, 30e-9)],
)
def test_series_chain_matches_exact_closed_form(rate, tau_p_s, tau_np_s):
    rng = np.random.default_rng(11)
    duration = 2e-3
    t = np.rint(poisson_times(rate, duration, rng)).astype(np.int64)
    chain = apply_nonparalyzable(apply_paralyzable(t, round(tau_p_s * 1e12)), round(tau_np_s * 1e12))
    measured = chain.size / duration
    assert measured == pytest.approx(series_measured_rate(rate, tau_p_s, tau_np_s), rel=0.02)


def test_patil_usman_hybrid_is_not_the_series_model():
    # Documented as empirical: at 10 MHz, 10 ns / 40 ns it is ~9 % below the exact series rate.
    exact = series_measured_rate(1e7, 10e-9, 40e-9)
    assert patil_usman_hybrid_rate(1e7, 10e-9, 40e-9) < 0.93 * exact


def test_unsorted_input_refused():
    with pytest.raises(ValueError, match="sorted"):
        apply_nonparalyzable(np.array([30, 0, 5, 10], dtype=np.int64), 10)


def test_float_timestamps_refused():
    with pytest.raises(TypeError, match="picoseconds"):
        apply_paralyzable(np.array([0.0, 1e-6]), 10)


def test_negative_dead_time_rejected():
    with pytest.raises(ValueError):
        apply_paralyzable(np.array([0, 1], dtype=np.int64), -1)
