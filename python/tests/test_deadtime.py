"""Dead-time filters against closed-form Poisson results."""

import numpy as np
import pytest

from dadl.closed_form import nonparalyzable_measured_rate, paralyzable_measured_rate
from dadl.deadtime import apply_nonparalyzable, apply_paralyzable
from dadl.sim import poisson_times

TRUE_RATE = 1e6  # per second
DURATION_S = 0.1
TAU_S = 500e-9  # n * tau = 0.5, so the two models differ by ~10%


@pytest.fixture(scope="module")
def arrivals():
    rng = np.random.default_rng(7)
    return np.rint(poisson_times(TRUE_RATE, DURATION_S, rng)).astype(np.int64)


def test_small_examples_match_rust_core():
    times = np.array([0, 5, 10, 25, 30], dtype=np.int64)
    assert apply_nonparalyzable(times, 10).tolist() == [0, 10, 25]
    assert apply_paralyzable(times, 10).tolist() == [0, 25]


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


def test_negative_dead_time_rejected():
    with pytest.raises(ValueError):
        apply_nonparalyzable(np.array([0, 1], dtype=np.int64), -1)
