"""Fixed-grid coincidences: same cases as the Rust core tests."""

import numpy as np

from dadl.coincidences import grid_coincidences


def a(*xs):
    return np.array(xs, dtype=np.int64)


def test_counts_shared_slots():
    assert grid_coincidences(a(0, 1000, 2500), a(100, 1900, 9000), 1000) == 2


def test_declared_delay_shifts_slots():
    assert grid_coincidences(a(50), a(950), 1000) == 1
    assert grid_coincidences(a(50), a(1050), 1000) == 0
    assert grid_coincidences(a(50), a(1050), 1000, delay_b_ps=1000) == 1


def test_several_events_in_one_slot_count_once():
    assert grid_coincidences(a(10, 20, 30), a(40), 1000) == 1


def test_negative_times_use_floor_division():
    assert grid_coincidences(a(-1), a(-999), 1000) == 1
    assert grid_coincidences(a(-1), a(0), 1000) == 0
