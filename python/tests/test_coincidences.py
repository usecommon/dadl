"""Fixed-grid coincidences (ordinary cases; edge cases are in tests/vectors/core.json)."""

import numpy as np
import pytest

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


def test_empty_inputs():
    assert grid_coincidences(a(), a(5), 1000) == 0


@pytest.mark.parametrize("width", [0.5, 1000.9, -5, 0])
def test_bad_widths_refused(width):
    with pytest.raises((TypeError, ValueError)):
        grid_coincidences(a(0), a(0), width)


def test_float_times_refused():
    with pytest.raises(TypeError):
        grid_coincidences(np.array([999.9]), a(0), 1000)


def test_delay_overflow_refused():
    with pytest.raises(OverflowError):
        grid_coincidences(a(2**63 - 6), a(0), 1000, delay_a_ps=-100)


def test_occupied_slots_equal_np_unique_on_sorted_and_unsorted_input():
    # Council round 9 (Pf2): sorted input is de-duplicated in one pass; the result is np.unique's.
    from dadl.coincidences import as_sorted, occupied_slots

    rng = np.random.default_rng(1)
    t = np.sort(rng.integers(0, 10**9, 50_000)).astype(np.int64)
    for x in (t, t[::-1].copy(), rng.permutation(t), np.repeat(t[:10], 3)):
        assert np.array_equal(occupied_slots(x, 1_000, 7), np.unique(np.floor_divide(x - 7, 1_000)))
    assert as_sorted(t) is t and np.array_equal(as_sorted(t[::-1]), t)
