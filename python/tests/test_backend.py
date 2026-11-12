"""The optional compiled core agrees with the Python reference beyond the shared vectors."""

import numpy as np
import pytest

from dadl import _backend, backend
from dadl.coincidences import grid_coincidences, occupied_slots
from dadl.deadtime import apply_nonparalyzable, apply_paralyzable
from dadl.histogram import cross_correlation, interarrival

needs_core = pytest.mark.skipif("rust" not in _backend.available(), reason="dadl-quantum-core not installed")


def test_backend_name():
    assert backend() in _backend.available()
    with _backend.use("python"):
        assert backend() == "python" and _backend.core() is None
    with pytest.raises(ValueError):
        with _backend.use("fortran"):
            pass


@needs_core
@pytest.mark.parametrize("seed", range(5))
def test_random_streams_agree(seed):
    rng = np.random.default_rng(seed)
    a = np.sort(rng.integers(-(10**9), 10**9, 5000))
    b = np.sort(rng.integers(-(10**9), 10**9, 5000))
    calls = [
        lambda: apply_nonparalyzable(a, 300_000),
        lambda: apply_paralyzable(a, 300_000),
        lambda: occupied_slots(a, 1000, 37),
        lambda: grid_coincidences(a, b, 100_000, 5, -7),
        lambda: interarrival(a, 1000, 2000),
        lambda: cross_correlation(a, b, 997, 1_000_000),
    ]
    for call in calls:
        with _backend.use("python"):
            ref = call()
        with _backend.use("rust"):
            got = call()
        np.testing.assert_array_equal(np.asarray(got), np.asarray(ref))


@needs_core
def test_non_contiguous_and_unsigned_inputs():
    a = np.arange(0, 2000, dtype=np.int64)
    with _backend.use("rust"):
        np.testing.assert_array_equal(apply_nonparalyzable(a[::2], 5), a[::2][::3])
        assert interarrival(np.arange(0, 50, 10, dtype=np.uint32), 10, 3).tolist() == [0, 4, 0]
        with pytest.raises(TypeError):
            apply_nonparalyzable(np.linspace(0, 1, 4), 1)


@needs_core
def test_errors_match():
    with _backend.use("rust"):
        with pytest.raises(ValueError, match="sorted"):
            cross_correlation(np.array([3, 1]), np.array([0]), 1, 5)
        with pytest.raises(OverflowError):
            occupied_slots(np.array([-(2**63) + 5]), 1000, 1000)
        with pytest.raises(ValueError):
            interarrival(np.array([1, 2]), 0, 3)
