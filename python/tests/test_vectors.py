"""Shared vectors in tests/vectors/core.json; the Rust core runs the same cases."""

import json
from pathlib import Path

import numpy as np
import pytest

from dadl.coincidences import grid_coincidences
from dadl.deadtime import apply_nonparalyzable, apply_paralyzable
from dadl.histogram import cross_correlation, interarrival

VECTORS = Path(__file__).resolve().parents[2] / "tests" / "vectors" / "core.json"
CASES = json.loads(VECTORS.read_text(encoding="utf-8"))["cases"] if VECTORS.exists() else []

ERRORS = {
    "Unsorted": (ValueError, "sorted"),
    "NegativeDeadTime": (ValueError, None),
    "NonPositiveWidth": (ValueError, None),
    "Overflow": (OverflowError, None),
}


def run(case):
    arr = lambda xs: np.array(xs, dtype=np.int64)  # noqa: E731
    op = case["op"]
    if op == "nonparalyzable":
        return apply_nonparalyzable(arr(case["times"]), case["tau"]).tolist()
    if op == "paralyzable":
        return apply_paralyzable(arr(case["times"]), case["tau"]).tolist()
    if op == "grid_coincidences":
        return grid_coincidences(arr(case["a"]), arr(case["b"]), case["width"], case["delay_a"], case["delay_b"])
    if op == "interarrival":
        return interarrival(arr(case["times"]), case["bin_width"], case["n_bins"]).tolist()
    if op == "cross_correlation":
        return cross_correlation(arr(case["a"]), arr(case["b"]), case["bin_width"], case["range"]).tolist()
    raise AssertionError(f"unknown op {op}")


@pytest.mark.skipif(not CASES, reason="shared vectors not present (e.g. running from the sdist)")
@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_shared_vector(case):
    expect = case["expect"]
    if isinstance(expect, dict):
        exc, match = ERRORS[expect["error"]]
        with pytest.raises(exc, match=match):
            run(case)
    else:
        assert run(case) == expect
