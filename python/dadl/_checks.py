"""Input checks shared by the per-tag functions.

The same contract holds in the Rust core (core/src/error.rs) and is pinned by
the shared test vectors in tests/vectors/:

- timestamps are integer picoseconds (int64); floats are refused, because a
  float array in seconds silently truncates to zero;
- inputs documented as sorted are checked, in every build;
- any arithmetic that would leave the int64 range is refused, never wrapped
  or saturated.
"""

from __future__ import annotations

import math
import numbers

import numpy as np

INT64_MIN = -(2**63)
INT64_MAX = 2**63 - 1


def timestamps(times: object, name: str = "times") -> np.ndarray:
    """Return `times` as a 1-D int64 array, refusing float or non-numeric input."""
    arr = np.asarray(times)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if arr.size == 0:
        return arr.astype(np.int64)
    if arr.dtype.kind not in "iu":
        raise TypeError(
            f"{name} must be integer picoseconds (int64), got dtype {arr.dtype}. "
            "Convert seconds with np.rint(t_s * 1e12).astype(np.int64)."
        )
    if arr.dtype.kind == "u" and arr.size and int(arr.max()) > INT64_MAX:
        raise ValueError(f"{name} exceeds the int64 range")
    return arr.astype(np.int64, copy=False)


def check_sorted(times: np.ndarray, name: str = "times") -> None:
    if times.size > 1:
        bad = np.flatnonzero(times[1:] < times[:-1])
        if bad.size:
            raise ValueError(f"{name} must be sorted ascending (first violation at index {int(bad[0]) + 1})")


def check_span(times: np.ndarray, name: str = "times") -> None:
    """Differences between any two tags must fit in int64."""
    if times.size > 1 and int(times.max()) - int(times.min()) > INT64_MAX:
        raise OverflowError(f"the time span of {name} exceeds the int64 range")


def integer(value: object, name: str, *, minimum: int | None = None) -> int:
    """An exact integer (not bool, not a float with a fractional part)."""
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        if isinstance(value, numbers.Real) and math.isfinite(float(value)) and float(value).is_integer():
            value = int(value)
        else:
            raise TypeError(f"{name} must be an integer, got {value!r}")
    value = int(value)
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {value}")
    return value


def finite(value: object, name: str, *, minimum: float | None = None, maximum: float | None = None) -> float:
    """A finite real number (not bool) within optional bounds."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a number, got {value!r}")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value}")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {value}")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} must be <= {maximum}, got {value}")
    return value
