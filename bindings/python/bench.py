"""Throughput of the per-tag functions on both backends (recorded, not gated).

    python bindings/python/bench.py [n_tags]

Needs dadl-quantum and, for the Rust column, dadl-quantum-core.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from dadl import _backend
from dadl.coincidences import grid_coincidences
from dadl.deadtime import apply_nonparalyzable
from dadl.histogram import cross_correlation, interarrival


def bench(fn, repeat: int = 3) -> float:
    best = float("inf")
    for _ in range(repeat):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    return best


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    rng = np.random.default_rng(0)
    span = n * 10**6  # 1 Mcps
    a = np.sort(rng.integers(0, span, n))
    b = np.sort(rng.integers(0, span, n))
    cases = {
        "apply_nonparalyzable (22 ns)": lambda: apply_nonparalyzable(a, 22_000),
        "grid_coincidences (2 ns)": lambda: grid_coincidences(a, b, 2_000),
        "interarrival (1 ns bins)": lambda: interarrival(a, 1_000, 10_000),
        "cross_correlation (+-20 ns)": lambda: cross_correlation(a, b, 20, 20_000),
    }
    print(f"{n:,} tags per channel, best of 3; Mtags/s")
    print(f"{'function':<32}{'python':>10}{'rust':>10}{'speed-up':>10}")
    for label, fn in cases.items():
        with _backend.use("python"):
            tp = bench(fn)
        row = f"{label:<32}{n / tp / 1e6:>10.1f}"
        if "rust" in _backend.available():
            with _backend.use("rust"):
                tr = bench(fn)
            row += f"{n / tr / 1e6:>10.1f}{tp / tr:>9.1f}x"
        print(row)


if __name__ == "__main__":
    main()
