"""Interval coverage of the calibrations on simulated data (the phase 2 gate).

Each replicate simulates one acquisition from known parameters with an
independent random stream, runs the calibration, and records whether each 95%
interval contains the truth. Run from the command line:

    python -m dadl.calibration.coverage dark --replicates 1000
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass

import numpy as np

from dadl.calibration.dark import calibrate_dark
from dadl.data import ChannelId
from dadl.sim import DetectorParams, simulate_dark
from dadl.stats import coverage_band

#: Simulated detector for the dark calibration gate: a Si SPAD module with
#: visible afterpulsing, 20 s per replicate (about 20,000 dark counts).
DARK_TRUTH = DetectorParams(
    efficiency=0.5, dark_rate_cps=1000.0, dead_time_ps=22_000, afterpulse_prob=0.03, afterpulse_tau_ps=50_000.0
)
DARK_DURATION_S = 20.0

#: Chance that a correctly calibrated method fails the gate, per calibration.
#: Split across the calibration's parameters (Bonferroni), so a calibration
#: with four parameters does not fail four times as often as one with one.
FALSE_FAIL = 0.01


@dataclass(frozen=True)
class CoverageResult:
    calibration: str
    replicates: int
    covered: dict[str, int]
    gate: tuple[float, float]

    def fraction(self, name: str) -> float:
        return self.covered[name] / self.replicates

    def passed(self) -> bool:
        lo, hi = self.gate
        return all(lo <= self.fraction(k) <= hi for k in self.covered)

    def __str__(self) -> str:
        lo, hi = self.gate
        lines = [f"{self.calibration}: {self.replicates} replicates, gate {lo:.1%}-{hi:.1%}"]
        for k in self.covered:
            f = self.fraction(k)
            lines.append(f"  {k:<20} {f:.1%}  {'ok' if lo <= f <= hi else 'FAIL'}")
        return "\n".join(lines)


def gate(replicates: int, n_params: int) -> tuple[float, float]:
    """Accepted coverage range: binomial band at confidence 1 - FALSE_FAIL / n_params."""
    return coverage_band(replicates, confidence=1 - FALSE_FAIL / n_params)


def dark_coverage(replicates: int = 1000, seed: int = 20261111) -> CoverageResult:
    truth = {
        "dark_rate_cps": DARK_TRUTH.dark_rate_cps,
        "afterpulse_prob": DARK_TRUTH.afterpulse_prob,
        "afterpulse_tau_ps": DARK_TRUTH.afterpulse_tau_ps,
        "dead_time_ps": float(DARK_TRUTH.dead_time_ps),
    }
    covered = dict.fromkeys(truth, 0)
    channel = ChannelId("sim", 0)
    for child in np.random.SeedSequence(seed).spawn(replicates):
        sim_rng, fit_seed = child.spawn(2)
        tags = simulate_dark(DARK_TRUTH, DARK_DURATION_S, np.random.default_rng(sim_rng))
        cal = calibrate_dark(tags, channel, backend="sim", seed=fit_seed.generate_state(1)[0], n_draws=10)
        for k, v in truth.items():
            iv = cal.estimates.get(k)
            if iv is not None and iv.contains(v):
                covered[k] += 1
    return CoverageResult("dark-afterpulse-deadtime", replicates, covered, gate(replicates, len(truth)))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m dadl.calibration.coverage")
    parser.add_argument("calibration", choices=["dark"])
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261111)
    args = parser.parse_args(argv)
    result = dark_coverage(args.replicates, args.seed)
    print(result)
    return 0 if result.passed() else 1


if __name__ == "__main__":
    sys.exit(main())
