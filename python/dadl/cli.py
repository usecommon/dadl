"""Command-line interface: `dadl quickstart <example>` and `dadl version`."""

from __future__ import annotations

import argparse
import sys

import numpy as np

from dadl import __version__
from dadl.analysis.chsh import Window, chsh
from dadl.closed_form import chsh_expected_grid
from dadl.sim import CHANNELS, DetectorParams, PairSource, simulate_chsh

#: Phase 1 keeps the example parameters in code; examples/chsh-spad/setup.yaml mirrors them.
CHSH_SPAD = {
    "seed": 2026,
    "duration_s": 1.0,
    "pair_rate": 1e5,
    "visibility": 0.95,
    "detector": {
        "efficiency": 0.5,
        "dark_rate": 500.0,
        "dead_time_ps": 22_000,
        "jitter_ps": 350.0,
        "afterpulse_prob": 0.005,
        "afterpulse_tau_ps": 50_000.0,
    },
    "window_ps": 2_000,
}


def quickstart_chsh_spad() -> int:
    p = CHSH_SPAD
    det = DetectorParams(**p["detector"])
    rng = np.random.default_rng(p["seed"])
    runs = simulate_chsh(
        PairSource(pair_rate=p["pair_rate"], visibility=p["visibility"]),
        {ch: det for ch in CHANNELS},
        p["duration_s"],
        rng,
    )
    result = chsh(runs, Window(width_ps=p["window_ps"], source="plan"))
    expected = chsh_expected_grid(
        pair_rate=p["pair_rate"],
        duration_s=p["duration_s"],
        efficiency=det.efficiency,
        dark_rate=det.dark_rate,
        jitter_sigma_ps=det.jitter_ps,
        width_ps=p["window_ps"],
        visibility=p["visibility"],
    )
    print(f"DADL {__version__} quickstart: chsh-spad (simulated data, seed {p['seed']})")
    print(f"  Source: {p['pair_rate']:.0e} pairs/s, visibility {p['visibility']}, {p['duration_s']} s per setting")
    print(
        f"  Detectors: eta {det.efficiency}, dark {det.dark_rate:.0f} cps, dead time {det.dead_time_ps / 1000:.0f} ns, "
        f"jitter {det.jitter_ps:.0f} ps, afterpulsing {det.afterpulse_prob:.1%}"
    )
    print(f"  Fixed window from the analysis plan: {p['window_ps'] / 1000:.1f} ns")
    print()
    print(f"  Raw S = {result.raw}")
    print(f"  Labels: {'; '.join(result.labels)}")
    print(f"  Closed-form expectation (ignores dead time and afterpulsing): S = {expected:.4f}")
    print()
    print("  Next: calibration fits and corrected values arrive in phase 2. See docs/design.md.")
    return 0


EXAMPLES = {"chsh-spad": quickstart_chsh_spad}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="dadl", description="Detector-Aware Data Layer for quantum optics experiments")
    sub = parser.add_subparsers(dest="command", required=True)
    qs = sub.add_parser("quickstart", help="run an example end to end on simulated data")
    qs.add_argument("example", choices=sorted(EXAMPLES))
    sub.add_parser("version", help="print the version")
    args = parser.parse_args(argv)
    if args.command == "version":
        print(__version__)
        return 0
    return EXAMPLES[args.example]()


if __name__ == "__main__":
    sys.exit(main())
