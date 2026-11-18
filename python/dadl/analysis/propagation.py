"""The `propagation` gate (phase 3): coverage of corrected values with calibration uncertainty.

    python -m dadl.analysis.propagation hom --replicates 1000 --jobs 4

Each replicate simulates its own calibration runs (dark acquisitions per channel), calibrates
them, writes and loads the artifacts as an analysis would, simulates the analysed acquisitions,
and checks whether the corrected value's 95% interval (statistical plus calibration, by joint
Monte Carlo) covers the injected truth. The accepted range is the binomial band of the phase 2
coverage gate (`dadl.calibration.coverage.gate`).

A control runs alongside: the same corrected value with the calibration uncertainty switched off
(point +- 1.96 statistical standard deviations). The gate (M3) requires the control to fall below
the band in at least one declared regime, or the regimes do not test propagation; it is reported
here, not gated per regime.

M2 ships the HOM regimes, run in the default test suite on a fast tier (200 replicates) as a
smoke test; M3 adds the CHSH regimes, the misspecified regime and the CI job.
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import tempfile
from dataclasses import dataclass

import numpy as np

from dadl.analysis.chsh import Window
from dadl.analysis.hom import hom_from_tags
from dadl.artifact import write_artifact
from dadl.calibration.coverage import CoverageResult, _children, gate
from dadl.calibration.dark import calibrate_dark
from dadl.closed_form import hom_visibility_expected
from dadl.data import ChannelId
from dadl.guardrails import Policy
from dadl.propagate import CalibrationSet, DatasetContext
from dadl.sim import DetectorParams, PairSource, simulate_dark, simulate_hom

UTC = dt.timezone.utc
CHANNELS = (ChannelId("sim", 0), ChannelId("sim", 1))
POLICY = Policy.default().allow("hom_visibility", "dead_time", "propagation gate: dead time from simulated dark runs")


@dataclass(frozen=True)
class HomRegime:
    source: PairSource
    detector: DetectorParams
    window_ps: int
    reference_delay_ps: float
    duration_s: float
    dark_duration_s: float
    #: Whether the dark calibration can identify the dead time here; if not (an SNSPD's 100 cps),
    #: the dark artifacts are not loaded and only the accidentals are corrected.
    dead_time_from_dark: bool

    @property
    def truth(self) -> float:
        """The two-point visibility the forward model targets: the source's, at the declared delays."""
        return hom_visibility_expected(self.source.visibility, self.reference_delay_ps, self.source.coherence_time_ps)


#: Declared before the first run (phase 3 plan, M3).
HOM_REGIMES: dict[str, HomRegime] = {
    # examples/hom-snspd: dark counts too few to identify the dead time; accidentals only.
    "hom-snspd": HomRegime(
        PairSource(pair_rate_per_s=2.0e5, visibility=0.92, coherence_time_ps=1.5),
        DetectorParams(efficiency=0.85, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=25.0),
        window_ps=200, reference_delay_ps=15.0, duration_s=0.5, dark_duration_s=2.0, dead_time_from_dark=False,
    ),
    # Bright background on SPADs: accidentals and dead time both matter, dead time from dark runs.
    "hom-bright": HomRegime(
        PairSource(pair_rate_per_s=2.0e6, visibility=0.92, coherence_time_ps=2.0),
        DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=60_000, jitter_sigma_ps=10.0),
        window_ps=1_000, reference_delay_ps=20.0, duration_s=0.2, dark_duration_s=1.0, dead_time_from_dark=True,
    ),
}

CORRECTED = "V corrected"
CONTROL = "V control"


def _replicate(regime: HomRegime, child: np.random.SeedSequence) -> tuple[bool, bool]:
    """(corrected interval covers the truth, control interval covers it) for one replicate."""
    cal_seq, run_seq, fit_seq = child.spawn(3)
    cal_rng, run_rng = np.random.default_rng(cal_seq), np.random.default_rng(run_seq)
    fit_seed = int(fit_seq.generate_state(1)[0])
    created = dt.datetime(2026, 11, 1, tzinfo=UTC)
    with tempfile.TemporaryDirectory() as tmp:
        for ch in CHANNELS:
            cal = calibrate_dark(simulate_dark(regime.detector, regime.dark_duration_s, cal_rng), ch, backend="sim", seed=fit_seed, n_draws=400)
            write_artifact(cal, tmp, accept_misfit=True, created=created)
        cset = CalibrationSet.load(
            tmp,
            analysis="hom_visibility",
            policy=POLICY,
            dataset=DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC)),
            kinds=None if regime.dead_time_from_dark else (),
        )
    runs = {
        name: simulate_hom(regime.source, regime.detector, regime.detector, delay, regime.duration_s, run_rng)
        for name, delay in (("dip", 0.0), ("reference", regime.reference_delay_ps))
    }
    window = Window(width_ps=regime.window_ps, source="plan")
    result, _ = hom_from_tags(runs, window, policy=POLICY, calibrations=cset, channels=CHANNELS, seed=fit_seed)
    c, truth = result.corrected, regime.truth
    half = 1.959964 * c.stderr_statistical
    return c.low <= truth <= c.high, c.value - half <= truth <= c.value + half


def hom_propagation(
    replicates: int = 1000, seed: int = 20261118, regime: str = "hom-bright", part: range | None = None
) -> CoverageResult:
    reg = HOM_REGIMES[regime]
    covered = {CORRECTED: 0, CONTROL: 0}
    for child in _children(seed, replicates, part, regime):
        ok, ok_control = _replicate(reg, child)
        covered[CORRECTED] += ok
        covered[CONTROL] += ok_control
    return CoverageResult(f"hom_visibility [{regime}]", replicates, covered, gate(replicates, 1))


def passed(result: CoverageResult) -> bool:
    """The corrected value is gated; the control is reported (the gate needs it low somewhere)."""
    lo, hi = result.gate
    return lo <= result.fraction(CORRECTED) <= hi


def run(replicates: int = 1000, seed: int = 20261118, regimes: list[str] | None = None) -> list[CoverageResult]:
    names = list(HOM_REGIMES) if regimes is None else regimes
    unknown = sorted(set(names) - set(HOM_REGIMES))
    if unknown:
        raise ValueError(f"unknown regime(s) {unknown}")
    return [hom_propagation(replicates, seed, r) for r in names]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m dadl.analysis.propagation")
    parser.add_argument("analysis", choices=["hom"])
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--regimes", nargs="+", help="run only these regimes (default: all)")
    args = parser.parse_args(argv)
    results = run(args.replicates, args.seed, args.regimes)
    for r in results:
        print(f"{r}\n  ({CONTROL}: calibration uncertainty off; reported, not gated per regime)")
    return 0 if all(passed(r) for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
