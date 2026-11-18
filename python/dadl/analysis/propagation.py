"""The `propagation` gate (phase 3): coverage of corrected values with calibration uncertainty.

    python -m dadl.analysis.propagation {hom,chsh,all} --replicates 1000 --jobs 4

Each replicate simulates its own calibration runs (dark acquisitions per channel and a pair run
per channel pair for the delay and jitter), calibrates them, writes and loads the artifacts as an
analysis would, simulates the analysed acquisitions, and checks whether the corrected value's 95%
interval (statistical plus calibration, by joint Monte Carlo) covers the injected truth. The
accepted range is the binomial band of the phase 2 coverage gate
(`dadl.calibration.coverage.gate`).

A control runs alongside: the same corrected value with the calibration uncertainty switched off
(point +- 1.96 statistical standard deviations). It must fall below the band in at least one
declared regime (`CONTROL_REGIMES`), or the regimes do not test propagation.

Regimes are declared before the first run (phase 3 plan, M3). Still to come in M3: the
misspecified regime, which needs the predictive check.
"""

from __future__ import annotations

import argparse
import datetime as dt
import multiprocessing
import os
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from dadl.analysis.chsh import Window, chsh
from dadl.analysis.hom import hom_from_tags
from dadl.artifact import write_artifact
from dadl.calibration.coverage import CoverageResult, _children, gate
from dadl.calibration.dark import calibrate_dark
from dadl.calibration.delay import calibrate_delay
from dadl.closed_form import chsh_ideal, hom_visibility_expected
from dadl.data import ChannelId
from dadl.guardrails import Policy
from dadl.propagate import MISFIT_ACCEPTED, CalibrationSet, DatasetContext
from dadl.sim import DetectorParams, PairSource, simulate_chsh, simulate_dark, simulate_hom, simulate_pairs

UTC = dt.timezone.utc
CREATED = dt.datetime(2026, 11, 1, tzinfo=UTC)
ACQUIRED = dt.datetime(2026, 11, 2, tzinfo=UTC)
HOM_CHANNELS = (ChannelId("sim", 0), ChannelId("sim", 1))
CHSH_CHANNELS = {name: ChannelId("sim", i) for i, name in enumerate(("A+", "A-", "B+", "B-"))}
CHSH_SETTINGS = {"a": 0.0, "a'": 45.0, "b": 22.5, "b'": 67.5}
#: Corrections the gate allows; misfit artifacts are accepted (labelled) and counted, since the
#: artifacts' predictive checks flag about 0.1-0.2% of clean calibrations each and a CHSH
#: replicate writes up to eight.
CORRECTIONS = ("dead_time", "accidental_subtraction", "window_capture", "misfit_calibration")


def _policy() -> Policy:
    policy = Policy.default()
    for analysis in ("hom_visibility", "chsh"):
        for correction in CORRECTIONS:
            policy = policy.allow(analysis, correction, "propagation gate: calibrations from simulated runs")
    return policy


POLICY = _policy()


@dataclass(frozen=True)
class Regime:
    analysis: str  # "hom" or "chsh"
    source: PairSource
    detector: DetectorParams
    window_ps: int
    duration_s: float
    dark_duration_s: float
    #: The pair run (no beamsplitter or analysers) that calibrates each channel pair's delay and jitter.
    pairs_duration_s: float
    #: Whether the dark calibration can identify the dead time; if not (100 cps on an SNSPD, 500
    #: cps on the example SPADs), the dark artifacts are not loaded and the dead time is not corrected.
    dead_time_from_dark: bool
    reference_delay_ps: float = 0.0

    @property
    def truth(self) -> float:
        """What the forward model targets: the source's two-point visibility at the declared delays
        (HOM), or S of the source's Werner state at the declared angles (CHSH)."""
        if self.analysis == "hom":
            return hom_visibility_expected(self.source.visibility, self.reference_delay_ps, self.source.coherence_time_ps)
        return chsh_ideal(self.source.visibility, CHSH_SETTINGS)


REGIMES: dict[str, Regime] = {
    # examples/hom-snspd: too few dark counts to identify the dead time.
    "hom-snspd": Regime(
        "hom",
        PairSource(pair_rate_per_s=2.0e5, visibility=0.92, coherence_time_ps=1.5),
        DetectorParams(efficiency=0.85, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=25.0),
        window_ps=200, duration_s=0.5, dark_duration_s=2.0, pairs_duration_s=0.1,
        dead_time_from_dark=False, reference_delay_ps=15.0,
    ),
    # Bright background on SPADs: accidentals and dead time both matter.
    "hom-bright": Regime(
        "hom",
        PairSource(pair_rate_per_s=2.0e6, visibility=0.92, coherence_time_ps=2.0),
        DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=60_000, jitter_sigma_ps=10.0),
        window_ps=1_000, duration_s=0.2, dark_duration_s=1.0, pairs_duration_s=0.05,
        dead_time_from_dark=True, reference_delay_ps=20.0,
    ),
    # examples/chsh-spad: 25k cps per channel, dead time negligible and not identifiable from 500
    # cps of dark counts; afterpulsing (0.5%) left in, unmodelled (uncorrelated, so in the
    # off-peak accidentals).
    "chsh-spad": Regime(
        "chsh",
        PairSource(pair_rate_per_s=1.0e5, visibility=0.95),
        DetectorParams(
            efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=350.0,
            afterpulse_prob=0.005, afterpulse_tau_ps=50_000,
        ),
        window_ps=2_000, duration_s=1.0, dark_duration_s=2.0, pairs_duration_s=0.2,
        dead_time_from_dark=False,
    ),
    # Dead-time selection (raw S about 1.6 sigma high at a source S of 1.98), with the dead time
    # from short dark runs: its uncertainty is comparable to the statistical error, so this is
    # the regime where switching it off must undercover.
    "chsh-dead-time": Regime(
        "chsh",
        PairSource(pair_rate_per_s=3.0e6, visibility=0.70),
        DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=50_000, jitter_sigma_ps=50.0),
        window_ps=1_000, duration_s=0.4, dark_duration_s=0.1, pairs_duration_s=0.02,
        dead_time_from_dark=True,
    ),
}
#: Regimes where the uncertainty-off control must fall below the band.
CONTROL_REGIMES = ("chsh-dead-time",)

CORRECTED = "corrected"
CONTROL = "control"
MISFITS = "misfit artifacts"
#: Largest share of replicates that may use a misfit artifact (the checks' nominal rate times up
#: to eight artifacts is about 1.2%).
MAX_MISFIT_SHARE = 0.02


def _calibrate(regime: Regime, rng: np.random.Generator, fit_seed: int, tmp: str) -> CalibrationSet:
    channels = HOM_CHANNELS if regime.analysis == "hom" else tuple(CHSH_CHANNELS.values())
    det = regime.detector
    if regime.dead_time_from_dark:
        for ch in channels:
            cal = calibrate_dark(simulate_dark(det, regime.dark_duration_s, rng), ch, backend="sim", seed=fit_seed, n_draws=400)
            write_artifact(cal, tmp, accept_misfit=True, created=CREATED)
    if regime.analysis == "hom":
        pairs = [HOM_CHANNELS]
    else:
        pairs = [(CHSH_CHANNELS["A" + x], CHSH_CHANNELS["B" + y]) for x in "+-" for y in "+-"]
    for a_id, b_id in pairs:
        a, b = simulate_pairs(regime.source.pair_rate_per_s, det, det, regime.pairs_duration_s, rng)
        delay = calibrate_delay(a, b, (a_id, b_id), backend="sim", seed=fit_seed, n_draws=400)
        write_artifact(delay, tmp, accept_misfit=True, created=CREATED)
    analysis = "hom_visibility" if regime.analysis == "hom" else "chsh"
    return CalibrationSet.load(tmp, analysis=analysis, policy=POLICY, dataset=DatasetContext(acquired=ACQUIRED))


def _replicate(regime: Regime, child: np.random.SeedSequence) -> tuple[bool, bool, bool]:
    """(corrected interval covers the truth, control interval covers it, a misfit artifact was used)."""
    cal_seq, run_seq, fit_seq = child.spawn(3)
    cal_rng, run_rng = np.random.default_rng(cal_seq), np.random.default_rng(run_seq)
    fit_seed = int(fit_seq.generate_state(1)[0])
    with tempfile.TemporaryDirectory() as tmp:
        cset = _calibrate(regime, cal_rng, fit_seed, tmp)
    window = Window(width_ps=regime.window_ps, source="plan")
    det = regime.detector
    if regime.analysis == "hom":
        runs = {
            name: simulate_hom(regime.source, det, det, delay, regime.duration_s, run_rng)
            for name, delay in (("dip", 0.0), ("reference", regime.reference_delay_ps))
        }
        result, _ = hom_from_tags(
            runs, window, policy=POLICY, calibrations=cset, channels=HOM_CHANNELS, seed=fit_seed,
            reference_delay_ps=regime.reference_delay_ps,
        )
    else:
        dets = dict.fromkeys(CHSH_CHANNELS, det)
        runs = simulate_chsh(regime.source, dets, regime.duration_s, run_rng, CHSH_SETTINGS)
        result = chsh(runs, window, POLICY, settings_deg=CHSH_SETTINGS, calibrations=cset, channels=CHSH_CHANNELS, seed=fit_seed)
    c, truth = result.corrected, regime.truth
    half = 1.959964 * c.stderr_statistical
    return c.low <= truth <= c.high, c.value - half <= truth <= c.value + half, MISFIT_ACCEPTED in cset.labels


def propagation(
    replicates: int = 1000, seed: int = 20261118, regime: str = "hom-bright", part: range | None = None
) -> CoverageResult:
    reg = REGIMES[regime]
    covered = {CORRECTED: 0, CONTROL: 0, MISFITS: 0}
    for child in _children(seed, replicates, part, regime):
        ok, ok_control, misfit = _replicate(reg, child)
        covered[CORRECTED] += ok
        covered[CONTROL] += ok_control
        covered[MISFITS] += misfit
    name = "hom_visibility" if reg.analysis == "hom" else "chsh"
    return CoverageResult(f"{name} [{regime}]", replicates, covered, gate(replicates, 1))


def regime_of(result: CoverageResult) -> str:
    return result.calibration.split("[", 1)[1].rstrip("]")


def passed(result: CoverageResult) -> bool:
    """The corrected value within the band; in a control regime, the control below it."""
    lo, hi = result.gate
    ok = lo <= result.fraction(CORRECTED) <= hi and result.fraction(MISFITS) <= MAX_MISFIT_SHARE
    if regime_of(result) in CONTROL_REGIMES:
        ok = ok and result.fraction(CONTROL) < lo
    return ok


def describe(result: CoverageResult) -> str:
    lo, hi = result.gate
    control = "must fall below the band" if regime_of(result) in CONTROL_REGIMES else "reported, not gated here"
    return "\n".join([
        f"{result.calibration}: {result.replicates} replicates, gate {lo:.1%}-{hi:.1%}  {'PASS' if passed(result) else 'FAIL'}",
        f"  {CORRECTED:<17} {result.fraction(CORRECTED):.1%}",
        f"  {CONTROL:<17} {result.fraction(CONTROL):.1%}  (calibration uncertainty off; {control})",
        f"  {MISFITS:<17} {result.fraction(MISFITS):.1%}  (at most {MAX_MISFIT_SHARE:.0%})",
    ])


#: Replicates per parallel task.
CHUNK = 50


def _unit(args: tuple[str, int, int, range | None]) -> CoverageResult:
    regime, replicates, seed, part = args
    return propagation(replicates, seed, regime, part)


def run(
    analysis: str = "all", replicates: int = 1000, seed: int = 20261118, jobs: int = 1, regimes: list[str] | None = None
) -> list[CoverageResult]:
    names = [r for r, reg in REGIMES.items() if analysis in ("all", reg.analysis)]
    if regimes is not None:
        unknown = sorted(set(regimes) - set(names))
        if unknown:
            raise ValueError(f"unknown regime(s) {unknown} for {analysis!r}")
        names = [r for r in names if r in regimes]
    if jobs <= 1:
        return [propagation(replicates, seed, r) for r in names]
    tasks = [(r, replicates, seed, range(lo, min(lo + CHUNK, replicates))) for r in names for lo in range(0, replicates, CHUNK)]
    # One BLAS thread per fresh worker (see dadl.calibration.coverage.run).
    single = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
    saved = {k: os.environ.get(k) for k in single}
    os.environ.update(single)
    try:
        with ProcessPoolExecutor(jobs, mp_context=multiprocessing.get_context("spawn")) as ex:
            parts = list(ex.map(_unit, tasks))
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    out = []
    for r in names:
        mine = [p for t, p in zip(tasks, parts) if t[0] == r]
        covered = {k: sum(p.covered[k] for p in mine) for k in (CORRECTED, CONTROL, MISFITS)}
        out.append(CoverageResult(mine[0].calibration, replicates, covered, mine[0].gate))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m dadl.analysis.propagation")
    parser.add_argument("analysis", choices=["hom", "chsh", "all"])
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--regimes", nargs="+", help="run only these regimes (default: all of the analysis)")
    args = parser.parse_args(argv)
    results = run(args.analysis, args.replicates, args.seed, args.jobs, args.regimes)
    for r in results:
        print(describe(r))
    return 0 if all(passed(r) for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
