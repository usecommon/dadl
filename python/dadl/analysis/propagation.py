"""The `propagation` gate (phase 3): coverage of corrected values with calibration uncertainty.

    python -m dadl.analysis.propagation {hom,chsh,all} --replicates 1000 --jobs 4
    python -m dadl.analysis.propagation chsh --regimes chsh-dead-time --shard 1/4 --out part1.json
    python -m dadl.analysis.propagation merge part*.json

A shard runs one contiguous slice of the replicates and writes its counts without gating; `merge`
sums the shards and applies the gate. Replicate i always uses child i of the seed, so shards and
jobs do not change the result.

Each replicate simulates its own calibration runs (dark acquisitions per channel and a pair run
per channel pair for the delay and jitter), calibrates them, writes and loads the artifacts as an
analysis would, simulates the analysed acquisitions, and checks whether the corrected value's 95%
interval (statistical plus calibration, by joint Monte Carlo) covers the injected truth. The
accepted range is the binomial band of the phase 2 coverage gate
(`dadl.calibration.coverage.gate`).

A control runs alongside: the same corrected value with the calibration uncertainty switched off
(point +- 1.96 statistical standard deviations). It must fall below the band in at least one
declared regime (`CONTROL_REGIMES`), or the regimes do not test propagation.

The forward model's predictive check (`dadl.forward_check`) runs in every replicate. In the clean
regimes it may flag at most `MAX_FLAG_SHARE` of replicates. The misspecified regime (`MISSPECIFIED`:
crosstalk between the two HOM detectors, which the model lacks, in the window) passes only if the
check flags at least `MIN_FLAG_SHARE` of its replicates; its coverage is reported, not gated
(it is expected to fail: that is what the check is for).

Regimes are declared before the first run (phase 3 plan, M3).
"""

from __future__ import annotations

import argparse
import json
import datetime as dt
import multiprocessing
import os
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
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
from dadl.forward_check import FORWARD_MISFIT_ACCEPTED
from dadl.guardrails import Policy
from dadl.propagate import MISFIT_ACCEPTED, CalibrationSet, DatasetContext
from dadl.sim import Crosstalk, DetectorParams, PairSource, simulate_chsh, simulate_dark, simulate_hom, simulate_pairs

UTC = dt.timezone.utc
CREATED = dt.datetime(2026, 11, 1, tzinfo=UTC)
ACQUIRED = dt.datetime(2026, 11, 2, tzinfo=UTC)
HOM_CHANNELS = (ChannelId("sim", 0), ChannelId("sim", 1))
CHSH_CHANNELS = {name: ChannelId("sim", i) for i, name in enumerate(("A+", "A-", "B+", "B-"))}
CHSH_SETTINGS = {"a": 0.0, "a'": 45.0, "b": 22.5, "b'": 67.5}
#: Corrections the gate allows; misfit artifacts are accepted (labelled) and counted, since the
#: artifacts' predictive checks flag about 0.1-0.2% of clean calibrations each and a CHSH
#: replicate writes up to eight.
CORRECTIONS = ("dead_time", "accidental_subtraction", "window_capture", "misfit_calibration", "misfit_forward_model")


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
    #: Crosstalk between the two HOM detectors ("A", "B"), in every acquisition that uses them;
    #: the forward model lacks it, so a regime with crosstalk is misspecified.
    crosstalk: tuple[Crosstalk, ...] = ()

    @property
    def misspecified(self) -> bool:
        return bool(self.crosstalk)

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
    # Misspecified: hom-bright with 0.1% crosstalk each way at 400 ps (inside the 1 ns window),
    # about 400 extra coincidences on a dip of about 4,000 and V several standard errors low.
    # The predictive check must flag it (at 0.02% it flagged 73% of trial runs, at 0.05% all).
    "hom-crosstalk": Regime(
        "hom",
        PairSource(pair_rate_per_s=2.0e6, visibility=0.92, coherence_time_ps=2.0),
        DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=60_000, jitter_sigma_ps=10.0),
        window_ps=1_000, duration_s=0.2, dark_duration_s=1.0, pairs_duration_s=0.05,
        dead_time_from_dark=True, reference_delay_ps=20.0,
        crosstalk=(Crosstalk("A", "B", 0.001, 400.0, 50.0), Crosstalk("B", "A", 0.001, 400.0, 50.0)),
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
FLAGGED = "forward misfit"
KEYS = (CORRECTED, CONTROL, MISFITS, FLAGGED)
#: Largest share of clean replicates the predictive check may flag (phase 3 plan, M3 acceptance).
MAX_FLAG_SHARE = 0.01
#: Smallest share of misspecified replicates it must flag.
MIN_FLAG_SHARE = 0.90
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
        a, b = simulate_pairs(regime.source.pair_rate_per_s, det, det, regime.pairs_duration_s, rng, _lower(regime.crosstalk))
        delay = calibrate_delay(a, b, (a_id, b_id), backend="sim", seed=fit_seed, n_draws=400)
        write_artifact(delay, tmp, accept_misfit=True, created=CREATED)
    analysis = "hom_visibility" if regime.analysis == "hom" else "chsh"
    return CalibrationSet.load(tmp, analysis=analysis, policy=POLICY, dataset=DatasetContext(acquired=ACQUIRED))


def _lower(links: tuple[Crosstalk, ...]) -> tuple[Crosstalk, ...]:
    """The HOM links with channels named as `simulate_pairs` names them ("a", "b")."""
    return tuple(Crosstalk(x.source.lower(), x.target.lower(), x.probability, x.delay_ps, x.jitter_sigma_ps) for x in links)


def _replicate(regime: Regime, child: np.random.SeedSequence) -> tuple[bool, bool, bool, bool]:
    """(corrected interval covers the truth, control interval covers it, a misfit artifact was used,
    the forward model's predictive check flagged the analysis)."""
    cal_seq, run_seq, fit_seq = child.spawn(3)
    cal_rng, run_rng = np.random.default_rng(cal_seq), np.random.default_rng(run_seq)
    fit_seed = int(fit_seq.generate_state(1)[0])
    with tempfile.TemporaryDirectory() as tmp:
        cset = _calibrate(regime, cal_rng, fit_seed, tmp)
    window = Window(width_ps=regime.window_ps, source="plan")
    det = regime.detector
    if regime.analysis == "hom":
        runs = {
            name: simulate_hom(regime.source, det, det, delay, regime.duration_s, run_rng, regime.crosstalk)
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
    covered = c.low <= truth <= c.high
    return covered, c.value - half <= truth <= c.value + half, MISFIT_ACCEPTED in cset.labels, FORWARD_MISFIT_ACCEPTED in result.labels


def propagation(
    replicates: int = 1000, seed: int = 20261118, regime: str = "hom-bright", part: range | None = None
) -> CoverageResult:
    reg = REGIMES[regime]
    covered = dict.fromkeys(KEYS, 0)
    for child in _children(seed, replicates, part, regime):
        for key, hit in zip(KEYS, _replicate(reg, child)):
            covered[key] += hit
    name = "hom_visibility" if reg.analysis == "hom" else "chsh"
    return CoverageResult(f"{name} [{regime}]", replicates, covered, gate(replicates, 1))


def regime_of(result: CoverageResult) -> str:
    return result.calibration.split("[", 1)[1].rstrip("]")


def passed(result: CoverageResult) -> bool:
    """Clean regimes: the corrected value within the band, few misfit artifacts and few flags; in
    a control regime, the control below the band. The misspecified regime: flagged enough."""
    if REGIMES[regime_of(result)].misspecified:
        return result.fraction(FLAGGED) >= MIN_FLAG_SHARE
    lo, hi = result.gate
    ok = lo <= result.fraction(CORRECTED) <= hi and result.fraction(MISFITS) <= MAX_MISFIT_SHARE
    ok = ok and result.fraction(FLAGGED) <= MAX_FLAG_SHARE
    if regime_of(result) in CONTROL_REGIMES:
        ok = ok and result.fraction(CONTROL) < lo
    return ok


def describe(result: CoverageResult) -> str:
    lo, hi = result.gate
    verdict = "PASS" if passed(result) else "FAIL"
    if REGIMES[regime_of(result)].misspecified:
        return "\n".join([
            f"{result.calibration}: {result.replicates} replicates, misspecified  {verdict}",
            f"  {FLAGGED:<17} {result.fraction(FLAGGED):.1%}  (predictive check; at least {MIN_FLAG_SHARE:.0%})",
            f"  {CORRECTED:<17} {result.fraction(CORRECTED):.1%}  (reported: the model lacks the crosstalk)",
            f"  {CONTROL:<17} {result.fraction(CONTROL):.1%}  (reported)",
            f"  {MISFITS:<17} {result.fraction(MISFITS):.1%}  (reported)",
        ])
    control = "must fall below the band" if regime_of(result) in CONTROL_REGIMES else "reported, not gated here"
    return "\n".join([
        f"{result.calibration}: {result.replicates} replicates, gate {lo:.1%}-{hi:.1%}  {verdict}",
        f"  {CORRECTED:<17} {result.fraction(CORRECTED):.1%}",
        f"  {CONTROL:<17} {result.fraction(CONTROL):.1%}  (calibration uncertainty off; {control})",
        f"  {MISFITS:<17} {result.fraction(MISFITS):.1%}  (at most {MAX_MISFIT_SHARE:.0%})",
        f"  {FLAGGED:<17} {result.fraction(FLAGGED):.1%}  (predictive check; at most {MAX_FLAG_SHARE:.0%})",
    ])


#: Replicates per parallel task.
CHUNK = 50


def _unit(args: tuple[str, int, int, range | None]) -> CoverageResult:
    regime, replicates, seed, part = args
    return propagation(replicates, seed, regime, part)


def shard_range(shard: str, replicates: int) -> range:
    """Shard "k/n" (1-based) as its contiguous slice of the replicates."""
    try:
        k, n = (int(x) for x in shard.split("/"))
    except ValueError:
        raise ValueError(f"shard must be k/n, got {shard!r}") from None
    if not 1 <= k <= n:
        raise ValueError(f"shard {shard!r}: need 1 <= k <= n")
    return range((k - 1) * replicates // n, k * replicates // n)


def run(
    analysis: str = "all",
    replicates: int = 1000,
    seed: int = 20261118,
    jobs: int = 1,
    regimes: list[str] | None = None,
    part: range | None = None,
    progress=None,
) -> list[CoverageResult]:
    """Every regime of `analysis`, or only `regimes`; with `part`, only those replicates (a shard:
    the counts are partial and the result is not gated). `progress` gets a line per finished chunk."""
    names = [r for r, reg in REGIMES.items() if analysis in ("all", reg.analysis)]
    if regimes is not None:
        unknown = sorted(set(regimes) - set(names))
        if unknown:
            raise ValueError(f"unknown regime(s) {unknown} for {analysis!r}")
        names = [r for r in names if r in regimes]
    span = part if part is not None else range(replicates)
    tasks = [(r, replicates, seed, range(lo, min(lo + CHUNK, span.stop))) for r in names for lo in range(span.start, span.stop, CHUNK)]
    start = time.monotonic()

    def report(done: int, task) -> None:
        if progress is not None:
            r, _, _, rng = task
            progress(f"{r}: replicates {rng.start}-{rng.stop - 1} done ({done}/{len(tasks)} chunks, {time.monotonic() - start:.0f} s)")

    if jobs <= 1:
        parts = []
        for i, t in enumerate(tasks):
            parts.append(_unit(t))
            report(i + 1, t)
    else:
        # One BLAS thread per fresh worker (see dadl.calibration.coverage.run).
        single = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")}
        saved = {k: os.environ.get(k) for k in single}
        os.environ.update(single)
        try:
            with ProcessPoolExecutor(jobs, mp_context=multiprocessing.get_context("spawn")) as ex:
                futures = {ex.submit(_unit, t): i for i, t in enumerate(tasks)}
                parts = [None] * len(tasks)
                for done, fut in enumerate(as_completed(futures), 1):
                    i = futures[fut]
                    parts[i] = fut.result()
                    report(done, tasks[i])
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
    out = []
    for r in names:
        mine = [p for t, p in zip(tasks, parts) if t[0] == r]
        covered = {k: sum(p.covered[k] for p in mine) for k in KEYS}
        out.append(CoverageResult(mine[0].calibration, replicates, covered, mine[0].gate))
    return out


def shard_record(results: list[CoverageResult], seed: int, part: range) -> dict:
    return {
        "seed": seed,
        "replicates": results[0].replicates,
        "range": [part.start, part.stop],
        "results": [{"calibration": r.calibration, "covered": r.covered} for r in results],
    }


def merge(records: list[dict]) -> list[CoverageResult]:
    """Sum shard records; refuses mixed seeds or totals, and a slice missing or run twice."""
    if not records:
        raise ValueError("no shards to merge")
    seeds = {r["seed"] for r in records}
    totals = {r["replicates"] for r in records}
    if len(seeds) != 1 or len(totals) != 1:
        raise ValueError(f"shards disagree: seeds {sorted(seeds)}, replicates {sorted(totals)}")
    replicates = totals.pop()
    out = []
    names = dict.fromkeys(res["calibration"] for r in records for res in r["results"])
    for name in names:
        mine = [(r["range"], res["covered"]) for r in records for res in r["results"] if res["calibration"] == name]
        seen = sorted(tuple(rng) for rng, _ in mine)
        expected_next = 0
        for lo, hi in seen:
            if lo != expected_next:
                raise ValueError(f"{name}: shards do not tile 0-{replicates} (gap or overlap at {lo})")
            expected_next = hi
        if expected_next != replicates:
            raise ValueError(f"{name}: shards cover replicates 0-{expected_next}, not all {replicates}")
        covered = {k: sum(c[k] for _, c in mine) for k in KEYS}
        out.append(CoverageResult(name, replicates, covered, gate(replicates, 1)))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m dadl.analysis.propagation")
    parser.add_argument("analysis", choices=["hom", "chsh", "all", "merge"])
    parser.add_argument("files", nargs="*", help="merge: the shard files")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--regimes", nargs="+", help="run only these regimes (default: all of the analysis)")
    parser.add_argument("--shard", help="k/n: run only the k-th of n slices of the replicates; not gated")
    parser.add_argument("--out", help="with --shard: write the partial counts here, for merge")
    args = parser.parse_args(argv)
    if args.analysis == "merge":
        if not args.files:
            parser.error("merge needs the shard files")
        results = merge([json.loads(open(f, encoding="utf-8").read()) for f in args.files])
    else:
        if args.files:
            parser.error("files are only for merge")
        if (args.shard is None) != (args.out is None):
            parser.error("--shard and --out go together")
        part = shard_range(args.shard, args.replicates) if args.shard else None
        log = lambda line: print(line, file=sys.stderr, flush=True)  # noqa: E731
        results = run(args.analysis, args.replicates, args.seed, args.jobs, args.regimes, part, log)
        if part is not None:
            with open(args.out, "w", encoding="utf-8") as fh:
                json.dump(shard_record(results, args.seed, part), fh, indent=1)
            for r in results:
                print(f"{r.calibration}: shard {args.shard} (replicates {part.start}-{part.stop - 1}) {r.covered}")
            return 0
    for r in results:
        print(describe(r))
    return 0 if all(passed(r) for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
