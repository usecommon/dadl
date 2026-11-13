"""Interval coverage of the calibrations on simulated data (the phase 2 gate).

Each replicate simulates one acquisition from known parameters with an
independent random stream, runs the calibration, and records whether each 95%
interval contains the truth. Run from the command line:

    python -m dadl.calibration.coverage {dark,pairs,dcr,all} --replicates 1000 [--jobs 4]

Each calibration is gated at more than one point of parameter space (a
"regime"), since a method can cover well at one point and fail at another:
the dark calibration's dead-time interval covered 96% at the original gate
point but 89% with fast afterpulsing (fixed in method v2).

`pairs` covers three calibrations from one simulated pair acquisition per
replicate (delay and jitter, accidentals, heralding efficiency), each with
its own gate. Truths for accidentals and heralding come from the closed forms
in `dadl.closed_form`.
"""

from __future__ import annotations

import argparse
import math
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

import numpy as np

from dadl.calibration.accidentals import calibrate_accidentals
from dadl.calibration.dark import calibrate_dark
from dadl.calibration.dcr import calibrate_dcr_temperature
from dadl.calibration.delay import calibrate_delay
from dadl.calibration.heralding import calibrate_heralding
from dadl.closed_form import accidentals_expected, arrhenius_dark_rate, heralding_ratio_expected
from dadl.data import ChannelId
from dadl.sim import DetectorParams, simulate_dark, simulate_pairs
from dadl.stats import coverage_band

#: Simulated detector for the dark calibration gate: a Si SPAD module with
#: visible afterpulsing, 20 s per replicate (about 20,000 dark counts).
DARK_TRUTH = DetectorParams(
    efficiency=0.5, dark_rate_cps=1000.0, dead_time_ps=22_000, afterpulse_prob=0.03, afterpulse_tau_ps=50_000.0
)
DARK_DURATION_S = 20.0

#: Dark-calibration regimes: (detector, seconds per replicate).
DARK_REGIMES: dict[str, tuple[DetectorParams, float]] = {
    "si-spad": (DARK_TRUTH, DARK_DURATION_S),
    # CMOS SPAD: strong, fast afterpulsing and a short dead time.
    "cmos-spad": (
        DetectorParams(dark_rate_cps=1000.0, dead_time_ps=10_000, afterpulse_prob=0.2, afterpulse_tau_ps=20_000.0),
        20.0,
    ),
}

#: Pair source for the delay, accidentals and heralding gates: a CW SPDC
#: source at 20,000 pairs/s, 2 s per replicate, a 2 ns fixed-grid window
#: with exact declared delays.
PAIR_RATE = 2e4
PAIR_A = DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, jitter_sigma_ps=350.0, delay_ps=0)
PAIR_B = DetectorParams(efficiency=0.4, dark_rate_cps=2000.0, jitter_sigma_ps=200.0, delay_ps=3_217)
PAIR_WIDTH_PS = 2_000
PAIR_DURATION_S = 2.0


@dataclass(frozen=True)
class PairRegime:
    rate: float
    a: DetectorParams
    b: DetectorParams
    width_ps: int
    duration_s: float
    #: Gates that only require coverage of at least 95% (discrete counts too low for an exact 95%).
    at_least: frozenset[str] = field(default_factory=frozenset)
    #: Calibrations not gated in this regime, with the reason recorded in the regime's comment.
    skip: frozenset[str] = field(default_factory=frozenset)


PAIR_REGIMES: dict[str, PairRegime] = {
    "cw-spad": PairRegime(PAIR_RATE, PAIR_A, PAIR_B, PAIR_WIDTH_PS, PAIR_DURATION_S),
    # Low-rate source on SNSPDs: few counts, a sharp peak, about 0.2 off-peak counts in all.
    # Declared before the first run: with so few counts any valid interval over-covers, so the
    # accidentals gate requires coverage of at least 95% only.
    "sparse-snspd": PairRegime(
        2e3,
        DetectorParams(efficiency=0.8, dark_rate_cps=50.0, dead_time_ps=30_000, jitter_sigma_ps=20.0, delay_ps=0),
        DetectorParams(efficiency=0.7, dark_rate_cps=50.0, dead_time_ps=30_000, jitter_sigma_ps=15.0, delay_ps=1_234),
        200,
        2.0,
        frozenset({"accidentals"}),
    ),
    # Si SPADs with afterpulsing (3% at 50 ns): afterpulses after true pairs must not reach the
    # off-peak windows (with the old 200 ns guard they biased the accidentals by +9%).
    "afterpulsing-spad": PairRegime(
        2e4,
        DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, jitter_sigma_ps=350.0, afterpulse_prob=0.03,
                       afterpulse_tau_ps=50_000.0, delay_ps=0),
        DetectorParams(efficiency=0.4, dark_rate_cps=2000.0, jitter_sigma_ps=200.0, afterpulse_prob=0.03,
                       afterpulse_tau_ps=50_000.0, delay_ps=3_217),
        2_000,
        2.0,
    ),
    # Bright background on the signal arm (8 Mcps recorded, 18% dead): dead time makes the
    # cross-correlation background step under the peak (a flat-background delay fit covered 52%)
    # and the accidentals are large. Heralding is not gated here: at this dead fraction misses
    # caused by dead time correlate with accidentals in the same slot, a limit the calibration
    # notes; its per-slot algebra is gated in the ideal regime below.
    "bright-background": PairRegime(
        2e5,
        DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, jitter_sigma_ps=350.0, delay_ps=0),
        DetectorParams(efficiency=0.4, dark_rate_cps=1e7, jitter_sigma_ps=200.0, delay_ps=3_217),
        2_000,
        0.5,
        skip=frozenset({"heralding-efficiency"}),
    ),
    # Broad timing response (multimode fibre, slow detectors): sigma about 1.8 ns.
    "broad-irf": PairRegime(
        2e4,
        DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, jitter_sigma_ps=1500.0, delay_ps=0),
        DetectorParams(efficiency=0.4, dark_rate_cps=2000.0, jitter_sigma_ps=1000.0, delay_ps=5_000),
        4_000,
        2.0,
    ),
}

#: DCR(T) gate: Arrhenius generation, 200 cps at 273 K, 0.5 eV, five
#: temperatures, a 10 s dark run at each.
DCR_REF_CPS = 200.0
DCR_EA_EV = 0.5
DCR_REF_K = 273.0
DCR_TEMPS_K = (253.0, 263.0, 273.0, 283.0, 293.0)
DCR_DURATION_S = 10.0

#: Chance that a correctly calibrated method fails the gate, per calibration.
#: Split across the calibration's parameters (Bonferroni), so a calibration
#: with four parameters does not fail four times as often as one with one.
#: 0.1%, not 1%: the gate runs on a fixed seed, so a chance failure would
#: block CI permanently, and with six calibrations at 1% each about 6% of
#: seeds fail somewhere. (Lowered from 1% after the delay/jitter gate read
#: 96.9% on seed 20261111 against a 96.8% limit, while three other seeds gave
#: 95.3%, 94.7% and 94.0%.)
FALSE_FAIL = 0.001


@dataclass(frozen=True)
class CoverageResult:
    calibration: str
    replicates: int
    covered: dict[str, int]
    gate: tuple[float, float]
    #: True when only the lower bound applies (coverage of at least 95%).
    at_least: bool = False

    def fraction(self, name: str) -> float:
        return self.covered[name] / self.replicates

    def passed(self) -> bool:
        lo, hi = self.gate
        return all(lo <= self.fraction(k) <= (1.0 if self.at_least else hi) for k in self.covered)

    def __str__(self) -> str:
        lo, hi = self.gate
        band = f"at least {lo:.1%}" if self.at_least else f"gate {lo:.1%}-{hi:.1%}"
        lines = [f"{self.calibration}: {self.replicates} replicates, {band}"]
        for k in self.covered:
            f = self.fraction(k)
            ok = lo <= f <= (1.0 if self.at_least else hi)
            lines.append(f"  {k:<20} {f:.1%}  {'ok' if ok else 'FAIL'}")
        return "\n".join(lines)


def _children(seed: int, replicates: int, part: range | None):
    """The SeedSequence of each replicate, child i of `seed`; `part` selects a subset of indices."""
    for i in range(replicates) if part is None else part:
        yield np.random.SeedSequence(seed, spawn_key=(i,))


def gate(replicates: int, n_params: int) -> tuple[float, float]:
    """Accepted coverage range: binomial band at confidence 1 - FALSE_FAIL / n_params."""
    return coverage_band(replicates, confidence=1 - FALSE_FAIL / n_params)


def dark_coverage(
    replicates: int = 1000, seed: int = 20261111, regime: str = "si-spad", part: range | None = None
) -> CoverageResult:
    det, duration_s = DARK_REGIMES[regime]
    truth = {
        "dark_rate_cps": det.dark_rate_cps,
        "afterpulse_prob": det.afterpulse_prob,
        "afterpulse_tau_ps": det.afterpulse_tau_ps,
        "dead_time_ps": float(det.dead_time_ps),
    }
    covered = dict.fromkeys(truth, 0)
    channel = ChannelId("sim", 0)
    for child in _children(seed, replicates, part):
        sim_rng, fit_seed = child.spawn(2)
        tags = simulate_dark(det, duration_s, np.random.default_rng(sim_rng))
        cal = calibrate_dark(tags, channel, backend="sim", seed=fit_seed.generate_state(1)[0], n_draws=10)
        for k, v in truth.items():
            iv = cal.estimates.get(k)
            if iv is not None and iv.contains(v):
                covered[k] += 1
    return CoverageResult(f"dark-afterpulse-deadtime [{regime}]", replicates, covered, gate(replicates, len(truth)))


def _covered(cal, truth: dict[str, float], covered: dict[str, int]) -> None:
    for k, v in truth.items():
        iv = cal.estimates.get(k)
        if iv is not None and iv.contains(v):
            covered[k] += 1


def pairs_truth(regime: str = "cw-spad") -> dict[str, dict[str, float]]:
    r = PAIR_REGIMES[regime]
    a, b = r.a, r.b
    ap_a = (a.afterpulse_prob, a.afterpulse_tau_ps * 1e-12)
    ap_b = (b.afterpulse_prob, b.afterpulse_tau_ps * 1e-12)
    rate_a = a.efficiency * r.rate + a.dark_rate_cps
    rate_b = b.efficiency * r.rate + b.dark_rate_cps
    jitter = math.hypot(a.jitter_sigma_ps, b.jitter_sigma_ps)
    return {
        "delay-jitter": {"delay_ps": float(b.delay_ps - a.delay_ps), "jitter_sigma_ps": jitter},
        "accidentals": {
            "accidental_coincidences": accidentals_expected(
                rate_a, rate_b, a.dead_time_ps * 1e-12, b.dead_time_ps * 1e-12, r.width_ps, r.duration_s, ap_a, ap_b
            )
        },
        "heralding-efficiency": {
            "heralding_efficiency": heralding_ratio_expected(
                r.rate, a.efficiency, b.efficiency, a.dark_rate_cps, b.dark_rate_cps,
                a.dead_time_ps * 1e-12, b.dead_time_ps * 1e-12, jitter, r.width_ps, ap_a, ap_b,
            )
        },
    }


def pairs_coverage(
    replicates: int = 1000, seed: int = 20261111, regime: str = "cw-spad", part: range | None = None
) -> list[CoverageResult]:
    """Delay/jitter, accidentals and heralding gates from shared simulated pair acquisitions."""
    r = PAIR_REGIMES[regime]
    truth = pairs_truth(regime)
    covered = {kind: dict.fromkeys(t, 0) for kind, t in truth.items() if kind not in r.skip}
    chs = (ChannelId("sim", 0), ChannelId("sim", 1))
    common = {"width_ps": r.width_ps, "backend": "sim", "n_draws": 10}
    for child in _children(seed, replicates, part):
        sim_rng, fit_seed = child.spawn(2)
        a, b = simulate_pairs(r.rate, r.a, r.b, r.duration_s, np.random.default_rng(sim_rng))
        s = int(fit_seed.generate_state(1)[0])
        delays = {"delay_a_ps": r.a.delay_ps, "delay_b_ps": r.b.delay_ps}
        if "delay-jitter" in covered:
            _covered(calibrate_delay(a, b, chs, backend="sim", seed=s, n_draws=10), truth["delay-jitter"], covered["delay-jitter"])
        # The accidentals measurand is this acquisition's: it follows the realised singles, so the
        # truth is n_a n_b W / T for the replicate (the long-run mean differs by the singles'
        # Poisson scatter, which outweighs the interval at high rates).
        acc_truth = {"accidental_coincidences": a.size * b.size * r.width_ps * 1e-12 / r.duration_s}
        if "accidentals" in covered:
            _covered(calibrate_accidentals(a, b, chs, seed=s, **delays, **common), acc_truth, covered["accidentals"])
        if "heralding-efficiency" in covered:
            her = calibrate_heralding(
                a, b, chs, delay_herald_ps=r.a.delay_ps, delay_signal_ps=r.b.delay_ps, seed=s, **common
            )
            _covered(her, truth["heralding-efficiency"], covered["heralding-efficiency"])
    return [
        CoverageResult(f"{kind} [{regime}]", replicates, c, gate(replicates, len(c)), kind in r.at_least)
        for kind, c in covered.items()
    ]


#: Ideal heralding regime: one herald per slot, its partner detected with probability ETA in the
#: same slot, and an accidental with probability P_ACC in every slot, independently; no dead
#: time. The exact truth is ETA, so this gates the per-slot algebra itself (method v2's
#: (C - lambda)/S reads ETA (1 - P_ACC), about three interval half-widths low here).
IDEAL_ETA = 0.5
IDEAL_P_ACC = 0.05
IDEAL_HERALDS = 20_000
IDEAL_SPACING = 2_000


def heralding_ideal_coverage(replicates: int = 1000, seed: int = 20261111, part: range | None = None) -> CoverageResult:
    width = 1_000
    covered = {"heralding_efficiency": 0}
    chs = (ChannelId("sim", 0), ChannelId("sim", 1))
    for child in _children(seed, replicates, part):
        sim_rng, fit_seed = child.spawn(2)
        rng = np.random.default_rng(sim_rng)
        # Heralds are sparse (one slot in IDEAL_SPACING) so that the shifted windows, which also
        # see other heralds' partners, read p within 1/(2 IDEAL_SPACING) of IDEAL_P_ACC; with one
        # herald in 20 slots they read 0.074 against 0.05, which the herald's own slot never sees.
        span = IDEAL_SPACING * IDEAL_HERALDS
        margin = 2_000
        slots = np.sort(rng.choice(span, IDEAL_HERALDS, replace=False)) + margin
        herald = slots * width + width // 2
        partner = herald[rng.random(slots.size) < IDEAL_ETA] + 10
        # Background: a Bernoulli process with probability IDEAL_P_ACC per slot (geometric gaps).
        n_slots = span + 2 * margin
        gaps = rng.geometric(IDEAL_P_ACC, int(n_slots * IDEAL_P_ACC * 1.2) + 1000)
        bg_slots = np.cumsum(gaps) - 1
        background = bg_slots[bg_slots < n_slots] * width + 300
        signal = np.sort(np.concatenate([partner, background]))
        cal = calibrate_heralding(herald, signal, chs, width_ps=width, guard_ps=200_000, backend="sim",
                                  seed=int(fit_seed.generate_state(1)[0]), n_draws=10)
        _covered(cal, {"heralding_efficiency": IDEAL_ETA}, covered)
    return CoverageResult("heralding-efficiency [ideal-accidentals]", replicates, covered, gate(replicates, 1))


def dcr_coverage(replicates: int = 1000, seed: int = 20261111, part: range | None = None) -> CoverageResult:
    """DCR(T) gate: dark calibrations at each temperature, then the Arrhenius fit."""
    truth = {"dcr_ref_cps": DCR_REF_CPS, "activation_energy_ev": DCR_EA_EV}
    covered = dict.fromkeys(truth, 0)
    channel = ChannelId("sim", 0)
    for child in _children(seed, replicates, part):
        streams = child.spawn(len(DCR_TEMPS_K) + 1)
        darks = []
        for t_k, stream in zip(DCR_TEMPS_K, streams):
            det = DetectorParams(dark_rate_cps=arrhenius_dark_rate(DCR_REF_CPS, DCR_EA_EV, t_k, DCR_REF_K))
            sim_rng, fit_seed = stream.spawn(2)
            tags = simulate_dark(det, DCR_DURATION_S, np.random.default_rng(sim_rng))
            darks.append(
                calibrate_dark(tags, channel, backend="sim", seed=int(fit_seed.generate_state(1)[0]), n_draws=10,
                               conditions={"temperature_K": t_k})
            )
        cal = calibrate_dcr_temperature(darks, reference_k=DCR_REF_K, seed=int(streams[-1].generate_state(1)[0]), n_draws=10)
        _covered(cal, truth, covered)
    return CoverageResult("dcr-temperature [five-temperatures]", replicates, covered, gate(replicates, len(truth)))


UNITS: dict[str, tuple[str, ...]] = {
    "dark": tuple(DARK_REGIMES),
    "pairs": (*PAIR_REGIMES, "ideal-accidentals"),
    "dcr": ("five-temperatures",),
}


#: Replicates per parallel task: small enough to keep every core busy to the end.
CHUNK = 125


def _unit(args: tuple[str, str, int, int, range | None]) -> list[CoverageResult]:
    name, regime, replicates, seed, part = args
    if name == "dark":
        return [dark_coverage(replicates, seed, regime, part)]
    if name == "pairs" and regime == "ideal-accidentals":
        return [heralding_ideal_coverage(replicates, seed, part)]
    if name == "pairs":
        return pairs_coverage(replicates, seed, regime, part)
    return [dcr_coverage(replicates, seed, part)]


def _merge(parts: list[CoverageResult], replicates: int) -> CoverageResult:
    first = parts[0]
    covered = {k: sum(p.covered[k] for p in parts) for k in first.covered}
    return CoverageResult(first.calibration, replicates, covered, first.gate, first.at_least)


def run(name: str, replicates: int = 1000, seed: int = 20261111, jobs: int = 1) -> list[CoverageResult]:
    """Every regime of one calibration family ("dark", "pairs", "dcr") or of all ("all").

    With jobs > 1 each regime's replicates are split into chunks run in parallel. Replicate i
    always uses child i of `seed`, so the result does not depend on `jobs`.
    """
    names = list(UNITS) if name == "all" else [name]
    if any(n not in UNITS for n in names):
        raise ValueError(f"unknown calibration {name!r}")
    units = [(n, r) for n in names for r in UNITS[n]]
    if jobs <= 1:
        return [res for n, r in units for res in _unit((n, r, replicates, seed, None))]
    tasks = [
        (n, r, replicates, seed, range(lo, min(lo + CHUNK, replicates)))
        for n, r in units
        for lo in range(0, replicates, CHUNK)
    ]
    # Fresh ("spawn") workers with one BLAS thread each: forked workers inherit the parent's
    # OpenBLAS thread pool on Linux, and several pools spinning on few cores slowed the CI gate
    # from minutes to over an hour.
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
    out: list[CoverageResult] = []
    for n, r in units:
        mine = [res for t, part in zip(tasks, parts) if (t[0], t[1]) == (n, r) for res in part]
        for cal in dict.fromkeys(res.calibration for res in mine):
            out.append(_merge([res for res in mine if res.calibration == cal], replicates))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m dadl.calibration.coverage")
    parser.add_argument("calibration", choices=["dark", "pairs", "dcr", "all"])
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument("--jobs", type=int, default=1, help="regimes run in parallel processes")
    args = parser.parse_args(argv)
    results = run(args.calibration, args.replicates, args.seed, args.jobs)
    for r in results:
        print(r)
    return 0 if all(r.passed() for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
