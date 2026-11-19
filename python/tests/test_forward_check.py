"""The forward model's predictive check (phase 3, M3). Its rates over 1,000 replicates are in the
`propagation` gate; these are the mechanics."""

from __future__ import annotations

import numpy as np
import pytest

from dadl.errors import Misspecified
from dadl.forward_check import FORWARD_MISFIT_ACCEPTED, OVERRIDE, PairData, enforce, run_check
from dadl.guardrails import Policy
from dadl.sim import Crosstalk, DetectorParams, PairSource, simulate_hom

SRC = PairSource(pair_rate_per_s=2.0e6, visibility=0.92, coherence_time_ps=2.0)
DET = DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=60_000, jitter_sigma_ps=10.0)
JITTER = 10.0 * 2**0.5
XT = (Crosstalk("A", "B", 0.001, 400.0, 50.0), Crosstalk("B", "A", 0.001, 400.0, 50.0))


def _pairs(seed: int, crosstalk=(), declared_offsets=None):
    rng = np.random.default_rng(seed)
    out = []
    for name, d in (("dip", 0.0), ("reference", 20.0)):
        run = simulate_hom(SRC, DET, DET, d, 0.2, rng, crosstalk)
        offsets = declared_offsets if declared_offsets is not None and d else ((d, -d) if d else (0.0,))
        out.append(PairData(name, run["A"], run["B"], 0.0, 1_000, 0.0, JITTER, offsets))
    return out


def test_clean_runs_pass_and_crosstalk_is_flagged():
    clean = run_check(_pairs(1))
    assert set(clean.pvalues) == {"dip", "reference"} and not clean.flagged
    bad = run_check(_pairs(1, XT))
    assert bad.flagged and bad.pvalue < 1e-6


def test_a_misdeclared_reference_delay_is_flagged():
    r = run_check(_pairs(2, declared_offsets=(0.0,)))
    assert r.flagged and r.worst() == "reference"


def test_enforce_refuses_unless_allowed():
    bad = run_check(_pairs(3, XT))
    with pytest.raises(Misspecified, match=OVERRIDE):
        enforce(bad, "hom_visibility", Policy.default())
    loose = Policy.default().allow("hom_visibility", OVERRIDE, "test")
    assert enforce(bad, "hom_visibility", loose) == (FORWARD_MISFIT_ACCEPTED,)
    assert enforce(run_check(_pairs(3)), "hom_visibility", Policy.default()) == ()


def test_sparse_histograms_are_skipped():
    empty = np.array([], dtype=np.int64)
    r = run_check([PairData("none", empty, np.arange(10, dtype=np.int64), 0.0, 1_000)])
    assert r.pvalues == {} and r.pvalue == 1.0 and not r.flagged


def test_a_peak_that_fills_the_window_is_fitted_from_the_outer_bins(tmp_path):
    # Council round 2: a 3% timing tail, 300 ps window, 80 ps reference delay. The background start
    # from the floor stalled (p about 1e-206 on the reference histogram at this seed); the start
    # from the outer bins reaches the fit (p about 0.4).
    from dadl.analysis.chsh import Window
    from dadl.analysis.hom import hom_from_tags
    from dadl.analysis.propagation import ACQUIRED, CREATED, HOM_CHANNELS, POLICY
    from dadl.artifact import write_artifact
    from dadl.calibration.delay import calibrate_delay
    from dadl.propagate import CalibrationSet, DatasetContext
    from dadl.sim import simulate_pairs

    def tail(tags, rng):
        out = {}
        for k, t in tags.items():
            t = np.asarray(t).copy()
            m = rng.random(t.size) < 0.03
            t[m] += np.round(rng.exponential(100.0, m.sum())).astype(np.int64)
            out[k] = np.unique(t)
        return out

    src = PairSource(pair_rate_per_s=2.0e5, visibility=0.92, coherence_time_ps=10.0)
    det = DetectorParams(efficiency=0.85, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=30.0)
    rng = np.random.default_rng(2007)
    a, b = simulate_pairs(src.pair_rate_per_s, det, det, 0.05, rng, ())
    tb = tail({"a": a, "b": b}, rng)
    cal = calibrate_delay(tb["a"], tb["b"], HOM_CHANNELS, backend="sim", seed=2007, n_draws=400)
    write_artifact(cal, tmp_path, accept_misfit=True, created=CREATED)
    cset = CalibrationSet.load(tmp_path, analysis="hom_visibility", policy=POLICY, dataset=DatasetContext(acquired=ACQUIRED))
    runs = {n: tail(simulate_hom(src, det, det, d, 0.5, rng, ()), rng) for n, d in (("dip", 0.0), ("reference", 80.0))}
    _, counts = hom_from_tags(
        runs, Window(width_ps=300, source="plan"), policy=POLICY, calibrations=cset, channels=HOM_CHANNELS,
        seed=2007, reference_delay_ps=80.0,
    )
    assert counts["forward_check_p"] > 0.1


@pytest.mark.slow
@pytest.mark.parametrize("seed", [36, 43])
def test_a_start_delay_tens_of_ps_off_does_not_stall_the_fit(tmp_path, seed):
    # Council round 3 (G3): with the start 14-27 ps off the peak, L-BFGS-B stopped on the
    # background step and a clean histogram read p 2e-6 (seed 36) and 7e-7 (seed 43); polished,
    # 0.69 and 0.82 here.
    import dataclasses
    import datetime as dt

    from dadl.analysis.chsh import SIGNS
    from dadl.artifact import write_artifact
    from dadl.calibration.delay import calibrate_delay
    from dadl.data import ChannelId
    from dadl.forward_check import delay_start, pair_pvalue
    from dadl.propagate import CalibrationSet
    from dadl.sim import simulate_chsh, simulate_pairs

    ids = {n: ChannelId("sim", i) for i, n in enumerate(("A+", "A-", "B+", "B-"))}
    det = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=350.0)
    dets = {n: dataclasses.replace(det, delay_ps=400 if n in ("A+", "B+") else 0) for n in ids}
    rng = np.random.default_rng(5000 + seed)
    for x in "+-":
        for y in "+-":
            a, b = simulate_pairs(2e5, dets["A" + x], dets["B" + y], 0.05, rng)
            cal = calibrate_delay(a, b, (ids["A" + x], ids["B" + y]), backend="sim", seed=1, n_draws=400)
            write_artifact(cal, tmp_path, accept_misfit=True, created=dt.datetime(2026, 11, 1, tzinfo=dt.timezone.utc))
    cset = CalibrationSet.load(tmp_path, analysis="chsh")
    settings = {"a": 0.0, "a'": 45.0, "b": 22.5, "b'": 67.5}
    runs = simulate_chsh(PairSource(pair_rate_per_s=2e5, visibility=0.70), dets, 1.0, rng, settings)
    assert ("a", "b'") in SIGNS
    d_cal, jitter = delay_start(cset, ids["A+"], ids["B-"])
    run = runs[("a", "b'")]
    pair = PairData("a,b' A+/B-", np.sort(run["A+"]), np.sort(run["B-"]), 0.0, 2000, d_cal, jitter)
    # The stalled fit read about 1e-6; the polished one is far above it. The threshold leaves room
    # for platform numerics (0.096 on Windows, Python 3.12, seed 43).
    assert pair_pvalue(pair) > 0.01
