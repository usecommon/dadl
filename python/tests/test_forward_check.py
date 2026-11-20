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


# Council round 7 (Q7): the check's statistics, below the propagation gate.


def test_grouped_pearson_degrees_of_freedom():
    from scipy import stats

    from dadl.forward_check import _grouped_pearson

    mu = np.full(8, 10.0)
    n = np.array([12, 8, 15, 5, 10, 11, 9, 10], dtype=float)
    assert _grouped_pearson(n, mu) == pytest.approx(stats.chi2.sf(np.sum((n - mu) ** 2 / mu), 7), rel=1e-12)


def test_bins_expecting_fewer_than_five_are_grouped():
    from scipy import stats

    from dadl.forward_check import MIN_EXPECTED, _grouped_pearson

    mu = np.full(6, 2.0)  # grouped in threes: 6 expected each
    n = np.array([0, 1, 0, 5, 4, 6], dtype=float)
    o, e = np.array([1.0, 15.0]), np.array([6.0, 6.0])
    assert MIN_EXPECTED == 5.0
    assert _grouped_pearson(n, mu) == pytest.approx(stats.chi2.sf(np.sum((o - e) ** 2 / e), 1), rel=1e-12)


def test_a_p_value_of_1e_4_is_flagged_and_bonferroni_applies():
    from dadl.forward_check import MISFIT_ALPHA, CheckResult, _grouped_pearson

    mu = np.full(10, 100.0)
    n = mu + np.array([25, -25, 25, -25, 20, -20, 0, 0, 0, 0])  # p = 1.3e-4
    p = _grouped_pearson(n, mu)
    assert 1e-6 < p < 1e-3
    assert CheckResult({"dip": p}, MISFIT_ALPHA).flagged
    assert CheckResult({"a": 0.01, "b": 0.5}, MISFIT_ALPHA).pvalue == pytest.approx(0.02)
    assert not CheckResult({"a": 0.0006, "b": 0.5}, MISFIT_ALPHA).flagged
    # Council round 8 (Fix 4): the check's share of the 1e-3 per-analysis budget is 4e-4.
    assert MISFIT_ALPHA == 4e-4
    assert CheckResult({"dip": 3.9e-4}, MISFIT_ALPHA).flagged and not CheckResult({"dip": 5e-4}, MISFIT_ALPHA).flagged


def _t3(t, period, res):
    """Times as a T3 file records them: sync count x period + floor(phase / res) x res, rounded to ps."""
    from dadl.adapters.ptu import _t3_to_ps

    nsync = np.floor(t / period).astype(np.int64)
    dtime = np.floor((t - nsync * period) / res).astype(np.int64)
    return _t3_to_ps(nsync, dtime, period * 1e-12, res * 1e-12)[0]


def _t3_pvalue(seed, res, period, cable, lattice):
    from dadl.forward_check import pair_pvalue
    from dadl.sim import simulate_pairs

    det = DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, dead_time_ps=22_000, jitter_sigma_ps=100.0)
    a, b = simulate_pairs(1e5, det, det, 0.5, np.random.default_rng(7000 + seed))
    b = b + int(cable)
    a, b = a[a >= 0], b[b >= 0]
    pair = PairData("A/B", _t3(a, period, res), _t3(b, period, res), float(int(cable)), 1_000, 0.0, 140.0,
                    lattice_ps=lattice)
    return pair_pvalue(pair)


def test_t3_tags_with_a_measured_sync_period_are_not_flagged():
    # Council round 7 (G1) left this an xfail: on T3 tags at 25 ps dtime and a 12500.3 ps sync the
    # check read p ~ 1e-100, binning 46 ps bins on the 1 ps step the tags show. Council round 8
    # (Fix 1): it bins on the dtime lattice the file's metadata gives.
    from dadl.forward_check import pair_pvalue
    from dadl.sim import simulate_pairs

    det = DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, dead_time_ps=22_000, jitter_sigma_ps=100.0)
    a, b = simulate_pairs(1e5, det, det, 0.5, np.random.default_rng(0))
    a, b = a[a >= 0], b[b >= 0]
    a, b = _t3(a, 12_500.3, 25), _t3(b, 12_500.3, 25)
    assert pair_pvalue(PairData("A/B", a, b, 0.0, 1_000, 0.0, 140.0, lattice_ps=25)) > 1e-3
    # Without the lattice (the tags show a 1 ps step) it misfires, as before.
    assert pair_pvalue(PairData("A/B", a, b, 0.0, 1_000, 0.0, 140.0)) < 1e-10


T3_CASES = [(25, 12_500.3, 30_000.0), (16, 12_345.6, 12_345.6), (16, 12_345.6, 0.0)]


@pytest.mark.parametrize("res, period, cable", T3_CASES)
def test_t3_lattice_check_one_seed(res, period, cable):
    # Council round 8 (Fix 1): 16 ps at 12345.6 ps, and cable delays that span one or two syncs.
    assert _t3_pvalue(0, res, period, cable, res) > 1e-3


@pytest.mark.slow
@pytest.mark.parametrize("res, period, cable", T3_CASES)
def test_t3_lattice_check_is_calibrated_over_seeds(res, period, cable):
    # r8-stat t3check2.py: p < 0.05 in 3% of 40 seeds per case; at most 4 of 40 allowed.
    ps = np.array([_t3_pvalue(seed, res, period, cable, res) for seed in range(40)])
    assert np.sum(ps < 0.05) <= 4


def test_t3_input_without_an_integer_lattice_skips_the_check():
    # Council round 8 (Fix 1 d): a dtime resolution that is not whole ps has no lattice to bin on.
    from dadl.analysis._twopoint import forward_check
    from dadl.data import ChannelId, RunTags, runs_lattice

    t3 = {"mode": "T3", "sync_period_ps": 12_500.3, "tick_ps": 26, "lattice_ps": None}
    runs = {"dip": RunTags({"A": np.arange(5), "B": np.arange(5)}, metadata=t3)}
    assert runs_lattice(runs) == 0
    assert runs_lattice({"dip": RunTags({}, metadata={**t3, "lattice_ps": 25})}) == 25
    assert runs_lattice({"dip": {"A": np.arange(3)}}) is None

    class _Set:
        members = ()

    p, labels = forward_check("hom_visibility", {"dip": (runs["dip"], (0.0,))}, ("A", "B"),
                              (ChannelId("sim", 0), ChannelId("sim", 1)), window_ps=1000, delays_ps={},
                              calibrations=_Set(), policy=None, lattice_ps=0)
    assert p is None and labels == ()


def test_the_alignment_puts_the_lattice_mid_bin_and_keeps_the_delay():
    # Council round 8 (Fix 1 c): the bin edges sit half a lattice step off the modal residue, and
    # the fitted offset is read back at the right place (a 40 ps T2-like lattice, delay 400 ps).
    from dadl.forward_check import _align, pooled_peak
    from dadl.sim import simulate_pairs

    det = DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, dead_time_ps=22_000, jitter_sigma_ps=60.0)
    a, b = simulate_pairs(2e5, det, det, 0.2, np.random.default_rng(5))
    a, b = np.sort(a[a >= 0] // 40 * 40), np.sort((b[b >= 0] + 400) // 40 * 40)
    pair = PairData("A/B", a, b, 400.0, 2_000, lattice_ps=40)
    shift, tooth = _align([pair], 400, 2_000, 40)
    assert tooth == 20 and (shift - 400) % 40 == 20
    fit = pooled_peak([pair])
    assert abs(fit.offset_ps) < 3 * fit.offset_se_ps + 1.0


def test_the_check_starts_from_the_artifact_in_either_order(tmp_path):
    # Council round 7 (Q2): a (B, A) artifact gives the delay of B after A with its sign flipped.
    from test_forward import CA, CB, SIGN_DET, SIGN_LATE, sign_sets

    from dadl.forward_check import delay_start
    from dadl.sim import simulate_pairs

    pair = simulate_pairs(2e5, SIGN_DET, SIGN_LATE, 0.05, np.random.default_rng(34))
    ab, ba = sign_sets(tmp_path, "hom_visibility", Policy.default(), {(CA, CB): pair})
    d_ab, j_ab = delay_start(ab, CA, CB)
    d_ba, j_ba = delay_start(ba, CA, CB)
    # The two fits of one run, mirrored, agree to well under a ps.
    assert d_ab == pytest.approx(400, abs=10) and d_ba == pytest.approx(d_ab, abs=1.0)
    assert j_ba == pytest.approx(j_ab, rel=0.01)
    assert delay_start(ba, CB, CA)[0] == pytest.approx(-d_ba, abs=1e-9)
