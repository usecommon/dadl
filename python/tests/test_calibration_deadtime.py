"""Series dead time (phase 3, M4): simulator, gap law, calibration, artifact and CLI.

Interval coverage over 1,000 replicates is the `coverage` gate (regimes series-spad and
passive-quench); these tests check recovery at a few points and the mechanics.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest
from scipy import stats

from dadl.artifact import load_calibration, write_artifact
from dadl.calibration.dark import calibrate_dark
from dadl.calibration.deadtime import KIND, calibrate_deadtime_series, gap_cdf, survival
from dadl.closed_form import series_measured_rate
from dadl.data import ChannelId
from dadl.deadtime import apply_nonparalyzable, apply_paralyzable
from dadl.errors import InsufficientCounts
from dadl.sim import Crosstalk, DetectorParams, detect, detect_channels, poisson_times, simulate_dark

CH = ChannelId("sim", 0)
UTC = dt.timezone.utc


def _det(tau_p: int, tau_np: int, **kw) -> DetectorParams:
    base = {"efficiency": 1.0, "dark_rate_cps": 1000.0, "jitter_sigma_ps": 0.0}
    return DetectorParams(dead_time_ps=tau_np, dead_time_p_ps=tau_p, **{**base, **kw})


def _runs(det: DetectorParams, rates, duration_s: float, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    return [detect(poisson_times(r, duration_s, rng), det, duration_s, rng) for r in rates]


def test_simulator_series_matches_the_filter_chain_and_closed_form():
    det = _det(20_000, 50_000, dark_rate_cps=0.0)
    rng = np.random.default_rng(1)
    t = poisson_times(1e7, 0.02, rng)
    got = detect(t, det, 0.02, np.random.default_rng(2))
    arrivals = np.sort(np.rint(t)).astype(np.int64)
    assert np.array_equal(got, apply_nonparalyzable(apply_paralyzable(arrivals, 20_000), 50_000))
    assert got.size / 0.02 == pytest.approx(series_measured_rate(1e7, 20e-9, 50e-9), rel=0.02)


def test_event_loops_apply_the_paralyzable_stage_too():
    # The afterpulse and crosstalk loops must give the filter chain's result when they add nothing.
    det = _det(20_000, 50_000, dark_rate_cps=0.0)
    t = poisson_times(1e7, 0.005, np.random.default_rng(3))
    chain = detect(t, det, 0.005, np.random.default_rng(4))
    looped = detect(t, DetectorParams(**{**det.__dict__, "afterpulse_prob": 1e-12}), 0.005, np.random.default_rng(4))
    assert np.array_equal(chain, looped)
    coupled = detect_channels({"a": t, "b": np.empty(0)}, {"a": det, "b": det}, 0.005, np.random.default_rng(4),
                              [Crosstalk("b", "a", 0.0)])
    assert np.array_equal(chain, coupled["a"])


def test_paralyzable_stage_must_not_exceed_the_dead_time():
    with pytest.raises(ValueError, match="dead_time_p_ps"):
        DetectorParams(dead_time_ps=10_000, dead_time_p_ps=20_000)


@pytest.mark.parametrize("rate,tau_p", [(5e6, 20_000.0), (2e7, 30_000.0), (5e6, 0.0)])
def test_gap_law_matches_simulated_gaps(rate, tau_p):
    tau_np = 50_000
    tags = _runs(_det(int(tau_p), tau_np, dark_rate_cps=0.0), [rate], 0.05, 5)[0]
    u = np.diff(tags) - tau_np
    n = rate / 1e12
    assert stats.kstest(u, lambda x: gap_cdf(x, n, tau_p)).pvalue > 0.001


def test_survival_is_continuous_where_the_exact_series_hands_over():
    n, tau = 2e-5, 30_000.0  # n tau = 0.6
    w = np.array([8 * tau - 1e-3, 8 * tau + 1e-3])
    s = survival(w, n, tau)
    assert s[0] == pytest.approx(s[1], rel=1e-6)
    # n tau > 1: the dominant mode comes from the Lambert W branch, still continuous.
    n2 = 5e-5
    assert survival(w, n2, tau)[0] == pytest.approx(survival(w, n2, tau)[1], rel=1e-5)


@pytest.mark.parametrize(
    "tau_p,tau_np,rates",
    [(20_000, 50_000, (2e6, 5e6, 1e7)), (40_000, 45_000, (1e6, 3e6, 6e6)), (5_000, 22_000, (5e6, 1.5e7, 3e7))],
)
def test_recovery_at_three_points(tau_p, tau_np, rates):
    cal = calibrate_deadtime_series(_runs(_det(tau_p, tau_np), rates, 0.01, tau_p), CH, backend="sim", seed=1, n_draws=200)
    assert cal.kind == KIND and cal.status == "ok"
    assert cal.estimates["dead_time_p_ps"].contains(tau_p)
    assert cal.estimates["dead_time_ps"].contains(tau_np)
    assert "dead_time_p_ps" not in cal.unidentified
    assert cal.conditions["lr_hard_edge"] > 3.84
    assert cal.samples.shape == (200, 2) and np.all(cal.samples[:, 0] <= cal.samples[:, 1])


def test_hard_edge_detector_gives_an_upper_limit():
    cal = calibrate_deadtime_series(_runs(_det(0, 50_000), (2e6, 5e6, 1e7), 0.01, 7), CH, backend="sim", seed=1, n_draws=200)
    iv = cal.estimates["dead_time_p_ps"]
    assert iv.low == 0.0 and iv.high > 0 and cal.unidentified == ("dead_time_p_ps",)
    assert cal.stand_in == ("dead_time_p_ps",)
    assert any("not rejected" in n for n in cal.notes)


def test_afterpulsing_is_held_at_the_dark_artifact():
    det = _det(30_000, 50_000, afterpulse_prob=0.03, afterpulse_tau_ps=30_000.0)
    rng = np.random.default_rng(9)
    dark = calibrate_dark(simulate_dark(_det(30_000, 50_000, dark_rate_cps=50_000.0, afterpulse_prob=0.03,
                                             afterpulse_tau_ps=30_000.0), 0.5, rng), CH, backend="sim", seed=1, n_draws=400)
    runs = _runs(det, (2e6, 5e6, 1e7), 0.01, 10)
    with_ap = calibrate_deadtime_series(runs, CH, afterpulse=dark, backend="sim", seed=1, n_draws=200)
    assert with_ap.estimates["dead_time_p_ps"].contains(30_000)
    assert with_ap.conditions["afterpulse_prob_held"] == pytest.approx(dark.estimates["afterpulse_prob"].value)
    assert any("held at the dark artifact" in n for n in with_ap.notes)


def test_the_dark_artifact_is_cited():
    # Council round 5 (V1): two dark artifacts from different dark runs gave identical inputs.
    det = _det(30_000, 50_000, dark_rate_cps=50_000.0, afterpulse_prob=0.03, afterpulse_tau_ps=30_000.0)
    runs = _runs(_det(30_000, 50_000, afterpulse_prob=0.03, afterpulse_tau_ps=30_000.0), (2e6, 5e6), 0.01, 10)
    cals = []
    for seed in (1, 2):
        dark = calibrate_dark(simulate_dark(det, 0.5, np.random.default_rng(seed)), CH, backend="sim", seed=1,
                              n_draws=400)
        cal = calibrate_deadtime_series(runs, CH, afterpulse=dark, backend="sim", seed=1, n_draws=10)
        assert set(dark.inputs) <= set(cal.inputs)
        assert cal.conditions["afterpulse_prob_sd_held"] > 0 and cal.conditions["afterpulse_tau_sd_held"] > 0
        cals.append(cal)
    assert cals[0].inputs != cals[1].inputs


def test_refusals():
    runs = _runs(_det(20_000, 50_000), (2e6, 5e6), 0.01, 11)
    with pytest.raises(ValueError, match="at least 2"):
        calibrate_deadtime_series(runs[:1], CH)
    with pytest.raises(InsufficientCounts, match="gaps"):
        calibrate_deadtime_series([runs[0][:100], runs[1]], CH)


def test_artifact_round_trip(tmp_path):
    cal = calibrate_deadtime_series(_runs(_det(20_000, 50_000), (2e6, 5e6, 1e7), 0.005, 12), CH, backend="sim", seed=1,
                                    n_draws=100)
    path = write_artifact(cal, tmp_path, created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    back = load_calibration(path)
    assert back.kind == "deadtime-series"
    assert back.estimates["dead_time_p_ps"].value == pytest.approx(cal.estimates["dead_time_p_ps"].value)
    again = write_artifact(cal, tmp_path / "again", created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    assert again.name == path.name


def test_cli_scan_and_calibrate(tmp_path, capsys):
    from dadl.cli import main

    scan, dark, out = tmp_path / "scan.npz", tmp_path / "dark.npz", tmp_path / "arts"
    assert main(["simulate", "chsh-spad", "--deadtime-scan", "--rates", "2e6,6e6", "--duration-s", "0.02",
                 "--out", str(scan)]) == 0
    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "5", "--out", str(dark)]) == 0
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(dark), "--channel", "A+", "--out", str(out)]) == 0
    assert main(["calibrate", "deadtime", "chsh-spad", "--tags", str(scan), "--channel", "A+", "--dark", str(out),
                 "--out", str(out)]) == 0
    text = capsys.readouterr().out
    assert "deadtime-series" in text and "settings [0, 1]" in text
    # Council round 5 (V1): the artifact cites the dark artifact it used.
    arts = [load_calibration(f) for f in sorted(out.glob("*.json"))]
    dark_ids = [f.stem for f, a in zip(sorted(out.glob("*.json")), arts) if a.kind == "dark-afterpulse-deadtime"]
    series = next(a for a in arts if a.kind == KIND)
    assert len(dark_ids) == 1
    assert any(i[0] == dark_ids[0] and i[1] == "previous-artifact" for i in series.inputs)
    assert main(["simulate", "chsh-spad", "--rates", "1e6,2e6", "--out", str(tmp_path / "x.npz")]) != 0
    assert "--rates applies to --deadtime-scan" in capsys.readouterr().err


def test_coverage_shards_merge_to_the_whole(tmp_path):
    import json

    from dadl.calibration.coverage import main, merge, run

    whole = run("deadtime", 4, seed=3, regimes=["series-spad"])
    files = []
    for k in (1, 2):
        out = tmp_path / f"s{k}.json"
        assert main(["deadtime", "--regimes", "series-spad", "--replicates", "4", "--seed", "3", "--jobs", "2",
                     "--shard", f"{k}/2", "--out", str(out)]) == 0
        files.append(json.loads(out.read_text()))
    merged = merge(files)
    assert [(r.calibration, r.covered, r.at_least) for r in merged] == [(r.calibration, r.covered, r.at_least) for r in whole]
    with pytest.raises(ValueError, match="not all"):
        merge(files[:1])


# Review council, round 1.


def _offset(runs, step_ps):
    return [r + i * step_ps for i, r in enumerate(runs)]


def test_interruption_check_takes_each_run_on_its_own_bounds():
    det = _det(20_000, 50_000)
    rng = np.random.default_rng(5000)
    runs = [detect(poisson_times(r, t, rng), det, t, rng) for r, t in zip((2e6, 5e6, 1e7), (0.02, 0.02, 0.015))]
    aligned = calibrate_deadtime_series(runs, CH, backend="sim", seed=1, n_draws=100)
    # The same runs taken 30 ms apart on one running clock.
    offset = calibrate_deadtime_series(_offset(runs, 30_000_000_000), CH, backend="sim", seed=1, n_draws=100)
    assert aligned.status == "ok" and offset.status == "ok"
    assert offset.ppc_pvalue == pytest.approx(aligned.ppc_pvalue)
    # A 20% hole cut out of one run is still an interruption.
    t = runs[1]
    holed = [runs[0], t[(t < 0.4 * t[-1]) | (t > 0.6 * t[-1])], runs[2]]
    cut = calibrate_deadtime_series(holed, CH, backend="sim", seed=1, n_draws=100)
    assert cut.status == "misfit" and any("far longer than its count rate allows" in n for n in cut.notes)


def _survival_exact(w: float, n: float, tau: float) -> float:
    import math
    from decimal import Decimal, localcontext

    with localcontext() as ctx:
        ctx.prec = 80
        a = Decimal(n) * (-Decimal(n) * Decimal(tau)).exp()
        w_, tau_ = Decimal(w), Decimal(tau)
        total, j = Decimal(0), 0
        while j * tau_ <= w_:
            total += (-a * (w_ - j * tau_)) ** j / Decimal(math.factorial(j))
            j += 1
        return float(total)


@pytest.mark.parametrize("x", [0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0])
def test_survival_tail_matches_the_exact_series(x):
    # Past 8 tau the second real mode has not decayed near n tau = 1 (the tail was 20-60% low).
    tau = 30_000.0
    n = x / tau
    w = np.array([1.0, 5.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0]) * tau
    exact = np.array([_survival_exact(v, n, tau) for v in w])
    assert np.max(np.abs(survival(w, n, tau) / exact - 1.0)) < 1e-6


def test_profile_is_smooth_with_afterpulsing():
    import dadl.calibration.deadtime as dtm

    p, tau_a = 0.05, 20_000.0
    runs = _runs(_det(20_000, 50_000, afterpulse_prob=p, afterpulse_tau_ps=tau_a), (2e6, 5e6, 1e7), 0.02, 777)
    gaps = [np.diff(r) for r in runs]
    m = min(g.min() for g in gaps)
    fit = dtm._Fit([dtm._Run(g.astype(float) - m, 1, tau_a) for g in gaps], p, tau_a)
    prof = np.array([fit.profile(t) for t in np.arange(18_000, 22_001, 100.0)])
    # The kernel's jump snapped to the grid made the profile jagged by units (sd 1.38).
    assert np.diff(prof, 2).std() < 0.05


def _afterpulse_artifact(p: float, tau_a: float, rel_sd: float):
    """A stand-in for a dark artifact: the fields `calibrate_deadtime_series` reads."""
    from types import SimpleNamespace

    from dadl.stats import Interval

    rng = np.random.default_rng(0)
    samples = np.column_stack([p * (1 + rel_sd * rng.standard_normal(400)),
                               tau_a * (1 + rel_sd * rng.standard_normal(400))])
    return SimpleNamespace(
        estimates={"afterpulse_prob": Interval(p, p * 0.9, p * 1.1), "afterpulse_tau_ps": Interval(tau_a, tau_a * 0.9, tau_a * 1.1)},
        stand_in=(), sample_params=("afterpulse_prob", "afterpulse_tau_ps"), samples=samples,
        kind="dark-afterpulse-deadtime", channels=(CH,), status="ok", ppc_pvalue=0.5, inputs=(),
    )


def _chi2_groups(counts: np.ndarray, prob: np.ndarray, min_expected: float = 20.0) -> float:
    """Pearson chi-square p of binned counts against bin probabilities, bins merged to >= 20 expected."""
    mu = prob * counts.sum()
    obs, exp, acc_o, acc_e = [], [], 0.0, 0.0
    for o, e in zip(counts, mu):
        acc_o += o
        acc_e += e
        if acc_e >= min_expected:
            obs.append(acc_o)
            exp.append(acc_e)
            acc_o = acc_e = 0.0
    obs[-1] += acc_o
    exp[-1] += acc_e
    o, e = np.array(obs), np.array(exp)
    return float(stats.chi2.sf(np.sum((o - e) ** 2 / e), o.size - 1))


STRESS = {"tau_p": 20_000, "tau_np": 50_000, "rate": 1e7, "p": 0.05, "tau_a": 200_000.0}


def test_pending_afterpulse_law_matches_the_simulator():
    # Council round 5 (I2): at 1e7 cps with p 5% and tau_a 200 ns, 6.8% of gaps start with an
    # afterpulse pending. The law with the pending term fits the simulated gaps; without it, not.
    import dadl.calibration.deadtime as dtm

    tp, tnp, rate, p, tau_a = STRESS.values()
    det = _det(tp, tnp, dark_rate_cps=0.0, afterpulse_prob=p, afterpulse_tau_ps=tau_a)
    g = np.diff(_runs(det, (rate,), 0.2, 5)[0]).astype(float) - tnp
    run = dtm._Run(g, 1, tau_a, float(tnp))
    plain = dtm._Run(g, 1, tau_a)
    n = rate / 1e12
    assert _chi2_groups(run.counts, run.probs(n, tp, p, tau_a)) > 1e-3
    assert _chi2_groups(plain.counts, plain.probs(n, tp, p, tau_a)) < 1e-6


def test_pending_fraction_matches_the_simulator():
    # The simulated fractions of gaps starting with a pending candidate (council round 5,
    # law_stress.txt and law_regime.txt: 2.3e7 and 1.3e7 gaps).
    import dadl.calibration.deadtime as dtm

    tp, tnp, rate, p, tau_a = STRESS.values()
    assert dtm.pending_fraction(rate / 1e12, tp, tnp, p, tau_a) == pytest.approx(0.0671, rel=0.03)
    assert dtm.pending_fraction(6.05e6 / 1e12, 40_000, 45_000, 0.02, 30_000.0) == pytest.approx(0.00284, rel=0.03)
    # The fit's own value (from the kernel it already evaluates) agrees with the fine grid.
    g = np.random.default_rng(0).exponential(1e5, 5_000)
    run = dtm._Run(g, 1, tau_a, float(tnp))
    x = run._afterpulse(rate / 1e12, tp, tau_a)[2]
    pi = dtm._pending(p, x, np.exp(-tnp / tau_a))
    assert pi == pytest.approx(dtm.pending_fraction(rate / 1e12, tp, tnp, p, tau_a), rel=1e-3)


def test_unmodelled_pending_afterpulses_are_noted():
    from dadl.calibration.coverage import DEADTIME_REGIMES
    from dadl.calibration.deadtime import UNMODELLED_Z

    tp, tnp, _, p, tau_a = STRESS.values()
    det = _det(tp, tnp, afterpulse_prob=p, afterpulse_tau_ps=tau_a)
    cal = calibrate_deadtime_series(_runs(det, (2e6, 5e6, 1e7), 0.09, 21), CH,
                                    afterpulse=_afterpulse_artifact(p, tau_a, 0.02), backend="sim", seed=1,
                                    n_draws=100)
    assert cal.conditions["gaps"] > 1e6
    assert cal.conditions["afterpulse_pending_max"] == pytest.approx(0.068, rel=0.05)
    r = cal.conditions["afterpulse_unmodelled"]
    assert r * np.sqrt(cal.conditions["gaps"]) > UNMODELLED_Z
    assert any("second pending afterpulse" in n for n in cal.notes)
    # The passive-quench coverage regime stays below the guard.
    reg = DEADTIME_REGIMES["passive-quench"].detector
    rng = np.random.default_rng(22)
    runs = [detect(poisson_times(r, 0.01, rng), reg, 0.01, rng) for r in DEADTIME_REGIMES["passive-quench"].rates_cps]
    pq = calibrate_deadtime_series(runs, CH, afterpulse=_afterpulse_artifact(reg.afterpulse_prob, reg.afterpulse_tau_ps, 0.02),
                                   backend="sim", seed=1, n_draws=100)
    assert pq.conditions["afterpulse_unmodelled"] * np.sqrt(pq.conditions["gaps"]) < UNMODELLED_Z
    assert not any("second pending afterpulse" in n for n in pq.notes)


@pytest.mark.parametrize(
    "tau_p,tau_np,rates,seed",
    [(10_000, 22_000, (5e5, 1e6, 2e6), 2), (50_000, 50_000, (2e6, 5e6, 1e7), 3)],
)
def test_tau_p_never_exceeds_tau_np(tau_p, tau_np, rates, seed):
    det = DetectorParams(efficiency=1.0, dark_rate_cps=500.0, dead_time_ps=tau_np, dead_time_p_ps=tau_p,
                         jitter_sigma_ps=0.0)
    rng = np.random.default_rng(seed)
    runs = [detect(poisson_times(r, 0.02, rng), det, 0.02, rng) for r in rates]
    cal = calibrate_deadtime_series(runs, CH, backend="sim", seed=1, n_draws=400)
    p_iv, np_iv = cal.estimates["dead_time_p_ps"], cal.estimates["dead_time_ps"]
    assert p_iv.value <= np_iv.value and p_iv.high <= np_iv.high
    assert np.all(cal.samples[:, 0] <= cal.samples[:, 1])


def test_afterpulse_uncertainty_widens_an_upper_limit():
    det = _det(0, 50_000, afterpulse_prob=0.03, afterpulse_tau_ps=30_000.0)
    runs = _runs(det, (2e6, 5e6, 1e7), 0.01, 27)  # tau_p MLE 0: the case that took no afterpulse uncertainty
    tight = calibrate_deadtime_series(runs, CH, afterpulse=_afterpulse_artifact(0.03, 30_000.0, 1e-9), backend="sim",
                                      seed=1, n_draws=100)
    loose = calibrate_deadtime_series(runs, CH, afterpulse=_afterpulse_artifact(0.03, 30_000.0, 0.5), backend="sim",
                                      seed=1, n_draws=100)
    assert loose.unidentified == ("dead_time_p_ps",) and tight.estimates["dead_time_p_ps"].value == 0.0
    assert any("afterpulse uncertainty adds" in n for n in loose.notes)
    assert loose.estimates["dead_time_p_ps"].high > tight.estimates["dead_time_p_ps"].high


def test_default_scan_rates_count_avalanches():
    from dadl.cli import example_path
    from dadl.config import load_setup
    from dadl.simsetup import deadtime_scan_rates

    setup = load_setup(example_path("chsh-spad"))
    rates = deadtime_scan_rates(setup)
    eff = min(d.efficiency for d in setup.detectors.values())
    tau = max(d.dead_time_ps for d in setup.detectors.values()) * 1e-12
    assert max(rates) * eff * tau == pytest.approx(0.5)


def test_coverage_merge_refuses_shards_with_different_gates(tmp_path):
    import json

    from dadl.calibration.coverage import main, merge

    files = []
    for k in (1, 2):
        out = tmp_path / f"s{k}.json"
        assert main(["deadtime", "--regimes", "series-spad", "--replicates", "2", "--seed", "3", "--shard", f"{k}/2",
                     "--out", str(out)]) == 0
        files.append(json.loads(out.read_text()))
    files[1]["results"][0]["gate"] = [0.5, 1.0]
    with pytest.raises(ValueError, match="disagree on the gate"):
        merge(files)


def test_survival_takes_an_integer_tau():
    # Council round 3 (N3): an int tau_p raised a casting error.
    w = np.array([1.0, 30_000.0, 90_000.0])
    assert np.array_equal(survival(w, 1e-6, 20_000), survival(w, 1e-6, 20_000.0))
    assert np.array_equal(gap_cdf(w, 1e-6, 20_000), gap_cdf(w, 1e-6, 20_000.0))


@pytest.mark.slow
@pytest.mark.parametrize("seed", [0, 1])
def test_a_run_past_the_count_rate_maximum_is_fitted_on_its_branch(seed):
    # Council round 3 (N2): the gap law depends on a run's rate through n exp(-n tau_p), which takes
    # each value twice. A 6e7 cps run (n tau_p 2.4) was fitted on the low branch (7e6 cps) and
    # tau_p was biased to tau_np: 4/4 intervals excluded 40 ns.
    det = _det(40_000, 45_000, dark_rate_cps=50_000.0, afterpulse_prob=0.02, afterpulse_tau_ps=30_000.0)
    rng = np.random.default_rng(seed)
    ap = calibrate_dark(simulate_dark(det, 1.0, rng), CH, backend="sim", seed=1, n_draws=400)
    runs = [detect(poisson_times(r, 0.02, rng), det, 0.02, rng) for r in (1e6, 3e6, 6e7)]
    cal = calibrate_deadtime_series(runs, CH, afterpulse=ap, seed=1, n_draws=10)
    assert cal.estimates["dead_time_p_ps"].contains(40_000.0)
    assert any("input rate 5" in n or "input rate 6" in n for n in cal.notes if n.startswith("run 2"))
    # Council round 4: the branch is decided once, and the note names the run past the maximum.
    assert any(n.startswith("run(s) [2]") and "count-rate maximum" in n for n in cal.notes)
    assert cal.conditions["branch_nll_gain"] > 0


def test_a_false_branch_switch_becomes_a_union():
    # Council round 5 (I1): on CI-seed passive-quench replicate 477 the fit past the maximum won by
    # 4.17 in NLL (margin 3) and switched to tau_p near 30 ns, excluding 40 ns. With margin 5 the
    # interval spans both fits.
    from dadl.calibration.coverage import DEADTIME_REGIMES, _children

    reg = DEADTIME_REGIMES["passive-quench"]
    child = list(_children(20261111, 478, range(477, 478), "passive-quench"))[0]
    sim_seq, fit_seq = child.spawn(2)
    rng = np.random.default_rng(sim_seq)
    seed = int(fit_seq.generate_state(1)[0])
    ap = calibrate_dark(simulate_dark(reg.detector, reg.dark_duration_s, rng), CH, backend="sim", seed=seed, n_draws=400)
    runs = [detect(poisson_times(r, reg.duration_s, rng), reg.detector, reg.duration_s, rng) for r in reg.rates_cps]
    cal = calibrate_deadtime_series(runs, CH, afterpulse=ap, backend="sim", seed=seed, n_draws=10)
    assert 3.0 < cal.conditions["branch_nll_gain"] <= 5.0
    assert any("fit about as well past the count-rate maximum" in n for n in cal.notes)
    assert cal.estimates["dead_time_p_ps"].contains(40_000.0)


def _passive_quench_scan(i: int):
    """The council round 4 scan (tau_p 40 ns, tau_np 45 ns, 1/3/6 Mcps), seeds [777, i, 1]."""
    det = _det(40_000, 45_000, dark_rate_cps=50_000.0, afterpulse_prob=0.02, afterpulse_tau_ps=30_000.0)
    rng = np.random.default_rng([777, i, 1])
    ap = calibrate_dark(simulate_dark(det, 1.0, rng), CH, backend="sim", seed=1, n_draws=400)
    runs = [detect(poisson_times(r, 0.02, rng), det, 0.02, rng) for r in (1e6, 3e6, 6e6)]
    return calibrate_deadtime_series(runs, CH, afterpulse=ap, backend="sim", seed=1, n_draws=10)


def _afterpulse_shift_ps(cal) -> float:
    note = next(n for n in cal.notes if n.startswith("afterpulse uncertainty adds"))
    return float(note.split()[3])


@pytest.mark.slow
def test_afterpulse_refits_hold_the_decided_branch():
    # Council round 4 (DT1): the sensitivity refits searched both branches and jumped to the
    # other minimum: a shift sd of 1.36e4 ps on this seed, and an interval of [16, 45] ns.
    cal = _passive_quench_scan(20)
    assert _afterpulse_shift_ps(cal) < 1_000.0
    assert cal.estimates["dead_time_p_ps"].low > 35_000.0


@pytest.mark.slow
@pytest.mark.parametrize("i", [36, 119])
def test_a_scan_below_the_maximum_stays_on_the_low_branch(i):
    # Council round 4 (DT2): with both branches searched, these scans put the 6 Mcps run at about
    # 9e7 cps with tau_p near 31 ns and a false pending-afterpulse note. The high branch
    # fits better by 2.1 and 1.5 in NLL, inside the margin: the point stays on the low branch and
    # the interval spans both fits.
    from dadl.calibration.deadtime import BRANCH_MARGIN

    cal = _passive_quench_scan(i)
    assert abs(cal.estimates["dead_time_p_ps"].value - 40_000.0) < 2_000.0
    assert 0 < cal.conditions["branch_nll_gain"] <= BRANCH_MARGIN
    assert any("fit about as well past the count-rate maximum" in n for n in cal.notes)
    assert not any("second pending afterpulse" in n for n in cal.notes)
    assert all("input rate 1e+08" not in n and "input rate 9" not in n for n in cal.notes if n.startswith("run 2"))


def test_series_refusals_of_inputs():
    # Council round 4 (S7).
    from types import SimpleNamespace

    runs = _runs(_det(20_000, 50_000), (2e6, 5e6), 0.01, 11)
    dup = np.sort(np.concatenate([runs[0], runs[0][:1]]))
    with pytest.raises(ValueError, match=r"duplicate timestamps in runs\[0\]"):
        calibrate_deadtime_series([dup, runs[1]], CH)
    with pytest.raises(ValueError, match="same acquisition"):
        calibrate_deadtime_series([runs[0], runs[0].copy()], CH)
    ap = _afterpulse_artifact(0.02, 30_000.0, 0.02)
    with pytest.raises(ValueError, match="dark artifact"):
        calibrate_deadtime_series(runs, CH, afterpulse=SimpleNamespace(**{**vars(ap), "kind": "delay-jitter"}))
    with pytest.raises(ValueError, match="for channel"):
        calibrate_deadtime_series(runs, CH, afterpulse=SimpleNamespace(**{**vars(ap), "channels": (ChannelId("sim", 1),)}))


def test_a_misfit_afterpulse_artifact_is_named():
    from types import SimpleNamespace

    runs = _runs(_det(20_000, 50_000, afterpulse_prob=0.02, afterpulse_tau_ps=30_000.0), (2e6, 5e6), 0.01, 11)
    ap = SimpleNamespace(**{**vars(_afterpulse_artifact(0.02, 30_000.0, 0.02)), "status": "misfit", "ppc_pvalue": 1e-4})
    cal = calibrate_deadtime_series(runs, CH, afterpulse=ap, backend="sim", seed=1, n_draws=10)
    assert any("afterpulse artifact for" in n and "misfit" in n for n in cal.notes)


def test_without_afterpulsing_a_rate_past_the_maximum_is_ambiguous():
    # Without afterpulsing both branches fit equally: the note gives both rates.
    runs = _runs(_det(20_000, 50_000), (5e6, 2e7, 1e8), 0.01, 3)
    cal = calibrate_deadtime_series(runs, CH, seed=1, n_draws=10)
    note = next(n for n in cal.notes if n.startswith("run 2"))
    assert " or " in note and "do not tell which" in note
