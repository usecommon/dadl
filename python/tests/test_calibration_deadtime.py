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
