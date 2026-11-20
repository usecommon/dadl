"""The calibration coverage gate's verdicts, regimes and shard merge (council round 7, Q1 and S1).

The 1,000-replicate runs are in CI; these check the machinery that turns counts into a verdict.
"""

from __future__ import annotations

import types

import pytest

from dadl.calibration import coverage
from dadl.calibration.coverage import (
    DEADTIME_REGIMES,
    MISFITS,
    CoverageResult,
    expected_results,
    gate,
    merge,
    shard_record,
)
from dadl.stats import Interval


def test_at_least_gates_only_the_lower_edge():
    band = gate(1000, 2)  # (0.924, 0.972)
    assert not CoverageResult("x", 1000, {"k": 923}, band, at_least=True).passed()
    assert CoverageResult("x", 1000, {"k": 949}, band, at_least=True).passed()
    assert CoverageResult("x", 1000, {"k": 1000}, band, at_least=True).passed()
    # Two-sided, 100% fails: the at-least rule is what lets it pass.
    assert not CoverageResult("x", 1000, {"k": 1000}, band).passed()
    assert CoverageResult("x", 1000, {"k": 950, MISFITS: 10}, band).passed()
    assert not CoverageResult("x", 1000, {"k": 950, MISFITS: 11}, band).passed()


def test_the_dead_time_results_and_their_rules():
    p, np_ = coverage.deadtime_coverage(1000, regime="series-spad", part=range(0))
    assert p.calibration == "deadtime-series [series-spad]" and not p.at_least
    assert np_.calibration == "deadtime-series dead_time_ps [series-spad]" and np_.at_least
    assert p.gate == np_.gate == gate(1000, 2)
    assert set(p.covered) == {"dead_time_p_ps", MISFITS} and set(np_.covered) == {"dead_time_ps"}


def test_the_dead_time_regimes():
    # Pinned: the truths and scans the gate's figures refer to (round 6: 95.5/99.9 and 94.6/99.0).
    table = {
        name: (r.detector.dead_time_p_ps, r.detector.dead_time_ps, r.detector.afterpulse_prob, r.rates_cps,
               r.duration_s, r.dark_duration_s)
        for name, r in DEADTIME_REGIMES.items()
    }
    assert table == {
        "series-spad": (20_000, 50_000, 0.0, (2e6, 5e6, 1e7), 0.02, 0.0),
        "passive-quench": (40_000, 45_000, 0.02, (1e6, 3e6, 6e6), 0.02, 1.0),
    }


def _fake_runs(monkeypatch, tau_np_interval):
    """Stub the calibrations so a replicate is cheap; record what the dead-time fit was given."""
    seen = []
    dark = object()

    def fake_series(runs, channel, *, afterpulse=None, **_):
        seen.append(afterpulse)
        est = {"dead_time_p_ps": Interval(20_000.0, 19_000.0, 21_000.0), "dead_time_ps": tau_np_interval}
        return types.SimpleNamespace(estimates=est, status="ok")

    monkeypatch.setattr(coverage, "calibrate_deadtime_series", fake_series)
    monkeypatch.setattr(coverage, "calibrate_dark", lambda *a, **k: dark)
    monkeypatch.setattr(coverage, "simulate_dark", lambda *a, **k: None)
    monkeypatch.setattr(coverage, "poisson_times", lambda *a, **k: None)
    monkeypatch.setattr(coverage, "detect", lambda *a, **k: None)
    return seen, dark


def test_a_planted_tau_np_miss_is_counted(monkeypatch):
    _fake_runs(monkeypatch, Interval(48_000.0, 47_000.0, 49_000.0))
    p, np_ = coverage.deadtime_coverage(3, regime="series-spad")
    assert p.covered["dead_time_p_ps"] == 3 and np_.covered["dead_time_ps"] == 0
    _fake_runs(monkeypatch, Interval(50_000.0, 49_000.0, 51_000.0))
    assert coverage.deadtime_coverage(3, regime="series-spad")[1].covered["dead_time_ps"] == 3


def test_passive_quench_holds_afterpulsing_at_a_dark_artifact(monkeypatch):
    seen, dark = _fake_runs(monkeypatch, Interval(45_000.0, 44_000.0, 46_000.0))
    coverage.deadtime_coverage(2, regime="passive-quench")
    assert seen == [dark, dark]
    seen.clear()
    coverage.deadtime_coverage(2, regime="series-spad")
    assert seen == [None, None]


def _shards(ranges, replicates=10):
    exp = expected_results("deadtime", ["series-spad"], replicates, 3)
    records = []
    for lo, hi in ranges:
        res = [CoverageResult(r.calibration, replicates, dict.fromkeys(r.covered, 0), r.gate, r.at_least)
               for r in exp.values()]
        records.append(shard_record(res, 3, range(lo, hi), "deadtime", ["series-spad"]))
    return records


def test_merge_tiles_the_replicates():
    merged = merge(_shards([(0, 5), (5, 10)]))
    assert [r.calibration for r in merged] == list(expected_results("deadtime", ["series-spad"], 10, 3))
    assert merged[1].at_least and merged[0].gate == gate(10, 2)
    with pytest.raises(ValueError, match="gap or overlap at 6"):
        merge(_shards([(0, 5), (6, 10)]))
    with pytest.raises(ValueError, match="gap or overlap at 4"):
        merge(_shards([(0, 6), (4, 10)]))
    with pytest.raises(ValueError, match="not all"):
        merge(_shards([(0, 5)]))


def test_merge_takes_the_gate_from_the_code_and_refuses_empty_or_partial_sets():
    # Council round 7 (S1): an empty shard set and one with a 0-100000% band both passed.
    with pytest.raises(ValueError, match="family and regimes"):
        merge([{"seed": 20261111, "replicates": 1000, "range": [0, 1000], "results": []}])
    good = _shards([(0, 10)])
    with pytest.raises(ValueError, match="no results"):
        merge([{**good[0], "results": []}])
    lax = [{**good[0], "results": [{**r, "gate": [0, 1000]} for r in good[0]["results"]]}]
    with pytest.raises(ValueError, match="not the code's"):
        merge(lax)
    with pytest.raises(ValueError, match="missing"):
        merge([{**good[0], "results": good[0]["results"][:1]}])
    with pytest.raises(ValueError, match="missing"):
        merge([{**good[0], "regimes": None}])
    with pytest.raises(ValueError, match="count"):
        merge([{**good[0], "results": [{**r, "covered": {}} for r in good[0]["results"]]}])
    with pytest.raises(ValueError, match="family or regimes"):
        merge([good[0], {**good[0], "family": "dark"}])


def test_the_reviewers_shard_files_are_refused(tmp_path):
    import json

    from dadl.analysis import propagation

    empty = {"seed": 20261111, "replicates": 1000, "range": [0, 1000], "results": []}
    lax = {"seed": 20261111, "replicates": 1000, "range": [0, 1000], "results": [
        {"calibration": "deadtime-series [series-spad]", "covered": {"dead_time_p_ps": 10, MISFITS: 0},
         "gate": [0, 1000], "at_least": False}]}
    for record in (empty, lax):
        path = tmp_path / "shard.json"
        path.write_text(json.dumps(record))
        with pytest.raises(ValueError):
            coverage.main(["merge", str(path)])
        with pytest.raises(ValueError):
            propagation.main(["merge", str(path)])
