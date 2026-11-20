"""Calibration sets (phase 3, M1): loading, refusals, labels and joint draws."""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from dadl.artifact import write_artifact
from dadl.calibration.base import Calibration
from dadl.calibration.dark import calibrate_dark
from dadl.data import ChannelId, tags_sha256
from dadl.errors import ConfigError, GuardrailRefusal, StaleCalibration
from dadl.guardrails import POLICY_MODIFIED, Policy
from dadl.propagate import (
    BOUND_ONLY,
    MISFIT_ACCEPTED,
    NO_VALIDITY,
    NOT_CHECKED,
    OUTSIDE_VALIDITY,
    SERIES_UNUSED,
    CalibrationSet,
    DatasetContext,
    Member,
    condition_mismatches,
)
from dadl.sim import DetectorParams, simulate_dark
from dadl.stats import Interval

UTC = dt.timezone.utc
CREATED = dt.datetime(2026, 11, 1, tzinfo=UTC)
A, B = ChannelId("sim", 1), ChannelId("sim", 2)
ANALYSIS = "hom_visibility"


def make_cal(channel=A, *, seed=0, rho=0.8, ppc=0.5, conditions=None, stand_in=(), inputs=None, kind="dark-afterpulse-deadtime"):
    """A calibration with two correlated posterior columns (rate, p) and known draws."""
    rng = np.random.default_rng(seed)
    z = rng.multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=2000)
    rate = 1000 + 10 * z[:, 0]
    p = np.clip(0.05 + 0.002 * z[:, 1], 0, 1)
    if "p" in stand_in:
        p = rng.uniform(0, 0.01, 2000)
    est = {
        "rate": Interval(float(np.median(rate)), *np.quantile(rate, [0.025, 0.975])),
        "p": Interval(float(np.median(p)), 0.0 if stand_in else float(np.quantile(p, 0.025)), float(np.quantile(p, 0.975))),
    }
    return Calibration(
        kind=kind,
        method="test",
        method_version="1",
        channels=(channel,),
        estimates=est,
        units={"rate": "s-1", "p": "1"},
        samples=np.column_stack([rate, p]),
        sample_params=("rate", "p"),
        converged=True,
        ppc_pvalue=ppc,
        inputs=inputs or ((f"{seed:064x}", "raw-tags", str(channel)),),
        backend="sim",
        conditions=conditions or {},
        unidentified=tuple(stand_in),
        stand_in=tuple(stand_in),
    )


def write(cal, d, **kw):
    kw.setdefault("created", CREATED)
    return write_artifact(cal, d, accept_misfit=True, **kw)


def ctx(**kw):
    kw.setdefault("acquired", dt.datetime(2026, 11, 2, tzinfo=UTC))
    return DatasetContext(**kw)


# -- loading and selection ----------------------------------------------------------------------


def test_loads_a_directory_and_selects_by_channel_and_kind(tmp_path):
    write(make_cal(A, seed=1), tmp_path)
    write(make_cal(B, seed=2), tmp_path)
    write(make_cal(A, seed=3, kind="delay-jitter"), tmp_path)
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), channels=[A], kinds=["dark-afterpulse-deadtime"])
    assert len(s.members) == 1 and s.get("dark-afterpulse-deadtime", A).channels == (A,)
    assert len(s.artifact_ids[0]) == 64
    assert s.labels == (NO_VALIDITY.format(aid=s.artifact_ids[0][:12]),)


def test_two_artifacts_of_one_kind_and_channel_are_ambiguous(tmp_path):
    write(make_cal(A, seed=1), tmp_path)
    write(make_cal(A, seed=2), tmp_path)
    with pytest.raises(ConfigError, match="pass the one to use"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())


# -- joint draws ---------------------------------------------------------------------------------


def test_columns_of_one_artifact_share_a_sample_index(tmp_path):
    cal = make_cal(A, seed=1, rho=0.9)
    write(cal, tmp_path)
    d = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx()).draws(1500, seed=7)
    rate, p = d.column("dark-afterpulse-deadtime", "rate", A), d.column("dark-afterpulse-deadtime", "p", A)
    # every drawn row is a stored row, so the joint structure is kept exactly
    stored = {tuple(r) for r in cal.samples}
    assert all((r, q) in stored for r, q in zip(rate, p))
    assert np.corrcoef(rate, p)[0, 1] == pytest.approx(np.corrcoef(cal.samples.T)[0, 1], abs=0.03)
    assert len(set(rate)) == 1500  # without replacement while n <= stored draws


def test_different_artifacts_are_drawn_independently(tmp_path):
    # Two channels whose stored draws came from the same generator (identical columns): drawing
    # them with a shared index would make the difference of the rates exactly zero (cycle 3.1).
    ca, cb = make_cal(A, seed=5), make_cal(B, seed=5)
    assert np.array_equal(ca.samples, cb.samples)
    write(ca, tmp_path)
    write(cb, tmp_path)
    d = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx()).draws(2000, seed=3)
    ra, rb = d.column("dark-afterpulse-deadtime", "rate", A), d.column("dark-afterpulse-deadtime", "rate", B)
    assert abs(np.corrcoef(ra, rb)[0, 1]) < 0.1
    assert np.std(ra - rb) == pytest.approx(np.sqrt(2) * np.std(ra), rel=0.1)


def test_draws_are_reproducible_and_seed_dependent(tmp_path):
    write(make_cal(A, seed=1), tmp_path)
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    a, b, c = s.draws(100, 1), s.draws(100, 1), s.draws(100, 2)
    col = lambda d: d.column("dark-afterpulse-deadtime", "rate", A)  # noqa: E731
    assert np.array_equal(col(a), col(b))
    assert not np.array_equal(col(a), col(c))


def test_more_draws_than_stored_resamples_with_replacement(tmp_path):
    write(make_cal(A, seed=1), tmp_path)
    d = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx()).draws(5000, seed=1)
    assert d.column("dark-afterpulse-deadtime", "rate", A).shape == (5000,)


def test_stand_in_columns_are_bounds_not_draws(tmp_path):
    write(make_cal(A, seed=1, stand_in=("p",)), tmp_path)
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    d = s.draws(100, seed=1)
    with pytest.raises(KeyError):
        d.column("dark-afterpulse-deadtime", "p", A)
    lo, hi = d.bound("dark-afterpulse-deadtime", "p", A)
    assert lo == 0 and hi == pytest.approx(s.get("dark-afterpulse-deadtime").estimates["p"].high)
    # The forward model says which parameter is a bound (council round 2); the set does not.
    assert BOUND_ONLY not in s.labels


# -- misfits -------------------------------------------------------------------------------------


def test_misfit_artifact_is_refused(tmp_path):
    write(make_cal(A, seed=1, ppc=1e-6), tmp_path)
    with pytest.raises(GuardrailRefusal, match="misfit"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())


def test_misfit_artifact_under_override_is_labelled(tmp_path):
    write(make_cal(A, seed=1, ppc=1e-6), tmp_path)
    policy = Policy.default().allow(ANALYSIS, "misfit_calibration", "only run available; re-acquiring Monday")
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), policy=policy)
    assert s.labels[0] == POLICY_MODIFIED and MISFIT_ACCEPTED in s.labels


def test_override_for_another_analysis_does_not_apply(tmp_path):
    write(make_cal(A, seed=1, ppc=1e-6), tmp_path)
    policy = Policy.default().allow("fringe_visibility", "misfit_calibration", "other analysis")
    with pytest.raises(GuardrailRefusal):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), policy=policy)


# -- disjointness --------------------------------------------------------------------------------


def test_calibration_on_the_analysed_tags_is_refused_without_override(tmp_path):
    tags = simulate_dark(DetectorParams(dark_rate_cps=1000.0, afterpulse_prob=0.0), 3.0, np.random.default_rng(1))
    cal = calibrate_dark(tags, A, backend="sim", seed=1, n_draws=200)
    write(cal, tmp_path)
    same = ctx(input_sha256=frozenset({tags_sha256(tags)}))
    with pytest.raises(GuardrailRefusal, match="disjoint"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=same)
    loose = Policy.default().allow(ANALYSIS, "calibration_outside_validity", "x").allow(ANALYSIS, "misfit_calibration", "x")
    with pytest.raises(GuardrailRefusal, match="disjoint"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=same, policy=loose)
    other = simulate_dark(DetectorParams(dark_rate_cps=1000.0), 1.0, np.random.default_rng(2))
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(input_sha256=frozenset({tags_sha256(other)})))
    assert s.get("dark-afterpulse-deadtime").estimates["dark_rate_cps"].contains(1000.0)


def test_shared_raw_file_hash_is_refused(tmp_path):
    h = "ab" * 32
    write(make_cal(A, seed=1, inputs=(("00" * 32, "raw-tags", "a"), (h, "raw-file", "setting0.ptu"))), tmp_path)
    with pytest.raises(GuardrailRefusal):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(input_sha256=frozenset({h})))


# -- validity range and conditions --------------------------------------------------------------


def test_outside_the_validity_range_is_stale(tmp_path):
    write(make_cal(A, seed=1), tmp_path, valid_until=dt.datetime(2026, 11, 3, tzinfo=UTC))
    CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    late = ctx(acquired=dt.datetime(2026, 11, 4, tzinfo=UTC))
    with pytest.raises(StaleCalibration, match="outside"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=late)
    with pytest.raises(StaleCalibration, match="acquisition time was not given"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=DatasetContext())


def test_outside_validity_under_override_is_labelled(tmp_path):
    write(make_cal(A, seed=1), tmp_path, valid_until=dt.datetime(2026, 11, 3, tzinfo=UTC))
    policy = Policy.default().allow(ANALYSIS, "calibration_outside_validity", "recalibrated next day, same settings")
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(acquired=dt.datetime(2026, 11, 4, tzinfo=UTC)), policy=policy)
    assert POLICY_MODIFIED in s.labels and OUTSIDE_VALIDITY in s.labels


def test_tagger_settings_must_match(tmp_path):
    write(make_cal(A, seed=1, conditions={"tagger_HWSync_Offset": 0.0, "tagger_HW_InpChan_TrgLevel": -100.0}), tmp_path)
    same = {"tagger_HWSync_Offset": 0.0, "tagger_HW_InpChan_TrgLevel": -100.0}
    CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(conditions=same))
    with pytest.raises(StaleCalibration, match="TrgLevel"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(conditions={**same, "tagger_HW_InpChan_TrgLevel": -120.0}))
    with pytest.raises(StaleCalibration, match="not given"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())


def test_condition_rules():
    assert condition_mismatches({"bias_V": 25.0}, {"bias_V": 25.0}) == []
    assert condition_mismatches({"bias_V": 25.0}, {"bias_V": 26.0})
    rng = {"temperature_min_K": 250.0, "temperature_max_K": 290.0}
    assert condition_mismatches(rng, {"temperature_K": 270.0}) == []
    assert condition_mismatches(rng, {"temperature_K": 300.0})
    assert condition_mismatches({"temperature_K": 253.0}, {"temperature_K": 253.4}) == []
    assert condition_mismatches({"temperature_K": 253.0}, {"temperature_K": 254.0})
    # diagnostics are not compared; a recorded temperature or bracket the dataset does not give is
    # a mismatch, as a bias is (council round 3: leaving the condition out passed the check)
    assert condition_mismatches({"gaps": 2e4}, {}) == []
    assert condition_mismatches({"temperature_K": 253.0}, {}) == ["temperature_K = 253 recorded, not given for the dataset"]
    assert condition_mismatches(rng, {}) == ["temperature_min_K = 250 recorded, not given for the dataset"]


def test_without_a_dataset_only_the_misfit_check_runs(tmp_path):
    write(make_cal(A, seed=1, conditions={"tagger_HWSync_Offset": 0.0}), tmp_path, valid_until=dt.datetime(2026, 11, 3, tzinfo=UTC))
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS)
    assert s.labels == (NOT_CHECKED,) and not s.dataset_checked


# -- council round 2 ------------------------------------------------------------------------------


def test_validity_starts_at_acquisition(tmp_path):
    # Calibration data taken 09:00, the artifact written 09:30 and valid to 18:00: a dataset taken
    # at 09:10 is covered (validity started at 09:30 before; council round 2).
    day = dt.datetime(2026, 11, 2, tzinfo=UTC)
    write(
        make_cal(A, seed=1), tmp_path, created=day.replace(hour=9, minute=30), acquired=day.replace(hour=9),
        valid_until=day.replace(hour=18),
    )
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(acquired=day.replace(hour=9, minute=10)))
    assert s.members[0].manifest["valid_for"]["from"].startswith("2026-11-02T09:00")
    with pytest.raises(StaleCalibration, match="outside"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(acquired=day.replace(hour=8, minute=50)))


def test_deadtime_series_artifacts_are_skipped_where_unused(tmp_path):
    write(make_cal(A, seed=1), tmp_path)
    write(make_cal(A, seed=2, kind="deadtime-series"), tmp_path)
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    assert [m.kind for m in s.members] == ["dark-afterpulse-deadtime"] and SERIES_UNUSED in s.labels
    # Asked for by kind, it is kept.
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), kinds=("deadtime-series",))
    assert [m.kind for m in s.members] == ["deadtime-series"] and SERIES_UNUSED not in s.labels


def test_a_set_loaded_without_a_dataset_labels_the_corrected_value(tmp_path):
    # Council round 2 (R4): a delay calibrated on the analysed tags, loaded without a
    # DatasetContext, gave a corrected value with neither refusal nor label.
    from dadl.analysis.chsh import Window
    from dadl.analysis.hom import hom_from_tags
    from dadl.calibration.delay import calibrate_delay
    from dadl.sim import PairSource, simulate_hom

    rng = np.random.default_rng(2)
    det = DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, dead_time_ps=30_000, jitter_sigma_ps=30.0)
    src = PairSource(pair_rate_per_s=5e5, visibility=0.9, coherence_time_ps=2.0)
    runs = {k: simulate_hom(src, det, det, d, 0.2, rng) for k, d in (("dip", 0.0), ("reference", 20.0))}
    ref = runs["reference"]
    write(calibrate_delay(np.sort(ref["A"]), np.sort(ref["B"]), (A, B), backend="sim", seed=1), tmp_path)
    policy = Policy.default().allow(ANALYSIS, "window_capture", "test")
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, policy=policy)
    r, _ = hom_from_tags(
        runs, Window(width_ps=1000, source="plan"), policy=policy, calibrations=s, channels=(A, B), reference_delay_ps=20.0
    )
    assert "window_capture" in r.corrections and NOT_CHECKED in r.corrected_labels


# -- council round 3 ------------------------------------------------------------------------------


def test_a_recorded_temperature_needs_the_dataset_s(tmp_path):
    # R2: dark artifacts recorded at 2.5 K loaded with no --condition and no label.
    write(make_cal(A, seed=1, conditions={"temperature_K": 2.5}), tmp_path)
    with pytest.raises(StaleCalibration, match="temperature_K = 2.5 recorded, not given for the dataset"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(conditions={"temperature_K": 2.6}))
    assert s.dataset_record() == {"acquired": "2026-11-02T00:00:00+00:00", "conditions": {"temperature_K": 2.6}}


def test_an_artifact_without_validity_is_labelled(tmp_path):
    write(make_cal(A, seed=1), tmp_path)
    write(make_cal(B, seed=2), tmp_path, valid_until=dt.datetime(2026, 11, 3, tzinfo=UTC))
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    no_range = [m.artifact_id[:12] for m in s.members if "valid_for" not in m.manifest]
    assert s.labels == (NO_VALIDITY.format(aid=no_range[0]),)


def test_ambiguity_does_not_depend_on_channel_order(tmp_path):
    # R4: one pair calibrated as (A, B) and again as (B, A); both were loaded, and the artifact ID
    # (so --created) chose which one was used.
    import dataclasses

    write(dataclasses.replace(make_cal(A, seed=1), channels=(A, B)), tmp_path)
    write(dataclasses.replace(make_cal(A, seed=2), channels=(B, A)), tmp_path)
    with pytest.raises(ConfigError, match="pass the one to use"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())


def test_another_detector_s_artifact_is_outside_its_validity(tmp_path):
    # R1: an artifact recorded for detector SNSPD-X9 was used for SNSPD-A1 on the same tagger input
    # (ChannelId equality is the input), getting around the refusal of the real detector's.
    old = ChannelId("HydraHarp-100", 0, role="A", detector="SNSPD-X9")
    new = ChannelId("HydraHarp-100", 0, role="A", detector="SNSPD-A1")
    write(make_cal(old, seed=1), tmp_path, valid_until=dt.datetime(2026, 11, 3, tzinfo=UTC))
    with pytest.raises(StaleCalibration, match="calibrated on detector SNSPD-X9, analysing SNSPD-A1"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), channels=[new])
    CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), channels=[old])
    CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), channels=[ChannelId("HydraHarp-100", 0)])
    policy = Policy.default().allow(ANALYSIS, "calibration_outside_validity", "same detector, renamed")
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(), channels=[new], policy=policy)
    assert OUTSIDE_VALIDITY in s.labels


def test_a_dead_time_longer_than_the_analysed_gaps_is_refused():
    # R1: another detector's dark artifact applied a 1 us dead time to a run whose gaps go down to
    # 30 ns.
    import dataclasses

    from dadl.forward import check_dead_time_gaps

    det = DetectorParams(efficiency=1.0, dark_rate_cps=50_000.0, dead_time_ps=1_000_000, jitter_sigma_ps=0.0)
    rng = np.random.default_rng(5)
    cal = calibrate_dark(simulate_dark(det, 2.0, rng), A, backend="sim", seed=1, n_draws=400)
    assert "dead_time_ps" not in cal.stand_in
    s = CalibrationSet(ANALYSIS, (Member("x" * 64, cal, {}),), ())
    fast = dataclasses.replace(det, dead_time_ps=30_000, dark_rate_cps=2e6)
    run = {"A": simulate_dark(fast, 0.05, rng)}
    with pytest.raises(StaleCalibration, match="analysed gaps are shorter, down to"):
        check_dead_time_gaps(s, {"A": A}, [run], analysis=ANALYSIS)
    policy = Policy.default().allow(ANALYSIS, "dead_time_contradicted", "checking")
    labels, counts = check_dead_time_gaps(s, {"A": A}, [run], analysis=ANALYSIS, policy=policy)
    rest = tuple(lab for lab in labels if not lab.startswith("the dark dead time of"))
    assert len(rest) == 1 and rest[0].startswith(CONTRADICTED) and 29_000 < counts["min_gap_ps"]["A"] < 31_000
    slow = {"A": simulate_dark(det, 1.0, rng)}
    assert check_dead_time_gaps(s, {"A": A}, [slow], analysis=ANALYSIS)[0] == ()


CONTRADICTED = "dead time contradicted by the analysed gaps, accepted: "


def _gap_check(cal_tau: int, true_tau: int, rate: float, rng, policy=None):
    """Calibrate a 1 s, 50 kcps dark run at cal_tau; check it on 0.05 s at rate with true_tau."""
    import dataclasses

    from dadl.forward import check_dead_time_gaps

    base = DetectorParams(efficiency=1.0, dark_rate_cps=50_000.0, dead_time_ps=cal_tau, jitter_sigma_ps=0.0)
    cal = calibrate_dark(simulate_dark(base, 1.0, rng), A, backend="sim", seed=1, n_draws=50)
    run = {"A": simulate_dark(dataclasses.replace(base, dead_time_ps=true_tau, dark_rate_cps=rate), 0.05, rng)}
    s = CalibrationSet(ANALYSIS, (Member("x" * 64, cal, {}),), ())
    return check_dead_time_gaps(s, {"A": A}, [run], analysis=ANALYSIS, policy=policy)


def test_the_gap_check_refuses_in_both_directions():
    # Council round 4 (DT4): the 0.9 factor passed a dead time up to about 16% too long, and a
    # dead time shorter than the detector's was never checked (22 ns on a 50 ns detector).
    rng = np.random.default_rng(3)
    with pytest.raises(StaleCalibration, match="longer than its upper end by"):
        _gap_check(22_000, 50_000, 1e6, rng)
    with pytest.raises(StaleCalibration, match="analysed gaps are shorter"):
        _gap_check(50_000, 46_000, 5e6, rng)
    for rate in (1e6, 5e6):
        labels, counts = _gap_check(50_000, 50_000, rate, rng)
        assert tuple(lab for lab in labels if not lab.startswith("the dark dead time of")) == () and 40_000 < counts["gap_check_floor_ps"]["A"] < 50_000
    policy = Policy.default().allow(ANALYSIS, "dead_time_contradicted", "checking")
    labels = _gap_check(22_000, 50_000, 1e6, rng, policy)[0]
    assert len(labels) == 1 and labels[0].startswith(CONTRADICTED) and "longer than its upper end" in labels[0]


def test_the_upward_gap_test_applies_to_an_upper_limit():
    from dadl.forward import check_dead_time_gaps

    import dataclasses

    det = DetectorParams(efficiency=1.0, dark_rate_cps=50_000.0, dead_time_ps=50_000, jitter_sigma_ps=0.0)
    cal = calibrate_dark(simulate_dark(det, 0.2, np.random.default_rng(4)), A, backend="sim", seed=1, n_draws=50)
    bound = dataclasses.replace(cal, unidentified=("dead_time_ps",), stand_in=("dead_time_ps",))
    s = CalibrationSet(ANALYSIS, (Member("x" * 64, bound, {}),), ())
    hi = bound.estimates["dead_time_ps"].high
    rng = np.random.default_rng(4)
    run = {"A": np.cumsum(int(hi) + 1_000_000 + rng.integers(0, 10_000, 5_000)).astype(np.int64)}
    with pytest.raises(StaleCalibration, match="dead time below"):
        check_dead_time_gaps(s, {"A": A}, [run], analysis=ANALYSIS)
    near = {"A": np.cumsum(int(hi) + rng.integers(1, 10_000, 5_000)).astype(np.int64)}
    labels, counts = check_dead_time_gaps(s, {"A": A}, [near], analysis=ANALYSIS)
    assert labels == () and counts["gap_check_floor_ps"]["A"] is None


def test_validity_without_from_starts_at_creation():
    # Council round 4 (S8): the schema allows valid_for without from (or to); _check raised KeyError.
    from dadl.propagate import _check

    cal = make_cal(A, seed=1)
    member = Member("x" * 64, cal, {"created": "2026-11-01T00:00:00Z", "valid_for": {"to": "2026-11-03T00:00:00Z"}})
    assert _check(member, ANALYSIS, Policy.default(), ctx(), {}) == []
    early = DatasetContext(acquired=dt.datetime(2026, 10, 1, tzinfo=UTC))
    with pytest.raises(StaleCalibration, match="outside 2026-11-01T00:00:00Z to 2026-11-03"):
        _check(member, ANALYSIS, Policy.default(), early, {})
    open_end = Member("x" * 64, cal, {"created": "2026-11-01T00:00:00Z", "valid_for": {"conditions": {}}})
    assert _check(open_end, ANALYSIS, Policy.default(), ctx(), {}) == []


def test_valid_for_conditions_are_enforced(tmp_path):
    # Council round 5 (V5): valid_for.conditions validated against the schema but were ignored.
    import json

    from dadl.artifact import artifact_id, canonical_json

    path = write(make_cal(A, seed=1), tmp_path / "orig")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["valid_for"] = {**manifest.get("valid_for", {}), "conditions": {"temperature_K": 300.0, "bias_V": 25.0}}
    aid = artifact_id(manifest)
    d = tmp_path / "mut"
    d.mkdir()
    (d / f"{aid}.json").write_bytes(canonical_json(manifest))
    (d / f"{aid}.posterior.npy").write_bytes(path.with_name(f"{path.stem}.posterior.npy").read_bytes())
    cold = ctx(conditions={"temperature_K": 250.0, "bias_V": 5.0})
    with pytest.raises(StaleCalibration, match="temperature_K = 300 recorded, 250"):
        CalibrationSet.load(d, analysis=ANALYSIS, dataset=cold)
    CalibrationSet.load(d, analysis=ANALYSIS, dataset=ctx(conditions={"temperature_K": 300.0, "bias_V": 25.0}))
    policy = Policy.default().allow(ANALYSIS, "calibration_outside_validity", "same detector, recalibrated at 250 K")
    s = CalibrationSet.load(d, analysis=ANALYSIS, dataset=cold, policy=policy)
    assert OUTSIDE_VALIDITY in s.labels and POLICY_MODIFIED in s.labels


# Council round 7.


def test_the_t3_sync_period_is_compared_to_1e_4():
    # G2: two acquisitions 10 Hz apart at 80 MHz were refused, and %g printed both as 1.25e-08.
    key = "tagger_MeasDesc_GlobalResolution"
    assert condition_mismatches({key: 1 / 79_999_870.0}, {key: 1 / 79_999_860.0}) == []
    (msg,) = condition_mismatches({key: 1 / 79_999_870.0}, {key: 1 / 76e6})
    assert "1.250002031e-08 recorded, 1.315789474e-08 for the dataset" in msg
    # Other settings stay exact.
    assert condition_mismatches({"tagger_HWSync_Offset": 1.0}, {"tagger_HWSync_Offset": 1.00001})


def test_a_delay_artifact_from_another_tagger_is_refused(tmp_path):
    # G6: the tagger's channel skew belongs to the instrument.
    from dadl.propagate import OTHER_TAGGER

    delay = make_cal(A, seed=1, kind="delay-jitter")
    write(delay, tmp_path, instrument={"model": "HydraHarp", "serial": "1111111"})
    other = ctx(instrument={"model": "HydraHarp", "serial": "2222222"})
    with pytest.raises(StaleCalibration, match="tagger serial 1111111 recorded, 2222222"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=other)
    policy = Policy.default().allow(ANALYSIS, "calibration_outside_validity", "skew checked against the old tagger")
    cset = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=other, policy=policy)
    assert OUTSIDE_VALIDITY in cset.labels and "skew checked against the old tagger" in policy.provenance()["overrides"][0]
    same = ctx(instrument={"model": "HydraHarp", "serial": "1111111"})
    assert OUTSIDE_VALIDITY not in CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=same).labels
    # Without the dataset's instrument nothing is compared.
    CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    # A dark artifact from another tagger is labelled only.
    d2 = tmp_path / "dark"
    write(make_cal(A, seed=2), d2, instrument={"model": "HydraHarp", "serial": "1111111"})
    labels = CalibrationSet.load(d2, analysis=ANALYSIS, dataset=other).labels
    assert OTHER_TAGGER in labels and OUTSIDE_VALIDITY not in labels


# --- Council round 9 ----------------------------------------------------------------------------

T2_DARK = {"tagger_mode": 2.0, "tagger_MeasDesc_GlobalResolution": 5e-12}
T3_RUN = {"tagger_mode": 3.0, "tagger_MeasDesc_GlobalResolution": 4.999937e-8, "tagger_MeasDesc_Resolution": 5e-12}


def test_a_t2_dark_artifact_covers_a_t3_dataset(tmp_path):
    # Council round 9 (Tw2): the T2 tick and the T3 sync period share one header key; a T2 dark
    # artifact (5 ps) was refused on a T3 run (50 ns), following DADL's own advice. Artifacts
    # without tagger_mode keep the old comparison.
    assert condition_mismatches(T2_DARK, T3_RUN) == []
    old = {"tagger_MeasDesc_GlobalResolution": 5e-12}
    assert condition_mismatches(old, T3_RUN) and condition_mismatches(T2_DARK, {**T3_RUN, "tagger_mode": 2.0})
    write(make_cal(A, seed=1, conditions=T2_DARK), tmp_path)
    CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(conditions=T3_RUN))


def test_a_changed_cw_t3_clock_is_labelled_not_refused(tmp_path):
    # Council round 9 (Tw3, the twin's hom_sync10): with source.type cw the sync is a clock only, so
    # a 20 to 10 MHz change labels dark and pair artifacts instead of refusing them.
    from dadl.propagate import T3_CLOCK_CHANGED

    write(make_cal(A, seed=1, conditions=T3_RUN, kind="heralding-efficiency"), tmp_path)
    run10 = {**T3_RUN, "tagger_MeasDesc_GlobalResolution": 9.999881e-8}
    with pytest.raises(StaleCalibration, match="tagger_MeasDesc_GlobalResolution = 4.999937e-08 recorded"):
        CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(conditions=run10))
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(conditions=run10, t3_cw=True))
    assert T3_CLOCK_CHANGED.format(a="4.999937e-08 s", b="9.999881e-08 s") in s.labels
    # A kind that may depend on the clock (deadtime-series) is still refused.
    series = tmp_path / "series"
    series.mkdir()
    write(make_cal(A, seed=2, conditions=T3_RUN, kind="deadtime-series"), series)
    with pytest.raises(StaleCalibration, match="tagger_MeasDesc_GlobalResolution"):
        CalibrationSet.load(series, analysis=ANALYSIS, dataset=ctx(conditions=run10, t3_cw=True),
                            kinds=["deadtime-series"])


def test_conditions_outside_the_compared_forms_are_labelled(tmp_path):
    # Council round 9 (Rd4): detector_temperature_K 233 against 298 passed silently, as did a lone
    # temperature_min_K.
    from dadl.propagate import NOT_COMPARED

    assert condition_mismatches({"temperature_min_K": 250.0}, {"temperature_K": 240.0})
    assert condition_mismatches({"temperature_min_K": 250.0}, {"temperature_K": 260.0}) == []
    assert condition_mismatches({"temperature_max_K": 250.0}, {"temperature_K": 260.0})
    assert condition_mismatches({"temperature_max_K": 250.0}, {}) == ["temperature_max_K = 250 recorded, not given for the dataset"]
    write(make_cal(A, seed=1, conditions={"detector_temperature_K": 233.0, "gaps": 2e4}), tmp_path)
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx(conditions={"detector_temperature_K": 298.0}))
    assert NOT_COMPARED.format(key="detector_temperature_K") in s.labels
    assert NOT_COMPARED.format(key="gaps") not in s.labels


def test_no_validity_labels_collapse_into_one(tmp_path):
    # Council round 9 (Tw8): CHSH repeated "has no validity range" once per artifact.
    from dadl.propagate import NO_VALIDITY_MANY, collapse_no_validity

    labels = (NO_VALIDITY.format(aid="aaa"), "other", NO_VALIDITY.format(aid="bbb"))
    assert collapse_no_validity(labels) == (NO_VALIDITY_MANY.format(aids="aaa, bbb"), "other")
    assert collapse_no_validity(labels[:2]) == labels[:2]
    write(make_cal(A, seed=1), tmp_path)
    write(make_cal(B, seed=2), tmp_path)
    s = CalibrationSet.load(tmp_path, analysis=ANALYSIS, dataset=ctx())
    ids = [m.artifact_id for m in s.members]
    assert s.labels_for(ids) == (NO_VALIDITY_MANY.format(aids=", ".join(a[:12] for a in ids)),)
