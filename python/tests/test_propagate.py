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
    NOT_CHECKED,
    OUTSIDE_VALIDITY,
    SERIES_UNUSED,
    CalibrationSet,
    DatasetContext,
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
    assert s.labels == ()


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
    # diagnostics and keys the dataset does not give are not compared
    assert condition_mismatches({"gaps": 2e4, "temperature_K": 253.0}, {}) == []


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
    assert "window_capture" in r.corrections and NOT_CHECKED in r.labels
