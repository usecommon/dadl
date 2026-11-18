"""Forward model and corrected HOM visibility with propagated calibrations (phase 3, M2)."""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from dadl.analysis.chsh import Window
from dadl.analysis.fringe import fringe_from_tags
from dadl.analysis.hom import hom_from_tags
from dadl.artifact import write_artifact
from dadl.calibration.dark import calibrate_dark
from dadl.data import ChannelId
from dadl.errors import GuardrailRefusal
from dadl.forward import DEAD_TIME_BOUND, FORWARD_MODEL, RunCounts, corrected_visibility, live_fraction
from dadl.guardrails import POLICY_MODIFIED, Policy
from dadl.propagate import CalibrationSet, DatasetContext
from dadl.sim import DetectorParams, PairSource, simulate_chsh_setting, simulate_dark, simulate_hom

UTC = dt.timezone.utc
CA, CB = ChannelId("sim", 1), ChannelId("sim", 2)
TRUE_V = 0.92
# High rates, so dead time and accidentals both bias the raw two-point visibility.
SRC = PairSource(pair_rate_per_s=2.0e6, visibility=TRUE_V, coherence_time_ps=2.0)
DET = DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=60_000, jitter_sigma_ps=10.0)
WINDOW = Window(width_ps=1_000, source="plan")
POLICY = Policy.default().allow("hom_visibility", "dead_time", "dead time calibrated from separate dark runs")


@pytest.fixture(scope="module")
def setup(tmp_path_factory):
    rng = np.random.default_rng(11)
    d = tmp_path_factory.mktemp("cal")
    for ch in (CA, CB):
        dark = simulate_dark(DET, 1.0, rng)
        write_artifact(calibrate_dark(dark, ch, backend="sim", seed=1, n_draws=400), d, created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    runs = {
        "dip": simulate_hom(SRC, DET, DET, 0.0, 0.5, rng),
        "reference": simulate_hom(SRC, DET, DET, 20.0, 0.5, rng),
    }
    ctx = DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC))
    cset = CalibrationSet.load(d, analysis="hom_visibility", dataset=ctx, policy=POLICY)
    return runs, cset


def test_corrected_hom_recovers_the_injected_visibility(setup):
    runs, cset = setup
    raw_only, _ = hom_from_tags(runs, WINDOW, subtract_accidentals=False)
    r, counts = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), seed=3)
    c = r.corrected
    # Raw is biased low by the accidentals: its exact interval excludes the truth, the corrected
    # one contains it. Dead time nearly cancels in this ratio (live fractions 0.914 and 0.913:
    # the reference run's higher coincidence rate offsets its higher singles), so the
    # calibration part of the error is small.
    assert r.raw.value == pytest.approx(raw_only.raw.value)
    assert r.raw.high < TRUE_V
    assert c.low <= TRUE_V <= c.high
    assert c.stderr_calibration > 0 and c.stderr_statistical > 0
    assert c.stderr == pytest.approx(np.hypot(c.stderr_statistical, c.stderr_calibration), rel=0.3)
    assert r.corrections == ("dead_time", "accidental_subtraction")
    assert r.artifact_ids == cset.artifact_ids and len(r.artifact_ids) == 2
    assert FORWARD_MODEL in r.labels and POLICY_MODIFIED in r.labels
    rec = r.as_record()
    assert rec["corrected_stderr_calibration"] == c.stderr_calibration
    assert "calibration" in rec["corrected_uncertainty"]


def test_dead_time_correction_needs_the_policy(setup):
    runs, cset = setup
    with pytest.raises(GuardrailRefusal, match="dead_time"):
        hom_from_tags(runs, WINDOW, calibrations=cset, channels=(CA, CB))


def test_stacked_corrections_are_refused(setup):
    runs, cset = setup
    with pytest.raises(GuardrailRefusal, match="twice"):
        hom_from_tags(runs, WINDOW, subtract_accidentals=True, policy=POLICY, calibrations=cset, channels=(CA, CB))


def test_without_calibrations_behaviour_is_unchanged(setup):
    runs, _ = setup
    r, _ = hom_from_tags(runs, WINDOW)
    assert r.corrected is not None and r.corrected.stderr_calibration is None
    assert r.corrections == () and r.artifact_ids == ()


def test_live_fraction_closed_form():
    # 1e5 cps on each side, 50 ns dead time: (1 - 0.005)^2, times 1 + r_c tau
    assert float(live_fraction(1e5, 1e5, 50e-9, 50e-9, 1e4)) == pytest.approx(0.995**2 * (1 + 1e4 * 50e-9))


def test_forward_model_closed_form_without_dead_time():
    # No calibrations: the forward model reduces to (C - A) per exposure; HOM V = 1 - s_dip / s_ref.
    high = RunCounts(10_000, 2_000, 20, 1.0, 10**5, 10**5)
    low = RunCounts(1_000, 2_000, 20, 1.0, 10**5, 10**5)
    c = corrected_visibility("hom_visibility", high, low, (CA, CB), calibrations=None, seed=1)
    assert c.estimate.value == pytest.approx(1 - (1_000 - 100) / (10_000 - 100))
    assert c.estimate.stderr_calibration == 0.0
    assert c.corrections == ("accidental_subtraction",)
    fringe = corrected_visibility(
        "bell_visibility", high, low, (CA, CB), calibrations=None, policy=Policy.default().allow("bell_visibility", "accidental_subtraction", "test"), seed=1
    )
    assert fringe.estimate.value == pytest.approx((9_900 - 900) / (9_900 + 900))


def test_fringe_correction_is_refused_under_the_default_policy():
    high = RunCounts(10_000, 2_000, 20, 1.0, 10**5, 10**5)
    low = RunCounts(1_000, 2_000, 20, 1.0, 10**5, 10**5)
    with pytest.raises(GuardrailRefusal):
        corrected_visibility("bell_visibility", high, low, (CA, CB), calibrations=None)


def _bounded_set(tmp_path, dark_rate_cps):
    """Dark artifacts whose dead time is stored only as an upper limit (too few short gaps)."""
    rng = np.random.default_rng(5)
    det = DetectorParams(efficiency=0.5, dark_rate_cps=dark_rate_cps, dead_time_ps=60_000)
    for ch in (CA, CB):
        cal = calibrate_dark(simulate_dark(det, 1.0, rng), ch, backend="sim", seed=1, n_draws=200)
        write_artifact(cal, tmp_path, accept_misfit=True, created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    policy = POLICY.allow("hom_visibility", "misfit_calibration", "test")
    cset = CalibrationSet.load(tmp_path, analysis="hom_visibility", dataset=DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC)), policy=policy)
    if not all("dead_time_ps" in m.calibration.unidentified for m in cset.members):
        pytest.skip("dead time was identified at this seed")
    return cset, policy, rng


def test_unidentified_dead_time_is_corrected_at_the_bound_and_reported_as_a_range(tmp_path):
    # 3,000 cps of dark counts bound the dead time below a few hundred ns; at 50k cps singles that
    # is a dead fraction of about 1%, an informative range.
    cset, policy, rng = _bounded_set(tmp_path, 3_000.0)
    src = PairSource(pair_rate_per_s=2.0e5, visibility=TRUE_V, coherence_time_ps=2.0)
    runs = {k: simulate_hom(src, DET, DET, d, 0.5, rng) for k, d in (("dip", 0.0), ("reference", 20.0))}
    r, counts = hom_from_tags(runs, WINDOW, policy=policy, calibrations=cset, channels=(CA, CB))
    assert DEAD_TIME_BOUND in r.labels
    assert counts["corrected_at_zero_dead_time"] != counts["corrected_at_bound_dead_time"]


def test_an_uninformative_dead_time_bound_is_refused(tmp_path):
    # 30 cps of dark counts bound the dead time only below about a millisecond: at 750k cps the
    # "range" would run from no correction to a negative live fraction.
    cset, policy, rng = _bounded_set(tmp_path, 30.0)
    runs = {k: simulate_hom(SRC, DET, DET, d, 0.1, rng) for k, d in (("dip", 0.0), ("reference", 20.0))}
    with pytest.raises(GuardrailRefusal, match="known only below"):
        hom_from_tags(runs, WINDOW, policy=policy, calibrations=cset, channels=(CA, CB))


# --- Fringe visibility from tags -------------------------------------------------------------

FRINGE_POLICY = (
    Policy.default()
    .allow("bell_visibility", "accidental_subtraction", "accidentals from each acquisition's off-peak bins")
    .allow("bell_visibility", "dead_time", "dead time calibrated from separate dark runs")
)
FRINGE_SRC = PairSource(pair_rate_per_s=2.0e6, visibility=0.95)
# A wide window, so the accidentals fill the fringe minimum (about 640 against 3,000 signal counts).
# Dead time alone raises the fringe (the detectors at the maximum are dead together), by less.
FRINGE_WINDOW = Window(width_ps=5_000, source="plan")


@pytest.fixture(scope="module")
def fringe_setup(tmp_path_factory):
    rng = np.random.default_rng(21)
    d = tmp_path_factory.mktemp("fringe-cal")
    for ch in (CA, CB):
        write_artifact(
            calibrate_dark(simulate_dark(DET, 1.0, rng), ch, backend="sim", seed=1, n_draws=400),
            d,
            created=dt.datetime(2026, 11, 1, tzinfo=UTC),
        )
    dets = {k: DET for k in ("A+", "A-", "B+", "B-")}

    def run(alpha, beta):
        tags = simulate_chsh_setting(FRINGE_SRC, alpha, beta, dets, 0.5, rng)
        return {"A": tags["A+"], "B": tags["B+"]}

    runs = {"max": run(0.0, 0.0), "min": run(0.0, 90.0)}
    cset = CalibrationSet.load(
        d, analysis="bell_visibility", dataset=DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC)), policy=FRINGE_POLICY
    )
    return runs, cset


def test_corrected_fringe_recovers_the_injected_visibility(fringe_setup):
    runs, cset = fringe_setup
    r, counts = fringe_from_tags(runs, FRINGE_WINDOW, policy=FRINGE_POLICY, calibrations=cset, channels=(CA, CB), seed=4)
    assert r.analysis == "bell_visibility"
    assert r.raw.high < FRINGE_SRC.visibility  # the accidentals fill the minimum
    assert r.corrected.low <= FRINGE_SRC.visibility <= r.corrected.high
    assert r.corrections == ("dead_time", "accidental_subtraction")
    assert r.artifact_ids == cset.artifact_ids
    assert FORWARD_MODEL in r.labels and POLICY_MODIFIED in r.labels
    assert counts["accidentals"] > 0


def test_fringe_raw_needs_no_calibrations_and_corrections_are_refused_by_default(fringe_setup):
    runs, cset = fringe_setup
    r, _ = fringe_from_tags(runs, FRINGE_WINDOW)
    assert r.corrected is None and r.raw.value < FRINGE_SRC.visibility
    with pytest.raises(GuardrailRefusal):
        fringe_from_tags(runs, FRINGE_WINDOW, calibrations=cset, channels=(CA, CB))
    with pytest.raises(GuardrailRefusal, match="twice"):
        fringe_from_tags(runs, FRINGE_WINDOW, subtract_accidentals=True, policy=FRINGE_POLICY, calibrations=cset, channels=(CA, CB))


def test_fringe_refuses_unequal_exposures(fringe_setup):
    runs, _ = fringe_setup
    cut = {"A": runs["min"]["A"][: runs["min"]["A"].size // 2], "B": runs["min"]["B"][: runs["min"]["B"].size // 2]}
    with pytest.raises(ValueError, match="same time"):
        fringe_from_tags({"max": runs["max"], "min": cut}, FRINGE_WINDOW)
