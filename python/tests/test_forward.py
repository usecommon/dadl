"""Forward model and corrected HOM visibility with propagated calibrations (phase 3, M2)."""

from __future__ import annotations

import datetime as dt
import json

import numpy as np
import pytest

from dadl.analysis.chsh import Window
from dadl.analysis.fringe import fringe_from_tags
from dadl.analysis.hom import hom_from_tags
from dadl.artifact import write_artifact
from dadl.calibration.dark import calibrate_dark
from dadl.data import ChannelId
from dadl.errors import GuardrailRefusal, Misspecified
from dadl.calibration.delay import calibrate_delay
from dadl.forward import (
    AFTERPULSING_UNMODELLED,
    CAPTURE_UNMODELLED,
    CORRECTED_ABOVE_ONE,
    DEAD_TIME_BOUND,
    FORWARD_MODEL,
    JITTER_BOUND,
    N_MC,
    RunCounts,
    corrected_visibility,
    grid_capture,
    live_fraction,
)
from dadl.forward_check import FORWARD_MISFIT_ACCEPTED
from dadl.guardrails import POLICY_MODIFIED, Policy
from dadl.propagate import CalibrationSet, DatasetContext
from dadl.sim import DetectorParams, PairSource, simulate_chsh_setting, simulate_dark, simulate_hom, simulate_pairs

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
    r, counts = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), seed=3, reference_delay_ps=20.0)
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


def test_the_record_carries_the_monte_carlo_seed_and_draws(setup):
    # Council round 5 (V3): the corrected interval depends on both, and neither was recorded.
    runs, cset = setup
    recs = []
    for seed in (0, 1):
        r, _ = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), seed=seed,
                             reference_delay_ps=20.0)
        recs.append(json.loads(r.as_record()["counts"])["monte_carlo"])
    assert recs[0] == {"seed": 0, "n_draws": N_MC} and recs[1] == {"seed": 1, "n_draws": N_MC}


def test_dead_time_correction_needs_the_policy(setup):
    runs, cset = setup
    with pytest.raises(GuardrailRefusal, match="dead_time"):
        hom_from_tags(runs, WINDOW, calibrations=cset, channels=(CA, CB), reference_delay_ps=20.0)


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
    r, counts = hom_from_tags(runs, WINDOW, policy=policy, calibrations=cset, channels=(CA, CB), reference_delay_ps=20.0)
    assert DEAD_TIME_BOUND in r.labels
    assert counts["corrected_at_zero_dead_time"] != counts["corrected_at_bound_dead_time"]


def test_an_uninformative_dead_time_bound_is_refused(tmp_path):
    # 30 cps of dark counts bound the dead time only below about a millisecond: at 750k cps the
    # "range" would run from no correction to a negative live fraction.
    cset, policy, rng = _bounded_set(tmp_path, 30.0)
    runs = {k: simulate_hom(SRC, DET, DET, d, 0.1, rng) for k, d in (("dip", 0.0), ("reference", 20.0))}
    with pytest.raises(GuardrailRefusal, match="known only below"):
        hom_from_tags(runs, WINDOW, policy=policy, calibrations=cset, channels=(CA, CB), reference_delay_ps=20.0)


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


# --- Window capture ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "offset,sigma,width", [(0, 35.4, 200), (15, 35.4, 200), (-40, 14, 1000), (300, 100, 200), (0, 5, 20), (500, 20, 10000)]
)
def test_grid_capture_closed_form_matches_the_integral(offset, sigma, width):
    # E[max(0, 1 - |D| / W)], D ~ N(offset, sigma), integrated on a fine grid around the peak. The
    # CLI's own quad integral read 2.4e-12 for (500, 20, 10000) (council round 4): it now calls this.
    from scipy import stats

    import dadl.cli

    d = np.linspace(offset - 12 * sigma, offset + 12 * sigma, 400_001)
    f = np.maximum(0.0, 1.0 - np.abs(d) / width) * stats.norm(offset, sigma).pdf(d)
    expected = float(np.sum(f) * (d[1] - d[0]))
    assert float(grid_capture(offset, sigma, width)) == pytest.approx(expected, rel=1e-6)
    assert not hasattr(dadl.cli, "_grid_capture")


def test_reference_delay_capture_is_modelled_with_a_delay_artifact(setup, tmp_path):
    runs, _ = setup
    rng = np.random.default_rng(31)
    a, b = simulate_pairs(SRC.pair_rate_per_s, DET, DET, 0.05, rng)
    write_artifact(
        calibrate_delay(a, b, (CA, CB), backend="sim", seed=1, n_draws=400), tmp_path, created=dt.datetime(2026, 11, 1, tzinfo=UTC)
    )
    ctx = DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC))
    policy = Policy.default().allow("hom_visibility", "window_capture", "delay and jitter from a pair run")
    cset = CalibrationSet.load(tmp_path, analysis="hom_visibility", dataset=ctx, policy=policy)
    kw = {"calibrations": cset, "channels": (CA, CB), "seed": 3}
    r, _ = hom_from_tags(runs, WINDOW, policy=policy, reference_delay_ps=20.0, **kw)
    assert r.corrections == ("accidental_subtraction", "window_capture")
    assert CAPTURE_UNMODELLED not in r.labels
    # The reference run's split pairs at +-20 ps lose more to the grid than the dip's: modelling it
    # raises V by C_dip/C_ref * (1 - capture ratio), here about 0.08 * 0.01.
    # No reference delay with calibrations is refused: the forward model needs the split (council
    # round 2: the reference was treated as unsplit and V read 10 sigma low, unlabelled).
    for missing in (None, 0.0):
        with pytest.raises(ValueError, match="reference_delay_ps"):
            hom_from_tags(runs, WINDOW, policy=policy, reference_delay_ps=missing, **kw)
    # Declaring about no delay for a run that has one is a misfit the predictive check sees (two
    # peaks at +-20 ps where the model has one); with the override the value can be compared.
    with pytest.raises(Misspecified, match="reference"):
        hom_from_tags(runs, WINDOW, policy=policy, reference_delay_ps=1.0, **kw)
    loose = policy.allow("hom_visibility", "misfit_forward_model", "comparing a misdeclared delay")
    r0, _ = hom_from_tags(runs, WINDOW, policy=loose, reference_delay_ps=1.0, **kw)
    assert FORWARD_MISFIT_ACCEPTED in r0.labels and FORWARD_MISFIT_ACCEPTED not in r.labels
    assert r.corrected.value - r0.corrected.value == pytest.approx(0.0008, abs=0.0003)
    # Not in the default HOM rule, like dead time (D10).
    with pytest.raises(GuardrailRefusal, match="window_capture"):
        hom_from_tags(runs, WINDOW, reference_delay_ps=20.0, **kw)


def test_reference_delay_without_a_delay_artifact_is_labelled(setup):
    runs, cset = setup
    r, _ = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), reference_delay_ps=20.0)
    assert CAPTURE_UNMODELLED in r.labels and "window_capture" not in r.corrections


# --- Council round 2 ----------------------------------------------------------------------------

def _member_set(analysis, cals):
    from types import MappingProxyType

    from dadl.propagate import Member

    members = tuple(Member(f"{i:064x}", c, MappingProxyType({})) for i, c in enumerate(cals, start=7))
    return CalibrationSet(analysis, members, ())


def _dark(ch, tau_ps=50_000.0, sd_ps=2_000.0, bound=None, seed=0, afterpulse=False):
    """A dark artifact holding the dead time only (or, with `afterpulse`, a fitted afterpulse
    probability as well); with `bound`, the dead time is a stand-in upper limit."""
    from dadl.calibration import dark
    from dadl.calibration.base import Calibration
    from dadl.stats import Interval

    rng = np.random.default_rng(seed)
    if bound is not None:
        col, est, flags = rng.uniform(0, bound, 2000), Interval(bound, 0.0, bound), ("dead_time_ps",)
    else:
        col, est, flags = rng.normal(tau_ps, sd_ps, 2000), Interval(tau_ps, tau_ps - 2 * sd_ps, tau_ps + 2 * sd_ps), ()
    estimates, units, samples, params = {"dead_time_ps": est}, {"dead_time_ps": "ps"}, [col], ["dead_time_ps"]
    if afterpulse:
        estimates["afterpulse_prob"] = Interval(0.01, 0.008, 0.012)
        units["afterpulse_prob"] = "1"
        samples.append(rng.normal(0.01, 0.001, 2000))
        params.append("afterpulse_prob")
    return Calibration(
        kind=dark.KIND, method=dark.METHOD, method_version=dark.METHOD_VERSION, channels=(ch,), estimates=estimates,
        units=units, samples=np.column_stack(samples), sample_params=tuple(params), converged=True, ppc_pvalue=0.5,
        inputs=(("x" * 64, "raw-tags", str(ch)),), backend="sim", unidentified=flags, stand_in=flags,
    )


def test_a_bounded_dead_time_gives_value_and_interval_at_the_bound():
    # Council round 2 (R1): the point was at zero and the interval at the bound, so with unequal
    # singles the value fell outside its own interval (0.9009 against [0.9045, 0.9072]).
    from dadl.forward import bound_range

    policy = Policy.default().allow("hom_visibility", "dead_time", "test")
    high = RunCounts(200_000, 4_000, 20, 1.0, 1_500_000, 1_500_000)
    low = RunCounts(20_000, 4_000, 20, 1.0, 1_200_000, 1_200_000)
    cset = _member_set("hom_visibility", [_dark(CA, bound=100_000), _dark(CB, bound=100_000)])
    c = corrected_visibility("hom_visibility", high, low, (CA, CB), calibrations=cset, policy=policy, seed=1)
    assert c.estimate.low <= c.estimate.value <= c.estimate.high
    ends = bound_range("hom_visibility", high, low, (CA, CB), calibrations=cset, policy=policy)
    assert c.estimate.value == pytest.approx(ends["at_bound"])
    assert DEAD_TIME_BOUND in c.labels


def test_a_threshold_crossed_only_at_the_bound_is_labelled():
    from dadl.forward import bound_range
    from dadl.guardrails import THRESHOLD_AFTER_CORRECTION, threshold_labels
    from dadl.result import Estimate

    policy = Policy.default().allow("hom_visibility", "dead_time", "test")
    high = RunCounts(400_000, 400, 20, 1.0, 1_500_000, 1_500_000)
    low = RunCounts(201_000, 400, 20, 1.0, 1_200_000, 1_200_000)
    cset = _member_set("hom_visibility", [_dark(CA, bound=100_000), _dark(CB, bound=100_000)])
    c = corrected_visibility("hom_visibility", high, low, (CA, CB), calibrations=cset, policy=policy, seed=1).estimate
    ends = bound_range("hom_visibility", high, low, (CA, CB), calibrations=cset, policy=policy)
    raw = Estimate(1 - 201_000 / 400_000, 0.0015, 0.4945, 0.5005)
    assert ends["at_zero"] < 0.5 < ends["at_bound"]
    assert THRESHOLD_AFTER_CORRECTION in threshold_labels("hom_visibility", raw, c, tuple(ends.values()))
    # The ends alone carry the crossing when the value does not.
    below = Estimate(0.499, 0.001, 0.497, 0.501)
    assert THRESHOLD_AFTER_CORRECTION in threshold_labels("hom_visibility", raw, below, (0.4975, 0.501))
    assert THRESHOLD_AFTER_CORRECTION not in threshold_labels("hom_visibility", raw, below)


def test_a_stand_in_jitter_is_used_at_its_bound(tmp_path):
    # Council round 2 (R3): an 80 ps tagger leaves the jitter a stand-in upper limit; every
    # corrected analysis raised KeyError on the missing jitter draws.
    rng = np.random.default_rng(11)
    t = np.sort(rng.uniform(0, 1e12, 20_000))
    a = (np.floor((t + 17 + rng.normal(0, 3, t.size)) / 80) * 80).astype(np.int64)
    b = (np.floor((t + 1030 + rng.normal(0, 3, t.size)) / 80) * 80).astype(np.int64)
    cal = calibrate_delay(np.sort(a), np.sort(b), (CA, CB), bin_width_ps=80, tick_ps=80, backend="sim")
    assert "jitter_sigma_ps" in cal.stand_in and cal.status == "ok"
    write_artifact(cal, tmp_path, created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    policy = Policy.default().allow("hom_visibility", "window_capture", "test")
    cset = CalibrationSet.load(
        tmp_path, analysis="hom_visibility", dataset=DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC)), policy=policy
    )
    high = RunCounts(200_000, 4_000, 20, 1.0, 1_000_000, 1_000_000, (20.0, -20.0))
    low = RunCounts(20_000, 4_000, 20, 1.0, 1_000_000, 1_000_000)
    c = corrected_visibility(
        "hom_visibility", high, low, (CA, CB), calibrations=cset, policy=policy, window_ps=800, declared_delay_ps=1013
    )
    assert np.isfinite(c.estimate.value) and "window_capture" in c.corrections
    # A window narrow against the bound: the capture at the bound differs from that at zero jitter.
    c = corrected_visibility(
        "hom_visibility", high, low, (CA, CB), calibrations=cset, policy=policy, window_ps=80, declared_delay_ps=1013
    )
    assert JITTER_BOUND in c.labels


def test_a_corrected_visibility_above_one_is_labelled_and_kept():
    high = RunCounts(500, 1_600, 200, 1.0, 10_000, 10_000)
    low = RunCounts(5, 1_600, 200, 1.0, 10_000, 10_000)
    c = corrected_visibility("hom_visibility", high, low, (CA, CB), calibrations=None)
    assert c.estimate.value > 1 and CORRECTED_ABOVE_ONE in c.labels


def test_only_the_artifacts_the_model_read_are_cited():
    # Council round 2 (R6): a deadtime-series artifact was cited though the model never read it,
    # and afterpulsing reported by a dark calibration went unlabelled (R7).
    import dataclasses

    from dadl.calibration import deadtime

    policy = Policy.default().allow("hom_visibility", "dead_time", "test")
    series = dataclasses.replace(
        _dark(CA, seed=3), kind=deadtime.KIND, method=deadtime.METHOD, method_version=deadtime.METHOD_VERSION
    )
    cset = _member_set("hom_visibility", [_dark(CA, seed=1), _dark(CB, seed=2), series])
    high = RunCounts(200_000, 4_000, 20, 1.0, 1_500_000, 1_500_000)
    low = RunCounts(20_000, 4_000, 20, 1.0, 1_200_000, 1_200_000)
    c = corrected_visibility("hom_visibility", high, low, (CA, CB), calibrations=cset, policy=policy, seed=1)
    assert len(c.artifact_ids) == 2 and cset.members[2].artifact_id not in c.artifact_ids
    assert AFTERPULSING_UNMODELLED not in c.labels
    fitted = _member_set("hom_visibility", [_dark(CA, seed=1, afterpulse=True), _dark(CB, seed=2)])
    c = corrected_visibility("hom_visibility", high, low, (CA, CB), calibrations=fitted, policy=policy, seed=1)
    assert AFTERPULSING_UNMODELLED in c.labels


def test_a_rate_drift_between_fringe_acquisitions_is_labelled():
    # Council round 2 (H1): +6% pair rate on the max run moved V across 1/sqrt(2) unlabelled.
    from dadl.analysis.fringe import SINGLES_DRIFT, TWO_POINT

    det = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=350.0)
    dets = dict.fromkeys(("A+", "A-", "B+", "B-"), det)
    for drift, labelled in ((0.0, False), (0.06, True)):
        rng = np.random.default_rng(100)
        rmax = simulate_chsh_setting(PairSource(pair_rate_per_s=2e5 * (1 + drift), visibility=0.70), 0.0, 0.0, dets, 2.0, rng)
        rmin = simulate_chsh_setting(PairSource(pair_rate_per_s=2e5, visibility=0.70), 0.0, 90.0, dets, 2.0, rng)
        runs = {"max": {"A": rmax["A+"], "B": rmax["B+"]}, "min": {"A": rmin["A+"], "B": rmin["B+"]}}
        r, counts = fringe_from_tags(runs, Window(width_ps=2000, source="plan"))
        assert TWO_POINT in r.labels
        assert (SINGLES_DRIFT in r.labels) == labelled
        assert counts["singles_ratio_a"] == pytest.approx(1 + drift, abs=0.02)


# --- Council round 4 ----------------------------------------------------------------------------

def test_a_partial_dark_set_gives_the_named_error():
    # Council round 4 (S2): a dark artifact for A only gave a bare KeyError (a CLI traceback).
    from dadl.forward import bound_range

    policy = Policy.default().allow("hom_visibility", "dead_time", "test")
    high = RunCounts(10_000, 500, 10, 1.0, 100_000, 100_000)
    low = RunCounts(3_000, 500, 10, 1.0, 100_000, 100_000)
    partial = _member_set("hom_visibility", [_dark(CA)])
    with pytest.raises(ValueError, match="none for B"):
        corrected_visibility("hom_visibility", high, low, (CA, CB), calibrations=partial, policy=policy)
    with pytest.raises(ValueError, match="none for B"):
        bound_range("hom_visibility", high, low, (CA, CB), calibrations=partial)


def test_hom_checks_the_dark_dead_times_against_the_analysed_gaps(setup):
    # Council round 4 (S4): no test failed with the HOM gap check removed. A 1 us dead time on
    # 60 ns detectors is refused downward, a 22 ns one upward; a matched one passes.
    from dadl.errors import StaleCalibration
    from dadl.propagate import OUTSIDE_VALIDITY

    runs, _ = setup
    kw = {"channels": (CA, CB), "reference_delay_ps": 20.0}
    long = _member_set("hom_visibility", [_dark(CA, 1e6), _dark(CB, 1e6, seed=1)])
    with pytest.raises(StaleCalibration, match="analysed gaps are shorter"):
        hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=long, **kw)
    loose = POLICY.allow("hom_visibility", "calibration_outside_validity", "checking")
    r, counts = hom_from_tags(runs, WINDOW, policy=loose, calibrations=long, **kw)
    assert OUTSIDE_VALIDITY in r.labels and 55_000 < counts["min_gap_ps"]["A"] < 61_000
    short = _member_set("hom_visibility", [_dark(CA, 22_000, 500), _dark(CB, 22_000, 500, seed=1)])
    with pytest.raises(StaleCalibration, match="longer than its upper end"):
        hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=short, **kw)
    matched = _member_set("hom_visibility", [_dark(CA, 60_000, 300), _dark(CB, 60_000, 300, seed=1)])
    r, counts = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=matched, **kw)
    assert OUTSIDE_VALIDITY not in r.labels and counts["gap_check_floor_ps"]["A"] == pytest.approx(58_200)
    # The forward model subtracts the accidentals (council round 4, S8: 0.0 was recorded).
    assert "accidental_subtraction" in r.corrections and counts["subtracted"] == 1.0


def test_hom_refuses_another_detector_s_artifacts(setup):
    # Council round 4 (S3): hom_from_tags took another detector's artifacts on the API path.
    from dadl.errors import StaleCalibration
    from dadl.propagate import OUTSIDE_VALIDITY

    runs, _ = setup
    old = [ChannelId(c.device, c.input, role=n, detector=f"X{n}") for n, c in (("A", CA), ("B", CB))]
    new = tuple(ChannelId(c.device, c.input, role=n, detector=f"NEW{n}") for n, c in (("A", CA), ("B", CB)))
    cset = _member_set("hom_visibility", [_dark(old[0], 60_000, 300), _dark(old[1], 60_000, 300, seed=1)])
    with pytest.raises(StaleCalibration, match="calibrated on detector XA, analysing NEWA"):
        hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=new, reference_delay_ps=20.0)
    loose = POLICY.allow("hom_visibility", "calibration_outside_validity", "same detectors, renamed")
    r, _ = hom_from_tags(runs, WINDOW, policy=loose, calibrations=cset, channels=new, reference_delay_ps=20.0)
    assert OUTSIDE_VALIDITY in r.labels
