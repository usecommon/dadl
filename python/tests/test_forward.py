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
    # Council round 7 (Q4): without a delay artifact the corrected V carries the unmodelled capture
    # of the split reference (labelled; mean z -0.75, 69 of 80 seeds covered), so this tested a 95%
    # event of a biased estimator at one seed. With the delay artifact the capture is modelled; the
    # truth is asserted within 3 standard errors (20 of 20 seeds, rng 11 + 1000 s).
    import dataclasses

    from dadl.propagate import Member

    runs, cset = setup
    a, b = simulate_pairs(SRC.pair_rate_per_s, DET, DET, 0.02, np.random.default_rng(12))
    delay = calibrate_delay(a, b, (CA, CB), backend="sim", seed=1, n_draws=400)
    full = dataclasses.replace(cset, members=cset.members + (Member("d" * 64, delay, {}),))
    policy = POLICY.allow("hom_visibility", "window_capture", "delay and jitter from a separate pair run")
    raw_only, _ = hom_from_tags(runs, WINDOW, subtract_accidentals=False)
    r, counts = hom_from_tags(runs, WINDOW, policy=policy, calibrations=full, channels=(CA, CB), seed=3,
                              reference_delay_ps=20.0)
    c = r.corrected
    # Raw is biased low by the accidentals. Dead time nearly cancels in this ratio (live fractions
    # 0.914 and 0.913: the reference run's higher coincidence rate offsets its higher singles), so
    # the calibration part of the error is small.
    assert r.raw.value == pytest.approx(raw_only.raw.value)
    assert r.raw.value < c.value
    assert abs(c.value - TRUE_V) < 3 * c.stderr
    assert CAPTURE_UNMODELLED not in r.all_labels
    assert c.stderr_calibration > 0 and c.stderr_statistical > 0
    assert c.stderr == pytest.approx(np.hypot(c.stderr_statistical, c.stderr_calibration), rel=0.3)
    assert r.corrections == ("dead_time", "accidental_subtraction", "window_capture")
    assert len(r.artifact_ids) == 3
    assert FORWARD_MODEL in r.all_labels and POLICY_MODIFIED in r.all_labels
    rec = r.as_record()
    assert rec["corrected_stderr_calibration"] == c.stderr_calibration
    assert "calibration" in rec["corrected_uncertainty"]
    # Without the delay artifact the split reference's capture is unmodelled, and labelled.
    r2, _ = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), seed=3,
                          reference_delay_ps=20.0)
    assert CAPTURE_UNMODELLED in r2.all_labels and r2.corrections == ("dead_time", "accidental_subtraction")
    assert np.isfinite(r2.corrected.value)


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
    # Unequal dead times (council round 7, Q3): a pair is lost to the shorter one's recovery only.
    assert float(live_fraction(1e5, 2e5, 50e-9, 20e-9, 1e4)) == pytest.approx(0.995 * 0.996 * (1 + 1e4 * 20e-9))


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


def _bounded_set(tmp_path, bound_ps):
    """Dark artifacts whose dead time is stored only as an upper limit (as a dark run with too few
    short gaps gives it), built directly: council round 7 (Q8) found the simulated version
    identified the dead time, and skipped, at 6 of 20 seeds."""
    for i, ch in enumerate((CA, CB)):
        write_artifact(_dark(ch, bound=bound_ps, seed=i), tmp_path, created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    cset = CalibrationSet.load(tmp_path, analysis="hom_visibility", dataset=DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC)), policy=POLICY)
    assert all("dead_time_ps" in m.calibration.unidentified for m in cset.members)
    return cset, POLICY, np.random.default_rng(5)


def test_unidentified_dead_time_is_corrected_at_the_bound_and_reported_as_a_range(tmp_path):
    # 3,000 cps of dark counts bound the dead time below a few hundred ns; at 50k cps singles that
    # is a dead fraction of about 1%, an informative range.
    cset, policy, rng = _bounded_set(tmp_path, 300_000.0)
    src = PairSource(pair_rate_per_s=2.0e5, visibility=TRUE_V, coherence_time_ps=2.0)
    runs = {k: simulate_hom(src, DET, DET, d, 0.5, rng) for k, d in (("dip", 0.0), ("reference", 20.0))}
    r, counts = hom_from_tags(runs, WINDOW, policy=policy, calibrations=cset, channels=(CA, CB), reference_delay_ps=20.0)
    assert DEAD_TIME_BOUND in r.all_labels
    assert counts["corrected_at_zero_dead_time"] != counts["corrected_at_bound_dead_time"]


def test_an_uninformative_dead_time_bound_is_refused(tmp_path):
    # 30 cps of dark counts bound the dead time only below about a millisecond: at 750k cps the
    # "range" would run from no correction to a negative live fraction.
    cset, policy, rng = _bounded_set(tmp_path, 1e9)
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
    assert FORWARD_MODEL in r.all_labels and POLICY_MODIFIED in r.all_labels
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
    assert CAPTURE_UNMODELLED not in r.all_labels
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
    assert FORWARD_MISFIT_ACCEPTED in r0.all_labels and FORWARD_MISFIT_ACCEPTED not in r.all_labels
    assert r.corrected.value - r0.corrected.value == pytest.approx(0.0008, abs=0.0003)
    # Not in the default HOM rule, like dead time (D10).
    with pytest.raises(GuardrailRefusal, match="window_capture"):
        hom_from_tags(runs, WINDOW, reference_delay_ps=20.0, **kw)


def test_reference_delay_without_a_delay_artifact_is_labelled(setup):
    runs, cset = setup
    r, _ = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), reference_delay_ps=20.0)
    assert CAPTURE_UNMODELLED in r.all_labels and "window_capture" not in r.corrections


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
    loose = POLICY.allow("hom_visibility", "dead_time_contradicted", "checking")
    # Under the override the gap check labels it, and the dead fraction (0.85) is then refused
    # (council round 6, R1 b).
    with pytest.raises(GuardrailRefusal, match="pair term and the unmodelled paralyzable stage are not trusted"):
        hom_from_tags(runs, WINDOW, policy=loose, calibrations=long, **kw)
    short = _member_set("hom_visibility", [_dark(CA, 22_000, 500), _dark(CB, 22_000, 500, seed=1)])
    with pytest.raises(StaleCalibration, match="longer than its upper end"):
        hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=short, **kw)
    matched = _member_set("hom_visibility", [_dark(CA, 60_000, 300), _dark(CB, 60_000, 300, seed=1)])
    r, counts = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=matched, **kw)
    assert OUTSIDE_VALIDITY not in r.all_labels and counts["gap_check_floor_ps"]["A"] == pytest.approx(58_200)
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
    assert OUTSIDE_VALIDITY in r.all_labels


# Review council, round 6.

#: hom-snspd dark runs (20 s at 100 cps, rng 10,000 + r, seed r) whose 30 ns dead time is reported
#: identified at microseconds (r6-council ident_hom20: 7 of 200, none covering).
SNSPD_FALSE_IDENT = (58, 64, 73, 99, 113, 157, 185)
#: chsh-spad dark runs (2 s at 500 cps) identified with three widths below zero, covering 22 ns.
CHSH_LOOSE_COVERING = (6, 8, 9, 10, 22, 38, 61, 74, 79, 116, 139, 144, 151, 160, 164, 166, 180, 192, 196, 197)


@pytest.fixture(scope="module")
def hom_snspd():
    from dadl.cli import example_path
    from dadl.config import load_setup
    from dadl.simsetup import simulate_setup

    setup = load_setup(example_path("hom-snspd"))
    return setup, simulate_setup(setup)


def _example_dark(setup, name: str, duration_s: float, r: int):
    return calibrate_dark(simulate_dark(setup.detectors[name], duration_s, np.random.default_rng(10_000 + r)),
                          setup.channel_id(name), backend="sim", seed=r)


def _snspd_hom(hom_snspd, r: int, runs=None):
    setup, tags = hom_snspd
    cal = _example_dark(setup, "A", 20.0, r)
    assert "dead_time_ps" not in cal.unidentified and cal.estimates["dead_time_ps"].low > 1e5
    other = dataclasses_replace_channel(cal, setup.channel_id("B"))
    cset = _member_set("hom_visibility", [cal, other])
    policy = Policy.default().allow("hom_visibility", "dead_time", "test")
    return hom_from_tags(runs or tags, setup.window, policy=policy, calibrations=cset,
                         channels=(setup.channel_id("A"), setup.channel_id("B")),
                         reference_delay_ps=setup.reference_delay_ps)


def _snspd_refused(hom_snspd, r: int) -> None:
    # Refused by the dead fraction (R1 b) at the referee's 85 kcps; the lower-end gap test only
    # labels (council round 6 addendum).
    with pytest.raises(GuardrailRefusal, match="pair term and the unmodelled paralyzable stage are not trusted"):
        _snspd_hom(hom_snspd, r)


def dataclasses_replace_channel(cal, channel):
    import dataclasses

    return dataclasses.replace(cal, channels=(channel,))


def test_a_falsely_identified_snspd_dead_time_is_refused(hom_snspd):
    # Council round 6 (R1): an SNSPD's 30 ns read 20 us, "identified", and gave V = 0.55 +- 6.2 with
    # exit 0. Its dead fraction at the analysed rate (1.7) is refused now.
    _snspd_refused(hom_snspd, SNSPD_FALSE_IDENT[0])


@pytest.mark.slow
@pytest.mark.parametrize("r", SNSPD_FALSE_IDENT[1:])
def test_every_falsely_identified_snspd_dead_time_is_refused(hom_snspd, r):
    _snspd_refused(hom_snspd, r)


def test_a_false_snspd_dead_time_at_a_low_rate_is_labelled(hom_snspd):
    # Council round 6 addendum: at 5 kcps the 26 us dead time is a dead fraction of 0.13, under the
    # refusal; the analysed gaps below its lower end label the corrected value instead.
    from dadl.forward import DEAD_TIME_BELOW_DARK

    _, tags = hom_snspd
    rng = np.random.default_rng(113)
    thin = {}
    for key, run in tags.items():
        span = max(int(a.max()) for a in run.values()) - min(int(a.min()) for a in run.values())
        thin[key] = {ch: np.sort(rng.choice(a, int(5_000 * span / 1e12), replace=False)) for ch, a in run.items()}
    r, counts = _snspd_hom(hom_snspd, 113, thin)
    head = DEAD_TIME_BELOW_DARK.split("{name}")[0]
    assert any(lab.startswith(head) for lab in r.corrected_labels)
    assert not any(lab.startswith(head) for lab in r.labels)
    assert counts["gap_check_floor_kind"]["A"] == "lower end (label only)"


def test_a_channel_never_carries_both_dead_time_labels(hom_snspd):
    # Council round 10 (Rg2): at about 4 kcps the false SNSPD identification of seed 58 carried
    # DEAD_TIME_BELOW_DARK ("likely too large") and DEAD_TIME_ABOVE_GAPS ("uses the dark value") on the
    # same channel: 4 dead-time labels for 2 channels (r10-reg p/r6/snspd.py, keep 0.045).
    from dadl.forward import DEAD_TIME_ABOVE_GAPS, DEAD_TIME_ABOVE_GAPS_BOUND, DEAD_TIME_BELOW_DARK

    setup, tags = hom_snspd
    rng = np.random.default_rng(1)
    thin = {k: {c: np.sort(v[c][rng.random(v[c].size) < 0.045]) for c in v} for k, v in tags.items()}
    r, _ = _snspd_hom(hom_snspd, 58, thin)
    heads = [t.split("{name}")[0] for t in (DEAD_TIME_BELOW_DARK, DEAD_TIME_ABOVE_GAPS, DEAD_TIME_ABOVE_GAPS_BOUND)]
    dead = [lab for lab in r.corrected_labels if any(lab.startswith(h) for h in heads)]
    assert dead
    for name in ("A", "B"):
        tag = f"{name} ({setup.channel_id(name)})"
        assert sum(tag in lab for lab in dead) <= 1, dead


def _chsh_spad_gap_check(r: int):
    from dadl.cli import example_path
    from dadl.config import load_setup
    from dadl.forward import check_dead_time_gaps
    from dadl.simsetup import simulate_setup

    setup = load_setup(example_path("chsh-spad"))
    cal = _example_dark(setup, "A+", 2.0, r)
    t = cal.estimates["dead_time_ps"]
    assert t.contains(setup.detectors["A+"].dead_time_ps) and t.value - 3 * (t.value - t.low) <= 0
    runs = simulate_setup(setup)
    labels, counts = check_dead_time_gaps(_member_set("chsh", [cal]), {"A+": setup.channel_id("A+")},
                                          list(runs.values()), analysis="chsh")
    # A loose point above some gaps is labelled DEAD_TIME_ABOVE_GAPS (council round 8), nothing else.
    assert tuple(lab for lab in labels if not lab.startswith("the dark dead time of")) == ()
    assert counts["gap_check_floor_kind"]["A+"] == ("lower end (label only)" if t.low > 0 else None)


def test_a_loose_covering_dark_dead_time_still_passes():
    # The lower-end test neither refuses nor labels loose identifications that cover the truth
    # (council: 20 of 22 chsh-spad 2 s ones).
    _chsh_spad_gap_check(CHSH_LOOSE_COVERING[0])


@pytest.mark.slow
def test_every_loose_covering_dark_dead_time_passes():
    for r in CHSH_LOOSE_COVERING[1:]:
        _chsh_spad_gap_check(r)


def test_an_identified_dead_time_above_the_dead_fraction_limit_is_refused():
    # Council round 6 (R1 b): a 22 ns dead time at 1e7 cps is a dead fraction of 0.22; it was
    # labelled only.
    run = RunCounts(coincidences=50_000, offpeak_sum=1_000, offpeak_shifts=10, duration_s=0.1,
                    singles_a=1_000_000, singles_b=1_000_000)
    cset = _member_set("hom_visibility", [_dark(CA, 22_000, 300), _dark(CB, 22_000, 300, seed=1)])
    with pytest.raises(GuardrailRefusal, match=r"rate x dead time = 0.22 \(limit 0.2\).*pair term and the unmodelled paralyzable stage"):
        corrected_visibility("hom_visibility", run, run, (CA, CB), calibrations=cset, policy=POLICY)
    slow = dataclasses_replace_run(run, 500_000)
    corrected_visibility("hom_visibility", slow, slow, (CA, CB), calibrations=cset, policy=POLICY)
    # The faster of the two runs sets the dead fraction (council round 7, Q3).
    with pytest.raises(GuardrailRefusal, match="rate x dead time = 0.22"):
        corrected_visibility("hom_visibility", run, slow, (CA, CB), calibrations=cset, policy=POLICY)
    with pytest.raises(GuardrailRefusal, match="rate x dead time = 0.22"):
        corrected_visibility("hom_visibility", slow, run, (CA, CB), calibrations=cset, policy=POLICY)


def dataclasses_replace_run(run, singles):
    import dataclasses

    return dataclasses.replace(run, singles_a=singles, singles_b=singles)


def test_calibrate_dark_warns_when_the_dead_time_is_far_above_the_declared(tmp_path, capsys):
    # Council round 6 (R1 c): the referee's channel A artifact (20 us for 30 ns) printed no warning.
    from dadl.cli import main

    dark = tmp_path / "dark.npz"
    assert main(["simulate", "hom-snspd", "--dark", "--out", str(dark)]) == 0
    assert main(["calibrate", "dark", "hom-snspd", "--tags", str(dark), "--channel", "A", "--out",
                 str(tmp_path / "arts")]) == 0
    text = capsys.readouterr().out
    assert "far above the setup's declared 30000 ps" in text and "usually a chance miss" in text


def test_corrected_labels_are_kept_apart_from_the_raw_ones(setup):
    # Council round 6 (R2, R3): labels of the corrected value were printed on the raw line, in one
    # flat list; forward-model V was not marked secondary.
    from dadl.forward import SECONDARY_V

    runs, cset = setup
    r, counts = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), seed=3,
                              reference_delay_ps=20.0)
    assert FORWARD_MODEL in r.corrected_labels and FORWARD_MODEL not in r.labels
    # Council round 10 (Ux2): "policy modified" is on the corrected value, not the raw one.
    assert SECONDARY_V in r.corrected_labels and POLICY_MODIFIED in r.corrected_labels
    assert POLICY_MODIFIED not in r.labels
    rec = r.as_record()
    assert FORWARD_MODEL in json.loads(rec["corrected_labels"]) and FORWARD_MODEL not in json.loads(rec["labels"])
    # R4: what the calibration part leaves out, from the artifacts.
    assert "tagger timebase and non-linearity" in counts["calibration_excludes"]


def test_threshold_labels_are_evaluated_per_value():
    # Council round 6 (R2): the referee's HOM run printed "the exact 95% interval contains the
    # threshold" beside raw V = 0.9186 [0.9162, 0.9209]; it was the corrected interval.
    from dadl.guardrails import NEAR_THRESHOLD_EXACT, corrected_threshold_labels, threshold_labels
    from dadl.result import Estimate

    raw = Estimate(0.9186, 0.0012, 0.9162, 0.9209)
    corrected = Estimate(0.5457, 6.1808, -1.1011, 2.7874, stderr_statistical=0.0066, stderr_calibration=6.1631)
    assert NEAR_THRESHOLD_EXACT not in threshold_labels("hom_visibility", raw)
    assert NEAR_THRESHOLD_EXACT in corrected_threshold_labels("hom_visibility", raw, corrected)


def test_fringe_forward_value_is_secondary():
    from dadl.forward import SECONDARY_V

    run = RunCounts(coincidences=10_000, offpeak_sum=100, offpeak_shifts=10, duration_s=1.0, singles_a=10**5,
                    singles_b=10**5)
    low = RunCounts(coincidences=500, offpeak_sum=100, offpeak_shifts=10, duration_s=1.0, singles_a=10**5,
                    singles_b=10**5)
    policy = Policy.default().allow("bell_visibility", "accidental_subtraction", "test")
    c = corrected_visibility("bell_visibility", run, low, (CA, CB), calibrations=None, policy=policy)
    assert SECONDARY_V in c.labels


def test_hom_and_fringe_without_dark_artifacts_say_dead_time_is_not_modelled():
    # Council round 10 (Ux12): CHSH said so, HOM and fringe did not.
    from dadl.forward import DEAD_TIME_UNMODELLED, chsh_model_label

    run = RunCounts(coincidences=10_000, offpeak_sum=100, offpeak_shifts=10, duration_s=1.0, singles_a=10**5,
                    singles_b=10**5)
    low = RunCounts(coincidences=500, offpeak_sum=100, offpeak_shifts=10, duration_s=1.0, singles_a=10**5,
                    singles_b=10**5)
    assert DEAD_TIME_UNMODELLED.split(" not")[0] in chsh_model_label(("accidental_subtraction",))
    for analysis in ("bell_visibility", "hom_visibility"):
        policy = Policy.default().allow(analysis, "accidental_subtraction", "test")
        c = corrected_visibility(analysis, run, low, (CA, CB), calibrations=None, policy=policy)
        assert DEAD_TIME_UNMODELLED in c.labels
    dark = _member_set("hom_visibility", [_dark(CA, 20_000), _dark(CB, 20_000, seed=1)])
    policy = Policy.default().allow("hom_visibility", "dead_time", "test")
    c = corrected_visibility("hom_visibility", run, low, (CA, CB), calibrations=dark, policy=policy)
    assert DEAD_TIME_UNMODELLED not in c.labels


def test_calibration_part_prints_two_significant_figures():
    # Council round 6 (R4): "calibration 0.0000" read as exactly known.
    from dadl.result import Estimate

    text = str(Estimate(0.9198, 0.0012, 0.9175, 0.9222, stderr_statistical=0.0012, stderr_calibration=1.3e-6))
    assert "stat 0.0012, calibration part 1.3e-06 (statistical part of the artifacts only)" in text
    # Council round 10 (Ux11): the corrected "±" is one sd, the interval 95%.
    assert text.startswith("0.9198 ± 0.0012 (1 sd) (95%: [")


def test_a_zero_calibration_part_prints_as_zero():
    # Council round 10 (Rg5 a): an exactly-zero part printed as "calibration 4.4e-16".
    from dadl.result import Estimate

    for cal in (0.0, 4.4e-16, 2.2e-16):
        text = str(Estimate(2.6947, 0.0105, 2.67, 2.71, stderr_statistical=0.0105, stderr_calibration=cal))
        assert "calibration part 0 (" in text
    assert "calibration part 6.3e-06" in str(Estimate(2.69, 0.01, 2.67, 2.71, 0.01, 6.3e-06))
    # A raw value keeps its exact-interval wording and no "(1 sd)".
    assert str(Estimate(0.92, 0.01, 0.90, 0.94)) == "0.9200 ± 0.0100 (exact 95%: [0.9000, 0.9400])"


def test_hom_refusals_say_hom_and_how_to_drop_the_dark_artifacts():
    # Council round 6 (R7): HOM refusals reused the Bell wording, with a Python-only remedy.
    from dadl.forward import WITHOUT_DARK

    with pytest.raises(GuardrailRefusal) as exc:
        Policy.default().check_correction("hom_visibility", "window_capture")
    assert "the headline V uses raw counts" in str(exc.value) and "Bell" not in str(exc.value)
    assert "--allow window_capture=REASON" in str(exc.value)
    assert "move the dark artifacts out of the calibrations folder" in WITHOUT_DARK
    run = RunCounts(coincidences=50_000, offpeak_sum=1_000, offpeak_shifts=10, duration_s=0.1,
                    singles_a=100_000, singles_b=100_000)
    loose = _member_set("hom_visibility", [_dark(CA, bound=1e7), _dark(CB, bound=1e7, seed=1)])
    with pytest.raises(GuardrailRefusal, match="move the dark artifacts out"):
        corrected_visibility("hom_visibility", run, run, (CA, CB), calibrations=loose, policy=POLICY)


def test_cli_hom_prints_the_corrected_labels_and_caveats(tmp_path, capsys):
    # Council round 6 (R2, R3, R4, R8) through the CLI: the corrected line is marked secondary, its
    # labels and exclusions are printed under it, and the declaration caveat is printed.
    from dadl.cli import main

    tags, her, arts = tmp_path / "tags.npz", tmp_path / "her.npz", tmp_path / "arts"
    assert main(["simulate", "hom-snspd", "--out", str(tags)]) == 0
    assert main(["simulate", "hom-snspd", "--heralding", "--duration-s", "0.2", "--out", str(her)]) == 0
    assert main(["calibrate", "pairs", "hom-snspd", "--tags", str(her), "--no-analysers", "--herald", "A",
                 "--signal", "B", "--out", str(arts)]) == 0
    capsys.readouterr()
    assert main(["analyze", "hom-snspd", "--tags", str(tags), "--calibrations", str(arts),
                 "--allow", "window_capture=delay and jitter from a pair run"]) == 0
    out = capsys.readouterr().out
    corrected = next(line for line in out.splitlines() if line.strip().startswith("Corrected (forward model"))
    assert "(secondary; not an indistinguishability claim)" in corrected
    raw = next(line for line in out.splitlines() if line.strip().startswith("Labels:"))
    assert "forward model" not in raw and "no validity range" not in raw
    assert "Labels of the corrected value:" in out and "Calibration part excludes: tagger timebase" in out
    assert "choosing them from the data biases V upward" in out


# --- Council round 7 (Q2): a delay-jitter artifact in the other channel order -------------------

#: B 400 ps late on 25 ps detectors: the declared window delay is +400 ps.
SIGN_DET = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=30_000, jitter_sigma_ps=25.0)
SIGN_LATE = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=30_000, jitter_sigma_ps=25.0, delay_ps=400)


def sign_sets(tmp_path, analysis, policy, pairs):
    """For each order, (A, B) and (B, A), a CalibrationSet of one delay-jitter artifact per channel
    pair, fitted from the same pair runs. `pairs` maps (id_a, id_b) to (tags_a, tags_b)."""
    created = dt.datetime(2026, 11, 1, tzinfo=UTC)
    ctx = DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC))
    out = {}
    for order in ("ab", "ba"):
        d = tmp_path / order
        d.mkdir()
        for (ia, ib), (a, b) in pairs.items():
            cal = (calibrate_delay(a, b, (ia, ib), backend="sim", seed=1, n_draws=400) if order == "ab"
                   else calibrate_delay(b, a, (ib, ia), backend="sim", seed=1, n_draws=400))
            write_artifact(cal, d, accept_misfit=True, created=created)
        out[order] = CalibrationSet.load(d, analysis=analysis, dataset=ctx, policy=policy)
    return out["ab"], out["ba"]


def assert_same_correction(ab, ba):
    """The corrected values of the two orders agree (only their Monte Carlo draws differ), and the
    predictive check ran on the analysed peak in both."""
    (r_ab, c_ab), (r_ba, c_ba) = ab, ba
    assert np.isfinite(r_ab.corrected.value)
    assert r_ba.corrected.value == pytest.approx(r_ab.corrected.value, abs=0.1 * r_ab.corrected.stderr)
    assert r_ba.corrected.stderr == pytest.approx(r_ab.corrected.stderr, rel=0.05)
    assert "forward_check_p" in c_ab and "forward_check_p" in c_ba
    assert c_ba["forward_check_p"] == pytest.approx(c_ab["forward_check_p"], rel=1e-6)


def test_fringe_takes_a_delay_artifact_in_either_channel_order(tmp_path):
    # Fringe V is a ratio of two runs' signals at one capture, so a capture computed at the wrong
    # delay cancels unless it vanishes: 10 ps detectors make a capture 700 ps off exactly 0.
    import dataclasses

    rng = np.random.default_rng(31)
    det, late = (dataclasses.replace(d, jitter_sigma_ps=10.0) for d in (SIGN_DET, SIGN_LATE))
    dets = {"A+": det, "A-": det, "B+": late, "B-": det}
    pair = simulate_pairs(2e5, det, late, 0.05, rng)
    policy = (Policy.default().allow("bell_visibility", "accidental_subtraction", "test")
              .allow("bell_visibility", "window_capture", "test"))
    sets = sign_sets(tmp_path, "bell_visibility", policy, {(CA, CB): pair})

    def run(alpha, beta):
        tags = simulate_chsh_setting(PairSource(pair_rate_per_s=2e5, visibility=0.95), alpha, beta, dets, 0.2, rng)
        return {"A": tags["A+"], "B": tags["B+"]}

    runs = {"max": run(0.0, 0.0), "min": run(0.0, 90.0)}
    window = Window(width_ps=200, source="plan", delays_ps={"B": 400})
    ab, ba = (fringe_from_tags(runs, window, policy=policy, calibrations=c, channels=(CA, CB), seed=2) for c in sets)
    assert ab[0].corrections == ("accidental_subtraction", "window_capture")
    assert_same_correction(ab, ba)
    c = ab[0].corrected
    assert abs(c.value - 0.95) < 4 * c.stderr


# --- Council round 7 (Q3): labels and refusals observed through the analysis path -------------


def _with_dataset(analysis, cals):
    """A set of these artifacts checked against a dataset context, as `CalibrationSet.load` makes it."""
    from types import MappingProxyType

    from dadl.propagate import Member

    ctx = DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC))
    members = tuple(Member(f"{i:064x}", c, MappingProxyType({})) for i, c in enumerate(cals, start=7))
    return CalibrationSet(analysis, members, (), True, ctx)


def test_fringe_checks_the_dark_dead_times_against_the_analysed_gaps(fringe_setup):
    # The fringe twin of the HOM test: refused when the gaps rule the dark dead time out, labelled
    # when only its lower end is below them; the check's p and the dataset are recorded.
    from dadl.errors import StaleCalibration
    from dadl.forward import DEAD_TIME_BELOW_DARK

    runs, _ = fringe_setup
    policy = FRINGE_POLICY.allow("bell_visibility", "window_capture", "test")
    a, b = simulate_pairs(2e5, DET, DET, 0.05, np.random.default_rng(22))
    delay = calibrate_delay(a, b, (CA, CB), backend="sim", seed=1, n_draws=400)
    kw = {"policy": policy, "channels": (CA, CB), "seed": 1}
    long = _with_dataset("bell_visibility", [_dark(CA, 1e6), _dark(CB, 1e6, seed=1), delay])
    with pytest.raises(StaleCalibration, match="analysed gaps are shorter"):
        fringe_from_tags(runs, FRINGE_WINDOW, calibrations=long, **kw)
    # 120 ns with a 95% lower end of 80 ns on 60 ns detectors: three widths below the point is at 0
    # (not tested), and the gaps fall below the lower end: a label, not a refusal.
    loose = _with_dataset("bell_visibility", [_dark(CA, 120_000, 20_000), _dark(CB, 120_000, 20_000, seed=1), delay])
    r, counts = fringe_from_tags(runs, FRINGE_WINDOW, calibrations=loose, **kw)
    head = DEAD_TIME_BELOW_DARK.split("{name}")[0]
    assert any(lab.startswith(head) for lab in r.corrected_labels)
    assert counts["gap_check_floor_kind"]["A"] == "lower end (label only)"
    assert 0.0 <= counts["forward_check_p"] <= 1.0
    assert counts["dataset"]["acquired"] == "2026-11-02T00:00:00+00:00"


def test_a_wide_dead_time_range_is_labelled_on_hom():
    # D3 range wider than the corrected value's error is labelled (the CHSH copy was the only one tested).
    from dadl.analysis._twopoint import forward_from_counts
    from dadl.forward import DEAD_TIME_RANGE_WIDE

    high = RunCounts(200_000, 4_000, 20, 1.0, 600_000, 600_000, (20.0, -20.0))
    low = RunCounts(20_000, 4_000, 20, 1.0, 600_000, 600_000)
    seen = {}
    for bound in (20_000.0, 60_000.0):
        cset = _member_set("hom_visibility", [_dark(CA, bound=bound), _dark(CB, bound=bound, seed=1)])
        fwd, extra = forward_from_counts("hom_visibility", high, low, (CA, CB), calibrations=cset, policy=POLICY,
                                         seed=1, window_ps=1_000, declared_delay_ps=0.0)
        width = abs(extra["corrected_at_bound_dead_time"] - extra["corrected_at_zero_dead_time"])
        seen[bound] = width / fwd.estimate.stderr
        assert (DEAD_TIME_RANGE_WIDE in fwd.labels) == (width > fwd.estimate.stderr)
    assert seen[20_000.0] < 1.0 < seen[60_000.0] < 10.0


def test_fringe_subtracts_each_acquisition_s_own_accidentals(fringe_setup, monkeypatch):
    # Council round 7 (Q3): the maximum's accidentals taken from the minimum's side windows survived.
    from dadl.analysis import fringe

    runs, cset = fringe_setup
    seen = []
    real = fringe.forward_from_counts

    def spy(analysis, high, low, *args, **kwargs):
        seen.append((high, low))
        return real(analysis, high, low, *args, **kwargs)

    monkeypatch.setattr(fringe, "forward_from_counts", spy)
    _, counts = fringe_from_tags(runs, FRINGE_WINDOW, policy=FRINGE_POLICY, calibrations=cset, channels=(CA, CB), seed=4)
    (high, low), = seen
    assert counts["accidentals_max"] != counts["accidentals"]
    assert high.accidentals == counts["accidentals_max"] and low.accidentals == counts["accidentals"]


def test_hom_subtracts_the_reference_s_own_accidentals(setup):
    # Council round 7 (Q3): the reference's accidentals taken from the dip's side windows survived
    # (equal singles in simulation); under drift they differ.
    from dadl.calibration.accidentals import shifted_coincidences, side_offsets

    runs, cset = setup
    _, counts = hom_from_tags(runs, WINDOW, policy=POLICY, calibrations=cset, channels=(CA, CB), seed=3,
                              reference_delay_ps=20.0)
    ref = shifted_coincidences(np.sort(runs["reference"]["A"]), np.sort(runs["reference"]["B"]), WINDOW.width_ps, 0, 0,
                               side_offsets(WINDOW.width_ps))
    assert counts["accidentals_reference"] == ref.sum() / ref.size
    assert counts["accidentals_reference"] != counts["accidentals"]


# --- Council round 8 ----------------------------------------------------------------------------


def _dark_counted(ch, tau_ps, sd_ps, gaps, rate_cps, seed=0):
    """`_dark` with the dark run's gap count and dark rate recorded, as calibrate_dark writes them."""
    import dataclasses

    from dadl.stats import Interval

    cal = _dark(ch, tau_ps, sd_ps, seed=seed)
    estimates = {**cal.estimates, "dark_rate_cps": Interval(rate_cps, 0.99 * rate_cps, 1.01 * rate_cps)}
    return dataclasses.replace(cal, estimates=estimates, units={**cal.units, "dark_rate_cps": "s-1"},
                               conditions={"gaps": float(gaps)})


def test_the_downward_floor_follows_the_dark_run_s_gap_law():
    # Council round 8 (St1, Fix 3): the floor is point - ln(K / 1e-4) / (N f), not three widths.
    from dadl.errors import StaleCalibration
    from dadl.forward import check_dead_time_gaps, gap_edge_excess, gap_floor_ps

    # r8-stat darktail extremes: a 263 ns point from 2,000 gaps at 20 kcps cannot be refused.
    assert gap_floor_ps(263_000.0, 2_000, 20_000.0, 4) <= 0
    assert gap_edge_excess(2) == pytest.approx(10.41, abs=0.01) and gap_edge_excess(4) == pytest.approx(11.18, abs=0.01)
    # A 1 us artifact from a 1 s, 50 kcps dark run on a channel with 30 ns gaps is still refused.
    rng = np.random.default_rng(3)
    runs = [{"A": np.cumsum(30_000 + rng.exponential(2e5, 20_000)).astype(np.int64)}]
    cset = _member_set("hom_visibility", [_dark_counted(CA, 1e6, 2_000, 50_000, 50_000.0)])
    with pytest.raises(StaleCalibration, match="analysed gaps are shorter"):
        check_dead_time_gaps(cset, {"A": CA}, runs, analysis="hom_visibility")
    loose = Policy.default().allow("hom_visibility", "dead_time_contradicted", "checking")
    _, counts = check_dead_time_gaps(cset, {"A": CA}, runs, analysis="hom_visibility", policy=loose)
    assert counts["gap_check_floor_kind"]["A"] == "exponential law"
    assert counts["gap_check_floor_ps"]["A"] == pytest.approx(gap_floor_ps(1e6, 50_000, 49_500.0, 1))


def test_a_dark_dead_time_above_the_analysed_gaps_is_labelled():
    # Council round 8 (Ph4, Fix 7): a dark point above many analysed gaps is labelled; the
    # correction keeps it. Council round 10 (Rg2): the label states how often it fires.
    from dadl.forward import DEAD_TIME_ABOVE_GAPS, check_dead_time_gaps

    rng = np.random.default_rng(4)
    runs = [{"A": np.cumsum(50_000 + rng.exponential(3e5, 20_000)).astype(np.int64)}]
    cset = _member_set("hom_visibility", [_dark_counted(CA, 52_000, 2_000, 2_000, 20_000.0)])
    labels, _ = check_dead_time_gaps(cset, {"A": CA}, runs, analysis="hom_visibility")
    head = DEAD_TIME_ABOVE_GAPS.split("{name}")[0]
    assert any(lab.startswith(head) and "52 ns" in lab and "fires on most analyses" in lab for lab in labels)
    below = _member_set("hom_visibility", [_dark_counted(CA, 49_000, 2_000, 2_000, 20_000.0)])
    labels, _ = check_dead_time_gaps(below, {"A": CA}, runs, analysis="hom_visibility")
    assert not any(lab.startswith(head) for lab in labels)


# --- Council round 9 ----------------------------------------------------------------------------


def _echo_channel(n_gaps, echoes, seed=0, edge_ps=1_000_000, mean_ps=1_000_000.0, echo_ps=50_000):
    """A channel with a 1 us edge and `echoes` clicks followed by a cable echo 50 ns later."""
    rng = np.random.default_rng(seed)
    tags = np.cumsum(edge_ps + rng.exponential(mean_ps, n_gaps)).astype(np.int64)
    pick = rng.choice(tags.size, echoes, replace=False)
    return np.sort(np.concatenate([tags, tags[pick] + echo_ps]))


@pytest.mark.parametrize("echoes", [2, 500])
def test_a_few_echoes_do_not_turn_off_the_upward_gap_test(echoes):
    # Council round 9 (Rd1): one echo set the shortest gap to the 50 ns dark bound, and more than
    # 100 made the spacing 0 (NaN): a 50 ns bound on a 1 us detector passed and S read 2.03 for
    # 1.98. The test reads gap j = max(10, 0.1% of the gaps) instead.
    from dadl.errors import StaleCalibration
    from dadl.forward import check_dead_time_gaps

    runs = [{"A": _echo_channel(1_000_000, echoes)}]
    cset = _member_set("hom_visibility", [_dark(CA, bound=50_000.0)])
    with pytest.raises(StaleCalibration, match="dead time below 50 ns of A"):
        check_dead_time_gaps(cset, {"A": CA}, runs, analysis="hom_visibility")


def test_the_upward_gap_test_keeps_the_round_4_power_and_fails_closed():
    # Council round 9 (Rd1): a 22 ns artifact on a 50 ns detector at 1 Mcps, 50,000 gaps, is
    # still refused; an all-equal gap run (zero spacing) is refused, not passed as NaN.
    from dadl.errors import StaleCalibration
    from dadl.forward import check_dead_time_gaps

    rng = np.random.default_rng(4)
    runs = [{"A": np.cumsum(50_000 + rng.exponential(1e6, 50_000)).astype(np.int64)}]
    cset = _member_set("hom_visibility", [_dark(CA, 22_000, 500)])
    with pytest.raises(StaleCalibration, match="longer than its upper end"):
        check_dead_time_gaps(cset, {"A": CA}, runs, analysis="hom_visibility")
    flat = [{"A": np.arange(1, 5_001, dtype=np.int64) * 1_000_000}]
    with pytest.raises(StaleCalibration, match="longer than its upper end by inf"):
        check_dead_time_gaps(_member_set("hom_visibility", [_dark(CA, bound=50_000.0)]), {"A": CA}, flat,
                             analysis="hom_visibility")


def _clean_excess_refusals(n_channels, n_gaps, seed, channels=4):
    from dadl.forward import gap_edge_excess, upward_gap_excess

    rng = np.random.default_rng(seed)
    refused = 0
    for _ in range(n_channels):
        gaps = np.sort(50_000 + rng.exponential(2e5, n_gaps))
        excess, j, m = upward_gap_excess(gaps, 50_000.0)
        refused += not excess <= gap_edge_excess(channels, j, m)
    return refused


def test_the_upward_threshold_is_the_exact_law_of_an_exponential_edge():
    # Council round 9 (Rd1): (j + 1) F.isf(alpha / K, 2 (j + 1), 2 M), 100((K / alpha)^(1/100) - 1)
    # at j = 0; the red team's g[k] with the g[0] threshold refused every clean channel.
    from scipy import stats

    from dadl.forward import gap_edge_excess

    assert gap_edge_excess(4, 0, 100) == pytest.approx(100 * ((4 / 1e-4) ** 0.01 - 1), rel=1e-9)
    assert gap_edge_excess(4, 10, 100) == pytest.approx(11 * stats.f.isf(2.5e-5, 22, 200), rel=1e-12)
    assert gap_edge_excess(2, 1000, 1001) > gap_edge_excess(4, 10, 100)
    assert _clean_excess_refusals(200, 2_000, seed=1) == 0


@pytest.mark.slow
def test_clean_exponential_edges_at_the_upper_end_are_not_refused():
    # Council round 9 (gapup.txt): 0 of 4,000 clean channels with the upper end at tau.
    assert _clean_excess_refusals(4_000, 2_000, seed=2) == 0


def test_the_gap_tests_have_their_own_override():
    # Council round 9 (Rd2): the validity override no longer accepts a contradicted dark dead time
    # (an echo-set 50 ns bound on a 1 us detector); dead_time_contradicted does, with the problems
    # in the label and the counts.
    from dadl.errors import StaleCalibration
    from dadl.forward import DEAD_TIME_CONTRADICTED, check_dead_time_gaps

    rng = np.random.default_rng(6)
    runs = [{"A": np.cumsum(1_000_000 + rng.exponential(1e6, 20_000)).astype(np.int64)}]
    cset = _member_set("hom_visibility", [_dark(CA, bound=50_000.0)])
    validity = Policy.default().allow("hom_visibility", "calibration_outside_validity", "taken before the artifact")
    with pytest.raises(StaleCalibration, match="--allow dead_time_contradicted=REASON"):
        check_dead_time_gaps(cset, {"A": CA}, runs, analysis="hom_visibility", policy=validity)
    loose = Policy.default().allow("hom_visibility", "dead_time_contradicted", "echo in the dark run")
    labels, counts = check_dead_time_gaps(cset, {"A": CA}, runs, analysis="hom_visibility", policy=loose)
    head = DEAD_TIME_CONTRADICTED.split("{problems}")[0]
    assert [lab for lab in labels if lab.startswith(head)] == [head + counts["dead_time_gap_problems"][0]]
    assert "dead time below 50 ns of A" in counts["dead_time_gap_problems"][0]


def test_a_dark_dead_time_far_below_the_declared_one_is_flagged():
    # Council round 9 (Rd1b): one rule for calibrate dark's warning and analyze's label.
    from dadl.forward import dead_time_far_below

    assert dead_time_far_below(_dark(CA, bound=50_000.0), 1_000_000) == 50_000.0
    assert dead_time_far_below(_dark(CA, 400_000, 1_000), 1_000_000) == 400_000.0
    assert dead_time_far_below(_dark(CA, 600_000, 1_000), 1_000_000) is None
    assert dead_time_far_below(_dark(CA, bound=50_000.0), 0) is None
