"""Corrected CHSH from the forward model (phase 3, M3)."""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

from dadl.analysis.chsh import Window, chsh
from dadl.artifact import write_artifact
from dadl.calibration.dark import calibrate_dark
from dadl.calibration.delay import calibrate_delay
from dadl.closed_form import chsh_ideal
from dadl.data import ChannelId
from dadl.errors import GuardrailRefusal
from dadl.forward import FORWARD_MODEL, chsh_model_label
from dadl.guardrails import POLICY_MODIFIED, Policy
from dadl.propagate import CalibrationSet, DatasetContext
from dadl.sim import DetectorParams, PairSource, simulate_chsh, simulate_dark, simulate_pairs

UTC = dt.timezone.utc
SETTINGS = {"a": 0.0, "a'": 45.0, "b": 22.5, "b'": 67.5}
IDS = {name: ChannelId("sim", i) for i, name in enumerate(("A+", "A-", "B+", "B-"))}
# Dead-time selection: detectors that fire together are dead together, so raw |E| reads high.
SRC = PairSource(pair_rate_per_s=3.0e6, visibility=0.70)
DET = DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=50_000, jitter_sigma_ps=50.0)
WINDOW = Window(width_ps=1_000, source="plan")
TRUTH = chsh_ideal(SRC.visibility, SETTINGS)


def _policy(*corrections: str) -> Policy:
    policy = Policy.default()
    for c in corrections:
        policy = policy.allow("chsh", c, "calibrated from separate runs")
    return policy


POLICY = _policy("dead_time", "accidental_subtraction", "window_capture")
CTX = DatasetContext(acquired=dt.datetime(2026, 11, 2, tzinfo=UTC))


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    rng = np.random.default_rng(1)
    d = tmp_path_factory.mktemp("chsh-cal")
    created = dt.datetime(2026, 11, 1, tzinfo=UTC)
    for cid in IDS.values():
        write_artifact(calibrate_dark(simulate_dark(DET, 1.0, rng), cid, backend="sim", seed=1, n_draws=400), d, accept_misfit=True, created=created)
    for x in "+-":
        for y in "+-":
            a, b = simulate_pairs(SRC.pair_rate_per_s, DET, DET, 0.02, rng)
            cal = calibrate_delay(a, b, (IDS["A" + x], IDS["B" + y]), backend="sim", seed=1, n_draws=400)
            write_artifact(cal, d, accept_misfit=True, created=created)
    runs = simulate_chsh(SRC, dict.fromkeys(IDS, DET), 0.2, rng, SETTINGS)
    return runs, d


@pytest.mark.slow
def test_corrected_s_removes_the_dead_time_selection(tmp_path):
    # Council round 7 (Q4): at 0.2 s the selection is about 1.6 raw standard errors, and "raw
    # excludes the truth" held at 7 of 20 seeds. At 2 s it is about 4.7 (20 seeds: 2.6 to 7.0).
    # The shift raw - corrected is checked against the dead-time model at the true dead time
    # (50 ns on every channel): the calibrated correction agrees within 3 se (20 seeds: |d| < 0.9).
    from test_forward import _dark

    rng = np.random.default_rng(1001)
    _calibrate(tmp_path, dict.fromkeys(IDS, DET), rng, dark_s=1.0, pairs_s=0.02, rate=SRC.pair_rate_per_s)
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=POLICY, dataset=CTX)
    runs = simulate_chsh(SRC, dict.fromkeys(IDS, DET), 2.0, rng, SETTINGS)
    r = chsh(runs, WINDOW, POLICY, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    raw_only = chsh(runs, WINDOW, settings_deg=SETTINGS)
    assert r.raw == raw_only.raw  # the headline is untouched
    assert r.raw.value > TRUTH  # dead-time selection raises raw S
    assert r.corrected.value < r.raw.value
    true = _with_darks(cset, [_dark(cid, 50_000, 1.0, seed=i) for i, cid in enumerate(IDS.values())])
    model = chsh(runs, WINDOW, POLICY, settings_deg=SETTINGS, calibrations=true, channels=IDS, seed=1)
    shift, expected = r.raw.value - r.corrected.value, r.raw.value - model.corrected.value
    assert abs(shift - expected) < 3 * r.corrected.stderr
    assert r.corrections == ("dead_time", "accidental_subtraction", "window_capture")
    assert len(r.artifact_ids) == 8
    for label in (FORWARD_MODEL, chsh_model_label(r.corrections), POLICY_MODIFIED):
        assert label in r.all_labels
    assert r.corrected.stderr_calibration > 0


def test_every_chsh_correction_needs_the_policy(data):
    runs, d = data
    # The default rule for CHSH allows none (Bell-type: raw only).
    for missing in ("dead_time", "accidental_subtraction", "window_capture"):
        policy = _policy(*(c for c in ("dead_time", "accidental_subtraction", "window_capture") if c != missing))
        cset = CalibrationSet.load(d, analysis="chsh", policy=policy, dataset=CTX)
        with pytest.raises(GuardrailRefusal, match=missing):
            chsh(runs, WINDOW, policy, calibrations=cset, channels=IDS)


def test_channels_and_complete_pairs_are_required(data, tmp_path):
    runs, d = data
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    with pytest.raises(ValueError, match="channels"):
        chsh(runs, WINDOW, POLICY, calibrations=cset)
    # Delay artifacts for only some channel pairs: refused rather than half-modelled.
    one = next(p for p in d.glob("*.json") if "delay-jitter" in p.read_text())
    for f in d.glob(one.name.split(".")[0] + "*"):
        (tmp_path / f.name).write_bytes(f.read_bytes())
    partial = CalibrationSet.load(tmp_path, analysis="chsh", policy=POLICY, dataset=CTX)
    with pytest.raises(ValueError, match="every channel pair"):
        chsh(runs, WINDOW, POLICY, calibrations=partial, channels=IDS)


# --- Council round 2 ----------------------------------------------------------------------------

def _calibrate(directory, dets, rng, *, dark_s=None, pairs_s=0.05, rate=2e5):
    created = dt.datetime(2026, 11, 1, tzinfo=UTC)
    if dark_s:
        for name, cid in IDS.items():
            cal = calibrate_dark(simulate_dark(dets[name], dark_s, rng), cid, backend="sim", seed=1, n_draws=400)
            write_artifact(cal, directory, accept_misfit=True, created=created)
    for x in "+-":
        for y in "+-":
            a, b = simulate_pairs(rate, dets["A" + x], dets["B" + y], pairs_s, rng)
            cal = calibrate_delay(a, b, (IDS["A" + x], IDS["B" + y]), backend="sim", seed=1, n_draws=400)
            write_artifact(cal, directory, accept_misfit=True, created=created)


def test_a_bounded_dead_time_gives_the_range_in_counts(tmp_path):
    # Council round 2 (B2): an SNSPD's dark run identifies no dead time; corrected S is at the bound
    # with the value at zero in the counts, and a range wider than the error is labelled.
    from dadl.calibration.delay import KIND as DELAY_KIND
    from dadl.forward import DEAD_TIME_BOUND, DEAD_TIME_RANGE_WIDE

    det = DetectorParams(efficiency=0.8, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=25.0)
    dets = dict.fromkeys(IDS, det)
    rng = np.random.default_rng(900)
    _calibrate(tmp_path, dets, rng, dark_s=20.0, pairs_s=0.2, rate=1e5)
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=POLICY, dataset=CTX)
    runs = simulate_chsh(PairSource(pair_rate_per_s=5e4, visibility=0.95), dets, 1.0, rng, SETTINGS)
    window = Window(width_ps=500, source="plan")
    r = chsh(runs, window, POLICY, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    c = r.corrected
    assert c.low <= c.value <= c.high
    assert c.value == pytest.approx(r.counts["corrected_at_bound_dead_time"])
    assert DEAD_TIME_BOUND in r.all_labels and DEAD_TIME_RANGE_WIDE in r.all_labels
    delays_only = CalibrationSet(cset.analysis, tuple(m for m in cset.members if m.kind == DELAY_KIND), cset.labels)
    r0 = chsh(runs, window, POLICY, settings_deg=SETTINGS, calibrations=delays_only, channels=IDS, seed=1)
    assert r.counts["corrected_at_zero_dead_time"] == pytest.approx(r0.corrected.value, abs=1e-12)
    assert "corrected_at_zero_dead_time" not in r0.counts


def test_a_stand_in_jitter_corrects_chsh(tmp_path):
    # Council round 2 (R3): a coarse tagger's jitter is an upper limit; CHSH raised KeyError.
    import dataclasses

    from dadl.forward import RunCounts, corrected_chsh
    from dadl.analysis.chsh import SIGNS

    rng = np.random.default_rng(11)
    t = np.sort(rng.uniform(0, 1e12, 20_000))
    a = np.sort((np.floor((t + 17 + rng.normal(0, 3, t.size)) / 80) * 80).astype(np.int64))
    b = np.sort((np.floor((t + 1030 + rng.normal(0, 3, t.size)) / 80) * 80).astype(np.int64))
    coarse = calibrate_delay(a, b, (IDS["A+"], IDS["B+"]), bin_width_ps=80, tick_ps=80, backend="sim")
    assert "jitter_sigma_ps" in coarse.stand_in
    for x in "+-":
        for y in "+-":
            cal = dataclasses.replace(coarse, channels=(IDS["A" + x], IDS["B" + y]))
            write_artifact(cal, tmp_path, created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    policy = _policy("accidental_subtraction", "window_capture")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=policy, dataset=CTX)
    run = {(x, y): RunCounts(5_000 if x == y else 1_000, 400, 20, 1.0, 100_000, 100_000) for x in "+-" for y in "+-"}
    c = corrected_chsh(
        {k: run for k in SIGNS}, IDS, SIGNS, calibrations=cset, policy=policy, window_ps=1_000,
        delays_ps={"B+": 1013.0, "B-": 1013.0},
    )
    assert np.isfinite(c.estimate.value) and "window_capture" in c.corrections
    # Council round 5 (V3): the Monte Carlo seed and draw count go into the result's counts.
    assert c.counts["monte_carlo"] == {"seed": 0, "n_draws": 4000}
    # Council round 7 (Q3): here the capture at the jitter bound differs from that at zero jitter,
    # and the corrected S says so.
    from dadl.forward import JITTER_BOUND

    assert JITTER_BOUND in c.labels


MISDECLARED_DET = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=150.0)
MISDECLARED_SRC = PairSource(pair_rate_per_s=2e5, visibility=0.70)


def _late(shift_ps):
    """A+ and B+ `shift_ps` later than A- and B-."""
    import dataclasses

    return {n: dataclasses.replace(MISDECLARED_DET, delay_ps=shift_ps if n in ("A+", "B+") else 0) for n in IDS}


def test_delays_that_miss_the_pair_peaks_are_labelled():
    # Council round 2 (B1): A+ and B+ 400 ps late, declared 0: raw S 2.20 for a 1.98 source, unlabelled.
    from dadl.analysis.chsh import CAPTURE_IMBALANCE

    truth = chsh_ideal(MISDECLARED_SRC.visibility, SETTINGS)
    rng = np.random.default_rng(7)
    runs = simulate_chsh(MISDECLARED_SRC, _late(400), 1.0, rng, SETTINGS)
    prefix = CAPTURE_IMBALANCE.split("{")[0]
    r = chsh(runs, Window(width_ps=1000, source="plan"), settings_deg=SETTINGS)
    assert r.raw.value - truth > 20 * r.raw.stderr
    assert any(lab.startswith(prefix) for lab in r.all_labels)
    assert abs(r.counts["capture_reweighted_s"] - truth) < 3 * r.raw.stderr
    assert r.counts["capture_offsets_ps"]["A+/B-"] == pytest.approx(-400, abs=15)
    right = chsh(runs, Window(width_ps=1000, source="plan", delays_ps={"A+": 400, "B+": 400}), settings_deg=SETTINGS)
    assert not any(lab.startswith(prefix) for lab in right.all_labels)


@pytest.mark.slow
def test_matched_delays_are_not_labelled():
    # chsh-spad physics (350 ps jitter, 2 ns window, afterpulsing): |raw - reweighted| stays well
    # under half a standard error.
    import math

    from dadl.analysis.chsh import CAPTURE_IMBALANCE, CAPTURE_LOW

    prefix = CAPTURE_IMBALANCE.split("{")[0]
    low = CAPTURE_LOW.split("{")[0]
    det = DetectorParams(
        efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=350.0, afterpulse_prob=0.005,
        afterpulse_tau_ps=50_000,
    )
    for seed in range(10):
        rng = np.random.default_rng(100 + seed)
        runs = simulate_chsh(PairSource(pair_rate_per_s=1e5, visibility=0.95), dict.fromkeys(IDS, det), 1.0, rng, SETTINGS)
        r = chsh(runs, Window(width_ps=2000, source="plan"), settings_deg=SETTINGS)
        assert not any(lab.startswith(prefix) for lab in r.all_labels)
        assert abs(r.raw.value - r.counts["capture_reweighted_s"]) < 0.25 * r.raw.stderr
    # Council round 3 (N1, G2): narrow peaks binned at W/30, and peaks as wide as the window fitted
    # only within it, gave offsets of +4 to +51 ps and false labels on matched delays.
    for pair_jitter, width, rate in ((30, 2000, 1e6), (30, 10_000, 1e6), (495, 800, 1e5)):
        det = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=pair_jitter / math.sqrt(2))
        for seed in range(3):
            runs = simulate_chsh(PairSource(pair_rate_per_s=rate, visibility=0.95), dict.fromkeys(IDS, det), 1.0,
                                 np.random.default_rng(200 + seed), SETTINGS)
            r = chsh(runs, Window(width_ps=width, source="plan"), settings_deg=SETTINGS)
            assert not any(lab.startswith((prefix, low)) for lab in r.all_labels), (pair_jitter, width, seed)
            assert max(abs(o) for o in r.counts["capture_offsets_ps"].values()) < 5.0 + 0.01 * pair_jitter
    # A window narrow against the jitter, delays matched: low capture, not a delay error (G5).
    runs = simulate_chsh(PairSource(pair_rate_per_s=1e5, visibility=0.95), dict.fromkeys(IDS, det), 1.0,
                         np.random.default_rng(200), SETTINGS)
    r = chsh(runs, Window(width_ps=600, source="plan"), settings_deg=SETTINGS)
    assert any(lab.startswith(low) for lab in r.all_labels)
    assert not any(lab.startswith(prefix) for lab in r.all_labels)


def test_a_delay_drift_since_calibration_is_refused(tmp_path):
    # Council round 2 (B1): delay artifacts taken before A+ and B+ drifted by 200 ps passed the
    # predictive check (it refits the delay); corrected S read 2.03-2.05 for 1.98.
    from dadl.errors import Misspecified
    from dadl.analysis.chsh import DELAY_DRIFT_ACCEPTED
    from dadl.forward_check import FORWARD_MISFIT_ACCEPTED

    rng = np.random.default_rng(3000)
    _calibrate(tmp_path, _late(0), rng)
    # Short pair runs leave some delay artifacts misfits; accepted here, the drift is the point.
    strict = _policy("accidental_subtraction", "window_capture", "misfit_calibration")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=strict, dataset=CTX)
    window = Window(width_ps=1000, source="plan")
    current = simulate_chsh(MISDECLARED_SRC, _late(0), 1.0, rng, SETTINGS)
    r = chsh(current, window, strict, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    assert FORWARD_MISFIT_ACCEPTED not in r.all_labels and DELAY_DRIFT_ACCEPTED not in r.all_labels
    drifted = simulate_chsh(MISDECLARED_SRC, _late(200), 1.0, rng, SETTINGS)
    with pytest.raises(Misspecified, match="calibrated delays"):
        chsh(drifted, window, strict, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    loose = strict.allow("chsh", "misfit_forward_model", "comparing a stale calibration")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=loose, dataset=CTX)
    r = chsh(drifted, window, loose, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    # Council round 8 (St7): the drift has its own label, apart from the predictive check's.
    assert DELAY_DRIFT_ACCEPTED in r.all_labels
    assert abs(r.counts["delay_drift_z"]) > 3.72
    # Council round 3: 100 ps at 150 ps jitter and a 1 ns window is refused too.
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=strict, dataset=CTX)
    drifted = simulate_chsh(MISDECLARED_SRC, _late(100), 1.0, rng, SETTINGS)
    with pytest.raises(Misspecified, match="calibrated delays"):
        chsh(drifted, window, strict, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)


@pytest.mark.slow
def test_current_delay_artifacts_are_not_refused(tmp_path_factory):
    # Council round 3 (G1): the drift threshold had no noise term and refused 2-65% of clean runs
    # (12/60 here: 150 ps jitter, 1 ns window, delays declared 0, artifacts current).
    from dadl.errors import Misspecified

    strict = _policy("accidental_subtraction", "window_capture", "misfit_calibration")
    window = Window(width_ps=1000, source="plan")
    zs = []
    for seed in range(12):
        rng = np.random.default_rng(5000 + seed)
        d = tmp_path_factory.mktemp(f"drift{seed}")
        _calibrate(d, _late(400), rng)
        cset = CalibrationSet.load(d, analysis="chsh", policy=strict, dataset=CTX)
        runs = simulate_chsh(MISDECLARED_SRC, _late(400), 1.0, rng, SETTINGS)
        try:
            r = chsh(runs, window, strict, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
        except Misspecified as exc:
            assert "calibrated delays" not in str(exc), seed
            continue
        zs.append(r.counts["delay_drift_z"])
    assert len(zs) >= 10 and max(abs(z) for z in zs) < 3.29


def test_a_corrected_s_without_dark_artifacts_says_so(tmp_path):
    # Council round 3 (R3): with no dark artifacts the label said dead time and window capture were
    # modelled, on a corrected S that kept the whole dead-time selection.
    from dadl.analysis.chsh import DEAD_TIME_UNCORRECTED
    from dadl.guardrails import DEAD_TIME_SELECTION

    det = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=1_000_000, jitter_sigma_ps=350.0)
    runs = simulate_chsh(SRC, dict.fromkeys(IDS, det), 0.2, np.random.default_rng(4), SETTINGS)
    policy = _policy("accidental_subtraction")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=policy, dataset=CTX)
    r = chsh(runs, Window(width_ps=2000, source="plan"), policy, settings_deg=SETTINGS, calibrations=cset, channels=IDS)
    assert r.corrections == ("accidental_subtraction",)
    label = chsh_model_label(r.corrections)
    assert label.startswith("accidentals modelled; dead time (no dark artifacts) and window capture") and label in r.all_labels
    assert DEAD_TIME_SELECTION in r.all_labels and DEAD_TIME_UNCORRECTED in r.all_labels
    assert r.counts["dataset"] == {"acquired": CTX.acquired.isoformat(), "conditions": {}}


# --- Council round 4 ----------------------------------------------------------------------------

def test_a_crosstalk_peak_fails_the_forward_check(tmp_path):
    # Council round 4 (S4): no test failed with the CHSH shape check removed.
    from dadl.errors import Misspecified
    from dadl.forward_check import FORWARD_MISFIT_ACCEPTED
    from dadl.sim import Crosstalk

    rng = np.random.default_rng(3000)
    _calibrate(tmp_path, _late(0), rng)
    strict = _policy("accidental_subtraction", "window_capture", "misfit_calibration")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=strict, dataset=CTX)
    runs = simulate_chsh(MISDECLARED_SRC, _late(0), 0.5, rng, SETTINGS, crosstalk=(Crosstalk("A+", "B+", 0.01, 600.0, 50.0),))
    window = Window(width_ps=1000, source="plan")
    with pytest.raises(Misspecified, match="does not fit the analysed coincidences"):
        chsh(runs, window, strict, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    loose = strict.allow("chsh", "misfit_forward_model", "checking")
    r = chsh(runs, window, loose, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    assert FORWARD_MISFIT_ACCEPTED in r.all_labels and r.counts["forward_check_p"] < 1e-3


def test_an_accidentals_only_corrected_s_runs_no_shape_check(tmp_path):
    # Council round 4 (S1): with no delay-jitter artifact CHSH still ran the shape check (HOM and
    # fringe skip it) and refused an IRF tail it cannot model without the calibrated jitter.
    det = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=150.0)
    rng = np.random.default_rng(11)
    runs = simulate_chsh(PairSource(pair_rate_per_s=2e5, visibility=0.7), dict.fromkeys(IDS, det), 1.0, rng, SETTINGS)
    for r in runs.values():
        for ch in ("B+", "B-"):
            t = r[ch].copy()
            late = rng.random(t.size) < 0.3
            t[late] += rng.exponential(400, late.sum()).astype(np.int64)
            r[ch] = np.unique(t)
    policy = _policy("accidental_subtraction")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=policy, dataset=CTX)
    r = chsh(runs, Window(1000, "plan"), policy, settings_deg=SETTINGS, calibrations=cset, channels=IDS)
    assert r.corrections == ("accidental_subtraction",) and "forward_check_p" not in r.counts


def _with_darks(cset, darks):
    """`cset`'s delay artifacts with these dark artifacts in place of its own."""
    from test_forward import _member_set

    from dadl.calibration.delay import KIND as DELAY_KIND

    delays = [m.calibration for m in cset.members if m.kind == DELAY_KIND]
    return _member_set("chsh", delays + darks)


def test_chsh_checks_the_dark_dead_times_against_the_analysed_gaps(data):
    # Council round 4 (S4, DT4): no test failed with the CHSH gap check removed. A 1 us dead time
    # on 50 ns detectors is refused downward, a 22 ns one upward.
    from test_forward import _dark

    from dadl.errors import StaleCalibration
    from dadl.forward import check_dead_time_gaps
    from dadl.propagate import OUTSIDE_VALIDITY

    runs, d = data
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    long = _with_darks(cset, [_dark(cid, 1e6, seed=i) for i, cid in enumerate(IDS.values())])
    with pytest.raises(StaleCalibration, match="analysed gaps are shorter"):
        chsh(runs, WINDOW, POLICY, calibrations=long, channels=IDS, seed=1)
    # Council round 9 (Rd2): the validity override no longer covers the gap tests.
    with pytest.raises(StaleCalibration, match="dead_time_contradicted"):
        check_dead_time_gaps(long, IDS, [runs[k] for k in runs], analysis="chsh",
                             policy=POLICY.allow("chsh", "calibration_outside_validity", "checking"))
    loose = POLICY.allow("chsh", "dead_time_contradicted", "checking")
    # Under the override the gap check passes it, labelled; its dead fraction (0.74) is then
    # refused (council round 6, R1 b).
    labels, counts = check_dead_time_gaps(long, IDS, [runs[k] for k in runs], analysis="chsh", policy=loose)
    # The 1 us point is also above the gaps: labelled (council round 8, DEAD_TIME_ABOVE_GAPS).
    rest = tuple(lab for lab in labels if not lab.startswith("the dark dead time of"))
    assert len(rest) == 1 and rest[0].startswith("dead time contradicted by the analysed gaps, accepted: dead time 1000 ns")
    assert OUTSIDE_VALIDITY not in labels and len(counts["dead_time_gap_problems"]) == 4
    assert 45_000 < counts["min_gap_ps"]["A+"] < 51_000
    with pytest.raises(GuardrailRefusal, match="pair term and the unmodelled paralyzable stage are not trusted"):
        chsh(runs, WINDOW, loose, calibrations=long, channels=IDS, seed=1)
    short = _with_darks(cset, [_dark(cid, 22_000, 500, seed=i) for i, cid in enumerate(IDS.values())])
    with pytest.raises(StaleCalibration, match="longer than its upper end"):
        chsh(runs, WINDOW, POLICY, calibrations=short, channels=IDS, seed=1)


def test_chsh_refuses_another_detector_s_artifacts(data):
    # Council round 4 (S3): the detector check ran only through CalibrationSet.load(channels=).
    import dataclasses

    from test_forward import _member_set

    from dadl.errors import StaleCalibration
    from dadl.propagate import OUTSIDE_VALIDITY

    runs, d = data
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    other = {n: ChannelId(c.device, c.input, role=n, detector=f"NEW-{n}") for n, c in IDS.items()}
    recorded = {n: ChannelId(c.device, c.input, role=n, detector=f"OLD-{n}") for n, c in IDS.items()}
    by_channel = {c: n for n, c in IDS.items()}
    cals = [dataclasses.replace(m.calibration, channels=tuple(recorded[by_channel[c]] for c in m.channels))
            for m in cset.members]
    old = _member_set("chsh", cals)
    with pytest.raises(StaleCalibration, match="calibrated on detector OLD-"):
        chsh(runs, WINDOW, POLICY, calibrations=old, channels=other, seed=1)
    loose = POLICY.allow("chsh", "calibration_outside_validity", "same detectors, renamed")
    assert OUTSIDE_VALIDITY in chsh(runs, WINDOW, loose, calibrations=old, channels=other, seed=1).all_labels
    assert OUTSIDE_VALIDITY not in chsh(runs, WINDOW, POLICY, calibrations=old, channels=recorded, seed=1).all_labels


def test_the_capture_balance_carries_the_peak_fits_noise():
    # Council round 4 (S4): no test failed with capture_sd_data set to 0. It shrinks as the run
    # grows (about 1/sqrt(duration)).
    import dataclasses

    det = dataclasses.replace(MISDECLARED_DET, jitter_sigma_ps=150.0)
    sds = []
    for duration in (0.25, 1.0):
        runs = simulate_chsh(MISDECLARED_SRC, dict.fromkeys(IDS, det), duration, np.random.default_rng(40), SETTINGS)
        sds.append(chsh(runs, WINDOW, settings_deg=SETTINGS).counts["capture_sd_data"])
    assert sds[0] > 0 and 0.3 < sds[1] / sds[0] < 0.7


def test_the_peak_fits_noise_can_decide_the_capture_label():
    # A shift of 0.6 statistical se, within 3.29 sd of the peak fits: not labelled; with no fit
    # noise it would be.
    from dadl.analysis.chsh import CAPTURE_IMBALANCE, capture_labels
    from dadl.result import Estimate

    raw = Estimate(2.0, 0.02)
    balance = {(x, y): (0.0, 300.0, 0.9, 5.0, 5.0) for x in "+-" for y in "+-"}
    assert capture_labels(raw, 2.012, balance, 0.005) == ()
    assert capture_labels(raw, 2.012, balance, 0.0)[0] == CAPTURE_IMBALANCE.format(dS=-0.012)


def test_labels_of_unused_artifacts_are_left_out(data, tmp_path):
    # Council round 6 (R2): the CHSH raw line carried "no validity range" for accidentals and
    # heralding artifacts the corrected value did not use.
    import shutil

    from dadl.calibration import calibrate_accidentals
    from dadl.calibration.delay import KIND as DELAY_KIND
    from dadl.propagate import NO_VALIDITY, NO_VALIDITY_MANY

    runs, d = data
    full = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    for m in full.members:
        if m.kind == DELAY_KIND:
            for f in d.glob(f"{m.artifact_id}*"):
                shutil.copy(f, tmp_path)
    a, b = simulate_pairs(SRC.pair_rate_per_s, DET, DET, 0.02, np.random.default_rng(9))
    unused = write_artifact(calibrate_accidentals(a, b, (IDS["A+"], IDS["B+"]), width_ps=WINDOW.width_ps, backend="sim",
                                                  seed=1), tmp_path, created=dt.datetime(2026, 11, 1, tzinfo=UTC))
    policy = _policy("accidental_subtraction", "window_capture")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=policy, dataset=CTX)
    r = chsh(runs, WINDOW, policy, calibrations=cset, channels=IDS, seed=1)
    gone = NO_VALIDITY.format(aid=unused.stem[:12])
    assert gone in cset.labels and gone not in r.all_labels
    # Council round 9 (Tw8): the used artifacts' "no validity range" is one label listing them.
    used = NO_VALIDITY_MANY.format(aids=", ".join(aid[:12] for aid in r.artifact_ids))
    assert len(r.artifact_ids) == 4 and used in r.corrected_labels and used not in r.labels
    assert FORWARD_MODEL in r.corrected_labels and FORWARD_MODEL not in r.labels


def test_chsh_takes_delay_artifacts_in_either_channel_order(tmp_path):
    # Council round 7 (Q2): B+ and B- 400 ps late, declared; artifacts written (B, A) give the
    # corrected S of artifacts written (A, B).
    from test_forward import SIGN_DET, SIGN_LATE, assert_same_correction, sign_sets

    rng = np.random.default_rng(32)
    dets = {"A+": SIGN_DET, "A-": SIGN_DET, "B+": SIGN_LATE, "B-": SIGN_LATE}
    pairs = {(IDS["A" + x], IDS["B" + y]): simulate_pairs(2e5, dets["A" + x], dets["B" + y], 0.05, rng)
             for x in "+-" for y in "+-"}
    policy = _policy("accidental_subtraction", "window_capture")
    sets = sign_sets(tmp_path, "chsh", policy, pairs)
    runs = simulate_chsh(PairSource(pair_rate_per_s=2e5, visibility=0.9), dets, 0.2, rng, SETTINGS)
    window = Window(width_ps=200, source="plan", delays_ps={"B+": 400, "B-": 400})
    ab, ba = (chsh(runs, window, policy, settings_deg=SETTINGS, calibrations=c, channels=IDS, seed=2) for c in sets)
    assert ab.corrections == ("accidental_subtraction", "window_capture")
    assert_same_correction((ab, ab.counts), (ba, ba.counts))



# --- Council round 7 (Q3): CHSH labels and refusals through chsh() --------------------------------


def test_chsh_labels_a_dark_dead_time_whose_lower_end_the_gaps_cross(data):
    # 100 ns with a 95% lower end of 66 ns on 50 ns detectors: three widths below the point is under
    # zero (not tested), the gaps fall below the lower end, and the dead fraction (0.08) is allowed.
    from test_forward import _dark

    from dadl.forward import DEAD_TIME_BELOW_DARK

    runs, d = data
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    loose = _with_darks(cset, [_dark(cid, 100_000, 17_000, seed=i) for i, cid in enumerate(IDS.values())])
    r = chsh(runs, WINDOW, POLICY, settings_deg=SETTINGS, calibrations=loose, channels=IDS, seed=1)
    head = DEAD_TIME_BELOW_DARK.split("{name}")[0]
    assert any(lab.startswith(head) for lab in r.corrected_labels)
    assert r.counts["gap_check_floor_kind"]["A+"] == "lower end (label only)"


def test_chsh_refuses_an_uninformative_dead_time_bound(data):
    # The CHSH twin of the HOM test: a bound of 1 us at 750 kcps is a dead fraction of 0.75.
    from test_forward import _dark

    runs, d = data
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    bounded = _with_darks(cset, [_dark(cid, bound=1e6, seed=i) for i, cid in enumerate(IDS.values())])
    with pytest.raises(GuardrailRefusal, match="known only below"):
        chsh(runs, WINDOW, POLICY, settings_deg=SETTINGS, calibrations=bounded, channels=IDS, seed=1)


def test_chsh_threshold_labels_belong_to_the_corrected_value(data):
    from dadl.guardrails import NEAR_THRESHOLD_EXACT, THRESHOLD_AFTER_CORRECTION

    _, d = data
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    # A source at S = 2.000 and a short run: the corrected interval holds the threshold at this seed
    # (a 95% event; the label rule below holds at any seed).
    src = PairSource(pair_rate_per_s=SRC.pair_rate_per_s, visibility=2**-0.5)
    runs = simulate_chsh(src, dict.fromkeys(IDS, DET), 0.02, np.random.default_rng(2), SETTINGS)
    r = chsh(runs, WINDOW, POLICY, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    c = r.corrected
    assert (NEAR_THRESHOLD_EXACT in r.corrected_labels) == (c.low <= 2.0 <= c.high)
    assert (THRESHOLD_AFTER_CORRECTION in r.corrected_labels) == (r.raw.value < 2.0 <= c.value)
    assert NEAR_THRESHOLD_EXACT not in r.labels
    assert c.low <= 2.0 <= c.high


# --- Council round 8 ----------------------------------------------------------------------------


def test_a_drift_below_the_refusal_threshold_is_not_refused(monkeypatch):
    # Council round 8 (Fix 4): the drift refuses at max(0.5 se, 3.72 sd) (two-sided 2e-4); a drift
    # at z 3.5, above the 0.5 se floor, is no longer refused. The capture label keeps 3.29 sd.
    import importlib

    import dadl.forward as F
    from dadl.errors import Misspecified
    from dadl.forward import OUTCOME_SIGNS
    from dadl.result import Estimate

    C = importlib.import_module("dadl.analysis.chsh")
    assert C.DRIFT_Z == 3.72 and C.CAPTURE_Z == 3.29
    monkeypatch.setattr(F, "_pair_capture", lambda *a, **k: (np.full(10, 0.9), 0.9, False))
    monkeypatch.setattr(C, "reweighted_s", lambda counts, caps: 2.0 if isinstance(caps[("+", "+")], float)
                        else np.full(10, 2.0))
    balance = dict.fromkeys(OUTCOME_SIGNS, (0.0, 100.0, 0.9, 1.0, 1.0))
    corrected = Estimate(2.0, 0.01, stderr_statistical=0.01)
    window = Window(width_ps=1000, source="plan")
    labels, counts = C._delay_drift(window, Policy.default(), None, IDS, 1, ({}, 2.035, balance, 0.01), corrected)
    assert labels == () and counts["delay_drift_z"] == pytest.approx(3.5)
    with pytest.raises(Misspecified, match="calibrated delays"):
        C._delay_drift(window, Policy.default(), None, IDS, 1, ({}, 2.038, balance, 0.01), corrected)
    loose = Policy.default().allow("chsh", "misfit_forward_model", "test")
    labels, _ = C._delay_drift(window, loose, None, IDS, 1, ({}, 2.038, balance, 0.01), corrected)
    assert labels == (C.DELAY_DRIFT_ACCEPTED,)


@pytest.mark.slow
def test_clean_chsh_spad_runs_pass_the_off_peak_check():
    # Council round 8 (Ph3, Fix 5): the off-peak check on the 16 histograms of clean chsh-spad
    # runs (afterpulsing 0.5% at 50 us included) does not flag, over 20 runs.
    from dadl.analysis._twopoint import OFFPEAK_ALPHA, offpeak_p
    from dadl.analysis.propagation import CHSH_SETTINGS, REGIMES
    from dadl.calibration.accidentals import shifted_coincidences, side_offsets

    reg = REGIMES["chsh-spad"]
    offsets = side_offsets(reg.window_ps)
    for seed in range(20):
        runs = simulate_chsh(reg.source, dict.fromkeys(("A+", "A-", "B+", "B-"), reg.detector), reg.duration_s,
                             np.random.default_rng(900 + seed), CHSH_SETTINGS)
        ps = [offpeak_p(shifted_coincidences(np.sort(r["A" + x]), np.sort(r["B" + y]), reg.window_ps, 0, 0, offsets),
                        offsets, 16) for r in runs.values() for x in "+-" for y in "+-"]
        assert min(ps) >= OFFPEAK_ALPHA / len(ps)


# --- Council round 9 -----------------------------------------------------------------------------


def _flat_sides(rng, n, mean=1.0):
    from dadl.calibration.accidentals import side_offsets

    offsets = side_offsets(2_000)
    return offsets, {f"s{i}": rng.poisson(mean, offsets.size) for i in range(n)}


def test_a_periodic_off_peak_side_among_16_is_refused():
    # Council round 9 (Pf1): at 20,000 fixed draws the doubled dispersion p could never reach
    # 2e-4 / 16; the r9-perf case (every 4th window, mean 1) passed CHSH unflagged.
    from dadl.analysis._twopoint import offpeak_check
    from dadl.errors import Misspecified

    rng = np.random.default_rng(3)
    offsets, sides = _flat_sides(rng, 15)
    bad = np.zeros(offsets.size, dtype=np.int64)
    bad[::4] = rng.poisson(4.0, bad[::4].size)
    sides = {"bad": bad, **sides}
    with pytest.raises(Misspecified, match="bad: p = "):
        offpeak_check("chsh", sides, offsets, None)
    p, labels, draws = offpeak_check("chsh", sides, offsets, POLICY.allow("chsh", "misfit_forward_model", "checking"))
    assert p < 2e-4 and labels and draws > 640_000
    # The flagged side runs to n_max; the sides after it stop once they cannot be flagged.
    assert draws < 700_000


def test_the_sequential_off_peak_p_is_reproducible_and_flat_sides_pass():
    # Council round 9 (Pf1): seeded from the counts, the same data give the same p; a flat side
    # stops after a few hundred draws with p = 20 / draws, not near the threshold.
    from dadl.analysis._twopoint import OFFPEAK_ALPHA, offpeak_check, sequential_dispersion_pvalue

    rng = np.random.default_rng(11)
    offsets, sides = _flat_sides(rng, 16, mean=1.3)
    first = [sequential_dispersion_pvalue(s, 16) for s in sides.values()]
    assert first == [sequential_dispersion_pvalue(s, 16) for s in sides.values()]
    for p, draws in first:
        assert p >= 1e-3 and (draws == 0 or p == pytest.approx(20 / draws) or draws == 640_000)
    p, labels, draws = offpeak_check("chsh", sides, offsets, None)
    assert p >= OFFPEAK_ALPHA and labels == () and 0 < draws < 16 * 20_000
    # The chi-square path above a mean of 5 per window draws nothing.
    assert sequential_dispersion_pvalue(rng.poisson(8.0, offsets.size), 16)[1] == 0


@pytest.mark.slow
def test_pulsed_accidentals_at_half_of_one_chsh_side_are_refused():
    # Council round 9 (Pf1, pulsed_chsh.py): one of 16 sides holds a pulsed share of 0.5 at a
    # 13157.9 ps period (50 ps jitter), mean 1.3 per window; refused in at least 28 of 30 seeds.
    from scipy.special import erf

    from dadl.analysis._twopoint import offpeak_check
    from dadl.calibration.accidentals import side_offsets
    from dadl.errors import Misspecified

    width, period, jit = 2_000, 13_157.9, 50.0
    offsets = side_offsets(width)
    w = np.zeros(offsets.size)
    for i, o in enumerate(offsets):
        lo, hi = o, o + width
        for n in range(int(np.floor((lo - 5 * jit) / period)), int(np.ceil((hi + 5 * jit) / period)) + 1):
            w[i] += 0.5 * (erf((hi - n * period) / (np.sqrt(2) * jit)) - erf((lo - n * period) / (np.sqrt(2) * jit)))
    lam = 1.3 * (0.5 + 0.5 * w / w.mean())
    refused = 0
    for seed in range(30):
        rng = np.random.default_rng(1000 + seed)
        sides = {f"s{i}": rng.poisson(1.3, offsets.size) for i in range(15)}
        sides["pulsed"] = rng.poisson(lam)
        try:
            offpeak_check("chsh", sides, offsets, None)
        except Misspecified:
            refused += 1
    assert refused >= 28


@pytest.mark.slow
def test_the_sequential_dispersion_p_is_valid_on_iid_poisson_sides():
    # Council round 9 (Pf1, null_seed.txt): 20,000 iid Poisson(1.3) sides at 16 sides; with the
    # seed from the counts every level is at or below nominal (seed 0 ran the 5% level at 6%).
    from dadl.analysis._twopoint import sequential_dispersion_pvalue
    from dadl.calibration.accidentals import side_offsets

    rng = np.random.default_rng(20261120)
    k = side_offsets(2_000).size
    ps = np.array([sequential_dispersion_pvalue(rng.poisson(1.3, k), 16)[0] for _ in range(20_000)])
    assert np.mean(ps <= 0.01) <= 0.012 and np.mean(ps <= 1e-3) <= 0.0015


def test_chsh_and_the_gap_check_do_not_modify_the_caller_s_tags(data):
    # Council round 9 (Pf3): as_sorted hands sorted input through without a copy.
    import copy

    from dadl.forward import check_dead_time_gaps

    runs, d = data
    before = copy.deepcopy(runs)
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    chsh(runs, WINDOW, POLICY, calibrations=cset, channels=IDS, seed=1)
    check_dead_time_gaps(cset, IDS, [runs[k] for k in runs], analysis="chsh", policy=POLICY)
    assert all(np.array_equal(runs[k][c], before[k][c]) for k in runs for c in runs[k])
