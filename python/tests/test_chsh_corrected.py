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


def test_corrected_s_removes_the_dead_time_selection(data):
    runs, d = data
    cset = CalibrationSet.load(d, analysis="chsh", policy=POLICY, dataset=CTX)
    r = chsh(runs, WINDOW, POLICY, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    raw_only = chsh(runs, WINDOW, settings_deg=SETTINGS)
    assert r.raw == raw_only.raw  # the headline is untouched
    assert r.raw.value - 1.96 * r.raw.stderr > TRUTH  # raw S excludes the source's S (dead-time selection)
    assert r.corrected.low <= TRUTH <= r.corrected.high
    assert r.corrected.value < r.raw.value
    assert r.corrections == ("dead_time", "accidental_subtraction", "window_capture")
    assert len(r.artifact_ids) == 8
    for label in (FORWARD_MODEL, chsh_model_label(r.corrections), POLICY_MODIFIED):
        assert label in r.labels
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
    assert DEAD_TIME_BOUND in r.labels and DEAD_TIME_RANGE_WIDE in r.labels
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
    assert any(lab.startswith(prefix) for lab in r.labels)
    assert abs(r.counts["capture_reweighted_s"] - truth) < 3 * r.raw.stderr
    assert r.counts["capture_offsets_ps"]["A+/B-"] == pytest.approx(-400, abs=15)
    right = chsh(runs, Window(width_ps=1000, source="plan", delays_ps={"A+": 400, "B+": 400}), settings_deg=SETTINGS)
    assert not any(lab.startswith(prefix) for lab in right.labels)


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
        assert not any(lab.startswith(prefix) for lab in r.labels)
        assert abs(r.raw.value - r.counts["capture_reweighted_s"]) < 0.25 * r.raw.stderr
    # Council round 3 (N1, G2): narrow peaks binned at W/30, and peaks as wide as the window fitted
    # only within it, gave offsets of +4 to +51 ps and false labels on matched delays.
    for pair_jitter, width, rate in ((30, 2000, 1e6), (30, 10_000, 1e6), (495, 800, 1e5)):
        det = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=pair_jitter / math.sqrt(2))
        for seed in range(3):
            runs = simulate_chsh(PairSource(pair_rate_per_s=rate, visibility=0.95), dict.fromkeys(IDS, det), 1.0,
                                 np.random.default_rng(200 + seed), SETTINGS)
            r = chsh(runs, Window(width_ps=width, source="plan"), settings_deg=SETTINGS)
            assert not any(lab.startswith((prefix, low)) for lab in r.labels), (pair_jitter, width, seed)
            assert max(abs(o) for o in r.counts["capture_offsets_ps"].values()) < 5.0 + 0.01 * pair_jitter
    # A window narrow against the jitter, delays matched: low capture, not a delay error (G5).
    runs = simulate_chsh(PairSource(pair_rate_per_s=1e5, visibility=0.95), dict.fromkeys(IDS, det), 1.0,
                         np.random.default_rng(200), SETTINGS)
    r = chsh(runs, Window(width_ps=600, source="plan"), settings_deg=SETTINGS)
    assert any(lab.startswith(low) for lab in r.labels)
    assert not any(lab.startswith(prefix) for lab in r.labels)


def test_a_delay_drift_since_calibration_is_refused(tmp_path):
    # Council round 2 (B1): delay artifacts taken before A+ and B+ drifted by 200 ps passed the
    # predictive check (it refits the delay); corrected S read 2.03-2.05 for 1.98.
    from dadl.errors import Misspecified
    from dadl.forward_check import FORWARD_MISFIT_ACCEPTED

    rng = np.random.default_rng(3000)
    _calibrate(tmp_path, _late(0), rng)
    # Short pair runs leave some delay artifacts misfits; accepted here, the drift is the point.
    strict = _policy("accidental_subtraction", "window_capture", "misfit_calibration")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=strict, dataset=CTX)
    window = Window(width_ps=1000, source="plan")
    current = simulate_chsh(MISDECLARED_SRC, _late(0), 1.0, rng, SETTINGS)
    r = chsh(current, window, strict, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    assert FORWARD_MISFIT_ACCEPTED not in r.labels
    drifted = simulate_chsh(MISDECLARED_SRC, _late(200), 1.0, rng, SETTINGS)
    with pytest.raises(Misspecified, match="calibrated delays"):
        chsh(drifted, window, strict, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    loose = strict.allow("chsh", "misfit_forward_model", "comparing a stale calibration")
    cset = CalibrationSet.load(tmp_path, analysis="chsh", policy=loose, dataset=CTX)
    r = chsh(drifted, window, loose, settings_deg=SETTINGS, calibrations=cset, channels=IDS, seed=1)
    assert FORWARD_MISFIT_ACCEPTED in r.labels
    assert abs(r.counts["delay_drift_z"]) > 3.29
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
    assert label.startswith("accidentals modelled; dead time (no dark artifacts) and window capture") and label in r.labels
    assert DEAD_TIME_SELECTION in r.labels and DEAD_TIME_UNCORRECTED in r.labels
    assert r.counts["dataset"] == {"acquired": CTX.acquired.isoformat(), "conditions": {}}
