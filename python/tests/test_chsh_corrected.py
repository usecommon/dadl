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
from dadl.forward import ACCIDENTALS_MODELLED, FORWARD_MODEL
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
    for label in (FORWARD_MODEL, ACCIDENTALS_MODELLED, POLICY_MODIFIED):
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
