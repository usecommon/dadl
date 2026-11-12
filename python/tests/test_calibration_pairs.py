"""Delay/jitter, accidentals, heralding efficiency and DCR(T) calibrations (milestone 2 of phase 2)."""

import json
import math
from pathlib import Path

import numpy as np
import pytest

from dadl.artifact import build_manifest
from dadl.calibration import (
    calibrate_accidentals,
    calibrate_dark,
    calibrate_dcr_temperature,
    calibrate_delay,
    calibrate_heralding,
)
from dadl.calibration.accidentals import shifted_coincidences, side_offsets
from dadl.calibration.coverage import (
    PAIR_A,
    PAIR_B,
    PAIR_DURATION_S,
    PAIR_RATE,
    PAIR_WIDTH_PS,
    dcr_coverage,
    pairs_coverage,
    pairs_truth,
)
from dadl.closed_form import arrhenius_dark_rate
from dadl.coincidences import grid_coincidences
from dadl.data import ChannelId
from dadl.errors import InsufficientCounts
from dadl.sim import DetectorParams, simulate_dark, simulate_pairs

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "schema" / "calibration-1.schema.json"
CHS = (ChannelId("sim", 0), ChannelId("sim", 1))


@pytest.fixture(scope="module")
def pairs():
    return simulate_pairs(PAIR_RATE, PAIR_A, PAIR_B, PAIR_DURATION_S, np.random.default_rng(7))


@pytest.fixture(scope="module")
def truth():
    return pairs_truth()


# --- delay and jitter ------------------------------------------------------------------------


def test_delay_recovers_injected_delay_and_jitter(pairs, truth):
    cal = calibrate_delay(*pairs, CHS, backend="sim", seed=1)
    for k, v in truth["delay-jitter"].items():
        assert cal.estimates[k].contains(v), (k, cal.estimates[k], v)
    assert cal.converged and cal.ppc_pvalue > 0.001
    assert cal.samples.shape == (2000, 2)
    assert cal.kind == "delay-jitter" and cal.channels == CHS


def test_delay_sign_is_b_minus_a(pairs):
    a, b = pairs
    cal = calibrate_delay(b, a, CHS[::-1], seed=1)
    assert cal.estimates["delay_ps"].contains(-(PAIR_B.delay_ps - PAIR_A.delay_ps))


def test_delay_flags_an_exponential_tail():
    # A detector with a long exponential IRF tail must fail the Gaussian predictive check.
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0, 2e12, 40_000))
    a = np.rint(t + rng.normal(0, 50, t.size)).astype(np.int64)
    b = np.sort(np.rint(t + 1000 + rng.normal(0, 50, t.size) + rng.exponential(300, t.size)).astype(np.int64))
    cal = calibrate_delay(a, b, CHS, seed=1)
    assert cal.ppc_pvalue < 0.01


def test_delay_needs_a_peak():
    rng = np.random.default_rng(0)
    a = np.sort(rng.integers(0, 10**12, 2000))
    b = np.sort(rng.integers(0, 10**12, 2000))
    with pytest.raises(InsufficientCounts):
        calibrate_delay(a, b, CHS)


def test_delay_peak_at_the_range_edge_is_not_converged(pairs):
    cal = calibrate_delay(*pairs, CHS, range_ps=3_500, seed=1)
    assert not cal.converged
    assert any("edge" in n for n in cal.notes)


# --- accidentals ---------------------------------------------------------------------------


def test_side_offsets():
    cw = side_offsets(1000, n_side=3, guard_ps=2500)
    assert cw.tolist() == [-8000, -7000, -6000, 3000, 4000, 5000]
    pulsed = side_offsets(1000, rep_period_ps=12_500, n_side=2, guard_ps=20_000)
    assert pulsed.tolist() == [-62_500, -50_000, 25_000, 37_500]
    # No window is the mirror image of another (mirrored windows are correlated).
    assert not set(cw.tolist()) & set((-cw).tolist())


def test_shifted_coincidences_match_direct_counts(pairs):
    a, b = pairs
    offsets = np.concatenate([side_offsets(PAIR_WIDTH_PS, n_side=5), [0, 12_345, -7_001]])
    got = shifted_coincidences(a, b, PAIR_WIDTH_PS, 0, PAIR_B.delay_ps, offsets)
    want = [grid_coincidences(a, b, PAIR_WIDTH_PS, 0, PAIR_B.delay_ps + int(o)) for o in offsets]
    assert got.tolist() == want


def test_accidentals_match_closed_form(pairs, truth):
    cal = calibrate_accidentals(*pairs, CHS, width_ps=PAIR_WIDTH_PS, delay_b_ps=PAIR_B.delay_ps, seed=1)
    assert cal.estimates["accidental_coincidences"].contains(truth["accidentals"]["accidental_coincidences"])
    assert cal.conditions["width_ps"] == PAIR_WIDTH_PS
    assert any("cross-check" in n for n in cal.notes)


def test_accidentals_large_sample_agrees_with_closed_form():
    # 20 s of data: the estimate must sit within 4% of R1 R2 width duration.
    from dadl.closed_form import accidentals_expected

    a, b = simulate_pairs(PAIR_RATE, PAIR_A, PAIR_B, 20.0, np.random.default_rng(11))
    cal = calibrate_accidentals(a, b, CHS, width_ps=PAIR_WIDTH_PS, delay_b_ps=PAIR_B.delay_ps, n_side=1000)
    exp = accidentals_expected(
        PAIR_A.efficiency * PAIR_RATE + PAIR_A.dark_rate_cps,
        PAIR_B.efficiency * PAIR_RATE + PAIR_B.dark_rate_cps,
        PAIR_A.dead_time_ps * 1e-12,
        PAIR_B.dead_time_ps * 1e-12,
        PAIR_WIDTH_PS,
        20.0,
    )
    assert cal.estimates["accidental_coincidences"].value == pytest.approx(exp, rel=0.04)


def test_periodic_background_fails_the_dispersion_check():
    # Background bunched at a 50 ns period, analysed with CW offsets: shifted counts are over-dispersed.
    rng = np.random.default_rng(5)
    pulses = np.arange(0, 2 * 10**12, 50_000)
    a = np.sort(rng.choice(pulses, 200_000, replace=False) + rng.integers(0, 300, 200_000))
    b = np.sort(rng.choice(pulses, 200_000, replace=False) + rng.integers(0, 300, 200_000))
    cal = calibrate_accidentals(a, b, CHS, width_ps=1000, n_side=60)
    assert cal.ppc_pvalue < 0.01
    pulsed = calibrate_accidentals(a, b, CHS, width_ps=1000, rep_period_ps=50_000, n_side=60)
    assert pulsed.estimates["accidental_coincidences"].value > 10 * cal.estimates["accidental_coincidences"].value


# --- heralding efficiency ------------------------------------------------------------------


def test_heralding_matches_closed_form(pairs, truth):
    cal = calibrate_heralding(
        *pairs, CHS, width_ps=PAIR_WIDTH_PS, delay_signal_ps=PAIR_B.delay_ps, backend="sim", seed=1
    )
    iv = cal.estimates["heralding_efficiency"]
    assert iv.contains(truth["heralding-efficiency"]["heralding_efficiency"])
    assert 0 < iv.low < iv.value < iv.high < 1
    lo, hi = np.quantile(cal.samples[:, 0], [0.025, 0.975])
    assert lo == pytest.approx(iv.low, abs=0.15 * (iv.high - iv.low))
    assert hi == pytest.approx(iv.high, abs=0.15 * (iv.high - iv.low))


def test_heralding_without_the_declared_delay_is_low(pairs):
    cal = calibrate_heralding(*pairs, CHS, width_ps=PAIR_WIDTH_PS)
    assert cal.estimates["heralding_efficiency"].value < 0.2


def test_heralding_needs_heralds():
    with pytest.raises(InsufficientCounts):
        calibrate_heralding(np.arange(50), np.arange(50), CHS, width_ps=1000)


# --- DCR(T) ----------------------------------------------------------------------------------


def _darks(temps, ea=0.5, seconds=10.0):
    out = []
    for t_k in temps:
        det = DetectorParams(dark_rate_cps=arrhenius_dark_rate(200.0, ea, t_k, 273.0))
        tags = simulate_dark(det, seconds, np.random.default_rng(int(t_k)))
        out.append(calibrate_dark(tags, CHS[0], backend="sim", seed=1, n_draws=50, conditions={"temperature_K": t_k}))
    return out


def test_dcr_recovers_arrhenius_parameters():
    cal = calibrate_dcr_temperature(_darks((253.0, 263.0, 273.0, 283.0, 293.0)), reference_k=273.0, seed=1)
    assert cal.estimates["dcr_ref_cps"].contains(200.0)
    assert cal.estimates["activation_energy_ev"].contains(0.5)
    assert cal.conditions["reference_K"] == 273.0
    assert len(cal.inputs) == 5 and cal.backend == "sim"


def test_dcr_default_reference_decorrelates():
    cal = calibrate_dcr_temperature(_darks((253.0, 273.0, 293.0)), seed=1)
    r = np.corrcoef(np.log(cal.samples[:, 0]), cal.samples[:, 1])[0, 1]
    assert abs(r) < 0.1
    assert 253.0 < cal.conditions["reference_K"] < 293.0


def test_dcr_refuses_bad_inputs():
    darks = _darks((253.0, 273.0, 293.0))
    with pytest.raises(ValueError, match="3 or more"):
        calibrate_dcr_temperature(darks[:2])
    other = calibrate_dark(
        simulate_dark(DetectorParams(dark_rate_cps=200.0), 10.0, np.random.default_rng(0)),
        CHS[1], seed=1, n_draws=50, conditions={"temperature_K": 300.0},
    )
    with pytest.raises(ValueError, match="same channel"):
        calibrate_dcr_temperature([*darks, other])
    no_temp = calibrate_dark(
        simulate_dark(DetectorParams(dark_rate_cps=200.0), 10.0, np.random.default_rng(0)), CHS[0], seed=1, n_draws=50
    )
    with pytest.raises(ValueError, match="temperature_K"):
        calibrate_dcr_temperature([*darks, no_temp])


# --- artifacts -------------------------------------------------------------------------------


@pytest.mark.skipif(not SCHEMA.exists(), reason="repository schema/ not present")
def test_every_new_kind_validates_against_the_schema(pairs):
    import datetime as dt

    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    validator = jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker())
    created = dt.datetime(2026, 11, 11, tzinfo=dt.timezone.utc)
    a, b = pairs
    cals = [
        calibrate_delay(a, b, CHS, seed=1, n_draws=20),
        calibrate_accidentals(a, b, CHS, width_ps=PAIR_WIDTH_PS, delay_b_ps=PAIR_B.delay_ps, seed=1, n_draws=20),
        calibrate_heralding(a, b, CHS, width_ps=PAIR_WIDTH_PS, delay_signal_ps=PAIR_B.delay_ps, seed=1, n_draws=20),
        calibrate_dcr_temperature(_darks((253.0, 273.0, 293.0)), seed=1, n_draws=20),
    ]
    for cal in cals:
        manifest, _ = build_manifest(cal, created=created)
        validator.validate(manifest)


# --- the phase 2 gate ------------------------------------------------------------------------


def test_closed_form_truths_are_sane(truth):
    eta = truth["heralding-efficiency"]["heralding_efficiency"]
    assert 0.3 * 0.7 < eta < PAIR_B.efficiency
    assert truth["delay-jitter"]["jitter_sigma_ps"] == pytest.approx(math.hypot(350, 200))


@pytest.mark.coverage
def test_pair_calibrations_coverage():
    results = pairs_coverage(1000)
    for r in results:
        print(r)
    assert all(r.passed() for r in results), "\n".join(map(str, results))


@pytest.mark.coverage
def test_dcr_calibration_coverage():
    result = dcr_coverage(1000)
    print(result)
    assert result.passed(), str(result)


def test_delay_is_unbiased_for_narrow_peaks():
    # Tags are whole ps: method v1 read bin [e, e + w) as continuous and was biased by +0.5 ps,
    # visible when the jitter is a few bin widths (here sigma = 35 ps, bins of 20 ps).
    det = DetectorParams(efficiency=0.85, dark_rate_cps=100, dead_time_ps=30_000, jitter_sigma_ps=25)
    est = []
    for seed in range(10):
        a, b = simulate_pairs(1e5, det, det, 1.0, np.random.default_rng(seed))
        est.append(calibrate_delay(a, b, CHS, n_draws=10).estimates["delay_ps"].value)
    assert abs(np.mean(est)) < 0.2


def test_heralds_sharing_a_slot_are_counted_once():
    # Two pairs per 2 ns slot, each herald with its own signal: eta is 1 per herald slot.
    t = np.arange(0, 10**10, 10**6)
    herald = np.sort(np.concatenate([t, t + 500]))
    signal = np.sort(np.concatenate([t + 100, t + 600]))
    cal = calibrate_heralding(herald, signal, CHS, width_ps=2000)
    assert cal.estimates["heralding_efficiency"].value == pytest.approx(1.0)
    assert any("shared a slot" in n for n in cal.notes)


def test_sparse_sharp_peak_still_gets_a_predictive_check_and_an_artifact(tmp_path):
    # Too few counts for a grouped chi-square: v2 returned NaN, which no artifact can hold.
    import datetime as dt

    from dadl.artifact import write_artifact

    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0, 1e12, 150))
    a = np.rint(t).astype(np.int64)
    b = np.sort(np.rint(t + rng.normal(0, 8, t.size)).astype(np.int64))
    cal = calibrate_delay(a, b, CHS, n_draws=10)
    assert 0.0 < cal.ppc_pvalue <= 1.0
    write_artifact(cal, tmp_path, created=dt.datetime(2026, 11, 12, tzinfo=dt.timezone.utc))
