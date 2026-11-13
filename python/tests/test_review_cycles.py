"""Findings of the outside-expert review cycles (docs/reviews/2026-11-review-cycles.md)."""

import math

import numpy as np
import pytest

from dadl.cli import main
from dadl.data import ChannelId
from dadl.errors import ConfigError
from dadl.sim import DetectorParams, simulate_pairs

CHS = (ChannelId("sim", 0), ChannelId("sim", 1))
SNSPD_A = DetectorParams(efficiency=0.85, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=25.0, delay_ps=0)
SNSPD_B = DetectorParams(efficiency=0.85, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=25.0, delay_ps=500)


def latch(tags: np.ndarray, start: float, fraction: float) -> np.ndarray:
    t0, span = int(tags[0]), int(tags[-1] - tags[0])
    lo, hi = t0 + int(start * span), t0 + int((start + fraction) * span)
    return tags[(tags < lo) | (tags > hi)]


def test_interruption_test_is_quiet_on_clean_streams_and_catches_a_latch():
    from dadl.stationarity import interruption_pvalue

    rng = np.random.default_rng(1)
    clean = [interruption_pvalue(simulate_pairs(2e4, SNSPD_A, SNSPD_B, 0.5, rng)[1])[0] for _ in range(20)]
    assert min(clean) > 1e-3
    _, b = simulate_pairs(2e4, SNSPD_A, SNSPD_B, 0.5, rng)
    assert interruption_pvalue(latch(b, 0.3, 0.02))[0] < 1e-6


def test_a_latched_signal_channel_marks_heralding_and_delay_misfit():
    from dadl.calibration import calibrate_delay, calibrate_heralding

    a, b = simulate_pairs(5e4, SNSPD_A, SNSPD_B, 0.5, np.random.default_rng(2))
    b = latch(b, 0.4, 0.1)
    her = calibrate_heralding(a, b, CHS, width_ps=200, delay_signal_ps=500, n_draws=10)
    assert her.status == "misfit" and any("gap" in n for n in her.notes)
    assert calibrate_delay(a, b, CHS, n_draws=10).status == "misfit"


def test_a_latch_in_the_hom_dip_run_is_labelled():
    from dadl.analysis.chsh import Window
    from dadl.analysis.hom import hom_from_tags
    from dadl.guardrails import INTERRUPTED

    rng = np.random.default_rng(3)
    runs = {}
    for name in ("dip", "reference"):
        a, b = simulate_pairs(5e4, SNSPD_A, SNSPD_B, 0.5, rng)
        runs[name] = {"A": a, "B": b}
    window = Window(width_ps=200, delays_ps={"B": 500}, source="plan")
    clean, _ = hom_from_tags(runs, window)
    assert INTERRUPTED not in clean.labels
    runs["dip"]["A"] = latch(runs["dip"]["A"], 0.5, 0.2)
    latched, _ = hom_from_tags(runs, window)
    assert INTERRUPTED in latched.labels


def test_heralding_closed_form_uses_the_conditional_live_fraction():
    # Statistician review: the unconditional live fraction read 1.3% low at 1 Mpairs/s, 20 ns.
    from dadl.calibration import calibrate_heralding
    from dadl.closed_form import heralding_ratio_expected

    a_det = DetectorParams(efficiency=0.8, dark_rate_cps=200.0, dead_time_ps=20_000, jitter_sigma_ps=50.0, delay_ps=0)
    b_det = DetectorParams(efficiency=0.8, dark_rate_cps=200.0, dead_time_ps=20_000, jitter_sigma_ps=50.0, delay_ps=1000)
    truth = heralding_ratio_expected(1e6, 0.8, 0.8, 200, 200, 2e-8, 2e-8, math.hypot(50, 50), 500)
    values = []
    for seed in range(4):
        a, b = simulate_pairs(1e6, a_det, b_det, 0.05, np.random.default_rng(seed))
        cal = calibrate_heralding(a, b, CHS, width_ps=500, delay_signal_ps=1000, dead_time_ps=20_000, n_draws=10)
        values.append(cal.estimates["heralding_efficiency"].value)
    assert abs(np.mean(values) / truth - 1) < 0.003
    assert any("live fraction given a recorded herald" in n for n in cal.notes)


def test_dcr_refuses_current_biased_detectors():
    from dadl.calibration import calibrate_dark, calibrate_dcr_temperature
    from dadl.sim import simulate_dark

    darks = []
    for i, t in enumerate((1.0, 1.5, 2.0)):
        det = DetectorParams(dark_rate_cps=100.0 * (i + 1), dead_time_ps=30_000)
        tags = simulate_dark(det, 2.0, np.random.default_rng(i))
        darks.append(calibrate_dark(tags, CHS[0], seed=0, conditions={"temperature_K": t, "bias_uA": 10.0 + i}))
    with pytest.raises(ValueError, match="current"):
        calibrate_dcr_temperature(darks)


def test_dark_only_dead_time_is_an_upper_limit_and_unidentified():
    from dadl.calibration import calibrate_dark
    from dadl.sim import simulate_dark

    tags = simulate_dark(SNSPD_A, 30.0, np.random.default_rng(4))
    cal = calibrate_dark(tags, CHS[0], seed=0)
    assert cal.estimates["dead_time_ps"].low == 0 and "dead_time_ps" in cal.unidentified
    assert "upper limit only" in str(cal)


def test_calibrate_dark_refuses_source_on_data(tmp_path, capsys):
    tags = tmp_path / "tags.npz"
    assert main(["simulate", "chsh-spad", "--out", str(tags)]) == 0
    capsys.readouterr()
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(tags), "--out", str(tmp_path / "a")]) == 1
    assert "source was on" in capsys.readouterr().err


def test_sim_backend_with_unmarked_tags_is_refused(tmp_path, capsys):
    from dadl.tagfile import read_tags, write_tags
    from dadl.config import load_setup
    from dadl.cli import example_path

    sim = tmp_path / "dark.npz"
    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "2", "--out", str(sim)]) == 0
    setup = load_setup(example_path("chsh-spad"))
    plain = write_tags(tmp_path / "recorded.npz", read_tags(sim, setup.channels, setup.setting_keys), setup.channels, simulated=False)
    capsys.readouterr()
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(plain), "--channel", "A+", "--out", str(tmp_path / "a")]) == 1
    assert "backend: file" in capsys.readouterr().err


def test_calibrate_pairs_refuses_each_artifact_on_its_own(tmp_path, capsys, monkeypatch):
    # SNSPD specialist review: a timing tail made the delay fit a misfit and stopped the command
    # before the accidentals and heralding artifacts were written.
    import dadl.calibration as cal_mod
    from dadl.calibration import delay as delay_mod

    real = delay_mod.calibrate_delay

    def misfit(*args, **kwargs):
        import dataclasses

        return dataclasses.replace(real(*args, **kwargs), ppc_pvalue=0.0)

    monkeypatch.setattr(cal_mod, "calibrate_delay", misfit)
    her = tmp_path / "her.npz"
    assert main(["simulate", "chsh-spad", "--heralding", "--duration-s", "1", "--out", str(her)]) == 0
    out = tmp_path / "art"
    argv = ["calibrate", "pairs", "chsh-spad", "--tags", str(her), "--herald", "A+", "--signal", "B+",
            "--no-analysers", "--out", str(out)]
    assert main(argv) == 1
    text = capsys.readouterr()
    assert "Window capture" in text.out and "1 of 3 artifacts not written" in text.err
    assert len(list(out.glob("*[0-9a-f].json"))) == 2


def test_picoharp_router_channel_four_is_read(tmp_path):
    from dadl.adapters import ptu, read_ptu
    from ptu_writer import write_ptu

    path = write_ptu(tmp_path / "f.ptu", ptu.PICOHARP_T2, [(10, "tag", 4), (20, "tag", 0)])
    assert len(read_ptu(path).tags) == 2


def test_unsupported_config_still_raises(tmp_path):
    from dadl.config import load_setup

    bad = tmp_path / "s.yaml"
    bad.write_text("schema_version: 1\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_setup(bad)


# --- Cycle 2 -----------------------------------------------------------------------------------


def test_bell_visibility_always_carries_its_caveats_and_checks_the_window():
    from dadl.analysis.chsh import FAIR_SAMPLING, Window
    from dadl.analysis.visibility import bell_visibility
    from dadl.errors import GuardrailRefusal
    from dadl.guardrails import SINGLE_BASIS, WINDOW_UNDECLARED

    r = bell_visibility(9000, 1000)
    assert {FAIR_SAMPLING, SINGLE_BASIS, WINDOW_UNDECLARED} <= set(r.labels)
    planned = bell_visibility(9000, 1000, window=Window(width_ps=2000, source="plan"))
    assert WINDOW_UNDECLARED not in planned.labels and planned.window["source"] == "plan"
    with pytest.raises(GuardrailRefusal):
        bell_visibility(9000, 1000, window=Window(width_ps=2000, source="data"))
    with pytest.raises(ValueError, match="whole numbers"):
        bell_visibility(90.5, 10.0)


def test_chsh_expectation_follows_the_declared_angles(tmp_path, capsys):
    from dadl.cli import example_path

    text = example_path("chsh-spad").read_text(encoding="utf-8")
    assert "settings_deg" in text
    alt = tmp_path / "alt.yaml"
    import re

    alt.write_text(re.sub(r"settings_deg:.*", 'settings_deg: {a: 0, "a\'": 45, b: -22.5, "b\'": 22.5}', text, count=1),
                   encoding="utf-8")
    assert main(["quickstart", str(alt)]) == 0
    out = capsys.readouterr().out
    expected = float(out.split("Closed-form expectation")[1].split("S = ")[1].split()[0])
    raw = float(out.split("Raw S = ")[1].split()[0])
    assert abs(expected) < 0.05 and abs(raw) < 0.1 and "warning" in out


def test_negative_s_beyond_minus_two_is_labelled():
    from dadl.guardrails import NEAR_THRESHOLD, threshold_labels
    from dadl.result import Estimate

    assert NEAR_THRESHOLD in threshold_labels("chsh", Estimate(-2.01, 0.02))
    assert threshold_labels("chsh", Estimate(-2.01, 0.02)) == threshold_labels("chsh", Estimate(2.01, 0.02))


def test_a_common_delay_offset_is_refused():
    from dadl.analysis.chsh import Window, chsh

    runs = {k: {ch: np.arange(0, 10**9, 10**6, dtype=np.int64) for ch in ("A+", "A-", "B+", "B-")}
            for k in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    window = Window(width_ps=2000, delays_ps={"A+": 300, "A-": 300, "B+": 500, "B-": 500}, source="plan")
    with pytest.raises(ValueError, match="reference channel"):
        chsh(runs, window)


def test_coarse_tick_delay_is_unbiased_with_the_triangle_kernel():
    from dadl.calibration import calibrate_delay

    rng = np.random.default_rng(7)
    covered = 0
    for _ in range(20):
        t = np.sort(rng.uniform(0, 5e11, 20_000))
        a = np.floor((t + rng.normal(0, 15, t.size)) / 80).astype(np.int64) * 80
        b = np.floor((t + 1220 + rng.normal(0, 15, t.size)) / 80).astype(np.int64) * 80
        cal = calibrate_delay(np.sort(a), np.sort(b), CHS, tick_ps=80, n_draws=10)
        covered += cal.estimates["delay_ps"].contains(1220.0)
    assert covered >= 16


def test_artifacts_round_trip_large_numbers_and_refuse_non_canonical_files(tmp_path):
    import datetime as dt
    import json

    from dadl.artifact import read_artifact, write_artifact
    from dadl.calibration import calibrate_dark
    from dadl.artifact import ArtifactError
    from dadl.sim import simulate_dark

    tags = simulate_dark(DetectorParams(dark_rate_cps=1000.0, dead_time_ps=22_000), 2.0, np.random.default_rng(1))
    cal = calibrate_dark(tags, CHS[0], conditions={"start_ps": 1.7e18})
    created = dt.datetime(2026, 11, 12, tzinfo=dt.timezone.utc)
    path = write_artifact(cal, tmp_path / "a", created=created)
    manifest, _ = read_artifact(path)
    assert manifest["conditions"]["start_ps"] == 1.7e18
    # Same data, no seed: the same artifact (the draw seed comes from the input hashes).
    again = write_artifact(calibrate_dark(tags, CHS[0], conditions={"start_ps": 1.7e18}), tmp_path / "b", created=created)
    assert again.name == path.name
    pretty = tmp_path / "c" / path.name
    pretty.parent.mkdir()
    pretty.write_text(json.dumps(json.loads(path.read_text()), indent=1), encoding="utf-8")
    (pretty.parent / path.name.replace(".json", ".posterior.npy")).write_bytes(
        path.with_name(path.name.replace(".json", ".posterior.npy")).read_bytes()
    )
    with pytest.raises(ArtifactError, match="canonical"):
        read_artifact(pretty)
    dup = path.read_bytes().replace(b'{"backend"', b'{"kind":"dcr-temperature","backend"', 1)
    pretty.write_bytes(dup)
    with pytest.raises(ArtifactError):
        read_artifact(pretty)


def test_bools_are_not_integers_and_csv_bom_is_accepted(tmp_path):
    from dadl import _checks
    from dadl.tagfile import read_tags

    with pytest.raises(TypeError):
        _checks.integer(True, "x")
    p = tmp_path / "t.csv"
    p.write_bytes("﻿setting,channel,timestamp_ps\n0,0,100\n".encode("utf-8"))
    assert read_tags(p)[("a", "b")]["A+"].tolist() == [100]
    p.write_text("setting,channel,timestamp_ps\n0,0,9223372036854775808\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        read_tags(p)


def test_conditional_laplace_draws_do_not_need_a_definite_covariance():
    # Found by the gates after cycle 2: a fit with an unidentified background crashed when drawing.
    from dadl.stats import Laplace

    cov = np.zeros((3, 3))
    cov[1:, 1:] = [[1.0, 0.2], [0.2, 2.0]]
    d = Laplace(np.zeros(3), cov).draws(np.random.default_rng(0), 500)
    assert d.shape == (500, 3) and np.all(d[:, 0] == 0)


# --- Cycle 3 -----------------------------------------------------------------------------------


def test_artifacts_from_different_data_under_one_seed_have_independent_draws():
    from dadl.calibration import calibrate_dark
    from dadl.sim import simulate_dark

    det = DetectorParams(dark_rate_cps=2000.0, dead_time_ps=22_000)
    a = calibrate_dark(simulate_dark(det, 5.0, np.random.default_rng(1)), CHS[0], seed=2026)
    b = calibrate_dark(simulate_dark(det, 5.0, np.random.default_rng(2)), CHS[1], seed=2026)
    r = np.corrcoef(a.samples[:, 0], b.samples[:, 0])[0, 1]
    assert abs(r) < 0.1


def test_jitter_far_below_the_tick_is_an_upper_limit_and_the_delay_still_covers():
    from dadl.calibration import calibrate_delay

    rng = np.random.default_rng(11)
    covered = 0
    for _ in range(15):
        t = np.sort(rng.uniform(0, 1e12, 20_000))
        off = int(rng.integers(0, 80))
        a = (np.floor((t + off + rng.normal(0, 3, t.size)) / 80) * 80).astype(np.int64)
        b = (np.floor((t + off + 1013 + rng.normal(0, 3, t.size)) / 80) * 80).astype(np.int64)
        cal = calibrate_delay(np.sort(a), np.sort(b), CHS, bin_width_ps=80, tick_ps=80, n_draws=50)
        assert "jitter_sigma_ps" in cal.unidentified and cal.estimates["jitter_sigma_ps"].low == 0
        assert "jitter_sigma_ps" in cal.stand_in
        covered += cal.estimates["delay_ps"].contains(1013.0)
    assert covered >= 12


def test_hom_compares_rates_when_the_acquisitions_differ_in_length():
    from dadl.analysis.chsh import Window
    from dadl.analysis.hom import hom_from_tags
    from dadl.sim import PairSource, simulate_hom

    src = PairSource(pair_rate_per_s=2e5, visibility=0.92, coherence_time_ps=1.5)
    rng = np.random.default_rng(5)
    dip = simulate_hom(src, SNSPD_A, SNSPD_A, 0.0, 0.5, rng)
    ref = simulate_hom(src, SNSPD_A, SNSPD_A, 15.0, 1.0, rng)
    result, counts = hom_from_tags({"dip": dip, "reference": ref}, Window(width_ps=200, source="plan"))
    assert abs(counts["duration_dip_s"] / counts["duration_reference_s"] - 0.5) < 0.01
    assert result.raw.low < 0.925 and result.raw.high > 0.91


def test_hom_from_tags_keeps_the_not_computed_label():
    from dadl.analysis.chsh import Window
    from dadl.analysis.hom import hom_from_tags
    from dadl.guardrails import CORRECTION_NOT_COMPUTED

    rng = np.random.default_rng(6)
    span = 10**12
    background = lambda n: np.sort(rng.integers(0, span, n)).astype(np.int64)  # noqa: E731
    runs = {k: {"A": background(400_000), "B": background(400_000)} for k in ("dip", "reference")}
    result, counts = hom_from_tags(runs, Window(width_ps=2000, source="plan"))
    assert counts["subtracted"] == 0.0 and CORRECTION_NOT_COMPUTED in result.labels


def test_limit_only_parameters_are_null_with_a_bound_in_the_record(tmp_path):
    import datetime as dt

    from dadl.artifact import build_manifest
    from dadl.calibration import calibrate_dark
    from dadl.sim import simulate_dark

    tags = simulate_dark(DetectorParams(dark_rate_cps=500.0, dead_time_ps=22_000), 2.0, np.random.default_rng(3))
    cal = calibrate_dark(tags, CHS[0])
    manifest, _ = build_manifest(cal, created=dt.datetime(2026, 11, 12, tzinfo=dt.timezone.utc))
    p = manifest["summary"]["afterpulse_prob"]
    assert p["value"] is None and p["bound"] == "upper" and p["unit"] == "1" and "definition" in p
    assert "afterpulse_prob" in manifest["posterior"]["stand_in"]
    assert "prior_dominated_params" not in manifest["checks"]
    assert {"numpy", "scipy", "python"} <= set(manifest["software"])


def test_hardware_calibrations_need_the_detector_serial_and_record_the_run(tmp_path, capsys):
    from dadl.artifact import read_artifact

    tags = tmp_path / "dark.npz"
    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "2", "--out", str(tags)]) == 0
    from dadl.cli import example_path

    base = example_path("chsh-spad").read_text(encoding="utf-8").replace("backend: sim", "backend: file")
    if "backend: file" not in base:
        base = base.replace("seed:", "backend: file\nseed:", 1)
    hw = tmp_path / "hw.yaml"
    hw.write_text(base + "hardware:\n  tagger: HydraHarp-1\n", encoding="utf-8")
    capsys.readouterr()
    argv = ["calibrate", "dark", str(hw), "--tags", str(tags), "--channel", "A+", "--out", str(tmp_path / "a")]
    assert main(argv) == 1 and "detector's serial" in capsys.readouterr().err
    hw.write_text(base + "hardware:\n  tagger: HydraHarp-1\n  detectors: {A+: SPCM-7}\n", encoding="utf-8")
    assert main(argv + ["--condition", "temperature_K=233.5", "--acquired", "2026-11-01T10:00:00Z"]) == 0
    (path,) = (tmp_path / "a").glob("*[0-9a-f].json")
    m = read_artifact(path)[0]
    assert m["conditions"]["temperature_K"] == 233.5 and m["conditions"]["duration_s"] > 1.5
    assert m["acquired"] == "2026-11-01T10:00:00Z"
    assert {i["role"] for i in m["inputs"]} == {"raw-tags", "raw-file", "config"}
    assert main(argv + ["--supersedes", "0" * 64]) == 1


# --- Cycle 4 -----------------------------------------------------------------------------------


@pytest.mark.parametrize("start", [0.0, 0.7])
def test_a_latch_at_the_start_or_end_of_a_run_is_caught(start):
    # QKD review: a trailing or leading latch left no internal gap and read status ok.
    from dadl.calibration import calibrate_heralding

    a, b = simulate_pairs(5e4, SNSPD_A, SNSPD_B, 0.5, np.random.default_rng(8))
    b = latch(b, start, 0.3)
    her = calibrate_heralding(a, b, CHS, width_ps=200, delay_signal_ps=500, n_draws=10)
    assert her.status == "misfit"


def test_dead_time_selection_is_labelled_and_ns_tags_are_refused(tmp_path, capsys):
    from dadl.analysis.chsh import Window, chsh
    from dadl.guardrails import DEAD_TIME_SELECTION
    from dadl.sim import PairSource, simulate_chsh

    det = DetectorParams(efficiency=0.5, dark_rate_cps=0.0, dead_time_ps=1_000_000, jitter_sigma_ps=350.0)
    runs = simulate_chsh(PairSource(pair_rate_per_s=1e6, visibility=0.7), {k: det for k in ("A+", "A-", "B+", "B-")},
                         0.05, np.random.default_rng(1))
    assert DEAD_TIME_SELECTION in chsh(runs, Window(width_ps=2000, source="plan")).labels

    from dadl.tagfile import read_tags, write_tags
    from dadl.config import load_setup
    from dadl.cli import example_path

    tags = tmp_path / "t.npz"
    assert main(["simulate", "chsh-spad", "--out", str(tags)]) == 0
    setup = load_setup(example_path("chsh-spad"))
    ns = {k: {c: v // 1000 for c, v in r.items()} for k, r in read_tags(tags, setup.channels).items()}
    path = write_tags(tmp_path / "ns.npz", ns, setup.channels, simulated=True)
    capsys.readouterr()
    assert main(["analyze", "chsh-spad", "--tags", str(path)]) == 1
    assert "picoseconds" in capsys.readouterr().err


def test_corrected_values_carry_a_label_and_tampered_records_are_refused(tmp_path):
    import datetime as dt
    import json

    from dadl.artifact import ArtifactError, artifact_id, canonical_json, read_artifact, write_artifact
    from dadl.calibration import calibrate_dark
    from dadl.guardrails import ACCIDENTALS_SUBTRACTED
    from dadl.sim import simulate_dark

    assert ACCIDENTALS_SUBTRACTED in hom_visibility_labels()
    tags = simulate_dark(DetectorParams(dark_rate_cps=1000.0, dead_time_ps=22_000), 2.0, np.random.default_rng(1))
    path = write_artifact(calibrate_dark(tags, CHS[0]), tmp_path, created=dt.datetime(2026, 11, 12, tzinfo=dt.timezone.utc))
    m = json.loads(path.read_text())
    m["checks"]["status"] = "misfit"
    raw = canonical_json(m)
    bad = tmp_path / "bad"
    bad.mkdir()
    new = bad / f"{artifact_id(m)}.json"
    new.write_bytes(raw)
    (bad / f"{artifact_id(m)}.posterior.npy").write_bytes(path.with_name(path.stem + ".posterior.npy").read_bytes())
    with pytest.raises(ArtifactError, match="status"):
        read_artifact(new)


def hom_visibility_labels():
    from dadl.analysis.visibility import hom_visibility

    return hom_visibility(1000.0, 120.0, accidentals=40.0).labels


def test_long_delays_fit_like_short_ones():
    # Numerical review, cycle 4: a relative Hessian step of 100 ps at a 1 us delay inflated the
    # interval 2-40x or refused the fit.
    from dadl.calibration import calibrate_delay

    widths = []
    for d in (1_000, 2_000_000):
        A = DetectorParams(efficiency=0.5, dark_rate_cps=200.0, jitter_sigma_ps=3.5, delay_ps=0)
        B = DetectorParams(efficiency=0.5, dark_rate_cps=200.0, jitter_sigma_ps=3.5, delay_ps=d)
        a, b = simulate_pairs(2e4, A, B, 0.5, np.random.default_rng(3))
        e = calibrate_delay(a, b, CHS, range_ps=d + 20_000, n_draws=10).estimates["delay_ps"]
        assert e.contains(d)
        widths.append(e.high - e.low)
    assert widths[1] == pytest.approx(widths[0], rel=0.05)


# --- Cycle 5 -----------------------------------------------------------------------------------


def test_duplicate_tags_are_refused_and_do_not_hide_dead_time_selection(tmp_path, capsys):
    from dadl.analysis.chsh import dead_time_selection

    runs = {("a", "b"): {ch: np.arange(0, 10**10, 10**6, dtype=np.int64) for ch in ("A+", "A-", "B+", "B-")}}
    with_dup = {k: {c: np.sort(np.append(v, v[0])) for c, v in r.items()} for k, r in runs.items()}
    counts = {("a", "b"): 1_000_000}
    assert dead_time_selection(with_dup, counts) == pytest.approx(dead_time_selection(runs, counts))
    tags = tmp_path / "t.csv"
    assert main(["simulate", "chsh-spad", "--out", str(tags)]) == 0
    lines = tags.read_text().splitlines()
    tags.write_text("\n".join(lines + [lines[2]]) + "\n")
    capsys.readouterr()
    assert main(["analyze", "chsh-spad", "--tags", str(tags)]) == 1
    assert "duplicate" in capsys.readouterr().err


def test_simulated_csv_stays_marked_and_analyze_labels_it(tmp_path, capsys):
    from dadl.guardrails import SIMULATED_INPUT
    from dadl.tagfile import tags_are_simulated

    tags = tmp_path / "t.csv"
    assert main(["simulate", "chsh-spad", "--out", str(tags)]) == 0
    assert tags_are_simulated(tags)
    capsys.readouterr()
    assert main(["analyze", "chsh-spad", "--tags", str(tags)]) == 0
    out = capsys.readouterr().out
    assert SIMULATED_INPUT in out and "Inputs:" in out


def test_policy_rule_subclasses_are_refused():
    from dadl.guardrails import DEFAULT_RULES, AnalysisRule, Policy

    class Lax(AnalysisRule):
        def is_at_least_as_strict_as(self, other):
            return True

    rules = dict(DEFAULT_RULES)
    rules["bell_visibility"] = Lax(frozenset({"accidental_subtraction"}), frozenset({"plan"}))
    with pytest.raises(TypeError):
        Policy(rules=rules)


def test_mismatched_dead_times_still_get_the_selection_label():
    from dadl.analysis.chsh import Window, chsh
    from dadl.guardrails import DEAD_TIME_SELECTION
    from dadl.sim import PairSource, simulate_chsh

    slow = DetectorParams(efficiency=0.5, dark_rate_cps=0.0, dead_time_ps=1_000_000, jitter_sigma_ps=350.0)
    fast = DetectorParams(efficiency=0.5, dark_rate_cps=0.0, dead_time_ps=20_000, jitter_sigma_ps=350.0)
    runs = simulate_chsh(PairSource(pair_rate_per_s=1e6, visibility=0.7),
                         {"A+": slow, "A-": fast, "B+": slow, "B-": fast}, 0.05, np.random.default_rng(2))
    assert DEAD_TIME_SELECTION in chsh(runs, Window(width_ps=2000, source="plan")).labels


def test_interruption_test_holds_its_level_at_low_counts():
    from dadl.stationarity import interruption_pvalue

    rng = np.random.default_rng(0)
    ps = np.array([interruption_pvalue(np.cumsum(rng.exponential(1e6, 22)).astype(np.int64))[0] for _ in range(5000)])
    assert np.mean(ps < 0.01) < 0.02


def test_percent_level_dark_drift_is_tolerated_and_recorded():
    from dadl.calibration.dark import _stationarity

    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(0, 1e12, 1_000_000))
    keep = rng.random(t.size) < 1.0 - 0.01 * (t / 1e12)  # 1% drift end to end
    tags = t[keep].astype(np.int64)
    st = _stationarity(tags, tags.size / 1e12, tags.size - 1, float(np.diff(tags).max()))
    assert st["p_drift"] > 0.01 and 0.98 < st["drift_ratio"] < 1.0
    keep = rng.random(t.size) < 1.0 - 0.2 * (t / 1e12)  # 20% drift
    tags = t[keep].astype(np.int64)
    assert _stationarity(tags, tags.size / 1e12, tags.size - 1, float(np.diff(tags).max()))["p_drift"] < 1e-6


# --- Cycle 6 -----------------------------------------------------------------------------------


def test_free_running_ingaas_afterpulsing_is_fitted():
    from dadl.calibration import calibrate_dark
    from dadl.sim import simulate_dark

    det = DetectorParams(dark_rate_cps=1000.0, dead_time_ps=10_000_000, afterpulse_prob=0.05,
                         afterpulse_tau_ps=20_000_000.0)
    cal = calibrate_dark(simulate_dark(det, 30.0, np.random.default_rng(4)), CHS[0], n_draws=50)
    assert cal.status == "ok"
    assert cal.estimates["afterpulse_prob"].contains(0.05) and cal.estimates["afterpulse_tau_ps"].contains(2e7)


def test_dcr_refuses_misfit_inputs_unless_accepted():
    import dataclasses

    from dadl.calibration import calibrate_dark, calibrate_dcr_temperature
    from dadl.sim import simulate_dark

    darks = []
    for i, t in enumerate((253.0, 263.0, 273.0)):
        det = DetectorParams(dark_rate_cps=100.0 * 2**i, dead_time_ps=22_000)
        darks.append(calibrate_dark(simulate_dark(det, 5.0, np.random.default_rng(i)), CHS[0],
                                    conditions={"temperature_K": t}))
    darks[0] = dataclasses.replace(darks[0], ppc_pvalue=0.0)
    with pytest.raises(ValueError, match="misfit"):
        calibrate_dcr_temperature(darks)
    cal = calibrate_dcr_temperature(darks, accept_misfit_inputs=True)
    assert any("misfit dark calibrations were accepted" in n for n in cal.notes)


def test_oversized_histograms_are_refused_not_allocated():
    from dadl.histogram import cross_correlation, interarrival

    t = np.arange(10, dtype=np.int64)
    with pytest.raises(ValueError, match="bins"):
        interarrival(t, 10, 2**40)
    with pytest.raises(ValueError, match="bins"):
        cross_correlation(t, t, 1, 2**40)


def test_records_keep_exact_intervals_and_the_grid_origin():
    from dadl.analysis.chsh import Window
    from dadl.analysis.visibility import bell_visibility

    rec = bell_visibility(30, 3, window=Window(width_ps=2000, source="plan")).as_record()
    assert rec["raw_low"] is not None and rec["raw_low"] < rec["raw_value"] < rec["raw_high"]
    assert rec["window_grid_origin_ps"] == 0


def test_short_runs_do_not_raise_the_dead_time_label():
    from dadl.analysis.chsh import Window, chsh
    from dadl.guardrails import DEAD_TIME_SELECTION
    from dadl.sim import PairSource, simulate_chsh

    det = DetectorParams(efficiency=0.9, dark_rate_cps=500.0, dead_time_ps=22_000, jitter_sigma_ps=350.0)
    labelled = 0
    for s in range(10):
        runs = simulate_chsh(PairSource(pair_rate_per_s=2e4, visibility=0.7),
                             {k: det for k in ("A+", "A-", "B+", "B-")}, 0.2, np.random.default_rng(s))
        labelled += DEAD_TIME_SELECTION in chsh(runs, Window(width_ps=2000, source="plan")).labels
    assert labelled == 0


def test_library_refuses_duplicate_timestamps():
    from dadl.analysis.chsh import Window, chsh

    t = np.arange(0, 10**9, 10**5, dtype=np.int64)
    runs = {k: {ch: t for ch in ("A+", "A-", "B+", "B-")} for k in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    runs[("a", "b")]["A+"] = np.sort(np.append(t, t[5]))
    with pytest.raises(ValueError, match="duplicate"):
        chsh(runs, Window(width_ps=2000, source="plan"))


def test_ci_coverage_matrix_runs_every_regime():
    import re
    from pathlib import Path

    from dadl.calibration.coverage import UNITS

    ci = Path(__file__).resolve().parents[2] / ".github" / "workflows" / "ci.yml"
    if not ci.exists():
        pytest.skip("repository CI file not present")
    args = re.findall(r'args: "([^"]+)"', ci.read_text(encoding="utf-8"))
    whole = {a.split()[0] for a in args if "--regimes" not in a}
    listed = {r for a in args if "--regimes" in a for r in a.split("--regimes")[1].split()}
    for family, regimes in UNITS.items():
        assert family in whole or set(regimes) <= listed, family


# --- Cycle 7 -----------------------------------------------------------------------------------


def test_chsh_guards_do_not_depend_on_tag_order():
    from dadl.analysis.chsh import Window, chsh
    from dadl.guardrails import DEAD_TIME_SELECTION
    from dadl.sim import PairSource, simulate_chsh

    det = DetectorParams(efficiency=0.5, dark_rate_cps=0.0, dead_time_ps=1_000_000, jitter_sigma_ps=350.0)
    rng = np.random.default_rng(3)
    runs = simulate_chsh(PairSource(pair_rate_per_s=1e6, visibility=0.7), {k: det for k in ("A+", "A-", "B+", "B-")},
                         0.05, rng)
    shuffled = {k: {c: rng.permutation(v) for c, v in r.items()} for k, r in runs.items()}
    w = Window(width_ps=2000, source="plan")
    assert DEAD_TIME_SELECTION in chsh(shuffled, w).labels
    assert chsh(shuffled, w).raw.value == chsh(runs, w).raw.value


def test_windows_need_a_declared_source_and_known_delay_channels(tmp_path):
    from dadl.analysis.chsh import Window, chsh

    with pytest.raises(TypeError):
        Window(width_ps=2000)  # type: ignore[call-arg]
    w = Window(width_ps=2000, source="plan", delays_ps={"Bminus": 3000})
    runs = {k: {ch: np.arange(0, 10**9, 10**5, dtype=np.int64) for ch in ("A+", "A-", "B+", "B-")}
            for k in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    with pytest.raises(ValueError, match="unknown channel"):
        chsh(runs, w)
    from dadl.cli import example_path
    from dadl.config import load_setup

    text = example_path("chsh-spad").read_text(encoding="utf-8").replace("delays_ps: {}", "delays_ps: {Bminus: 3000}")
    bad = tmp_path / "s.yaml"
    bad.write_text(text, encoding="utf-8")
    with pytest.raises(ConfigError, match="unknown channel"):
        load_setup(bad)


def test_python_artifact_path_matches_the_cli(tmp_path):
    import dataclasses
    import datetime as dt

    from dadl.artifact import ArtifactError, load_calibration, write_artifact
    from dadl.calibration import calibrate_dark, calibrate_dcr_temperature
    from dadl.sim import simulate_dark

    when = dt.datetime(2026, 11, 12, tzinfo=dt.timezone.utc)
    paths = []
    for i, t in enumerate((253.0, 263.0, 273.0)):
        det = DetectorParams(dark_rate_cps=100.0 * 2**i, dead_time_ps=22_000)
        cal = calibrate_dark(simulate_dark(det, 5.0, np.random.default_rng(i)), CHS[0], conditions={"temperature_K": t})
        assert cal.backend == "sim"  # device "sim" makes the default backend "sim"
        paths.append(write_artifact(cal, tmp_path, created=when))
    loaded = [load_calibration(p) for p in paths]
    assert calibrate_dcr_temperature(loaded).kind == "dcr-temperature"
    with pytest.raises(ArtifactError, match="misfit"):
        write_artifact(dataclasses.replace(loaded[0], ppc_pvalue=0.0), tmp_path / "x", created=when)


def test_estimates_format_like_numbers():
    from dadl.result import Estimate

    assert f"{Estimate(2.69412, 0.01049):.3f}" == "2.694 ± 0.010"


# --- Cycle 8 -----------------------------------------------------------------------------------


def test_numpy_conditions_and_named_inputs_are_written(tmp_path):
    import datetime as dt

    from dadl.artifact import read_artifact, write_artifact
    from dadl.calibration import calibrate_dark
    from dadl.sim import simulate_dark

    tags = simulate_dark(DetectorParams(dark_rate_cps=1000.0, dead_time_ps=22_000), 2.0, np.random.default_rng(1))
    cal = calibrate_dark(tags, CHS[0], conditions={"temperature_K": np.mean([249.0, 251.0])})
    m = read_artifact(write_artifact(cal, tmp_path, created=dt.datetime(2026, 11, 12, tzinfo=dt.timezone.utc)))[0]
    assert m["conditions"]["temperature_K"] == 250.0
    assert m["inputs"][0]["name"] == "sim:0"
    with pytest.raises(ValueError, match="gaps"):
        calibrate_dark(tags, CHS[0], conditions={"gaps": 1.0})


def test_clock_resets_and_broken_npz_are_refused(tmp_path):
    from dadl.tagfile import read_tags

    p = tmp_path / "t.csv"
    p.write_text("setting,channel,timestamp_ps\n0,0,100\n0,0,200\n0,0,50\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="back in time"):
        read_tags(p)
    z = tmp_path / "t.npz"
    z.write_bytes(b"PK\x03\x04broken")
    with pytest.raises(ConfigError, match="npz"):
        read_tags(z)


def test_stale_declared_delay_marks_heralding_misfit_in_the_library():
    from dadl.calibration import calibrate_heralding

    a_det = DetectorParams(efficiency=0.5, dark_rate_cps=500.0, jitter_sigma_ps=350.0, delay_ps=0)
    b_det = DetectorParams(efficiency=0.4, dark_rate_cps=500.0, jitter_sigma_ps=350.0, delay_ps=3000)
    a, b = simulate_pairs(2e4, a_det, b_det, 1.0, np.random.default_rng(0))
    assert calibrate_heralding(a, b, CHS, width_ps=2000, delay_signal_ps=3000, n_draws=10).status == "ok"
    stale = calibrate_heralding(a, b, CHS, width_ps=2000, delay_signal_ps=0, n_draws=10)
    assert stale.status == "misfit" and any("pair peak" in n for n in stale.notes)


def test_acquired_needs_an_offset_and_records_are_json(tmp_path, capsys):
    import json

    from dadl.analysis.visibility import bell_visibility

    rec = bell_visibility(9000, 1000).as_record()
    labels = json.loads(rec["labels"])
    assert len(labels) == 3 and rec["simulated"] is False and rec["dadl_version"]
    tags = tmp_path / "dark.npz"
    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "2", "--out", str(tags)]) == 0
    capsys.readouterr()
    argv = ["calibrate", "dark", "chsh-spad", "--tags", str(tags), "--channel", "A+", "--out", str(tmp_path / "a")]
    assert main(argv + ["--acquired", "2026-11-01T09:00"]) == 1
    assert "UTC offset" in capsys.readouterr().err
