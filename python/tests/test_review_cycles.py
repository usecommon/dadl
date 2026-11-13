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
    plain = write_tags(tmp_path / "recorded.npz", read_tags(sim, setup.channels, setup.setting_keys), setup.channels)
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
