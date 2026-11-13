"""HOM: simulator, closed form, setup files and the hom-snspd example (milestone 3 of phase 2)."""

from pathlib import Path

import numpy as np
import pytest

from dadl.analysis.chsh import Window
from dadl.analysis.hom import hom_from_tags
from dadl.cli import example_path, main
from dadl.closed_form import hom_split_probability, hom_visibility_expected
from dadl.coincidences import grid_coincidences
from dadl.config import load_setup, parse_setup, read_yaml
from dadl.errors import ConfigError
from dadl.sim import DetectorParams, PairSource, simulate_hom

ROOT_EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "hom-snspd" / "setup.yaml"
SRC = PairSource(pair_rate_per_s=2e5, visibility=0.9, coherence_time_ps=2.0)
DET = DetectorParams(efficiency=0.8, dark_rate_cps=0.0, dead_time_ps=30_000, jitter_sigma_ps=10)


def test_closed_form():
    assert hom_split_probability(1.0, 0.0, 1.0) == 0.0
    assert hom_split_probability(0.9, 1e3, 1.0) == pytest.approx(0.5)
    assert hom_visibility_expected(0.9, 20.0, 2.0) == pytest.approx(0.9)


def test_simulated_dip_matches_the_closed_form():
    rng = np.random.default_rng(1)
    c = {}
    for delay in (0.0, 2.0, 20.0):
        r = simulate_hom(SRC, DET, DET, delay, 2.0, rng)
        c[delay] = grid_coincidences(r["A"], r["B"], 500)
    for delay in (0.0, 2.0):
        ratio = hom_split_probability(0.9, delay, 2.0) / hom_split_probability(0.9, 20.0, 2.0)
        assert c[delay] / c[20.0] == pytest.approx(ratio, abs=4 * np.sqrt(c[delay]) / c[20.0])


def test_bunched_photons_are_counted_once():
    # Perfect interference: every pair leaves by one port, so singles fall to about
    # eta(2 - eta) / 2 per pair per detector and only accidentals (about 1 here) remain.
    src = PairSource(pair_rate_per_s=1e5, visibility=1.0, coherence_time_ps=2.0)
    r = simulate_hom(src, DET, DET, 0.0, 1.0, np.random.default_rng(2))
    assert grid_coincidences(r["A"], r["B"], 500) <= 8
    eta = DET.efficiency
    assert r["A"].size == pytest.approx(1e5 * eta * (2 - eta) / 2, rel=0.02)


def test_hom_from_tags_reports_raw_first_and_subtracts_off_peak_accidentals():
    noisy = DetectorParams(efficiency=0.8, dark_rate_cps=2e5, dead_time_ps=30_000, jitter_sigma_ps=10)
    rng = np.random.default_rng(3)
    runs = {"dip": simulate_hom(SRC, noisy, noisy, 0.0, 1.0, rng), "reference": simulate_hom(SRC, noisy, noisy, 40.0, 1.0, rng)}
    result, counts = hom_from_tags(runs, Window(width_ps=1000, source="plan"))
    assert counts["accidentals"] > 0
    assert result.corrected.value > result.raw.value
    assert result.corrected.value == pytest.approx(0.9, abs=4 * result.corrected.stderr)
    assert result.window["width_ps"] == 1000
    with pytest.raises(ValueError, match="reference"):
        hom_from_tags({"dip": runs["dip"]}, Window(width_ps=1000, source="plan"))


def test_example_loads():
    setup = load_setup(example_path("hom-snspd"))
    assert setup.analysis_type == "hom" and setup.channels == {0: "A", 1: "B"}
    assert setup.setting_keys == ("dip", "reference") and setup.reference_delay_ps == 15


@pytest.mark.skipif(not ROOT_EXAMPLE.exists(), reason="repository examples/ not present")
def test_repository_example_matches_bundled_copy():
    assert ROOT_EXAMPLE.read_bytes() == example_path("hom-snspd").read_bytes()


@pytest.mark.parametrize(
    "change, match",
    [
        (lambda d: d["analysis"].pop("reference_delay_ps"), "reference_delay_ps"),
        (lambda d: d["analysis"].update(reference_delay_ps=0), "away from the dip"),
        (lambda d: d["analysis"].update(settings_deg={"a": 0}), "chsh setups"),
        (lambda d: d["detectors"].pop("B"), "missing"),
        (lambda d: d.update(channels={0: "A", 1: "A+"}), "names must be one of"),
        (lambda d: d["source"].update(coherence_time_ps=0), "coherence_time_ps"),
    ],
)
def test_bad_hom_setups_are_refused(change, match):
    data = read_yaml(example_path("hom-snspd"))
    change(data)
    with pytest.raises(ConfigError, match=match):
        parse_setup(data)


def test_quickstart_simulate_analyze(tmp_path, capsys):
    assert main(["quickstart", "hom-snspd"]) == 0
    out = capsys.readouterr().out
    assert "Raw V = 0.9" in out and "Closed-form expectation" in out and "Accidentals subtracted" in out
    tags = tmp_path / "hom.npz"
    assert main(["simulate", "hom-snspd", "--out", str(tags)]) == 0
    assert main(["analyze", "hom-snspd", "--tags", str(tags)]) == 0
    out = capsys.readouterr().out
    assert "at the dip" in out and "Raw V = 0.9" in out
    assert main(["init", str(tmp_path / "x"), "--example", "hom-snspd"]) == 0


def test_accidentals_above_c_min_are_not_subtracted_and_it_is_said(capsys):
    # Independent dense streams with every on-peak coincidence removed: C_min = 0 while the
    # off-peak windows still see accidentals, so subtraction would give V > 1 and is skipped.
    rng = np.random.default_rng(0)
    a = np.sort(rng.choice(10**10, 200_000, replace=False))
    b = np.sort(rng.choice(10**10, 200_000, replace=False))
    b = b[~np.isin(b // 200, a // 200)]
    runs = {"dip": {"A": a, "B": b}, "reference": {"A": a, "B": np.sort(np.concatenate([b, a + 50]))}}
    result, counts = hom_from_tags(runs, Window(width_ps=200, source="plan"))
    assert counts["c_min"] == 0 and counts["accidentals"] > 0
    assert result.corrected is None and counts["subtracted"] == 0.0
