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
    assert "Raw V = 0.9" in out and "Closed-form expectation" in out
    # Phase 3: the corrected line comes from the forward model with simulated dark calibrations;
    # 100 cps of dark counts only bound an SNSPD's dead time, and the quickstart says so.
    assert "Corrected (forward model: accidental subtraction, then window capture)" in out and "Dead time not corrected" in out
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


def test_analyze_with_calibrations(tmp_path, capsys):
    # A HOM setup whose detectors have enough dark counts to identify the dead time.
    setup = tmp_path / "setup.yaml"
    setup.write_text(example_path("hom-snspd").read_text(encoding="utf-8").replace("dark_rate_cps: 100", "dark_rate_cps: 20000"), encoding="utf-8")
    dark, tags, art = tmp_path / "dark.npz", tmp_path / "hom.npz", tmp_path / "art"
    assert main(["simulate", str(setup), "--dark", "--duration-s", "2", "--out", str(dark)]) == 0
    assert main(["calibrate", "dark", str(setup), "--tags", str(dark), "--out", str(art)]) == 0
    her = tmp_path / "pairs.npz"
    assert main(["simulate", str(setup), "--heralding", "--duration-s", "0.2", "--out", str(her)]) == 0
    assert main(["calibrate", "pairs", str(setup), "--tags", str(her), "--no-analysers", "--herald", "A", "--signal", "B", "--out", str(art)]) == 0
    assert main(["simulate", str(setup), "--out", str(tags)]) == 0
    capsys.readouterr()
    # Dead time is not in the default HOM rule (decision D10): refused without --allow.
    assert main(["analyze", str(setup), "--tags", str(tags), "--calibrations", str(art)]) == 1
    assert "dead_time" in capsys.readouterr().err
    args = ["analyze", str(setup), "--tags", str(tags), "--calibrations", str(art), "--allow", "dead_time=calibrated from dark runs"]
    assert main(args) == 1  # the delay-jitter artifact brings window capture, also an explicit override
    assert "window_capture" in capsys.readouterr().err
    assert main([*args, "--allow", "window_capture=delay and jitter from a pair run"]) == 0
    out = capsys.readouterr().out
    assert "Corrected (forward model: dead time, then accidental subtraction, then window capture): V = 0.9" in out
    assert "Calibrations: " in out and "(modified)" in out
    # --allow needs a reason.
    assert main(["analyze", str(setup), "--tags", str(tags), "--allow", "dead_time"]) == 1
    # Council round 4 (S2): a dark artifact for A only gave a KeyError traceback; now the reason.
    one = tmp_path / "one"
    assert main(["calibrate", "dark", str(setup), "--tags", str(dark), "--channel", "A", "--out", str(one)]) == 0
    capsys.readouterr()
    assert main(["analyze", str(setup), "--tags", str(tags), "--calibrations", str(one), "--allow",
                 "dead_time=calibrated from dark runs"]) == 1
    err = capsys.readouterr().err
    assert "every channel needs one" in err and "Traceback" not in err


def test_analyze_chsh_with_calibrations_needs_the_policy_and_every_pair(tmp_path, capsys):
    tags, dark, her = tmp_path / "chsh.npz", tmp_path / "dark.npz", tmp_path / "her.npz"
    assert main(["simulate", "chsh-spad", "--out", str(tags)]) == 0
    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "2", "--out", str(dark)]) == 0
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(dark), "--out", str(tmp_path / "dark-art")]) in (0, 1)
    capsys.readouterr()
    # The default CHSH rule allows no correction (Bell-type: raw only).
    assert main(["analyze", "chsh-spad", "--tags", str(tags), "--calibrations", str(tmp_path / "dark-art")]) == 1
    assert "dead_time" in capsys.readouterr().err
    # A delay-jitter artifact for one channel pair of four: refused rather than half-modelled.
    assert main(["simulate", "chsh-spad", "--heralding", "--duration-s", "0.5", "--out", str(her)]) == 0
    art = tmp_path / "pairs-art"
    assert main(["calibrate", "pairs", "chsh-spad", "--tags", str(her), "--no-analysers", "--herald", "A+", "--signal", "B+", "--out", str(art)]) == 0
    allow = ["--allow", "accidental_subtraction=off-peak bins", "--allow", "window_capture=pair runs"]
    capsys.readouterr()
    assert main(["analyze", "chsh-spad", "--tags", str(tags), "--calibrations", str(art), *allow]) == 1
    assert "every channel pair" in capsys.readouterr().err


def test_analyze_chsh_prints_the_dead_time_range(tmp_path, capsys):
    # Council round 3 (R5): with every dead time a bound the CHSH CLI printed the label "counts give
    # the value at zero" but not the range.
    tags, dark, art = tmp_path / "chsh.npz", tmp_path / "dark.npz", tmp_path / "dark-art"
    assert main(["simulate", "chsh-spad", "--out", str(tags)]) == 0
    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "2", "--out", str(dark)]) == 0
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(dark), "--out", str(art)]) in (0, 1)
    capsys.readouterr()
    allow = ["--allow", "dead_time=dark runs", "--allow", "accidental_subtraction=off-peak bins"]
    assert main(["analyze", "chsh-spad", "--tags", str(tags), "--calibrations", str(art), *allow]) == 0
    out = capsys.readouterr().out
    assert "Dead time known only as an upper limit: S = " in out and " at zero, " in out and " at the bound" in out


def _pulsed(period, n_pulses, rng, p_single=0.01, p_pair=0.0, jitter=50.0):
    """A pulsed source: per pulse, A and B each click with p_single, and a pair with p_pair."""
    k = np.arange(n_pulses)
    a, b, pair = (rng.random(n_pulses) < p for p in (p_single, p_single, p_pair))
    out = {}
    for name, hit in (("A", a | pair), ("B", b | pair)):
        t = k[hit] * period + rng.normal(0, jitter, int(hit.sum())) + 1000
        out[name] = np.unique(t.astype(np.int64))
    return out


def test_pulsed_t3_input_needs_a_declared_cw_source():
    # Council round 7 (G3): at 76 MHz the side windows sample random phases of the pulse train and
    # the accidentals read about W / period low; the corrected V was 0.962 for a true 1, unlabelled.
    from dadl.data import RunTags
    from dadl.errors import GuardrailRefusal
    from dadl.guardrails import T3_INPUT

    period = 13157.9
    rng = np.random.default_rng(2)
    plain = {"dip": _pulsed(period, 2_000_000, rng), "reference": _pulsed(period, 2_000_000, rng, p_pair=0.002)}
    meta = {"mode": "T3", "sync_period_ps": period, "tick_ps": 2}
    runs = {k: RunTags(v, metadata=meta) for k, v in plain.items()}
    w = Window(1000, "plan")
    with pytest.raises(GuardrailRefusal, match="Declare `source.type: cw`"):
        hom_from_tags(runs, w, subtract_accidentals=True)
    result, counts = hom_from_tags(runs, w)
    assert result.corrected is None and T3_INPUT in result.labels
    # Declared CW (the sync a clock only), the off-peak windows still show the pulse train: refused
    # (council round 8, Ph3), and under the override computed and labelled. Here it is wrong.
    from dadl.analysis._twopoint import OFFPEAK_ACCEPTED
    from dadl.errors import Misspecified
    from dadl.guardrails import Policy

    with pytest.raises(Misspecified, match="off-peak windows are not flat"):
        hom_from_tags(runs, w, source_type="cw")
    loose = Policy.default().allow("hom_visibility", "misfit_forward_model", "test")
    declared, _ = hom_from_tags(runs, w, source_type="cw", policy=loose)
    assert declared.corrected is not None and T3_INPUT in declared.labels
    assert OFFPEAK_ACCEPTED in declared.corrected_labels
    assert counts["accidentals"] < 0.5 * counts["c_min"]
    # T2 input (no sync period) is not labelled T3, and the off-peak check refuses it the same way.
    t2 = {k: RunTags(v, metadata={"mode": "T2", "tick_ps": 1}) for k, v in plain.items()}
    for r in (t2, plain):
        with pytest.raises(Misspecified, match="off-peak windows are not flat"):
            hom_from_tags(r, w)
        res, _ = hom_from_tags(r, w, policy=loose)
        assert res.corrected is not None and T3_INPUT not in res.labels


def test_hom_takes_a_delay_artifact_in_either_channel_order(tmp_path):
    # Council round 7 (Q2): every artifact in the tests was in analysis order, every window delay 0;
    # a sign error in the capture moved V by about 15 se unnoticed.
    from test_forward import CA, CB, SIGN_DET, SIGN_LATE, assert_same_correction, sign_sets

    from dadl.guardrails import Policy
    from dadl.sim import simulate_pairs

    rng = np.random.default_rng(30)
    policy = Policy.default().allow("hom_visibility", "window_capture", "test")
    sets = sign_sets(tmp_path, "hom_visibility", policy, {(CA, CB): simulate_pairs(2e5, SIGN_DET, SIGN_LATE, 0.05, rng)})
    src = PairSource(pair_rate_per_s=2e5, visibility=0.92, coherence_time_ps=1.5)
    runs = {k: simulate_hom(src, SIGN_DET, SIGN_LATE, d, 0.2, rng) for k, d in (("dip", 0.0), ("reference", 15.0))}
    window = Window(width_ps=200, source="plan", delays_ps={"B": 400})
    ab, ba = (hom_from_tags(runs, window, policy=policy, calibrations=c, channels=(CA, CB), seed=2,
                            reference_delay_ps=15.0) for c in sets)
    assert ab[0].corrections == ("accidental_subtraction", "window_capture")
    assert_same_correction(ab, ba)
    c = ab[0].corrected
    assert abs(c.value - 0.92) < 4 * c.stderr


def test_the_predictive_check_looks_at_the_declared_delay(tmp_path):
    # Council round 7 (Q2): crosstalk 100 ps beside the pair peak at the declared +400 ps is
    # flagged; checked at -400 ps (a sign error) the window holds no peak and passes.
    from test_forward import CA, CB, SIGN_DET, SIGN_LATE, sign_sets

    from dadl.errors import Misspecified
    from dadl.guardrails import Policy
    from dadl.sim import Crosstalk, simulate_pairs

    rng = np.random.default_rng(33)
    policy = Policy.default().allow("hom_visibility", "window_capture", "test")
    ab, _ = sign_sets(tmp_path, "hom_visibility", policy, {(CA, CB): simulate_pairs(2e5, SIGN_DET, SIGN_LATE, 0.05, rng)})
    src = PairSource(pair_rate_per_s=2e5, visibility=0.92, coherence_time_ps=1.5)
    xt = (Crosstalk("A", "B", 0.02, 500.0, 20.0),)
    runs = {k: simulate_hom(src, SIGN_DET, SIGN_LATE, d, 0.2, rng, xt) for k, d in (("dip", 0.0), ("reference", 15.0))}
    window = Window(width_ps=200, source="plan", delays_ps={"B": 400})
    with pytest.raises(Misspecified):
        hom_from_tags(runs, window, policy=policy, calibrations=ab, channels=(CA, CB), seed=2, reference_delay_ps=15.0)


# --- Council round 8: pulsed light on T2 and plain tags (Ph3, Fix 5) ----------------------------


def _pulsed_plain(period, n_pulses, rng, jitter=50.0):
    """A HOM pair of acquisitions from a pulsed source as plain tags (no T3 metadata): 1% singles
    per pulse on each detector, and same-pulse pairs at 0.2% per pulse in the reference only."""
    def run(p_pair):
        k = np.arange(n_pulses)
        a, b, pr = (rng.random(n_pulses) < p for p in (0.01, 0.01, p_pair))
        out = {}
        for name, hit in (("A", a | pr), ("B", b | pr)):
            t = k[hit] * period + rng.normal(0, jitter, int(hit.sum())) + 1000
            out[name] = np.unique(t.astype(np.int64))
        return out

    return {"dip": run(0.0), "reference": run(0.002)}


def test_pulsed_light_on_plain_tags_is_refused_off_peak():
    # r8-council pulsed_t2.py: at a 13157.9 ps period the side windows sample random phases of the
    # pulse train, the accidentals read low and V came out 0.96 for a true 1, unlabelled.
    from dadl.analysis._twopoint import OFFPEAK_ACCEPTED
    from dadl.errors import Misspecified
    from dadl.guardrails import Policy

    w = Window(1000, "plan")
    runs = _pulsed_plain(13_157.9, 2_000_000, np.random.default_rng(2))
    with pytest.raises(Misspecified, match="off-peak windows are not flat"):
        hom_from_tags(runs, w)
    raw, counts = hom_from_tags(runs, w, subtract_accidentals=False)
    assert raw.corrected is None and "offpeak_p" not in counts
    loose = Policy.default().allow("hom_visibility", "misfit_forward_model", "test")
    res, counts = hom_from_tags(runs, w, policy=loose)
    assert OFFPEAK_ACCEPTED in res.corrected_labels and counts["offpeak_p"] < 2e-4


def test_pulsed_light_at_a_period_on_the_window_grid_is_not_refused():
    # Every shift lands on a pulse: the off-peak estimate is right and the windows are flat.
    runs = _pulsed_plain(12_500.0, 2_000_000, np.random.default_rng(2))
    _, counts = hom_from_tags(runs, Window(1000, "plan"))
    assert counts["offpeak_p"] > 0.01


@pytest.mark.slow
def test_the_off_peak_check_is_calibrated_on_clean_hom_snspd():
    # Council round 8 (Fix 5): off-peak p below 0.01 in at most 2% of 200 clean hom-snspd analyses
    # (council: 0.25% at the Bonferroni over the two acquisitions).
    from dadl.analysis.propagation import REGIMES

    reg = REGIMES["hom-snspd"]
    rng = np.random.default_rng(31)
    w = Window(reg.window_ps, "plan")
    ps = []
    for _ in range(200):
        runs = {k: simulate_hom(reg.source, reg.detector, reg.detector, d, reg.duration_s, rng)
                for k, d in (("dip", 0.0), ("reference", reg.reference_delay_ps))}
        ps.append(hom_from_tags(runs, w)[1]["offpeak_p"])
    assert np.mean(np.array(ps) < 0.01) <= 0.02
