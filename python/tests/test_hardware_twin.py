"""DADL's hardware twin: simulated PTU files and a simulated live tagger with faults."""

import numpy as np
import pytest
from ptu_writer import encode, encode_t3, write_ptu as write_fixture

from dadl.adapters import Adapter, Faults, PtuFile, SimulatedTagger, check_accounting, ptu
from dadl.adapters.ptu_write import DEVICES, MAX_INPUT, T3_RECORD_TYPES, encode_records, encode_t3_records, write_ptu
from dadl.cli import example_path, main
from dadl.config import load_setup
from dadl.errors import ConfigError
from dadl.tagfile import read_tags, tags_are_simulated, write_ptu_dir

FAMILIES = [ptu.PICOHARP_T2, ptu.HYDRAHARP1_T2, ptu.HYDRAHARP2_T2, ptu.MULTIHARP_T2]


@pytest.mark.parametrize("rtype", FAMILIES)
def test_vectorized_encoder_matches_the_record_by_record_writer(rtype):
    rng = np.random.default_rng(rtype & 0xFF)
    _, wrap = ptu.FORMATS[rtype]
    far = (1 << 25) + 7 if rtype == ptu.HYDRAHARP2_T2 else 9  # an HH2 gap needing two count records
    ticks = np.sort(np.concatenate([rng.integers(0, 40 * wrap, 2000), [50 * wrap, 50 * wrap, (50 + far) * wrap + 3]]))
    marker = rng.random(ticks.size) < 0.05
    inputs = np.where(marker, rng.integers(1, 15, ticks.size), rng.integers(0, 3, ticks.size))
    events = [(int(t), "marker" if m else "tag", int(i)) for t, i, m in zip(ticks, inputs, marker)]
    np.testing.assert_array_equal(encode_records(rtype, ticks, inputs, marker), np.array(encode(rtype, events), dtype=np.uint32))


@pytest.mark.parametrize("rtype", [ptu.PICOHARP_T3, ptu.HYDRAHARP1_T3, ptu.HYDRAHARP2_T3, ptu.MULTIHARP_T3])
def test_vectorized_t3_encoder_matches_the_record_by_record_writer(rtype):
    rng = np.random.default_rng(rtype & 0xFFFFF)
    family, wrap = ptu.T3_FORMATS[rtype]
    far = 1023 + 7 if family == "hh2" else 9  # an HH2 gap needing two count records
    nsync = np.sort(np.concatenate([rng.integers(0, 40 * wrap, 2000), [50 * wrap, 50 * wrap, (50 + far) * wrap + 3]]))
    marker = rng.random(nsync.size) < 0.05
    inputs = np.where(marker, rng.integers(1, 15, nsync.size), rng.integers(1, 5, nsync.size))
    dtime = np.where(marker, 0, rng.integers(0, 4096, nsync.size))
    events = [(int(s), "marker" if m else "tag", int(i), int(d)) for s, i, m, d in zip(nsync, inputs, marker, dtime)]
    np.testing.assert_array_equal(encode_t3_records(rtype, nsync, dtime, inputs, marker),
                                  np.array(encode_t3(rtype, events), dtype=np.uint32))


@pytest.mark.parametrize("device", sorted(DEVICES))
def test_written_t3_file_reads_back_quantized(tmp_path, device):
    _, tick, _ = DEVICES[device]
    period = 12_500
    rng = np.random.default_rng(1)
    tags = {inp: np.sort(rng.choice(10**10, 3000, replace=False)) for inp in range(1, min(MAX_INPUT[device], 2) + 1)}
    reader = PtuFile(write_ptu(tmp_path / "f.ptu", tags, device=device, sync_period_ps=period, markers=[(10**9, 3)]))
    got = reader.read()
    assert reader.mode == "T3" and reader.record_type == T3_RECORD_TYPES[device]
    assert got.metadata["simulated"] is True and got.metadata["exact_ps"] is True
    by_input = {ch.input: got[ch] for ch in got.channels}
    for inp, t in tags.items():
        np.testing.assert_array_equal(by_input[inp], t // period * period + t % period // tick * tick)
    np.testing.assert_array_equal(by_input[0], np.unique(np.concatenate(list(tags.values())) // period) * period)
    assert [(e.time_ps, e.value) for e in reader.events() if e.kind == "marker"] == [(10**9 // period * period, 3)]


def test_t3_writer_refuses_sync_tags_and_long_periods(tmp_path):
    with pytest.raises(ValueError, match="input 0"):
        write_ptu(tmp_path / "a.ptu", {0: np.array([1, 2]), 1: np.array([5])}, sync_period_ps=12_500)
    with pytest.raises(ValueError, match="dtime range"):
        write_ptu(tmp_path / "b.ptu", {1: np.array([5])}, device="picoharp", sync_period_ps=100_000)


@pytest.mark.parametrize("device", sorted(DEVICES))
def test_written_file_reads_back_quantized_and_marked_simulated(tmp_path, device):
    _, tick, hw = DEVICES[device]
    rng = np.random.default_rng(0)
    tags = {0: np.sort(rng.integers(0, 10**12, 5000)), 1: np.sort(rng.integers(0, 10**12, 5000))}
    reader = PtuFile(write_ptu(tmp_path / "f.ptu", tags, device=device, markers=[(5 * 10**11, 3)]))
    got = reader.read()
    assert got.metadata["simulated"] is True
    assert reader.describe().model == hw and reader.config()["CreatorSW_Name"] == "DADL"
    for ch in got.channels:
        np.testing.assert_array_equal(got[ch], tags[ch.input] // tick * tick)
    assert [e.value for e in reader.events() if e.kind == "marker"] == [3]


def test_recordings_are_not_marked_simulated(tmp_path):
    path = write_fixture(tmp_path / "setting0.ptu", ptu.HYDRAHARP2_T2, [(10, "tag", 1)])
    assert PtuFile(path).read().metadata["simulated"] is False
    assert not tags_are_simulated(tmp_path)


def test_negative_tags_are_refused(tmp_path):
    from dadl.errors import AdapterError

    with pytest.raises(AdapterError, match="negative"):
        write_ptu(tmp_path / "f.ptu", {1: np.array([-5, 10])})


def test_ptu_folder_round_trip_gives_the_same_result(tmp_path, capsys):
    assert main(["simulate", "chsh-spad", "--out", str(tmp_path / "t.npz")]) == 0
    assert main(["simulate", "chsh-spad", "--out", str(tmp_path / "ptu"), "--format", "ptu"]) == 0
    assert tags_are_simulated(tmp_path / "ptu")
    capsys.readouterr()
    assert main(["analyze", "chsh-spad", "--tags", str(tmp_path / "t.npz")]) == 0
    direct = capsys.readouterr().out.split("Raw S = ")[1].split()[0]
    assert main(["analyze", "chsh-spad", "--tags", str(tmp_path / "ptu")]) == 0
    via_ptu = capsys.readouterr().out.split("Raw S = ")[1].split()[0]
    assert direct == via_ptu


def test_t3_ptu_folder_reads_back_the_simulated_tags(tmp_path):
    """A HydraHarp T3 folder (1 ps dtime, 12.5 ns sync) returns the simulated tags exactly."""
    from dadl.sim import DetectorParams, PairSource, simulate_hom
    from dadl.tagfile import ptu_file, tag_tick_ps, tagger_settings

    det = DetectorParams(efficiency=0.8, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=25.0)
    rng = np.random.default_rng(9)
    runs = {key: simulate_hom(PairSource(pair_rate_per_s=1e5), det, det, delay, 0.05, rng)
            for key, delay in (("dip", 0.0), ("reference", 50.0))}
    runs = {k: {ch: t[t >= 0] for ch, t in run.items()} for k, run in runs.items()}
    for i, key in enumerate(("dip", "reference")):
        write_ptu(ptu_file(tmp_path, i), {1: runs[key]["A"], 2: runs[key]["B"]}, sync_period_ps=12_500)
    # Input 0 (derived syncs) is left out unless mapped; the detector inputs read back exactly.
    got = read_tags(tmp_path, {1: "A", 2: "B"}, ("dip", "reference"))
    for key in runs:
        for ch in ("A", "B"):
            np.testing.assert_array_equal(got[key][ch], runs[key][ch])
    assert tag_tick_ps(tmp_path) == 1
    assert tagger_settings(tmp_path)["tagger_MeasDesc_Resolution"] == pytest.approx(1e-12)


def test_calibration_from_simulated_ptu_is_recorded_as_sim(tmp_path, capsys):
    import json

    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "5", "--out", str(tmp_path / "d"), "--format", "ptu"]) == 0
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(tmp_path / "d"), "--channel", "A+", "--out", str(tmp_path / "art")]) == 0
    (manifest,) = (tmp_path / "art").glob("*[0-9a-f].json")
    assert json.loads(manifest.read_text())["backend"] == "sim"


def test_ptu_folder_errors(tmp_path):
    setup = load_setup(example_path("hom-snspd"))
    runs = {"dip": {"A": np.array([1, 2]), "B": np.array([3])}, "reference": {"A": np.array([5]), "B": np.array([6])}}
    write_ptu_dir(tmp_path, runs, setup.channels, setup.setting_keys)
    got = read_tags(tmp_path, setup.channels, setup.setting_keys)
    assert got["dip"]["A"].tolist() == [1, 2] and got["reference"]["B"].tolist() == [6]
    with pytest.raises(ConfigError, match="not in the setup's channel map"):
        read_tags(tmp_path, {0: "A"}, setup.setting_keys)
    (tmp_path / "setting1.ptu").unlink()
    with pytest.raises(ConfigError, match="missing setting1.ptu"):
        read_tags(tmp_path, setup.channels, setup.setting_keys)


# --- simulated live tagger ---------------------------------------------------------------------


@pytest.fixture(scope="module")
def hom_setup():
    return load_setup(example_path("hom-snspd"))


def test_live_tagger_without_faults_delivers_everything(hom_setup):
    dev = SimulatedTagger(hom_setup, "reference", batch_s=0.25)
    assert isinstance(dev, Adapter) and dev.capabilities().live
    batches = list(dev.stream())
    assert len(batches) == 4 and all(b.metadata["simulated"] for b in batches)
    tags = dev.read()
    truth = dev.truth()
    for ch in tags.channels:
        np.testing.assert_array_equal(tags[ch], truth[ch.input])
    assert dev.events() == ()


def test_every_lost_tag_is_reported(hom_setup):
    faults = Faults(
        buffer_overflows=((0.10, 0.02), (0.48, 0.05)),  # the second starts in dropped batch 1
        dropped_batches=(1,),
        config_changes=((0.6, "HWInpChan_TrgLevel(0)", -60),),
    )
    dev = SimulatedTagger(hom_setup, "dip", batch_s=0.25, faults=faults)
    tags = dev.read()
    kinds = [e.kind for e in dev.events()]
    assert kinds.count("buffer-overflow") == 2 and kinds.count("dropped-batch") == 1 and "config-change" in kinds
    assert all(v == 0 for v in check_accounting(tags, dev.events(), dev.truth()).values())
    lost = sum(e.value for e in dev.events() if e.kind != "config-change")
    assert lost > 0 and sum(tags[c].size for c in tags.channels) + lost == sum(a.size for a in dev.truth().values())


def test_silent_loss_is_detected(hom_setup):
    dev = SimulatedTagger(hom_setup, "dip", faults=Faults(dropped_batches=(2,)))
    tags = dev.read()
    unreported = tuple(e for e in dev.events() if e.kind != "dropped-batch")
    assert any(v > 0 for v in check_accounting(tags, unreported, dev.truth()).values())


def test_overlapping_overflows_are_refused(hom_setup):
    with pytest.raises(ValueError, match="overlap"):
        SimulatedTagger(hom_setup, faults=Faults(buffer_overflows=((0.1, 0.2), (0.2, 0.1))))


@pytest.mark.parametrize("device", ["hydraharp", "picoharp", "multiharp", "timeharp260p"])
def test_delay_fit_is_unbiased_on_every_device_tick(device):
    # Found by running the existing calibrations on twin data: method v2 assumed 1 ps ticks
    # and read 5 ps MultiHarp delays 2 ps high.
    from dadl.calibration import calibrate_delay
    from dadl.calibration.delay import tag_tick
    from dadl.data import ChannelId
    from dadl.sim import DetectorParams, simulate_pairs

    _, tick, _ = DEVICES[device]
    det = DetectorParams(efficiency=0.85, dark_rate_cps=100, dead_time_ps=30_000, jitter_sigma_ps=25)
    est, notes = [], set()
    for seed in range(8):
        a, b = simulate_pairs(1e5, det, det, 1.0, np.random.default_rng(seed))
        a, b = a // tick * tick, b // tick * tick
        assert tag_tick(a, b) == tick
        cal = calibrate_delay(a, b, (ChannelId("A", 0), ChannelId("B", 1)), n_draws=10)
        est.append(cal.estimates["delay_ps"].value)
        notes.update(cal.notes)
        assert cal.conditions["tick_ps"] == tick
    assert abs(np.mean(est)) < 0.3
    assert any("raised" in n for n in notes) == (20 % tick != 0)


def test_accounting_compares_times_not_counts(hom_setup):
    # Time-tagging review W11: count-only accounting passed lost-plus-duplicated, shifted, swapped
    # and extra data. Every one of those must now show up.
    from dadl.data import ChannelId, TagSet

    dev = SimulatedTagger(hom_setup, "dip")
    good = dev.read()
    truth = dev.truth()
    assert all(v == 0 for v in check_accounting(good, (), truth).values())
    (c0, c1) = good.channels

    def variant(arr0, arr1):
        return TagSet({c0: np.sort(arr0), c1: np.sort(arr1)})

    a, b = good[c0], good[c1]
    cases = {
        "lost and duplicated": variant(np.concatenate([a[500:], a[:500][::-1][:0], np.repeat(a[500:1000], 1)]), b),
        "shifted": variant(a + 1_000_000, b),
        "swapped": variant(b, a),
    }
    for name, ts in cases.items():
        assert sum(check_accounting(ts, (), truth).values()) > 0, name
    extra = TagSet({c0: a, c1: b, ChannelId(c0.device, 7): np.array([5])})
    assert check_accounting(extra, (), truth)[7] == 1


def test_presets_refuse_inputs_the_device_lacks(tmp_path):
    with pytest.raises(ValueError, match="do not exist"):
        write_ptu(tmp_path / "f.ptu", {0: np.array([1]), 3: np.array([2])}, device="picoharp")
    write_ptu(tmp_path / "g.ptu", {0: np.array([1]), 1: np.array([2])}, device="picoharp")


def test_declared_tick_overrides_a_spoiled_gcd():
    # Time-tagging review Q6: an 80 ps tick with one channel offset by 30 ps made the GCD 10 ps.
    from dadl.calibration import calibrate_delay
    from dadl.data import ChannelId
    from dadl.sim import DetectorParams, simulate_pairs

    det = DetectorParams(efficiency=0.85, dark_rate_cps=100, dead_time_ps=30_000, jitter_sigma_ps=60)
    a, b = simulate_pairs(1e5, det, det, 1.0, np.random.default_rng(2))
    a, b = a // 80 * 80, b // 80 * 80 + 30
    cal = calibrate_delay(a, b, (ChannelId("A", 0), ChannelId("B", 1)), tick_ps=80, n_draws=10)
    assert cal.conditions["tick_ps"] == 80 and any("tick grid" in n for n in cal.notes)


# --- Council round 7: T3 input (tick, sync period, pulsed sources, dtime range, divider) ---------


def _patch_float(path, key, value):
    """Overwrite a float header entry in place (as the reviewers' realhead.py did)."""
    import struct
    from pathlib import Path

    data = bytearray(Path(path).read_bytes())
    i = data.index(key.encode().ljust(32, b"\0"))
    struct.pack_into("<d", data, i + 40, value)
    Path(path).write_bytes(bytes(data))


def _t3_folder(d, sync_hz=(None, None), device="timeharp260p", n=2000):
    """Two T3 files (25 ps dtime on the TimeHarp 260 P) at a 12.5 ns sync, the period then set to
    1 / sync rate as an instrument measuring it would write it."""
    from dadl.tagfile import ptu_file

    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    for i, hz in enumerate(sync_hz):
        a, b = np.sort(rng.integers(0, 10**10, n)), np.sort(rng.integers(0, 10**10, n))
        p = write_ptu(ptu_file(d, i), {1: a, 2: b}, device=device, sync_period_ps=12_500)
        if hz is not None:
            _patch_float(p, "MeasDesc_GlobalResolution", 1.0 / hz)
    return d


def test_the_t3_tick_is_the_dtime_resolution(tmp_path):
    # G1: a measured (non-whole) sync period left a 1 ps tick, so tau_np widened by 1 ps on 25 ps data.
    from dadl.tagfile import tag_tick_ps

    whole = _t3_folder(tmp_path / "whole")
    assert tag_tick_ps(whole) == 25 and PtuFile(whole / "setting0.ptu").tick_ps == 25
    measured = _t3_folder(tmp_path / "measured", (79_999_870.0, 79_999_870.0))
    # The dtime resolution, plus 1 ps for the rounding to whole ps.
    assert tag_tick_ps(measured) == 26
    tags = read_tags(measured, {1: "A", 2: "B"}, ("dip", "reference"))
    assert tags["dip"].metadata["tick_ps"] == 26
    assert tags["dip"].metadata["sync_period_ps"] == pytest.approx(1e12 / 79_999_870.0)


def test_a_drifting_sync_period_is_one_setting_and_another_rate_is_not(tmp_path):
    # G2: two files 10 Hz apart at 80 MHz were refused as different tagger settings.
    from dadl.propagate import condition_mismatches
    from dadl.tagfile import tagger_settings

    a = tagger_settings(_t3_folder(tmp_path / "a", (79_999_870.0, 79_999_860.0)))
    b = tagger_settings(_t3_folder(tmp_path / "b", (79_999_880.0, 79_999_880.0)))
    assert condition_mismatches(a, b) == []
    with pytest.raises(ConfigError, match=r"1.250002031e-08 in setting0.ptu, 1.315789474e-08 in setting1.ptu"):
        tagger_settings(_t3_folder(tmp_path / "c", (79_999_870.0, 76e6)))
    (msg,) = condition_mismatches(a, tagger_settings(_t3_folder(tmp_path / "d", (76e6, 76e6))))
    assert "1.250002031e-08" in msg


def test_a_t3_dtime_range_shorter_than_the_sync_period_is_refused(tmp_path):
    # G4: a HydraHarp at 1 ps covers 32.8 ns; with a 100 ns sync, late photons are never recorded.
    from dadl.errors import AdapterError
    from dadl.tagfile import ptu_file

    rng = np.random.default_rng(0)
    t = np.sort(rng.integers(0, 10**10, 5000))
    p = write_ptu(ptu_file(tmp_path, 0), {1: t, 2: t + 500}, device="hydraharp", sync_period_ps=12_500)
    _patch_float(p, "MeasDesc_GlobalResolution", 100e-9)
    reader = PtuFile(p)
    tagset = reader.read()
    assert "t3-range" in [e.kind for e in reader.events()]
    assert "dtime range 32.77 ns < sync period 100 ns" in tagset.metadata["note"]
    with pytest.raises(AdapterError, match="coarser"):
        read_tags(tmp_path, {1: "A", 2: "B"}, ("dip",))


def _t3_hom(d, extra=(), declared=None):
    """A HOM experiment as two HydraHarp T3 files (1 ps dtime, 12.5 ns sync, marked simulated),
    detectors on inputs 1 and 2, and its setup.yaml (source.type as `declared`)."""
    from pathlib import Path

    from dadl.sim import DetectorParams, PairSource, simulate_hom
    from dadl.tagfile import ptu_file

    det = DetectorParams(efficiency=0.8, dark_rate_cps=2e4, dead_time_ps=30_000, jitter_sigma_ps=25.0)
    rng = np.random.default_rng(9)
    for i, delay in enumerate((0.0, 15.0)):
        run = simulate_hom(PairSource(pair_rate_per_s=1e5, coherence_time_ps=1.5), det, det, delay, 0.05, rng)
        ev = []
        for inp, t in ((1, run["A"]), (2, run["B"])):
            t = t[t >= 0]
            ev += [(int(x) // 12_500, "tag", inp, int(x) % 12_500) for x in t]
        write_fixture(ptu_file(d, i), ptu.HYDRAHARP2_T3, ev, sync_period_s=12.5e-9,
                      extra=[("DADL_Simulated", "bool", True), *extra])
    text = Path(example_path("hom-snspd")).read_text().replace("channels: {0: A, 1: B}", "channels: {1: A, 2: B}")
    text = text.replace("  type: spdc-cw\n", f"  type: {declared}\n" if declared else "")
    (d / "setup.yaml").write_text(text)
    return d / "setup.yaml"


def test_t3_hom_corrects_only_from_a_declared_cw_source(tmp_path, capsys):
    # G3: a pulsed T3 HOM read V 0.962 for 1, unlabelled; the source must be declared CW.
    from dadl.guardrails import T3_INPUT

    setup = _t3_hom(tmp_path)
    assert main(["analyze", str(setup), "--tags", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert T3_INPUT in out and "Accidentals not subtracted" in out
    _t3_hom(tmp_path, declared="cw")
    assert main(["analyze", str(setup), "--tags", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert T3_INPUT in out and "Accidentals subtracted" in out


def test_a_divided_t3_sync_is_refused_for_dark_and_labelled_in_analyses(tmp_path, capsys):
    # G5: whether a divided T3 file's period counts divided syncs is unverified (M9).
    setup = _t3_hom(tmp_path, extra=[("HWSync_Divider", "int", 2)], declared="cw")
    assert main(["calibrate", "dark", str(setup), "--tags", str(tmp_path), "--out", str(tmp_path / "art")]) == 1
    assert "T3 sync divider convention is unverified" in capsys.readouterr().err
    assert main(["calibrate", "deadtime", str(setup), "--tags", str(tmp_path), "--out", str(tmp_path / "art")]) == 1
    assert "T3 sync divider convention is unverified" in capsys.readouterr().err
    assert main(["analyze", str(setup), "--tags", str(tmp_path)]) == 0
    assert "T3 sync divider 2: the divider convention is unverified" in capsys.readouterr().out


# --- Council round 8: the T3 lattice for histograms (Fix 1) ------------------------------------


def _t3_pairs_folder(tmp_path, period=13_157.9, written=13_158):
    """A HOM setup on inputs 1 and 2 and a pair run written as TimeHarp 260 P T3 files (25 ps
    dtime) at a whole-ps sync, then given the measured period `period` as an instrument writes it.
    Returns (setup path, folder)."""
    from dadl.sim import DetectorParams, simulate_pairs
    from dadl.tagfile import ptu_file

    setup = tmp_path / "setup.yaml"
    text = example_path("hom-snspd").read_text(encoding="utf-8")
    setup.write_text(text.replace("channels: {0: A, 1: B}", "channels: {1: A, 2: B}")
                     .replace("jitter_sigma_ps: 25", "jitter_sigma_ps: 60"), encoding="utf-8")
    det = DetectorParams(efficiency=0.85, dark_rate_cps=100.0, dead_time_ps=30_000, jitter_sigma_ps=60.0)
    a, b = simulate_pairs(2e5, det, det, 0.1, np.random.default_rng(8))
    a, b = a[a >= 0], b[b >= 0]
    folder = tmp_path / "pairs"
    folder.mkdir()
    for i in (0, 1):
        p = write_ptu(ptu_file(folder, i), {1: a, 2: b}, device="timeharp260p", sync_period_ps=written)
        _patch_float(p, "MeasDesc_GlobalResolution", period * 1e-12)
    return setup, folder


def test_t3_lattice_and_tick_are_kept_apart(tmp_path):
    from dadl.tagfile import tag_lattice_ps, tag_tick_ps

    _, folder = _t3_pairs_folder(tmp_path)
    assert tag_lattice_ps(folder) == 25 and tag_tick_ps(folder) == 26
    tags = read_tags(folder, {1: "A", 2: "B"}, ("dip", "reference"))
    assert tags["dip"].metadata["lattice_ps"] == 25 and tags["dip"].metadata["tick_ps"] == 26


def test_calibrate_pairs_on_t3_with_a_measured_period_covers_the_delay(tmp_path, capsys):
    # Council round 8 (Ph1, St3): calibrate pairs passed the 26 ps gap bound as the tick and the
    # delay read -12 ps with its interval excluding 0, every artifact a misfit.
    from dadl.artifact import load_calibration

    setup, folder = _t3_pairs_folder(tmp_path)
    art = tmp_path / "art"
    assert main(["calibrate", "pairs", str(setup), "--tags", str(folder), "--no-analysers", "--herald", "A",
                 "--signal", "B", "--out", str(art)]) == 0
    out = capsys.readouterr().out
    assert "tick grid" not in out
    (delay,) = [c for c in map(load_calibration, art.glob("*.json")) if c.kind == "delay-jitter"]
    assert delay.estimates["delay_ps"].contains(0.0) and delay.status == "ok"
    assert delay.conditions["tick_ps"] == 25.0


def test_t3_delay_at_the_lattice_is_calibrated():
    # Council round 8 (t3delay.py): 25 ps T3 at a 12500.3 ps sync, the delay fit at the lattice.
    from dadl.adapters.ptu import _t3_to_ps
    from dadl.calibration.delay import calibrate_delay
    from dadl.data import ChannelId
    from dadl.sim import DetectorParams, simulate_pairs

    det = DetectorParams(efficiency=0.5, dark_rate_cps=2000.0, dead_time_ps=22_000, jitter_sigma_ps=100.0)

    def t3(t, period=12_500.3, res=25):
        nsync = np.floor(t / period).astype(np.int64)
        dtime = np.floor((t - nsync * period) / res).astype(np.int64)
        return _t3_to_ps(nsync, dtime, period * 1e-12, res * 1e-12)[0]

    ids = (ChannelId("x", 0), ChannelId("x", 1))
    misfit = covered = 0
    for seed in range(20):
        a, b = simulate_pairs(1e5, det, det, 0.2, np.random.default_rng(500 + seed))
        a, b = a[a >= 0], b[b >= 0]
        cal = calibrate_delay(t3(a), t3(b), ids, backend="sim", seed=1, n_draws=200, tick_ps=25,
                              t3_period_ps=12_500.3)
        misfit += cal.ppc_pvalue < 1e-3
        covered += cal.estimates["jitter_sigma_ps"].contains(100 * 2**0.5)
        assert not any("tick grid" in n for n in cal.notes)
    assert misfit == 0 and covered >= 19
