"""DADL's hardware twin: simulated PTU files and a simulated live tagger with faults."""

import numpy as np
import pytest
from ptu_writer import encode, write_ptu as write_fixture

from dadl.adapters import Adapter, Faults, PtuFile, SimulatedTagger, check_accounting, ptu
from dadl.adapters.ptu_write import DEVICES, encode_records, write_ptu
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
