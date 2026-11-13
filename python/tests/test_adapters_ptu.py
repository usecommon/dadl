"""PicoQuant PTU (T2) adapter, against fixtures written by tests/ptu_writer.py."""

import numpy as np
import pytest
from ptu_writer import write_ptu

from dadl.adapters import Adapter, PtuFile, read_ptu
from dadl.adapters import ptu
from dadl.data import ChannelId
from dadl.errors import AdapterError

HH_TYPES = [ptu.HYDRAHARP1_T2, ptu.HYDRAHARP2_T2, ptu.TIMEHARP260N_T2, ptu.TIMEHARP260P_T2, ptu.MULTIHARP_T2]


def _events(rng, inputs, n, span):
    ev = []
    truth = {}
    for inp in inputs:
        ticks = np.sort(rng.choice(span, n, replace=False))
        truth[inp] = ticks
        ev += [(int(t), "tag", inp) for t in ticks]
    return ev, truth


@pytest.mark.parametrize("rtype", [ptu.PICOHARP_T2, *HH_TYPES])
def test_round_trip_across_many_rollovers(tmp_path, rtype):
    rng = np.random.default_rng(rtype & 0xFFFF)
    _, wrap = ptu.FORMATS[rtype]
    inputs = (0, 1) if rtype == ptu.PICOHARP_T2 else (0, 1, 3)
    ev, truth = _events(rng, inputs, 2000, 40 * wrap)
    res = 4e-12 if rtype == ptu.PICOHARP_T2 else 1e-12
    path = write_ptu(tmp_path / "f.ptu", rtype, ev, resolution_s=res)
    reader = PtuFile(path)
    tags = reader.read()
    scale = round(res * 1e12)
    assert sorted(ch.input for ch in tags.channels) == sorted(inputs)
    for ch in tags.channels:
        np.testing.assert_array_equal(tags[ch], truth[ch.input] * scale)
    counts = reader.record_counts()
    assert counts["tags"] == len(ev)
    assert counts["tags"] + counts["overflow_records"] + counts["markers"] + counts["unknown"] == counts["total"]
    (roll,) = [e for e in reader.events() if e.kind == "rollover"]
    assert roll.value >= 39
    assert tags.metadata["exact_ps"] is True
    # Duration from the data (last tag plus one tick), not the requested acquisition time.
    assert tags.duration_ps == max(int(tags[c][-1]) for c in tags.channels) + scale


def test_hh2_single_and_multi_overflow_records_agree(tmp_path):
    rng = np.random.default_rng(1)
    ev, truth = _events(rng, (1,), 500, 300 * ptu.HH2_WRAP)
    multi = read_ptu(write_ptu(tmp_path / "m.ptu", ptu.HYDRAHARP2_T2, ev))
    single = read_ptu(write_ptu(tmp_path / "s.ptu", ptu.HYDRAHARP2_T2, ev, hh2_single_overflows=True))
    (ch,) = multi.channels
    np.testing.assert_array_equal(multi[ch], truth[1])
    np.testing.assert_array_equal(single[ch], truth[1])


def test_markers_are_events_not_tags(tmp_path):
    ev = [(100, "tag", 1), (5_000, "marker", 2), (9_000, "tag", 1), (12_000, "marker", 1)]
    reader = PtuFile(write_ptu(tmp_path / "f.ptu", ptu.HYDRAHARP2_T2, ev))
    tags = reader.read()
    assert [int(x) for x in tags[ChannelId("HydraHarp-1234567", 1)]] == [100, 9_000]
    markers = [(e.time_ps, e.value) for e in reader.events() if e.kind == "marker"]
    assert markers == [(5_000, 2), (12_000, 1)]


def test_header_is_provenance(tmp_path):
    reader = PtuFile(write_ptu(tmp_path / "f.ptu", ptu.MULTIHARP_T2, [(1, "tag", 1)]))
    cfg = reader.config()
    assert cfg["HWInpChan_TrgLevel(1)"] == -60 and cfg["UsrHeadName"] == "test fixture"
    info = reader.describe()
    assert info.model == "HydraHarp" and info.serial == "1234567" and info.resolution_ps == pytest.approx(1.0)
    assert isinstance(reader, Adapter)
    assert reader.capabilities().replay and not reader.capabilities().live


def test_non_integer_resolution_rounds_and_says_so(tmp_path):
    tags = read_ptu(write_ptu(tmp_path / "f.ptu", ptu.HYDRAHARP2_T2, [(3, "tag", 1)], resolution_s=2.5e-12))
    (ch,) = tags.channels
    assert int(tags[ch][0]) == 8  # 7.5 ps rounds to even
    assert tags.metadata["exact_ps"] is False and "rounded" in tags.metadata["note"]


def test_truncated_or_padded_files_are_refused(tmp_path):
    ev = [(10, "tag", 1), (20, "tag", 1)]
    with pytest.raises(AdapterError, match="truncated"):
        read_ptu(write_ptu(tmp_path / "a.ptu", ptu.HYDRAHARP2_T2, ev, n_records=3))
    with pytest.raises(AdapterError, match="trailing"):
        read_ptu(write_ptu(tmp_path / "b.ptu", ptu.HYDRAHARP2_T2, ev, trailing=b"\x00" * 4))


def test_t3_and_unknown_formats_are_refused(tmp_path):
    p = write_ptu(tmp_path / "t3.ptu", ptu.HYDRAHARP2_T2, [], extra=[])
    data = p.read_bytes().replace((ptu.HYDRAHARP2_T2).to_bytes(8, "little", signed=True), (0x01010304).to_bytes(8, "little"))
    p.write_bytes(data)
    with pytest.raises(AdapterError, match="T3"):
        PtuFile(p)
    (tmp_path / "x.ptu").write_bytes(b"NOTAPTU!" + b"\x00" * 64)
    with pytest.raises(AdapterError, match="magic"):
        PtuFile(tmp_path / "x.ptu")
    (tmp_path / "y.ptu").write_bytes(ptu.MAGIC + b"1.0.00\x00\x00" + b"\x00" * 20)
    with pytest.raises(AdapterError, match="truncated"):
        PtuFile(tmp_path / "y.ptu")


def test_out_of_order_records_are_refused(tmp_path):
    p = write_ptu(tmp_path / "f.ptu", ptu.HYDRAHARP2_T2, [(10, "tag", 1), (20, "tag", 1)])
    data = bytearray(p.read_bytes())
    a, b = data[-8:-4], data[-4:]
    data[-8:] = b + a
    p.write_bytes(bytes(data))
    with pytest.raises(AdapterError, match="time order"):
        read_ptu(p)


@pytest.mark.parametrize("rtype", [ptu.PICOHARP_T2, ptu.HYDRAHARP1_T2, ptu.HYDRAHARP2_T2, ptu.MULTIHARP_T2])
def test_times_agree_with_ptufile(tmp_path, rtype):
    """Independent cross-check, run only where the third-party `ptufile` package is installed.

    Channel numbers are not compared: DADL follows PicoQuant's demo code (HydraHarp
    sync = 0, input k = k + 1), which ptufile does not.
    """
    ptufile = pytest.importorskip("ptufile")
    rng = np.random.default_rng(3)
    _, wrap = ptu.FORMATS[rtype]
    ev, _ = _events(rng, (0, 1), 300, 20 * wrap)
    path = write_ptu(tmp_path / "f.ptu", rtype, ev)
    ours = read_ptu(path)
    with ptufile.PtuFile(path) as f:
        rec = f.decode_records()
    theirs = np.sort(rec["time"][rec["channel"] >= 0].astype(np.int64))
    mine = np.sort(np.concatenate([ours[c] for c in ours.channels]))
    np.testing.assert_array_equal(mine, theirs)


@pytest.mark.parametrize("code, lost", [(0, False), (1, False), (2, False), (3, True), (-2, True), (-13, True), (-99, True)])
def test_stop_reason_marks_data_loss(tmp_path, code, lost):
    # Expert review W5: a FIFO overrun or lost clock ended the acquisition early with data lost;
    # before that review, such files were read as complete.
    path = write_ptu(tmp_path / "f.ptu", ptu.HYDRAHARP2_T2, [(10, "tag", 1), (20, "tag", 1)],
                     acquisition_ms=360_000_000.0, extra=[("TTResult_StopReason", "int", code)])
    reader = PtuFile(path)
    tags = reader.read()
    assert tags.metadata["incomplete"] is lost
    assert any(e.kind == ("data-loss" if lost else "stop-reason") for e in reader.events())
    assert tags.duration_ps == 21  # from the data, not the requested 100 hours


def test_incomplete_files_are_refused_in_a_folder(tmp_path):
    from dadl.errors import ConfigError
    from dadl.tagfile import read_tags

    write_ptu(tmp_path / "setting0.ptu", ptu.HYDRAHARP2_T2, [(10, "tag", 1)], extra=[("TTResult_StopReason", "int", -2)])
    with pytest.raises(ConfigError, match="FIFO overrun"):
        read_tags(tmp_path, {1: "A"}, ("s",))


def test_picoharp_channels_beyond_its_inputs_are_refused(tmp_path):
    path = write_ptu(tmp_path / "f.ptu", ptu.PICOHARP_T2, [(10, "tag", 5)])
    with pytest.raises(AdapterError, match="valid inputs"):
        read_ptu(path)
