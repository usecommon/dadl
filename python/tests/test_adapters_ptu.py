"""PicoQuant PTU (T2 and T3) adapter, against fixtures written by tests/ptu_writer.py."""

from fractions import Fraction

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


def test_unknown_formats_are_refused(tmp_path):
    p = write_ptu(tmp_path / "t3.ptu", ptu.HYDRAHARP2_T2, [], extra=[])
    # 0x00010308: a T3 code on a hardware byte no supported device uses.
    data = p.read_bytes().replace((ptu.HYDRAHARP2_T2).to_bytes(8, "little", signed=True), (0x00010308).to_bytes(8, "little"))
    p.write_bytes(data)
    with pytest.raises(AdapterError, match="unsupported PTU record type"):
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


# --- T3 mode ---------------------------------------------------------------------------------

T3_TYPES = [ptu.PICOHARP_T3, ptu.HYDRAHARP1_T3, ptu.HYDRAHARP2_T3, ptu.TIMEHARP260N_T3, ptu.TIMEHARP260P_T3,
            ptu.MULTIHARP_T3]
PERIOD_PS = 12_500  # an 80 MHz sync


def _t3_events(rng, inputs, n, n_sync, max_dtime):
    """Per input, n distinct sync periods with a random dtime each; events in time order."""
    ev, truth = [], {}
    for inp in inputs:
        nsync = np.sort(rng.choice(n_sync, n, replace=False))
        dtime = rng.integers(0, max_dtime, n)
        truth[inp] = (nsync, dtime)
        ev += [(int(s), "tag", inp, int(d)) for s, d in zip(nsync, dtime)]
    ev.sort(key=lambda e: (e[0], e[3]))
    return ev, truth


@pytest.mark.parametrize("rtype", T3_TYPES)
def test_t3_round_trip_across_many_rollovers(tmp_path, rtype):
    rng = np.random.default_rng(rtype & 0xFFFFF)
    _, wrap = ptu.T3_FORMATS[rtype]
    picoharp = rtype == ptu.PICOHARP_T3
    res_ps = 4 if picoharp else 1
    inputs = (1, 2, 4)  # PicoHarp: router channels 1-4; HydraHarp family: record channels 0, 1, 3
    ev, truth = _t3_events(rng, inputs, 1500, 40 * wrap, PERIOD_PS // res_ps)
    path = write_ptu(tmp_path / "f.ptu", rtype, ev, resolution_s=res_ps * 1e-12, sync_period_s=PERIOD_PS * 1e-12)
    reader = PtuFile(path)
    assert reader.mode == "T3" and reader.sync_period_s == pytest.approx(PERIOD_PS * 1e-12)
    tags = reader.read()
    assert sorted(ch.input for ch in tags.channels) == [0, *inputs]
    by_input = {ch.input: tags[ch] for ch in tags.channels}
    for inp, (nsync, dtime) in truth.items():
        np.testing.assert_array_equal(by_input[inp], nsync * PERIOD_PS + dtime * res_ps)
    # The sync channel: one tag per sync period that started a recorded measurement.
    syncs = np.unique(np.concatenate([nsync for nsync, _ in truth.values()]))
    np.testing.assert_array_equal(by_input[0], syncs * PERIOD_PS)
    counts = reader.record_counts()
    assert counts["tags"] == len(ev) and counts["sync_tags"] == syncs.size
    assert counts["tags"] + counts["overflow_records"] + counts["markers"] + counts["unknown"] == counts["total"]
    (roll,) = [e for e in reader.events() if e.kind == "rollover"]
    assert roll.value >= 39
    (sync_event,) = [e for e in reader.events() if e.kind == "t3-sync"]
    assert sync_event.value == syncs.size
    assert tags.metadata["exact_ps"] is True and tags.metadata["mode"] == "T3"
    assert tags.metadata["sync_periods"] == int(syncs.max()) + 1
    assert "T3" in tags.clock
    assert tags.duration_ps == max(int(a[-1]) for a in by_input.values()) + res_ps


def test_t3_hh2_single_and_multi_overflow_records_agree(tmp_path):
    rng = np.random.default_rng(2)
    # Gaps of thousands of wraps need several count records (at most 1023 wraps each).
    ev, truth = _t3_events(rng, (1,), 300, 3000 * ptu.HH_T3_WRAP, 1000)
    multi = PtuFile(write_ptu(tmp_path / "m.ptu", ptu.HYDRAHARP2_T3, ev))
    single = PtuFile(write_ptu(tmp_path / "s.ptu", ptu.HYDRAHARP2_T3, ev, hh2_single_overflows=True))
    a, b = multi.read(), single.read()
    assert multi.record_counts()["overflow_records"] < single.record_counts()["overflow_records"]
    nsync, dtime = truth[1]
    for tags in (a, b):
        (ch,) = [c for c in tags.channels if c.input == 1]
        np.testing.assert_array_equal(tags[ch], nsync * PERIOD_PS + dtime)


@pytest.mark.parametrize("rtype", [ptu.PICOHARP_T3, ptu.HYDRAHARP1_T3, ptu.MULTIHARP_T3])
def test_t3_markers_are_events_at_their_sync(tmp_path, rtype):
    wrap = ptu.T3_FORMATS[rtype][1]
    ev = [(3, "tag", 1, 100), (5, "marker", 2, 0), (wrap + 7, "tag", 1, 50), (3 * wrap + 1, "marker", 9, 0)]
    res_ps = 4 if rtype == ptu.PICOHARP_T3 else 1
    reader = PtuFile(write_ptu(tmp_path / "f.ptu", rtype, ev, resolution_s=res_ps * 1e-12))
    tags = reader.read()
    (ch,) = [c for c in tags.channels if c.input == 1]
    assert tags[ch].tolist() == [3 * PERIOD_PS + 100 * res_ps, (wrap + 7) * PERIOD_PS + 50 * res_ps]
    markers = [(e.time_ps, e.value) for e in reader.events() if e.kind == "marker"]
    assert markers == [(5 * PERIOD_PS, 2), ((3 * wrap + 1) * PERIOD_PS, 9)]
    counts = reader.record_counts()
    assert counts["markers"] == 2 and counts["tags"] == 2
    assert counts["tags"] + counts["overflow_records"] + counts["markers"] + counts["unknown"] == counts["total"]
    # The last marker is the last event: the duration includes it.
    assert tags.duration_ps == (3 * wrap + 1) * PERIOD_PS + res_ps


def test_t3_non_integer_sync_period_rounds_and_says_so(tmp_path):
    period_s = 1 / 76.2e6  # 13123.359... ps
    rng = np.random.default_rng(5)
    # Late syncs (about 10 hours in) check that the period's fraction keeps sub-ps precision.
    nsync = np.sort(rng.choice(np.arange(2_700_000_000_000, 2_700_000_100_000), 200, replace=False))
    dtime = rng.integers(0, 13_000, nsync.size)
    ev = [(int(s), "tag", 1, int(d)) for s, d in zip(nsync, dtime)]
    tags = read_ptu(write_ptu(tmp_path / "f.ptu", ptu.HYDRAHARP2_T3, ev, sync_period_s=period_s))
    (ch,) = [c for c in tags.channels if c.input == 1]
    p, r = Fraction(period_s) * 10**12, Fraction(1e-12) * 10**12
    expected = [round(int(s) * p + int(d) * r) for s, d in zip(nsync, dtime)]
    assert tags[ch].tolist() == expected
    assert tags.metadata["exact_ps"] is False and "rounded" in tags.metadata["note"]


@pytest.mark.parametrize("channel", [0, 5])
def test_t3_picoharp_illegal_channels_are_refused(tmp_path, channel):
    path = write_ptu(tmp_path / "f.ptu", ptu.PICOHARP_T3, [(1, "tag", channel, 10)], resolution_s=4e-12)
    with pytest.raises(AdapterError, match="valid inputs are 1-4"):
        read_ptu(path)


def test_t3_dtime_beyond_the_sync_period_is_refused(tmp_path):
    # A dtime of 5 ns after a 1 ns sync period lands after the next period's tag: out of order.
    ev = [(0, "tag", 1, 5_000), (1, "tag", 1, 0)]
    with pytest.raises(AdapterError, match="time order.*sync period"):
        read_ptu(write_ptu(tmp_path / "f.ptu", ptu.HYDRAHARP2_T3, ev, sync_period_s=1e-9))


def test_t3_unknown_special_records_are_counted_not_dropped(tmp_path):
    p = write_ptu(tmp_path / "f.ptu", ptu.HYDRAHARP2_T3, [(1, "tag", 1, 5), (2, "marker", 1, 0), (3, "tag", 1, 6)])
    marker = ((1 << 31) | (1 << 25) | 2).to_bytes(4, "little")
    channel0 = ((1 << 31) | 2).to_bytes(4, "little")  # special on channel 0: no meaning in T3
    data = p.read_bytes()
    assert data.count(marker) == 1
    p.write_bytes(data.replace(marker, channel0))
    reader = PtuFile(p)
    tags = reader.read()
    counts = reader.record_counts()
    assert counts == {"tags": 2, "overflow_records": 0, "markers": 0, "unknown": 1, "total": 3, "sync_tags": 2}
    assert [e.value for e in reader.events() if e.kind == "unknown-special"] == [1]
    assert sum(tags[c].size for c in tags.channels if c.input != 0) == 2


def test_t3_header_resolutions_are_checked(tmp_path):
    with pytest.raises(AdapterError, match="sync period"):
        PtuFile(write_ptu(tmp_path / "a.ptu", ptu.HYDRAHARP2_T3, [], sync_period_s=10.0))
    with pytest.raises(AdapterError, match="MeasDesc_Resolution"):
        PtuFile(write_ptu(tmp_path / "b.ptu", ptu.HYDRAHARP2_T3, [], resolution_s=1e-3))


@pytest.mark.parametrize("rtype", [ptu.PICOHARP_T3, ptu.HYDRAHARP1_T3, ptu.HYDRAHARP2_T3, ptu.MULTIHARP_T3])
def test_t3_times_agree_with_ptufile(tmp_path, rtype):
    """Independent cross-check, run only where the third-party `ptufile` package is installed.

    ptufile reports T3 `time` in sync periods and `dtime` in resolution units; detector
    records only (DADL's sync channel is derived, and channel numbers differ, as in T2).
    """
    ptufile = pytest.importorskip("ptufile")
    rng = np.random.default_rng(4)
    wrap = ptu.T3_FORMATS[rtype][1]
    res_ps = 4 if rtype == ptu.PICOHARP_T3 else 1
    ev, _ = _t3_events(rng, (1, 2), 300, 20 * wrap, PERIOD_PS // res_ps)
    path = write_ptu(tmp_path / "f.ptu", rtype, ev, resolution_s=res_ps * 1e-12)
    ours = read_ptu(path)
    with ptufile.PtuFile(path) as f:
        rec = f.decode_records()
    keep = rec["channel"] >= 0
    theirs = np.sort(rec["time"][keep].astype(np.int64) * PERIOD_PS + rec["dtime"][keep].astype(np.int64) * res_ps)
    mine = np.sort(np.concatenate([ours[c] for c in ours.channels if c.input != 0]))
    np.testing.assert_array_equal(mine, theirs)
