"""PicoQuant PTU files, T2 and T3 mode, read without vendor software.

A PTU file is the 8-byte magic `PQTTTR\\0\\0`, an 8-byte version, a list of
tagged header entries ending with `Header_End`, then 32-bit records.
Supported record types. T2, where every event carries its absolute time:

| Device | Record type | Tick wrap per overflow |
| --- | --- | --- |
| PicoHarp 300 | 0x00010203 | 210698240 |
| HydraHarp, format v1 | 0x00010204 | 33552000 |
| HydraHarp v2, TimeHarp 260 N/P, generic T2 (MultiHarp, PicoHarp 330) | 0x01010204, 0x00010205, 0x00010206, 0x00010207 | 33554432 x count |

T3, where every event carries the number of the sync period it fell in
(nsync) and its start-stop time after that sync (dtime):

| Device | Record type | Fields, most significant first (bits) | nsync wrap per overflow |
| --- | --- | --- | --- |
| PicoHarp 300 | 0x00010303 | channel 4, dtime 12, nsync 16 | 65536 |
| HydraHarp, format v1 | 0x00010304 | special 1, channel 6, dtime 15, nsync 10 | 1024 |
| HydraHarp v2, TimeHarp 260 N/P, generic T3 (MultiHarp, PicoHarp 330) | 0x01010304, 0x00010305, 0x00010306, 0x00010307 | as HydraHarp v1 | 1024 x count |

Record types, bit layouts and overflow rules follow PicoQuant's demo code
(github.com/PicoQuant/PicoQuant-Time-Tagged-File-Format-Demos: PTU/C/ptudemo.cc
and PTU/Python/Read_PTU.py). A T3 overflow record (PicoHarp: channel 15 with
dtime 0; HydraHarp family: special, channel 63) adds one wrap, or on v2-family
files `nsync` wraps (0 meaning 1). Markers are PicoHarp channel 15 with a
non-zero dtime (the dtime field holds the marker bits) and HydraHarp-family
special records on channels 1-15.

T3 times: a T3 file stores no absolute times. The absolute time of a T3 event
is nsync x sync period + dtime x resolution, with the sync period from
`MeasDesc_GlobalResolution` and the dtime unit from `MeasDesc_Resolution`.
This assumes the sync period in the header held for the whole acquisition: a
sync source that drifts by a fraction f moves a tag at time t by about f x t.
Markers carry no dtime and are placed at their sync. The sync input is a
channel of its own (input 0, as in T2), but a T3 file records only the syncs
that started a recorded start-stop measurement: the sync channel holds one
tag per distinct nsync of a tag record, not every sync pulse, a `t3-sync`
event says so, and `sync_periods` in the metadata counts the sync periods
the data span. Its singles rate is therefore not the sync rate.

Channels: PicoHarp inputs keep their numbers (T2: 0 and 1; T3 detector
records 1-4, with the sync as 0). On the HydraHarp family the sync input is 0
and detector input k (numbered from 0 in the record) is k + 1, matching
PicoQuant's own T2 numbering; T3 files use the same numbers, so an input
keeps its number in either mode.

T3 files come from a sync (usually pulsed) source. DADL models CW light: in a T3 file the sync is
used as a clock only. Corrected values are refused on T3 input unless the setup declares
`source.type: cw`. Pulsed sources are not supported on any input (T2, T3 or plain tags): pulsed
accidentals sit in side peaks at multiples of the repetition period. On T2 and plain tags DADL
detects them through the off-peak windows (not flat) and refuses the accidental-subtracted values;
a period that is a whole multiple of the window grid is not detected, but then the off-peak
estimate is also right.

T3 tick and lattice: a recorded gap between two T3 tags is within one dtime resolution of the
true gap, whatever the sync period, so the error bound of a T3 gap is the dtime resolution
(`tick_ps`; plus 1 ps when the tags were rounded to whole ps), not the common divisor of the tags.
It is used only for dead-time intervals. Histograms of pair delays use `lattice_ps`, the dtime
resolution itself (council round 8): within one sync the difference of two tags is an exact
multiple of it, and binning a 25 ps lattice in 26 ps bins biased the delay by -12 ps. A
resolution that is not a whole number of ps has no integer lattice (`lattice_ps` None). A T3 file whose dtime
range (2^bits x resolution) is shorter than the sync period records no photon late in each
period: a `t3-range` event and a note say so, and `dadl.tagfile.read_tags` refuses it. The sync
divider convention in T3 (`HWSync_Divider` > 1) is unverified; the metadata record it
(`sync_divider`).

Times: T2 ticks are multiplied by `MeasDesc_GlobalResolution` to give ps.
When that (in T3: both the sync period and the dtime resolution) is a whole
number of ps (1, 4, 5, 80 ps, a 12500 ps sync period, ...) the conversion is
exact integer arithmetic; otherwise tags are rounded to the nearest ps and a
note says so.

The whole record body is decoded in memory: about 100 bytes per record at
peak in T2 (a 160 MB file took 4.2 GB; with markers about 215 bytes per
record) and up to about 170 bytes per record in T3, which adds a sync tag per
occupied sync period.
Chunked reading for very large files is planned. Metadata checks read the
header only.

Data loss: PicoQuant reports how an acquisition ended in
`TTResult_StopReason`. Codes for lost data (3 error, -2 FIFO overrun, -9 queue
overrun, -10 transfer failure, -12 reference clock lost, -13 sync lost, and
any other negative code) become a `data-loss` event and mark the TagSet
`incomplete`; 0 (time over), 1 (manual stop) and 2 (count overflow stop) end a
complete acquisition. The duration comes from the data (last record time),
not from `MeasDesc_AcquisitionTime`, which is the requested time.

Every record is accounted for: tags, rollovers, markers and unrecognised
special records add up to `TTResult_NumberOfRecords`, and a file shorter or
longer than the header says is refused. Rollovers and unknown records are
reported as summary `Event`s, markers one by one.
"""

from __future__ import annotations

import math
import struct
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy as np

from dadl.adapters import Capabilities, DeviceInfo, Event
from dadl.data import ChannelId, TagSet
from dadl.errors import AdapterError

MAGIC = b"PQTTTR\x00\x00"
MAX_HEADER_ENTRIES = 100_000

TY_EMPTY8 = 0xFFFF0008
TY_BOOL8 = 0x00000008
TY_INT8 = 0x10000008
TY_BITSET64 = 0x11000008
TY_COLOR8 = 0x12000008
TY_FLOAT8 = 0x20000008
TY_TDATETIME = 0x21000008
TY_FLOAT8_ARRAY = 0x2001FFFF
TY_ANSI_STRING = 0x4001FFFF
TY_WIDE_STRING = 0x4002FFFF
TY_BINARY_BLOB = 0xFFFFFFFF

PICOHARP_T2 = 0x00010203
HYDRAHARP1_T2 = 0x00010204
HYDRAHARP2_T2 = 0x01010204
TIMEHARP260N_T2 = 0x00010205
TIMEHARP260P_T2 = 0x00010206
MULTIHARP_T2 = 0x00010207
PICOHARP_T3 = 0x00010303
HYDRAHARP1_T3 = 0x00010304
HYDRAHARP2_T3 = 0x01010304
TIMEHARP260N_T3 = 0x00010305
TIMEHARP260P_T3 = 0x00010306
MULTIHARP_T3 = 0x00010307
#: TTResult_StopReason values that end a complete acquisition; any other value means data were lost.
COMPLETE_STOP_REASONS = {0: "time over", 1: "manual stop", 2: "count overflow stop"}
STOP_REASON_NAMES = {3: "error", -2: "FIFO overrun", -9: "queue overrun", -10: "data transfer failure",
                     -12: "reference clock lost", -13: "sync lost"}
#: PicoHarp 300 T2 channel values: 0 is the sync input (Ch0), 1-4 are Ch1 through a PHR router
#: (channel 1 without one). PicoQuant's T2 demo code treats channel > 4 as illegal.
PICOHARP_MAX_INPUT = 4
#: PicoHarp 300 T3 detector channels are 1-4; PicoQuant's T3 demo code treats 0 and > 4 as illegal.
PICOHARP_T3_MIN_INPUT = 1
#: Sync periods (MeasDesc_GlobalResolution of a T3 file) accepted: 1 ps to 1 s.
T3_PERIOD_RANGE_S = (1e-12, 1.0)

#: Width of the T3 dtime field per family (bits).
T3_DTIME_BITS = {"picoharp": 12, "hh1": 15, "hh2": 15}

#: Header bytes read first; a longer header is read on in steps.
HEADER_CHUNK = 1 << 16
#: More inputs than any PicoQuant device has; a larger HW_InpChannels is a corrupt header.
MAX_HEADER_INPUTS = 256

PICOHARP_WRAP = 210_698_240
HYDRAHARP1_WRAP = 33_552_000
HH2_WRAP = 33_554_432
PICOHARP_T3_WRAP = 65_536
HH_T3_WRAP = 1_024

#: Record type -> (family, wrap). "hh1" counts one wrap per overflow record;
#: "hh2" counts the record's 25-bit value (0 meaning 1).
FORMATS = {
    PICOHARP_T2: ("picoharp", PICOHARP_WRAP),
    HYDRAHARP1_T2: ("hh1", HYDRAHARP1_WRAP),
    HYDRAHARP2_T2: ("hh2", HH2_WRAP),
    TIMEHARP260N_T2: ("hh2", HH2_WRAP),
    TIMEHARP260P_T2: ("hh2", HH2_WRAP),
    MULTIHARP_T2: ("hh2", HH2_WRAP),
}
#: T3 record type -> (family, nsync wrap): the same overflow rules per family as T2, with
#: "hh2" counting the record's 10-bit nsync value (0 meaning 1).
T3_FORMATS = {
    PICOHARP_T3: ("picoharp", PICOHARP_T3_WRAP),
    HYDRAHARP1_T3: ("hh1", HH_T3_WRAP),
    HYDRAHARP2_T3: ("hh2", HH_T3_WRAP),
    TIMEHARP260N_T3: ("hh2", HH_T3_WRAP),
    TIMEHARP260P_T3: ("hh2", HH_T3_WRAP),
    MULTIHARP_T3: ("hh2", HH_T3_WRAP),
}
T3_TYPES = frozenset(T3_FORMATS)

_INT64_MAX = int(np.iinfo(np.int64).max)


def parse_header(data: bytes) -> tuple[dict[str, Any], int]:
    """Header entries as {name: value} (indexed entries as "name(i)"), and the offset of the first record."""
    if data[:8] != MAGIC:
        raise AdapterError("not a PTU file (bad magic); only PicoQuant .ptu files are supported")
    tags: dict[str, Any] = {"__version__": data[8:16].rstrip(b"\x00").decode("ascii", "replace")}
    pos = 16
    for _ in range(MAX_HEADER_ENTRIES):
        if pos + 48 > len(data):
            raise AdapterError("PTU header is truncated (no Header_End)")
        ident, idx, typ = struct.unpack_from("<32siI", data, pos)
        raw = data[pos + 40 : pos + 48]
        pos += 48
        name = ident.split(b"\x00", 1)[0].decode("ascii", "replace")
        key = name if idx < 0 else f"{name}({idx})"
        if typ in (TY_ANSI_STRING, TY_WIDE_STRING, TY_FLOAT8_ARRAY, TY_BINARY_BLOB):
            (length,) = struct.unpack("<q", raw)
            if length < 0 or pos + length > len(data):
                raise AdapterError(f"PTU header entry {key} claims {length} bytes beyond the end of the file")
            blob = data[pos : pos + length]
            pos += length
            if typ == TY_ANSI_STRING:
                value: Any = blob.split(b"\x00", 1)[0].decode("latin-1")
            elif typ == TY_WIDE_STRING:
                value = blob.decode("utf-16-le", "replace").split("\x00", 1)[0]
            elif typ == TY_FLOAT8_ARRAY:
                value = np.frombuffer(blob[: len(blob) // 8 * 8], dtype="<f8").tolist()
            else:
                value = f"<{length} bytes>"
        elif typ in (TY_FLOAT8, TY_TDATETIME):
            (value,) = struct.unpack("<d", raw)
        elif typ == TY_BOOL8:
            value = struct.unpack("<q", raw)[0] != 0
        elif typ in (TY_INT8, TY_BITSET64, TY_COLOR8):
            (value,) = struct.unpack("<q", raw)
        elif typ == TY_EMPTY8:
            value = None
        else:
            raise AdapterError(f"PTU header entry {key} has unknown type 0x{typ:08X}")
        if name == "Header_End":
            return tags, pos
        tags[key] = value
    raise AdapterError("PTU header has too many entries")


def _unwrap(incr: np.ndarray, wrap: int) -> np.ndarray:
    """Running overflow offset; refuses sums beyond int64."""
    # A crafted file of maximal overflow records wrapped the int64 sum silently (QKD review, cycle 4).
    if incr.size and float(np.sum(incr, dtype=np.float64)) + float(wrap) >= float(_INT64_MAX):
        raise AdapterError("PTU overflow records exceed the int64 tick range")
    return np.cumsum(incr, dtype=np.int64)


def _decode(records: np.ndarray, family: str, wrap: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """T2: ticks and input numbers of every tag, plus the non-tag records."""
    rec = records.astype(np.uint64)
    if family == "picoharp":
        ch = (rec >> np.uint64(28)).astype(np.int64)
        tt = (rec & np.uint64(0x0FFFFFFF)).astype(np.int64)
        special = ch == 0xF
        low = tt & 0xF
        overflow = special & (low == 0)
        marker = special & (low != 0)
        unknown = np.zeros(rec.size, dtype=bool)
        incr = np.where(overflow, wrap, 0)
        is_tag = ~special
        tag_input = ch
        marker_bits = low
    else:
        sp = (rec >> np.uint64(31)).astype(bool)
        ch = ((rec >> np.uint64(25)) & np.uint64(0x3F)).astype(np.int64)
        tt = (rec & np.uint64(0x1FFFFFF)).astype(np.int64)
        overflow = sp & (ch == 0x3F)
        sync = sp & (ch == 0)
        marker = sp & (ch >= 1) & (ch <= 15)
        unknown = sp & ~overflow & ~sync & ~marker
        if family == "hh1":
            incr = np.where(overflow, wrap, 0)
        else:
            incr = np.where(overflow, wrap * np.where(tt == 0, 1, tt), 0)
        is_tag = ~sp | sync
        tag_input = np.where(sync, 0, ch + 1)
        marker_bits = ch
    time = _unwrap(incr, wrap) + tt
    extra = {
        "rollovers": int(np.sum(incr) // wrap) if incr.size else 0,
        "overflow_records": int(overflow.sum()),
        "marker_times": time[marker],
        "marker_bits": marker_bits[marker],
        "unknown": int(unknown.sum()),
    }
    if int(is_tag.sum() + overflow.sum() + marker.sum() + unknown.sum()) != rec.size:
        raise AssertionError("PTU record accounting does not add up")  # pragma: no cover
    return time[is_tag], tag_input[is_tag], extra


def _decode_t3(
    records: np.ndarray, family: str, wrap: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """T3: unwrapped nsync, dtime and input number of every tag, plus the non-tag records."""
    rec = records.astype(np.uint64)
    if family == "picoharp":
        ch = (rec >> np.uint64(28)).astype(np.int64)
        dtime = ((rec >> np.uint64(16)) & np.uint64(0xFFF)).astype(np.int64)
        nsync = (rec & np.uint64(0xFFFF)).astype(np.int64)
        special = ch == 0xF
        overflow = special & (dtime == 0)
        marker = special & (dtime != 0)
        unknown = np.zeros(rec.size, dtype=bool)
        incr = np.where(overflow, wrap, 0)
        is_tag = ~special
        tag_input = ch
        marker_bits = dtime
    else:
        sp = (rec >> np.uint64(31)).astype(bool)
        ch = ((rec >> np.uint64(25)) & np.uint64(0x3F)).astype(np.int64)
        dtime = ((rec >> np.uint64(10)) & np.uint64(0x7FFF)).astype(np.int64)
        nsync = (rec & np.uint64(0x3FF)).astype(np.int64)
        overflow = sp & (ch == 0x3F)
        marker = sp & (ch >= 1) & (ch <= 15)
        # Special records on channel 0 or 16-62 have no meaning in T3; they are counted, not read.
        unknown = sp & ~overflow & ~marker
        if family == "hh1":
            incr = np.where(overflow, wrap, 0)
        else:
            incr = np.where(overflow, wrap * np.where(nsync == 0, 1, nsync), 0)
        is_tag = ~sp
        tag_input = ch + 1
        marker_bits = ch
    n = _unwrap(incr, wrap) + nsync
    used = is_tag | marker
    extra = {
        "rollovers": int(np.sum(incr) // wrap) if incr.size else 0,
        "overflow_records": int(overflow.sum()),
        "marker_nsync": n[marker],
        "marker_bits": marker_bits[marker],
        "unknown": int(unknown.sum()),
        "sync_periods": int(n[used].max()) + 1 if used.any() else 0,
    }
    if int(is_tag.sum() + overflow.sum() + marker.sum() + unknown.sum()) != rec.size:
        raise AssertionError("PTU record accounting does not add up")  # pragma: no cover
    return n[is_tag], dtime[is_tag], tag_input[is_tag], extra


def _whole_ps(seconds: float) -> int | None:
    """`seconds` as a whole number of ps, or None when it is not one (to 1e-9 relative)."""
    res_ps = seconds * 1e12
    whole = round(res_ps)
    if whole >= 1 and abs(res_ps - whole) <= 1e-9 * res_ps:
        return int(whole)
    return None


def _to_ps(ticks: np.ndarray, resolution_s: float) -> tuple[np.ndarray, bool]:
    whole = _whole_ps(resolution_s)
    if whole is not None:
        if ticks.size and int(ticks.max()) > _INT64_MAX // whole:
            raise AdapterError("PTU times exceed the int64 ps range")
        return ticks * np.int64(whole), True
    return np.rint(ticks.astype(np.float64) * (resolution_s * 1e12)).astype(np.int64), False


def _t3_to_ps(nsync: np.ndarray, dtime: np.ndarray, period_s: float, resolution_s: float) -> tuple[np.ndarray, bool]:
    """nsync x period + dtime x resolution in ps: exact when both are whole ps, else rounded."""
    period, res = _whole_ps(period_s), _whole_ps(resolution_s)
    if period is not None and res is not None:
        if nsync.size and int(nsync.max()) * period + int(dtime.max()) * res > _INT64_MAX:
            raise AdapterError("PTU times exceed the int64 ps range")
        return nsync * np.int64(period) + dtime * np.int64(res), True
    period_ps, res_ps = period_s * 1e12, resolution_s * 1e12
    if nsync.size and float(nsync.max()) * period_ps + float(dtime.max()) * res_ps >= 0.99 * _INT64_MAX:
        raise AdapterError("PTU times exceed the int64 ps range")
    # The header's period exactly (the float64 product period_s * 1e12 is off by up to 1e-12 ps,
    # several ps after 1e12 syncs), split into whole ps (integer arithmetic) and a fraction.
    exact_period = Fraction(period_s) * 10**12
    base = math.floor(exact_period)
    frac = float(exact_period - base)
    rest = np.rint(nsync.astype(np.float64) * frac + dtime.astype(np.float64) * res_ps).astype(np.int64)
    return nsync * np.int64(base) + rest, False


def _resolution(header: Mapping[str, Any], key: str, path: Path, lo: float, hi: float, what: str) -> float:
    res = header.get(key)
    if not isinstance(res, float) or not lo <= res <= hi:
        raise AdapterError(f"{path}: {key} is missing or outside {what}")
    return res


class PtuFile:
    """Adapter for one PicoQuant PTU file in T2 or T3 mode."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        # Only the header is read here; `read` loads the records. Metadata checks (simulated
        # flag, tick, serial) used to read whole files several times each (tagger engineer
        # review, cycle 2).
        try:
            with self.path.open("rb") as fh:
                head = fh.read(HEADER_CHUNK)
                while True:
                    try:
                        self.header, self._offset = parse_header(head)
                        break
                    except AdapterError as exc:
                        more = fh.read(len(head))
                        if not more or "truncated" not in str(exc) and "beyond the end" not in str(exc):
                            raise
                        head += more
        except OSError as exc:
            raise AdapterError(f"cannot read {self.path}: {exc}") from exc
        rtype = self.header.get("TTResultFormat_TTTRRecType")
        if not isinstance(rtype, int) or isinstance(rtype, bool):
            raise AdapterError(f"{self.path}: TTResultFormat_TTTRRecType must be an integer, got {rtype!r}")
        if rtype in FORMATS:
            self.mode = "T2"
        elif rtype in T3_FORMATS:
            self.mode = "T3"
        else:
            raise AdapterError(f"{self.path}: unsupported PTU record type {rtype!r}")
        self.record_type = int(rtype)
        #: The time unit of a tag: the T2 tick, or the T3 dtime resolution.
        self.resolution_s: float
        #: T3 only: the sync period (MeasDesc_GlobalResolution in a T3 file); None in T2.
        self.sync_period_s: float | None = None
        if self.mode == "T2":
            self.resolution_s = _resolution(self.header, "MeasDesc_GlobalResolution", self.path, 1e-15, 1e-6,
                                            "1 fs - 1 us")
        else:
            self.resolution_s = _resolution(self.header, "MeasDesc_Resolution", self.path, 1e-15, 1e-6,
                                            "1 fs - 1 us")
            self.sync_period_s = _resolution(self.header, "MeasDesc_GlobalResolution", self.path,
                                             *T3_PERIOD_RANGE_S, "1 ps - 1 s (the T3 sync period)")
        self._events: tuple[Event, ...] = ()
        self._counts: dict[str, int] = {}

    @property
    def tick_ps(self) -> int | None:
        """The tick of the tags in whole ps (council round 7, G1). T2: `MeasDesc_GlobalResolution`
        when it is a whole number of ps, else None. T3: the dtime resolution rounded up to whole
        ps, plus 1 when the tags are rounded to whole ps (each by up to 0.5 ps)."""
        if self.mode == "T2":
            return _whole_ps(self.resolution_s)
        assert self.sync_period_s is not None
        res = _whole_ps(self.resolution_s)
        exact = res is not None and _whole_ps(self.sync_period_s) is not None
        tick = res if res is not None else math.ceil(self.resolution_s * 1e12)
        return int(tick) + (0 if exact else 1)

    @property
    def lattice_ps(self) -> int | None:
        """The lattice of pair delays in whole ps, for histograms (council round 8). T2:
        `MeasDesc_GlobalResolution`; T3: the dtime resolution. None when it is not a whole number
        of ps. Unlike `tick_ps` it has no +1 for rounding: within one sync two T3 tags share the
        rounded sync part, so their difference is an exact multiple of the resolution."""
        return _whole_ps(self.resolution_s)

    @property
    def sync_divider(self) -> int | None:
        """T3 only: `HWSync_Divider` when above 1, else None."""
        divider = self.header.get("HWSync_Divider")
        if self.mode == "T3" and isinstance(divider, int) and not isinstance(divider, bool) and divider > 1:
            return int(divider)
        return None

    def dtime_range_short(self) -> str | None:
        """T3 only: a note when the dtime range (2^bits x resolution) is shorter than the sync
        period, else None (council round 7, G4)."""
        if self.mode != "T3":
            return None
        assert self.sync_period_s is not None
        span = (1 << T3_DTIME_BITS[self.family]) * self.resolution_s
        if span >= self.sync_period_s * (1 - 1e-9):
            return None
        return (f"dtime range {span * 1e9:.4g} ns < sync period {self.sync_period_s * 1e9:.4g} ns: photons late in "
                "each period are not recorded")

    @property
    def family(self) -> str:
        return (FORMATS if self.mode == "T2" else T3_FORMATS)[self.record_type][0]

    def describe(self) -> DeviceInfo:
        n_inputs = self.header.get("HW_InpChannels")
        if self.family == "picoharp":
            channels = tuple(range(PICOHARP_MAX_INPUT + 1))
        elif isinstance(n_inputs, int) and 0 < n_inputs <= MAX_HEADER_INPUTS:
            # As the header reports it; check against the inputs that carry tags after `read`.
            channels = tuple(range(n_inputs))
        else:
            channels = ()
        return DeviceInfo(
            model=str(self.header.get("HW_Type", "PicoQuant")),
            serial=_str_or_none(self.header.get("HW_SerialNo")),
            firmware=_str_or_none(self.header.get("HW_Version")),
            channels=channels,
            resolution_ps=self.resolution_s * 1e12,
        )

    def capabilities(self) -> Capabilities:
        return Capabilities(live=False, replay=True)

    def config(self) -> Mapping[str, Any]:
        """Scalar header entries, as provenance (strings, numbers and booleans)."""
        return {k: v for k, v in self.header.items() if isinstance(v, (str, int, float, bool))}

    def events(self) -> tuple[Event, ...]:
        return self._events

    def record_counts(self) -> dict[str, int]:
        """How every record was used (filled by `read`)."""
        return dict(self._counts)

    def _read_t2(self, records: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool, dict[str, Any]]:
        family, wrap = FORMATS[self.record_type]
        ticks, inputs, extra = _decode(records, family, wrap)
        ps, exact = _to_ps(ticks, self.resolution_s)
        marker_ps, _ = _to_ps(extra["marker_times"], self.resolution_s)
        if family == "picoharp" and inputs.size and int(inputs.max()) > PICOHARP_MAX_INPUT:
            raise AdapterError(
                f"{self.path}: PicoHarp T2 record on input {int(inputs.max())}; valid inputs are 0-{PICOHARP_MAX_INPUT}"
            )
        extra["tag_records"] = int(ticks.size)
        return ps, inputs, marker_ps, exact, extra

    def _read_t3(self, records: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool, dict[str, Any]]:
        family, wrap = T3_FORMATS[self.record_type]
        assert self.sync_period_s is not None
        nsync, dtime, inputs, extra = _decode_t3(records, family, wrap)
        if family == "picoharp" and inputs.size:
            lo, hi = int(inputs.min()), int(inputs.max())
            if lo < PICOHARP_T3_MIN_INPUT or hi > PICOHARP_MAX_INPUT:
                bad = lo if lo < PICOHARP_T3_MIN_INPUT else hi
                raise AdapterError(
                    f"{self.path}: PicoHarp T3 record on input {bad}; valid inputs are "
                    f"{PICOHARP_T3_MIN_INPUT}-{PICOHARP_MAX_INPUT}"
                )
        ps, exact = _t3_to_ps(nsync, dtime, self.sync_period_s, self.resolution_s)
        markers = extra["marker_nsync"]
        marker_ps, _ = _t3_to_ps(markers, np.zeros_like(markers), self.sync_period_s, self.resolution_s)
        # The syncs that started a recorded measurement are the sync channel (input 0), as in T2.
        syncs = np.unique(nsync)
        sync_ps, _ = _t3_to_ps(syncs, np.zeros_like(syncs), self.sync_period_s, self.resolution_s)
        extra["tag_records"] = int(nsync.size)
        extra["sync_tags"] = int(syncs.size)
        ps = np.concatenate([ps, sync_ps])
        inputs = np.concatenate([inputs, np.zeros(syncs.size, dtype=np.int64)])
        return ps, inputs, marker_ps, exact, extra

    def read(self) -> TagSet:
        n_records = self.header.get("TTResult_NumberOfRecords")
        if not isinstance(n_records, int) or n_records < 0:
            raise AdapterError(f"{self.path}: TTResult_NumberOfRecords is missing")
        body_bytes = self.path.stat().st_size - self._offset
        if body_bytes != 4 * n_records:
            raise AdapterError(
                f"{self.path}: header says {n_records} records ({4 * n_records} bytes) but the file holds "
                f"{body_bytes} bytes of records; the file is truncated or has trailing data"
            )
        records = np.fromfile(self.path, dtype="<u4", count=n_records, offset=self._offset)
        reader = self._read_t2 if self.mode == "T2" else self._read_t3
        ps, inputs, marker_ps, exact, extra = reader(records)
        family = self.family

        info = self.describe()
        device = info.model if info.serial is None else f"{info.model}-{info.serial}"
        tags: dict[ChannelId, np.ndarray] = {}
        for inp in np.unique(inputs):
            arr = ps[inputs == inp]
            # Records are written in time order; refuse rather than silently sort if not. In T3
            # this also catches a dtime beyond the header's sync period.
            if arr.size > 1 and np.any(arr[1:] < arr[:-1]):
                hint = " (in T3: a dtime longer than the sync period in the header?)" if self.mode == "T3" else ""
                raise AdapterError(f"{self.path}: input {int(inp)} is not in time order{hint}")
            tags[ChannelId(device, int(inp))] = arr

        events = [Event("rollover", value=extra["rollovers"], detail={"records": extra["overflow_records"]})]
        events += [Event("marker", time_ps=int(t), value=int(b)) for t, b in zip(marker_ps, extra["marker_bits"])]
        if extra["unknown"]:
            events.append(Event("unknown-special", value=extra["unknown"]))
        if self.mode == "T3":
            events.append(Event("t3-sync", value=extra["sync_tags"], detail={
                "note": "input 0 holds only the syncs that started a recorded measurement, not every sync",
                "sync_periods": extra["sync_periods"],
            }))
        short_range = self.dtime_range_short()
        if short_range:
            events.append(Event("t3-range", detail={"note": short_range}))
        stop = self.header.get("TTResult_StopReason")
        incomplete = False
        if stop is None:
            # QKD review, cycle 4: without a stop reason, completeness is unknown, not assured.
            events.append(Event("stop-reason", detail={"reason": "not recorded: completeness unknown"}))
        if stop is not None:
            if not isinstance(stop, int) or isinstance(stop, bool):
                raise AdapterError(f"{self.path}: TTResult_StopReason must be an integer, got {stop!r}")
            if stop in COMPLETE_STOP_REASONS:
                events.append(Event("stop-reason", value=stop, detail={"reason": COMPLETE_STOP_REASONS[stop]}))
            else:
                incomplete = True
                events.append(Event("data-loss", value=stop, detail={
                    "reason": STOP_REASON_NAMES.get(stop, "unrecognised stop reason"),
                    "source": "TTResult_StopReason",
                }))
        divider = self.header.get("HWSync_Divider")
        if isinstance(divider, int) and divider > 1 and any(ch.input == 0 for ch in tags) and family != "picoharp":
            events.append(Event("sync-divider", value=divider,
                                detail={"note": "input 0 (sync) is divided: it records one event in every " f"{divider}"}))
        self._events = tuple(events)
        self._counts = {
            "tags": extra["tag_records"],
            "overflow_records": extra["overflow_records"],
            "markers": int(marker_ps.size),
            "unknown": extra["unknown"],
            "total": int(records.size),
        }
        if self.mode == "T3":
            # Derived from the tag records, not records themselves: outside the accounting sum.
            self._counts["sync_tags"] = extra["sync_tags"]

        # Duration from the data: the requested acquisition time is not the time acquired.
        last = max([int(a[-1]) for a in tags.values() if a.size] + [int(t) for t in marker_ps], default=-1)
        duration = last + max(1, int(round(info.resolution_ps))) if last >= 0 else None
        metadata: dict[str, Any] = {
            "source": "ptu",
            "file": self.path.name,
            "record_type": f"0x{self.record_type:08X}",
            "resolution_ps": info.resolution_ps,
            "exact_ps": exact,
            # Written by DADL's hardware twin (dadl.adapters.ptu_write), not by an instrument.
            "simulated": self.header.get("DADL_Simulated") is True or self.header.get("CreatorSW_Name") == "DADL",
            "incomplete": incomplete,
            "inputs_with_tags": sorted(ch.input for ch in tags),
        }
        if self.mode == "T3":
            assert self.sync_period_s is not None
            metadata.update({
                "mode": "T3",
                "sync_period_ps": self.sync_period_s * 1e12,
                "sync_periods": extra["sync_periods"],
                "time_model": "nsync x sync period + dtime x resolution; the header's sync period is assumed "
                              "constant over the acquisition",
                "sync_tags": "input 0 holds only the syncs that started a recorded measurement",
            })
            if self.sync_divider is not None:
                metadata["sync_divider"] = self.sync_divider
        metadata["tick_ps"] = self.tick_ps
        metadata["lattice_ps"] = self.lattice_ps
        if incomplete:
            metadata["data_loss"] = STOP_REASON_NAMES.get(stop, f"stop reason {stop}")
        if not exact:
            what = "the tick is" if self.mode == "T2" else "the sync period or dtime resolution is"
            metadata["note"] = f"{what} not a whole number of ps; tags were rounded to the nearest ps"
        if short_range:
            metadata["note"] = "; ".join(n for n in (metadata.get("note"), short_range) if n)
        return TagSet(tags, duration_ps=duration, clock=f"{device} {self.mode} clock", metadata=metadata)


def _str_or_none(v: Any) -> str | None:
    return None if v in (None, "") else str(v)


def read_ptu(path: str | Path) -> TagSet:
    """Read a PTU file (T2 or T3) into a TagSet (see `PtuFile` for events and provenance)."""
    return PtuFile(path).read()
