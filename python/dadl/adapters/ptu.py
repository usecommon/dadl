"""PicoQuant PTU files, T2 mode, read without vendor software.

A PTU file is the 8-byte magic `PQTTTR\\0\\0`, an 8-byte version, a list of
tagged header entries ending with `Header_End`, then 32-bit records.
Supported record types (T2: every event carries its absolute time):

| Device | Record type | Tick wrap per overflow |
| --- | --- | --- |
| PicoHarp 300 | 0x00010203 | 210698240 |
| HydraHarp, format v1 | 0x00010204 | 33552000 |
| HydraHarp v2, TimeHarp 260 N/P, generic T2 (MultiHarp, PicoHarp 330) | 0x01010204, 0x00010205, 0x00010206, 0x00010207 | 33554432 x count |

T3 files (times relative to a sync pulse) are refused with a clear error.

Channels: PicoHarp inputs keep their numbers (0 and 1). On the HydraHarp
family the sync input is 0 and detector input k (numbered from 0 in the
record) is k + 1, matching PicoQuant's own numbering.

Times: ticks are multiplied by `MeasDesc_GlobalResolution` to give ps. When
that is a whole number of ps (1, 4, 5, 80 ps, ...) the conversion is exact
integer arithmetic; otherwise tags are rounded to the nearest ps and a note
says so.

The whole record body is decoded in memory: about 100 bytes per record at
peak (a 160 MB file took 4.2 GB; with markers about 215 bytes per record).
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

import struct
from collections.abc import Mapping
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
#: TTResult_StopReason values that end a complete acquisition; any other value means data were lost.
COMPLETE_STOP_REASONS = {0: "time over", 1: "manual stop", 2: "count overflow stop"}
STOP_REASON_NAMES = {3: "error", -2: "FIFO overrun", -9: "queue overrun", -10: "data transfer failure",
                     -12: "reference clock lost", -13: "sync lost"}
#: PicoHarp 300 T2 channel values: 0 is the sync input (Ch0), 1-4 are Ch1 through a PHR router
#: (channel 1 without one). PicoQuant's T2 demo code treats channel > 4 as illegal.
PICOHARP_MAX_INPUT = 4

#: Header bytes read first; a longer header is read on in steps.
HEADER_CHUNK = 1 << 16
#: More inputs than any PicoQuant device has; a larger HW_InpChannels is a corrupt header.
MAX_HEADER_INPUTS = 256

T3_TYPES = {0x00010303, 0x00010304, 0x01010304, 0x00010305, 0x00010306, 0x00010307}

PICOHARP_WRAP = 210_698_240
HYDRAHARP1_WRAP = 33_552_000
HH2_WRAP = 33_554_432

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


def _decode(records: np.ndarray, family: str, wrap: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Ticks and input numbers of every tag, plus the non-tag records."""
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
    # A crafted file of maximal overflow records wrapped the int64 sum silently (QKD review, cycle 4).
    if incr.size and float(np.sum(incr, dtype=np.float64)) + float(wrap) >= float(np.iinfo(np.int64).max):
        raise AdapterError("PTU overflow records exceed the int64 tick range")
    offset = np.cumsum(incr, dtype=np.int64)
    time = offset + tt
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


def _to_ps(ticks: np.ndarray, resolution_s: float) -> tuple[np.ndarray, bool]:
    res_ps = resolution_s * 1e12
    whole = round(res_ps)
    if whole >= 1 and abs(res_ps - whole) <= 1e-9 * res_ps:
        if ticks.size and int(ticks.max()) > np.iinfo(np.int64).max // whole:
            raise AdapterError("PTU times exceed the int64 ps range")
        return ticks * np.int64(whole), True
    return np.rint(ticks.astype(np.float64) * res_ps).astype(np.int64), False


class PtuFile:
    """Adapter for one PicoQuant PTU file in T2 mode."""

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
        if rtype in T3_TYPES:
            raise AdapterError(f"{self.path}: T3 mode (record type 0x{rtype:08X}) is not supported; record in T2 mode")
        if rtype not in FORMATS:
            raise AdapterError(f"{self.path}: unsupported PTU record type {rtype!r}")
        self.record_type = int(rtype)
        res = self.header.get("MeasDesc_GlobalResolution")
        if not isinstance(res, float) or not res > 0:
            raise AdapterError(f"{self.path}: MeasDesc_GlobalResolution is missing or not positive")
        self.resolution_s = res
        self._events: tuple[Event, ...] = ()
        self._counts: dict[str, int] = {}

    def describe(self) -> DeviceInfo:
        family = FORMATS[self.record_type][0]
        n_inputs = self.header.get("HW_InpChannels")
        if family == "picoharp":
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
        family, wrap = FORMATS[self.record_type]
        ticks, inputs, extra = _decode(records, family, wrap)
        ps, exact = _to_ps(ticks, self.resolution_s)
        marker_ps, _ = _to_ps(extra["marker_times"], self.resolution_s)

        if family == "picoharp" and inputs.size and int(inputs.max()) > PICOHARP_MAX_INPUT:
            raise AdapterError(
                f"{self.path}: PicoHarp T2 record on input {int(inputs.max())}; valid inputs are 0-{PICOHARP_MAX_INPUT}"
            )
        info = self.describe()
        device = info.model if info.serial is None else f"{info.model}-{info.serial}"
        tags: dict[ChannelId, np.ndarray] = {}
        for inp in np.unique(inputs):
            arr = ps[inputs == inp]
            # T2 records are written in time order; refuse rather than silently sort if not.
            if arr.size > 1 and np.any(arr[1:] < arr[:-1]):
                raise AdapterError(f"{self.path}: input {int(inp)} is not in time order")
            tags[ChannelId(device, int(inp))] = arr

        events = [Event("rollover", value=extra["rollovers"], detail={"records": extra["overflow_records"]})]
        events += [Event("marker", time_ps=int(t), value=int(b)) for t, b in zip(marker_ps, extra["marker_bits"])]
        if extra["unknown"]:
            events.append(Event("unknown-special", value=extra["unknown"]))
        stop = self.header.get("TTResult_StopReason")
        incomplete = False
        if stop is None:
            # QKD review, cycle 4: without a stop reason, completeness is unknown, not assured.
            events.append(Event("stop-reason", detail={"reason": "not recorded: completeness unknown"}))
        if stop is not None:
            stop = int(stop)
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
            "tags": int(ticks.size),
            "overflow_records": extra["overflow_records"],
            "markers": int(marker_ps.size),
            "unknown": extra["unknown"],
            "total": int(records.size),
        }

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
        if incomplete:
            metadata["data_loss"] = STOP_REASON_NAMES.get(stop, f"stop reason {stop}")
        if not exact:
            metadata["note"] = "the tick is not a whole number of ps; tags were rounded to the nearest ps"
        return TagSet(tags, duration_ps=duration, clock=f"{device} T2 clock", metadata=metadata)


def _str_or_none(v: Any) -> str | None:
    return None if v in (None, "") else str(v)


def read_ptu(path: str | Path) -> TagSet:
    """Read a PTU T2 file into a TagSet (see `PtuFile` for events and provenance)."""
    return PtuFile(path).read()
