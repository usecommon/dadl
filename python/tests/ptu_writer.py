"""A minimal PTU writer (T2 and T3), used only to make test fixtures for dadl.adapters.ptu.

It follows PicoQuant's published file layout: the magic, a version, tagged
header entries ending with Header_End, then 32-bit records with overflow
records inserted whenever the time field (T2) or sync counter (T3) would wrap.
The T3 records are packed field by field as PicoQuant's demo code reads them
(PTU/C/ptudemo.cc, ProcessPHT3 and ProcessHHT3), independently of the reader.
"""

from __future__ import annotations

import struct
from pathlib import Path

from dadl.adapters import ptu

TY = {
    "int": ptu.TY_INT8,
    "float": ptu.TY_FLOAT8,
    "bool": ptu.TY_BOOL8,
    "str": ptu.TY_ANSI_STRING,
    "wstr": ptu.TY_WIDE_STRING,
    "empty": ptu.TY_EMPTY8,
}


def _entry(name: str, kind: str, value, idx: int = -1) -> bytes:
    ident = name.encode("ascii").ljust(32, b"\x00")
    typ = TY[kind]
    if kind == "str":
        blob = value.encode("latin-1") + b"\x00"
        blob = blob.ljust((len(blob) + 7) // 8 * 8, b"\x00")
        return struct.pack("<32siIq", ident, idx, typ, len(blob)) + blob
    if kind == "wstr":
        blob = (value + "\x00").encode("utf-16-le")
        blob = blob.ljust((len(blob) + 7) // 8 * 8, b"\x00")
        return struct.pack("<32siIq", ident, idx, typ, len(blob)) + blob
    if kind == "float":
        return struct.pack("<32siId", ident, idx, typ, value)
    if kind == "empty":
        return struct.pack("<32siIq", ident, idx, typ, 0)
    return struct.pack("<32siIq", ident, idx, typ, int(value))


def encode(record_type: int, events: list[tuple[int, str, int]], *, hh2_single_overflows: bool = False) -> list[int]:
    """T2 records for events (tick, "tag" | "marker", input or marker bits), in time order."""
    family, wrap = ptu.FORMATS[record_type]
    out: list[int] = []
    ofl = 0
    for tick, kind, value in sorted(events, key=lambda e: e[0]):
        k = (tick - ofl) // wrap
        if k:
            if family == "hh2" and not hh2_single_overflows:
                while k:
                    n = min(k, (1 << 25) - 1)
                    out.append((1 << 31) | (0x3F << 25) | n)
                    ofl += n * wrap
                    k -= n
            else:
                for _ in range(k):
                    out.append((0xF << 28) if family == "picoharp" else (1 << 31) | (0x3F << 25) | 1)
                    ofl += wrap
        t = tick - ofl
        if family == "picoharp":
            if kind == "tag":
                out.append((value << 28) | t)
            else:
                out.append((0xF << 28) | (t & ~0xF) | value)
        elif kind == "tag":
            out.append((1 << 31) | t if value == 0 else ((value - 1) << 25) | t)
        else:
            out.append((1 << 31) | (value << 25) | t)
    return out


def encode_t3(
    record_type: int, events: list[tuple[int, str, int, int]], *, hh2_single_overflows: bool = False
) -> list[int]:
    """T3 records for events (nsync, "tag" | "marker", input or marker bits, dtime), in nsync order.

    Inputs use the reader's numbering: PicoHarp channels as recorded (1-4 are valid);
    HydraHarp-family input k >= 1 is record channel k - 1. Markers ignore dtime.
    """
    family, wrap = ptu.T3_FORMATS[record_type]
    out: list[int] = []
    ofl = 0
    for nsync, kind, value, dtime in sorted(events, key=lambda e: e[0]):
        k = (nsync - ofl) // wrap
        if k:
            if family == "hh2" and not hh2_single_overflows:
                while k:
                    n = min(k, (1 << 10) - 1)
                    out.append((1 << 31) | (0x3F << 25) | n)
                    ofl += n * wrap
                    k -= n
            else:
                for _ in range(k):
                    # PicoHarp: channel 15 with dtime 0. HydraHarp v1 counts one wrap whatever nsync holds.
                    out.append((0xF << 28) if family == "picoharp" else (1 << 31) | (0x3F << 25) | 1)
                    ofl += wrap
        n = nsync - ofl
        if family == "picoharp":
            if kind == "tag":
                assert 0 <= dtime < (1 << 12)
                out.append((value << 28) | (dtime << 16) | n)
            else:
                out.append((0xF << 28) | (value << 16) | n)
        elif kind == "tag":
            assert value >= 1 and 0 <= dtime < (1 << 15)
            out.append(((value - 1) << 25) | (dtime << 10) | n)
        else:
            out.append((1 << 31) | (value << 25) | n)
    return out


def write_ptu(
    path: str | Path,
    record_type: int,
    events: list[tuple],
    *,
    resolution_s: float = 1e-12,
    sync_period_s: float = 12.5e-9,
    acquisition_ms: float | None = 1000.0,
    n_records: int | None = None,
    trailing: bytes = b"",
    extra: list[tuple[str, str, object]] = (),
    hh2_single_overflows: bool = False,
) -> Path:
    """A T2 fixture (events of 3) or, for a T3 record type, a T3 fixture (events of 4).

    In T3, `resolution_s` is the dtime unit and `sync_period_s` the sync period
    (MeasDesc_GlobalResolution); in T2 both header resolutions are `resolution_s`.
    """
    t3 = record_type in ptu.T3_TYPES
    if t3:
        records = encode_t3(record_type, events, hh2_single_overflows=hh2_single_overflows)
    else:
        records = encode(record_type, events, hh2_single_overflows=hh2_single_overflows)
    picoharp = record_type in (ptu.PICOHARP_T2, ptu.PICOHARP_T3)
    header = [
        _entry("File_GUID", "str", "{00000000-0000-0000-0000-000000000000}"),
        _entry("HW_Type", "str", "PicoHarp 300" if picoharp else "HydraHarp"),
        _entry("HW_SerialNo", "str", "1234567"),
        _entry("HW_InpChannels", "int", 5),
        _entry("HW_Markers", "int", 4),
        _entry("UsrHeadName", "wstr", "test fixture"),
        _entry("MeasDesc_GlobalResolution", "float", sync_period_s if t3 else resolution_s),
        _entry("MeasDesc_Resolution", "float", resolution_s),
        _entry("Measurement_Mode", "int", 3 if t3 else 2),
        _entry("Measurement_SubMode", "int", 0),
        _entry("TTResultFormat_TTTRRecType", "int", record_type),
        _entry("TTResultFormat_BitsPerRecord", "int", 32),
        _entry("TTResult_NumberOfRecords", "int", len(records) if n_records is None else n_records),
        _entry("TTResult_StopReason", "int", 0),
        _entry("HWSync_Offset", "empty", None),
        _entry("HWInpChan_TrgLevel", "int", -50, idx=0),
        _entry("HWInpChan_TrgLevel", "int", -60, idx=1),
    ]
    if acquisition_ms is not None:
        header.append(_entry("MeasDesc_AcquisitionTime", "float", acquisition_ms))
    header += [_entry(n, k, v) for n, k, v in extra]
    header.append(_entry("Header_End", "empty", None))
    body = struct.pack(f"<{len(records)}I", *records)
    path = Path(path)
    path.write_bytes(ptu.MAGIC + b"1.0.00\x00\x00" + b"".join(header) + body + trailing)
    return path
