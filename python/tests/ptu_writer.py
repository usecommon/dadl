"""A minimal PTU T2 writer, used only to make test fixtures for dadl.adapters.ptu.

It follows PicoQuant's published file layout: the magic, a version, tagged
header entries ending with Header_End, then 32-bit T2 records with overflow
records inserted whenever the time field would wrap.
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
    """Records for events (tick, "tag" | "marker", input or marker bits), in time order."""
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


def write_ptu(
    path: str | Path,
    record_type: int,
    events: list[tuple[int, str, int]],
    *,
    resolution_s: float = 1e-12,
    acquisition_ms: float | None = 1000.0,
    n_records: int | None = None,
    trailing: bytes = b"",
    extra: list[tuple[str, str, object]] = (),
    hh2_single_overflows: bool = False,
) -> Path:
    records = encode(record_type, events, hh2_single_overflows=hh2_single_overflows)
    header = [
        _entry("File_GUID", "str", "{00000000-0000-0000-0000-000000000000}"),
        _entry("HW_Type", "str", "HydraHarp" if record_type != ptu.PICOHARP_T2 else "PicoHarp 300"),
        _entry("HW_SerialNo", "str", "1234567"),
        _entry("HW_InpChannels", "int", 5),
        _entry("HW_Markers", "int", 4),
        _entry("UsrHeadName", "wstr", "test fixture"),
        _entry("MeasDesc_GlobalResolution", "float", resolution_s),
        _entry("MeasDesc_Resolution", "float", resolution_s),
        _entry("Measurement_Mode", "int", 2),
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
