"""Write time tags as PicoQuant PTU files (T2 mode): the file half of DADL's hardware twin.

`write_ptu` turns tags (int64 ps per input) into the file a PicoQuant
instrument would have written: the PTU header, then 32-bit T2 records with
the device's own tick, timer width and rollover records. Simulated data then
enters DADL through the same adapter (`dadl.adapters.ptu`) as a recording.

Tags are quantized to the device tick (floor), as the instrument does; with
the HydraHarp's 1 ps tick, simulated ps tags pass through unchanged.

Every file this module writes says so: the header carries
`CreatorSW_Name = "DADL"` and `DADL_Simulated = True`, and the reader marks
such files `simulated` in their metadata, so simulated data cannot pass for
a recording.

The test suite keeps an independent, record-by-record writer
(python/tests/ptu_writer.py); the two must produce identical bytes.
"""

from __future__ import annotations

import struct
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from dadl import __version__
from dadl.adapters import ptu
from dadl.errors import AdapterError

#: Device presets: record type, tick (ps), HW_Type.
DEVICES: dict[str, tuple[int, int, str]] = {
    "picoharp": (ptu.PICOHARP_T2, 4, "PicoHarp 300"),
    "hydraharp": (ptu.HYDRAHARP2_T2, 1, "HydraHarp"),
    "timeharp260p": (ptu.TIMEHARP260P_T2, 25, "TimeHarp 260 P"),
    "multiharp": (ptu.MULTIHARP_T2, 5, "MultiHarp"),
}

_HH2_MAX_COUNT = (1 << 25) - 1


def _entry(name: str, typ: int, value: Any, idx: int = -1) -> bytes:
    ident = name.encode("ascii")
    if len(ident) > 31:
        raise ValueError(f"PTU tag name {name!r} is longer than 31 characters")
    ident = ident.ljust(32, b"\x00")
    if typ == ptu.TY_ANSI_STRING:
        blob = str(value).encode("latin-1") + b"\x00"
        blob = blob.ljust((len(blob) + 7) // 8 * 8, b"\x00")
        return struct.pack("<32siIq", ident, idx, typ, len(blob)) + blob
    if typ == ptu.TY_FLOAT8:
        return struct.pack("<32siId", ident, idx, typ, float(value))
    if typ == ptu.TY_EMPTY8:
        return struct.pack("<32siIq", ident, idx, typ, 0)
    if typ == ptu.TY_BOOL8:
        return struct.pack("<32siIq", ident, idx, typ, 1 if value else 0)
    return struct.pack("<32siIq", ident, idx, typ, int(value))


def encode_records(
    record_type: int,
    ticks: np.ndarray,
    inputs: np.ndarray,
    is_marker: np.ndarray | None = None,
) -> np.ndarray:
    """32-bit T2 records for events sorted by tick: tags (inputs) and markers (marker bits).

    Overflow records are inserted wherever the timer would wrap, exactly as
    the instrument does: one per wrap for PicoHarp and HydraHarp v1, one
    record carrying the wrap count (up to 2^25 - 1) for the HydraHarp v2 family.
    """
    family, wrap = ptu.FORMATS[record_type]
    ticks = np.asarray(ticks, dtype=np.int64)
    inputs = np.asarray(inputs, dtype=np.int64)
    marker = np.zeros(ticks.size, dtype=bool) if is_marker is None else np.asarray(is_marker, dtype=bool)
    if ticks.size and (np.any(np.diff(ticks) < 0) or ticks[0] < 0):
        raise ValueError("ticks must be non-negative and sorted")
    q = ticks // wrap
    delta = np.diff(q, prepend=0)
    t = (ticks - q * wrap).astype(np.uint32)

    if family == "picoharp":
        if np.any(inputs[~marker] > 0xE) or np.any((inputs[marker] < 1) | (inputs[marker] > 0xF)):
            raise ValueError("PicoHarp T2 inputs must be 0..14 and marker bits 1..15")
        event = np.where(
            marker,
            (np.uint32(0xF) << np.uint32(28)) | (t & np.uint32(0x0FFFFFF0)) | inputs.astype(np.uint32),
            (inputs.astype(np.uint32) << np.uint32(28)) | t,
        )
        n_ofl = delta
        ofl_value = np.uint32(0xF) << np.uint32(28)
    else:
        if np.any(inputs[~marker] > 63) or np.any((inputs[marker] < 1) | (inputs[marker] > 15)):
            raise ValueError("HydraHarp-family inputs must be 0..63 and marker bits 1..15")
        special = np.uint32(1) << np.uint32(31)
        chan = np.where(marker, inputs, np.where(inputs == 0, 0, inputs - 1)).astype(np.uint32)
        is_special = marker | (inputs == 0)
        event = np.where(is_special, special, np.uint32(0)) | (chan << np.uint32(25)) | t
        if family == "hh2":
            n_ofl = -(-delta // _HH2_MAX_COUNT)
        else:
            n_ofl = delta
        ofl_value = special | (np.uint32(0x3F) << np.uint32(25))

    total = int(ticks.size + n_ofl.sum())
    out = np.empty(total, dtype=np.uint32)
    # Each event i is preceded by n_ofl[i] overflow records.
    event_pos = np.arange(ticks.size) + np.cumsum(n_ofl)
    out[event_pos] = event
    is_ofl = np.ones(total, dtype=bool)
    is_ofl[event_pos] = False
    if family == "hh2":
        # Counts per overflow record: full chunks of _HH2_MAX_COUNT, then the remainder.
        idx = np.repeat(np.arange(ticks.size), n_ofl)
        first = np.repeat(np.cumsum(n_ofl) - n_ofl, n_ofl)
        k = np.arange(idx.size) - first  # chunk number within its event
        remaining = delta[idx] - k * _HH2_MAX_COUNT
        counts = np.minimum(remaining, _HH2_MAX_COUNT).astype(np.uint32)
        out[is_ofl] = ofl_value | counts
    elif family == "hh1":
        out[is_ofl] = ofl_value | np.uint32(1)
    else:
        out[is_ofl] = ofl_value
    return out


def write_ptu(
    path: str | Path,
    tags: Mapping[int, np.ndarray],
    *,
    device: str = "hydraharp",
    markers: Sequence[tuple[int, int]] = (),
    duration_ps: int | None = None,
    serial: str = "SIM-0001",
    trigger_levels_mv: Mapping[int, int] | None = None,
    comment: str = "DADL simulation, not a recording",
) -> Path:
    """Write `tags` ({input: sorted int64 ps}) and `markers` ([(time_ps, bits)]) as a T2 PTU file."""
    if device not in DEVICES:
        raise ValueError(f"device must be one of {sorted(DEVICES)}, got {device!r}")
    record_type, tick_ps, hw_type = DEVICES[device]
    parts_t, parts_i, parts_m = [], [], []
    for inp, arr in sorted(tags.items()):
        arr = np.asarray(arr, dtype=np.int64)
        if arr.size and arr[0] < 0:
            raise AdapterError("PTU times cannot be negative; shift the tags to start at 0 or later")
        parts_t.append(arr // tick_ps)
        parts_i.append(np.full(arr.size, int(inp), dtype=np.int64))
        parts_m.append(np.zeros(arr.size, dtype=bool))
    if markers:
        m = np.asarray(markers, dtype=np.int64).reshape(-1, 2)
        parts_t.append(m[:, 0] // tick_ps)
        parts_i.append(m[:, 1])
        parts_m.append(np.ones(len(m), dtype=bool))
    ticks = np.concatenate(parts_t) if parts_t else np.empty(0, dtype=np.int64)
    inputs = np.concatenate(parts_i) if parts_i else np.empty(0, dtype=np.int64)
    marker = np.concatenate(parts_m) if parts_m else np.empty(0, dtype=bool)
    order = np.argsort(ticks, kind="stable")
    records = encode_records(record_type, ticks[order], inputs[order], marker[order])

    n_inputs = max([int(k) for k in tags] + [0]) + 1
    if duration_ps is None:
        duration_ps = int(ticks.max() + 1) * tick_ps if ticks.size else 0
    header = [
        _entry("File_GUID", ptu.TY_ANSI_STRING, "{00000000-0000-0000-0000-000000000000}"),
        _entry("File_Comment", ptu.TY_ANSI_STRING, comment),
        _entry("CreatorSW_Name", ptu.TY_ANSI_STRING, "DADL"),
        _entry("CreatorSW_Version", ptu.TY_ANSI_STRING, __version__),
        _entry("DADL_Simulated", ptu.TY_BOOL8, True),
        _entry("HW_Type", ptu.TY_ANSI_STRING, hw_type),
        _entry("HW_SerialNo", ptu.TY_ANSI_STRING, serial),
        _entry("HW_InpChannels", ptu.TY_INT8, n_inputs),
        _entry("Measurement_Mode", ptu.TY_INT8, 2),
        _entry("Measurement_SubMode", ptu.TY_INT8, 0),
        _entry("MeasDesc_Resolution", ptu.TY_FLOAT8, tick_ps * 1e-12),
        _entry("MeasDesc_GlobalResolution", ptu.TY_FLOAT8, tick_ps * 1e-12),
        _entry("MeasDesc_AcquisitionTime", ptu.TY_FLOAT8, duration_ps / 1e9),
        _entry("TTResultFormat_TTTRRecType", ptu.TY_INT8, record_type),
        _entry("TTResultFormat_BitsPerRecord", ptu.TY_INT8, 32),
        _entry("TTResult_NumberOfRecords", ptu.TY_INT8, int(records.size)),
        _entry("TTResult_StopReason", ptu.TY_INT8, 0),
    ]
    for inp, level in sorted((trigger_levels_mv or {}).items()):
        header.append(_entry("HWInpChan_TrgLevel", ptu.TY_INT8, int(level), idx=int(inp)))
    header.append(_entry("Header_End", ptu.TY_EMPTY8, None))
    path = Path(path)
    path.write_bytes(ptu.MAGIC + b"1.0.00\x00\x00" + b"".join(header) + records.astype("<u4").tobytes())
    return path
