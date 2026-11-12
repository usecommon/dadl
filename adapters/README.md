# Adapters

An adapter turns a vendor's files or device stream into DADL `TagSet`s:
int64 picoseconds per channel on a declared clock, with every non-tag record
reported as a typed event.

| Adapter | Status | Module |
| --- | --- | --- |
| PicoQuant PTU, T2 mode (PicoHarp 300, HydraHarp v1/v2, TimeHarp 260 N/P, MultiHarp) | Available on `main` (phase 2) | `dadl.adapters.ptu` |
| PicoQuant T3 mode | Planned | |
| Swabian Instruments, qutools, IDQ | Seeking contributors | |

```python
from dadl.adapters import PtuFile

reader = PtuFile("run.ptu")
tags = reader.read()          # TagSet: {ChannelId(device, input): int64 ps}
reader.events()               # rollovers, markers, unknown records, stop reason
reader.config()               # scalar header entries, kept as provenance
```

The PTU reader is DADL's own (no vendor software). Its tests write PTU
fixtures with `python/tests/ptu_writer.py` and, when the third-party
`ptufile` package is installed, check that both readers give the same times.

Other formats: export to the generic `.npz` / `.csv` tag file that
`dadl analyze` reads (integer columns `setting`, `channel`, `timestamp_ps`).

Writing an adapter: see [TEMPLATE.md](TEMPLATE.md). Proprietary vendor SDKs
are never committed to this repository.

See `docs/design.md`, "Adapter contract", for the Rust trait this mirrors.
