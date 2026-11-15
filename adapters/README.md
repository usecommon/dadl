# Adapters

An adapter turns a vendor's files or device stream into DADL `TagSet`s:
int64 picoseconds per channel on a declared clock, with every non-tag record
reported as a typed event.

| Adapter | Status | Module |
| --- | --- | --- |
| PicoQuant PTU, T2 mode (PicoHarp 300, HydraHarp v1/v2, TimeHarp 260 N/P, MultiHarp) | Available (0.1.0) | `dadl.adapters.ptu` |
| PicoQuant PTU, T3 mode (the same devices) | Available (unreleased) | `dadl.adapters.ptu` |
| Simulated tagger (hardware twin): live batches with injectable buffer overflows, dropped batches and config changes | Available (0.1.0) | `dadl.adapters.simdevice` |
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

T3 files store each event's sync count (nsync) and start-stop time (dtime),
not an absolute time. The reader returns nsync x sync period + dtime x
resolution in ps, from the header's `MeasDesc_GlobalResolution` (the sync
period) and `MeasDesc_Resolution`: exact when both are whole ps, rounded
with a note otherwise, and only as good as the assumption that the sync
period held for the whole run. Input 0 holds the syncs that started a
recorded measurement, not every sync pulse (a `t3-sync` event says so).
`write_ptu(..., sync_period_ps=12_500)` writes simulated T3 files.

## Hardware twin

With no instrument at hand, DADL can stand in for one:

```bash
dadl simulate chsh-spad --out run/ --format ptu --device hydraharp   # setting0.ptu ... setting3.ptu
# devices: hydraharp, multiharp (5 ps), multiharp150n (80 ps), picoharp and timeharp260p
# (two inputs: the HOM example only)
dadl analyze chsh-spad --tags run/                                    # read through the PTU adapter
```

The files are real T2 PTU files (header, device tick, timer rollovers), so
simulated data enters through the same adapter as a recording, with known
ground truth. They carry `CreatorSW_Name = DADL` and `DADL_Simulated = True`;
DADL marks them simulated, and calibrations from them are recorded with
`backend: sim`.

`SimulatedTagger` emulates a live device: it streams batches and injects
faults, each reported as an event. `check_accounting` compares the delivered
tag times with the simulator's truth: every tag is either delivered or inside
a window reported lost, with no extras or duplicates. A real device has no
truth to compare with; there the rule is record accounting plus the vendor's
loss flags (see [TEMPLATE.md](TEMPLATE.md), rule 7).

PTU files stopped by a FIFO overrun, transfer failure or lost clock
(`TTResult_StopReason`) are marked incomplete and refused by the CLI, and the
duration is taken from the data, not the requested acquisition time.

```python
from dadl import load_setup
from dadl.adapters import Faults, SimulatedTagger, check_accounting

dev = SimulatedTagger(load_setup("setup.yaml"), faults=Faults(dropped_batches=(3,)))
for batch in dev.stream(): ...            # TagSet per 0.1 s batch
check_accounting(dev.read(), dev.events(), dev.truth())   # {input: 0, ...}
```

A twin reproduces only the behaviour we model. Instrument quirks still need
one real recording.

Other formats: export to the generic `.npz` / `.csv` tag file that
`dadl analyze` reads (integer columns `setting`, `channel`, `timestamp_ps`).

Writing an adapter: see [TEMPLATE.md](TEMPLATE.md). Proprietary vendor SDKs
are never committed to this repository.

See `docs/design.md`, "Adapter contract", for the Rust trait this mirrors.
