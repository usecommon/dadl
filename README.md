# DADL

**Detector-Aware Data Layer for quantum optics experiments.**

DADL is an open-source starter kit for single-photon experiments. It is
designed to keep a versioned model of every detector channel (efficiency, dark
counts, dead time, afterpulsing, crosstalk, delay, jitter) and carry that
model's uncertainty into Bell S, visibility and fidelity, always reporting raw
and corrected values side by side and refusing corrections that would make a
result look better than the data supports.

**Read the overview:** [usecommon.github.io/dadl](https://usecommon.github.io/dadl/),
a short technical report on the problem DADL solves, how it is used in an
experiment, its guardrails and its roadmap.

> **Status: pre-alpha, phase 1 (partial), phase 2 in progress.** Anything may change in 0.0.x. The
> canonical phase plan is in [docs/design.md](docs/design.md#mvp-scope-and-milestones).

## What works today (on `main`; 0.0.3 is the latest release)

- A seeded time-tag simulator with injected efficiency, dark counts,
  non-paralyzable dead time, jitter, delay and cascading afterpulsing.
- **Calibration (phase 2, on `main`, unreleased):** dark rate, afterpulse
  probability and time constant, and dead time fitted from a dark acquisition,
  with 95% intervals whose coverage is checked in CI at a few simulated
  regimes (see "Known limits"), written as
  content-addressed artifacts:
  `dadl simulate chsh-spad --dark --out dark.npz`, then
  `dadl calibrate dark chsh-spad --tags dark.npz --out artifacts/`. From a
  pair-source acquisition, `dadl calibrate pairs` fits the relative delay and
  combined jitter, off-peak accidentals in the declared window and the
  heralding efficiency; `calibrate_dcr_temperature` fits DCR(T) (Arrhenius)
  from dark calibrations at several temperatures.
- `setup.yaml` drives everything: copy it with `dadl init`, edit it, run it.
- Raw CHSH S from fixed-grid coincidences with a planned window, from
  simulated tags or **your own tags** (`dadl analyze`).
- Bell fringe and HOM dip visibility, with accidental subtraction only where
  the guardrail policy allows it. `dadl quickstart hom-snspd` simulates a HOM
  dip on two SNSPDs and analyses it from the tags.
- **PicoQuant PTU files** (T2: PicoHarp, HydraHarp, TimeHarp 260, MultiHarp)
  read by DADL's own adapter (phase 2, on `main`).
- Guardrails that refuse invalid corrections, and record the policy
  fingerprint, overrides and window on every result.
- **Hardware twin** (phase 2, on `main`): `dadl simulate SETUP --out DIR
  --format ptu` writes the simulation as PicoQuant PTU files (marked
  simulated), and `SimulatedTagger` emulates a live tagger with buffer
  overflows, dropped batches and config changes. See `adapters/README.md`.
- **Optional compiled core** (phase 2, on `main`): with `dadl-quantum-core`
  installed, dead time, coincidences and histograms run in Rust (about 50x
  faster); `dadl.backend()` says which is in use.

Not yet: calibration uncertainty in results, corrected CHSH, signing,
the explain engine, MCP.

Uncertainties today are statistical only. Analyses (S, visibility) report a
Gaussian 1σ error, plus an exact 95% interval for visibility at low counts.
Calibrations report 95% intervals that need not be symmetric, validated by
simulation; each artifact's `uncertainty` block lists the components not
included (tagger timebase, temperature sensor, connectors, declared delays).

### Known limits

The calibrations were reviewed by outside-expert reviewers (AI agents, not
yet people; [register](docs/reviews/2026-11-expert-assumptions.md)). Their
limits, stated plainly:

- **CW sources only.** `source.type` must be `spdc-cw`. Pulsed sources and
  multi-pair emission are not modelled; heralding efficiency reads low when
  the mean pair number per slot is not small.
- **Single-exponential afterpulsing.** Detectors with a multi-component
  afterpulse kernel get the fast component's probability, not the total.
- **Hard-edge, non-paralyzable dead time.** Valid for actively quenched SPAD
  modules. The reported value is the effective channel dead time: the longer
  of detector and tagger dead time, and afterpulses inside it are not seen.
  Soft recovery (SNSPDs, passive quenching) is not modelled; the series
  dead-time model is planned.
- **Afterpulse time constant near the detection threshold** covers somewhat
  under 95% conditional on detection.
- **DCR(T)** is Arrhenius only, valid inside the calibrated temperature range
  at one bias, for SPADs.
- **Rate dependence.** Delay, jitter and efficiency are valid at the singles
  rates recorded in the artifact only.
- **Heralding efficiency** needs a heralding run with no analysers in the
  path (`--no-analysers`); on Bell data it measures the analyser projection.
- **Bell results** assume fair sampling and fair coincidence, with static
  settings: they do not close detection, locality or freedom-of-choice
  loopholes (nor the memory loophole, where a local model uses earlier
  trials). "Coincidences per singles" is the fraction of one arm's clicks
  that post-selection kept, about the other arm's efficiency times the window
  capture; closing the detection loophole needs about 0.83 (CHSH) or 2/3
  (Eberhard). CHSH needs both outputs of each analyser (A+, A−, B+, B−); the
  two-detector, 16-measurement scheme with rotating polarisers is not
  supported.
- **Dead time can raise S.** A recorded click makes the partner arm's dead
  state depend on the outcome, so |E| grows with η²Rτ: a V = 0.7 source
  (S = 1.98) read S = 2.14 with a 1 µs dead time at 3·10⁶ pairs/s. Results
  where coincidence rate × dead time exceeds 1% carry a "dead-time selection"
  label; tags whose gaps are impossible for the declared dead time (often a
  ns export read as ps) are refused.
- **Coverage gates** run 1,000 simulated replicates per calibration at a few
  regimes and catch gross failures (coverage below about 93%), not small ones.
- **Gaussian timing response.** The delay fit models a Gaussian IRF (with the
  tagger tick's quantization modelled exactly) on a two-level background. A detector with an exponential timing tail (most SNSPDs,
  many SPADs) is marked `misfit` and the delay artifact refused; the other
  pair calibrations are checked on their own.
- **SNSPDs.** `family` is a label only. Recovery, timing tails, bias-current
  dependence and latching signatures are not modelled; dark data cannot
  identify an SNSPD's dead time (reported as an upper limit, marked
  unidentified); DCR(T) refuses current-biased detectors. A latch or other
  interruption in any acquisition is detected (longest gap against the count
  rate) and marks calibrations `misfit` and analyses "interrupted".
- **Heralding efficiency is grid-specific.** It counts pairs in the same fixed
  slot, so it includes the window capture (about 1 − E|Δt|/W) and the signal
  arm's live fraction given a recorded herald. `calibrate pairs` prints the
  capture implied by the fitted jitter.
- **Off-peak checks have limited power.** A side-window excess smaller than
  about twice the accidentals interval can pass the check unseen.
- **Corrected visibility** is refused when the accidentals exceed the minimum
  count (or the corrected value would exceed 1), which removes the datasets
  nearest V = 1; report the raw value with its exact interval there.
- **HOM** is two-fold and two-point: the visibility at the declared dip delay
  against a reference delay, on a 50:50 beamsplitter, with detectors that do
  not resolve photon number. Four-fold heralded HOM, delay-scan fits and
  splitter imbalance (which caps V at 2RT/(R² + T²)) are not handled. The
  reference delay must be small against the window. Acquisitions of different
  length are compared as rates; singles cannot normalise drift (they fall at
  the dip by design).
- **Whole-file PTU reading.** A PTU file is decoded in memory at about 100
  bytes per record (215 with many markers): a 4 GB file needs about 100 GB.
  Chunked reading is planned. Unrecognised special records and a divided sync
  input mapped as a detector are refused.
- **Tagger settings** (input offsets, trigger levels, dead-time settings, the
  tick) are recorded in each calibration's `conditions` from PTU headers;
  they are not yet compared when the calibration is used.
- **Double clicks** in CHSH (both outputs of one side in a slot) count in two
  outcome combinations, which pulls E toward 0.
- **Reproducibility.** The same tags, setup and code give the same artifact
  ID when `created` is fixed (`--created`); posterior draws are seeded from
  the setup's `seed` mixed with the input hashes.
- **Records, not certificates.** An artifact names the tagger, input and
  detector serial, the input file and setup hashes, the acquisition time
  and the conditions you give (`--condition temperature_K=233`), but it is
  not a calibration record in the ISO/IEC 17025 sense: uncertainty is
  statistical only, and traceability, validity (`--valid-until`) and
  supersession (`--supersedes`) are recorded, not enforced. The ID gives
  integrity, not authenticity: signing is planned. Parameters reported only
  as a bound have `value: null` and `bound: "upper"`, and draw columns that
  are stand-ins rather than a posterior are listed in `posterior.stand_in`.
- A failed predictive check marks the calibration `misfit`; the CLI refuses to
  write it unless `--accept-misfit` is given.

## What DADL will not do

- Make low-cost detectors "research-grade". Lost efficiency and jitter cannot
  be recovered in software.
- Close the detection loophole in software. Efficiency correction is
  equivalent to assuming fair sampling.
- Correct data in a QKD security path, or detect attacks on detectors. DADL
  has no blinding, after-gate, time-shift or efficiency-mismatch detection
  and does not support gated detectors; nothing it produces belongs in a QKD
  security argument.
- Verify, in phase 1, that a window declared as planned really was fixed
  before acquisition, or that the tags you pass in are unfiltered. It records
  what you declare; signed analysis plans arrive in phase 3.

## Quickstart

Requires Python 3.10 or later. (On macOS the system `python3` is 3.9, and pip
then reports "No matching distribution found for dadl-quantum"; install a
newer Python first.) The PyPI package is `dadl-quantum`; the import name and
the command are both `dadl`.

```bash
pip install dadl-quantum
dadl quickstart
```

This loads the bundled `chsh-spad` setup: a polarization-entangled SPDC
source on four detectors with known imperfections. It prints raw S with a
fixed, pre-declared window next to the closed-form expectation.

To change the experiment, copy the setup, edit it and run your copy:

```bash
dadl init my-experiment
dadl quickstart my-experiment/setup.yaml
```

## Using your own data

Export your time tags as a table with three integer columns, in `.npz`
(arrays named `setting`, `channel`, `timestamp_ps`) or `.csv` (header
`setting,channel,timestamp_ps`):

| Column | Meaning |
| --- | --- |
| `setting` | CHSH setting pair, 0–3, in the order (a,b), (a,b′), (a′,b), (a′,b′) |
| `channel` | Detector input number, mapped to A+, A−, B+, B− by `channels:` in `setup.yaml` |
| `timestamp_ps` | int64 time tag in picoseconds |

Set the window width and per-channel delays in `setup.yaml` *before* looking
at the data, then:

```bash
dadl analyze my-experiment/setup.yaml --tags tags.npz
```

`dadl simulate my-experiment/setup.yaml --out tags.npz` writes simulated tags
in the same format, so you can check the round trip first.

From Python:

```python
import numpy as np
import dadl

setup = dadl.load_setup("my-experiment/setup.yaml")
runs = dadl.read_tags("tags.npz", setup.channels, setup.setting_keys)   # or build the dict yourself:
# runs[("a", "b")] = {"A+": np.int64 array of ps, "A-": ..., "B+": ..., "B-": ...}
result = dadl.chsh(runs, setup.window, settings_deg=setup.settings_deg)  # the planned window and delays
print(result)               # raw S, labels
print(result.as_record())   # flat record with labels, policy fingerprint, window

print(dadl.bell_visibility(1450, 550))            # raw fringe visibility
print(dadl.hom_visibility(1000, 120, accidentals=40))  # dip visibility, corrected value labelled
```

PicoQuant `.ptu` files (T2 mode) read directly: `from dadl.adapters import read_ptu`, then
`read_ptu("run.ptu")` returns int64 ps per channel, with rollovers and markers reported as events
(see `adapters/README.md`). From the command line, pass a folder of `setting0.ptu`,
`setting1.ptu`, ... as `--tags`. On HydraHarp, MultiHarp and TimeHarp 260 files, input 0 is the
sync input and detector inputs are numbered from 1, as in PicoQuant's software; map those
numbers under `channels:`. Other formats: export to the generic tag file.

Before calibrating recorded data, edit the copied `setup.yaml`:

```yaml
backend: file                 # not `sim`: calibrations of recorded data are refused under `sim`
hardware:
  tagger: HydraHarp-1234567   # the tagger's serial, recorded as the channel's device
  detectors: {A+: SPCM-0042}  # optional detector serials
```

The fixed coincidence grid starts at t = 0 of the tag file: slot k holds
times t with k·W ≤ t − delay < (k + 1)·W, and the record says so
(`grid_origin_ps`). One channel must have delay 0, so a common offset cannot
move the grid. Shifting all tags (for example re-exporting from the first
tag) still moves pairs across slot edges: each shift changes S within its
statistical error, but choosing the best of ten shifts doubles the chance of
a spurious S > 2 + 1.96σ. Keep the export that was planned.

## Developing DADL

From a clone of this repository:

```bash
pip install -e "./python[test]"
pytest python/tests
cargo test --workspace
pip install ./bindings/python    # optional: the compiled core (needs Rust)
```

The Rust core is pinned to the toolchain in `rust-toolchain.toml`; its
minimum supported Rust version is 1.74. The pure-Python package is the
reference implementation; `bindings/python` (PyO3, built with maturin) puts
the Rust core behind the same functions as the optional package
`dadl-quantum-core`. Both backends are held to the shared cases in
`tests/vectors/core.json`, and `python bindings/python/bench.py` compares
their throughput (about 50x for dead time, coincidences and histograms).

## Repository layout

| Path | Contents |
| --- | --- |
| `core/` | Rust core (`dadl` crate): dead time, fixed-grid coincidences, histograms |
| `bindings/python/` | PyO3 bindings: the optional `dadl-quantum-core` package (module `dadl_core`) |
| `python/` | Python package `dadl`: simulator, analyses, guardrails, config, tag files, CLI, tests |
| `tests/vectors/` | Test cases shared by the Rust and Python suites |
| `schema/` | JSON Schemas for calibration artifacts (`calibration/1`, draft) and `setup.yaml` (`setup-1`), with a validated example |
| `examples/` | Example setups, starting with `chsh-spad` (also bundled in the package) |
| `adapters/` | Vendor adapters (phase 2 onward) |
| `mcp/` | MCP agent surface (v0.2+, seeking contributors) |
| `docs/` | Design doc, review records and the GitHub Pages site (`docs/index.html`) |

## Citing

See [CITATION.cff](CITATION.cff). Please cite the version you used and keep
DADL's labels with any number you publish.

## Licence

Dual-licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT licence ([LICENSE-MIT](LICENSE-MIT))

at your option. Unless you state otherwise, any contribution you submit for
inclusion is dual-licensed as above, without additional terms.
