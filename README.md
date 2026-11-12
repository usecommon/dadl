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

## What works today (v0.0.3)

- A seeded time-tag simulator with injected efficiency, dark counts,
  non-paralyzable dead time, jitter, delay and cascading afterpulsing.
- **Calibration (phase 2, on `main`, unreleased):** dark rate, afterpulse
  probability and time constant, and dead time fitted from a dark acquisition,
  with 95% intervals whose coverage is checked in CI, written as
  content-addressed artifacts:
  `dadl simulate chsh-spad --dark --out dark.npz`, then
  `dadl calibrate dark chsh-spad --tags dark.npz --out artifacts/`.
- `setup.yaml` drives everything: copy it with `dadl init`, edit it, run it.
- Raw CHSH S from fixed-grid coincidences with a planned window, from
  simulated tags or **your own tags** (`dadl analyze`).
- Bell fringe and HOM dip visibility, with accidental subtraction only where
  the guardrail policy allows it.
- Guardrails that refuse invalid corrections, and record the policy
  fingerprint, overrides and window on every result.

Not yet: the other calibrations, calibration uncertainty in results, corrected CHSH, signing,
vendor adapters, the explain engine, MCP. Uncertainties today are statistical
only (Gaussian, 1σ).

## What DADL will not do

- Make low-cost detectors "research-grade". Lost efficiency and jitter cannot
  be recovered in software.
- Close the detection loophole in software. Efficiency correction is
  equivalent to assuming fair sampling.
- Correct data in a QKD security path. There, DADL only keeps monitoring
  records, which are not proof of security.
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
runs = dadl.read_tags("tags.npz", setup.channels)   # or build the dict yourself:
# runs[("a", "b")] = {"A+": np.int64 array of ps, "A-": ..., "B+": ..., "B-": ...}
result = dadl.chsh(runs, dadl.Window(width_ps=2000, source="plan"))
print(result)               # raw S, labels
print(result.as_record())   # flat record with labels, policy fingerprint, window

print(dadl.bell_visibility(1450, 550))            # raw fringe visibility
print(dadl.hom_visibility(1000, 120, accidentals=40))  # dip visibility, corrected value labelled
```

Vendor formats (PTU, Swabian, Photon-HDF5) come as adapters in phase 2.

## Developing DADL

From a clone of this repository:

```bash
pip install -e "./python[test]"
pytest python/tests
cargo test --workspace
```

The Rust core is pinned to the toolchain in `rust-toolchain.toml`; its
minimum supported Rust version is 1.74. The Python package is the reference
implementation in phase 1, and both are held to the shared cases in
`tests/vectors/core.json`. PyO3 bindings to the Rust core arrive in phase 2.

## Repository layout

| Path | Contents |
| --- | --- |
| `core/` | Rust core (`dadl` crate): dead time, fixed-grid coincidences, histograms |
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
