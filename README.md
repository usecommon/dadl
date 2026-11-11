# DADL

**Detector-Aware Data Layer for quantum optics experiments.**

DADL is an open-source starter kit for single-photon experiments. It keeps a
versioned model of every detector channel (efficiency, dark counts, dead time,
afterpulsing, crosstalk, delay, jitter) and carries that model's uncertainty
into Bell S, visibility and fidelity. It always reports raw and corrected
values side by side, and it refuses corrections that would make a result look
better than the data supports.

**Read the overview:** [usecommon.github.io/dadl](https://usecommon.github.io/dadl/),
a short technical report on the problem DADL solves, how it is used in an
experiment, its guardrails and its roadmap.

> **Status: pre-alpha, phase 1.** The simulator, closed-form tests and the
> guardrail regression cases are in place. Calibration fits, vendor adapters,
> the explain engine and the MCP agent surface come in later phases. See
> [docs/design.md](docs/design.md).

## What DADL will not do

- Make low-cost detectors "research-grade". Lost efficiency and jitter cannot
  be recovered in software.
- Close the detection loophole in software. Efficiency correction is
  equivalent to assuming fair sampling.
- Correct data in a QKD security path. There, DADL only keeps monitoring
  records, which are not proof of security.

## Quickstart (phase 1)

Requires Python 3.10 or later. The PyPI package is `dadl-quantum`; the import
name and the command are both `dadl`.

```bash
pip install dadl-quantum
dadl quickstart chsh-spad
```

To work on DADL itself, from a clone of this repository:

```bash
pip install -e "./python[test]"
dadl quickstart chsh-spad
pytest python/tests
```

`dadl quickstart chsh-spad` simulates a polarization-entangled SPDC source on
four detectors with known imperfections, computes CHSH S with fixed,
pre-declared windows, and prints the raw value next to the closed-form
expectation.

The Rust core builds and tests separately:

```bash
cargo test --workspace
```

## Repository layout

| Path | Contents |
| --- | --- |
| `core/` | Rust core (`dadl` crate): per-tag loops such as dead time, fixed-grid coincidences, histograms |
| `python/` | Python package `dadl`: simulator, analyses, guardrails, CLI, tests |
| `schema/` | JSON Schema for calibration artifacts (`calibration/1`) |
| `examples/` | Example setups, starting with `chsh-spad` |
| `adapters/` | Vendor adapters (phase 2 onward) |
| `mcp/` | MCP agent surface (phase 4) |
| `docs/` | Design doc and the GitHub Pages site (`docs/index.html`) |

The Python package is pure NumPy in phase 1. PyO3 bindings to the Rust core
arrive in phase 2, built as abi3 wheels with maturin.

## Licence

Dual-licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT licence ([LICENSE-MIT](LICENSE-MIT))

at your option. Unless you state otherwise, any contribution you submit for
inclusion is dual-licensed as above, without additional terms.
