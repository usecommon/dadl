# Contributing to DADL

Thanks for helping. DADL is developed in public on GitHub.

## Ground rules

- **Raw is always preserved.** No API may return a corrected value without
  its raw counterpart and the provenance it used.
- **Guardrails are code, not advice.** Any change that loosens a guardrail
  must keep the regression tests in `python/tests/test_guardrails.py` passing
  and needs maintainer review. Once the project has two or more maintainers,
  that means two approvals. Until then, the change stays open as a public
  pull request for at least 7 days before merge, with the physics argument in
  its description.
- **Physics changes need a test.** Prefer a closed-form or published-result
  test over a tolerance tuned to make CI green.
- **Rust and Python agree.** A change to dead time, coincidences or
  histograms updates `tests/vectors/core.json`, which both suites run.
- **Docs say what works today.** Planned features are written as planned.

## Development setup

- Python 3.10 or later: `pip install -e "./python[test]"`, then
  `pytest python/tests`.
- Rust: the toolchain is pinned in `rust-toolchain.toml`;
  `cargo test --workspace` and `cargo clippy --workspace --all-targets -- -D warnings`.

CI runs both on Windows, Linux and macOS for every pull request, and builds
and tests the sdist and wheel.

Record user-visible changes in `CHANGELOG.md`.

## Gates

Two statistical gates check that 95% intervals cover the injected truth over
1,000 simulated replicates per regime (binomial 99.9% band). They are slow and
are not part of `pytest`:

- calibrations: `python -m dadl.calibration.coverage FAMILY --replicates 1000 --jobs N`
  (FAMILY is dark, pairs, dcr, deadtime or all; CI seed 20261111);
- corrected values: `python -m dadl.analysis.propagation hom|chsh --replicates 1000 --jobs N`
  (CI seed 20261118).

`--jobs` runs replicates in parallel processes; the result does not depend on
it. `--regimes` picks regimes. Slow regimes run in CI as shards
(`--shard k/n --out FILE`), summed and gated by `merge FILES`.

Run the gate of a calibration or analysis whose statistics you change, on the
CI seed, before opening the pull request. A change to a pinned calibration
module (`MODULES` in `python/tests/test_method_versions.py`), or to shared code
they call (`calibration/base.py`), also needs a new `METHOD_VERSION` for each
affected method, unless outputs are byte-identical: re-pin the hash and add `notes`. Then
`python python/tests/test_method_versions.py` regenerates the table that
`test_method_versions.py` checks. Never loosen a gate to make it pass.

Quick gate runs (100-200 replicates) judge corrected coverage only. The control rule of
chsh-dead-time needs about 400 replicates: its control covers about 87.5%, and at 100 replicates
it lands inside the band about one run in three. A FAIL on the control rule in a quick run is not
evidence of a problem; re-run at 400 or more.

The CI job `gates` passes only when every coverage and propagation job
passed. Once it is a required check on `main`, a failed gate blocks the merge.

## Licence of contributions

Unless you state otherwise, any contribution you submit is dual-licensed
under MIT OR Apache-2.0, without additional terms or conditions.
