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

## Licence of contributions

Unless you state otherwise, any contribution you submit is dual-licensed
under MIT OR Apache-2.0, without additional terms or conditions.
