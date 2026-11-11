# Contributing to DADL

Thanks for helping. DADL is developed in public on GitHub.

## Ground rules

- **Raw is always preserved.** No API may return a corrected value without
  its raw counterpart and the calibration artifacts it used.
- **Guardrails are code, not advice.** Any change that loosens a guardrail,
  the approval flow or signing needs approval from two maintainers and must
  keep the guardrail regression tests in `python/tests/test_guardrails.py`
  passing.
- **Physics changes need a test.** Prefer a closed-form or published-result
  test over a tolerance tuned to make CI green.

## Development setup

- Python 3.10 or later: `pip install -e "./python[test]"`, then
  `pytest python/tests`.
- Rust stable: `cargo test --workspace`.

CI runs both on Windows, Linux and macOS for every pull request.

## Licence of contributions

Unless you state otherwise, any contribution you submit is dual-licensed
under MIT OR Apache-2.0, without additional terms or conditions.
