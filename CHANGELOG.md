# Changelog

All notable changes. DADL is pre-alpha: in 0.0.x anything may change between
releases. The `calibration/1` schema is a draft and is versioned separately
from the packages.

## Unreleased

### Phase 2, milestone 1: data model and the dark calibration

- `dadl calibrate dark` and `dadl.calibration.calibrate_dark`: dark rate,
  afterpulse probability and time constant, and dead time, fitted jointly from
  one dark acquisition's inter-arrival gaps (maximum likelihood, Laplace
  posterior; an upper limit on afterpulsing when none is detected), with a
  posterior predictive check.
- `dadl.artifact`: content-addressed calibration artifacts (RFC 8785 canonical
  manifest, SHA-256 ID, `.npy` posterior table), verified on read. JCS test
  vectors in `tests/vectors/jcs.json`.
- `dadl.data` (`ChannelId`, `TagSet`, `Dataset`), `dadl.models`
  (`DetectorModel`), `dadl.stats` (Wilson, Clopper–Pearson, Beta and Gamma
  posteriors, Laplace).
- `schema/setup-1.schema.json` for `setup.yaml`. `calibration/1` (still a
  draft) gains `summary`, `notes` and `unidentified_params`, and the `npy`
  posterior format.
- `dadl simulate --dark` writes dark acquisitions for calibration.
- Simulator: afterpulsing is self-exciting (afterpulses cascade and every
  avalanche starts dead time), matching the closed form
  `closed_form.dark_afterpulse_rate`. Simulated tags for a given seed differ
  from 0.0.3; the quickstart prints S = 2.6941 instead of 2.6890.
- CI: a coverage job runs the 1,000-replicate gate on every pull request.
- New dependency: SciPy (>= 1.10).

### Phase 2, milestone 2: delay, jitter, accidentals, heralding, DCR(T)

- `dadl.histogram.cross_correlation` and `dadl::histogram::cross_correlation`
  (Rust): histogram of delays b - a within a range, held to new shared vectors.
- `calibrate_delay` (kind `delay-jitter`): relative delay and combined jitter
  from a Gaussian peak plus flat background fitted to the cross-correlation
  (Poisson maximum likelihood, Laplace posterior, grouped chi-square check).
- `calibrate_accidentals` (kind `accidentals`, new in `calibration/1`):
  expected accidentals in the declared fixed-grid window from 200 shifted
  windows beyond a 200 ns guard (CW) or at the side peaks of a pulsed pump;
  Gamma posterior and an off-peak dispersion check. R1·R2·Δt is a
  cross-check only.
- `calibrate_heralding` (kind `heralding-efficiency`): (C - A) / S with the
  accidentals' Gamma-Poisson uncertainty and a Beta posterior; an exact,
  deterministic interval.
- `calibrate_dcr_temperature` (kind `dcr-temperature`, new in
  `calibration/1`): Arrhenius fit of dark calibrations taken at several
  temperatures (`conditions["temperature_K"]`).
- `dadl calibrate pairs SETUP --tags F --herald A+ --signal B+ --out DIR`
  writes the delay/jitter, accidentals and heralding artifacts, using the
  setup's declared window and delays, and warns when the fitted delay
  disagrees with the declared one (it never moves the window).
- `simulate_pairs` and closed forms `arrhenius_dark_rate`,
  `accidentals_expected` and `heralding_ratio_expected`.
- Coverage gate: `python -m dadl.calibration.coverage {dark,pairs,dcr,all}`;
  the CI job now runs all six calibrations.

- Release: the crate is now published to crates.io by the release workflow with
  Trusted Publishing, in a `crates` environment, instead of by hand with a token.

## 0.0.3 (2026-11-11)

Applies the third review ([docs/reviews/2026-11-council.md](docs/reviews/2026-11-council.md)).

### Breaking

- `two_photon_visibility(..., analysis=...)` is replaced by
  `bell_visibility()` and `hom_visibility()`; the analysis type can no longer
  be chosen by the caller. HOM visibility now uses the dip convention
  (C_max − C_min) / C_max.
- Simulator and closed-form parameters are unit-suffixed to match
  `setup.yaml`: `dark_rate_cps`, `jitter_sigma_ps`, `pair_rate_per_s`.
- `load_setup` returns a typed `Setup` and refuses unknown or duplicate keys,
  wrong types (including `1.0e5`, which YAML reads as a string) and
  non-default policies.
- Rust: every public function returns `Result<_, DadlError>` instead of
  panicking; parameters are named with their units.
- `closed_form.series_hybrid_measured_rate` is replaced by
  `series_measured_rate` (the exact series model); the old formula is kept as
  `patil_usman_hybrid_rate`, documented as empirical.

### Fixed

- Corrected visibility error now propagates the accidentals' variance and
  uses the correct derivatives (the old error was about 10–30% too small).
- Corrected visibility above 1, negative accidentals and non-finite counts
  are refused.
- Unsorted input to dead-time filters and histograms is refused in every
  build, in Rust and Python (release builds used to return wrong answers).
- Overflow near the int64 limits is an error instead of a saturated or
  wrapped value (Rust used to report a false coincidence).
- Float timestamps are refused with a hint to convert to int64 picoseconds.
- The default policy fingerprint is pinned and covers the thresholds; rule
  tables are read-only; looser rules must go through `Policy.allow`, which
  validates the correction name.
- Licence texts ship in the wheel and sdist; the sdist test suite passes.

### Added

- `dadl quickstart [example|setup.yaml]` reads the YAML; `dadl init`,
  `dadl simulate`, `dadl analyze` (own tags in `.npz`/`.csv`); `--version`;
  `python -m dadl`; errors exit 1 with a message.
- Results record `policy_fingerprint`, `overrides` (with reasons) and
  `window`; threshold labels for values within 1.96σ of a threshold and for
  zero-σ results at a boundary; CHSH gets threshold labels.
- Shared Rust/Python test vectors (`tests/vectors/core.json`), a 200-seed
  pull study of the CHSH error, an exact series dead-time test.
- CI: release-mode tests, clippy, MSRV 1.74, sdist/wheel install test,
  Dependabot, a release check that the tag matches both versions and is on
  main, and separate build and publish jobs.
- CHANGELOG, CITATION.cff, a validated example calibration artifact.

## 0.0.2

First published release: `dadl-quantum` on PyPI and `dadl` on crates.io.
