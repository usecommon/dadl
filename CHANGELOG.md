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
  the CI job now runs all six calibrations. The chance that a correct
  calibration fails its gate is lowered from 1% to 0.1% (band 92.6–97.1% for
  one parameter, 92.3–97.3% for four): the seed is fixed, so a chance failure
  blocks CI permanently. The delay/jitter gate read 96.9% for jitter on seed
  20261111, above the old 96.8% limit; three other seeds gave 95.3%, 94.7%
  and 94.0%. The accidentals gate first read 91.7% because mirrored off-peak
  windows are correlated; the windows no longer mirror each other (94.0%).

### Phase 2, milestone 3: PTU adapter and the hom-snspd example

- `dadl.adapters`: the `Adapter` protocol (Python subset of the Rust trait),
  `DeviceInfo`, `Capabilities` and typed `Event`s; `adapters/TEMPLATE.md`
  for contributors.
- `dadl.adapters.ptu` (`PtuFile`, `read_ptu`): DADL's own reader for
  PicoQuant PTU files in T2 mode (PicoHarp 300, HydraHarp v1/v2, TimeHarp
  260 N/P, MultiHarp). Exact integer ps conversion when the tick is whole
  ps; rollovers handled (including multi-rollover records); markers,
  rollovers and unknown records reported as events; every record accounted
  for; truncated files, trailing data, out-of-order records and T3 files are
  refused. Tests write fixtures with `python/tests/ptu_writer.py` and, when
  `ptufile` is installed, check that both readers give the same times.
- `analysis.type: hom` in `setup.yaml` (channels A and B; acquisitions "dip"
  and "reference" at `analysis.reference_delay_ps`; `source.coherence_time_ps`),
  `simulate_hom`, `dadl.analysis.hom.hom_from_tags` and the `hom-snspd`
  example: `dadl quickstart hom-snspd`, `dadl simulate` and `dadl analyze`
  work for HOM setups. The MVP tagger is now PicoQuant PTU instead of Swabian.
- Tag files take named settings (`read_tags(..., settings=...)`): HOM uses
  0 = dip, 1 = reference.
- Fix: the delay/jitter fit (method `xcorr-gaussian-mle` v2) treats each
  histogram bin as the whole-ps delays it holds. Version 1 was biased by
  +0.5 ps, visible only when the jitter is a few bin widths.
- `AdapterError` for malformed or unsupported vendor data.

### Phase 2, milestone 4: Python bindings to the Rust core

- `bindings/python`: a PyO3 crate (its own Cargo workspace) built with
  maturin into the optional package `dadl-quantum-core` (module
  `dadl_core`, abi3 wheels for Python 3.10+). It exposes dead time,
  occupied slots, fixed-grid coincidences, inter-arrival and
  cross-correlation histograms over NumPy int64 arrays, reading contiguous
  arrays in place. The GIL is held while the core runs, so no other Python
  thread can write to an array the core is reading.
- `dadl.backend()` and `DADL_BACKEND` (`python` / `rust`): `dadl` uses the
  compiled core when it is installed, else the pure-Python reference.
  Argument types are checked in Python for both; the shared vectors and
  random-stream comparisons run on both backends.
- `bindings/python/bench.py`: on 2,000,000 tags per channel (Apple M-series)
  the core is 61x faster for the non-paralyzable filter, 47x for fixed-grid
  coincidences, 48x for inter-arrival histograms and 5x for
  cross-correlation (already vectorized in NumPy). Recorded, not gated.
- CI: a bindings job builds the core on Windows, Linux and macOS and runs
  every test with `DADL_BACKEND=rust`. Release: abi3 wheels (manylinux2014
  x86-64 and aarch64, macOS universal2, Windows x64) and an sdist for
  `dadl-quantum-core`, published with Trusted Publishing.

### Hardware twin

- `dadl simulate SETUP --out DIR --format ptu [--device hydraharp|multiharp|picoharp|timeharp260p]`
  writes each simulated acquisition as a T2 PTU file (`setting0.ptu`, ...),
  with the device's tick, timer rollovers and header. `--tags DIR` reads such
  a folder in `analyze` and `calibrate`. The CHSH quickstart gives the same S
  through PTU files as through a tag file.
- `dadl.adapters.ptu_write.write_ptu`: a vectorized PTU writer, byte-for-byte
  equal to the independent writer in the test suite. Files it writes carry
  `CreatorSW_Name = DADL` and `DADL_Simulated = True`; the reader marks them
  `simulated`, and calibrations from them are recorded as `backend: sim`.
- `dadl.adapters.SimulatedTagger`: a live tagger emulated from a setup,
  streaming batches and injecting buffer overflows, dropped batches and
  configuration changes as typed events; `check_accounting` verifies that
  tags delivered plus tags reported lost equal tags generated.
- `dadl.simsetup` (`simulate_setup`, `simulate_setup_dark`), moved out of the CLI.
- Fix found by running the existing calibrations on twin data: the delay fit
  (method `xcorr-gaussian-mle` v3) finds the tagger tick from the data, rounds
  the bin width up to whole ticks, and models each bin accordingly. Version 2
  assumed 1 ps ticks and read 5 ps MultiHarp delays 2 ps high.
- Coverage analysis ([review record](docs/reviews/2026-11-phase2-review.md#follow-up-test-coverage-and-untested-regimes-2026-11-12)):
  the dark calibration was probed at six detector regimes. Its dead-time
  interval (method `interarrival-mle` v2) now uses exact quantiles of the
  shortest gap in whole ticks; v1 covered 89% with fast afterpulsing, v2 96%.
  A sparse, sharp coincidence peak no longer gives a NaN predictive check
  (which blocked saving the artifact); HOM says when accidentals are not
  subtracted.

### Twenty-round review of phase 2 ([docs/reviews/2026-11-phase2-review.md](docs/reviews/2026-11-phase2-review.md))

- Bindings hold the GIL while reading NumPy buffers (releasing it allowed a
  data race with other Python threads).
- Heralding efficiency counts heralds per slot (method `herald-ratio` v2):
  two pairs in one slot no longer halve it.
- Artifacts are written atomically; a manifest without a valid posterior
  entry raises `ArtifactError`.
- `dadl-quantum` uses `dadl-quantum-core` only when the versions match.
- CI builds the release wheels of `dadl-quantum-core` on every pull request.
- Notes: the dark rate includes undetected afterpulses (bounded by the upper
  limit); PTU files are read whole.

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
