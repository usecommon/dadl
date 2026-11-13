# Changelog

All notable changes. DADL is pre-alpha: in 0.x anything may change between
releases. The `calibration/1` artifact schema is frozen from 0.1.0 (new
optional fields may be added; changes to existing ones need `calibration/2`)
and is versioned separately from the packages.

## 0.1.0 (2026-11-13)

Phase 2: calibrations with simulation-checked 95% intervals, content-addressed
artifacts, PicoQuant PTU input, a hardware twin and the optional compiled core,
after an outside-expert review of the assumptions and ten further review
cycles. Pre-alpha: read the README's "Known limits".

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
- New dependency: SciPy (>= 1.13).

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
- Coverage gates run at more than one regime: the dark calibration at a Si
  SPAD and a CMOS SPAD (20% afterpulsing), the pair calibrations at a CW SPAD
  source, a sparse low-rate SNSPD source and a broad (about 1.8 ns) IRF:
  12 gates, run in parallel (`--jobs`) in the existing CI check. In the sparse
  regime about 0.2 off-peak counts are expected in all, so the accidentals
  gate there requires coverage of at least the band's lower edge, 92.6% (declared
  before the first run).
- Fix found by the broad-IRF regime: the delay fit crashed on broad peaks
  (overflow, or a singular Hessian) because its starting width came from a
  fixed +-50 bins. The start now widens until it holds +-4 sigma, the fit is
  bounded, and when the background fits to zero the delay and jitter are
  reported conditional on it, with a note, instead of failing.
- Calibration probe campaign (three rounds; [review record](docs/reviews/2026-11-phase2-review.md#calibration-probe-campaign-2026-11-12)):
  - dark calibration v4: the dead-time interval averages the shortest-gap law
    over the uncertainty of f(0) and runs up to m when afterpulsing is not
    detected (it covered as little as 0% with a few hundred gaps); afterpulse
    intervals come from the profile likelihood below 5,000 gaps; the minimum
    is 20 gaps instead of 200, so cold DCR(T) runs calibrate; at low dark
    rates the dead time is reported as the bound it is, not as 0 ps;
  - delay: narrow peaks are refitted with finer bins (92% to 94-95%);
  - accidentals: an exact dispersion check for sparse counts (false alarms
    from 4.5% to 0.1%).
  - the coverage gates split each regime's replicates across processes
    (results identical to a serial run) and the profile intervals compute to
    the precision they need: all 12 gates in 95 s on 4 local cores, after the
    low-count fixes had pushed the CI job to 36 minutes.

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

### Outside-expert review of the assumptions ([docs/reviews/2026-11-expert-assumptions.md](docs/reviews/2026-11-expert-assumptions.md))

Five reviewers (detector physicist, statistician, metrologist, time-tagging
engineer, quantum optics experimentalist) rated 13 assumptions wrong and 14
questionable. Changes:

- **Heralding efficiency** (method `herald-ratio` v3): per-slot estimator
  η = (q − p)/(1 − p) with the accidental probability p propagated exactly;
  unbiased up to 20% accidentals (the old ratio read η(1 − p)).
  `calibrate pairs` refuses CHSH and HOM setups unless `--no-analysers`
  confirms a heralding run, and records delays, duration, singles rates and
  analysers in `conditions`. `dadl simulate --heralding` writes such a run.
- **Accidentals** (v2): the off-peak guard is at least 10 µs and ten afterpulse
  time constants past the dead time (afterpulses after true pairs had biased
  it +9% to +165%); a trend test joins the dispersion check.
- **Misfit status**: a predictive check below p = 0.001 sets
  `checks.status: misfit`; `dadl calibrate` refuses to write it unless
  `--accept-misfit` is given.
- **PTU data loss**: non-complete `TTResult_StopReason` codes become
  `data-loss` events and `incomplete` metadata, refused by the CLI; the
  duration comes from the data; PicoHarp channel values above 4 (sync 0, router inputs 1-4) are refused; the
  sync divider is reported.
- **Channel identity**: artifacts record the tagger serial and input, with the
  role and detector serial alongside; hardware setups need `hardware.tagger`.
- `source.type` is validated (`spdc-cw` only).
- **Exact visibility intervals** (Clopper–Pearson) with `NEAR_THRESHOLD_EXACT`
  and `LOW_COUNTS` labels; the delta method covered 63–93% near V = 1.
- **Dark calibration** (v5): "effective channel dead time"; up to three
  short gaps far below all others are dropped as double triggers (a single
  spurious gap near the dead time cannot be told apart and is noted, not
  dropped); profile intervals keyed on the afterpulse count;
  posterior draws reweighted to agree with the intervals; a stationarity check
  (split-half rate, longest gap) flags drift and latches as misfit.
- **DCR(T)**: refuses runs at different bias; records the temperature range;
  notes the SPAD-only validity and few-gap runs.
- **Delay** (v4): the tick comes from the adapter (the tag GCD is a
  cross-check); a two-level background models the dead-time step (bright
  background coverage 52% → 98%).
- `check_accounting` compares tag times, not counts; twin presets enforce the
  real input counts; a MultiHarp 150 N (80 ps) preset.
- Bell results carry "fair sampling and fair coincidence assumed", static
  settings and loophole labels, and coincidences per singles; HOM results are
  labelled as a two-point visibility at the declared delay.
- Artifacts gain an `uncertainty` block (statistical, 95%, components not
  included).
- Coverage gates: new regimes with afterpulsing SPADs, a bright background and
  an ideal heralding check against the true η; accidentals are judged against
  each replicate's own expectation.
- Docs: README uncertainty wording and "Known limits"; design.md notes on the
  dead-time model and misfit refusal; "signed" becomes "content-addressed,
  signing planned".

### Ten review cycles ([docs/reviews/2026-11-review-cycles.md](docs/reviews/2026-11-review-cycles.md))

Ten further cycles of AI-agent role-play review (three reviewers each, 30 in
all) found 47 problems rated wrong, all fixed or stated as limits; the
register lists every finding and its status. Method versions in this release:
dark `interarrival-mle` 7, delay `xcorr-gaussian-mle` 8, heralding
`herald-ratio` 6, accidentals `offpeak-grid` 3, DCR(T) `arrhenius-wls` 2; a
test now fails if a method's code changes without a new version.

**Breaking changes**

- `Window(width_ps, source, ...)`: `source` is required (no default "plan");
  `delays_ps` is read-only, and delays for unknown channels are refused.
- `write_tags(..., simulated=...)` is required.
- `write_artifact` refuses misfit calibrations unless `accept_misfit=True`.
- CHSH and HOM refuse duplicate timestamps and need one channel at delay 0;
  `bell_visibility` and `hom_visibility` need whole counts.
- `hom_from_tags` compares rates across acquisitions of different length;
  `hom_visibility(..., exposure_ratio=)`.
- Tag files must be in time order per channel (a clock reset is refused, not
  sorted).
- Corrected visibility above 1 is left out and labelled instead of raising.
- Artifact manifests: `summary` entries may have `value: null` with `bound`,
  plus `definition`; `posterior.stand_in`; named `inputs`; `acquired` and
  `instrument`; software versions; UCUM units. Non-converged fits are misfits.
- `as_record()` writes JSON cells and adds `dadl_version` and `simulated`.
- `dadl simulate` refuses to overwrite without `--force`; `--device` needs
  `--format ptu`; `--duration-s` applies to `--dark` and `--heralding` only.
- SciPy 1.13 or later.

**New**

- Labels: interrupted (latch, gaps at the start, middle or end of a run),
  dead-time selection, output imbalance, simulated input, accidentals
  subtracted, correction not computed, window not given, single basis.
- `calibrate pairs`: per-artifact refusal, window capture, coincidences to
  accidentals and pair rate (secondary), peak-placement check, declared dead
  time and afterpulse time constant passed; `calibrate dark`: per-channel
  handling, refusal of source-on and unit-error data.
- Record options: `--condition`, `--acquired`, `--created`, `--valid-until`,
  `--supersedes`; `load_calibration`; `dadl-quantum[core]` extra.
- Delay fit: exact tick model, two-level background, multi-start and Newton
  polishing, jitter upper limit below the tick; dark fit: time constants up to
  1 ms (free-running InGaAs); heralding: closed form with the conditional live
  fraction.
- Coverage gates: 10 regimes and 42 checks plus a misfit-share limit, split
  across CI jobs; `--regimes`.
- Input safety: YAML alias expansion and npz decompression limits, PTU header
  checks, histogram size limit in both backends, GIL declared for
  free-threaded Python.

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
