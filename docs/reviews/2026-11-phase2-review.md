# Fourth review: twenty rounds on phase 2 (2026-11-12)

DADL still has no outside reviewers. This review was run by the agent that wrote phase 2, in twenty
rounds, each through one expert lens, over everything merged for phase 2 (milestones 1–4) and the
release path. Suspected defects were confirmed with a probe or a failing test before being fixed.
A self-review is weaker than an independent one: a human physicist's and statistician's review
remains the top open item.

## Rounds

| # | Lens | Scope |
| --- | --- | --- |
| 1 | Statistician | Dark/afterpulse/dead-time fit: likelihood, gradient, Laplace posterior, dead-time interval, upper limit |
| 2 | Statistician | Delay/jitter, accidentals, heralding and DCR(T) posteriors |
| 3 | Statistician | Coverage gate: design, truths, seeds, the two rule changes |
| 4 | Physicist | Simulator: self-exciting dead time and afterpulsing, dark counts, jitter |
| 5 | Physicist | HOM: beamsplitter model, bunching on non-number-resolving detectors, accidentals |
| 6 | Rust engineer | `cross_correlation`, overflow and sorting contract, vectors |
| 7 | FFI engineer | PyO3 bindings: borrowing, GIL, error mapping, ABI |
| 8 | Instrument engineer | PTU T2 format: record layouts, rollovers, markers, units |
| 9 | Adapter contract | Silent drops, provenance, declared capabilities, template |
| 10 | Data integrity | Artifacts: JCS, hashes, write semantics, read verification |
| 11 | Configuration | `setup.yaml` loader against `setup-1` schema; HOM and CHSH branches |
| 12 | CLI and UX | `quickstart`, `simulate`, `analyze`, `calibrate` for both analysis types |
| 13 | Python API | Backend selection, exports, error types, version coupling |
| 14 | Performance | Python reference hot paths, memory, the benchmark |
| 15 | Metrology | Units, naming, "as detected" efficiency, what each estimate is |
| 16 | Release engineer | Release workflow, wheel matrix, Trusted Publishing, version checks |
| 17 | CI engineer | Job coverage, backends, platforms, oldest dependencies |
| 18 | Security | Parsers (YAML, PTU, npz, artifacts), resource limits |
| 19 | Documentation | README, design, site, CHANGELOG against the code |
| 20 | Test quality | What the tests would and would not catch |

## Defects fixed

| Round | Finding | Evidence | Fix |
| --- | --- | --- | --- |
| 7 | The bindings released the GIL while borrowing a NumPy buffer, so another Python thread could write to the array while Rust read it (a data race) | rust-numpy's borrow tracking covers Rust borrows only | The GIL is held while the core runs |
| 2 | Heralding efficiency counted heralds per detection but coincidences per slot; two pairs in one slot halved eta | Probe: two pairs per slot, every herald with its partner, reported 0.5 instead of 1 | Heralds counted per slot (method `herald-ratio` v2), with a note when slots are shared; test added |
| 10 | A crash while writing an artifact left a partial file, and every later write of that ID then failed with "exists with different content" | Code reading | Write to a temporary file and rename atomically; test checks no temporary files remain |
| 10 | A manifest without a valid `posterior` entry raised `KeyError` instead of `ArtifactError` | Probe | Checked and refused with `ArtifactError`; test added |
| 13 | `dadl-quantum` used any installed `dadl-quantum-core`, even from another release | Code reading | Versions must match; otherwise the core is ignored with a warning (or refused under `DADL_BACKEND=rust`); test added |
| 16 | The `dadl-quantum-core` wheel builds ran only on a release tag, so the first test of the release path would have been the release | Workflow reading | CI builds the same four wheels on every pull request and tests the native ones against the shared vectors |
| 1 | With afterpulsing not detected, the dark rate silently included any afterpulses below the detection threshold | Model reading | The artifact notes say so; the bias is bounded by the reported upper limit on `afterpulse_prob` |
| 8 | Docs claimed the PTU `HW_InpChannels` count was "sync plus inputs" without a source, and did not state that the whole file is read into memory | Code reading | Comment corrected; memory use (about 50 bytes per record) stated |

## Checked and found sound

- **Dark fit.** The analytic gradient matches the model in all three coordinates. The dead-time
  interval follows the exponential law of the shortest gap. The upper limit is the largest over
  a grid of afterpulse time constants, so it is conservative.
- **Cross-correlation.** The Rust and Python implementations agree on the 11 shared vectors and
  on randomised streams. All index arithmetic is within range after the `2 * range_ps` overflow check.
- **PTU.** Times and rollovers agree exactly with the independent `ptufile` reader for PicoHarp,
  HydraHarp v1 and v2 and MultiHarp fixtures. Truncation, trailing data and out-of-order records
  are refused.
- **Parsers.** YAML uses a safe loader with a 1 MB limit and duplicate-key refusal. `.npz` and
  `.npy` files are loaded with `allow_pickle=False`. PTU header lengths are checked against the file size.
- **DCR(T) coverage.** It was on the low side on the pre-set seed (DCR_ref 93.9%, Ea 93.4%), so a
  weighted-least-squares bias was suspected. Two independent seeds gave 95.4% / 94.6% and
  94.8% / 95.2%, and pooled over 3,000 replicates coverage is 94.7% and 94.4%: chance, not bias.
- **HOM.** Bunched photons on one detector are recorded at most once, and singles fall as
  eta(2 − eta)/2 per pair (tested). The simulated dip matches the closed form.

## Recorded, not fixed (deferred)

| Round | Item | Why deferred |
| --- | --- | --- |
| 3 | Each gate tests one parameter point per calibration. Low counts, strong afterpulsing and broad IRFs are not covered | Adds CI time; plan a second regime per calibration in phase 3 |
| 3 | The gate rule changed twice after seeing results (milestones 1 and 2), both disclosed | Future changes should be decided before running |
| 8 | PTU files are read whole; T3 mode is refused | Chunked reading and T3 are planned |
| 11 | The `setup-1` schema allows `reference_delay_ps` in a CHSH setup and A/B names in either type; the loader refuses both | Express with `if`/`then` in the schema later |
| 2 | The delay fit assumes a Gaussian IRF; a tail shows as a failing predictive check, not a better fit | The design's IRF has an exponential tail; fit it in phase 3 |
| 19 | No replay file from real PicoQuant hardware | Needs a recording from a lab |

## Follow-up: test coverage and untested regimes (2026-11-12)

**Code coverage.** Line and branch coverage of `python/dadl`, both backends combined, is 89%
(2,661 statements, 676 branches). Most of the rest is refusal paths: bad input that must raise a
clear error. That excludes `calibration/coverage.py`, whose 1,000-replicate gates run in their
own CI job. The Rust core and bindings are not measured yet.

**Regime probe.** The gates test one parameter point per calibration. The dark calibration was
probed at six points, 1,000 replicates each (95% intervals; band about 92.6–97.3%):

| Regime | Dark rate | Afterpulse p | τa | Dead time |
| --- | --- | --- | --- | --- |
| Gate point (1 kcps, p 3%, 20 s) | 94.5% | 94.5% | 94.6% | 95.8% |
| CMOS SPAD (p 20%, τa 20 ns, 10 ns dead time) | 94.5% | 95.3% | 94.7% | 96.2% (was 89.0%) |
| Low counts (2 s) | 94.5% | 94.8% | 94.0% | 94.2% |
| High rate (100 kcps) | 95.2% | 93.7% | 94.9% | 95.8% |
| Long afterpulse (τa 1 µs) | 94.6% | 94.5% | 94.2% | 95.0% |
| SNSPD (100 cps, no afterpulsing) | 94.6% | 99.9% (one-sided limit, conservative by design) | not identified | 95.2% |

**Defects found and fixed.**

| Finding | Evidence | Fix |
| --- | --- | --- |
| Dead-time interval covered 89.0% with fast afterpulsing | The shortest gap equals the dead time exactly in 7% of runs, and the continuous interval ended just below it | Exact discrete quantiles of the shortest gap in whole ticks (method `interarrival-mle` v2); 96.2% |
| A sparse, sharp coincidence peak gave a NaN predictive check, so `calibrate pairs` could not save the artifact | Probe: 150 pairs, 8 ps jitter, no background | Parametric-bootstrap deviance check when there are too few bin groups for a chi-square; test added |
| HOM dropped the accidental-subtracted value without saying why when accidentals exceed C_min | Probe | `counts["subtracted"]`, and the CLI says so |
| The delay fit read 5 ps MultiHarp delays 2 ps high (found by running existing calibrations on hardware-twin data) | Probe | Tick-aware bins (`xcorr-gaussian-mle` v3) |

**What is still needed for higher coverage**, most important first:

1. **A second regime in each CI gate.** Done: the dark gate also runs a CMOS-SPAD regime, and the
   pair gates a sparse low-rate SNSPD regime and a broad-IRF regime (12 gates in all, run in
   parallel in CI). All pass on the pre-set seed 20261111. The broad-IRF regime first crashed the
   delay fit (a starting width from a fixed +-50 bins, an unbounded fit); fixed and tested.
2. **Misspecification tests per calibration.** Only the dark fit (two afterpulse time constants) and
   the delay fit (an exponential IRF tail) are tested against data their model cannot describe.
   Heralding needs multi-pair emission, accidentals need drifting rates, and DCR(T) needs a
   tunnelling term.
3. **A pulsed-source simulator.** The pulsed-pump path of the accidentals calibration (side peaks
   at the repetition period) has only a synthetic test.
4. **Rust coverage.** Run `cargo llvm-cov` on the core and the bindings in CI.
5. **Refusal paths.** About 150 uncovered lines are input errors in `config`, `tagfile`, `artifact`
   and `ptu`. A parametrised table of bad inputs would cover most of them cheaply.
6. **PTU edge cases.** Header types not yet exercised (float arrays, bitsets, dates), T3 mode, and
   one real recording per device family.

## Calibration probe campaign (2026-11-12)

Three rounds of regime and misspecification probes, 300–1,000 replicates per case, until a round
found nothing significant in a well-specified case. Each finding was fixed and tested; the 12 CI
gates were rerun after every change and pass on seed 20261111 with their earlier numbers, except
where a fix was meant to change them (noted).

| Round | Cases | Found | Fix | Before → after |
| --- | --- | --- | --- | --- |
| 1 | 15 | Dead-time interval with few gaps and afterpulsing (method treated f(0) as known, and the dark-count law as certain when afterpulsing was not detected) | Exact tick law averaged over the posterior of f(0); interval up to m when afterpulsing is not ruled out | 0–85% → 96–98% (50–200 gaps) |
| 1 | | Afterpulse p and τa intervals at low counts (detection-dependent switch; Gaussian shape) | Profile-likelihood intervals below 5,000 gaps; conservative grid limit when undetected | p 90.9% → 94.8–99.6%; τa 82–91% → 89–95% |
| 1 | | DCR(T) refused cold or dim scans (200-gap minimum) | Minimum 20 gaps | 100% errors → 93–97.5% coverage |
| 2 | 19 | Narrow peaks (5 ps jitter, 20 ps bins) | Refit with bins of about a third of the jitter | 92% → 94–95% |
| 2 | | Accidentals dispersion check gave false alarms with sparse counts | Exact conditional test below 5 counts per window | 4.5% → 0.1% of good runs flagged |
| 3 | 19 | Low-rate dead time reported as 0 ps (point estimate clipped) | Maximum-likelihood value m, with a note that the data only bound it above | value 0 → m |

Round 3 (fresh seed, 600 replicates per case): every well-specified case inside the 99.9% band,
every predictive check flagging 0.2–1.7% of good runs (nominal 1%), and wrong-model DCR(T) data (an
unmodelled 50 cps tunnelling floor) flagged in 100% of runs. Remaining limits: below about 30 gaps
the afterpulse time constant is rarely identified and covers about 88% when it is; one-sided limits
at p = 0 and sparse accidentals over-cover by design.
