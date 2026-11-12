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
- **HOM.** Bunched photons on one detector are recorded at most once, and singles fall as
  eta(2 − eta)/2 per pair (tested). The simulated dip matches the closed form.

## Recorded, not fixed (deferred)

| Round | Item | Why deferred |
| --- | --- | --- |
| 3 | Each gate tests one parameter point per calibration. Low counts, strong afterpulsing and broad IRFs are not covered | Adds CI time; plan a second regime per calibration in phase 3 |
| 3 | The gate rule changed twice after seeing results (milestones 1 and 2), both disclosed | Future changes should be decided before running |
| 2 | DCR(T) coverage on the pre-set seed is on the low side (93.9% and 93.4%, gate 92.4–97.2%) | See the independent-seed check below |
| 8 | PTU files are read whole; T3 mode is refused | Chunked reading and T3 are planned |
| 11 | The `setup-1` schema allows `reference_delay_ps` in a CHSH setup and A/B names in either type; the loader refuses both | Express with `if`/`then` in the schema later |
| 2 | The delay fit assumes a Gaussian IRF; a tail shows as a failing predictive check, not a better fit | The design's IRF has an exponential tail; fit it in phase 3 |
| 19 | No replay file from real PicoQuant hardware | Needs a recording from a lab |
