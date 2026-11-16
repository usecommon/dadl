# Phase 3 plan: inference and plans

Nov 15, 2026 · decisions taken Nov 16, 2026 · builds on 0.1.1

Phase 3 in [design.md](design.md#mvp-scope-and-milestones) delivers "forward model, joint Monte
Carlo, corrected CHSH and HOM; signed analysis plans with trial counting; streaming O(n)
coincidences", gated by `guardrails`. This document turns that row, and the items the design and
the review records defer to phase 3, into milestones of one pull request each. Nothing here is
code yet. Decisions that need the owner are marked **Decision** and collected at the end.

## Starting point

What 0.1.1 gives phase 3: six calibrations as `Calibration` objects with 2,000 posterior
draws, `stand_in` columns, `unidentified` parameters, `conditions` (including `tagger_*` from PTU
headers) and a misfit status; content-addressed artifacts (`dadl.artifact`, `calibration/1`
frozen, with an unused optional `signature` block); Laplace and conjugate posteriors reweighted to
match the profile intervals (W10); analyses that report raw values, or raw plus an
accidentals-subtracted secondary value, with statistical errors only; a self-declared
`Window.source` and `grid_origin_ps` fixed at 0; the `coverage` CI job (10 regimes, 42 checks);
and the series dead-time closed form, tested against a simulated chain but never fitted. There
is no HDF5 code, no signing, and no `dadl verify` beyond self-consistency on read.

In flight elsewhere, not planned here: published results as executable tests (accidental
subtraction, coincidence-time loophole, NIST comparison), and PTU T3 reading plus crosstalk
simulation. Phase 3 consumes them: the published-results tests join the `guardrails` job once
merged (M10), and the crosstalk simulator becomes a misspecification regime (M3).

## The phase 3 gate

Phase 2 was gated by numbers; phase 3 should be too. The design names one job, `guardrails`. A
refusal suite alone does not test the new numbers phase 3 produces (corrected values with
calibration uncertainty), so the gate has three parts, all required on every pull request:

1. **`guardrails`** (design): every regression case refuses or labels as specified, on Windows,
   Linux and macOS. The catalogue is in M10; it includes window shopping under a signed plan:
   a result whose window, delays, grid origin or calibration set differ from the signed plan is
   refused, and every analysis run under a plan increments its trial count.
2. **`propagation`** (new): 1,000 replicates per regime of simulated CHSH and HOM runs, each
   with its own simulated calibration runs. The corrected value's 95% interval (statistical plus
   calibration, by joint Monte Carlo) must cover the injected truth within the same binomial band
   rule as phase 2 (92.6–97.1% for one quantity), on a seed fixed in advance; the raw interval is
   checked against its own expectation. A **control** must fail: the same run with calibration
   uncertainty switched off must fall below the band in at least one declared regime, otherwise
   the regime does not test propagation.
3. **`coverage`** (extended): the phase 2 job plus the series dead-time regimes (M4) and, if
   NUTS is adopted (M5), the same gates run with NUTS posteriors on the fast tier.

**Decision D1 (taken: three-part gate).** Adopt the three-part gate, or keep the design's `guardrails` job
only. If adopted, design.md's phase table row changes to name all three.

## Milestones

One PR each. Dependencies are on merged milestones. Rough order: M0 and M4 can start at once;
M1 → M2 → M3 is the critical path; M6 → M7 needs only M1.

| # | Milestone | Depends on |
| --- | --- | --- |
| M0 | Inference library benchmark (NumPyro, PyMC, Stan) and the decision record | none |
| M1 | Calibration sets: loading, bracketing, joint draws | none |
| M2 | Forward model and corrected HOM and fringe visibility | M1 |
| M3 | Corrected CHSH and the `propagation` gate | M2 |
| M4 | Series dead-time model: simulator and calibration | none |
| M5 | NUTS backend for joint fits | M0, M2 |
| M6 | Signing primitives, trust file and `dadl verify` | M1 |
| M7 | Signed analysis plans, ledger and trial counting | M6 |
| M8 | Photon-HDF5 export with phconvert validation | M1 (M6 for signatures in the file) |
| M9 | Streaming O(n) coincidences and chunked PTU reading | none |
| M10 | `guardrails` job, docs and release 0.2.0 | all |

### M0. Inference library benchmark

Design: "the probabilistic-programming library (NumPyro, PyMC or Stan) is chosen then, after an
install benchmark on Windows" (open question in design.md).

- Add `bench/ppl/` (not shipped in the wheel): one script per library fitting the same reference
  model, the dark inter-arrival likelihood of `models.DarkAfterpulseModel` at 20,000 gaps plus
  the M2 HOM joint model on a small histogram.
- CI workflow `ppl-bench` (manual trigger plus one run per PR to this milestone) on
  windows-latest, ubuntu-latest, macos-latest, Python 3.10 and 3.13, in a clean venv. Record:
  `pip install` success and time with no compiler on PATH, installed size, import time, time to
  first sample (compile), effective samples per second, R-hat, agreement of posterior quantiles
  with the existing Laplace/profile intervals.
- Candidates: NumPyro (JAX CPU wheels on Windows), PyMC (PyTensor; C compiler for speed),
  CmdStanPy (needs a CmdStan toolchain install).
- Output: `docs/reviews/2026-xx-ppl-benchmark.md` with the table and a recommendation.
- Acceptance: all three measured on all three OSes or the failure recorded; licences checked
  (NumPyro Apache-2.0, PyMC Apache-2.0, Stan BSD-3; all compatible with MIT OR Apache-2.0).

**Decision D2 (taken: NumPyro preferred, M0 confirms).** Library choice, after M0. Recommendation going in: NumPyro, as an optional
extra `dadl-quantum[mcmc]`, if JAX installs from wheels on Windows without a compiler; PyMC
second. Stan's separate toolchain install fails the "clean install on a lab PC" test unless the
benchmark shows otherwise. Whatever wins stays optional: the base install keeps NumPy and SciPy
only, and Laplace stays the default for quickstart.

### M1. Calibration sets

The bridge from recorded artifacts to analyses. New module `python/dadl/propagate.py`
(`CalibrationSet`, `draws`).

- `CalibrationSet.load(paths | directory, channels, kinds)`: loads artifacts with
  `load_calibration`, verified on read; refuses misfit artifacts unless the policy allows it
  (label "policy modified").
- Joint draws: columns from one artifact always share a sample index (design: "always drawn
  together"); different artifacts are drawn independently by permuting indices with a seeded
  `SeedSequence` keyed on the artifact IDs (same rule as `draws_rng`, cycle 3.1).
- Refuses to propagate `stand_in` columns as a posterior (base.py comment, cycle 3.3). A
  parameter known only as an upper limit enters as a bound: the corrected value is computed at
  the bound and at zero and reported as a range, with a label. **Decision D3** below.
- Bracketing (design, "Calibration choice"): refuse an artifact whose `valid_for` time range or
  recorded `conditions` do not bracket the dataset, raising `StaleCalibration`; overridable via
  `Policy.allow(analysis, "calibration_outside_validity", reason)`. This also closes the open
  item that `tagger_*` conditions are recorded but not compared at use (cycle 2.8).
- Disjointness: refuse an artifact whose `raw-tags` input hashes overlap the analysed tags'
  hashes (design: calibration data disjoint from analysis data; expert review Q4, "refuse delays
  calibrated on the analysed tags", still open).
- Files: `propagate.py`, `errors.py` (add `StaleCalibration` from the design's error table),
  `guardrails.py` (new correction names), `result.py` (`artifact_ids` filled).
- Tests: `test_propagate.py`: shared index within an artifact, independence across artifacts
  (correlation of a difference, as cycle 3.1), stand-in refusal, bracketing and disjointness
  refusals, override labelling.

### M2. Forward model, corrected HOM and fringe visibility

New module `python/dadl/forward.py` implementing the design's "Forward model specification" for
the two-channel correlation histogram, on reduced data only (histograms and counts; no per-tag
Python loops).

- `expected_histogram(params, bins)`: μ_ab(k) with the Gaussian IRF (the tail waits for the IRF
  item below), the delay, the off-peak accidental term Â_ab from the accidentals artifact, and
  dead-time loss as the conditional live fraction already used in
  `closed_form.heralding_ratio_expected`. Afterpulse and crosstalk terms are parameters fixed at
  zero unless an artifact supplies them; if a channel's dark artifact reports afterpulsing above
  its upper limit and no kernel is propagated, the result carries a label.
- Corrected value = the metric fitted under the forward model with nuisance draws from M1, never
  a chain of subtractions. The joint Monte Carlo pushes each draw set through the fit; the
  reported interval combines the statistical and calibration parts, and both are reported
  separately.
- `hom_visibility` and `bell_visibility` gain `calibrations: CalibrationSet | None`. With none,
  behaviour is unchanged (raw, plus the existing off-peak subtraction where allowed). With one,
  the corrected value comes from the forward model and replaces the subtraction path; asking for
  both is refused ("single model, no double counting").
- `result.py`: `Estimate` gains optional `stderr_calibration`; `UNCERTAINTY_KIND` becomes per
  estimate. `LabelledValue.corrections` lists corrections in the order applied (design,
  "Minimum metadata").
- CLI: `dadl analyze setup.yaml --tags ... --calibrations DIR` prints raw first, then corrected
  with its artifact IDs. The quickstart's "Corrected values arrive in phase 3" line (cli.py:192)
  becomes a corrected line computed from simulated calibration runs.
- Tests: closed-form HOM and Werner fringe with known accidentals and dead time reproduced to
  within Monte Carlo error; threshold-crossed-after-correction label on the 45% → 82% case;
  refusal of stacked corrections; the M3 `propagation` regimes for HOM run here on a fast tier
  (200 replicates) as a smoke test.

### M3. Corrected CHSH and the `propagation` gate

- `chsh(..., calibrations=)`: corrected S as a labelled secondary value from the forward model
  per setting pair, with "fair sampling assumed" when efficiency enters and "accidentals
  modelled; secondary, not a Bell claim" always. The headline stays raw (design guardrail table);
  `THRESHOLD_AFTER_CORRECTION` applies.
- Dead-time selection (cycle 4.10) is modelled, not only labelled: the forward model includes
  each channel's dead time, so corrected S removes the selection bias at the declared rate;
  raw S keeps the label.
- `python/dadl/analysis/propagation.py` (the gate, like `calibration/coverage.py`):
  `python -m dadl.analysis.propagation {hom,chsh,all} --replicates 1000 --jobs N`. Regimes,
  declared before the first run: CW SPAD CHSH at the chsh-spad parameters; high-rate CHSH with
  dead-time selection; HOM on SNSPDs (hom-snspd); bright-background HOM; and one misspecified
  regime (crosstalk injected by the in-flight crosstalk simulator, or afterpulsing not
  propagated), where the pass criterion is that the predictive check flags misfit in at least
  90% of replicates instead of returning confident wrong intervals (design, "Validation").
- CI: `propagation` job alongside `coverage`, split across runners under the same 50-minute
  limit.
- Acceptance: all regimes within band; the uncertainty-off control falls below the band in the
  declared regime(s); misfit share of clean runs at most 1%.

### M4. Series dead-time model

Design: "Phase 2 implements a hard-edge non-paralyzable effective channel dead time only; the
series model below is phase 3"; README limit; expert review Q2.

- Simulator: `DetectorParams` gains `dead_time_p_ps` (paralyzable stage, default 0); `sim.py`
  applies paralyzable then non-paralyzable, as `test_series_chain_matches_exact_closed_form`
  does. The Rust core already has both filters; add a series test vector to `tests/vectors/`.
- Calibration: `calibrate_deadtime_series` fitting τ_p and τ_np in the time domain from an
  illuminated inter-arrival histogram at two or more known rates (dark data cannot identify τ_p:
  paralysis needs arrivals during the dead time). Likelihood: gap density of the series chain for
  a Poisson input, with afterpulsing held at the dark artifact's draws. Model comparison with the
  hard-edge model; τ_p reported as an upper limit when not identified.
- `simulate --deadtime-scan` writes the illuminated runs; `dadl calibrate deadtime`.
- Schema: **Decision D4**, new kind name.
- Tests: closed form via the simulated chain (exists), recovery at three (τ_p, τ_np) points,
  `coverage` regimes `series-spad` and `passive-quench` (1,000 replicates each, band as phase 2),
  method version pinned in `method_versions.json`.
- Stays out: soft (continuous) recovery of SNSPDs and passive-quench efficiency recovery curves.
  The README limit is reworded, not removed.

### M5. NUTS backend

- `python/dadl/inference/` with `laplace.py` (moved from `stats.py`, unchanged numbers) and
  `nuts.py` (the M0 library, imported lazily). Selected by `--sampler {laplace,nuts}`; default
  `laplace`.
- First users: the M2/M3 joint fits at low counts, and the dark fit below 5,000 gaps where draws
  needed reweighting (W10). Not every calibration moves; each move bumps the method version.
- `checks` records R-hat, bulk and tail ESS and divergences; R-hat > 1.01 or any divergence makes
  the artifact a misfit. Draws thinned to `N_DRAWS`.
- Reproducibility: NUTS output is not bit-reproducible across platforms (design, "Content
  addressing"), so artifact IDs are stable on one stack only; `dadl verify` (M6) compares
  statistically.
- Tests: Laplace and NUTS intervals agree where Laplace is valid (20,000-gap dark fit); `coverage`
  fast tier with `--sampler nuts`; optional-dependency tests skip cleanly without the extra.
- CI: one job with the `mcmc` extra on Windows, Linux and macOS.

### M6. Signing primitives and `dadl verify`

Design, "Signing": payload `"dadl-calibration-v1\0"` + JCS bytes, base64url `sig`, `key_id` =
SHA-256 of the public key, per-chain `seq`, an `allowed_signers`-style trust file; each site
holds its own key.

- `python/dadl/signing.py`: ed25519 sign and verify, key generation (`dadl key new`, labelled
  "unverified local dev key" by default), trust file parsing (key IDs, identities, validity,
  revocation). **Decision D5** on the crypto dependency.
- The artifact ID stays the hash of the manifest without `signature`, so signing an existing
  artifact does not change its ID and old IDs stay valid. `seq` and the signed-over ID go into the
  signature block: a new optional field, allowed under the frozen-schema rule.
- `dadl sign ARTIFACT --key K` and `dadl verify PATH...`: integrity (exists today), signature
  against the trust file, chain checks (`supersedes`, monotonic `seq`, no forks), and with
  `--recompute` re-runs the calibration from its hashed inputs when they are available and
  checks agreement: each summary value within its interval's Monte Carlo tolerance, not byte
  equality. Artifacts whose inputs are missing are reported "unverifiable", not failed.
- Manifest-to-ID test vectors (design: "planned") added to `tests/vectors/`, checked by Python
  and, for JCS and SHA-256, by the Rust core.
- Tests: tampered manifest, wrong key, revoked key, rollback (`seq` lower than seen), fork, and a
  recompute that disagrees. All join the `guardrails` catalogue.

### M7. Signed analysis plans, ledger and trial counting

Closes the `source: plan` gap (design "What phase 1 cannot check"; third review "Window
`source: plan` is self-declared"; cycle 5.4 window shopping; cycle 1.18 unrecorded grid phase).

- `python/dadl/plan.py`, schema `schema/plan-1.schema.json`. A plan holds: analysis type, window
  width, declared delays and `grid_origin_ps`, the delay calibration artifact ID (Bell timing
  from a separate, frozen segment), sample size and stopping rule, exclusion rules (setting-blind
  telemetry only), guardrail policy fingerprint, and the setup hash. Signed with M6.
- Ledger: an append-only, hash-chained JSON-lines file per project (`plans/ledger.jsonl`), each
  entry signed. Entries: plan registered, analysis run (inputs, artifact set, result hash),
  plan superseded. `dadl plan new|sign|show|log`.
- Analyses take `plan=` instead of a bare `Window` for Bell types. `Window(source="plan")` without
  a verified plan is labelled "plan not verified" (or refused: **Decision D6**). A result whose
  window, delays, grid origin, calibration set or policy differ from the plan is refused.
- Trial counting: every analysis under a plan appends an entry, and the result reports how many
  runs, windows and artifact sets were tried under that plan (design, "Analysis-plan ledger").
  Combining Bell results across runs with a test martingale needs CH/Eberhard and stays out; a
  result under a plan with more than one trial carries a "N trials under this plan" label and no
  combined significance.
- Exclusion: excluded fraction and the estimate on all data reported when the plan has
  exclusion rules.
- What a signature cannot show: that the plan existed before acquisition. The ledger orders the
  plan before the data's own timestamps (PTU `file_created_local`, `acquired`) and refuses a plan
  registered after them, but a lab controlling its own clock can backdate. Stated as a limit;
  trusted timestamps (RFC 3161) stay out (design, "Security design").
- Quickstart and `init` create and sign a plan with the dev key, so the starter kit shows the
  path without extra steps.
- Tests: window shopping (best of 72 windows and delays from cycle 5.4) under a signed plan is
  refused and counted; an edited plan fails verification; a ledger with a removed entry fails the
  chain check; a delay artifact computed from the analysed tags is refused (with M1).

### M8. Photon-HDF5 export

Design, "Data model": write under `/user/dadl_calibration`, fill `/setup/detectors` with posterior
medians in SI units, `timestamps_unit` = 1e-12 s for int64 ps tags, validated with phconvert.
Promised for phase 2, marked planned in cycle 8.13.

- `python/dadl/hdf5.py` behind the optional extra `dadl-quantum[hdf5]` (h5py; phconvert in the
  test extra only). `dadl export SETUP --tags ... --calibrations DIR --out file.h5`.
- Writes `/photon_data` (timestamps, detectors, `timestamps_specs`), `/setup`, `/identity`,
  `/provenance`, and `/user/dadl_calibration/{artifacts,channels,events,analysis}` holding each
  artifact's canonical manifest bytes and posterior table unchanged, so IDs and signatures verify
  from the file.
- Reader `read_photon_hdf5` for round trips; `raw-tags` hashes cover datasets, so writing the
  calibration group into the same file does not change them (test).
- Tests: phconvert's validator passes on every exported example; round trip of tags and
  artifacts; the artifact verifies by ID and signature from inside the file.
- **Decision D7** on proposing `/calibration` upstream now or after a user exists.

### M9. Streaming O(n) coincidences and chunked PTU reading

The design table puts streaming coincidences in phase 3. PTU files are read whole at about 100
bytes per record (README limit; cycle 2.10).

- Rust core: `grid_coincidences` and `occupied_slots` over chunks with carried state at chunk
  edges; Python binding taking an iterator of chunks. The pure-Python reference stays the
  oracle; shared vectors include chunk boundaries inside a slot and inside a dead time.
- `adapters/ptu.py`: chunked T2 decoding with rollover state carried across chunks; whole-file
  and chunked reads must agree record for record. (T3 is in flight elsewhere; coordinate so the
  chunked decoder covers both.)
- Acceptance: memory bounded by the chunk size (test with a 1 GB synthetic file in the bindings
  job, recorded, not gated); results identical to the whole-file path on all examples.

### M10. `guardrails` job, docs and release 0.2.0

- `python/tests/guardrails/` collects the regression catalogue, run as the `guardrails` CI job on
  three OSes. Each case names the result it must refuse or label:
  1. accidental subtraction raising Bell visibility past 1/√2 (45% → 82%);
  2. coincidence-time loophole: data-dependent or event-centred window;
  3. window shopping under a signed plan (M7) and a plan edited after signing;
  4. delays calibrated on the analysed tags (M1, Q4);
  5. calibration outside its validity or conditions, and tagger settings that differ (M1);
  6. stacked corrections (subtraction plus forward model) (M2);
  7. threshold crossed only after correction, for S, fringe and HOM visibility (M2, M3);
  8. stand-in draws propagated as a posterior (M1);
  9. tampered, rolled-back or forked artifacts and ledgers (M6, M7);
  10. the published-results tests from the in-flight work, once merged.
- Docs: design.md phase table (phase 3 Done, gate named), "What works today", the forward model
  and signing sections reconciled with what was built; README "Known limits" updated (corrected
  values, series model, plans and what signatures do not prove); CHANGELOG 0.2.0.
- A review cycle in the format of `docs/reviews/2026-11-review-cycles.md`, before release, with
  at least a Bayesian computation, a Bell-test red team and a security lens.

### Smaller items assigned to phase 3

From the review records; each can ride with the milestone named or go to 0.2.x.

| Item | Source | Where |
| --- | --- | --- |
| IRF with exponential tail (Gaussian core plus tail) in the delay fit and forward model | phase 2 review round 2; README "Gaussian timing response" | Own PR after M2, **Decision D8** |
| `grid_origin_ps` declared rather than fixed at 0 | cycle 1.18 | M7 |
| Pair-labelled simulation truth for the heralding gate | cycle 7.14; summary of cycles | M3 (the propagation simulator labels pairs) |
| Rotating seeds for the coverage gates | expert review Q11 | M10, if the owner agrees to change the gate rule before running |
| Stored draws gated, not only intervals | cycle 7.17 | M1 (the draws are now used) |

## Decisions

The owner took all nine on Nov 16, 2026, each as recommended below. D2 is a
preference going into M0; the benchmark still has to confirm it.

| # | Question | Options | Decided |
| --- | --- | --- | --- |
| D1 | What gates phase 3 | `guardrails` only (design) / plus `propagation` and extended `coverage` | Three-part gate: corrected intervals are the new numbers and need a coverage test |
| D2 | MCMC library | NumPyro / PyMC / Stan, after M0 | NumPyro as an optional extra if it installs from wheels on Windows; Laplace stays default |
| D3 | Upper-limit parameters (unidentified dead time, undetected afterpulsing) in corrections | Range at bound and zero / uniform prior / refuse | Range plus a label; never a uniform prior presented as a posterior |
| D4 | Schema for the series dead-time artifact | New kind `deadtime-series` in `calibration/1` (enum widened) / reuse `dark-afterpulse-deadtime` with a new method / `calibration/2` | New kind in `calibration/1`, documented as the one permitted kind of additive change; 0.1.x readers refuse it, which is acceptable |
| D5 | Crypto dependency for ed25519 | `cryptography` (optional extra) / PyNaCl / shell out to `ssh-keygen -Y` or minisign | `cryptography` as the `sign` extra; accept SSH signatures in the trust file later |
| D6 | `source: "plan"` without a verified plan, for CHSH | Refuse / label "plan not verified" | Label in 0.2.0, refuse in 0.3.0 after one release of warning; the quickstart uses a signed dev plan from day one |
| D7 | Upstream Photon-HDF5 `/calibration` proposal | Propose with M8 / ship under `/user` and propose after a user exists | Ship under `/user`; propose after one external user or the hardware pilot |
| D8 | IRF tail in phase 3 | In phase 3 / 0.2.x | In phase 3: without it SNSPD delay artifacts are refused, so hom-snspd cannot use corrected values on real data |
| D9 | Should calibration artifacts be signed in phase 3 (design: signing v0.2+) | Plans only / plans and artifacts | Both: same primitive, little extra code; HSM, TPM and QKD attestation stay out |

## Risks

| Risk | Impact | Mitigation |
| --- | --- | --- |
| Corrected values read as the result | Overcorrection in print, the project's main risk | Raw headline unchanged; corrected values labelled, after raw, with artifact IDs; threshold-after-correction refusals in the gate |
| Joint Monte Carlo runtime in CI | `propagation` job over the time limit | Reduced-data forward model only; fast tier per PR, full tier split across runners; budget measured in M2 before regimes are fixed |
| Propagation gate passes trivially | False confidence | Uncertainty-off control that must fail; misspecified regime that must be flagged |
| JAX or PyMC install problems on Windows | Users without MCMC | Optional extra; Laplace default; M0 before any dependency is added |
| Signatures mistaken for proof of timing or honesty | Overclaimed plans | README and result text say what is shown (unchanged plan, ordering in the ledger) and what is not (time of creation, honest data) |
| Series dead-time not identifiable at lab rates | Wide or prior-driven τ_p | Requires a rate scan; τ_p reported as a bound when not identified; prior-sensitivity check |
| Schema changes after the freeze | Broken readers | Only additive fields; D4 decided before M4 merges |
| Scope against team size | Phase slips | M0, M4, M9 are independent and can be deferred to 0.2.x without blocking the gate; M1–M3, M6, M7, M10 are the gate's minimum |

## Out of phase 3

- Explain engine and counterfactual attribution (phase 4).
- MCP server, approval tokens, audit-log writer process, monitor daemon, change-point detection.
- State-space drift filter between calibrations, and `drift-filter-state` artifacts.
- CH/Eberhard, test martingales and combined significance across runs; tomography; g(2).
- Pulsed sources and multi-pair statistics beyond the CW term; rate-dependent efficiency.
- Crosstalk calibration (the simulation is in flight; the calibration is v0.2+), multi-component
  afterpulse kernels, surviving-afterpulse likelihood (cycle 7.13), soft SNSPD recovery.
- HSM, TPM or OS-keystore signer process; RFC 3161 timestamps; QKD attestation.
- Swabian, qutools and IDQ adapters; Arrow C Stream adapters; plugin hash pinning.
- Operator features: more than four inputs per setup, `--json`, history and latest-artifact lookup
  (cycle 8.18).
- In flight elsewhere: published-results executable tests, PTU T3 reading, crosstalk simulation.
