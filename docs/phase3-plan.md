# Phase 3 plan: inference and plans

Nov 15, 2026 · decisions taken Nov 16, 2026 · builds on 0.1.1

Phase 3 in [design.md](design.md#mvp-scope-and-milestones) delivers "forward model, joint Monte
Carlo, corrected CHSH and HOM; signed analysis plans with trial counting; streaming O(n)
coincidences", gated by `guardrails`. This document turns that row, and the items the design and
the review records defer to phase 3, into milestones of one pull request each. M0-M3 are merged
and M4 is in review (#26) (CHANGELOG, Unreleased); the rest is plan. Decisions that need the owner are marked **Decision** and collected at the end.

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
   rule as phase 2 (92.6–97.1% for one quantity), on a seed fixed in advance. A **control** must fail: the same run with calibration
   uncertainty switched off must fall below the band in at least one declared regime, otherwise
   the regime does not test propagation.
3. **`coverage`** (extended): the phase 2 job plus the series dead-time regimes (M4) and, if
   NUTS is adopted (M5b), the same gates run with NUTS posteriors on the fast tier.

**Decision D1 (taken: three-part gate).** Adopt the three-part gate, or keep the design's `guardrails` job
only. If adopted, design.md's phase table row changes to name all three.

## Milestones

One PR each. Dependencies are on merged milestones. Rough order: M0 and M4 can start at once;
M1 → M2 → M3 is the critical path; M6 → M7 needs only M1. Council round 8 split M5 into M5a and
M5b and moved most of M10's smaller items into M2b, M3b and M4b.

| # | Milestone | Depends on |
| --- | --- | --- |
| M0 | Inference library benchmark (NumPyro, PyMC, Stan) and the decision record | none |
| M1 | Calibration sets: loading, bracketing, joint draws | none |
| M2 | Forward model and corrected HOM and fringe visibility | M1 |
| M2b | IRF with an exponential tail (D8), with the common-vs-confined misfit test | M2 |
| M3 | Corrected CHSH and the `propagation` gate | M2 |
| M3b | Corrected-CHSH additions: efficiency normalisation, D3 union, pair-labelled heralding truth | M3 |
| M4 | Series dead-time model: simulator and calibration | none |
| M4b | Pre-registered dead-time gates and regimes, coarse and Fano checks | M4 |
| M5a | Series dead time at use in the forward model | M2, M4 |
| M5b | NUTS backend for joint fits | M0, M5a, the B1 decision |
| M6 | Signing primitives, trust file and `dadl verify` | M1 |
| M7 | Signed analysis plans, ledger and trial counting | M6 |
| M8 | Photon-HDF5 export with phconvert validation | M1, M6 |
| M9 | Streaming O(n) coincidences and chunked PTU reading | none (its divided-T3 check needs a real file: external, O10) |
| M10 | `guardrails` job, docs and release 0.2.0 | the gate minimum (Risks), plus the optional milestones that merged |

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

**Decision D2 (taken: NumPyro, confirmed by M0 in [the benchmark record](reviews/2026-11-ppl-benchmark.md)).** Library choice, after M0. Recommendation going in: NumPyro, as an optional
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
- CI: `propagation` job alongside `coverage`, split across runners under the 40-45-minute step
  timeouts.
- Acceptance: all regimes within band; the uncertainty-off control falls below the band in the
  declared regime(s); misfit artifacts at most 2% and forward flags at most 1% of clean runs.

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
  `coverage` regimes `series-spad` and `passive-quench` (1,000 replicates each; τ_p band as phase 2;
  τ_np at least 95%, one-sided, see `coverage.deadtime_coverage`),
  method version pinned in `method_versions.json`.
- Stays out: soft (continuous) recovery of SNSPDs and passive-quench efficiency recovery curves.
  The README limit is reworded, not removed.

### M5a. Series dead time at use

The forward model reads the dark dead time only; the `deadtime-series` artifacts that M4 writes
have no consumer, and user-facing text says so ("not yet supported; see Known limits"). Council
round 8 moved the entry items that do not need NUTS here.

- Forward-model live fraction from an identified `deadtime-series` artifact (m/n from the series
  closed form; M4 council round 1, D4), with a propagation regime that uses it.
- Minimum method version at use (council round 5, V4): a table (kind, method) -> version in
  `propagate` that refuses or labels older artifacts, with a test, before anything reads
  `deadtime-series` artifacts.
- Shared row index (B5): a series tau_p's draws lose their correlation (about 0.3) with the dark
  afterpulse draws they were propagated from; `CalibrationSet.draws` needs a shared row index
  across a series artifact and its cited dark artifact, with the reviewer's correlation test.
- Scan range at use (T6): check the analysed rate against a series artifact's recorded scan
  rates, renaming `scan_rate_low_cps` / `scan_rate_high_cps` to the `_min_`/`_max_` bracket
  convention when the check is wired.
- The dead time at the analysed gaps (council round 8, Ph4 and Ph5, one change): truncate every
  dead-time draw, point and bound at the analysed shortest gap plus the D5 allowance, use the
  additive pair live fraction (1 - fa)(1 - fb) + r_c min tau, and redesign the control regime.
  Truncation alone gave 92.5% with a 93.0% control; both together are nearly unbiased (mean z
  +0.13) but make chsh-dead-time cover with its control (94.0% / 94.0%), so the control regime
  must change to one whose dead-time uncertainty survives truncation (owner item O8).
- The dead fraction above 20% (R1 b), refused today, is in scope once the series model is read.
- Whether M5a joins the gate minimum is owner item O6 (council: yes).

### M5b. NUTS backend

- `python/dadl/inference/` with `laplace.py` (moved from `stats.py`, unchanged numbers) and
  `nuts.py` (the M0 library, imported lazily). Selected by `--sampler {laplace,nuts}`; default
  `laplace`.
- First users, pending the B1 decision below: the M2/M3 forward model has no likelihood of its own,
  so "joint fits at low counts" exist only if B1 adds one; otherwise the dark fit below 5,000 gaps
  where draws needed reweighting (W10). Not every calibration moves; each move bumps the method
  version.
- `checks` records R-hat, bulk and tail ESS and divergences; R-hat > 1.01 or any divergence makes
  the artifact a misfit. Draws thinned to `N_DRAWS`.
- Reproducibility: NUTS output is not bit-reproducible across platforms (design, "Content
  addressing"), so artifact IDs are stable on one stack only; `dadl verify` (M6) compares
  statistically.
- Tests: Laplace and NUTS intervals agree where Laplace is valid (20,000-gap dark fit); `coverage`
  fast tier with `--sampler nuts`; optional-dependency tests skip cleanly without the extra.
- CI: one job with the `mcmc` extra on Windows, Linux and macOS.
- Entry items (council round 6; V4, B5 and T6 moved to M5a in round 8):
  - **What NUTS samples (B1; owner item).** Name the model NUTS samples (likelihood, parameters,
    priors): a new joint likelihood (the dip or fringe histogram with the IRF, or calibration and
    counts together), or state "calibration fits only". The M2/M3 forward model has no likelihood
    of its own (closed form with Monte Carlo draws). Re-run `bench/ppl` on the chosen model.
  - **Starts and modes (B2).** Chains from overdispersed starts (the Laplace covariance scaled by
    4-9, plus one chain per known branch or mode), so multimodality shows in R-hat; branch
    enumeration and unions stay outside the sampler; the start rule is recorded in `checks`.
  - **Boundaries (B3).** Dead time stays conditioned on the shortest gap m with exact draws;
    tau_p = tau_np sigmoid(z) and smooth transforms for the fit boxes, with proper priors.
    Divergences are recorded, but the "any divergence makes a misfit" rule applies only after the
    reparametrization is tested on passive-quench and series-spad.
  - **Priors (B4).** Stored draws are posteriors under flat priors in the fit's coordinates (log
    rate, logit p, log tau_a) inside the fit's box; rank checks confirm this (round 6). Near p = 0
    that prior is improper without the box. Declare proper priors per calibration, record them in
    the artifact, and bump the method version where draws move.
  - **Stored draws gated (B7, review item 7.17).** A rank-statistic test per propagated
    calibration (as the round-6 `sbc_dark.py` probe), run with `--sampler nuts`.

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
bytes per record in T2 and up to about 170 in T3 (README limit; cycle 2.10; council round 7, S3).

Entry item (council round 7, G5): check the T3 sync divider convention on one real divided T3
file (the nsync step against `TTResult_SyncRate`). Until then `calibrate dark` and `calibrate
deadtime` refuse T3 files with `HWSync_Divider` > 1 and analyses label them. The file is an
external dependency (lab data, owner item O10); without it M9 ships the refusal only.

- Rust core: `grid_coincidences` and `occupied_slots` over chunks with carried state at chunk
  edges; Python binding taking an iterator of chunks. The pure-Python reference stays the
  oracle; shared vectors include chunk boundaries inside a slot and inside a dead time.
- `adapters/ptu.py`: chunked T2 decoding with rollover state carried across chunks; whole-file
  and chunked reads must agree record for record. (T3 is in flight elsewhere; coordinate so the
  chunked decoder covers both.)
- Acceptance: memory bounded by the chunk size (test with a 1 GB synthetic file in the bindings
  job, recorded, not gated); results identical to the whole-file path on all examples.

### M2b. IRF tail (D8)

- IRF with an exponential tail (Gaussian core plus tail) in the delay fit and the forward model
  (phase 2 review round 2; README "Gaussian timing response"). Without it SNSPD delay artifacts
  are refused, so hom-snspd cannot use corrected values on real data (D8).
- The common-vs-confined misfit test for CHSH: a timing tail common to all channels does not bias
  S, so the predictive check need not refuse it (M4 council round 2, B5).
- Whether M2b is in the gate minimum is for the owner (D8 says "in phase 3").

### M3b. Corrected-CHSH additions

- Opt-in per-detector efficiency normalisation for corrected S: a new correction name, a policy
  line and a gate regime (M4 council round 2, B4).
- D3 union interval (zero to the bound) for corrected values at a dead-time bound, and a
  chsh-dead-time regime whose dark runs identify the dead time or that gates the range (M4 council
  round 2, B2, G; pre-registered with O8).
- Pair-labelled simulation truth for the heralding gate (review item 7.14): the propagation
  simulator does not label pairs yet.
- `dadl simulate --heralding --pair A-,B+`: heralding runs for the other three CHSH channel pairs
  (council round 6, R6).
- May slip to 0.2.x (owner item O13).

### M4b. Dead-time gates and regimes

Every change here alters a gate or a regime, so all are pre-registered in one owner sitting
before any run (owner item O8).

- Rotating seeds for the coverage gates (expert review Q11), with a second seed (round 8, St2).
- tau_np dead-time gate with power (M4 council round 1, G1).
- A coverage regime with an unidentified tau_p (council round 6, T7 c); a regime past the
  count-rate maximum with the user declaration and the per-run branch flip search (M4 council
  rounds 3 and 4); a series regime at tau_p 5 ns / tau_np 22 ns (round 7, Q5).
- Regimes closer to real labs (council round 8, Ph6): per-channel spread (efficiency +-20%, dead
  time +-10%), tagger jitter of 10-20 ps after the dead time, an SNSPD at 2-5 Mcps per channel,
  and a SPAD with an exponential IRF tail (recording its refusal rate).
- Coarse grouped predictive check and the count-variance (Fano) check for the series dead time
  (M4 council round 1, S5 and D1; each changes the misfit share, so the gate is re-run).
- The 0.5 s passive-quench probe and the second-order pending-afterpulse term, if a regime needs
  it (M4 council round 5, I2).
- Gates with M5a (owner item O13).

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
  values, series model, plans and what signatures do not prove); the site (`docs/index.html`),
  which still describes 0.1.1 as released; CHANGELOG 0.2.0.
- A review cycle in the format of `docs/reviews/2026-11-review-cycles.md`, before release, with
  at least a Bayesian computation, a Bell-test red team and a security lens.

### Smaller items assigned to phase 3

From the review records; each can ride with the milestone named or go to 0.2.x. Council round 8
moved most rows from M10 to M2b, M3b, M4b and M5a.

| Item | Source | Where |
| --- | --- | --- |
| IRF with exponential tail (Gaussian core plus tail) in the delay fit and forward model | phase 2 review round 2; README "Gaussian timing response" | M2b, **Decision D8** |
| Common-vs-confined misfit test for CHSH: a timing tail common to all channels does not bias S, so the predictive check need not refuse it | M4 council round 2 (B5) | M2b |
| `grid_origin_ps` declared rather than fixed at 0 | cycle 1.18 | M7 |
| Pair-labelled simulation truth for the heralding gate | cycle 7.14; summary of cycles | M3b (M3 merged without it; the propagation simulator does not label pairs) |
| Rotating seeds for the coverage gates, with a second seed | expert review Q11; council round 8 (St2) | M4b (owner item O8) |
| τ_np dead-time gate with power: two-sided around its exact expected coverage, or a goodness-of-fit gate on the gap - τ_np histogram, with a random sub-tick dead-time phase in the regime | M4 council round 1 (G1) | M4b (owner item O8) |
| Coarse grouped predictive check, Bonferroni-combined with the fine one, for the series dead time | M4 council round 1 (S5) | M4b (changes the misfit share; gate re-run) |
| Count-variance (Fano) stationarity check for dead-time scans | M4 council round 1 (D1) | M4b (needs the model's renewal Fano factor; gate re-run) |
| Opt-in per-detector efficiency normalisation for corrected S (new correction name, policy line and a gate regime) | M4 council round 2 (B4) | M3b; 0.2.x if not done |
| D3 union interval (zero to the bound) for corrected values at a dead-time bound, and a chsh-dead-time regime whose dark runs identify the dead time or that gates the range | M4 council round 2 (B2, G) | M3b (regime pre-registered with O8) |
| Forward model live fraction m/n from the series closed form when the set holds an identified `deadtime-series` artifact | M4 council round 1 (D4) | M5a (propagation gate re-run) |
| The dead time at the analysed gaps (truncation, the additive pair term and a new control regime, as one change) | council round 8 (Ph4, Ph5) | M5a (owner item O8 for the regime) |
| Stored draws gated, not only intervals (a rank-statistic test per propagated calibration) | cycle 7.17; council round 6 (B7) | M5b entry item (moved from M1; 7.17 stays open with that target) |
| `dadl simulate --heralding --pair A-,B+`: heralding runs for the other three CHSH channel pairs | council round 6 (R6) | M3b |
| A coverage regime with an unidentified tau_p, to gate the upper limit's level (one-sided about 97.5%) | council round 6 (T7 c) | M4b (owner item O8) |
| Coverage regime with a dead-time scan past the count-rate maximum (1/tau_p) | M4 council round 3 (N2) | M4b (owner item O8) |
| User declaration of scan runs past the count-rate maximum | M4 council round 4 (DT3) | M4b, with the past-maximum regime |
| Per-run branch flip search to restore the passive-quench gate's speed | M4 council round 4 | M4b, checked on the past-maximum regime |
| A 0.5 s passive-quench probe (0.5 s coverage 92.5%, matching the no-afterpulse control's 93.1%: check the profile's calibration at 1e6 gaps) | M4 council round 5 (I2) | M4b |
| The second-order pending-afterpulse term (two pending, or one beside the record's own), if a regime needs it | M4 council round 5 (I2) | M4b |
| Minimum method version at use: a table (kind, method) -> version in `propagate` that refuses or labels older artifacts, with a test | M4 council round 5 (V4) | M5a entry item, before anything reads `deadtime-series` artifacts |
| A series coverage regime at tau_p 5 ns / tau_np 22 ns, 5-30 Mcps (the third point of `test_recovery_at_three_points`, in no regime: 3 of 20 seeds missed there) | M4 council round 7 (Q5) | M4b (owner item O8) |
| Gate regimes closer to real labs: per-channel spread, tagger jitter after the dead time, an SNSPD at Mcps, a SPAD IRF tail | council round 8 (Ph6) | M4b (owner item O8) |
| Pulsed sources: side-peak accidentals at multiples of the repetition period, so T3 input can be corrected without declaring the source CW | M4 council round 7 (G3) | v0.2+ (with the pulsed-source item below); round 8 refuses pulsed light seen off-peak on any input |
| The predictive check on T3 tags with a measured sync period | M4 council round 7 (G1, Fix 1) | Done in council round 8: histograms bin on the dtime lattice (`lattice_ps`) |

## Decisions

The owner took all nine on Nov 16, 2026 (D10, D11 later, during M2), each as recommended below.
D2 is a preference going into M0; the benchmark still has to confirm it.

| # | Question | Options | Decided |
| --- | --- | --- | --- |
| D1 | What gates phase 3 | `guardrails` only (design) / plus `propagation` and extended `coverage` | Three-part gate: corrected intervals are the new numbers and need a coverage test |
| D2 | MCMC library | NumPyro / PyMC / Stan, after M0 | NumPyro as an optional extra if it installs from wheels on Windows; Laplace stays default |
| D3 | Upper-limit parameters (unidentified dead time, undetected afterpulsing) in corrections | Range at bound and zero / uniform prior / refuse | Range plus a label; never a uniform prior presented as a posterior. **As built:** the value and its interval are at the bound; the two ends (zero and the bound) are points in `counts`, labelled when they differ by more than the error. The union interval is M3b. |
| D4 | Schema for the series dead-time artifact | New kind `deadtime-series` in `calibration/1` (enum widened) / reuse `dark-afterpulse-deadtime` with a new method / `calibration/2` | New kind in `calibration/1`, documented as the one permitted kind of additive change; 0.1.x readers refuse it, which is acceptable |
| D5 | Crypto dependency for ed25519 | `cryptography` (optional extra) / PyNaCl / shell out to `ssh-keygen -Y` or minisign | `cryptography` as the `sign` extra; accept SSH signatures in the trust file later |
| D6 | `source: "plan"` without a verified plan, for CHSH | Refuse / label "plan not verified" | Label in 0.2.0, refuse in 0.3.0 after one release of warning; the quickstart uses a signed dev plan from day one |
| D7 | Upstream Photon-HDF5 `/calibration` proposal | Propose with M8 / ship under `/user` and propose after a user exists | Ship under `/user`; propose after one external user or the hardware pilot |
| D8 | IRF tail in phase 3 | In phase 3 / 0.2.x | In phase 3: without it SNSPD delay artifacts are refused, so hom-snspd cannot use corrected values on real data |
| D9 | Should calibration artifacts be signed in phase 3 (design: signing v0.2+) | Plans only / plans and artifacts | Both: same primitive, little extra code; HSM, TPM and QKD attestation stay out |
| D10 | Dead time in the default HOM correction rule (taken during M2) | Keep an explicit override / add to the default rule | Keep the override: `Policy.allow("hom_visibility", "dead_time", reason)`, labelled "policy modified"; the pinned default policy fingerprint is unchanged |
| D11 | Window capture (fixed-grid loss at a HOM reference delay, from a delay-jitter artifact) in the default HOM rule (taken during M2) | Keep an explicit override / add to the default rule | Keep the override, as D10: `Policy.allow("hom_visibility", "window_capture", reason)`; the default policy fingerprint is unchanged |

## Owner items

Consolidated from council rounds 4-8. Each needs the owner; the council's recommendation is in
the last column.

| # | Item | Source | Status | Recommendation |
| --- | --- | --- | --- | --- |
| O1 | Make "Statistical gates" (which now needs the slow-test job) a required check on `main` | r4 C1, r8 Pl5 | Open | Edit ruleset 27075146 |
| O2 | Release requires gates on the tagged SHA (or a RELEASING policy); forbid direct pushes to main | r7 S6 | Open | Require it in `release.yml`'s `check` job |
| O3 | Dev version (still 0.1.1; delay is at method 10, deadtime at 7) | r5 V2, r6, r7 | Open | 0.2.0.dev0 before more artifacts |
| O4 | Calibration schema `$id` | r5 V7 | Open | With O3 |
| O5 | D4 decision text ("0.1.x readers refuse") | r5 V6 | Open | Reword as r5 L2 |
| O6 | M5 scope (B1), now M5a / M5b; whether M5a joins the gate minimum | r6 B1, r8 Pl1 | Open | Split; M5a in the gate minimum |
| O7 | Accept the behaviour changes in release notes: refusal above 20% dead fraction; T3 needs `source.type: cw`; `t3-range` refusal; off-peak refusal of pulsed light on T2 and plain tags; forward check at 4e-4 and drift at 3.72 sd; the 1e-3 per-analysis refusal budget as the stated standard. Round 9: on T3 only `source.type: cw` declares the clock (the templates' `spdc-cw` is refused); the gap tests have their own override `dead_time_contradicted` (`calibration_outside_validity` no longer covers them); the off-peak check can now refuse in CHSH; the upward gap test tolerates a few short gaps | r6, r7, r8, r9 | Partly in CHANGELOG | Accept |
| O8 | Pre-register gate and regime changes in one sitting: rotating seeds with a second seed, the tau_np gate, the D3 union with chsh-dead-time, past-maximum with declaration and flip search, unidentified tau_p, 5/22 ns, regimes closer to real labs (Ph6), a control regime that survives the M5a truncation | r1-r8 | Open | One sitting before any M4b or M5a run |
| O9 | Purge the bench binaries from history | r7 S5 | Removed from the tree | Leave (default) |
| O10 | A real divided T3 file | r7 G5 | Open | Source it, or ship the refusal |
| O11 | Withdrawn deferral (pending afterpulse moved into M4) | r5 item 4 | Plan updated | Close |
| O12 | Propagation gate refusal tolerance: replace "any refusal fails" with at most 5 refused per 1,000 (binomial at 1e-3, `FALSE_FAIL`), refused counted as not covered; the misfit (at most 18) and flag caps in the same binomial form; check the p < 0.01 share too | r8 St2, St5 | New | Adopt (a gate rule, so the owner's) |
| O13 | Which of M2b, M3b and M4b may slip to 0.2.x | r8 Pl2 | New | M4b gates with M5a; M3b may slip |
| O14 | PR #26 merge method, and whether fixes to merged milestones go in their own PRs or reviews are capped at a number of rounds | r8 Pl6 | New | Owner's choice (the reviewer advises a non-squash merge so the round commits stay traceable) |
| O15 | Whether the example and `dadl init` setups keep `source.type: spdc-cw`, or ship `cw` with a comment for T3 users (changing a template changes what new users declare) | r9 Rd3 | New | Keep `spdc-cw`; the T3 refusal text explains |
| O16 | Extend M9 to chunked analyses: coincidences, shifted counts, cross-correlation histograms and gap statistics over time blocks, plus the per-pair difference extraction (Pf8). Plan text and milestone scope | r9 Pf4, Pf8 | New | Adopt into M9; README limit until then |
| O17 | The analytic gradient of the peak fit (Pf7) in the M10 performance pass, with its acceptance run (forward-check p within optimiser tolerance, flagged counts equal at 200 replicates) as the entry test | r9 Pf7 | New | Schedule in M10 |
| O18 | CI budget: the gate jobs are 477-494 of about 700-790 job-minutes, and the 20-job concurrency limit sets the wall time. Buy concurrency, or run the propagation gate nightly and on release rather than on every push | r9 Pf5 | New | Owner's choice |

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
| Scope against team size | Phase slips | M0, M4, M9 are independent and can be deferred to 0.2.x without blocking the gate; M1–M3, M6, M7, M10 are the gate's minimum. M5a is proposed for it (owner item O6); which of M2b, M3b and M4b may slip is owner item O13 |

## Out of phase 3

- Explain engine and counterfactual attribution (phase 4).
- MCP server, approval tokens, audit-log writer process, monitor daemon, change-point detection.
- State-space drift filter between calibrations, and `drift-filter-state` artifacts.
- CH/Eberhard, test martingales and combined significance across runs; tomography; g(2).
- Pulsed sources and multi-pair statistics beyond the CW term; rate-dependent efficiency.
- Crosstalk calibration (the simulation is in flight; the calibration is v0.2+), multi-component
  afterpulse kernels, surviving-afterpulse likelihood (cycle 7.13), soft SNSPD recovery.
- The reverse series order for actively quenched modules (non-paralyzable hold-off, then a
  paralyzable reset window; M4 council round 1, D2).
- HSM, TPM or OS-keystore signer process; RFC 3161 timestamps; QKD attestation.
- Swabian, qutools and IDQ adapters; Arrow C Stream adapters; plugin hash pinning.
- Operator features: more than four inputs per setup, `--json`, history and latest-artifact lookup
  (cycle 8.18).
- In flight elsewhere: published-results executable tests, PTU T3 reading, crosstalk simulation.
