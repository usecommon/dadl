# Phase 3 review council (2026-11)

Phase 3 (milestones M0-M4 on branch `m4/series-deadtime`) went through ten review rounds before merge. Each
round had three reviewers, each through a different expert lens, working read-only against the branch's HEAD and
backing every finding with their own probe scripts through the package's functions. A council of three seats
(P = physicist, M = maintainer, U = user advocate) then re-checked the findings numerically and ruled on each: a
finding is adopted with two of three votes, as a fix (ADOPT-FIX), a stated limit (ADOPT-LIMIT), a deferral to a
named milestone (DEFER), an owner decision, or rejected. The next round started only after the adopted fixes were
committed, one commit per round. From round 8 on, every analysis-time refusal shares one budget: on clean data
within the tested regimes DADL refuses an analysis with probability at most 1 in 1,000, all checks together, and a
label that fires on more than 1 in 1,000 clean analyses states its rate.

Reviewers and council seats are AI agents playing outside experts, not people; a human domain review is still
open. The rulings and the reviewer reports are kept verbatim in [phase3/](phase3/) (paths scrubbed; probe
directories and data files not copied).

An id is cited as "r{N} {Id}" or "council round N, {Id}" (for example r8 Pl6, council round 9, Rd1). Within a
round the ids are unique; a council row that merges findings lists both ids ("S6 + D8"). Severity is the
reviewer's: **wrong** (wrong numbers or misleading in realistic use), **questionable** (needs evidence or a stated
limit), **minor** (documentation, polish); round 10's merge review also used **blocker**. "Where fixed" gives
the round's commit; each section names the method versions that commit moved (deadtime = `series-gap-mle`,
delay = `xcorr-gaussian-mle`).

## Round 1: M4 series dead time and gates

Lenses: point-process statistician (S), detector physicist (D), gate statistician (G). Ruling: [phase3/r1-council.md](phase3/r1-council.md); reports: [r1-stat.md](phase3/r1-stat.md), [r1-phys.md](phase3/r1-phys.md), [r1-gate.md](phase3/r1-gate.md).

Fixes in 666df7b; method versions: deadtime 1 -> 2.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| S1 | wrong | Interruption check treats separate scan runs as one acquisition: unequal durations or one running clock give status misfit on good data | ADOPT-FIX | 666df7b |
| S2 | wrong | `survival` drops the second real mode past 8 tau_p; tail 20-80% low near n tau_p = 1; tau_p undercovers | ADOPT-FIX | 666df7b |
| S3 | wrong | Afterpulse kernel jump snapped to the grid: profile NLL jagged by units, intervals partly grid noise | ADOPT-FIX | 666df7b |
| S4 | wrong | Afterpulses pending at the next gap not modelled: tau_p biased low (-1.3 sd at p 5%, tau_a 200 ns, 1e7 cps), status ok | ADOPT-FIX (guard) + DEFER (model) | 666df7b; deferred (plan) |
| S5 | questionable | Pearson check on thousands of fine groups has little power against localized misfits | ADOPT-LIMIT + DEFER | limit (below); deferred (plan) |
| S6 + D8 | minor | tau_p point (and upper limit) can exceed tau_np | ADOPT-FIX | 666df7b |
| S7 | minor | With tau_p MLE = 0 the afterpulse uncertainty is not added to the upper limit | ADOPT-FIX | 666df7b |
| D1 | questionable | Slow intensity modulation biases tau_p low and passes the check; tolerance unstated | ADOPT-LIMIT + DEFER | limit (below); deferred (plan) |
| D2 | questionable | Series order (P then NP) is not an actively quenched SPAD's (NP hold-off then paralyzable reset window); regime comment says it is | ADOPT-LIMIT | limit (below) |
| D3 | questionable | README says gradual passive-quench recovery "shows as a misfit"; often it is absorbed as tau_p close to tau_np | ADOPT-LIMIT | limit (below) |
| D4 | questionable | Forward model ignores an identified tau_p: live fraction high by up to 2.8-3.7% per channel at the 20% gate | ADOPT-LIMIT + DEFER | limit (below); deferred (plan) |
| D5 | minor | tau_np interval excludes the truth on a real tagger (tagger jitter after the dead time) | ADOPT-LIMIT | limit (below) |
| D6 | minor | Default scan rates count photons, not avalanches; with efficiency 0.5 the scan only reaches n tau = 0.25 and tau_p goes unidentified | ADOPT-FIX | 666df7b |
| D7 | minor | Unrecorded avalanches restart tau_p but never afterpulse or emit crosstalk | ADOPT-LIMIT | limit (below) |
| G1 | questionable | tau_np at-least gate: stated reason wrong (about 1.5 gaps per tick, not tens of thousands); simulator's integer phase over-covers; gate has little power | ADOPT-FIX (text) + DEFER (gate) | 666df7b; deferred (plan) |
| G2 | minor | Plan's gate rules drifted from code (raw-interval check, 2% misfit cap, tau_np one-sided) | ADOPT-FIX | 666df7b |
| G3 | minor | `coverage.merge` takes `gate`/`at_least` from the first shard unchecked | ADOPT-FIX | 666df7b |
| G4 | minor | One-sided rule reuses the two-sided lower edge: false-fail 1.6e-4, not 5e-4 | ADOPT-LIMIT | limit (below) |

Limits of round 1 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | deadtime.py module docstring, end of "Predictive check" |
| L2 | README "Hard-edge dead time" bullet; deadtime.py docstring after "under steady light"; `calibrate deadtime` help |
| L3 | README, replacing "Soft (continuous) recovery ... shows as a misfit." |
| L4 | README, replacing "The analyses' forward model uses the dark calibration's non-paralyzable dead time only."; same sentence in CHANGELOG "Stays out" |
| L5 | README; `DEFINITIONS["dead_time_ps"]` gains "...; the shortest recorded gap, so it includes the tagger's timing jitter" |
| L6 | sim.py module docstring, README |
| L7 | coverage.py, beside the at-least rule in `CoverageResult._ok` or `deadtime_coverage` |

## Round 2: corrected values: CHSH, HOM/fringe, propagation

Lenses: Bell-test experimentalist (B), HOM/fringe experimentalist (H), propagation statistician (R). Ruling: [phase3/r2-council.md](phase3/r2-council.md); reports: [r2-bell.md](phase3/r2-bell.md), [r2-hom.md](phase3/r2-hom.md), [r2-prop.md](phase3/r2-prop.md).

Fixes in 98bbc25; method versions: none changed.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| B1 | wrong | Delay errors on A+ and B+ (misdeclared, or drift after a delay calibration) raise raw and corrected S past 2 (2.20 at 400 ps for 1.98) with no label; the predictive check refits the delay freely | ADOPT-FIX | 98bbc25 |
| B2 + R2 | wrong | CHSH at a dead-time bound: no value at zero, label promises one; D3 not implemented for CHSH | ADOPT-FIX | 98bbc25 |
| R1 | wrong | Two-point at a bound: point at zero, interval at the bound; point can sit outside its interval; threshold crossing at the bound unlabelled | ADOPT-FIX | 98bbc25 |
| R3 | wrong | Stand-in jitter (coarse tagger, legitimate ok artifact) crashes every corrected analysis with KeyError | ADOPT-FIX | 98bbc25 |
| H1 | wrong | Fringe: a pump drift between max and min moves V across 1/sqrt(2), unlabelled; singles show it | ADOPT-FIX | 98bbc25 |
| H2 | wrong | `hom_from_tags(calibrations=...)` without `reference_delay_ps` treats the reference as unsplit and drops the capture label (V 10 sigma low, unlabelled) | ADOPT-FIX | 98bbc25 |
| H3 | questionable | Predictive-check fit stalls from the floor background: p ~1e-206 on clean-enough data, refused as Misspecified | ADOPT-FIX | 98bbc25 |
| B3 | questionable | Corrected point (shortest-gap dead time) sits about 0.4 sd low in its own interval; CHANGELOG gives no sign | ADOPT-LIMIT | limit (below) |
| B4 | questionable | Unequal efficiencies: corrected S labelled "detector effects modelled" keeps the whole imbalance bias (z +13.5, 0/6 covered) | ADOPT-FIX (label) + DEFER | 98bbc25; deferred (plan) |
| B5 | questionable | Predictive check refuses Bell runs with mild common IRF tails that do not bias S; only a blanket override | ADOPT-LIMIT + DEFER | limit (below); deferred (plan) |
| B6 + H5 | minor | README stale: "Not yet: ... corrected CHSH", "statistical only", crosstalk "no analysis flags it", `--calibrations` HOM only, python/README quickstart line | ADOPT-FIX | 98bbc25 |
| H4 | questionable | Forward-model corrected V above 1 reported unlabelled; README says refused | ADOPT-FIX | 98bbc25 |
| H6 | minor | CLI "V reads low" warning printed even when window capture is modelled | ADOPT-FIX | 98bbc25 |
| H7 | minor | TWO_POINT does not mention a reference delay inside the coherence time | ADOPT-FIX | 98bbc25 |
| R4 | questionable | Without a DatasetContext the use-time refusals are skipped silently; slices of the analysed run pass disjointness | ADOPT-FIX + ADOPT-LIMIT | 98bbc25; limit (below) |
| R5 | questionable | Validity starts at write time, not acquisition; a before/after pair cannot be combined | ADOPT-FIX + ADOPT-LIMIT | 98bbc25; limit (below) |
| R6 | questionable | Unused deadtime-series artifacts cited (and can refuse an analysis); BOUND_ONLY set by any stand-in column (afterpulse_prob on every quiet dark) | ADOPT-FIX | 98bbc25 |
| R7 | questionable | M2's required "afterpulsing not propagated" label missing | ADOPT-FIX | 98bbc25 |
| R8 | minor | Gate stores 400 draws, users get 2,000; Monte Carlo error unstated | ADOPT-LIMIT | limit (below) |
| G (from B2's evidence) | - | chsh-dead-time corrects at least one channel at its bound in 40/40 replicates; the control regime's "calibration uncertainty" is partly a fixed bound | ADOPT-LIMIT + DEFER | limit (below); deferred (plan) |

Limits of round 2 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | forward.py docstrings of `_dead_times` and `_channel_dead_times`; CHANGELOG M3 bullet |
| L2 | README "Unequal output efficiencies" bullet; `corrected_chsh` docstring |
| L3 | README Bell limits, after the corrected-CHSH paragraph |
| L4 | propagate.py module docstring, "Disjointness" |
| L5 | propagate.py docstring, README calibration paragraph |
| L6 | propagation.py module docstring |
| L7 | propagation.py, comment on chsh-dead-time |

## Round 3: goodness of fit, guardrail red team, numerics

Lenses: goodness-of-fit specialist (G), guardrail red team (R), numerical analyst (N). Ruling: [phase3/r3-council.md](phase3/r3-council.md); reports: [r3-gof.md](phase3/r3-gof.md), [r3-red.md](phase3/r3-red.md), [r3-num.md](phase3/r3-num.md).

Fixes in ede8b76; method versions: deadtime 3.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| G1 + G2 + N1 (+ G4) | wrong | CHSH capture balance and delay drift: narrow peaks binned too coarsely (N1), wide peaks fitted in too small a region (G2), drift refusal with no noise term (G1); false raw labels and false refusals (3-65%) on clean, correctly calibrated data | ADOPT-FIX | ede8b76 |
| G1 (label part) | wrong | Label, not refuse, when the shift exceeds 0.5 se but is within the noise | REJECT | not adopted |
| G3 | wrong | Predictive check's L-BFGS-B fit stalls when the start delay is tens of ps off; refuses clean data (2/60) | ADOPT-FIX | ede8b76 |
| G4 | minor | Data-side fit pulls the ±400 ps offsets toward the declared delay (drift statistic +0.0011 on the null) | ADOPT-LIMIT | limit (below) |
| G5 | minor | `CAPTURE_IMBALANCE` says "declared delays miss the pair peaks" when only the window is narrow against the jitter | ADOPT-FIX | ede8b76 |
| R1 | wrong | Another detector's dark artifacts accepted on the same tagger input (identity ignores `detector`), getting around the bound refusal; a 1 us dead time applied to a run with 30 ns gaps | ADOPT-FIX | ede8b76 |
| R2 | wrong | Bracketing fails open: no `valid_for` = valid forever unlabelled; recorded temperature or bracket not given for the dataset is skipped; result does not record the dataset context | ADOPT-FIX | ede8b76 |
| R3 | wrong | Corrected CHSH says "dead time and window capture modelled" when neither was applied | ADOPT-FIX | ede8b76 |
| R4 | questionable | One-per-kind-and-channel check depends on channel order; user can pick the artifact by `--created` | ADOPT-FIX | ede8b76 |
| R5 | minor | CHSH CLI never prints the D3 range | ADOPT-FIX | ede8b76 |
| R6 | minor | Docs say validity and tagger settings are not enforced; code enforces them | ADOPT-FIX | ede8b76 |
| N2 | wrong | Series gap law: Newton stays on the low-rate branch; a scan past the count-rate maximum biases tau_p to tau_np (4/4 intervals exclude the truth) | ADOPT-FIX | ede8b76 |
| N3 | minor | `survival`/`gap_cdf` crash on an int tau_p | ADOPT-FIX | ede8b76 |

Limits of round 3 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | chsh.py `_delay_drift` docstring; README corrected-CHSH limits bullet, after the B1 numbers |
| L2 | chsh.py `_delay_drift` docstring |

## Round 4: dead-time statistics, research software, CI

Lenses: dead-time statistician (DT), research software engineer (S), CI engineer (C). Ruling: [phase3/r4-council.md](phase3/r4-council.md); reports: [r4-dt.md](phase3/r4-dt.md), [r4-rse.md](phase3/r4-rse.md), [r4-ci.md](phase3/r4-ci.md).

Fixes in ef8bc0c; method versions: deadtime 4.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| DT1 | wrong | Afterpulse-sensitivity refit searches both branches and lands on the other minimum (shift sd up to 1.5e4 ps; "identified" interval becomes most of [0, tau_np]) | ADOPT-FIX | ef8bc0c |
| DT2 | wrong | Both-branch MLE puts a below-maximum scan on the high branch in about 5% (tau_p about 31 ns, rate 9e7, false pending note) | ADOPT-FIX | ef8bc0c |
| DT3 | questionable | Joined interval is a patch: widens on both sides (p90 12.9 vs 4.7 ns), its note fires exactly when the held fit is rejected | ADOPT-FIX | ef8bc0c |
| DT4a | questionable | Gap check's 0.9 factor has no basis: dead times up to about 16% too long pass | ADOPT-FIX | ef8bc0c |
| DT4b | questionable | A calibrated dead time shorter than the detector's is never checked (22 ns artifact on a 50 ns detector passes, +2.7% live fraction) | ADOPT-FIX | ef8bc0c |
| DT5 | minor | Docstring and CHANGELOG justify the join by "high branch can fit better by up to a few NLL" (Asimov: worse by 2.5-5.6) | ADOPT-FIX | ef8bc0c |
| S1 | wrong | CHSH runs the shape check with no delay-jitter artifact; HOM/fringe skip it; refuses accidentals-only corrected S | ADOPT-FIX | ef8bc0c |
| S2 | questionable | Partial dark set in HOM/fringe gives a bare `KeyError` (CLI traceback) | ADOPT-FIX | ef8bc0c |
| S3 | questionable | Detector-serial check (round 3 R1) runs only through `CalibrationSet.load(channels=)`; API paths skip it | ADOPT-FIX | ef8bc0c |
| S4 | questionable | Tests survive removal of the CHSH shape check, the gap check in all three paths, and `capture_sd_data` | ADOPT-FIX | ef8bc0c |
| S5 | questionable | `cli._grid_capture` (quad) returns about 0 for a narrow peak in a wide window: false "window misses the pair peak", pair rate "high by 1/capture" | ADOPT-FIX | ef8bc0c |
| S6 | minor | Gap check off without record when the floor is <= 0 | ADOPT-FIX (counts only) | ef8bc0c |
| S7 | minor | `calibrate_deadtime_series`: duplicate timestamp gives a numpy error; same run twice accepted; afterpulse artifact kind/channel/misfit unchecked | ADOPT-FIX | ef8bc0c |
| S8 | minor | Counts keys and docs drift: `subtracted` 0.0 with subtraction; counts keys undocumented; CLI help "HOM: corrected" only; schema `valid_for` "not enforced" and `KeyError 'from'`; CHSH quickstart's fixed dark note | ADOPT-FIX | ef8bc0c |
| C1 | questionable | No gate is a required check; a failed shard makes the merge job skipped (counts as passing) | ADOPT-FIX + owner item | ef8bc0c; owner item |
| C2 | questionable | chsh-spad propagation estimated 27-41 min of 45 | ADOPT-FIX | ef8bc0c |
| C3 | minor | Propagation shards leave a worker idle (chunks 50/50/25) | ADOPT-FIX | ef8bc0c |
| C4 | minor | Coverage gate prints nothing until done | ADOPT-FIX | ef8bc0c |
| C5 | minor | Comments understate passive-quench cost (4.6x, 19 s a replicate, not 3x / 8 s) | ADOPT-FIX | ef8bc0c |
| C6 | minor | CI-structure test checks only coverage regimes | ADOPT-FIX | ef8bc0c |
| C7 | minor | CONTRIBUTING has no gate instructions | ADOPT-FIX | ef8bc0c |
| C8 | minor | No `concurrency:` group; about 51 jobs per push queue | ADOPT-FIX | ef8bc0c |
| (council) | - | Per-run flip search in place of the both-branch grid (about 2.5x faster gate) | DEFER | deferred (plan) |

Limits of round 4 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | forward.py `check_dead_time_gaps` docstring; README use-time checks paragraph |
| L2 | deadtime.py module docstring, "Two rates per gap law", replacing the round-3 text |
| L3 | CHANGELOG round-4 entry |

## Round 5: interval coverage, provenance, claims

Lenses: interval statistician (I), provenance archivist (V), claims auditor (C). Ruling: [phase3/r5-council.md](phase3/r5-council.md); reports: [r5-int.md](phase3/r5-int.md), [r5-prov.md](phase3/r5-prov.md), [r5-claims.md](phase3/r5-claims.md).

Fixes in 7bcc32b; method versions: deadtime 5.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| I1 | wrong | `BRANCH_MARGIN` 3 lets about 0.5% of below-maximum scans switch to tau_p about 30 ns, every one a miss | ADOPT-FIX | 7bcc32b |
| I2 | wrong | Pending afterpulses unmodelled: tau_p bias -170 ps that does not shrink (88% at 0.5 s); guard fixed-fraction and its note understates tenfold | ADOPT-FIX (model) | 7bcc32b |
| I3 | minor | Coverage history reads as drift; record the miss split and switch count | ADOPT-FIX (docs) | 7bcc32b |
| V1 | wrong | `deadtime-series` does not cite its dark artifact | ADOPT-FIX | 7bcc32b |
| V2 | questionable | Package still says 0.1.1; `valid_for.from` rule changed under the same version | ADOPT-LIMIT + owner item | limit (below); owner item |
| V3 | questionable | Corrected records lack the Monte Carlo seed and draw count | ADOPT-FIX | 7bcc32b |
| V4 | questionable | Method versions never checked at use | DEFER (M5) | deferred (plan) |
| V5 | questionable | `valid_for.conditions` validates but is ignored | ADOPT-FIX | 7bcc32b |
| V6 | minor | "readers from 0.1.x refuse" is false | ADOPT-LIMIT | limit (below) |
| V7 | minor | Widened schema keeps the v0.1.0 `$id` | ADOPT-LIMIT + owner item | limit (below); owner item |
| C1 | questionable | Current passive-quench coverage figure recorded nowhere | ADOPT-FIX (docs) | 7bcc32b |
| C2 | questionable | Gap-check floor quoted only for a favourable case | ADOPT-LIMIT | limit (below) |
| C3 | questionable | README "What works today (0.1.0)" lists unreleased features | ADOPT-LIMIT | limit (below) |
| C4 | minor | design.md phase 3 "Not started" | ADOPT-LIMIT | limit (below) |
| C5 | minor | Coverage `--jobs` help stale; `run()` families omit deadtime | ADOPT-FIX (docs) | 7bcc32b |
| C6 | minor | CI comment "pushes to main all run" is false | ADOPT-FIX (docs) | 7bcc32b |
| C7 | minor | CONTRIBUTING calls every calibration module pinned | ADOPT-FIX (docs) | 7bcc32b |
| C8 | minor | `--rates` default holds for four settings only | ADOPT-FIX (docs) | 7bcc32b |
| C9 | minor | Counts-key table wrong about raw/corrected keys | ADOPT-FIX (docs) | 7bcc32b |
| C10 | minor | README says fringe has a `--calibrations` path | ADOPT-FIX (docs) | 7bcc32b |
| C11 | minor | Plan and site stale (a-d) | ADOPT-FIX (a-c) + DEFER (d, M10) | 7bcc32b; deferred (plan) |
| C12 | minor | Docstrings/help miss refusals and limits | ADOPT-FIX (docs) | 7bcc32b |
| C13 | minor | design.md says afterpulsing is fitted jointly | ADOPT-FIX (docs) | 7bcc32b |

Limits of round 5 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | CHANGELOG Unreleased, README "Calibration artifacts" |
| L2 | schema `kind.description`, both copies |
| L3 | schema top-level description, both copies |
| L4 | README:382-387 and forward.py:253-258, same words |
| L5 | README |
| L6 | design.md:621 |
| L7 | deadtime.py docstring "Afterpulsing", README dead-time section |

## Round 6: Bayesian computation, journal referee, metrology

Lenses: journal referee (R), metrologist (T), Bayesian computation (B). Ruling: [phase3/r6-council.md](phase3/r6-council.md); reports: [r6-ref.md](phase3/r6-ref.md), [r6-metro.md](phase3/r6-metro.md), [r6-bayes.md](phase3/r6-bayes.md).

Fixes in 0ec2a0b (the addendum in c7d3cf2); method versions: deadtime 6 (dark unchanged at 7).

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| R1 | wrong | Falsely identified dark dead time (SNSPD: 7/200, 0/7 cover) reaches corrected HOM unrefused | ADOPT-FIX | 0ec2a0b |
| R2 | questionable | Corrected-value labels printed on the raw line, one flat `labels` list | ADOPT-FIX | 0ec2a0b |
| R3 | questionable | Forward-model HOM/fringe not marked secondary | ADOPT-FIX | 0ec2a0b |
| R4 | questionable | "Calibration uncertainty" does not say what it excludes; prints 0.0000 | ADOPT-FIX + ADOPT-LIMIT | 0ec2a0b; limit (below) |
| R5 | questionable | README insufficient for a methods section | ADOPT-FIX (docs) | 0ec2a0b |
| R6 | minor | `--heralding` writes A+/B+ only; four pairs needed for CHSH capture | ADOPT-LIMIT + DEFER (`--pair`, M10) | limit (below); deferred (plan) |
| R7 | minor | HOM refusals reuse CHSH wording; Python-only remedy | ADOPT-FIX (text) | 0ec2a0b |
| R8 | minor | HOM prints no selection caveat | ADOPT-FIX | 0ec2a0b |
| T1 | questionable | tau_np interval zero width; misses an off-grid truth every time | ADOPT-FIX | 0ec2a0b |
| T2 | questionable | Dark artifact at other bias/temperature used silently by `calibrate deadtime` | ADOPT-FIX | 0ec2a0b |
| T3 | questionable | No readable uncertainty budget for tau_p | ADOPT-FIX | 0ec2a0b |
| T4 | minor | "likelihood ratio" printed for 2 ln LR | ADOPT-FIX | 0ec2a0b |
| T5 | minor | tau_np definition mixes quantity and estimator | ADOPT-FIX (text) | 0ec2a0b |
| T6 | minor | Validity range not machine-readable | ADOPT-FIX (record) + DEFER (check, M5) | 0ec2a0b; deferred (plan) |
| T7 | minor | Certificate lines: no level, too many digits, limit level | ADOPT-FIX (a, b) + ADOPT-LIMIT (c) | 0ec2a0b; limit (below) |
| T8 | minor | "input rate ... cps" is the avalanche-candidate rate; `_sd_held` unit | ADOPT-FIX | 0ec2a0b |
| B1 | questionable | M5's named NUTS users have no likelihood | ADOPT-FIX (plan) | 0ec2a0b |
| B2 | questionable | Multimodality invisible to R-hat from MLE starts | ADOPT-FIX (plan) | 0ec2a0b |
| B3 | questionable | Hard boundaries / support parameters vs "any divergence" | ADOPT-FIX (plan) | 0ec2a0b |
| B4 | questionable | Stored draws imply flat priors in fit coordinates (improper in logit p) | ADOPT-LIMIT + DEFER (M5) | limit (below); deferred (plan) |
| B5 | questionable | Series tau_p draws lose correlation (about 0.3) with dark afterpulse draws | DEFER (M5 entry) | deferred (plan) |
| B6 | minor | Union intervals drawn as one split normal (18% of draws in the gap) | ADOPT-FIX | 0ec2a0b |
| B7 | minor | "Stored draws gated" (M1) not done; no draw-level test | DEFER (M5) + plan fix | deferred (plan) |
| B8 | minor | Benchmark speed summary selective | ADOPT-FIX (docs) | 0ec2a0b |
| R1-a' | - | The lower-end floor refuses nominal 5% misses per channel (15-21% of chsh-dead-time replicates) | ADOPT-FIX: withdraw R1(a)'s refusal | c7d3cf2 |
| R1-est | - | Change the dark estimator (do not identify here, or widen) | REJECT | not adopted |
| R1-regime | - | Keep the refusal and change the regime (longer dark runs) | REJECT | not adopted |

Limits of round 6 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | README "Calibration artifacts", near the dead-time text |
| L2 | README "Raw values ... Corrected values" paragraph, and the docstring of `corrected_visibility` / `corrected_chsh` |
| L3 | dark.py is pinned: README dead-time section and the `deadtime.py` docstring "Fit" |
| L4 | README "For a CHSH setup", and the `calibrate pairs` / quickstart printed hint |
| L5 | `deadtime.py` docstring "Intervals" and the unidentified note |
| L6 | README "Calibration draws", and `docs/phase3-plan.md` M5 |

## Round 7: time tagger, security, test quality

Lenses: time-tagger engineer (G), security (S), test quality (Q). Ruling: [phase3/r7-council.md](phase3/r7-council.md); reports: [r7-tag.md](phase3/r7-tag.md), [r7-sec.md](phase3/r7-sec.md), [r7-test.md](phase3/r7-test.md).

Fixes in 8d7c50a; method versions: deadtime 7.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| G1 | wrong | T3 tick lost; round-6 tau_np widening is 1 ps (10/20 cover at res 25) | ADOPT-FIX | 8d7c50a |
| G2 | wrong (conditional on the header convention in #1) | T3 `MeasDesc_GlobalResolution` checked at rtol 1e-9; `%g` hides the difference | ADOPT-FIX | 8d7c50a |
| G3 | wrong | Pulsed T3 HOM: side-window accidentals low by about W/period; corrected V 0.962 vs 1, no label | ADOPT-FIX + ADOPT-LIMIT + DEFER (pulsed side peaks, v0.2+ as plan:410) | 8d7c50a; limit (below); deferred (plan) |
| G4 | questionable | No check that 2^bits x res covers the T3 sync period | ADOPT-FIX | 8d7c50a |
| G5 | questionable | T3 sync divider convention unverified | ADOPT-FIX (guard) + DEFER (verification, M9) | 8d7c50a; deferred (plan) |
| G6 | questionable | Tagger model and serial not checked at use time | ADOPT-FIX | 8d7c50a |
| S1 | questionable | Gate merge passes on empty results; band read from the shard file | ADOPT-FIX | 8d7c50a |
| S2 | minor | Crafted `.npy` header: MemoryError escapes `ArtifactError` | ADOPT-FIX | 8d7c50a |
| S3 | minor | PTU memory figure holds for T2 only | ADOPT-FIX (docs) | 8d7c50a |
| S4 | minor | `ppl-bench.yml` checkouts persist the token | ADOPT-FIX | 8d7c50a |
| S5 | minor | Committed Mach-O binaries `bench/ppl/models/{dark,delay}` | ADOPT-FIX (remove in a normal commit) + OWNER (history) | 8d7c50a; owner item |
| S6 | minor | A tag on a main commit with no completed gates can release | OWNER | owner item |
| Q1 | questionable | Gate pass/fail logic and coverage-merge tiling untested (10 survivors) | ADOPT-FIX (tests) | 8d7c50a |
| Q2 | questionable | Channel order and declared-delay sign untested (V 0.910 -> 0.835 unnoticed) | ADOPT-FIX (tests) | 8d7c50a |
| Q3 | questionable | Wired labels and refusals never observed through the analysis path (10 survivors) | ADOPT-FIX (tests) | 8d7c50a |
| Q4 | questionable | Two tests pass only at their seed (5/20, 69/80) | ADOPT-FIX (tests) | 8d7c50a |
| Q5 | questionable | `test_recovery_at_three_points`: three 95% intervals at one seed (15/20); 5 ns point ungated | ADOPT-LIMIT + DEFER (M10) | limit (below); deferred (plan) |
| Q6 | questionable | Deadtime interval logic defended by a golden test whose oracle is the code | ADOPT-FIX (tests) | 8d7c50a |
| Q7 | minor | Forward-check statistics unprotected below the gate | ADOPT-FIX (tests) | 8d7c50a |
| Q8 | minor | `pytest.skip` on a seed premise; `_dark_artifact_for` channel/multiple untested | ADOPT-FIX (tests) | 8d7c50a |
| Q9 | minor | Slow tests in the default tier; `slow` marker never deselected | ADOPT-FIX (CI) | 8d7c50a |

Limits of round 7 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | README line 18 after "PicoQuant PTU (T2 and T3) input", and the `adapters/ptu.py` module docstring T3 paragraph |
| L2 | adapters/ptu.py:54-56 |
| L3 | docstring of `test_recovery_at_three_points`, `test_tau_p_never_exceeds_tau_np`, `test_afterpulse_is_held_at_the_dark_artifact` |

## Round 8: physics, statistical design and the plan as a whole

Lenses: physics as a whole (Ph), statistical design as a whole (St), plan (Pl). Ruling: [phase3/r8-council.md](phase3/r8-council.md); reports: [r8-phys.md](phase3/r8-phys.md), [r8-stat.md](phase3/r8-stat.md), [r8-plan.md](phase3/r8-plan.md).

Fixes in 2e15eba; method versions: delay 8 -> 9.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| Ph1 / St3 | wrong | T3 tick (res + 1) used as the histogram lattice: delay -12 ps, jitter 0/10 covered, every artifact a misfit | ADOPT-FIX | 2e15eba |
| Ph2 | wrong (open item from round 7) | Forward check (and `pooled_peak`) misfire on T3: bin on the dtime lattice | ADOPT-FIX | 2e15eba |
| Ph3 | wrong | Pulsed light on T2 or plain tags: V 10-40 sd off, unlabelled | ADOPT-FIX + ADOPT-LIMIT | 2e15eba; limit (below) |
| Ph4 | wrong | Dark dead times above the analysed shortest gap (400/400 replicates) | ADOPT-LIMIT + DEFER (M5a, together with Ph5 and a redesigned control regime) | limit (below); deferred (plan) |
| Ph5 | questionable | Pair live fraction under-corrects at second order | ADOPT-LIMIT + DEFER (M5a, only together with Ph4) | limit (below); deferred (plan) |
| Ph6 | questionable | Gate regimes narrower than real labs | DEFER (M4b / M3b regime sitting, owner) | deferred (plan); owner item |
| Ph7 | minor | The 20% dead-fraction refusal gives the wrong reason | ADOPT-FIX (wording) | 2e15eba |
| St1 | wrong | Downward gap floor on the wrong scale (2.5e-4 to 3.7e-4 per channel) | ADOPT-FIX | 2e15eba |
| St2 | wrong | The gate aborts on any refusal; the CI pass is seed luck | ADOPT-FIX (record refusals by type, verdict unchanged) + OWNER (tolerance, second seed) | 2e15eba; owner item |
| St4 | questionable | No joint budget; labels and refusals at mixed alphas | ADOPT-FIX | 2e15eba |
| St5 | questionable | `delay.py` dof groups-5 runs hot (9.5% at p<0.05) | ADOPT-FIX (method 8 -> 9) | 2e15eba |
| St6 | minor | Upward gap test's "1e-4 per channel" wrong | ADOPT-FIX | 2e15eba |
| St7 | minor | Drift and predictive-check events share a label and a key | ADOPT-FIX | 2e15eba |
| Pl1 | questionable | M5's only consumer of the series artifact is optional | ADOPT-FIX (plan text: M5a / M5b; neutral strings) + OWNER (M5a in the gate minimum) | 2e15eba; owner item |
| Pl2 | questionable | M10 overloaded | ADOPT-FIX (plan text: M2b, M3b, M4b, M10) + OWNER (which may slip) | 2e15eba; owner item |
| Pl3 | questionable | Homeless deferrals (pair labels, D8, site) | ADOPT-FIX (plan text) | 2e15eba |
| Pl4 | questionable | Owner items scattered and dropped | ADOPT-FIX (plan section) | 2e15eba |
| Pl5 | questionable | Slow-test job not blocking | ADOPT-FIX (`gates` needs) + OWNER (required check) | 2e15eba; owner item |
| Pl6 | questionable | PR #26 description stale; two PRs under one title | ADOPT-FIX (title and body) + OWNER (merge method) | 2e15eba; owner item |
| Pl7 | minor | Decision texts do not match what was built | ADOPT-FIX (D3 as-built note, lines 8 and 200) + OWNER (D4 text) | 2e15eba; owner item |
| Pl8 | minor | Dependency table drift | ADOPT-FIX | 2e15eba |

Limits of round 8 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | README, new "Refusals and their rates" section; `guardrails.py` module docstring |
| L2 | README line 18 limit and `adapters/ptu.py` T3 paragraph, extending round 7's L1 |
| L3 | README dead-time limits; `forward.py` docstrings of `live_fraction` and `_channel_dead_times` |

## Round 9: red team 2, hardware twin, performance

Lenses: red team 2 (Rd), hardware twin (Tw), performance (Pf). Ruling: [phase3/r9-council.md](phase3/r9-council.md); reports: [r9-red.md](phase3/r9-red.md), [r9-twin.md](phase3/r9-twin.md), [r9-perf.md](phase3/r9-perf.md).

Fixes in e8c1aba; method versions: delay 10; deadtime 7 re-pinned (outputs byte-identical).

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| Rd1 | wrong | One echo turns off the upward gap test; S 2.029 for 1.980, unrefused | ADOPT-FIX (council form, not the red team's as written) + ADOPT-LIMIT | e8c1aba; limit (below) |
| Rd1b | wrong | No analysis-time cross-check of the dark dead time against the setup's declared `dead_time_ps`; D3 says "the truth lies between" | ADOPT-FIX (label) | e8c1aba |
| Rd2 | questionable | `calibration_outside_validity` masks the gap refusals; reason dropped | ADOPT-FIX (own override name, problems recorded) | e8c1aba |
| Rd3 | questionable | `spdc-cw` counts as the T3 clock declaration; templates bypass G3 | ADOPT-FIX + ADOPT-LIMIT (README sentence) | e8c1aba; limit (below) |
| Rd4 | minor | Unrecognised condition names recorded but never compared | ADOPT-FIX (label, help text) | e8c1aba |
| Rd5 | minor | Forward check skipped on T3, note not printed | ADOPT-FIX | e8c1aba |
| Tw1 | wrong | T3 delay biased by pairs straddling a sync edge, status ok | ADOPT-FIX (delay "9" -> "10") | e8c1aba |
| Tw2 | questionable | T2 dark artifacts refused on T3 analyses (shared key) | ADOPT-FIX | e8c1aba |
| Tw3 | questionable | A CW T3 clock change refuses every artifact | ADOPT-FIX (label under the explicit `cw` declaration) | e8c1aba |
| Tw4 | questionable | Refusals point CLI users to `Policy.allow`; override names and the 0.5 K rule undocumented | ADOPT-FIX | e8c1aba |
| Tw5 | questionable | SNSPD refusal suggests a remedy no command delivers | ADOPT-FIX (wording) | e8c1aba |
| Tw6 | minor | One dark PTU file cannot be calibrated as given | ADOPT-FIX | e8c1aba |
| Tw7 | minor | README calls serials optional | ADOPT-FIX | e8c1aba |
| Tw8 | minor | First-run path not signposted; NO_VALIDITY repeated 8x | ADOPT-FIX | e8c1aba |
| Pf1 | wrong | Off-peak dispersion half cannot refuse in CHSH; doubled cost | ADOPT-FIX (sequential MC, counts-seeded) + ADOPT-LIMIT (residual power) | e8c1aba; limit (below) |
| Pf1(b) | wrong | Null-statistic cache in `accidentals.py` | REJECT (superseded by Fix 1; it would re-pin accidentals for no gain) | not adopted |
| Pf2 | questionable | `occupied_slots` uses `np.unique` on sorted input | ADOPT-FIX (from `perf.patch`) | e8c1aba |
| Pf3 | questionable | Repeated sorts and copies; 3x peak memory | ADOPT-FIX (from `perf.patch`) | e8c1aba |
| Pf4 | questionable | Real datasets do not fit in memory; not stated | ADOPT-LIMIT + DEFER (M9, chunked analyses) | limit (below); deferred (plan) |
| Pf5 | questionable | chsh-spad shards near timeout; unit jobs without timeouts | ADOPT-FIX (timeouts; Fix 1 and Fix 9 cut the cost) | e8c1aba |
| Pf6 | minor | Series dead-time constants recomputed per call | ADOPT-FIX (re-pin deadtime sha256, version stays "7", output-unchanged note) | e8c1aba |
| Pf7 | minor | Peak fit uses finite-difference gradients | DEFER (M10 performance pass; not bit-identical, needs the stated acceptance test) | deferred (plan) |
| Pf8 | minor | Repeated full-array passes | ADOPT-FIX (single `np.quantile`, memoised grid counts: bit-identical) + DEFER (difference extraction to M9) | e8c1aba; deferred (plan) |

Limits of round 9 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | README dead-time limits; `forward.check_dead_time_gaps` docstring |
| L2 | README T3 paragraph, replacing the sentence at README:94-98, extending round 8's L2 |
| L3 | README Known limits; CHANGELOG |
| L4 | README, after round 8's L2; `_twopoint.offpeak_check` docstring |

## Round 10: regressions, first-time user, merge readiness

Lenses: interactions and regressions (Rg), first-time user (Ux), merge readiness (Mg). Ruling: [phase3/r10-council.md](phase3/r10-council.md); reports: [r10-reg.md](phase3/r10-reg.md), [r10-ux.md](phase3/r10-ux.md), [r10-merge.md](phase3/r10-merge.md).

Fixes in the round-10 commit; method versions: deadtime 8; delay stays 10.

| id | severity | finding | ruling | where fixed |
| --- | --- | --- | --- | --- |
| Rg1 | questionable | tau_p cap ends at m while tau_np is widened a tick; coverage 7-8 / 60 on coarse ticks, fully paralyzable | ADOPT-FIX (variant B) + ADOPT-LIMIT | the round-10 commit; limit (below) |
| Rg2 | questionable | `DEAD_TIME_ABOVE_GAPS` breaks the label rule; duplicates `DEAD_TIME_BELOW_DARK` | ADOPT-FIX | the round-10 commit |
| Rg3 | minor | `dark_*` keys will be labelled "not compared" once M5a loads series artifacts | DEFER (M5a entry items) | deferred (plan) |
| Rg4 | minor | 100-replicate quick gate cannot judge the chsh-dead-time control rule | ADOPT-FIX (describe text) + ADOPT-LIMIT; `passed` unchanged (a gate rule, O12) | the round-10 commit; limit (below) |
| Rg5 | minor | "calibration 4.4e-16"; identified bound printed [0, ...] | ADOPT-FIX | the round-10 commit |
| Ux1 | wrong (misleading in realistic use) | analyze never prints the overrides or reasons | ADOPT-FIX | the round-10 commit |
| Ux2 | wrong (misquote risk) | raw value labelled "policy modified" for correction-only overrides | ADOPT-FIX | the round-10 commit |
| Ux3 | questionable | CHSH dead end on missing delay-jitter artifacts | ADOPT-FIX | the round-10 commit |
| Ux4 | questionable | policy refusals one per run | ADOPT-FIX (one pre-flight, same check function) | the round-10 commit |
| Ux5 | questionable | condition refusal offers only the override | ADOPT-FIX | the round-10 commit |
| Ux6 | questionable | flat label list, no severity order | ADOPT-FIX (one per line, overrides first, fixed excludes once) + DEFER (severity grouping, M10) | the round-10 commit; deferred (plan) |
| Ux7 | questionable | no label when the predictive check did not run | ADOPT-FIX | the round-10 commit |
| Ux8 | minor | README cites a constant name the output never shows | ADOPT-FIX (README quotes opening words) + DEFER (printed codes, M10) | the round-10 commit; deferred (plan) |
| Ux9 | minor | no single label and refusal table; `--allow` help lists one name | ADOPT-FIX (table, help list, sync test) | the round-10 commit |
| Ux10 | minor | "dead fraction of 2.66" | ADOPT-FIX | the round-10 commit |
| Ux11 | minor | corrected "±" unlabelled; calibration part reads as total | ADOPT-FIX | the round-10 commit |
| Ux12 | minor | HOM without dark artifacts: no "dead time not modelled" label | ADOPT-FIX | the round-10 commit |
| Ux13 | minor | error above header in redirected logs | ADOPT-FIX | the round-10 commit |
| Mg1 | blocker | replicate hash fails on Linux; required check red | ADOPT-FIX (fixture with tolerance; not the rounded hash) | the round-10 commit |
| Mg2 | blocker (process) | "Statistical gates" not a required check | Owner (O1, moved before merge) | owner item |
| Mg3 | must decide | linear history: squash or rebase | Owner (O14 reworded) | owner item |
| Mg4 | should, before merge | council rulings only in a scratch directory | ADOPT-FIX (record designed below) | the round-10 commit |
| Mg5 | should | PR title and body stale | ADOPT-FIX (draft, updated for round 10) | PR title and body, by the merge auditor (draft in r10-merge.md) |
| Mg6 | minor | no consolidated behaviour-changes list | ADOPT-FIX | the round-10 commit |
| Mg7 | minor | artifact compatibility of method bumps unstated | ADOPT-FIX + ADOPT-LIMIT | the round-10 commit; limit (below) |
| Mg8 | minor | "this branch starts it" | ADOPT-FIX | the round-10 commit |
| Mg9 | minor | CONTRIBUTING omits the byte-identical re-pin | ADOPT-FIX | the round-10 commit |
| Mg10 | minor | python/README omits L4 | ADOPT-FIX | the round-10 commit |
| Mg11 | minor | plan header "rounds 4-8"; O7 status stale | ADOPT-FIX | the round-10 commit |

Limits of round 10 (the wording is in the ruling; it lives where the ruling put it):

| id | where the wording lives |
| --- | --- |
| L1 | CONTRIBUTING, Gates; `propagation.describe` for control regimes below 400 replicates |
| L2 | README dead-time limits, after round 1's L3; `calibrate_deadtime_series` docstring |
| L3 | CHANGELOG, top of Unreleased; README calibration artifacts |
| L4 | python/README.md, the off-peak sentence |

## Owner items

Consolidated O1-O21 as the round-10 council left them (also in `docs/phase3-plan.md`, "Owner items").
**Before merge**

| # | Item | Source | Council recommendation |
|---|---|---|---|
| O14 | Merge method: linear history is required, so squash or rebase-merge (a merge commit needs a ruleset bypass). Review cap: this round is the last | r8 Pl6, r10 Mg3 | Rebase-merge, so each round's commit stays traceable from the record (Fix 7). Squash is acceptable now that the record carries the ids |
| O1 | Make "Statistical gates" (and the slow-test job it needs) required on `main` (ruleset 27075146). Until then, check the gates by hand | r4 C1, r8 Pl5, r10 Mg2 | Do it before merging #26, not after: cheap, and it protects this merge |
| O19 | Publishing the review record in the repo: the AI reviewer reports and rulings verbatim (about 600 KB, paths scrubbed) under `docs/reviews/phase3/` | r10 Mg4 | Publish, as phase 2 did; the summary states that reviewers are AI agents |
| O20 | Only if Fix 1's Linux run deviates by more than 1e-8 relative: accept per-platform fixtures (macOS and Linux, from a run whose gates are green) instead of one tolerance | r10 Mg1 | Accept, with the observed deviation recorded in the test |

**After merge, in this order**

| # | Item | Source | Council recommendation |
|---|---|---|---|
| O3 | Dev version 0.2.0.dev0 (still 0.1.1; delay at 10, deadtime at 8) | r5 V2, r6, r7 | First, before more artifacts are written from main |
| O4 | Calibration schema `$id` | r5 V7 | With O3 |
| O2 | Release requires gates on the tagged SHA; forbid direct pushes to main | r7 S6 | In `release.yml`'s `check` job |
| O8 | Pre-register gate and regime changes in one sitting (rotating seeds, tau_np gate, D3 union, past-maximum, unidentified tau_p, 5/22 ns, realistic regimes, a control regime that survives M5a); add a coarse-tick (T3) series regime, so Rg1's case is gated | r1-r8, r10 Rg1 | One sitting before any M4b or M5a run |
| O12 | Propagation gate tolerance: at most 5 refused per 1,000, binomial misfit and flag caps; the quick-run control rule (Rg4) | r8 St2, St5, r10 Rg4 | Adopt (a gate rule) |
| O6 | M5 split into M5a / M5b; M5a in the gate minimum; Rg3 is an M5a entry item | r6 B1, r8 Pl1, r10 Rg3 | Split; M5a in the gate minimum |
| O13 | Which of M2b, M3b, M4b may slip to 0.2.x | r8 Pl2 | M4b gates with M5a; M3b may slip |
| O16 | M9 extended to chunked analyses (Pf4, Pf8) | r9 Pf4, Pf8 | Adopt into M9; README limit until then |
| O17 | Pf7 analytic gradient in M10, with its acceptance run as entry test | r9 Pf7 | Schedule in M10 |
| O21 | M10 output pass: label severity groups (Ux6) and printed label codes (Ux8) | r10 Ux6, Ux8 | Schedule in M10 |
| O18 | CI budget: buy concurrency, or run the propagation gate nightly and on release | r9 Pf5 | Owner's choice |
| O15 | Templates keep `spdc-cw` or ship `cw` | r9 Rd3 | Keep `spdc-cw` |
| O10 | A real divided T3 file | r7 G5 | Source it, or ship the refusal |
| O5 | D4 decision text ("0.1.x readers refuse") | r5 V6 | Reword as r5 L2 |
| O7 | Accept every behaviour change in the release notes (the CHANGELOG list from Fix 8, round 10 included) | r6-r10 | Accept at release |
| O9 | Purge the bench binaries from history | r7 S5 | Leave (default) |
| O11 | Withdrawn deferral (pending afterpulse into M4) | r5 item 4 | Close |
