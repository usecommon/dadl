> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 4 (dead-time statistics, research software engineering, CI)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: DT = dead-time statistician (r4-dt.md, rows 1-5 as DT1-DT5), S = research software engineer
(r4-rse.md, rows 1-8 as S1-S8), C = CI engineer (r4-ci.md, C1-C8). Council probes and the prototype:
`scratchpad/council/r4-council/` (`proto/` is a `git archive HEAD` copy of `python/` with the branch rule
applied; `branch_rule.diff` is the change; `w1.py` the end-to-end probe, `gaprule.py` the gap-check probe,
`prof.py` the timing; logs `e2e_*.log`, `cov_*.log`, JSON `e2e_*.json`).

## Checked numerically

All runs on branch m4/series-deadtime at 926560b with the repo venv, through the package's own functions,
one BLAS thread. "Prototype" = HEAD plus `branch_rule.diff` (deadtime.py only).

- **DT1 (refit jumps branch): confirmed, and fixed by the prototype.** `w1.py pq 0,20,24,25` at HEAD
  reproduces the reviewer's numbers exactly: afterpulse shift 410 ps, **1.36e4, 1.53e4, 8.4e3 ps** (intervals
  [16.0, 45.0], [8.1, 45.0], [21.3, 45.0] ns). Prototype (refits held on each run's decided branch): 410,
  **412, 439, 309 ps**. Over 160 pq replicates the largest shift is 512 ps (HEAD: up to 1.5e4; s30 at HEAD up
  to 6.25e3, prototype 200 ps). Replicate 20 is then not covered ([40.43, 45.0] ns, truth 40): it was covered
  at HEAD only through the spurious 13.6 ns sd, and is one of the expected 5%.
- **DT2 + DT3 (both-branch search on a below-maximum scan; joined interval): confirmed** from the reviewer's
  stored replicates (`r4-dt/summ.py`): 8/160 pq and 7/96 s30 global MLEs put the top run at about 9e7-1e8
  cps; the held-minus-global NLL is 0 in all other null runs and at most **2.05 (pq), 2.06 (s30)**; past the
  maximum (1/3/60 Mcps) it is 0, 0.89, 1.44, 1.47, 1.58, 1.99, then 5.4 median, up to 13.4. Asimov (reviewer's
  `asimov.py`) has the high branch worse by 5.6 and 2.5 NLL in expectation. The join's p90 widths 12.9 (pq)
  and 5.4 ns (s30) against 4.7 and 2.2 ns held: confirmed.
- **The adopted branch rule (Fix 1), verified end to end** through `calibrate_deadtime_series` (afterpulse
  widening and tau_np cap included; `w1.py`, seeds [777, i, k] as the reviewer's):

  | scan | HEAD covered, median / p90 width | prototype covered, median / p90 width | prototype decisions |
  |---|---|---|---|
  | pq 1/3/6 Mcps, tau_p 40 ns, 160 reps | (4 reps only: see DT1) | **153/160**, 4.73 / **5.07 ns** | 152 low, 8 union (8 covered), 0 switched |
  | s30 1/3/10 Mcps, tau_p 30 ns, 96 reps | 90/96, 2.26 / **15.2 ns** | **90/96**, 2.22 / **2.41 ns** | 89 low, 7 union (7 covered) |
  | past 1/3/60 Mcps, tau_p 40 ns, 48 reps | 48/48, 6.34 / 7.98 ns | **47/48**, 3.53 / 6.49 ns | 36 switched (36 covered), 11 union (11), 1 low (missed: global MLE on the low branch too) |

  Point means: pq 39.65 ns (HEAD global 39.21), s30 29.79 (29.53). The false "afterpulse still pending" note
  (W2) on pq: **0/160** (reviewer: 1.5% at HEAD); s30 7/96 against 14/96 at HEAD (its true pending fraction,
  0.46%, sits just under the 0.5% threshold); past 36/48 (true: the 6e7 run has pending afterpulses).
- **CI regimes with the prototype (coverage module, CI seed 20261111, `DADL_BACKEND=python`):**
  series-spad 1,000 replicates: tau_p **95.5%** (gate 92.4-97.2%), misfit 0.3%, tau_np 99.2%: ok (p = 0, so
  the rule does not run there; unchanged). passive-quench, first 300 replicates (`cov_pq_proto.log`, `cov_pq_head.log`): tau_p **94.7%** (HEAD 95.7%;
  gate at 300 replicates 90.0-98.7%; the 3 replicates of difference are HEAD covering through the DT1
  inflation), misfit 0.0%, tau_np 98.0% (both): ok. CPU 3,638 s against 4,995 s at HEAD (12.1 vs 16.7 s a
  replicate, 1.37x faster).
- **Speed.** One pq calibration (`prof.py`, one core): HEAD 14.5 s CPU (fit 9.2, interval 3.3, refits 1.9);
  prototype 10.1 s (fit 7.9, interval 1.2, refits 0.9): about 1.4x faster (gate run above: 1.37x), not the round-2 speed (about 4.2 s). The
  remaining cost is the both-branch grid search (7.7 s of which 1.9 s is the held fit it now reuses). A
  per-run flip search (Brent with one run moved to the other branch, 1.8 s in place of 5.7 s; `prof2.py`)
  would restore most of it but only tests single-run flips: not adopted now (DEFER below).
- **DT4 (gap check): confirmed, and the replacement rule checked** (`gaprule.py`). Null, 150 dark
  calibrations per detector x 4 analysed rates (2e4-5e6 cps; plain 50 kcps 1 s, with 2% afterpulsing,
  0.1 s, and a series tau_p 40 / tau_np 45 ns detector): threshold at the 3-width floor itself (no 0.9):
  **0/2,400 refused**; upward test (below) **0/2,400 refused**, largest excess 7.0 (series, 5 Mcps; the
  statistic is Exp(1), so 9 gives about 1.2e-4 per channel). Misses: calibrated 50.7 ns on a 46 ns detector
  is now refused downward (1,122 gaps below the floor; passed with 0.9); 22.2 ns on a 50 ns detector at 1 Mcps
  refused upward (excess 1,339); 41.7 on 50 ns refused (431); 22.1 on 50 ns at 2e4 cps not (1,000 gaps,
  excess 0.8: no power at low rates, stated).
- **S1, S2, S4-S8: confirmed by rerunning the reviewer's probes** at HEAD: `p5` (CHSH refused `Misspecified`
  with an empty calibration folder, p = 2.2e-9); `p1` (`KeyError` from `corrected_visibility` and
  `bound_range`; CHSH gives the named `ValueError`; duplicate timestamp -> "Geometric sequence cannot include
  zero"; another channel's dark artifact accepted); `p2` (API `hom_from_tags` with another detector's
  artifacts accepted, `load(channels=)` refuses; `counts['subtracted'] = 0.0` with accidental subtraction);
  `p4` (`KeyError 'from'`); `p6` (floor -60.8 ns, check off). `cli._grid_capture(500, 20, 10000)` =
  2.4e-12 against `forward.grid_capture` 0.95; equal at (0, 100, 1000). Schema text "Recorded only: DADL does
  not yet enforce it" at calibration-1.schema.json:90. S4 mutations accepted on the reviewer's evidence
  (`mutations*.txt`), not rerun.
- **C1: confirmed** (`gh api .../rulesets/27075146`: 11 required contexts, none a gate; no `if: always()` on
  the merge jobs, ci.yml:249, :339). **C2, C3, C5, C6, C8: confirmed by code** (chsh-spad unsharded at
  ci.yml:283 with timeout 45; `CHUNK = 50` fixed at propagation.py:293/330; no `concurrency:`). Timings
  accepted from the reviewer's logs, not rerun. Propagation regimes all calibrate the delay (propagation.py:
  202), so S1's guard does not change gate results.
- **Method versions.** Fix 1 and Fix 6 edit `calibration/deadtime.py` (pinned, now "3"): **bump to "4"**
  (`METHOD_VERSION`, `tests/method_versions.json` version and sha256). `forward.py`, `chsh.py`, `cli.py`,
  `propagate.py` are not pinned. No other pinned module changes.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| DT1 | Afterpulse-sensitivity refit searches both branches and lands on the other minimum (shift sd up to 1.5e4 ps; "identified" interval becomes most of [0, tau_np]) | Y/Y/Y | ADOPT-FIX | Fix 1(c). Verified: max 512 ps over 160 pq |
| DT2 | Both-branch MLE puts a below-maximum scan on the high branch in about 5% (tau_p about 31 ns, rate 9e7, false pending note) | Y/Y/Y | ADOPT-FIX | Fix 1 (branch decided once; point, rates and notes from the decided fit). Verified 0/160 false pending notes |
| DT3 | Joined interval is a patch: widens on both sides (p90 12.9 vs 4.7 ns), its note fires exactly when the held fit is rejected | Y/Y/Y | ADOPT-FIX | Fix 1. User declaration of runs past the maximum: DEFER to M10 with the past-maximum coverage regime (data decide correctly in 47/48 with the margin; a declaration needs a CLI/schema field and its own regime) |
| DT4a | Gap check's 0.9 factor has no basis: dead times up to about 16% too long pass | Y/Y/Y | ADOPT-FIX | Fix 6(a): drop `GAP_FRACTION`; the 3-width floor alone (0/2,400 null refusals). State the remaining blind zone (L1) |
| DT4b | A calibrated dead time shorter than the detector's is never checked (22 ns artifact on a 50 ns detector passes, +2.7% live fraction) | Y/Y/Y | ADOPT-FIX | Fix 6(b): upward test, also for bounds. 0/2,400 null; refuses the reviewer's case |
| DT5 | Docstring and CHANGELOG justify the join by "high branch can fit better by up to a few NLL" (Asimov: worse by 2.5-5.6) | Y/Y/Y | ADOPT-FIX | Fix 1 (rewrite) and Fix 9 (CHANGELOG) |
| S1 | CHSH runs the shape check with no delay-jitter artifact; HOM/fringe skip it; refuses accidentals-only corrected S | Y/Y/Y | ADOPT-FIX | Fix 2 |
| S2 | Partial dark set in HOM/fringe gives a bare `KeyError` (CLI traceback) | Y/Y/Y | ADOPT-FIX | Fix 3 |
| S3 | Detector-serial check (round 3 R1) runs only through `CalibrationSet.load(channels=)`; API paths skip it | Y/Y/Y | ADOPT-FIX | Fix 4. Propagation re-run |
| S4 | Tests survive removal of the CHSH shape check, the gap check in all three paths, and `capture_sd_data` | Y/Y/Y | ADOPT-FIX | Fix 5 (tests) |
| S5 | `cli._grid_capture` (quad) returns about 0 for a narrow peak in a wide window: false "window misses the pair peak", pair rate "high by 1/capture" | Y/Y/Y | ADOPT-FIX | Fix 7(a) |
| S6 | Gap check off without record when the floor is <= 0 | Y/Y/N (counts) ; N/N/Y (label) | ADOPT-FIX (counts only) | Fix 6(c). No label: the upward test still runs on such artifacts, and a label would fire on every short dark run |
| S7 | `calibrate_deadtime_series`: duplicate timestamp gives a numpy error; same run twice accepted; afterpulse artifact kind/channel/misfit unchecked | Y/Y/Y | ADOPT-FIX | Fix 1(e) (inside the deadtime version 4 change) |
| S8 | Counts keys and docs drift: `subtracted` 0.0 with subtraction; counts keys undocumented; CLI help "HOM: corrected" only; schema `valid_for` "not enforced" and `KeyError 'from'`; CHSH quickstart's fixed dark note | Y/Y/Y | ADOPT-FIX | Fix 7(b)-(f), Fix 9 |
| C1 | No gate is a required check; a failed shard makes the merge job skipped (counts as passing) | Y/Y/Y | ADOPT-FIX + owner item | Fix 8(a) aggregator `gates` job; owner: make `gates` required in the `main` ruleset |
| C2 | chsh-spad propagation estimated 27-41 min of 45 | Y/Y/Y | ADOPT-FIX | Fix 8(b): through the shard/merge path, 4 shards |
| C3 | Propagation shards leave a worker idle (chunks 50/50/25) | Y/Y/Y | ADOPT-FIX | Fix 8(c) |
| C4 | Coverage gate prints nothing until done | Y/Y/Y | ADOPT-FIX | Fix 8(d) |
| C5 | Comments understate passive-quench cost (4.6x, 19 s a replicate, not 3x / 8 s) | Y/Y/Y | ADOPT-FIX | Fix 8(e): re-measure after Fix 1 and write the measured number |
| C6 | CI-structure test checks only coverage regimes | Y/Y/Y | ADOPT-FIX | Fix 8(f) |
| C7 | CONTRIBUTING has no gate instructions | Y/Y/Y | ADOPT-FIX | Fix 9 |
| C8 | No `concurrency:` group; about 51 jobs per push queue | Y/Y/Y | ADOPT-FIX | Fix 8(g), cancel in progress for pull requests only |
| (council) | Per-run flip search in place of the both-branch grid (about 2.5x faster gate) | N/Y/N (now) | DEFER | M10, with the past-maximum regime: it must be checked against multi-run flips there before it replaces the exact profile |

### Limit wordings

- **L1** (forward.py `check_dead_time_gaps` docstring; README use-time checks paragraph): "The check refuses
  a dark dead time that the analysed gaps contradict in either direction: gaps shorter than its floor (the
  point minus three lower half-widths, about 7% below the point for a 1 s dark run at 50 kcps), or a shortest
  analysed gap more than 9 / (edge density) above its upper end (chance about 1e-4 per channel). A dead time
  between the floor and the truth passes (live-fraction error up to about 1.5% at 5 Mcps for a 1 s dark run);
  at low analysed rates (about 1,000 gaps) the upward test has no power."
- **L2** (deadtime.py module docstring, "Two rates per gap law", replacing the round-3 text): "Each run is
  taken below the count-rate maximum (n tau_p < 1, the designed scan) unless the fit with each run on its
  better branch is better by more than `BRANCH_MARGIN` = 3 in NLL (a likelihood ratio of about 20); then
  those runs are fitted past it. Within the margin the interval spans both fits and a note names the runs.
  Below the maximum the high branch fits worse in expectation (Asimov: 2.5-5.6 NLL for the CI and CHANGELOG
  scans) and wins by at most 2.1 on 256 simulated scans; past it (1/3/60 Mcps) the margin is passed in 36 of
  48 and the union covers the rest but one. Refits for the afterpulse sensitivity hold each run's branch.
  Past the maximum, coverage is checked by probes only (no CI regime; M10)."
- **L3** (CHANGELOG round-4 entry): "deadtime method version 4: the round-3 joined interval is removed;
  intervals on below-maximum scans are narrower (pq p90 5.1 ns from 12.9; s30 2.4 from 15.2), tau_p and the
  fitted rates no longer jump to the high branch (about 5% of scans in round 3)."

## Implementation order (adopted fixes)

**Pinned method version: `deadtime` "3" -> "4"** (Fix 1; `METHOD_VERSION` and `tests/method_versions.json`).
**Gates to re-run:** coverage deadtime **passive-quench (all 12 shards) and series-spad** on the CI seed
after Fix 1; propagation **all five regimes** at 1,000 replicates on the CI seed after Fixes 2-6 (Fix 4 and
Fix 6 run in every corrected path; expect coverage and flag shares unchanged, no gap or detector refusals on
clean replicates). Then the CI timing comment (C5) from the measured shard times.

1. **Series branch rule (DT1, DT2, DT3, DT5, S7; deadtime version 4; prototype `r4-council/branch_rule.diff`).**
   (a) `calibrate_deadtime_series`: fit with every run on the low branch (`both_branches = False`), its
   profile interval and rates. When p > 0, also the both-branch fit, reusing the held minimum and nll0
   (`_fit_tau(..., held_min=, nll0=)`, which skips its internal held search). Let `past` = runs whose
   both-branch rate has n tau_p > 1 at its MLE, and gain = held NLL - both NLL. No `past`: keep the held fit.
   gain > `BRANCH_MARGIN` (3.0): use the both-branch fit (tau, interval, rates, lr), note "run(s) {past}
   fitted past the count-rate maximum (n tau_p > 1): that fit is better by {gain:.1f} in NLL". Otherwise:
   held point and rates, interval = union with the both-branch profile interval, note "run(s) {past} fit
   about as well past the count-rate maximum ({gain:.1f} in NLL, margin 3): the interval spans that fit (..)
   and the one below it (..)". Record `branch_nll_gain` in conditions. Remove the joined held interval.
   (b) Pending fraction, rate notes and tau_np's edge density use the decided rates (they already read
   `rates`). The p = 0 "either rate" note stays (only then is the both-branch re-evaluation run).
   (c) Sensitivity refits: `other.both_branches = False`, `other._n_cache = list(rates)` (decided rates).
   (d) Docstrings: L2 in the module docstring; `_fit_tau` and `best_n` docstrings say when each mode is used.
   (e) S7: refuse a zero gap ("duplicate timestamps in runs[i]", as `chsh.check_no_duplicates`); refuse two
   runs with the same `tags_sha256`; `afterpulse.kind == "dark-afterpulse-deadtime"` and `channel in
   afterpulse.channels` else `ValueError`; a misfit afterpulse artifact adds a note naming it (no policy on
   this path).
   Tests: pq seed `[777, 20, 1]` gives an afterpulse shift < 1 ns; pq seeds with a spurious high-branch
   global minimum (from `r4-dt/cov_pq.json`, e.g. the 8 with a 9e7 rate) give a low-branch point within 2 ns of
   40 and no pending note; the round-3 N2 past scan (1/3/60 Mcps) 2 seeds: interval covers 40 ns and the
   note names run 2; the S7 refusals. `METHOD_VERSION = "4"`, method_versions.json.
2. **CHSH shape check (S1).** `chsh._forward_check`: return `(None, ())` when `delay_start` gives no jitter for
   the pairs (`corrected_chsh` already needs all four or none); leave `forward_check_p` out of the counts then.
   Test: `r4-rse/p5` physics with accidentals only is not refused.
3. **Partial dark set (S2).** `forward._dead_times` and `bound_range` go through `_channel_dead_times` with
   `{"A": a, "B": b}` (one helper, the named `ValueError`). Test: `r4-rse/p1` cases raise `ValueError`; CLI
   prints the reason, no traceback.
4. **Detector check on API paths (S3).** In the corrected CHSH, HOM and fringe paths run
   `detector_mismatches` for each used member against the analysis `channels` (through `_outside`, the same
   `calibration_outside_validity` override and `OUTSIDE_VALIDITY` label). Test: `r4-rse/p2` via
   `hom_from_tags` and `chsh` refused; with the override, labelled.
5. **Tests for round-3 features (S4).** (a) `chsh()` with delay artifacts and a crosstalk-like extra peak
   raises `Misspecified`, and with the override carries `FORWARD_MISFIT_ACCEPTED`; (b) the 1 us-on-30 ns
   physics through `chsh()` and `hom_from_tags()` (StaleCalibration; `min_gap_ps` in counts), plus the
   upward case of Fix 6; (c) `capture_sd_data > 0` and smaller for a 4x longer run, and one case where
   3.29 sd_data > 0.5 se decides the label.
6. **Gap check (DT4, S6).** `forward.check_dead_time_gaps`: (a) remove `GAP_FRACTION`; count gaps below the
   floor. (b) Upward: with the analysed gaps sorted, when at least 101, f_edge = 100 / (g[100] - g[0]); refuse
   when (g[0] - tau.high) x f_edge > `GAP_EDGE_EXCESS` = 9, message "dead time {tau} ns of {ch}: the shortest
   analysed gap is {g0} ns, longer than its upper end by {x:.0f} mean edge spacings"; applies to bounds too
   (an upper limit below every analysed gap is contradicted). (c) Counts: `gap_check_floor_ps` per channel
   (None when the floor is <= 0: downward not checked). L1 in the docstring. Tests: `r4-dt/gap_miss.py` cases
   (46 vs 50.7 refused; 22 on 50 ns refused); a matched detector at 1 and 5 Mcps not refused.
7. **CLI and docs consistency (S5, S8).** (a) Delete `cli._grid_capture`; call `float(forward.grid_capture(...))`;
   test (500, 20, 10000) prints capture 0.95. (b) HOM forward path: `subtracted` from `corrections`.
   (c) cli.py module docstring and `--calibrations` help name CHSH. (d) schema `valid_for` description:
   "Enforced when the calibration is used (from acquisition, or creation, to `to`)"; `_check` reads
   `valid_for.get("from")` falling back to acquired/created; test with `{"to": ...}`. (e) CHSH quickstart dark
   note: "not attempted in the quickstart". (f) README: a short table of counts keys (`capture_offsets_ps`,
   `capture_sd_data`, `delay_drift_diff`, `delay_drift_z`, `min_gap_ps`, `gap_check_floor_ps`, `dataset`,
   `forward_check_p`).
8. **CI (C1-C6, C8).** (a) Job `gates`: `needs:` every coverage, coverage-merge, propagation,
   propagation-merge job (and python-minimum, bindings, core-wheels if intended as gates), `if: always()`,
   one step failing unless every `needs.*.result == 'success'`. (b) chsh-spad into `propagation-shard` with
   4 shards and its merge. (c) `propagation.run`: `chunk = max(1, min(CHUNK, -(-len(span) // jobs)))`.
   (d) `coverage.run`: one stderr line per finished chunk (as propagation). (e) ci.yml comment and the
   926560b statement: the measured per-replicate CPU after Fix 1. (f) `test_ci_coverage_matrix_runs_every_regime`
   also checks `propagation.REGIMES` and that each sharded regime has a merge job in `needs:`.
   (g) `concurrency: {group: ci-${{ github.ref }}, cancel-in-progress: ${{ github.event_name == 'pull_request' }}}`.
9. **Docs.** CONTRIBUTING "Gates" paragraph (C7: both commands, `--jobs`, shard/merge, CI seeds 20261111 and
   20261118, when a change needs them, method-version bump and `test_method_versions.py`, that `gates`
   blocks merges once required). CHANGELOG round-4 entries with L3. L1, L2 in README. docs/phase3-plan.md
   smaller-items rows: "User declaration of scan runs past the count-rate maximum | M4 council round 4 (DT3)
   | M10, with the past-maximum regime" and "Per-run branch flip search to restore the passive-quench gate's
   speed | M4 council round 4 | M10, checked on the past-maximum regime".

## Owner items

- Add the new `gates` job (Fix 8a) as a required status check in the `main` ruleset (27075146). Until then a
  red coverage or propagation gate does not block a merge.
- Agree the past-maximum coverage regime for M10 (already deferred in round 3), now also covering the branch
  margin and the deferred declaration and flip search.
