> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Merge readiness: PR #26 (m4/series-deadtime, HEAD e8c1aba)

## Verdict

**NO-GO at e8c1aba.** One test fails on Linux. Everything else is ready. Expected to be GO after one small test fix, a green re-run, and a new PR body.

CI run 44800651345 (CI on e8c1aba, 19:16-20:07 UTC) concluded **failure**. Three jobs failed:
- "Python 3.12, slow tests (ubuntu-latest)"
- "Build sdist and wheel, test the installed package" (a **required** check, so GitHub itself blocks the merge)
- "Statistical gates" (it fails because python-slow failed)

The cause is one test, `test_the_performance_changes_leave_the_gate_replicates_unchanged`, added in round 9 (python/tests/test_propagation.py:225-248). It fails in all 5 regimes on ubuntu. For example, hom-snspd gives `3367e4c94b6f0ae3` against the pinned `8d494923702d38aa`. The two Linux jobs produce identical hashes, and locally on macOS the test passes (`pytest -k "performance_changes and (hom-snspd or hom-bright)" -m slow`: 2 passed).

So the pin is a SHA-256 of full-precision `repr` floats taken on one platform. It is bit-exact on macOS and not portable to Linux libm/BLAS. It is not evidence of a value change:
- Every coverage and propagation gate job and both sharded merges passed on ubuntu.
- The deadtime figures equal round 8/9 exactly (series-spad tau_p 95.5%, tau_np 99.9%; passive-quench 94.6%, 99.0%).
- The propagation regimes show 0 refused, 0.0% forward misfit, drift and off-peak accepted (hom-crosstalk misfit artifacts 0.3%).
- The PPL benchmark run passed.

The round-9 claim "three gate replicates per regime hash as before" was therefore verified on macOS only.

Fix (must, before merge), either:
(a) hash the values rounded to about 10 significant digits (or compare the numbers with rtol 1e-9 against pinned values), or
(b) pin per platform.
Re-run CI, and merge only on a fully green run, "Statistical gates" included.

## Findings

| # | Severity | Finding | Evidence | Where | Suggested fix |
|---|---|---|---|---|---|
| 1 | blocker | CI red on HEAD: a cross-platform hash pin (see Verdict). It also fails a required check | run 44800651345; job logs saved in r10-merge/slow.log, pkg.log | python/tests/test_propagation.py:221-248 | Make the pin tolerant of last-digit float differences; re-run |
| 2 | blocker (process) | "Statistical gates" and the slow-test job are not required checks. Ruleset "main" (27075146) requires 11 contexts (Rust x3, MSRV, Python 3.10/3.12 x3, package). GitHub would let this PR merge with the gates red or still pending | `gh api repos/<owner>/dadl/rulesets/27075146` | ruleset | Owner item O1. Make it required, or check the gates by hand before merging |
| 3 | must decide | The merge method in O14 is constrained. The ruleset has `required_linear_history`, so "Create a merge commit" is refused. The plan's "non-squash merge" recommendation means **rebase-merge** here: 13 commits, branch 0 behind main. The owner role has `bypass_mode: always`, so a bypass could also push a merge commit | ruleset rules: `required_linear_history`, `update`, bypass RepositoryRole 5 | docs/phase3-plan.md:487 (O14) | Choose squash (as #20-#25) or rebase. Reword O14 as "rebase-merge (linear history is required)" |
| 4 | should, before merge | The council rulings exist only in a temporary scratch directory. The plan's Source column and the CHANGELOG cite "r8 Pl6", "r9 Rd3", "council round 6 addendum" (53 citations in the plan) with no record in the repo. Phase 2 published its reviews (`docs/reviews/2026-11-council.md`, `2026-11-review-cycles.md`). The scratch directory is lost on reboot | `ls docs/reviews`: no phase 3 council file; rulings in `.../scratchpad/council/r1..r9-council.md` | docs/reviews/ | Add `docs/reviews/2026-11-phase3-council.md` (rounds 1-9: adopted, limited, deferred, rejected, owner items), or at least copy the rulings somewhere durable now |
| 5 | should | The PR title and body are stale. They describe the M4 commit only (method 1, local coverage 99.2/97.3, tau_np "at least 92.4%", passive-quench in 4 shards). They miss 9 rounds that touch M2/M3 analyses, T3, the delay method (8 to 10), new refusals and the overrides | `gh pr view 26` | PR | Replace with the draft below |
| 6 | minor | The CHANGELOG is accurate but organised by review round: about 350 lines under "M4". The user-noticeable behaviour changes are scattered as inline "Behaviour change" tags, with no single list. O7 asks the owner to accept them "in release notes" | CHANGELOG.md:104-460 | CHANGELOG.md | Add a short "Behaviour changes since 0.1.1" list at the top of Unreleased (the eight items in the PR draft). This can wait for M10, but is cheap now |
| 7 | minor | Artifact compatibility is not stated for method bumps. delay goes from 8 (0.1.1) to 10. No minimum-version check runs at use (deferred to M5a, plan:442), so 0.1.x delay artifacts load silently. The impact is small: T3 is new in Unreleased, and Fix 4 keeps T2 outputs identical to the last digit. But the round 8 predictive check now uses (groups-1) degrees of freedom, so an old artifact's misfit flag was computed differently. deadtime-series 1 to 7 never shipped | `grep METHOD_VERSION`; `git show v0.1.1:python/tests/method_versions.json` (delay 8, no deadtime) | CHANGELOG.md:8-13 | One sentence: "delay artifacts from 0.1.x (method 8) are still accepted; refit them if their misfit flag matters; no method version is checked at use until M5a" |
| 8 | minor | "this branch starts it at acquisition" will read wrongly once merged | README.md:364; CHANGELOG.md:10-11 (both new in this PR) | README, CHANGELOG | "unreleased versions start it at acquisition" |
| 9 | minor | CONTRIBUTING says a change to a pinned module "also needs a new METHOD_VERSION". It omits the byte-identical re-pin exception that round 9 used for deadtime 7 and that test_method_versions.py:7 implements (`notes`) | CONTRIBUTING.md (Gates section); python/tests/method_versions.json deadtime `notes` | CONTRIBUTING.md | Add "unless outputs are byte-identical: re-pin the hash and add `notes`" |
| 10 | minor | The PyPI README (python/README.md) says that on T2/plain tags DADL "detects" pulsed light through the off-peak windows and refuses. It omits round 9's L4 power limit: reliable only at about half or more of the off-peak counts; smaller shares pass unlabelled, with V and \|E\| low. It defers to the root README's known limits | python/README.md:24-29 vs README.md:104-110 | python/README.md | Add "when pulsed accidentals are a large share (see known limits)" |
| 11 | minor | The plan's owner-items header says "Consolidated from council rounds 4-8", but the table includes round 9 (O15-O18). O7's status "Partly in CHANGELOG" is outdated: every item is now tagged in the CHANGELOG, though acceptance is still open | docs/phase3-plan.md:469, 480 | plan | "rounds 4-9"; status "In CHANGELOG; acceptance open" |

## Checked and found sound

- Working tree clean. Branch even with origin and 0 commits behind main. Mergeable.
- No stray files. No scratch output, logs or large files (none over 100 KB among added or modified files). The CmdStan binaries `bench/ppl/models/{dark,delay}` are deleted and gitignored (O9 history purge left as the default). No scratch or tmp paths leak into the repo.
- `python/tests/method_versions.json` matches the code: accidentals 3, dark 7, dcr 2, heralding 6 (sha unchanged from v0.1.1); delay 10; deadtime 7 with the re-pin note. `test_method_versions` passes locally.
- The CHANGELOG records every round 1-9 behaviour change listed in O7 and the r9 ruling (off-peak refusal in CHSH, `dead_time_contradicted`, `spdc-cw` refused on T3, the gap-test tolerance, the 20% dead fraction, `t3-range`, the refusal budget) and every method version step. The memory limit (L3) appears in the CHANGELOG and README:330.
- README: the override table includes `dead_time_contradicted` (README:480-483). The counts keys `offpeak_draws` and `dead_time_gap_problems` are documented. The refusal budget section, T3 `cw` text, L2/L4 and "unreleased" markers are in place. The CLI matches: `dadl calibrate deadtime --help`, `--allow` lists the name (guardrails.py:107). Version 0.1.1 is stated as intentional (O3).
- docs/phase3-plan.md: the M2b/M3b/M4b/M5a/M5b milestones, dependencies and smaller-items routing are consistent. O1-O18 are present, each with source, status and recommendation, and each is actionable (O1 names ruleset 27075146).
- CI structure: the `gates` job uses `if: always()` and fails on any non-success among the python-slow, coverage, shard and merge jobs, so a failed shard cannot pass as a skipped merge.

## Before merge vs after

Before merge:
1. Fix the hash-pin test (finding 1) and get a fully green run, "Statistical gates" included. Check the gates by hand, since finding 2 means GitHub won't require them.
2. Replace the PR title and body (draft below).
3. Owner picks squash or rebase (finding 3 / O14).
4. Preserve the council rulings (finding 4). Strongly advised: the scratch directory is temporary.
5. Optional and cheap: findings 6-11 as one docs commit. A docs-only commit still re-runs the roughly 80-minute CI.

After merge (owner items, not blockers): O1 required checks (do this first; it protects the next PR); O3/O4 dev version 0.2.0.dev0 and schema `$id` before more artifacts are written from main; O2 release gating; O5; O6/O13 milestone scope; O7 acceptance at release; O8 pre-registration sitting before any M4b/M5a run; O10 a real divided T3 file; O12 refusal tolerance; O15-O18.

## Replacement PR title and body (draft)

Title: M4 series dead time, plus nine review rounds across the corrected analyses

Body:

M4 of docs/phase3-plan.md, plus the fixes from nine review-council rounds. The rounds touched the M2/M3 corrected analyses, T3 input and the delay calibration, not only M4. Each round is one commit, and the CHANGELOG (Unreleased, "Phase 3, M4") records them round by round.

What M4 adds
- Simulator: a paralyzable stage (`dead_time_p_ps`, default 0) before the non-paralyzable dead time, in every event path. The shared Rust/Python vectors gain series cases.
- Calibration: kind `deadtime-series` (D4), method `series-gap-mle` 7. It fits tau_p and tau_np from illuminated runs at two or more rates, using the exact series-chain gap law and pending afterpulses. tau_p is an upper limit when the hard-edge model is not rejected.
- CLI: `dadl simulate SETUP --deadtime-scan` and `dadl calibrate deadtime`.
- Coverage regimes series-spad and passive-quench (1,000 replicates; passive-quench in 12 CI shards with a merge step).

What the review rounds changed beyond M4
- Corrected HOM, fringe and CHSH: an off-peak flatness check (sequential Monte Carlo), dead-time gap tests against the analysed data, the 20% dead-fraction refusal, delay-drift noise, and labels per value.
- Analysis-time refusals share one budget: at most 1e-3 per clean analysis.
- T3: the dtime lattice and the sync period in the delay fit and the predictive check. The sync source declaration is checked. T2 and T3 calibrations are compared by mode.
- Calibration sets: validity from acquisition, tagger and detector checks, and recorded conditions labelled when not compared.
- Performance: sorted input is not re-sorted, and CHSH grid counts are memoised. Values are unchanged.

Behaviour changes users will notice
- New refusals: identified dark dead time above 20% dead fraction; T3 without an explicit `source.type: cw` (the templates' `spdc-cw` no longer counts); T3 range and divided T3 in `calibrate dark`/`deadtime`; non-flat off-peak windows (pulsed light), now in CHSH too; analysed gaps that contradict the dark dead time.
- A new override, `dead_time_contradicted`. `calibration_outside_validity` no longer covers the gap tests.
- Refusal hints give the CLI form (`--allow NAME=REASON`).
- Method versions: delay `xcorr-gaussian-mle` 8 -> 10 (T3 lattice and sync edges; T2 outputs unchanged). deadtime `series-gap-mle` is new at 7. accidentals, dark, dcr and heralding are unchanged.
- Artifacts: still `calibration/1`, with the new kind `deadtime-series`. 0.1.x delay artifacts still load (no version check at use until M5a). `valid_for` now starts at acquisition.
- Known limit: analyses hold all tags in memory (about 80-150 GB for a 10-minute 1 Mcps CHSH run; chunked analysis in M9).
- The package still reports 0.1.1 (owner item O3).

Gates (CI seed, 1,000 replicates)
- Deadtime coverage: tau_p/tau_np 95.5/99.9% (series-spad), 94.6/99.0% (passive-quench).
- Delay coverage is identical to round 8.
- Propagation: every clean regime covers, with no refusals. hom-crosstalk is flagged in 100% of replicates. The "Statistical gates" job aggregates every coverage and propagation job.
- Unit, slow, minimum-dependency, bindings, wheels and package jobs: green on this commit (update after the re-run).

Owner items (docs/phase3-plan.md, "Owner items")
- Before or at merge: O14, the merge method. main requires linear history, so the choice is squash or rebase.
- After merge: O1 make "Statistical gates" a required check; O3/O4 dev version and schema `$id`; O2 release gating; O5; O6/O13 milestone scope; O7 accept the behaviour changes in the release notes; O8 pre-register gate changes; O10 a real divided T3 file; O12 refusal tolerance; O15 template source type; O16 chunked analyses in M9; O17 analytic gradient in M10; O18 CI budget.

Stays out: soft SNSPD recovery and passive-quench efficiency curves. The forward model still uses the dark non-paralyzable dead time until M5a.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
