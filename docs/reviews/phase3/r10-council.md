> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 10 (final: regressions, first-time user, merge readiness)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: **Rg** = r10-reg.md rows 1-5, **Ux** = r10-ux.md rows 1-13, **Mg** = r10-merge.md rows 1-11.
Council probes and outputs are in `scratchpad/council/r10-council/`:
- `capfix.py` (variants `head`, `A` = profile runs to m + tick, `B` = the council fix below), outputs `capfix_*.txt`
- `hashpin.py`, outputs `hash_*.txt` and `vals_*.json`: the round-9 replicate hash at full precision and rounded,
  under 1 and 8 BLAS threads, plus the flattened values (the form of the fixture below)

Round 8's refusal budget (1e-3 per clean analysis) and label rule (a label firing on more than 1e-3 of clean
analyses states its rate) stand. No ruling below adds a statistical refusal. This is the last review round: after
the fixes, the merge auditor's check is a fully green CI run plus the spot checks named in the order below.

## Checked numerically

All runs at HEAD e8c1aba, repo venv, one BLAS thread unless stated.

### Rg1 (confirmed): the tau_p cap ends below the truth on coarse ticks

`capfix.py`, tau_p = tau_np = 22000 + U(0, tick) ps, tags floored to the tick, 2/5/10 Mcps, 20 ms. Seeds 500+ are
the reviewer's; seeds 9000+ are fresh.

| tick (ps) | seeds | HEAD tau_p covered | fix B tau_p covered | tau_np covered (both) | median tau_p width HEAD / B |
|---|---|---|---|---|---|
| 25 | 500+, 20 | 2 / 20 | 19 / 20 | 20 / 20 | 1452 / 1477 |
| 25 | 9000+, 60 | **8 / 60** | **60 / 60** | 60 / 60 | 1011 / 1036 |
| 80 | 500+, 20 | 4 / 20 | 19 / 20 | 20 / 20 | 1139 / 1219 |
| 80 | 9000+, 60 | **7 / 60** | **60 / 60** | 60 / 60 | 921 / 1001 |
| 1 (control) | 9000+, 60 | 59 / 60 | 59 / 60 (same miss, s9) | 60 / 60 | 1076 / 1077 |
| 80, tau_p = 0.8 tau_np | 500+, 20 | 18 / 20 | 18 / 20 | 20 / 20 | (unchanged) |

Every HEAD miss has the upper end exactly at m (22000) with the truth up to one tick above. Variant A (profile
to m + tick) and B give the same intervals at tick 25. **B is adopted**: it does not evaluate the binned likelihood
above the shortest gap, and it leaves every fit whose profile stops below m unchanged. The width costs about one
tick, only when the profile reaches m. The invariants tau_p <= tau_np (value, upper end, draws) were asserted on
every run and hold.

### Mg1 (confirmed): the replicate hash is a single-platform pin

- On macOS (`hash_*.txt`), the full-precision hashes equal the pins (hom-snspd `8d494923702d38aa`, chsh-spad
  `5aa8e5cff09ddac9`) with 1 and with 8 BLAS threads. The CI log (`r10-merge/slow.log`:282-397) shows all five
  regimes differing on ubuntu, with both Linux jobs agreeing with each other. That is a platform difference in
  libm/BLAS last bits, not a value change: every gate job on the same run passed.
- **The merge auditor's option (a), hashing values rounded to 10 significant figures, is rejected.** A value whose
  11th digit straddles a rounding point still flips the hash. chsh-spad carries 159 floats in 3 replicates. With
  last-bit noise of 1e-14 relative, the chance that at least one rounds differently is about 159 x 2e-4, or 3%, per
  regime per run: a flaky test. Option (b) alone (per-platform hashes) still hides how large the difference is.
- **Council form:** compare the flattened values against a pinned JSON fixture, with integers, strings, booleans
  and None compared exactly, and floats with a relative tolerance. A refactor bug of the kind the test guards
  against (de-duplication, sort, memo key) changes counts (integers, exact) or values far above 1e-9.

### Rg2, Rg4, Rg5, Ux1, Ux2, Ux3 (confirmed by code read or rerun)

- Rg2: `forward.py:476-481` adds `DEAD_TIME_ABOVE_GAPS` with no rate text, independently of `DEAD_TIME_BELOW_DARK`
  (`:466`) on the same channel. With a bound (`used = upper`), "a non-paralyzable dead time cannot be" is said of an
  upper limit.
- Rg4: binomial chance that the chsh-dead-time control (true 87.5%) falls below the band low end: 0.37 at 100
  replicates (band 87%), 0.70 at 200 (89%), 0.98 at 400 (91%).
- Rg5(b): `Interval(16e6, 281121, 3.2e7)` prints `16000000 (95% interval [0, 32000000])` (`stats.py:36-50`).
- Ux1: `cli.py:79` prints only the fingerprint and "(modified)"; the override list is printed only in the
  quickstart (`cli.py:295`).
- Ux2: `guardrails.labels_for` and `propagate.py:481-482` add `POLICY_MODIFIED` whenever `policy.is_modified`,
  including on the raw value.
- Ux3: `forward.py:899` gives no remedy.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| Rg1 | tau_p cap ends at m while tau_np is widened a tick; coverage 7-8 / 60 on coarse ticks, fully paralyzable | Y/Y/Y | ADOPT-FIX (variant B) + ADOPT-LIMIT | Fix 2; L2 |
| Rg2 | `DEAD_TIME_ABOVE_GAPS` breaks the label rule; duplicates `DEAD_TIME_BELOW_DARK` | Y/Y/Y | ADOPT-FIX | Fix 3 |
| Rg3 | `dark_*` keys will be labelled "not compared" once M5a loads series artifacts | -/Y/Y | DEFER (M5a entry items) | plan |
| Rg4 | 100-replicate quick gate cannot judge the chsh-dead-time control rule | -/Y/Y | ADOPT-FIX (describe text) + ADOPT-LIMIT; `passed` unchanged (a gate rule, O12) | Fix 6; L1 |
| Rg5 | "calibration 4.4e-16"; identified bound printed [0, ...] | Y/Y/Y | ADOPT-FIX | Fix 3 |
| Ux1 | analyze never prints the overrides or reasons | Y/Y/Y | ADOPT-FIX | Fix 4 |
| Ux2 | raw value labelled "policy modified" for correction-only overrides | Y/Y/Y | ADOPT-FIX | Fix 4 |
| Ux3 | CHSH dead end on missing delay-jitter artifacts | Y/Y/Y | ADOPT-FIX | Fix 4 |
| Ux4 | policy refusals one per run | -/Y/Y | ADOPT-FIX (one pre-flight, same check function) | Fix 4 |
| Ux5 | condition refusal offers only the override | -/Y/Y | ADOPT-FIX | Fix 4 |
| Ux6 | flat label list, no severity order | Y/N/Y on three severity groups; Y/Y/Y on one label per line | ADOPT-FIX (one per line, overrides first, fixed excludes once) + DEFER (severity grouping, M10) | Fix 5 |
| Ux7 | no label when the predictive check did not run | Y/Y/Y | ADOPT-FIX | Fix 4 |
| Ux8 | README cites a constant name the output never shows | -/Y/Y | ADOPT-FIX (README quotes opening words) + DEFER (printed codes, M10) | Fix 5 |
| Ux9 | no single label and refusal table; `--allow` help lists one name | -/Y/Y | ADOPT-FIX (table, help list, sync test) | Fix 5 |
| Ux10 | "dead fraction of 2.66" | Y/Y/Y | ADOPT-FIX | Fix 4 |
| Ux11 | corrected "±" unlabelled; calibration part reads as total | Y/Y/Y | ADOPT-FIX | Fix 4 |
| Ux12 | HOM without dark artifacts: no "dead time not modelled" label | Y/Y/Y | ADOPT-FIX | Fix 4 |
| Ux13 | error above header in redirected logs | -/Y/Y | ADOPT-FIX | Fix 4 |
| Mg1 | replicate hash fails on Linux; required check red | Y/Y/Y | ADOPT-FIX (fixture with tolerance; not the rounded hash) | Fix 1 |
| Mg2 | "Statistical gates" not a required check | -/Y/Y | Owner (O1, moved before merge) | owner |
| Mg3 | linear history: squash or rebase | -/Y/Y | Owner (O14 reworded) | owner |
| Mg4 | council rulings only in a scratch directory | Y/Y/Y | ADOPT-FIX (record designed below) | Fix 7 |
| Mg5 | PR title and body stale | -/Y/Y | ADOPT-FIX (draft, updated for round 10) | Fix 9 |
| Mg6 | no consolidated behaviour-changes list | -/Y/Y | ADOPT-FIX | Fix 8 |
| Mg7 | artifact compatibility of method bumps unstated | Y/Y/Y | ADOPT-FIX + ADOPT-LIMIT | Fix 8; L3 |
| Mg8 | "this branch starts it" | -/Y/Y | ADOPT-FIX | Fix 8 |
| Mg9 | CONTRIBUTING omits the byte-identical re-pin | -/Y/- | ADOPT-FIX | Fix 8 |
| Mg10 | python/README omits L4 | Y/Y/Y | ADOPT-FIX | Fix 8; L4 |
| Mg11 | plan header "rounds 4-8"; O7 status stale | -/Y/- | ADOPT-FIX | Fix 8 |

**Votes on Ux6.** P and U: a label that bounds the number must not sit beside bookkeeping. M votes N on the
three severity groups now: that needs a severity on every label constant (new API surface the gates and tests do
not cover) the day before merge. All three accept one label per line with the user's overrides printed first and
the constant excludes list printed once. The grouping goes to M10 with the printed label codes (Ux8).

**Votes on Rg3.** M and U: unreachable today (`SERIES_UNUSED`), and M5a rewrites that path. P abstains.

**Seats on items marked "-".** P abstains on CLI-text, CI and process items.

## Limit wordings

- **L1** (CONTRIBUTING, Gates; `propagation.describe` for control regimes below 400 replicates): "Quick gate runs
  (100-200 replicates) judge corrected coverage only. The control rule of chsh-dead-time needs about 400
  replicates: its control covers about 87.5%, and at 100 replicates it lands inside the band about one run in
  three. A FAIL on the control rule in a quick run is not evidence of a problem; re-run at 400 or more."
- **L2** (README dead-time limits, after round 1's L3; `calibrate_deadtime_series` docstring): "On a coarse tag
  tick (PTU T3, about 25-80 ps), a dead time close to fully paralyzable has tau_p's upper end at tau_np's upper end,
  one tick above the shortest recorded gap, since either stage can lie up to a tick above it. The interval is then
  about one tick wider than the profile alone (council round 10, Rg1: coverage 60 / 60 against 7-8 / 60 before)."
- **L3** (CHANGELOG, top of Unreleased; README calibration artifacts): "delay artifacts from 0.1.x (method 8)
  still load. Their misfit flag was computed before round 8's (groups - 1) degrees of freedom; refit them if it
  matters. No method version is checked at use until M5a."
- **L4** (python/README.md, the off-peak sentence): "... and refuses when pulsed accidentals are a large share of
  the off-peak counts (about half or more in a typical CHSH run; smaller shares pass unlabelled, with V and |E|
  reading low; see Known limits)."

## Implementation order

1. **Fix 1 (Mg1; tests/test_propagation.py, new `tests/data/r8_gate_replicates.json`). Blocker.**
   - On macOS at HEAD, first confirm the five full-precision hashes equal `R8_REPLICATE_HASHES`, then dump the
     flattened `(raw, corrected, counts without offpeak_p/offpeak_draws, misfit)` values, three replicates per
     regime, as `hashpin.py` does: a list of [path, value]. Record the generating command in the test docstring.
   - The test compares paths exactly, non-floats exactly, and floats with
     `math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)`, with NaN equal to NaN. On failure it lists the five worst
     paths with their relative differences. Remove the hash dict.
   - Rg1 does not touch this: the gate regimes never load deadtime-series artifacts.
   - If the Linux run deviates above 1e-9: when the worst is below 1e-8, set `rel_tol` to 1e-7 with a comment
     giving the observed figure. Above that, stop: that is an optimiser-path difference, and O20 applies.
   - Test of the test: a fixture with one count changed by 1, and one float changed by 1e-6 relative, fails.
2. **Fix 2 (Rg1; calibration/deadtime.py, README, CHANGELOG).** Before `cap_high = min(est_p.high, est_np.high)`:

   ```python
   if high >= tau_hi:  # the profile reached the shortest gap: tau_p, like tau_np, may sit a tick above it (r10 Rg1)
       est_p = Interval(est_p.value, est_p.low, max(est_p.high, est_np.high))
       stat_high = max(stat_high, est_np.high)
   ```

   - The draws need no change: they are already capped by the widened `np_draws`.
   - **deadtime `series-gap-mle` 7 -> 8**; regenerate `tests/method_versions.json`, drop the round-9 `notes`.
   - Tests:
     - fast: one fresh seed of the `capfix.py` case (tick 25, tau_p = tau_np off grid) covers;
     - slow: 20 seeds at tick 25 and at tick 80 cover in at least 18 / 20 each (council 60 / 60);
     - control: at tick 1, and at tick 80 with tau_p = 0.8 tau_np, the interval equals method 7's;
     - tau_p <= tau_np (value, upper end, every draw) in all of them.
   - Mutation harness: add the mutant "drop the new block"; the slow test must kill it.
3. **Fix 3 (Rg2, Rg5; forward.py, stats.py, result.py).**
   - Rg2:
     - Skip `DEAD_TIME_ABOVE_GAPS` for a channel that already has `DEAD_TIME_BELOW_DARK` or a gap problem.
     - New text: "the dark dead time of {name} ({x}) is longer than {k} of the analysed gaps (shortest {g}); this
       fires on most analyses whose dark run is short against the analysed one (council round 8, 400 of 400
       chsh-dead-time replicates). The correction uses the dark value and its uncertainty (see Known limits)".
     - With a bound: "the dark run bounds the dead time of {name} below {x}; {k} analysed gaps are shorter, so the
       bound is loose; the correction uses the bound (see Known limits)".
     - Tests: r6 SNSPD seed 58 never carries both labels on one channel; the existing round-8 test matches the new
       opening words.
   - Rg5(a): print the calibration part as `0` when it is below 1e-12 of |value| (or exactly 0).
   - Rg5(b): `Interval.__str__` never rounds a positive end to 0. Raise `d` until every nonzero end shows two
     significant figures. Test: `Interval(16e6, 281121, 3.2e7)` shows a nonzero lower end; `Interval(1, 0.98, 1.02)`
     is unchanged.
4. **Fix 4 (Ux1-5, Ux7, Ux10-13; cli.py, guardrails.py, propagate.py, forward.py).**
   - Ux1: after the policy line, print `Overrides:` with one line per allowed name: NAME, reason, and for
     `calibration_outside_validity` / `dead_time_contradicted` the artifact ID and the condition or problem
     overridden (from the result's labels and counts). Reuse the quickstart's line. README:475 "recorded" becomes
     "printed, and kept in the result's `overrides`".
   - Ux2: when every override only adds allowed corrections, the raw value's label is "policy loosened for the
     corrected value only; raw value unchanged". `POLICY_MODIFIED` stays on the corrected value. When no corrected
     value is produced, print "--allow had no effect: no corrected value". Test: raw S under `--allow dead_time`
     with no calibrations carries no "policy modified".
   - Ux3: append to the `forward.py:899` message: "take a heralding run for each missing pair, or move the
     delay-jitter artifacts (.json and .posterior.npy) out of the folder to correct without window capture". In
     the CHSH quickstart epilogue, write the pairs artifact to its own folder and say why.
   - Ux4: before the first refusal, `cmd_analyze` asks the same `Policy.check_correction` about every correction
     the calibration set implies, and refuses once listing every missing `--allow NAME=REASON`. Test: the user
     advocate's HOM folder gives one refusal naming dead_time and window_capture.
   - Ux5: "not given for the dataset: pass --condition NAME=<value at acquisition>; if it is unknown, --allow
     calibration_outside_validity=REASON".
   - Ux7: corrected label "no predictive check (needs delay-jitter artifacts): crosstalk and peak-shape misfits are
     not tested", for CHSH, HOM and fringe without delay-jitter artifacts.
   - Ux10: "rate x dead time = 2.66 (limit 0.2)".
   - Ux11: "± 0.0105 (1 sd)"; "calibration part 6.3e-06 (statistical part of the artifacts only)".
   - Ux12: the CHSH "dead time (no dark artifacts) not modelled" label on HOM and fringe.
   - Ux13: flush stdout before writing the error.
5. **Fix 5 (Ux6, Ux8, Ux9; cli.py, README).**
   - Labels printed one per line, indented, with the overrides first. The fixed "Calibration part excludes" list
     printed once, after them.
   - README table "Labels and refusals": the printed opening words, refusal or label, raw or corrected, what to do,
     override name. Include the flag correspondence `--accept-misfit` / `misfit_calibration` and
     `--allow-condition-mismatch` / `calibration_outside_validity`.
   - `--allow` help lists all `KNOWN_CORRECTIONS` (9). README:196, :514 quote the label's opening words instead of
     the constant name.
   - Test: every name in `KNOWN_CORRECTIONS` appears in the README table and in `dadl analyze --help`.
6. **Fix 6 (Rg4; analysis/propagation.py, CONTRIBUTING).** In `describe`, for a control regime below 400
   replicates: "control rule not decidable below about 400 replicates". L1 in CONTRIBUTING. `passed` unchanged.
7. **Fix 7 (Mg4; docs/reviews/). Form of the record:**
   - `docs/reviews/2026-11-phase3-council.md`, the summary, in the form of `2026-11-review-cycles.md`:
     - opening paragraph: AI reviewers, seats, the two-of-three rule, the refusal budget;
     - per round: lenses, then a table: id, severity, finding (one line), ruling, where fixed (commit and method
       version); then the round's limits by id, with the file where each wording now lives;
     - closing: the owner items table, as in the plan.
   - `docs/reviews/phase3/`: verbatim `rN-council.md` (rounds 1-10) and the reviewer reports `rN-<lens>.md` that
     define the ids. Probe directories are not copied; a header line in each file says so.
   - Replace every absolute scratch or home path with `<scratch>/` (r1-gate.md contains a home path).
   - Id scheme: an id is resolved as "r{N} {Id}" or "council round N, {Id}". Within each round the ids are unique.
   - Test (`tests/test_review_record.py`):
     - every `r\d+ [A-Z][A-Za-z]*\d+` and `council round \d+[^)]*?, [A-Z][A-Za-z]*\d+` cited in the plan,
       CHANGELOG, README and `python/dadl` appears as a row id in that round's section of the summary;
     - no file under `docs/reviews/phase3/` contains a home or scratch path.
8. **Fix 8 (Mg6-Mg11, round 10 CHANGELOG; docs).**
   - CHANGELOG: "Behaviour changes since 0.1.1" at the top of Unreleased (the PR draft's eight items, plus round
     10's: deadtime 8, overrides printed, all needed `--allow` names in one refusal, raw values no longer marked
     "policy modified"), L3, and a round 10 section.
   - README:364 and CHANGELOG:10-11: "unreleased versions start it at acquisition".
   - CONTRIBUTING: "unless outputs are byte-identical: re-pin the hash and add `notes`".
   - python/README: L4.
   - Plan:
     - header "rounds 4-10";
     - O7 status "In CHANGELOG; acceptance open";
     - Rg3 in the M5a entry items;
     - Ux6 grouping and Ux8 codes in M10;
     - the owner table below.
9. **Gates and merge check.**
   - Re-run `coverage deadtime`, 1,000 replicates, CI seed (series-spad, passive-quench). tau_p coverage may only
     rise against 95.5% / 94.6%, and tau_np must be equal (99.9% / 99.0%). Record both in CHANGELOG.
   - The propagation gate and `coverage delay` run in CI and must equal the record. No value change is expected:
     Rg1 is outside analyses; the rest is labels and display.
   - Re-run `test_published_results.py`; report any changed label text.
   - Re-run the r7 `mutate.py` harness with the Rg1 mutant.
   - The full CI must be green, "Statistical gates" included. Then the merge auditor checks the Linux tolerance
     outcome (Fix 1), confirms `method_versions.json` shows deadtime 8, and replaces the PR title and body with
     the r10-merge draft. In the draft: deadtime `series-gap-mle` is new at 8, the gate figures come from the
     re-run, and the round count is ten.

**Method versions.**
- deadtime `series-gap-mle` 7 -> 8.
- delay stays 10.
- accidentals 3, dark 7, dcr 2 and heralding 6 are untouched, with sha256 unchanged.

## Final owner items (consolidated O1-O18, plus round 10)

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
