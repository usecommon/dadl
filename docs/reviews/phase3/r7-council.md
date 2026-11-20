> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 7 (time tagger, security, test quality)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: G = time-tagger engineer (r7-tag.md rows 1-6 as G1-G6), S = security (r7-sec.md rows 1-6 as
S1-S6), Q = test quality (r7-test.md rows 1-9 as Q1-Q9). Council probes and outputs:
`scratchpad/council/r7-council/` (copies of the reviewers' `t3tick.py`, `realhead.py`, `pulsed.py`; outputs
`t3tick_25_12500.3.txt`, `t3tick_25_12500.txt`, `realhead.txt`, `pulsed.txt`).

## Checked numerically

All runs on HEAD c7d3cf2 with the repo venv, through the package's own functions, one BLAS thread.

### G1: the T3 tick is lost (confirmed, wrong)

`t3tick.py` (hard edge, true tau_np = 22000 + U(0, res), 2 and 5 Mcps, 0.01 s, tags through `_t3_to_ps`),
20 seeds each:

| res, period | tick seen by `tag_tick` | `tag_tick_ps` recorded | tau_np covered | covered with +-res |
|---|---|---|---|---|
| 25 ps, 12500.3 ps | 1 | 1.0 | **10/20** | 20/20 |
| 25 ps, 12500 ps | 25 | 25.0 | 20/20 | 20/20 |

The reviewer's 29/60 and 59/60 reproduce in kind. The intervals are 14 ps wide where the recorded gaps are
quantized at 25 ps (seed 2: [22012, 22026] for a truth of 22006.5). Code read: `deadtime.py:642` takes
`tag_tick(*tags)` (gcd of offsets), which a non-whole or non-multiple period reduces to 1; `tagfile.tag_tick_ps`
(tagfile.py:287-305) returns None or gcd(res, period) in T3. Physics: each T3 time is floored to the dtime bin
within its sync, so any recorded gap is within one dtime resolution of the true gap, whatever the period; the
round-6 law holds with tick = res (plus 1 ps when the reader rounded tags to whole ps).

Related, not reported: `delay.py:220` and `forward_check.py:237,304` also take `tag_tick` from the tags, so on
such T3 data they bin at 1 ps multiples on a 25 ps comb. `calibrate_delay` already accepts `tick_ps`; the fix
passes it (below). Whether the forward check then misfires is unmeasured; Fix 1 adds a test.

### G2: the T3 sync period checked as a setting at rtol 1e-9 (confirmed)

`realhead.py`: two files of one folder at 79,999,870 and 79,999,860 Hz give ConfigError "the tagger settings
differ between files (tagger_MeasDesc_GlobalResolution)"; across folders 10 Hz apart,
`condition_mismatches` returns `tagger_MeasDesc_GlobalResolution = 1.25e-08 recorded, 1.25e-08 for the
dataset` (`propagate.py:128-131`, `EQUAL_RTOL = 1e-9` at :64, `%g` format). The council did not settle whether
PicoQuant writes the period from the measured sync rate; the ruling does not depend on it. In T3 the value is
the source's repetition period (a property of the laser, which drifts at the ppm level), not a tagger setting
that places tags, and a 1e-9 equality on it is wrong under either convention. A changed period moves
cross-sync times by f x gap only (round-7 G "sound" bullet 3).

### G3: pulsed T3 HOM with CW accidentals (confirmed, wrong)

`pulsed.py` (rng 2, 2e7 pulses, W = 1000 ps): at 12500 ps the side estimate is 1488 against c_min 1484 and the
corrected value is not computed; at 13157.9 ps the side estimate is **372.4** against c_min **1858** and the
corrected V is **0.9622 +- 0.0011** with only the standard two-point label (true accidental-free V = 1).
`side_offsets` (calibration/accidentals.py:65-78) already takes `rep_period_ps`, but nothing passes it, and
an integer step drifts by k x (period fraction) over the 10 us guard. `config.py:228` refuses `source.type`
other than CW in setup files only; the tags path has no such tie.

### S1: sharded gate merges fail open (confirmed by code read)

`coverage.merge` (coverage.py:595-623): an empty `results` list yields `out = []`; the gate is taken from the
shard record (`first["gate"]`, :618, :621). `main` then exits 0 with nothing printed. The reviewer's
`empty.json` and `lax.json` outputs agree with this reading.

### Q2-Q4, Q8: spot checks

- Q3(a): `chsh.py:459` does include `gap_labels`; the finding is that removing it is not caught (mutant
  `h_drop_gap_labels`), not a code bug.
- Q4: `fl_chsh_cov.txt` 5/20, `fl_hom80.txt` 69/80, `fl_dt_recovery.txt` 15/20, `fl_bounded_skip.txt` 14/20
  as reported. `test_chsh_corrected.py:65` asserts a 1.96-se exclusion whose expected excess is about 1.6 se.
- Mutation totals in `r7-test/results.txt`: 84 killed (5 by the source hash only), 48 survived, 53
  behavioural survivors; consistent.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| G1 | T3 tick lost; round-6 tau_np widening is 1 ps (10/20 cover at res 25) | Y/Y/Y | ADOPT-FIX | Fix 1 |
| G2 | T3 `MeasDesc_GlobalResolution` checked at rtol 1e-9; `%g` hides the difference | Y/Y/Y | ADOPT-FIX | Fix 2 |
| G3 | Pulsed T3 HOM: side-window accidentals low by about W/period; corrected V 0.962 vs 1, no label | Y/Y/Y | ADOPT-FIX + ADOPT-LIMIT + DEFER (pulsed side peaks, v0.2+ as plan:410) | Fix 3; wording L1 |
| G4 | No check that 2^bits x res covers the T3 sync period | Y/Y/Y | ADOPT-FIX | Fix 4 |
| G5 | T3 sync divider convention unverified | Y/Y/Y (guard); Y/Y/N (verification deferred) | ADOPT-FIX (guard) + DEFER (verification, M9) | Fix 4; M9 entry item: one real divided T3 file |
| G6 | Tagger model and serial not checked at use time | Y/Y/Y | ADOPT-FIX | Fix 5 |
| S1 | Gate merge passes on empty results; band read from the shard file | Y/Y/Y | ADOPT-FIX | Fix 6 |
| S2 | Crafted `.npy` header: MemoryError escapes `ArtifactError` | N/Y/Y | ADOPT-FIX | Fix 7 |
| S3 | PTU memory figure holds for T2 only | Y/Y/Y | ADOPT-FIX (docs) | Fix 8, wording L2 |
| S4 | `ppl-bench.yml` checkouts persist the token | -/Y/Y | ADOPT-FIX | Fix 8 |
| S5 | Committed Mach-O binaries `bench/ppl/models/{dark,delay}` | -/Y/Y | ADOPT-FIX (remove in a normal commit) + OWNER (history) | Fix 8; owner item 2 |
| S6 | A tag on a main commit with no completed gates can release | Y/Y/Y (risk real) | OWNER | Owner item 1 (release policy) |
| Q1 | Gate pass/fail logic and coverage-merge tiling untested (10 survivors) | Y/Y/Y | ADOPT-FIX (tests) | Fix 9 |
| Q2 | Channel order and declared-delay sign untested (V 0.910 -> 0.835 unnoticed) | Y/Y/Y | ADOPT-FIX (tests) | Fix 10 |
| Q3 | Wired labels and refusals never observed through the analysis path (10 survivors) | Y/Y/Y | ADOPT-FIX (tests) | Fix 11 |
| Q4 | Two tests pass only at their seed (5/20, 69/80) | Y/Y/Y | ADOPT-FIX (tests) | Fix 12 |
| Q5 | `test_recovery_at_three_points`: three 95% intervals at one seed (15/20); 5 ns point ungated | Y/Y/Y (limit); Y/N/Y (regime at M10) | ADOPT-LIMIT + DEFER (M10) | Wording L3; M10 item: a 5 ns / 22 ns series regime |
| Q6 | Deadtime interval logic defended by a golden test whose oracle is the code | Y/Y/Y | ADOPT-FIX (tests) | Fix 13 |
| Q7 | Forward-check statistics unprotected below the gate | Y/Y/N | ADOPT-FIX (tests) | Fix 14 |
| Q8 | `pytest.skip` on a seed premise; `_dark_artifact_for` channel/multiple untested | Y/Y/Y | ADOPT-FIX (tests) | Fix 15 |
| Q9 | Slow tests in the default tier; `slow` marker never deselected | N/Y/Y | ADOPT-FIX (CI) | Fix 16 |

Votes on G1. P: a 14 ps interval on 25 ps-quantized gaps is the round-6 defect again, on exactly the hardware
(T3 at a laser's measured period) where it is the normal case. M: the fix is a tick passed in, not a new law;
the simulator's 1 ps grid path is unchanged, so the coverage gate figures should not move. U: yes.

Votes on G3. P: a corrected V 35 sd below the truth (0.962 against 1, sd 0.0011) printed with no label is wrong;
pulsed accidentals sit at k x period and DADL does not model them. M: refuse rather than half-model (an integer
`rep_period_ps` drifts over the guard; the pulsed route is v0.2+). U: refuse with a message that says what to
declare; never silent.

Votes on S2. P: not a physics issue, no vote of substance (N: low priority). M: the `ArtifactError` contract and
the "skip unreadable" path in `calibrate deadtime --dark` are both broken by it; the fix is a few lines. U: a
traceback from a downloaded artifact folder is a bad experience.

Votes on Q7. U: no user-visible effect today (N). P, M: alpha and Bonferroni define the misfit rate the gate
relies on; two cheap tests.

Votes on Q9. P: test tiering is not a physics question (N). M, U: 7.5 min serial per matrix job; the marker
already exists.

## Limit wordings

- **L1** (README line 18 after "PicoQuant PTU (T2 and T3) input", and the `adapters/ptu.py` module docstring T3
  paragraph): "T3 files come from a sync (usually pulsed) source. DADL models CW light: in a T3 file the sync is
  used as a clock only. Corrected values are refused on T3 input unless the setup declares `source.type: cw`;
  pulsed accidentals sit in side peaks at multiples of the repetition period, which DADL does not yet model."
- **L2** (adapters/ptu.py:54-56): "The whole record body is decoded in memory: about 100 bytes per record at
  peak in T2 (a 160 MB file took 4.2 GB; with markers about 215 bytes per record) and up to about 170 bytes per
  record in T3, which adds a sync tag per occupied sync period."
- **L3** (docstring of `test_recovery_at_three_points`, `test_tau_p_never_exceeds_tau_np`,
  `test_afterpulse_is_held_at_the_dark_artifact`): "Smoke test at one seed: it asserts 95% intervals once, so it
  is not a coverage check (about 86% of seeds pass all three points). Coverage is gated by
  `coverage deadtime` in the series-spad and passive-quench regimes; the 5 ns / 22 ns point is in neither (M10)."

## Implementation order (adopted fixes)

1. **Fix 1 (G1; data.py, tagfile.py, adapters/ptu.py, calibration/deadtime.py, cli.py).** (a) `PtuFile`
   metadata gets `tick_ps`: T2 the whole-ps `MeasDesc_GlobalResolution`; T3 `ceil(MeasDesc_Resolution in ps)`,
   plus 1 when `exact_ps` is False (each tag rounded by up to 0.5 ps). (b) `tag_tick_ps(path)` returns that
   value for T3 (the dtime resolution), not gcd(res, period) or None. (c) `calibrate_deadtime_series(...,
   tick_ps=None)`: `tick = max(tick_ps, tag_tick(*tags))` when given; record `tag_tick_ps` as used and
   `tag_tick_source` ("declared" or "tags"). The CLI `calibrate deadtime` passes `tag_tick_ps(folder)`.
   (d) The CLI `calibrate pairs` passes the same value to `calibrate_delay(tick_ps=...)` (delay.py unchanged;
   its existing note fires if tags are off the grid). Tests: `t3tick.py`'s generator as a fixture,
   res 25 ps, period 12500.3 ps, 20 seeds (slow-marked; one seed fast): tau_np covers in >= 19/20 with the
   tick passed (council: 20/20) and `tag_tick_ps` recorded 25 (26 with rounding); `tag_tick_ps` on a written
   T3 folder with a non-whole period returns the dtime resolution; `forward_check.pair_pvalue` on a T3
   pair at res 25, period 12500.3 is not flagged (if it is, report to the council before loosening anything).
2. **Fix 2 (G2; tagfile.py, propagate.py).** Per-key tolerance: `tagger_MeasDesc_GlobalResolution` compares at
   rtol 1e-4 (in T2 the ticks are discrete, so this is still exact; in T3 it admits a laser's drift and
   refuses a different repetition rate); the same tolerance in `tagger_settings`' within-folder comparison
   (replace `settings != found` by a per-key `isclose`). Keep the key name (old artifacts stay checkable).
   Print every mismatching value with `:.10g` in `condition_mismatches` and name both values in the folder
   ConfigError. Tests: `realhead.py`'s 10 Hz-apart folders pass both checks; 80 vs 76 MHz is refused; the
   message shows `1.250002031e-08`.
3. **Fix 3 (G3; adapters/ptu.py, tagfile.read_tags, analysis/hom.py, fringe.py, chsh.py, cli.py, README).**
   `read_tags` carries `sync_period_ps` from T3 metadata into the TagSet metadata it returns. When any analysed
   run has it: refuse the corrected value (`GuardrailRefusal`: "T3 input from a sync source: the side-window
   accidentals assume CW light; pulsed accidentals sit in side peaks. Declare `source.type: cw` in the setup
   if the source is CW and the sync is a clock only") unless the setup or call declares `source_type="cw"`
   explicitly; the raw result gets the label "T3 input: DADL models CW light (see README)". Do not route
   `rep_period_ps` into `side_offsets` now (non-integer periods drift; pulsed model deferred to v0.2+). L1.
   Tests: the `pulsed.py` 76 MHz case (reduced to 2e6 pulses) refuses the corrected V; with `source_type="cw"`
   declared it computes; a T2 file is unaffected.
4. **Fix 4 (G4, G5; adapters/ptu.py, tagfile.read_tags, cli.py).** (a) In T3, when `2**dtime_bits x res <
   period`, add a `t3-range` event and a metadata note ("dtime range X ns < sync period Y ns: photons late in
   each period are not recorded"); `read_tags` refuses such files (`AdapterError`, remedy: a coarser
   resolution). (b) T3 with `HWSync_Divider > 1`: `calibrate dark`, `calibrate deadtime` and `calibrate dcr`
   refuse it ("the T3 sync divider convention is unverified; rates may be off by the divider"), and analyses
   label it. Tests: `edges.py`'s 100 ns HydraHarp file is refused by `read_tags` with the event recorded by
   `PtuFile`; a divider-2 T3 file is refused by `calibrate dark` and labelled by `analyze`.
5. **Fix 5 (G6; propagate.py).** In `_check`, when the artifact's `instrument` and the dataset's both carry
   `model` (and `serial`), a difference is a condition mismatch for `delay-jitter` and `deadtime-series`
   kinds (refused, overridable with the existing `--allow-condition-mismatch REASON`), and a label
   ("calibration taken on another tagger") for other kinds. Test: a delay artifact from serial X refused on
   serial Y, passes with the override and a recorded reason; a dark artifact is labelled only.
6. **Fix 6 (S1; calibration/coverage.py, analysis/propagation.py, ci.yml).** (a) Recompute each gate in
   `merge` from the code (as propagation's `gate(replicates, 1)`), ignoring the shard's band; keep the shard's
   band only to refuse a disagreement with the code. (b) Shards record the requested family and regimes;
   `merge` refuses unless the merged calibration names equal the set the code produces for them. (c) An empty
   merge (no results) raises. Same three in `propagation.merge`. Tests: `empty.json` and `lax.json` are
   refused by both modules; kill `v_merge_no_tile` (Fix 9).
7. **Fix 7 (S2; artifact.py, cli.py).** Before `np.load`, parse the header with `np.lib.format.read_magic` and
   `read_array_header_1_0/2_0`; require dtype `<f8`, 2-D shape `(n_samples, len(params))` from the manifest,
   `fortran_order` False, and file size = header + rows x cols x 8; check the size before `read_bytes`.
   Catch `MemoryError` as `ArtifactError`. `_dark_artifact_for` then skips the file as intended. Test: the
   reviewer's `npy.py huge-header` table gives `ArtifactError` from `read_artifact`, `CalibrationSet.load`, and
   is skipped by `_dark_artifact_for`.
8. **Fix 8 (S3-S5; docs and repo hygiene).** L2 in ptu.py; `persist-credentials: false` on both
   `ppl-bench.yml` checkouts; `git rm bench/ppl/models/dark bench/ppl/models/delay` and `.gitignore`
   `bench/ppl/models/*` with `!bench/ppl/models/*.stan`.
9. **Fix 9 (Q1; tests/test_propagation.py, tests/test_calibration_coverage.py new).** Table tests on synthetic
   `CoverageResult`s: clean regime at 94% passes; control regime with control inside the band fails
   (`p_no_control_gate`, `p_control_half_wide`); flags at 2% in a clean regime fail (`p_max_flag_1`); misfits
   at 3% fail (`p_max_misfit_1`); misspecified regime at 89% flagged fails (`p_min_flag_0`). Coverage:
   `at_least` at 94.9% fails and at 100% passes (`v_np_not_at_least`, `v_at_least_two_sided`); a
   `deadtime_coverage` run with a planted miss counts it (`v_np_always_covered`); coverage `merge` with a gap
   and an overlap raises (`v_merge_no_tile`); a regime table test pins truth and the passive-quench dark
   artifact held (`v_dark_ap_held`). Replace the restating `passed()` test.
10. **Fix 10 (Q2; tests/test_hom.py, test_forward.py, test_chsh_corrected.py).** One HOM, one fringe and one
    CHSH test with a (B, A) delay-jitter artifact, B 400 ps late, declared delay 400 ps, asserting the
    corrected value equals the (A, B) run's (as `p_sign.py`). Must kill `f_capture_sign`,
    `f_pair_capture_sign`, `c_delay_start_sign`, `m_declared_swap`, `r_declared_swap`, `t_declared_swap`.
11. **Fix 11 (Q3; tests).** (a) `chsh()` with a loose dark set below the 20% dead fraction asserts
    `DEAD_TIME_BELOW_DARK` in `corrected_labels` (`h_drop_gap_labels`); (b) CHSH twin of
    `test_an_uninformative_dead_time_bound_is_refused` (`f_chsh_no_loose_bound`); (c)-(e) one assertion each
    for CHSH threshold labels, `JITTER_BOUND`, HOM/fringe `DEAD_TIME_RANGE_WIDE` (`h_no_threshold_labels`,
    `f_chsh_no_jitter_bound`, `t_range_wide_x10`); (f) fringe twin of
    `test_hom_checks_the_dark_dead_times_against_the_analysed_gaps` asserting labels, `forward_check_p` and
    `dataset` (`r_drop_use_labels`, `r_no_check_p`, `t_no_dataset`); (g) unequal dead times in
    `test_live_fraction_closed_form` and unequal run rates in the dead-fraction test (`f_vis_rate_min`,
    `f_live_max`); plus unequal dip/reference singles (`m_ref_side_from_dip`, `r_high_side_min`).
12. **Fix 12 (Q4; tests).** (a) `test_corrected_s_removes_the_dead_time_selection`: assert `raw.value > TRUTH`
    and the shift `raw - corrected` within 3 se of its expectation from the dead-time model; raise duration until
    the selection is >= 4 se if the expectation is too small to assert. (b)
    `test_corrected_hom_recovers_the_injected_visibility`: add the delay artifact (capture modelled), or move
    the coverage assertion to a slow 20-seed sweep with >= 17/20; keep one-seed sanity (finite, label
    present). Acceptance: each re-run over 20 seeds passes >= 19/20 (record in the PR).
13. **Fix 13 (Q6; tests/test_calibration_deadtime.py).** Behavioural oracles: 2 ln LR at both profile ends of
    tau_p equals 3.84 to 0.05 (`d_profile_90`); `est_p.high <= est_np.high` on a scan where the profile reaches
    the shortest gap (the tau_p cap mutant); 999 gaps refused (`MIN_GAPS` mutant); tau_np draws spread >= one
    tick (tick-jitter mutant); edge density with and without afterpulsing differ (afterpulse-term mutant);
    a union scan whose interval spans both branches (`d_no_union`, `d_union_no_widen`); afterpulse widening
    widens (`d_no_ap_widen`); a gap shifted by one ps moves tau_np by one ps (`d_shift_m_off1`).
    `test_the_tau_p_budget_is_recorded` stays, docstring "characterisation test (values from one run)".
14. **Fix 14 (Q7; tests/test_forward_check.py).** A known histogram whose p lies in (1e-6, 1e-3) is flagged
    (`MISFIT_ALPHA`); `CheckResult({"a": 0.01, "b": 0.5}).pvalue == 0.02` (Bonferroni); one bin set with
    expectations between 1 and 5 is grouped (`MIN_EXPECTED`); the dof of a known grouped histogram (dof mutant).
15. **Fix 15 (Q8; tests).** `_bounded_set` built with `_dark(bound=...)` deterministically, no `pytest.skip`;
    a CLI test with two dark artifacts in `--dark` (two channels, and two for one channel) asserting the right
    one is taken and the duplicate refused (`l_dark_any_channel`, `l_dark_first_of_many`); `calibrate
    deadtime` with mixed time units refused (`l_deadtime_no_time_units`).
16. **Fix 16 (Q9; tests, ci.yml).** Mark the four named tests and any test over 15 s `slow`; matrix jobs run
    `-m "not slow and not coverage"`; one Linux job runs `-m slow`. No test is dropped from CI.
17. **Method versions and gates.** deadtime `METHOD_VERSION` "6" -> "7" (Fix 1 changes `deadtime.py`; one
    bump), regenerate `tests/method_versions.json`. `delay.py`, `dark.py`, `accidentals.py` untouched (stay "8",
    "7", "3"). Re-run: `coverage deadtime --replicates 1000` on the CI seed, both regimes (simulator tags are
    on the 1 ps grid and no tick is declared: figures must equal round 6's exactly; any change is a bug);
    the propagation gate and `test_published_results.py` (Fixes 2, 3, 5 must not refuse or relabel any
    simulated or published input: published inputs are T2 or simulated); both merges via the CI shard path
    after Fix 6 (same verdicts). Re-run the r7-test harness (`mutate.py`) on the new HEAD and report the
    survivors named in Fixes 9-15 as killed; equivalent survivors listed in r7-test.md need no test.

## Owner items

1. **Release requires gates (S6).** Whether `release.yml`'s `check` job must require a successful "Statistical
   gates" check run on `$GITHUB_SHA` (council recommends it: `gh api repos/{repo}/commits/$GITHUB_SHA/check-runs`),
   or the policy "tag only commits whose main CI completed" is documented in RELEASING; and whether the ruleset
   forbids direct pushes to main. Release workflow and rulesets are the owner's.
2. **Committed binaries in history (S5).** Fix 8 removes them from the tree; purging them from history is a
   history rewrite and the owner's call (the council does not recommend it: 4.8 MB, never shipped).
3. **Behaviour change (G3, G4).** T3 users lose corrected values unless they declare `source.type: cw`, and T3
   files with a dtime range shorter than the sync period are refused; release notes should say so.
4. **Carried:** round-6 owner items 1-3 (dev version still 0.1.1 while deadtime goes 6 -> 7; schema `$id`; D4).
