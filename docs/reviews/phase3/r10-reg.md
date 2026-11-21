> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Interactions and regressions (round 10)

HEAD e8c1aba, repo venv, one BLAS thread. Probes copied and adapted in `scratchpad/council/r10-reg/p/`
(subdirs `r1`, `r2`, `r4`, `r4rse`, `r6`, `r7`, `r8`, `r9c`, `r9red`, `r9tw`, `dt`); gate logs in
`r10-reg/gates/`; the chsh-spad drift scan in `r10-reg/drift/`; unit suite log `r10-reg/unit.log`
(659 passed, 8 skipped, `-m "not slow and not coverage"`). The repo tree is unchanged (`git status` clean).

37 adopted fixes from rounds 1-9 were re-run against HEAD with the original probe or the council's probe.
36 still hold. One holds, but a later fix undercuts it (finding 1). The other findings are conflicts
between fixes, a label that breaks round 8's own label rule, and the gate's quick version.

## Findings

| # | Severity | Finding | Evidence | Where | Suggested fix |
|---|---|---|---|---|---|
| 1 | questionable | **Round 1's tau_p cap and the round 6/7 tick widening conflict: on coarse-tick data, a nearly fully paralyzable detector's tau_p interval misses the truth.** Round 6 (T1) and round 7 (G1) widened tau_np by one tick on each side, because the true dead time can sit up to one tick above the shortest recorded gap m. tau_p's profile and its cap (round 1 S6) still stop at `tau_hi = m`. When tau_p is within a tick of tau_np and the truth is off the grid, the tau_p interval ends below the truth. | `p/r1/cap_offgrid.py` (tau_p = tau_np = 22000 + U(0, tick), tags floored to the tick, 2/5/10 Mcps, 20 ms, 20 seeds). **tick 25: tau_p covered 2/20; tick 80: 4/20.** tau_np covered 20/20 in both. Typical miss: truth 22016, tau_p [19705, 22000], tau_np [21975, 22025]. Controls: tick 1 gives 19/20; tick 80 with tau_p = 0.8 tau_np gives 18/20. CI regimes are unaffected (1 ps grid). README L3 (round 1) says passive quenching "close to fully paralyzable" fits as tau_p near tau_np, which is exactly this case on PTU T3 or other coarse taggers. | `calibration/deadtime.py:686` (`tau_hi = float(m)  # tau_p <= tau_np`); the widening at `:789` | Let the tau_p profile and cap run to m + tick, i.e. to the widened tau_np upper end. Alternatively, when the profile reaches tau_hi, set high = min(high + tick, est_np.high). Bump the deadtime method version. Test: `cap_offgrid.py 25 20` covers in at least 18/20. |
| 2 | questionable | **`DEAD_TIME_ABOVE_GAPS` (round 8, Fix 7) breaks round 8's own label standard, and it duplicates `DEAD_TIME_BELOW_DARK` (round 6).** Round 8's budget says a label whose clean-data rate exceeds 1e-3 must state that rate in its text; README:390-391 says the labels' "text says how often". This label fires on most clean analyses but says nothing about how often. It also fires next to `DEAD_TIME_BELOW_DARK` on the same channel: the latter says the correction is "likely too large", the former says it "uses the dark value". | Fires on clean data: chsh-dead-time 400/400 replicates (r8-council W4); `m_clean100` (r9 red team's echo-free control) on all 4 channels; the r9 twin's clean SNSPD HOM (T3, dark bound 2.33 us vs analysed edge 40 ns) on both channels. Same channel, both labels: r6 SNSPD seed 113 and seed 58 at about 4 kcps (`p/r6/snspd.py`) give 4 dead-time labels for 2 channels. | `forward.py:77-80`, `:476-481` | Put its clean-data rate in the text ("fires on most analyses with a short dark run"). Do not add it for a channel that already carries `DEAD_TIME_BELOW_DARK` or a gap problem. |
| 3 | minor | **Latent conflict: round 6's T2 key names against round 9's Rd4 label.** Round 6 named the dark run's conditions `dark_temperature_K` and `dark_bias_V` so the use-time checks would ignore them. Round 9's `uncompared_conditions` labels any recorded key containing "temperature". Once M5a loads deadtime-series artifacts, every use would print "condition dark_temperature_K recorded, not compared". Unreachable today: analyses skip series artifacts (`SERIES_UNUSED`). | `p/dt/t2v1.py`: `uncompared_conditions(series.conditions, {"temperature_K": 233})` returns `['dark_temperature_K']`. | `propagate.py:207-212`; `calibration/deadtime.py:654` | Exclude `dark_*` keys in `uncompared_conditions`, or add this to the M5a entry items. |
| 4 | minor | **The gate's quick version cannot decide the control rule.** At 100 replicates the band is 87-100%. chsh-dead-time's control covers at about 87.5%, so the regime prints FAIL whenever control lands at 87% or above. | Seed 777001, 100 replicates: **chsh-dead-time FAIL** (corrected 98.0%, control 93.0%). The same seed at 400 replicates: **PASS** (corrected 96.8%, control 87.5%, band 91.0-98.25%). Binomial chance that control falls below the band, at a true control of 87.5%: 0.37 (100 replicates), 0.77 (200), 0.98 (400). | `analysis/propagation.py` `passed` / `describe` | In `describe`, print "control rule not decidable below about 400 replicates" when the replicates are fewer. Or document in CONTRIBUTING that quick runs judge corrected coverage only. |
| 5 | minor | **Display regressions from round 6 (Fix 7, Fix 10).** (a) A calibration part that is exactly zero now prints as floating-point noise: "calibration 4.4e-16" or "2.2e-16". The fix replaced "0.0000" with a number that implies an uncertainty. (b) The 2-significant-figure rounding prints an identified dark dead time (lower end 281 ns) as "16000000 (95% interval [0, 16000000])". That reads as an unidentified upper limit, the very distinction round 6 R1 relies on. | (a) `m_darkonly` with `dead_time_contradicted`; `r2-prop/p1_bound.py`, `p8_threshold.py`, `p6_chsh_bound.py` reruns. (b) `calibrate_dark` on the round 6 SNSPD seed 58: `low` = 281121, printed as 0. | `result.py` / `base.py` `__str__` (round 6 Fix 7, Fix 10) | (a) Print 0 when the part is below 1e-12 of the value. (b) Never round a positive interval end to 0; widen the precision of the ends instead. |

Not a finding, but worth watching: on the chsh-spad gate regime at seed 777001, 1 of 400 replicates was
"drift accepted" (replicate 62: z 3.77 against 3.72, diff 0.0056 against 0.5 se = 0.0053). Over the 400,
drift z has mean -0.02, sd 1.015 and KS p 0.55, so the statistic is calibrated. One event where about 0.08
are expected has chance about 8%. It stays inside the 1% flag cap.

## Gates on a non-CI seed (777001)

- **Propagation, all regimes, 100 replicates:**
  - hom-snspd: corrected 94.0%, control 94.0%, misfit artifacts 1.0%; PASS.
  - hom-bright: corrected 96.0%, control 96.0%; PASS.
  - hom-crosstalk: flagged 100%, corrected 29%; PASS.
  - chsh-spad: corrected 95.0%, drift accepted 1.0%; PASS.
  - chsh-dead-time: corrected 98.0%, control 93.0%; FAIL on the control rule (finding 4).
  - Zero refusals in every regime.
- **chsh-dead-time follow-up, 400 replicates:** corrected 96.8%, control 87.5%, misfit artifacts 0.2%; PASS.
- **chsh-spad, 400 replicates through the gate path (`drift/`):** covered 382/400 (95.5%); forward-check
  p < 0.01 in 1/400; off-peak p < 0.01 in 3/400.
- **Coverage deadtime, 200 replicates:**
  - series-spad: tau_p 95.5%, misfit status 0.0%, tau_np 100%.
  - passive-quench: tau_p 94.5%, misfit status 0.0%, tau_np 99.5%.
  - Band 89.0-99.5%; all ok.

## Checked and found sound (fix re-run against HEAD; probe in `p/`)

- **R1 S1:** the per-run interruption check gives ppc 0.49 ok for same-origin runs, runs offset on one
  clock, and 20/20/15 ms durations. A 20% hole still gives status misfit.
- **R1 S2:** `survival` against an 80-digit series at n tau 0.2-3.0, out to 40 tau: worst relative error
  4.4e-8.
- **R3 N3:** an int tau gives the same result as a float tau.
- **R1 S3:** the profile is smooth with the pending term on: sd of the second difference 0.007 (0.009 at
  h/10).
- **R1 S6:** value, high and every draw keep tau_p <= tau_np, with and without a declared tick (but see
  finding 1).
- **R2 R1 and B2:** the point sits at the bound and inside its interval. `p8_threshold` now carries "threshold
  crossed only after correction". CHSH reports the bound and zero ends.
- **R2 R3:** stand-in jitter no longer raises KeyError (0.9025).
- **R2 B1 and R3 G1/G2:**
  - Misdeclared delays (400 and 200 ps) carry `CAPTURE_IMBALANCE`; correct declarations do not.
  - Drift after calibration of 100-400 ps: refused 16/16 under the round 8 `DRIFT_Z` 3.72.
- **R2 H2:** with calibrations and no reference delay, `ValueError`.
- **R2 H1:** a +6% pump drift carries `SINGLES_DRIFT`.
- **R3 N2, R4 DT1-DT3, R5 I1:**
  - Afterpulse shift of 331-459 ps on pq seeds 0, 20, 24, 25 (round 4: 309-439).
  - Past-maximum scan covered 23/24 (round 4: 47/48), with the decisions split between switched and union
    as before.
  - The replicate-477 union test passes in the unit suite.
- **R4 DT4, R8 St1, R9 Rd1 (gap tests through `check_dead_time_gaps`):**
  - Refused: 22 ns on a 50 ns detector at 1 Mcps (upward); 50 ns dark on 46/44/40 ns detectors (downward);
    1 us on 30 ns at 3e7 cps.
  - Matched null: 0/80 refused.
  - Passes as stated in the docstring: 22 on 50 ns at 10 kcps.
  - r8-stat's seed 4242 replicates 69 and 360 run to completion.
- **R4 S1, S2, S3, S7:**
  - S1: a CHSH analysis with accidentals only is not refused, and `forward_check_p` is absent.
  - S2: a partial dark set gives a named `ValueError` on all three paths.
  - S3: the API detector check refuses detector X0 when analysing NEW0. It needs `dead_time_contradicted`
    allowed first, since the gap test now fires earlier on that probe.
  - S7: duplicate timestamps and a dark artifact for the wrong channel are refused.
- **R5 I2:** the pending law at 3 fresh seeds of each of the two stress cases (2e6 gaps):
  - Chi-square p with the pending term: 0.17-0.67.
  - Without it at 1e7 cps: p 1e-20 to 1e-27.
  - At 6e6 cps / 90 ns the term is not needed at 2.3e6 gaps (p 0.18-0.51 without it), as round 5 found.
- **R5 V1:** two dark artifacts give different series inputs.
- **R5 V5:** `valid_for.conditions` is enforced (temperature recorded, not given for the dataset: refused).
- **R6 R1(b) and Fix 1':** the 7 SNSPD false identifications are refused at 85 kcps by the 20% rule. At
  about 4 kcps they pass, labelled `DEAD_TIME_BELOW_DARK` (floor kind "lower end (label only)", still
  reachable after round 8's exponential floor).
- **R6 T1:** off-grid tau_np covered 98/100 (tick 25), 60/60 (tick 80), 100/100 (tick 4).
- **R6 T2:** dark artifact conditions differing from the scan are refused; with the override, a note.
- **R7 G1:** T3 at 25 ps / 12500.3 ps covered 20/20 with tick 25 or 26 passed. Without a tick, 4/20, so the
  API needs the tick, as ruled.
- **R7 G2:** sync rates 10 Hz apart are accepted; 80 vs 76 MHz is refused; values print with 10 digits.
- **R7 S1:** coverage and propagation merges fail closed:
  - refused: empty results, a lax band, an unknown regime, and old shards without a family;
  - a full-tiling shard below the band gives exit 1;
  - a passing shard gives exit 0.
- **R7 G3 and R9 Rd3:** T3 with `spdc-cw` is refused with the new text; `source_type=None` gives a raw value
  only.
- **R8 Ph1 and R9 Tw1:**
  - Delay at lattice 25 with the T3 period: 0/10 misfit, jitter covered 10/10.
  - `res80.py` at 3 fresh seeds: every period, including 49999.37 and 12500.3, gives a delay interval
    containing 4372 with status ok.
  - The twin's T3 HOM folder calibrated by the HEAD CLI: delay 4371.95 [4371.57, 4372.33], p 0.80.
- **R9 Tw2:** a T2 dark artifact is accepted on the twin's T3 analysis.
- **R8 Ph3 and R9 Pf1 (off-peak check):**
  - Pulsed HOM (T2/plain tags) refused at 12500.37, 13157.9 and 12345.6 ps; not at 12500.0.
  - CHSH (16 sides, via `offpeak_check`): pulsed share 0.5 and 1.0 refused 15/15 at both periods; 0.3
    refused 2/15 and 11/15 (L4 says 8-17 of 30).
  - The perf probe's every-4th-window side is refused.
  - Null, 3,000 sixteen-side analyses: 0 refused, adjusted p < 0.01 in 0.9%.
  - A refused analysis costs 3.7-4.7 s.
- **R9 Rd1, Rd1b, Rd2 (CLI on the red team's data):**
  - `m_echo_few`, `m_echo100` and `m_darkonly` are refused on 4/4 channels.
  - `calibration_outside_validity` no longer unlocks the refusal.
  - `dead_time_contradicted` gives the accepted label with the problems, plus the "below half the declared"
    label.
  - `m_clean100` passes: S 1.9937 [1.9875, 2.0038].
- **Thresholds across modules:** they agree with README "Refusals and their rates" and the round 8 budget:
  - forward check 4e-4; drift 3.72 (refuse) and capture label 3.29;
  - off-peak 2e-4; gap tests 1e-4 + 1e-4;
  - calibration-time `MISFIT_ALPHA` 1e-3, a separate class.
- **Checks still reachable after later rounds:**
  - the 3-width floor (artifacts without `gaps` only);
  - the lower-end label (SNSPD, chsh-dead-time);
  - the p = 0 "either rate" note;
  - `DEAD_FRACTION_HIGH`, kept for old records as documented.
