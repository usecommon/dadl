> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 8 (physics as a whole, statistical design as a whole, plan)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: Ph = r8-phys.md rows 1-7, St = r8-stat.md rows 1-7, Pl = r8-plan.md rows P1-P8. Council probes
and outputs are in `scratchpad/council/r8-council/`:
- `w4.py`, `w4run.sh`, `w4jobs.txt`, `w4sum.py`; raw `w4/*.jsonl`, summary `w4sum.txt`
- `bud.py`, `budrun.sh`, `budjobs.txt`, `budsum.py`; raw `bud/*.jsonl`, summary `budsum.txt`
- copies of `t3delay.py` and `pulsed_t2.py`, with outputs `t3delay.txt` and `pulsed_t2.txt`

## Checked numerically

All runs used HEAD 8d7c50a with the repo venv, one BLAS thread, and the package's own functions.

### W4 (Ph4, with Ph5): truncating the dark dead times at the analysed shortest gap

`w4.py` runs the propagation gate's exact replicate path: `_children`, `P._calibrate` and `P.POLICY`, with
`calibration_outside_validity` also allowed so that nothing aborts. Each replicate is analysed four times with
`forward._channel_dead_times` / `forward.live_fraction` patched:
- **ship**: as shipped.
- **trunc**: every draw, point and bound is clipped at g_min + 1 ps per channel (the simulator's lattice; the
  simulator has no tagger jitter after the dead time).
- **add**: the additive pair form (1-fa)(1-fb) + r_c·min tau.
- **both**: trunc and add together.

CI seed 20261118, replicates 0-399 (chsh-dead-time) and 0-199 (hom-bright):

| regime | variant | corrected | control | mean z (se) | rms z | se / se_stat |
|---|---|---|---|---|---|---|
| chsh-dead-time (400) | ship | **95.5%** | **87.5%** | -0.43 (0.06) | 1.25 | 1.18 |
| | trunc | 92.5% | 93.0% | +0.42 (0.05) | 1.10 | 1.10 |
| | add | 96.3% | 80.2% | -0.99 (0.06) | 1.60 | 1.23 |
| | both | 94.0% | **94.0%** | +0.13 (0.05) | 1.03 | 1.12 |
| hom-bright (200) | all four | 95.0% | 95.0-95.5% | -0.05 to +0.04 | 1.02 | 1.00 |

The band is 92.6-97.1% at 1,000 replicates; the control must fall below 92.6%.

- Truncation changes S in 400/400 chsh-dead-time replicates, by +0.85 stat se on average. This confirms Ph4's
  48/48: the dark points sit above the analysed edge on essentially every channel.
- **Truncation alone fails the gate.** Corrected sits at 92.5%, at the band's lower edge (expected to fail about
  half of 1,000-replicate runs). The control is 93.0%, not below the band. The +0.42 z bias left behind is
  Ph5's pair-term under-correction, which today's over-long dark dead times cancel (ship: -0.43).
- **The additive form alone is worse than shipped**: z -0.99, control 80%. It must not land without truncation.
- **Both together** are nearly unbiased (z +0.13, rms 1.03). Physically, this is the right model. But the
  regime then pins tau to about 1 ps from 1.2M analysed gaps, so the dead time's uncertainty no longer matters.
  The control covers as well as the corrected value (94.0% / 94.0%). **chsh-dead-time stops being a control
  regime**, and the gate as declared fails on the control rule.
- hom-bright does not move in any variant.

Conclusion: W4 is right physics, but it cannot be adopted on its own. It needs Ph5 in the same change, and a
redesigned control regime: one whose dead-time uncertainty survives truncation, for example a low analysed rate
where the gap edge is not sharp, or a series artifact at use. Redesigning a regime is a gate change for the owner.

### Refusal budget on clean data (St1, St4, St5, Ph3)

`bud.py` uses the gate path with a permissive policy and fresh seed 777. Runs: chsh-spad 400, hom-snspd 400,
hom-crosstalk 200 (misspecified, for power). It records every statistic, and the delay ppc under both dof rules.

| check | chsh-spad (400) | hom-snspd (400) | hom-crosstalk (200, power) |
|---|---|---|---|
| delay ppc, groups-5 (shipped): p<0.05 / p<0.01 / p<1e-3 | 5.6% / 1.25% / 1 of 1,600 | **9.5%** / 1.25% / 0 of 400 | 13.5% / 6.5% / 1 |
| delay ppc, groups-1: same | 3.4% / 0.56% / 0 | 2.5% / 0.25% / 0 | 8.0% / 3.0% / 1 |
| forward check adj p: <4e-4 / <1e-3 / <0.01 | 0 / 0 / 0.50% | 0 / 0 / 0 | **200/200 / 200/200** / 100% |
| drift z: mean, sd, KS p | +0.06, 0.98, 0.53 | n/a | n/a |
| drift: expected refusal rate at 3.29 / 3.72 sd (each run's own max(0.5 se, k sd)); 0.5 se binds | 5.0e-4 / **2.0e-4**; 96% | n/a | n/a |
| off-peak p, Bonferroni over acquisitions (16 / 2 / 2), exact dispersion at n_sim 20,000: <2e-4 / <0.01 / <0.05 | 0 / 0.25% / 4.0% | 0 / 0.25% / 2.25% | 0 / 1.5% / 7.0% |
| off-peak p per acquisition: <0.05 / <0.01 | 4.7% / 1.05% | 3.9% / 0.25% | 6.3% / 2.0% |
| coverage, corrected / control | 93.5 / 94.0% | 94.5 / 94.8% | 29.5% (misspecified) |

- The off-peak check is calibrated on clean data, including chsh-spad's afterpulsing (0.5%, 50 us), and is
  conservative in the tail. The mean side count is 1.3 per window. At that count, today's `offpeak_pvalue`
  (exact dispersion at n_sim 2,000) cannot go below about 1e-3 after doubling, so it **cannot fire at
  alpha 2e-4**. The analysis-time check needs n_sim >= 20,000.
- Pulsed T2 (`pulsed_t2.txt`, Ph3 reproduced): at 12500.37 / 13157.9 / 12345.6 ps the side estimate is
  813 / 372 / 173 against about 1,870 accidentals. Corrected V is 0.973 / 0.962 / 0.957 +- 0.0011 (truth 1) with
  no label, and the off-peak check gives p = 0 in all three. At 12500.0 ps it gives p = 0.45, and there the
  estimate is right.
- Forward check at alpha 4e-4 keeps full power on the misspecified regime (200/200; MIN_FLAG_SHARE is 90%).
- Delay dof (St5) reproduced: groups-5 runs hot at p<0.05 (9.5% in hom-snspd, binomial p about 1e-4 against 5%).
  groups-1 is conservative, as `forward_check` already is.
- Downward gap floor (St1). On r8-stat's 12,000 chsh-dead-time dark artifacts (`r8-stat/dt/`), the floor
  point - ln(K/alpha)/(gaps x dark rate), with K = 4, alpha = 1e-4 and 1/(N f) = 25 ns, sits 265 ns below the
  point. It is **above the truth in 0/12,000** artifacts (0.6 expected at alpha/K) and above zero in 3. The
  3-width floor misses in 3/12,000. At 1 s and 50 kcps the same law gives a floor 4.2 ns below the point, so it
  keeps its power where the dark data have it.

### T3 tick as histogram lattice (Ph1, St3), reproduced

`t3delay.txt`: 25 ps T3 at a 12500.3 ps period, 10 seeds per tick.

| tick passed | misfit | jitter covered (truth 141.4) | delay median |
|---|---|---|---|
| 26 (what the CLI passes today, `tag_tick_ps`) | 10/10 | **0/10** (136.7) | **-12.1 ps** |
| 25 | 0/10 | 10/10 | -0.1 |
| none (1 ps from the tags) | 10/10 | 10/10 | — |

At period 12500.0, tick 26 is still 10/10 misfit, and ticks 25 and none are both clean. Code read:
`ptu.py:416-426` returns res + 1 when the period is not a whole number of ps; `cli.py:812` passes it to
`calibrate_delay`; `forward_check.py:237,304` take `tag_tick` from the tags.

## The refusal budget (one standard for every check)

**Standard.** On clean data within the declared regimes, an analysis is refused at analysis time with
probability at most **1e-3 per analysis**. This holds for all checks together, whatever the number of channels,
pairs or acquisitions: each check corrects for its own multiplicity with Bonferroni. It replaces round 6's
"1e-4 per channel for a downward refusal", which was stated for one test and was not met (St1).

| check | kind | share of the 1e-3 | threshold | measured on clean data |
|---|---|---|---|---|
| forward-model predictive check (`Misspecified`, overridable `misfit_forward_model`) | refuse | 4e-4 | `forward_check.MISFIT_ALPHA` 1e-3 -> 4e-4, Bonferroni over histograms | 0 / 800 below 4e-4 (r8-council), 0 / 2,000 below 1e-3 (r8-stat); power kept (200/200) |
| CHSH delay drift (`Misspecified`, overridable) | refuse | 2e-4 | max(0.5 se, 3.72 sd); `DRIFT_Z` 3.29 -> 3.72 (two-sided 2e-4) | expected 2.0e-4 from the measured z law |
| off-peak flatness (new, Ph3; `Misspecified`, overridable) | refuse | 2e-4 | `offpeak` p < 2e-4 / n_acquisitions, exact dispersion at n_sim >= 20,000 | 0 / 800 |
| dead-time gap tests, downward + upward (`StaleCalibration`) | refuse | 1e-4 + 1e-4 | downward: floor = point - ln(K / 1e-4) / (gaps x dark rate); upward: excess > 100((K / 1e-4)^(1/100) - 1) (10.4 at K = 2, 11.2 at K = 4) | downward 0 / 12,000 artifacts; upward 0 / 2,400 channels (r8-stat) |
| consequence rules (dead fraction > 20%, loose bound, T3 sync source, t3-range, condition mismatch) | refuse | none (deterministic, not tests) | unchanged | 0 in every declared clean regime; the gate records them |

**Calibration time (per artifact, a separate class).** The predictive check refuses to write a misfit artifact
at `MISFIT_ALPHA` 1e-3 per artifact (`calibration/base.py`, shared by the pinned methods; unchanged). A user who
gets this refusal re-takes a cheap calibration run; no analysed data are lost. With the delay dof fix, a
four-pair CHSH set (8 artifacts) meets it in at most about 0.8% of sets (measured: 1/400 replicates at groups-5,
0/400 at groups-1). Stated as such (limit L1).

**Labels.** Labels are not refusals. Any label may fire at up to 5% per channel or check on clean data. Where
its clean-data rate exceeds 1e-3, the label text states it (as `DEAD_TIME_BELOW_DARK` already does).
`INTERRUPTED` (1e-3) and `CAPTURE_IMBALANCE` (keeps 3.29 sd, as a new constant `CAPTURE_Z`) stay labels.

Verified (the two tables above): on clean chsh-spad and hom-snspd, the analysis-time checks under this budget
refused 0 of 800 analyses. The statistics are calibrated or conservative (drift z sd 0.98; forward and off-peak
tails at or below nominal), so the budget holds by construction rather than by a count that 800 runs cannot
resolve at 1e-3.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| Ph1 / St3 | T3 tick (res + 1) used as the histogram lattice: delay -12 ps, jitter 0/10 covered, every artifact a misfit | Y/Y/Y | ADOPT-FIX | Fix 1 |
| Ph2 | Forward check (and `pooled_peak`) misfire on T3: bin on the dtime lattice | Y/Y/Y | ADOPT-FIX | Fix 1 |
| Ph3 | Pulsed light on T2 or plain tags: V 10-40 sd off, unlabelled | Y/Y/Y | ADOPT-FIX + ADOPT-LIMIT | Fix 5; L2 |
| Ph4 | Dark dead times above the analysed shortest gap (400/400 replicates) | Y/N/N (fix now); Y/Y/Y (limit + defer) | ADOPT-LIMIT + DEFER (M5a, together with Ph5 and a redesigned control regime) | Fix 7 (label); L3; M5a entry item |
| Ph5 | Pair live fraction under-corrects at second order | Y/N/N (fix now); Y/Y/Y (defer with Ph4) | ADOPT-LIMIT + DEFER (M5a, only together with Ph4) | L3 |
| Ph6 | Gate regimes narrower than real labs | Y/Y/Y | DEFER (M4b / M3b regime sitting, owner) | O8 |
| Ph7 | The 20% dead-fraction refusal gives the wrong reason | Y/Y/Y | ADOPT-FIX (wording) | Fix 7 |
| St1 | Downward gap floor on the wrong scale (2.5e-4 to 3.7e-4 per channel) | Y/Y/Y | ADOPT-FIX | Fix 3 |
| St2 | The gate aborts on any refusal; the CI pass is seed luck | Y/Y/Y (record); owner (tolerance) | ADOPT-FIX (record refusals by type, verdict unchanged) + OWNER (tolerance, second seed) | Fix 6; O12 |
| St4 | No joint budget; labels and refusals at mixed alphas | Y/Y/Y | ADOPT-FIX | Fix 4 (budget above); L1 |
| St5 | `delay.py` dof groups-5 runs hot (9.5% at p<0.05) | Y/Y/Y | ADOPT-FIX (method 8 -> 9) | Fix 2 |
| St6 | Upward gap test's "1e-4 per channel" wrong | Y/Y/Y | ADOPT-FIX | Fix 3 |
| St7 | Drift and predictive-check events share a label and a key | Y/Y/Y | ADOPT-FIX | Fix 6 |
| Pl1 | M5's only consumer of the series artifact is optional | Y/Y/Y | ADOPT-FIX (plan text: M5a / M5b; neutral strings) + OWNER (M5a in the gate minimum) | Fix 8; O6 |
| Pl2 | M10 overloaded | Y/Y/Y | ADOPT-FIX (plan text: M2b, M3b, M4b, M10) + OWNER (which may slip) | Fix 8; O13 |
| Pl3 | Homeless deferrals (pair labels, D8, site) | Y/Y/Y | ADOPT-FIX (plan text) | Fix 8 |
| Pl4 | Owner items scattered and dropped | -/Y/Y | ADOPT-FIX (plan section) | Fix 8; consolidated list below |
| Pl5 | Slow-test job not blocking | Y/Y/Y | ADOPT-FIX (`gates` needs) + OWNER (required check) | Fix 9; O1 |
| Pl6 | PR #26 description stale; two PRs under one title | -/Y/Y | ADOPT-FIX (title and body) + OWNER (merge method) | Fix 9; O14 |
| Pl7 | Decision texts do not match what was built | Y/Y/Y | ADOPT-FIX (D3 as-built note, lines 8 and 200) + OWNER (D4 text) | Fix 8; O5 |
| Pl8 | Dependency table drift | -/Y/Y | ADOPT-FIX | Fix 8 |

**Votes on Ph4/Ph5.**
- P: a dead time longer than every analysed gap is excluded by the definition of a non-paralyzable dead time.
  "both" is nearly unbiased (z +0.13, rms 1.03), so I vote to fix now with the regime change.
- M: the measurements show that the fix cannot land alone. Truncation without Ph5 gives 92.5% and a 93.0%
  control, which fails the gate. The pair kills the declared control regime (94.0% / 94.0%). Changing the
  control regime is a gate change that the owner pre-registers (O8). Today's corrected value covers
  (95.5%, round 6: 95.5%). Defer to M5a as one change: truncation, the additive form, and a control regime whose
  dead-time uncertainty survives truncation.
- U: the shipped interval covers, so the user is not misled about uncertainty. The user is misled about the
  dead time. A label now saying so, plus the M5a route, is enough.

**Votes on St2.**
- All three: recording refusals by type instead of aborting is pure measurement, and it is adopted.
- M: whether a clean regime may pass with some refusals, and how many, changes a gate's pass condition. Since
  round 1 the regimes and gate rules have been the owner's to declare, so the owner must agree to the tolerance.
  Until then "any refusal fails a clean regime" stays (same verdict as the abort today). Fix 3 makes refusals
  in the gate's policy about 1e-4 per replicate, so the seed-luck risk largely goes away with it.
- P, U: agree.

**Votes on Pl4.** P abstains (process).

## Limit wordings

- **L1** (README, new "Refusals and their rates" section; `guardrails.py` module docstring): "On clean data within
  the tested regimes, DADL refuses an analysis at analysis time with probability at most 1 in 1,000, all checks
  together: the forward-model check 4e-4, the CHSH delay-drift check 2e-4, the off-peak flatness check 2e-4,
  and the dead-time gap tests 2e-4. Each check corrects for its own number of channels or acquisitions.
  Separately, a calibration whose model check fails (p < 0.001 per artifact) is not written; a CHSH calibration
  set of eight artifacts meets this in under 1% of sets, and the remedy is to re-take that calibration run.
  Labels are not refusals: some fire on a few percent of clean analyses, and their text says how often."
- **L2** (README line 18 limit and `adapters/ptu.py` T3 paragraph, extending round 7's L1): "Pulsed sources are
  not supported on any input (T2, T3 or plain tags): pulsed accidentals sit in side peaks at multiples of the
  repetition period. On T2 and plain tags DADL detects them through the off-peak windows (not flat) and
  refuses the accidental-subtracted values; a period that is a whole multiple of the window grid is not
  detected, but then the off-peak estimate is also right."
- **L3** (README dead-time limits; `forward.py` docstrings of `live_fraction` and `_channel_dead_times`): "The
  forward model takes each channel's non-paralyzable dead time from the dark artifact. At high analysed rates
  the analysed gaps measure it more precisely, and often show it to be shorter than the dark value (label
  `DEAD_TIME_ABOVE_GAPS`). The pair live fraction is first order: the term r_c·tau·(fa + fb) is not included. In
  the chsh-dead-time regime the two effects cancel: the corrected S covers at 95.5% with mean z -0.43 statistical
  se. Neither effect is corrected alone, because each alone biases the value more (council round 8). Reading the
  dead time from the analysed gaps, together with the second-order term, is milestone M5a."

## Implementation order (adopted fixes)

1. **Fix 1 (Ph1, St3, Ph2; adapters/ptu.py, tagfile.py, data.py, cli.py, forward_check.py, analysis paths).**
   - (a) `PtuFile.lattice_ps`. T3: the dtime resolution when it is a whole number of ps, else None. T2: the
     whole-ps global resolution. `tick_ps` keeps round 7's meaning (error bound, res + 1 with rounding) and is
     used only for dead-time intervals. `tagfile.tag_lattice_ps(path)` mirrors `tag_tick_ps`.
   - (b) `calibrate pairs` passes `tick_ps=tag_lattice_ps(folder)` to `calibrate_delay`.
   - (c) `read_tags` carries `lattice_ps` in the TagSet metadata. `PairData.lattice_ps` (default None) is filled
     by the analysis path. `pair_pvalue` and `pooled_peak` use it when given, else `tag_tick`. They align the
     bin edges half a lattice step off the modal residue of (b - a - shift) mod lattice in the window, so that
     +-1 ps rounding never crosses an edge.
   - (d) T3 with a resolution that is not whole ps: skip the forward check with a note ("no integer lattice").
     Do not bin at 1 ps.
   - (e) The delay off-grid note does not fire on T3 input for residues of at most 1 ps plus (period mod res).
     This is a delay.py change; it lands with Fix 2's bump.
   - Tests:
     - The strict xfail becomes a pass at lattice 25 (12500.3 ps).
     - Add 16 ps / 12345.6 ps, and a cable delay that spans syncs (r8-stat `t3check2.py` cases): p<0.05 in at
       most 4 of 40 seeds each (slow-marked; one seed fast).
     - `t3delay.py` at lattice 25: 0/10 misfit and jitter covered in at least 19/20 seeds.
     - A T3 folder at 25 ps / 13157.9 ps through `calibrate pairs` gives a delay interval that contains 0.
     - `tag_lattice_ps` returns 25 on that folder and `tag_tick_ps` still returns 26.
2. **Fix 2 (St5, Fix 1e; calibration/delay.py).** `_grouped_pearson(..., n_params)` call with 1 (groups - 1, as
   forward_check). Add the T3 note condition. `METHOD_VERSION` "8" -> "9"; regenerate `tests/method_versions.json`.
   Test: 400 hom-snspd pair runs at groups-1 give p<0.05 in at most 6% (council: 2.5%) (slow).
3. **Fix 3 (St1, St6; forward.py, README).**
   - Downward floor: point - ln(K / `GAP_DOWN_ALPHA`) / (gaps x dark_rate.low), with K the analysed channels
     and `GAP_DOWN_ALPHA` = 1e-4. `gaps` comes from the artifact's `conditions["gaps"]`; an artifact without it
     keeps the 3-width floor. Floor kind "exponential law". The count rule is unchanged. The lower-end label
     (round 6 Fix 1') stays when the floor is at or below zero.
   - Upward: `GAP_EDGE_EXCESS` = 100((K / 1e-4)^(1/100) - 1) per analysis.
   - Docstring and README L1 numbers as in the budget table. No dark.py change ("7").
   - Tests:
     - The two r8-stat refusals (chsh-dead-time seed 4242, replicates 69 and 360) run to completion.
     - A 1 us artifact on 30 ns gaps with 1 s / 50 kcps dark data is still refused.
     - The floor for the darktail extremes (point 263 ns, gaps 2,000, 20 kcps) is at or below zero.
4. **Fix 4 (St4; forward_check.py, analysis/chsh.py, README, guardrails docstring).**
   - `forward_check.MISFIT_ALPHA` 1e-3 -> 4e-4.
   - `DRIFT_Z` 3.29 -> 3.72 for the refusal, and a new `CAPTURE_Z` = 3.29 for the `CAPTURE_IMBALANCE` label
     (unchanged behaviour).
   - L1 in the README.
   - Tests: Fix 14 of round 7 updated to the new alpha. A drift at z 3.5 above the 0.5 se floor is labelled
     only, not refused.
5. **Fix 5 (Ph3; analysis/_twopoint.py, hom.py, fringe.py, chsh.py, README).**
   - Each analysis computes, on each acquisition's shifted counts it already has, p = min(1,
     2 min(dispersion_pvalue(counts, n_sim=20000), trend_pvalue)). It flags when the smallest p < 2e-4 /
     (number of acquisitions or pairs).
   - Flagged: refuse the accidental-subtracted raw and corrected values with `Misspecified`: "the off-peak
     windows are not flat (pulsed light or a periodic background): the accidentals cannot be read off-peak".
     Overridable with `misfit_forward_model`, and then labelled `OFFPEAK_ACCEPTED`.
   - `accidentals.py` unchanged ("3").
   - Tests:
     - `pulsed_t2.py` at 13157.9 ps (2e6 pulses) is refused.
     - At 12500.0 ps it is not refused.
     - The 20 covering chsh-spad replicates are not flagged.
     - Off-peak p on 200 clean hom-snspd analyses is below 0.01 in at most 2% (slow).
6. **Fix 6 (St2, St7; analysis/propagation.py, analysis/chsh.py).**
   - `_replicate` catches `DadlError` and returns a refusal record (type and message). `propagation` counts
     refusals under "refused: <type>", and a refused replicate counts as not covered.
   - `passed` is unchanged in effect: any refusal fails a clean regime until the owner sets a tolerance (O12).
   - The drift's override label becomes `DELAY_DRIFT_ACCEPTED`, counted under "drift accepted". "forward
     misfit" now counts only the predictive check. A new "off-peak accepted" key is added. The 1% cap applies
     to their sum, so the verdict is unchanged.
   - `describe` prints each key.
   - Tests:
     - A planted `StaleCalibration` in one replicate is counted and fails the regime, without a traceback.
     - Merge tiling keeps the refusal keys.
7. **Fix 7 (Ph4 limit, Ph7, Pl1 strings; forward.py, propagate.py, cli.py, README).**
   - New label `DEAD_TIME_ABOVE_GAPS`: "the dark dead time of {name} ({x}) is longer than {k} of the analysed
     gaps (shortest {g}); a non-paralyzable dead time cannot be. The correction uses the dark value and its
     uncertainty (see Known limits)". It fires under the count rule `GAP_COUNT` and goes in
     `corrected_labels`.
   - L3.
   - Reword the 20% refusal as Ph7 suggests.
   - Replace "(M5)" / "until M5" in user-facing strings with "not yet supported; see Known limits".
   - Test: a chsh-dead-time replicate carries the label.
8. **Fix 8 (Pl1-Pl4, Pl7, Pl8; docs/phase3-plan.md).**
   - Split M5:
     - M5a, series dead time at use: V4, T6, B5, the live fraction from `deadtime-series`, W4 truncation
       together with the Ph5 pair term, and a redesigned control regime. It depends on M2 and M4.
     - M5b, NUTS: B1-B4, B7. It depends on M0, M5a and the B1 decision.
   - Re-split M10 into M2b (D8 IRF tail), M3b (efficiency normalisation, D3 union, heralding with pair-labelled
     simulation truth: retarget the plan:358 row and review item 7.14), M4b (the pre-registered dead-time gates
     and regimes, coarse and Fano checks, the 0.5 s probe, Ph6), and M10 (guardrails job, docs including the
     site `docs/index.html`, review cycle, release).
   - Add an "Owner items" section with the list below.
   - Add an "As built" note to D3.
   - Line 8: "M0-M3 merged, M4 in review (#26)". Qualify line 200.
   - Dependencies: M8 depends on M1 and M6. M10 depends on the gate minimum plus the optional milestones that
     merged. The T3 divided file is an external dependency (O10).
9. **Fix 9 (Pl5, Pl6; ci.yml, PR #26).**
   - Add the slow-test job to `gates`' `needs`.
   - Rewrite the PR #26 title ("M4 series dead time, council fixes to M1-M3 corrected values, T3 and CI") and
     body: one section per area, current method versions and gate figures, a "Behaviour changes" list (O7), the
     owner items, and a commit-by-round guide.
   - Merge only after the HEAD run is green.
10. **Method versions and gates.**
    - **Pinned methods.** delay "8" -> "9" (Fix 2; one bump covers Fix 1e). dark "7", deadtime "7",
      accidentals "3", dcr "2" and heralding "6" are untouched; their sha256 must not change.
    - **Re-run `coverage delay`** at 1,000 replicates on the CI seed. Estimates and intervals are unchanged; only
      the ppc p moves, so coverage must equal the current figures exactly, and the misfit share falls.
    - **Re-run the propagation gate**, all five regimes at 1,000 replicates via the CI shards. Expected:
      - corrected and control figures identical to the current record (the changes move only checks and labels,
        not values);
      - zero refusals;
      - "forward misfit" at or below today's;
      - hom-crosstalk flagged at least 90% (council: 100% at 4e-4).
      Any coverage change is a bug.
    - **Also re-run:**
      - `coverage deadtime`: deadtime is untouched, but the T3 plumbing touches `tag_tick_ps` callers, so the
        figures must equal round 7's exactly;
      - `test_published_results.py`: Fix 5 must not refuse a published input; report any label or refusal,
        do not suppress it;
      - the unit suite, and the r7-test `mutate.py` harness on the new HEAD.

## Owner items (consolidated, rounds 4-8; proposed for the plan's new section)

| # | Item | Source | Status at 8d7c50a | Council recommendation |
|---|---|---|---|---|
| O1 | Make "Statistical gates" and the slow-test job required checks on `main` | r4 C1, r8 Pl5 | Open (11 contexts, neither among them) | Edit ruleset 27075146 |
| O2 | Release requires gates on the tagged SHA (or a RELEASING policy); forbid direct pushes to main | r7 S6 | Open | Require it in `release.yml`'s `check` job |
| O3 | Dev version (still 0.1.1; delay goes to 9, deadtime is 7) | r5 V2, carried r6, r7 | Open | 0.2.0.dev0 before more artifacts |
| O4 | Calibration schema `$id` | r5 V7 | Open | With O3 |
| O5 | D4 decision text ("0.1.x readers refuse") | r5 V6 | Open | Reword as r5 L2 |
| O6 | M5 scope (B1), now M5a / M5b; whether M5a joins the gate minimum | r6 B1, r8 Pl1 | Open | Split; M5a in the gate minimum |
| O7 | Accept the behaviour changes in release notes: R1(b) refusal above 20% dead fraction; T3 needs `source.type: cw`; `t3-range` refusal; **new:** off-peak refusal of pulsed light on T2/tags; forward check at 4e-4 and drift at 3.72 sd; the 1e-3 per-analysis budget as the stated standard | r6 item 1, r7 item 3, r8 | Partly in CHANGELOG | Accept |
| O8 | Pre-register gate and regime changes in one sitting: rotating seeds (with a second seed, St2), tau_np gate (r1 G1), D3 union with chsh-dead-time (r2), past-maximum with declaration and flip search (r3 N2, r4), unidentified tau_p (r6 T7c), 5/22 ns (r7 Q5), per-channel spread / tagger jitter / SNSPD at Mcps / SPAD tail (r8 Ph6), **a control regime that survives W4 truncation** (r8 Ph4) | r1-r8 | Open | One sitting before any M4b/M5a run |
| O9 | Purge the bench binaries from history | r7 S5 | Removed from tree | Leave (default) |
| O10 | A real divided T3 file | r7 G5 | Open | Source it, or ship the refusal |
| O11 | Withdrawn deferral (pending afterpulse moved into M4) | r5 item 4 | Plan updated | Close |
| O12 | Propagation gate refusal tolerance: replace "any refusal fails" with "at most 5 refused per 1,000 (binomial at 1e-3, FALSE_FAIL), refused counted as not covered". Give the misfit (8e-3 per replicate: at most 18) and flag caps the same binomial form. Check the p<0.01 share too (St5) | r8 St2 | New | Adopt; this is a gate rule, so it is the owner's |
| O13 | Which of M2b / M3b / M4b may slip to 0.2.x | r8 Pl2 | New | M4b gates with M5a; M3b may slip |
| O14 | PR #26 merge method. The reviewer advises a non-squash merge so the round commits stay traceable; the owner's practice for phase 3 has been squash merges (#20-#25). Also: fixes to merged milestones in their own PRs, or a review-round cap | r8 Pl6 | New | Owner's choice; the council does not decide it |
