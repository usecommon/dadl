> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Phase 3 statistical design as a whole (multiplicity, joint false-refusal rate, gate coherence, T3 tick)

HEAD 8d7c50a, repo venv, one BLAS thread. Probes and raw outputs are in `scratchpad/council/r8-stat/`:
- `joint.py`: the gate's own calibration and analysis path (`propagation._calibrate`/`_replicate` logic),
  with every check made non-fatal by a permissive policy. Each check's statistic is recorded, and the
  default-policy verdict is recomputed from that statistic. Outputs are in `out/*.jsonl`, summarised by
  `summarize.py`.
- Runs on clean data: chsh-spad 400 at seed 4242 and 400 at the CI seed 20261118; chsh-dead-time 400,
  hom-bright 400 and hom-snspd 600 (example duration 1 s) at seed 4242.
- `darktail.py`: 12,000 chsh-dead-time dark artifacts, in `dt/`.
- `delayppc.py`: 6,000 clean delay calibrations, in `dl/`.
- T3 probes: `t3check.py`, `t3check2.py`, `t3delay.py` and the end-to-end `e2e/`.
- Gate reproduction: `rep69.py`, `rep154.py`, and shard 3/16 of seed 4242.

## Joint false-refusal rate on clean data (default policy with the corrections allowed, no overrides)

| check (refuses unless overridden) | nominal | chsh-spad (example) 800 | chsh-dead-time 400 | hom-snspd (example) 600 | hom-bright 400 |
|---|---|---|---|---|---|
| calibration misfit, any artifact (write refused by `calibrate`; load refused without `misfit_calibration`) | 1e-3 per artifact | 6/800 (0.75%, 95% 0.28-1.6%) | 2/400 | 1/600 (per-artifact 2.5e-3 from 4,000, row 5) | 1/400 |
| forward-model predictive check (`Misspecified`) | 1e-3 per analysis | 0; adj. p<0.01 in 0.4% | 0; 0.5% | 0; 0.33% | 0; 1.0% |
| CHSH delay drift (`Misspecified`) | "3.29 sd" | 1/800 observed; **5.1e-4 expected** from the recorded z (threshold median 3.49 sd) | 0 (threshold 5.9 sd) | n/a | n/a |
| dark dead time contradicted by gaps, downward (`StaleCalibration`) | "1e-4 per channel" | n/a | **2/400**; 2.5e-4 per channel from 12,000 artifacts (row 1) | n/a | 0 |
| upward gap test, dead fraction > 20%, loose bound | ~1e-4 / 0 / 0 | n/a | 0 (excess never > 0) | n/a | 0 |
| **any refusal per clean analysis** | | **~0.8%** (analysis time alone ~5e-4) | **~1-1.5%** (analysis time alone ~1e-3) | ~0.3% | ~0.3-0.5% |
| labels firing on clean data | | INTERRUPTED 2/800 | DEAD_TIME_BELOW_DARK on 62/400 (15.5%) of analyses | none beyond fixed text | none |

The 1-in-10^4 budget is not met jointly. It was stated per channel, for the downward gap test only (rounds 4
and 6), and even there it is missed (row 1). Every other refusing check runs at 1e-3 per artifact or per
analysis, and a CHSH analysis multiplies them: 4-8 artifacts, 4 channels, plus the forward and drift checks.
A clean CHSH workflow therefore meets at least one refusal about once in 70-120 tries, and the analysis-time
refusals alone run at 5e-4 to 1e-3. No document states a per-analysis budget.

## Findings

| # | Severity | Finding | Evidence | Where | Suggested fix |
|---|---|---|---|---|---|
| 1 | wrong | The downward gap check refuses clean dark artifacts at 2.5e-4 to 3.7e-4 per channel, about 1-1.5e-3 per four-channel analysis. The documented rate is "about 1 in 10^4". The 3-width floor assumes the dead time's lower tail is exponential with the interval's own width as its scale. That fails when the dark fit falsely detects afterpulsing on a detector with none: the interval turns two-sided and narrow (B-: [124.6, 128.5] ns for a 50 ns truth). It also fails on plain shortest-gap tail events (263 ns point, [192, 263] ns). | `darktail.py`: floor above the truth in 3/12,000 artifacts (95% 5e-5-7.3e-4); 2 of the 3 had afterpulsing reported, which happens in 7/4,441 identified artifacts. `joint.py` chsh-dead-time seed 4242: 2/400 analyses refused, replicates 69 and 360 ("dead time 74.04 ns of B- ... 11884 analysed gaps are shorter"; "dead time 128 ns ... 60440"). Pooled 5/13,600 per channel; P(>=5 at 1e-4) = 0.013. | forward.py:272-276, 359-377; README:419-425 | Take the floor from the pivot the dark data actually give, not from the fitted interval's width. The shortest dark gap minus tau is Exp(N f) for N dark gaps at rate f, and afterpulsing only raises the edge density, so this is conservative: floor = point - ln(K/alpha)/(N f), with K the analysed channels and alpha the per-analysis budget. At 0.1 s and 20 kcps this is point - 265 ns, so the test never fires there, which is correct: such a dark run cannot exclude anything at 1e-4. At 1 s and 50 kcps it is point - 4.2 ns, close to today's "7% below". Re-run chsh-dead-time and hom-bright. |
| 2 | wrong | The propagation gate has no tolerance for a refusal: any refused clean replicate raises out of `_replicate` and aborts the gate. The POLICY allows the misfit overrides but not `calibration_outside_validity`, and the dead-fraction and loose-bound refusals are not overridable. With finding 1's ~1e-3 per replicate, a 1,000-replicate chsh-dead-time gate aborts on about 63% of seeds. The CI seed's pass is luck (about 37%), not evidence. The gate also never measures the default-policy refusal rate: forward and drift misfits are allowed and only counted, against a 1% cap (10x alpha). | `python -m dadl.analysis.propagation chsh --regimes chsh-dead-time --replicates 400 --seed 4242 --shard 3/16`: traceback `StaleCalibration ... 74.04 ns of B-`, exit 1. `rep69.py` 4242 360: same through `_replicate`. | analysis/propagation.py:212-237, 63, 166-176 | Catch `DadlError` per replicate and count it under a key "refused", split by type. Gate that share against a declared per-analysis budget with a binomial bound (FALSE_FAIL), and decide whether refused replicates count as not covered or leave the denominator (the former is conservative). Give the gate's misfit and flag caps the same binomial treatment. Run one regime on a second, rotating seed. |
| 3 | wrong | Round-7 Fix 1(d) breaks T3 delay calibration at a measured sync period. `calibrate pairs` passes `tag_tick_ps`, which for a non-whole period is the gap-law bound res + 1 = 26 ps. The delay fit then bins a 25 ps lattice in 26 ps bins: every clean artifact is a misfit, and the jitter and delay are biased outside their intervals. The note blames "an offset between channels?". | e2e, identical photons, TimeHarp 260 P T3 at 25 ps, header period 12500.0 vs 12500.3 ps: ppc p 0.571 vs **9.9e-55 (no artifact written)**; jitter 141.6 [138.8, 144.4] vs **137.0 [134.3, 139.7]**; delay 2.1 vs **-10.3 [-14.1, -6.4]** (truths 141.4 and 0). `t3delay.py`, 10 seeds each: tick 26 gives 10/10 misfit and 0/10 jitter covered; tick 25 gives 0/10 misfit and 10/10 covered; no tick (1 ps from the tags) gives 10/10 misfit. | cli.py:812; adapters/ptu.py:416-426; tagfile.py:302-318 | Keep two quantities apart. The lattice step (the exact dtime resolution, 25) is used for binning and the quantization model. The error bound (resolution + 1 when rounded) is used only for the gap-law widening. For example, `PtuFile.lattice_ps` = res, with `tick_ps` keeping its round-7 meaning for deadtime, or deadtime adds the +1 itself from `exact_ps`. Pass the lattice step to `calibrate_delay`, `pair_pvalue` and `pooled_peak`. Test: the e2e pair above gives p > 1e-3 and jitter covered in >= 19/20 seeds. |
| 4 | questionable | Labels and refusals do not follow one error-rate standard, and there is no joint budget. Round 6 declared 1e-4 "the standard for a downward refusal", but the forward check (1e-3), the drift refusal (5e-4 per chsh-spad analysis), calibration misfits (1e-3 per artifact, x4-8) and the INTERRUPTED label (1e-3) all sit 5-10x looser. The drift refusal tells the user "a drift since calibration", but on clean data it fires on a 7 ps apparent shift at 350 ps jitter. Its statistic is well calibrated; the 0.5 se consequence floor binds in 96-98% of chsh-spad runs and sits at only 3.3-3.5 sd. | Joint table above. Drift z: mean -0.03 to +0.03, sd 1.00-1.02, KS p 0.45-0.75 (1,200 runs). `rep154.py` chsh-spad 4242 154: `Misspecified ... A-/B+: calibrated 5, analysed -2 ps; ... z = 3.8`. | analysis/chsh.py:73, 521; forward_check.py:41; calibration/base.py:19 | State one per-analysis false-refusal budget in the README (for example 1e-3 at analysis time) and split it across checks (alpha spending). Forward check 4e-4, drift 2e-4 (for example, refuse at max(0.5 se, 3.7 sd) or at 0.5 se plus 3.29 sd), gap tests 2e-4 shared over the channels (finding 1), interruption 2e-4. Document calibration-misfit refusals separately as a per-artifact rate with its CHSH multiplicity (about 0.5-1% of four-pair calibration sets). |
| 5 | questionable | The two predictive checks of the same Gaussian peak use opposite degrees-of-freedom rules. `delay.py` uses groups - 5, which is anti-conservative between the Chernoff-Lehmann bounds; `forward_check` uses groups - 1, and its own docstring records that groups - params gave about 3x the nominal rate. In the SNSPD example regime the delay misfit rate is 2.5x nominal, and phase 3 multiplies it by four delay artifacts per CHSH. | `delayppc.py`, hom-snspd pair runs (0.1 s), 4,000 artifacts: p<0.05 6.3% (binomial p 1e-4), p<0.01 1.2%, **p<0.001 10/4,000 = 2.5e-3** (p 0.008). chsh-spad, 2,000 artifacts: 5.1%, 0.65%, 0.05% (nominal). Same direction in `joint.py`: hom-snspd p<0.05 7.8% of 600. | calibration/delay.py:168, 403 (unchanged since v0.1.1; feeds phase 3) | Use groups - 1 in `delay.py` as `forward_check` does (method version bump; re-run the delay coverage gate). Have the gates check the p<0.01 share, not only p<0.001, so a check that runs 2-3x hot is visible at 1,000 replicates. |
| 6 | minor | The upward gap test's stated "chance about 1e-4 per channel" is the wrong number in both directions. Against a matched upper end the exact tail is (1 + 9/100)^-100 = 1.8e-4. In practice the upper end is the shortest dark gap, which is never below the truth, so on clean data the excess is never positive and the rate is far below 1e-4. | 2,400 checked channels (chsh-dead-time, hom-bright): excess > 0 in 0. | forward.py:276, README:424 | Write "at most about 2e-4 per channel (1.09^-100), far less in practice because the upper end is the shortest dark gap". |
| 7 | minor | The gate's "forward misfit" key merges two checks: drift refusals carry the same `FORWARD_MISFIT_ACCEPTED` label, so the 1% cap cannot tell them apart and the CHANGELOG reads it as the predictive check's rate. | `joint.py` chsh-spad 4242: the one labelled replicate (154) is a drift event, with forward-check adjusted p >= 0.001. | analysis/chsh.py:523; analysis/propagation.py:237 | Count drift and predictive-check events under separate keys, or give the drift its own label. |

## The open item from round 7: `pair_pvalue` on T3 tags at a 1 ps tick (strict xfail)

Diagnosis. Within one sync period a T3 time is `n P + r floor(phi / r)`, so the recorded difference of two
tags lies on the lattice r·Z, plus a constant coset offset (kP mod r) when the pair spans k syncs, plus up to
1 ps from rounding. For CW light the two phases are independent and uniform, so the recorded difference is
the true one convolved with a triangle of half-width r and evaluated on the lattice. `_ShapeModel`'s
box-shifted integral with `tick` = r is exactly that model, provided the bins are whole multiples of r and
aligned to the lattice. Today the tick is read from the tags. For a non-whole period that gives 1 ps, so the
bins are int(sigma/3) = 46 ps wide over a 25 ps lattice: each bin holds 1 or 2 teeth, and p is about 1e-100.
Passing the gap-law bound of 26 ps makes it worse (p about 1e-150). Measured in `t3check.py`, 3 seeds each:

| res, period | tick from tags (today) | tick = res | tick = res + 1 (round-7 `tick_ps`) |
|---|---|---|---|
| 25, 12500.3 | median p 6e-98 | **0.42** | 1e-144 |
| 25, 12500.0 (exact) | 0.23 | 0.23 | 2e-156 |
| 16, 12345.6 | 4e-28 | **0.44** | 2e-54 |
| 4, 12500.3 | 3e-4 | **0.46** | 1.5e-3 |

Null check with tick = res, 40 seeds each (`t3check2.py`): p<0.05 in 3%, 3%, 3% and 3% of runs for
(25, 12500.3, no cable delay), (25, 12500.3, 30 ns cable delay), (16, 12345.6, one-period cable delay) and
(16, 12345.6, no cable delay). That is nominal, including pairs that span sync periods.

Proposed fix (statistically exact for CW, no loosening):
1. Carry the lattice step separately from the gap-law bound (finding 3). `PairData` gets `tick_ps`, filled by
   the analysis path from the TagSet metadata of T3 input (as `t3_input` already reads it). `pair_pvalue`
   and `pooled_peak` use `tick = tick_ps` when given, else `tag_tick(a, b)`.
2. Align the bin edges to the lattice phase: take the modal residue of (b - a - shift) mod tick inside the
   window, and put the edges just below it. This absorbs input offsets that are not multiples of r, and the
   kP mod r coset of a cable delay longer than the period. Rounding adds +0 or +1 ps, which stays inside a
   bin.
3. Limit to state: when the window's pairs split between two sync spans (|declared delay mod P| comparable
   to sigma), the histogram mixes two cosets offset by (P mod r). The exact model is the two-coset sum
   weighted by the spanning share. It is negligible for sigma >> r, and can be checked as a sigma ~ r case
   in the test.
4. Turn the xfail into a pass at tick 25, and add the 16 ps / 12345.6 ps case and a cable delay that spans
   syncs. Apply the same tick to `calibrate_delay` through the CLI (finding 3) and to `capture_balance`'s
   `pooled_peak`. The CHSH drift statistic on T3 input is unmeasured until then.

## Checked and found sound

- The delay-drift statistic is well calibrated on clean data: z ~ N(0,1) over 1,200 CHSH runs in 3
  regimes/seeds. The capture-balance shift divided by its sd_data has sd 0.99-1.05.
- The forward check is conservative, as designed: adjusted p<0.01 in 0.25-1.0% of analyses against a
  nominal 1%. No clean analysis is flagged in 2,000.
- `DEAD_TIME_BELOW_DARK` is an accurate label. With about 1.2M analysed gaps it fires exactly when the dark
  lower end misses (62/400 chsh-dead-time analyses, as 5% per channel compounds over 4 channels). Coverage
  given the label is 93.5% (n = 62), against 97.0% without it.
- Coverage on a seed other than CI's (4242) sits in the band in every regime: chsh-spad 93.5%
  (CI seed 94.0%; pooled 93.75% of 800, p = 0.10 against 95%, worth watching), chsh-dead-time 96.5% with
  control 86.0% below the band, hom-bright 95.0%, hom-snspd 95.3%. The coverage bands themselves show no
  sign of having been tuned to the CI seed. The seed dependence found is finding 2, which is about
  refusals.
- INTERRUPTED fires on 2/800 clean chsh-spad analyses, consistent with its alpha of 1e-3. No analysis-time
  refusal in 1,000 clean HOM analyses.
- The one-sided rules are coherent: the control is required below the band's lower edge, and `at_least`
  keeps the two-sided lower edge (false-fail rate below the budget).
