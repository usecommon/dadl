> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 3 (goodness of fit, guardrail red team, numerics)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: G = goodness-of-fit specialist (r3-gof.md), R = guardrail red team (r3-red.md, its rows 1-6
as R1-R6), N = numerical analyst (r3-num.md). Duplicates are merged under one id. Council probes:
`scratchpad/council/r3-council/` (`proto.py` is the prototype of the merged CHSH fix; outputs `out_*.txt`).

## Checked numerically

All runs on branch m4/series-deadtime at 98bbc25 with the repo venv, through the package's own functions
(the prototype monkeypatches `chsh.capture_balance` and `chsh._delay_drift` only).

- **G1 + G2 + N1 (CHSH capture balance and delay drift): confirmed, and the merged fix verified.** The
  reviewers' numbers stand (rerun of `r3-num/p10_drift.py` before the patch: refusals on clean data with
  30 ps jitter and a 10 ns window; `r3-gof/out_p1b_*`: 2-65% false drift refusals). Three causes, one check:
  (a) the pooled fit bins at W/30 when no jitter is given, so a narrow peak sits in one or two bins (N1);
  (b) it fits only |delay - declared| < W, so a peak as wide as the window trades its offset against the
  background (G2); (c) the drift threshold is an effect size with no noise term (G1).
  Prototype (`proto.py`): fit the pooled peak iteratively on a region half = W + 5 sigma with bins
  max(tick, sigma/3) (at most 4,000 bins), starting from the current coarse fit, polished, up to 3 passes until
  the bin width is stable; data-side sd of the reweighted S by the delta method (offset se sigma/sqrt(A), jitter
  se sigma/sqrt(2A), A the fitted peak area); artifact-side sd as the sd of the reweighted S over the delay
  artifacts' capture draws (n = 200); refuse when |diff| > max(0.5 stat se, 3.29 sqrt(sd_art^2 + sd_data^2));
  raw label when |shift| > max(0.5 raw se, 3.29 sd_data) or a capture < `MIN_CAPTURE`.
  - **False drift refusals on clean data (current artifacts), `drift.py`:** 150 ps / 1 ns, delays declared:
    **0/60** (was 2/60); same, declared 0: **0/60** (was 12/60); 350 ps / 2 ns declared: **0/60** (was 5/60);
    declared 0: **0/60** (was 6/60); 350 ps / 1 ns: 0/20 (was 7/20); 350 ps / 0.8 ns: **0/20** (was 13/20);
    21 ps / 1 ns: 0/20. Narrow peaks (`n1drift.py`, N1's physics, 30 ps, 10 ns, 1e6/s, 4 s): **0/4** (was 3/4),
    offsets 0.0 +- 0.1 ps (were +4 to +17). Null z = diff/sd: mean +0.13 to +0.22, sd 0.88-0.95 (slightly
    conservative). Total 0/264 at the 0.001 design.
  - **Power (refused by the drift check, plus by the predictive check):** 150 ps / 1 ns: 100 ps 20/20,
    200 ps 20/20, 400 ps 19/20 + 1. 350 ps / 2 ns: 200 ps 17/20 + 3, 400 ps 12/20 + 8, declared 0 with 200 ps
    14/20 + 6. 21 ps / 1 ns: 100 and 200 ps 20/20 (z -231, -582). Every 200-400 ps run is refused. The
    predictive-check refusals at 400 ps are G3's stall (z of the drift statistic -23 there), so after Fix 2
    they become drift refusals. 100 ps at 350 ps / 2 ns: 0/20 (mean shift 0.41 stat se, below the 0.5 se floor
    by design; L1). Sub-floor 50 ps at 150 ps / 1 ns: 2/30 (was 11/30).
  - **Raw label on matched delays, `raw.py`:** pair jitter 495 ps (G2's case), W 0.8 ns: **0/20** (was 9/20),
    offsets sd **2.7 ps** (was 51); W 1 and 1.5 ns: 0/20; W 0.6 ns: 20/20, all through `MIN_CAPTURE` alone
    (capture 0.43, the case G5 renames). Narrow peaks: 15 ps / 2 ns 0/10, 30 ps / 2 ns 0/10, 30 ps / 10 ns
    0/6, offsets within 0.2 ps. Misdeclared (A+, B+ 400 ps late, W 2 ns): 20/20 labelled, (raw - rew)/se +3.95.
  - **G4 rechecked after the wider fit:** null drift statistic mean +0.0005 (pooled over 240 runs, about
    2.9 se of the mean; 0.1 stat se, a tenth of the effective threshold). Smaller than G4's +0.0011; harmless.
- **G3 (predictive-check fit stalls): confirmed.** `r3-council/p13_fc.py zero 0 36 43`: worst histogram p
  2.0e-6 and 7.2e-7 unpolished, **0.687 and 0.82** polished (NLL lower by 24 and 26); these are the 2/60
  "check refusals" in `out_drift_350_2000_zero_0_60.txt`.
- **N2 (series fit on the wrong rate branch): confirmed.** `r3-num/p14_branch_nll.py` rerun: at tau_p = 40 ns
  the 6e7 cps run's NLL is 30.3 (seed 0) and 33.4 (seed 1) lower on the high branch (n 6.0e7) than at
  Newton's low-branch result (7.1e6); the low-branch profile is lowest at the shortest gap (45 ns = tau_np),
  so tau_p is biased to tau_np. Not rerun: the 4/4 coverage failures of `p13_branch_ap.py` (consistent).
- **N3 (`survival` on an int tau): confirmed.** `survival(np.array([1.0]), 1e-6, 20000)` raises
  `UFuncTypeError` (cast float64 to int64).
- **R1-R6: confirmed by code.** `data.py:34-35` (`role`, `detector` with `compare=False`); `propagate._check`
  has no identity or gap check. `propagate.condition_mismatches`: `temperature_K` and `_min_/_max_` checked
  only `if ... in dataset`, `bias*`/`tagger_*` refused when absent; `_check`: `if valid:` (no `valid_for`, no
  check, no label). `forward.py:572` `labels = [FORWARD_MODEL, ACCIDENTALS_MODELLED]` whatever `corrections`
  holds (our round-2 Fix 9d wording). `propagate.py` `key = (m.kind, m.channels)` (ordered tuple).
  `cli.py:286` and `:954` `_print_forward(..., {}, "S")`. `cli.py:398` "recorded, not enforced";
  README.md:251 "not yet compared", :264 "recorded, not enforced"; docs/design.md:343 same.
  Not rerun (accepted on the reviewer's CLI evidence): R1's V 0.9219 vs 0.9202, R3's corrected S 2.1486.
- **Method versions.** N2 and N3 edit `calibration/deadtime.py` (pinned, `deadtime` version "2"). **Bump
  to "3"** (`METHOD_VERSION`) and update `tests/method_versions.json` (version and sha256). No other
  pinned module is edited.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| G1 + G2 + N1 (+ G4) | CHSH capture balance and delay drift: narrow peaks binned too coarsely (N1), wide peaks fitted in too small a region (G2), drift refusal with no noise term (G1); false raw labels and false refusals (3-65%) on clean, correctly calibrated data | Y/Y/Y | ADOPT-FIX | Fix 1 (one coherent fix; verified 0/264 false refusals, 200-400 ps drift always refused). Propagation re-run chsh-spad, chsh-dead-time |
| G1 (label part) | Label, not refuse, when the shift exceeds 0.5 se but is within the noise | N/N/Y | REJECT | That noise is already in the corrected interval's calibration part (capture draws); a label would fire on most clean runs with short delay calibrations. Record `delay_drift_z` and `delay_drift_diff` in counts instead (Fix 1) |
| G3 | Predictive check's L-BFGS-B fit stalls when the start delay is tens of ps off; refuses clean data (2/60) | Y/Y/Y | ADOPT-FIX | Fix 2. Propagation re-run all five regimes |
| G4 | Data-side fit pulls the ±400 ps offsets toward the declared delay (drift statistic +0.0011 on the null) | Y/N/N (fix) ; Y/Y/Y (state) | ADOPT-LIMIT | After Fix 1 the bias is +0.0005 (0.1 stat se). L2 in `_delay_drift`'s docstring; the Fix 1 null test tolerates it |
| G5 | `CAPTURE_IMBALANCE` says "declared delays miss the pair peaks" when only the window is narrow against the jitter | Y/Y/Y | ADOPT-FIX | Fix 1 (separate `CAPTURE_LOW` text) |
| R1 | Another detector's dark artifacts accepted on the same tagger input (identity ignores `detector`), getting around the bound refusal; a 1 us dead time applied to a run with 30 ns gaps | Y/Y/Y | ADOPT-FIX | Fix 4 (serial check, gap check). Propagation re-run (all regimes run the gap check) |
| R2 | Bracketing fails open: no `valid_for` = valid forever unlabelled; recorded temperature or bracket not given for the dataset is skipped; result does not record the dataset context | Y/Y/Y | ADOPT-FIX | Fix 5 |
| R3 | Corrected CHSH says "dead time and window capture modelled" when neither was applied | Y/Y/Y | ADOPT-FIX | Fix 3 |
| R4 | One-per-kind-and-channel check depends on channel order; user can pick the artifact by `--created` | Y/Y/Y | ADOPT-FIX | Fix 5 (`frozenset`) |
| R5 | CHSH CLI never prints the D3 range | Y/Y/Y | ADOPT-FIX | Fix 6 |
| R6 | Docs say validity and tagger settings are not enforced; code enforces them | Y/Y/Y | ADOPT-FIX | Fix 7 |
| N2 | Series gap law: Newton stays on the low-rate branch; a scan past the count-rate maximum biases tau_p to tau_np (4/4 intervals exclude the truth) | Y/Y/Y | ADOPT-FIX | Fix 8. **deadtime method version 2 -> 3**; coverage re-run series-spad, passive-quench. A coverage regime with a run past 1/tau_p: DEFER to M10 (regimes are declared before the first run; owner's agreement, next to rotating seeds); a unit test now |
| N3 | `survival`/`gap_cdf` crash on an int tau_p | Y/Y/Y | ADOPT-FIX | Fix 8 (same version bump) |

### Limit wordings

- **L1** (chsh.py `_delay_drift` docstring; README corrected-CHSH limits bullet, after the B1 numbers):
  "The drift check refuses only a drift that moves S by more than half its statistical error and more than
  3.29 standard deviations of the check's own noise (the delay artifacts' capture draws and the analysed
  peak fit). A smaller drift is not refused (100 ps at 350 ps jitter and a 2 ns window: 0.41 se, 0/20
  refused); one within the artifacts' noise is covered by the calibration part of the interval. Longer delay
  calibrations make the check sharper."
- **L2** (chsh.py `_delay_drift` docstring): "On clean data the statistic has a small positive mean (+0.0005,
  about 0.1 statistical se, from the peak fit at offsets of +-400 ps); it is well inside the threshold."
- **G5 label** (chsh.py, new `CAPTURE_LOW`): "window narrow against the pair jitter: capture {c:.2f} on
  {pair} with the peak on the declared delay (offsets in counts); S is post-selected more strongly on the
  timing". `CAPTURE_IMBALANCE` keeps its text and fires only on the shift rule or on a low capture whose
  fitted offset exceeds its sigma/2.
- **R3 label** (forward.py, replacing the constant `ACCIDENTALS_MODELLED` text, built from `corrections`):
  "{modelled} modelled; {not} not modelled; detector efficiencies not (an output imbalance biases corrected S
  as it does raw S); secondary, not a Bell claim (the headline is raw S)", e.g. "accidentals modelled; dead
  time (no dark artifacts) and window capture (no delay-jitter artifacts) not modelled; ...".
  `DEAD_TIME_UNCORRECTED` (added when the raw value carries `DEAD_TIME_SELECTION` and no dark artifact is
  used): "dead time not modelled: the dead-time selection bias of the raw value stays in the corrected value".
- **R2 labels** (propagate.py): `NO_VALIDITY` "calibration {aid} has no validity range: not checked against
  the acquisition time"; condition mismatch text for a missing temperature or bracket: "{key} = {v:g}
  recorded, not given for the dataset" (same as `bias*`).
- **R6 docs**: cli.py:398 "--valid-until: end of the intended validity (ISO 8601); an analysis acquired
  outside it is refused unless calibration_outside_validity is allowed"; README.md:249-251 and docs/design.md:343
  "tagger settings (`tagger_*`) and bias are compared when the calibration is used: a mismatch, or a setting
  recorded but not given for the dataset, is refused (overridable)"; README.md:263-264 "validity is enforced
  (from acquisition, or creation without one, to `--valid-until`); an artifact with no validity range is
  labelled; supersession is recorded, not enforced".

## Implementation order (adopted fixes)

**Pinned method version: `deadtime` "2" -> "3"** (Fix 8; `calibration/deadtime.py` `METHOD_VERSION` and
`tests/method_versions.json` version and sha256). No other pinned module changes.
**Gates:** coverage, deadtime regimes series-spad and passive-quench (all shards, CI seed), after Fix 8.
Propagation, all five regimes (hom-snspd, hom-bright, hom-crosstalk, chsh-spad, chsh-dead-time) at 1,000
replicates on the CI seed, once after Fixes 1-5: Fix 2 changes predictive-check p-values everywhere, Fix 4's
gap check runs in every corrected path, Fix 1 changes the CHSH flag sources. Expect coverage unchanged, CHSH
drift flags 0-0.3%, flag shares in clean regimes <= 1%, hom-crosstalk still >= 90%.

1. **CHSH capture balance and drift (G1, G2, N1, G5; prototype `r3-council/proto.py`).**
   (a) `forward_check.pooled_peak(pairs) -> PeakFit | None`, `PeakFit(offset_ps, sigma_ps, offset_se_ps,
   sigma_se_ps)` (update `capture_balance`; nothing else calls it). Keep the current coarse fit as pass 0 (and
   its `MIN_COUNTS` check over ±W); then up to 3 passes: sigma clipped to [tick/2, 4W], half = ceil((W + 5 sigma)
   / tick) ticks, width = max(tick, floor(sigma/3), ceil(2 half / 4000)) in ticks, `fit_peak(..., polish=True)`
   from the previous (d, sigma); stop when the width does not change. se: sigma / sqrt(A) and sigma / sqrt(2A),
   A = exp(theta[1]). Docstring: the region extends 5 sigma beyond the window so the background is fitted
   outside the peak; bins resolve the peak.
   (b) `chsh.py`: `_setting_outcome_counts(runs, window)` computed once; `reweighted_s` takes those counts
   (no recount per call). `capture_balance` returns per pair `(offset, sigma, capture, offset_se, sigma_se)` and
   `sd_data` (delta method: shift each pair's offset by its se, then its sigma, recompute S, sum squares).
   Counts: `capture_offsets_ps`, `capture_reweighted_s`, new `capture_sd_data`.
   (c) Raw label: `CAPTURE_IMBALANCE` when `|shift| > max(0.5 raw.stderr, 3.29 sd_data)`, or a capture <
   `MIN_CAPTURE` with |offset| > sigma/2; `CAPTURE_LOW` (text above) when a capture < `MIN_CAPTURE` otherwise.
   Module docstring updated (both conditions; 3.29 = two-sided 0.001).
   (d) `_delay_drift`: `_pair_capture(..., n=400, seed)`; `sd_art` = sd over draws of S reweighted by each
   draw's four captures; `sd = hypot(sd_art, sd_data)`; refuse (or label under the override, as now) when
   `|diff| > max(0.5 stat se, 3.29 sd)`; message adds "z = {diff/sd:.1f}". Counts: `delay_drift_diff`,
   `delay_drift_z`. L1, L2 in the docstring.
   Tests (`test_chsh_corrected.py`): `test_matched_delays_are_not_labelled` gains a 30 ps pair jitter with 2 and
   10 ns windows and a 495 ps jitter at 0.8 ns (no label; offsets within 5 ps); W 0.6 ns at 495 ps gives
   `CAPTURE_LOW`, not `CAPTURE_IMBALANCE`; a slow-marked null test: 12 seeds of the `r3-council/drift.py`
   physics (150 ps, 1 ns, declared 0, current artifacts) with no drift refusal; the existing 200 ps drift test
   still refuses, and 100 ps at 150 ps / 1 ns refuses.
2. **Predictive-check stall (G3).** `forward_check.pair_pvalue`: after the L-BFGS-B starts, when the grouped
   p < 0.01 refit with `polish=True` (and add a start at the histogram's highest bin) and keep the better NLL;
   polishing only lowers the NLL, so the check stays the same test with a better optimiser and costs little.
   Remove "the predictive check tolerates that" from `fit_peak`'s docstring. Test: `r3-council/p13_fc.py`
   seeds 36 and 43 (zero) give p > 0.1 via `pair_pvalue` (fixture by seed, slow-marked).
3. **Corrected CHSH labels (R3).** `forward.corrected_chsh` (and the HOM/fringe forward paths if they carry the
   same constant): build the label from `corrections` (R3 wording); add `DEAD_TIME_UNCORRECTED` in
   `chsh._corrected` when the raw labels hold `DEAD_TIME_SELECTION` and "dead_time" is not in corrections.
   Keep `ACCIDENTALS_MODELLED` exported (now a format string). Test: `r3-red/chsh/p_dts.py` physics with an
   empty calibration folder: label says "dead time ... not modelled" and carries `DEAD_TIME_UNCORRECTED`.
4. **Wrong detector's artifacts (R1).** (a) `propagate._check`: when the artifact's channel record and the
   setup's `ChannelId` (passed through `CalibrationSet.load(..., channels=...)`, which the CLI and the
   analysis paths already know) both give `detector` and they differ, refuse with `StaleCalibration`
   ("calibrated on detector X9, analysing A1"), overridable by `calibration_outside_validity` (then
   `OUTSIDE_VALIDITY`). `ChannelId` equality is unchanged (L-note in its docstring: identity is the input; the
   detector is checked at use). (b) New `forward.check_dead_time_gaps(calibrations, channels, tags)`, called
   from the corrected HOM, fringe and CHSH paths with the analysed tags: for each used dark artifact whose dead
   time is identified (not a bound), count the analysed gaps shorter than 0.9 x its 95% lower end; when more
   than max(10, 1e-3 x gaps) refuse with `StaleCalibration` ("dead time {tau} ns of {ch}: {k} analysed gaps
   are shorter, down to {min} ns"), overridable as above; record `min_gap_ps` per channel in counts. Bounds are
   not checked (gaps below an upper limit are consistent with it). Tests: `r3-red` hw_new/hw_old physics
   (serial mismatch refused; with the override labelled); a 1 us identified dark artifact on a 30 ns-gap run
   refused; the gate regimes' clean replicates pass (checked by the propagation re-run).
5. **Bracketing and ambiguity (R2, R4).** `condition_mismatches`: a recorded `temperature_K`, or a
   `_min_/_max_` bracket, whose dataset value is not given is a mismatch "recorded, not given for the
   dataset" (as `bias*`). `_check`: no `valid_for` adds `NO_VALIDITY` (wording above). `LabelledValue.counts`
   gains `dataset` {acquired, conditions} when a DatasetContext was given (CLI `_calibration_set` passes it
   through). Ambiguity key `(m.kind, frozenset(m.channels))`. Tests: the `r3-red/mix` set with no
   `--condition` refuses (temperature recorded); no `valid_for` labelled; `r3-red/amb` raises `ConfigError`.
6. **CLI D3 range for CHSH (R5).** cli.py:286 and :954: `_print_forward(result, result.counts, "S")`. Test:
   CLI on a bounded-dark CHSH set prints "at zero ... at the bound".
7. **Docs (R6, limits).** R6 wordings; L1 in README; CHANGELOG entries for Fixes 1-8 and the deadtime method
   version 3; docs/phase3-plan.md smaller-items row: "Coverage regime with a dead-time scan past the count-rate
   maximum (1/tau_p) | M4 council round 3 (N2) | M10, next to rotating seeds (owner's agreement)".
8. **Series rate branch and int tau (N2, N3; deadtime method version 3).** `deadtime._Fit.best_n`: when tau
   > 0, run Newton also from the other branch's rate, n' = -W_k(-n tau e^{-n tau}) / tau with k = -1 if
   n tau < 1 else 0 (`scipy.special.lambertw`, real part), and keep the lower NLL; cache the winning rate per
   run (warm start keeps its branch). When afterpulsing is zero or unidentified and the two NLLs agree within
   0.5, the run's note says "input rate {n1:.4g} or {n2:.4g} cps (fitted; past the count-rate maximum the
   gaps do not tell which)". `survival`, `gap_cdf`: `tau = float(tau)` (and `n = float(n)`) on entry.
   `METHOD_VERSION = "3"`; `tests/method_versions.json`. Tests: `r3-num/p13_branch_ap.py`'s scan (1e6, 3e6,
   6e7 cps, tau_p 40 ns, 2% afterpulsing), 2 seeds: interval covers 40 ns; `survival` with an int tau equals
   the float result. Then the coverage re-run (series-spad, passive-quench).
