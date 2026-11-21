> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 2 (corrected values: CHSH, HOM/fringe, propagation)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: B = Bell-test experimentalist (r2-bell.md), H = HOM/fringe experimentalist (r2-hom.md),
R = propagation statistician (r2-prop.md). Duplicates are merged under one id. Council probes:
`scratchpad/council/r2-council/`.

## Checked numerically

All runs on branch m4/series-deadtime at 666df7b with the repo venv, through the package's own functions.

- **B1 (per-channel delay errors raise S, unflagged): confirmed.** `r2-bell/p9_misdeclared.py` rerun:
  A+ and B+ 400 ps late, declared 0: raw S 2.1994 +- 0.0087 (truth 1.9799), labels fair-sampling and static
  only; 200 ps: 2.0500. `p10_stale.py` rerun (delay artifacts calibrated before a 100-400 ps drift, strict
  policy): 16/16 pass the predictive check; corrected S 1.985-2.003 (100 ps), 2.033-2.051 (200 ps), 2.181-2.207
  (400 ps). Council prototype of the fix (`r2-council/capbal.py`): fit the peak of each of the four channel
  pairs with `forward_check._ShapeModel` on the histogram pooled over settings, compute
  `forward.grid_capture(offset, sigma, W)`, and recompute S with each combination's counts divided by its
  capture. Shift 0/100/200/400 ps: (raw - reweighted)/se = +0.08 / +1.99 / +7.47 / +23.7, reweighted S
  1.983 / 1.987 / 1.983 / 1.995 (truth 1.980); fitted offsets within 2 ps of the injected ones. Clean chsh-spad
  physics (350 ps jitter, 2 ns window, afterpulsing), 10 seeds: |raw - reweighted| <= 0.15 se. A 0.5 se rule
  separates them.
- **B2 + R2 (CHSH at a dead-time bound: no D3 range) and R1 (two-point: point at zero, interval at the
  bound): confirmed.** `r2-prop/p1_bound.py`: unequal singles, `0.9009 +- 0.0007 (95%: [0.9045, 0.9072])`,
  point outside its interval; at_zero 0.9009, at_bound 0.9059. Code: `forward._dead_times` returns
  `(np.full(n, hi), lo=0, hi)` and the caller unpacks the second element as the point.
  `p8_threshold.py`: raw 0.4975, corrected point 0.4975, interval [0.5196, 0.5245], no "threshold crossed only
  after correction". `p6_chsh_bound.py`: S at the bound 2.6379 [2.6310, 2.6445] with the label, no dead time
  2.6453 (outside), and no `corrected_at_zero_dead_time` in counts (`grep at_zero` finds only _twopoint.py and
  cli.py). `r2-bell/p3_unid.py`: 34/60 dark artifacts of the chsh-dead-time regime are bounds (median 60 ns,
  50-80 ns, truth 50 ns). Council `r2-council/d3_union.py`, 40 replicates of chsh-dead-time on the CI seed: **every
  replicate (40/40) corrects at least one channel at its bound**; interval at the bound covers 37/40, at zero
  7/40, the union of the two 38/40 and is 1.74x wider (median). The union interval would likely over-cover
  the 92.6-97.1% band in this regime, so it is not adopted now (see ruling).
- **R3 (stand-in jitter crashes corrections): confirmed.** `p2c_held.py`: a real `calibrate_delay` on 80 ps
  ticks, 3 ps jitter, gives status ok with `stand_in ('jitter_sigma_ps',)`; `corrected_visibility` then raises
  `KeyError ('delay-jitter', 'jitter_sigma_ps', ...)`.
- **H1 (fringe drift moves V across 1/sqrt(2), unlabelled): confirmed.** `p5_fringe_drift.py`: raw V 0.6989 /
  0.7070 / 0.7147 at 0 / +3% / +6% pair-rate drift on the max run (6 seeds each, all six runs at +6% between
  0.7128 and 0.7165); singles ratio max/min 1.0003 / 1.0317 / 1.0593. Labels: fair sampling and single-basis only.
- **H2 (no reference delay with calibrations drops the capture label): confirmed.** `p8_undeclared.py`: V
  0.8661 +- 0.0051 (truth 0.92) both ways; label present when declared (100), absent when None.
- **H3 (predictive-check fit stalls from the floor background): confirmed for the reviewer's seed.**
  `p2c.py` (seed 2007): package p 6.2e-206 for the reference histogram; its three L-BFGS-B starts stop after
  10-13 iterations ("ABNORMAL" twice) with log background -7.09 against 0.16 and an NLL about 400 units
  higher; five starts with the background from the outer bins all reach p 0.403. `p2d.py` on 40 other seeds
  flagged nothing (the reviewer saw 2 in 100), so it is rare but real. chsh-spad (pair jitter ~495 ps, W 2 ns)
  has no "far" bins at all, hom-snspd about two, so both gate regimes start from the floor or near it.
- **R6 (unused or mislabelled artifacts): confirmed.** `p9_unused.py`: a deadtime-series artifact is cited (3
  IDs). `p7_boundlabel.py`: two dark artifacts with identified dead time and stand-in `afterpulse_prob` give a
  corrected HOM labelled "parameter known only as an upper limit" with no D3 range.
- **R4, R5, R7, H4, H6: confirmed by code.** `propagate._check`: `if dataset is None: return labels`;
  `artifact.py:253` `{"from": manifest["created"], ...}` though `manifest["acquired"]` exists when given;
  `grep -i afterpuls` over forward.py, _twopoint.py, hom/chsh/fringe.py, forward_check.py: nothing;
  `corrected_visibility` has no above-1 check (README.md:197-199 says refused); `cli.py:121-126` prints the
  "V reads low" warning whatever `result.corrections` holds.
- **Gate bands.** `coverage.gate(1000, 1)` = (0.926, 0.971).
- **Method versions.** No adopted fix edits a module in `python/dadl/calibration/` (the pinned set:
  accidentals, dark, dcr, deadtime, delay, heralding). No version bump; the coverage gate is not affected.
- Not rerun (accepted on the reviewers' evidence and code reading): B3 (point offset; matches r1-gate's 0.37
  sd), B4 (efficiencies; the raw-S figure 2.033 is also in chsh.py's own comment), B5 (IRF tail refusals), H7,
  R8.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| B1 | Delay errors on A+ and B+ (misdeclared, or drift after a delay calibration) raise raw and corrected S past 2 (2.20 at 400 ps for 1.98) with no label; the predictive check refits the delay freely | Y/Y/Y | ADOPT-FIX | Fix 3 (raw: capture-balance label) and Fix 5 (calibrated: refuse when the artifacts' capture disagrees with the analysed peaks). Propagation re-run chsh-spad, chsh-dead-time |
| B2 + R2 | CHSH at a dead-time bound: no value at zero, label promises one; D3 not implemented for CHSH | Y/Y/Y | ADOPT-FIX | Fix 1 (range in counts, both ends in the threshold labels, wide-range label). Union interval: DEFER to M10 (owner's agreement; chsh-dead-time would over-cover) |
| R1 | Two-point at a bound: point at zero, interval at the bound; point can sit outside its interval; threshold crossing at the bound unlabelled | Y/Y/Y | ADOPT-FIX | Fix 1 (point at the bound) |
| R3 | Stand-in jitter (coarse tagger, legitimate ok artifact) crashes every corrected analysis with KeyError | Y/Y/Y | ADOPT-FIX | Fix 2 |
| H1 | Fringe: a pump drift between max and min moves V across 1/sqrt(2), unlabelled; singles show it | Y/Y/Y (labels); N/N/Y (refuse in the Bell rule) | ADOPT-FIX | Fix 7 (labels only; refusal would change the default policy) |
| H2 | `hom_from_tags(calibrations=...)` without `reference_delay_ps` treats the reference as unsplit and drops the capture label (V 10 sigma low, unlabelled) | Y/Y/Y | ADOPT-FIX | Fix 6 (ValueError) |
| H3 | Predictive-check fit stalls from the floor background: p ~1e-206 on clean-enough data, refused as Misspecified | Y/Y/Y (extra start); N/N/N ("not checked" on non-convergence: could hide misfits) | ADOPT-FIX | Fix 4. Propagation re-run all five regimes (p-values can change wherever the floor start was used) |
| B3 | Corrected point (shortest-gap dead time) sits about 0.4 sd low in its own interval; CHANGELOG gives no sign | Y/N/N (point from draws); Y/Y/Y (state it) | ADOPT-LIMIT | Wording L1 and the CHANGELOG sign. The point stays the artifact's reported value (consistent with the dark calibration and with Fix 1) |
| B4 | Unequal efficiencies: corrected S labelled "detector effects modelled" keeps the whole imbalance bias (z +13.5, 0/6 covered) | Y/Y/Y (label); Y/N/N (normalisation now) | ADOPT-FIX (label) + DEFER | Fix 9 (label wording), L2. Opt-in per-detector efficiency normalisation: M10 smaller items (0.2.x if not done; new correction name, policy line and a gate regime) |
| B5 | Predictive check refuses Bell runs with mild common IRF tails that do not bias S; only a blanket override | Y/Y/Y (limit); Y/N/N (common-vs-confined test now) | ADOPT-LIMIT + DEFER | L3 (README). The common-vs-confined misfit test rides with the IRF-tail PR (D8) |
| B6 + H5 | README stale: "Not yet: ... corrected CHSH", "statistical only", crosstalk "no analysis flags it", `--calibrations` HOM only, python/README quickstart line | Y/Y/Y | ADOPT-FIX | Fix 11 |
| H4 | Forward-model corrected V above 1 reported unlabelled; README says refused | Y/Y/Y | ADOPT-FIX | Fix 8 (label; keep the unclipped value) and README wording |
| H6 | CLI "V reads low" warning printed even when window capture is modelled | N/Y/Y | ADOPT-FIX | Fix 9 |
| H7 | TWO_POINT does not mention a reference delay inside the coherence time | Y/Y/Y | ADOPT-FIX | Fix 9 (label wording; warning when setup.yaml gives `coherence_time_ps`) |
| R4 | Without a DatasetContext the use-time refusals are skipped silently; slices of the analysed run pass disjointness | N/Y/Y | ADOPT-FIX + ADOPT-LIMIT | Fix 9 (label), L4 |
| R5 | Validity starts at write time, not acquisition; a before/after pair cannot be combined | Y/Y/Y (start); Y/N/N (pair now) | ADOPT-FIX + ADOPT-LIMIT | Fix 10 (`from` = acquired when recorded), L5. A drift-aware pair is the drift filter: stays "Out of phase 3" |
| R6 | Unused deadtime-series artifacts cited (and can refuse an analysis); BOUND_ONLY set by any stand-in column (afterpulse_prob on every quiet dark) | Y/Y/Y | ADOPT-FIX | Fix 9 |
| R7 | M2's required "afterpulsing not propagated" label missing | Y/Y/Y | ADOPT-FIX | Fix 9 |
| R8 | Gate stores 400 draws, users get 2,000; Monte Carlo error unstated | N/Y/Y (limit); N/N/N (gate at 2,000: measured no difference) | ADOPT-LIMIT | L6 |
| G (from B2's evidence) | chsh-dead-time corrects at least one channel at its bound in 40/40 replicates; the control regime's "calibration uncertainty" is partly a fixed bound | Y/Y/Y (state); Y/N/N (lengthen the dark run now) | ADOPT-LIMIT + DEFER | L7. Regime change (dark run long enough to identify the dead time, or a regime gating the D3 range) to M10 next to "Rotating seeds" (regimes were declared before the first run; owner's agreement) |

### Limit wordings

- **L1** (forward.py docstrings of `_dead_times` and `_channel_dead_times`; CHANGELOG M3 bullet): "The corrected
  point uses each dark artifact's reported dead time, the shortest gap: an upper edge, high by about one mean
  gap divided by the number of gaps. The interval comes from the posterior draws, so the corrected point sits
  low in it (about 0.4 sd for S in chsh-dead-time)." CHANGELOG: "the corrected S 0.26 sigma low (+-0.29)".
- **L2** (README "Unequal output efficiencies" bullet; `corrected_chsh` docstring): "Corrected S does not
  correct detector efficiencies: an output imbalance biases it as it does raw S (0.6/0.4 on both sides: corrected
  S 2.033 for a 1.980 source, 13 sigma high). The 'output imbalance' label applies to both values."
- **L3** (README Bell limits, after the corrected-CHSH paragraph): "The predictive check fits a Gaussian peak.
  A timing tail of a few percent (SPAD diffusion, SNSPD) fails it on a Bell run's 16 histograms even when the
  delay artifacts pass (3% of clicks with a 300 ps tail: 6 of 8 runs refused), although a tail common to all
  channels does not bias S. Until the IRF tail (D8) lands, real-detector Bell runs will usually need
  `--allow misfit_forward_model=...`, which also stops the check flagging crosstalk."
- **L4** (propagate.py module docstring, "Disjointness"): "Disjointness compares whole tag-array and file
  hashes: a calibration made on a slice or a copy with any change of the analysed run is not caught."
- **L5** (propagate.py docstring, README calibration paragraph): "One artifact per kind and channel: with
  calibrations before and after the run, pass the one to use; the drift between them is not propagated."
- **L6** (propagation.py module docstring): "The gate writes 400 draws per artifact (users get 2,000). With
  N_MC = 4,000 joint draws the interval ends have a Monte Carlo sd of about 4% of the joint standard error,
  the same at 400 and 2,000 rows."
- **L7** (propagation.py, comment on chsh-dead-time): "About 57% of its 0.1 s dark artifacts give the dead
  time as a bound (median 60 ns for 50 ns), so nearly every replicate corrects at least one channel at its
  bound: the regime tests the corrected value at the bound and the identified channels' draws, not the D3
  range from zero to the bound."
- **D3 label** (forward.py `DEAD_TIME_BOUND`, replacing the current text): "dead time known only as an upper
  limit: value and interval at the bound; counts give the value at zero (the truth lies between)".
- **forward_check.py docstring** ("What it does not check"): "For CHSH the calibrated delays are compared
  with the analysed peaks through the window capture (`chsh`); for HOM and fringe a drift between calibration
  and run is not seen."

## Implementation order (adopted fixes)

No pinned method version changes (no file in `python/dadl/calibration/` is edited). The coverage gate is
not affected. **The propagation gate needs a re-run of all five regimes** (hom-snspd, hom-bright,
hom-crosstalk, chsh-spad, chsh-dead-time) at 1,000 replicates on the CI seed: Fix 4 changes the predictive
check's p-values wherever the floor start was used (every regime may be touched; chsh-spad and hom-snspd
most), and Fix 5 adds a flag source in both CHSH regimes. At 10-40 min per regime locally with all cores
this is about 1-3 hours: affordable; run it once after Fixes 1-10 are in. Fixes 1, 2 and 6-10 change labels,
counts or the two-point point at a bound only (no gated quantity), so alone they would not need it.

1. **D3 made consistent (R1, B2 + R2).** `forward._dead_times`: return the bound as the point
   (`out.append((np.full(n, hi), hi, hi))` in s) so value, interval and label agree. New
   `forward.chsh_bound_range(counts, channels, signs, *, calibrations, window_ps, delays_ps)`: the point S
   (point counts and point captures) with bounded channels at 0 and at the bound, or None; `chsh._corrected`
   adds `corrected_at_zero_dead_time` and `corrected_at_bound_dead_time` to the result's counts (the CLI
   already prints them). `guardrails.threshold_labels` gains `ends: Sequence[float] = ()` (mirrored with the
   value for CHSH) and labels "threshold crossed only after correction" when `raw.value < threshold <=
   max(corrected.value, *ends)`; hom, fringe and chsh pass the two ends. New label `DEAD_TIME_RANGE_WIDE`
   ("the range from zero to the dead-time bound moves the corrected value by more than its standard error:
   quote the range") when `|at_bound - at_zero| > corrected.stderr`. Replace `DEAD_TIME_BOUND` text (above);
   cli.py:898 "(the corrected value and interval use the bound)". Tests (`test_forward.py`,
   `test_chsh_corrected.py`): with a bounded dark artifact and unequal singles, `low <= value <= high` (two-point
   and CHSH); CHSH counts carry both keys and at_zero equals S with no dark artifacts; the `p8_threshold` case
   gets "threshold crossed only after correction"; the `r2-bell/p5_snspd.py` case carries the wide-range label.
2. **Stand-in jitter (R3).** `forward._capture` and `_pair_capture`: when `"jitter_sigma_ps"` is in
   `cal.stand_in` (or `unidentified`), use `np.full(n, hi)` as the jitter draws and `hi` as the point, and also
   compute the point capture at jitter 0; when they differ by more than 1e-3 (relative) in any run add the
   label "delay jitter known only as an upper limit: window capture at the bound". Tests: `r2-prop/p2c_held.py`
   as a regression (HOM and a four-pair CHSH) returns a finite value with no KeyError.
3. **Raw capture-balance check (B1, raw path).** `forward_check.py`: factor the fit out of `pair_pvalue` into
   `fit_peak(hist, edges, tick, d0, sigma0, offsets) -> (theta, mu)` (pair_pvalue unchanged in result
   except for Fix 4) and add `pooled_peak(pairs) -> (offset_ps, sigma_ps)` that sums the histograms of
   `PairData` sharing a channel pair and declared delay and fits once. `analysis/chsh.py`: new
   `capture_balance(runs, window) -> (s_reweighted, {pair: (offset, sigma, capture)})` using `pooled_peak` per
   outcome pair (pooled over the four settings) and `forward.grid_capture`; S recomputed with each
   combination's raw count divided by its capture. In `chsh()`, record `counts["capture_offsets_ps"]` and
   `counts["capture_reweighted_s"]`; add label `CAPTURE_IMBALANCE` ("unequal window capture across outcome
   pairs: the declared delays miss the pair peaks (offsets in counts), shifting S by about {dS:+.3f}; declare
   the delays from a pair calibration") when `|S_raw - S_reweighted| > 0.5 raw.stderr` or any capture < 0.5;
   skip pairs below `MIN_COUNTS`. Tests: the `p9_misdeclared` physics with A+, B+ 400 ps late: label present,
   reweighted S within 3 se of 1.980; declared correctly: no label; 10 seeds of the chsh-spad physics: no label.
4. **Predictive-check start (H3).** `forward_check.pair_pvalue`: always add one start with the background from
   the outer bins, `max(median(hist[:3] + hist[-3:]) / width, floor)`, to the loop of starts (keep the best NLL).
   Do not treat non-convergence as "not checked". Test: the seed-2007 reference histogram of `r2-hom/p2c.py`
   (stored as a small fixture of the two tag arrays or regenerated by seed) gives p > 0.1.
5. **Calibrated delays against the analysed peaks (B1, corrected path).** In `chsh._corrected`, with
   delay-jitter artifacts: compare S reweighted by the data captures (Fix 3) with S reweighted by the
   artifacts' point captures (`_pair_capture(...)[1]`); when they differ by more than 0.5 x the corrected
   value's `stderr_statistical`, raise `Misspecified` ("the calibrated delays do not match the analysed pair
   peaks (A+/B+: calibrated ..., analysed ... ps): a drift since calibration changes the window capture") unless
   the policy allows `misfit_forward_model`, then label `FORWARD_MISFIT_ACCEPTED` (so the gate counts it as a
   flag). Update the forward_check.py docstring (above). Tests: `p10_stale` physics at 200 ps drift with the
   strict policy raises `Misspecified`; with the override carries the label; with current artifacts passes.
6. **HOM reference delay required with calibrations (H2).** `hom_from_tags`: when `calibrations` is given and
   `reference_delay_ps` is None or 0, raise `ValueError("with calibrations, pass reference_delay_ps (setup.yaml's
   analysis.reference_delay_ps): the forward model needs the reference's split")`. Fix the docstring ("without
   one" -> "required with calibrations"). Test: raises; propagation and CLI paths unchanged.
7. **Fringe drift labels (H1).** `analysis/fringe.py`: always add `fringe.TWO_POINT` ("two-point visibility at
   the declared settings: not normalised for rate drift between the max and min acquisitions"); per channel,
   compare singles rates between the two acquisitions with a Poisson z
   (`(r1 - r2) / sqrt(n1/T1^2 + n2/T2^2)`), record `singles_ratio_a`, `singles_ratio_b` in counts, and add
   `SINGLES_DRIFT` ("singles differ between the max and min acquisitions beyond 4 sigma: a rate drift or a
   setting-dependent loss; V moves with it (+3% on the max run read 0.707 for 0.699)") when either |z| > 4.
   Docstring: the singles test assumes setting-independent singles (Werner state, no analyser-dependent loss).
   Tests: the `p5_fringe_drift` physics at +6% carries SINGLES_DRIFT; at 0% it does not.
8. **Corrected V above 1 (H4).** `forward.corrected_visibility`: when the point exceeds 1, add
   `CORRECTED_ABOVE_ONE` ("corrected value above 1: the accidentals fluctuate above the minimum count; quote the
   interval"); value unclipped. README.md:197-199: "The subtraction path refuses ...; the forward model keeps
   the unclipped value and labels it." Test: `r2-hom/p7_above1.py`'s HOM case (C_dip 5, A 8) carries the label.
9. **Labels and citations.** (a) R6: `Corrected.artifact_ids` lists only the members the forward model read
   (dark artifacts of the analysed channels, delay-jitter artifacts of the pairs used); `CalibrationSet.load`
   with `kinds=None` skips `deadtime-series` for hom_visibility, bell_visibility and chsh and adds the label
   "deadtime-series artifacts present but not used: the forward model takes the dark dead time until M5"; drop
   the set-level `BOUND_ONLY` (the forward model's `DEAD_TIME_BOUND` and Fix 2's jitter label say which
   parameter); keep the constant exported for compatibility. Test: `p9_unused` cites 2 IDs; `p7_boundlabel` has
   no bound label. (b) R7: forward model labels `AFTERPULSING_UNMODELLED` ("afterpulsing reported by a dark
   calibration is not in the forward model (no kernel propagated)") when a used dark artifact has
   `afterpulse_prob` neither stand-in nor unidentified. (c) R4: `CalibrationSet` records `dataset_checked`;
   loading without a DatasetContext adds "calibrations not checked against the dataset (no DatasetContext:
   disjointness, validity and conditions not tested)"; L4 in the docstring. Test: `r2-prop/p12_nodataset.py`.
   (d) B4: `ACCIDENTALS_MODELLED` becomes "accidentals, dead time and window capture modelled; detector
   efficiencies not (an output imbalance biases corrected S as it does raw S); secondary, not a Bell claim (the
   headline is raw S)". (e) H7: `hom.TWO_POINT` adds "or the reference delay is not many coherence times"; the
   CLI warns when setup.yaml gives `source.coherence_time_ps` and `reference_delay_ps < 3 x` it. (f) H6:
   `cli._print_hom_extras` prints the reference-delay warning only when `"window_capture"` is not in
   `result.corrections` (else "the raw V reads low; the corrected value models the capture").
10. **Validity from acquisition (R5).** `artifact.build_manifest`: `valid_for.from` = `manifest["acquired"]`
    when recorded, else `created`. L5 in propagate.py and README. Test: dark acquired 09:00, created 09:30,
    valid to 18:00 accepts a dataset at 09:10 (`r2-prop/p3_validity.py`). Artifact IDs change only for artifacts
    with both `acquired` and `valid_until`; update any pinned ID that has both.
11. **Docs (B6 + H5, limits, deferrals).** README.md:57-60 (corrected CHSH and calibration uncertainty are
    done; uncertainties of corrected values include calibration), :88-95 (crosstalk flagged only with a
    delay-jitter artifact, unflagged without), a CHSH `--calibrations` paragraph (`--allow dead_time=...
    accidental_subtraction=... window_capture=...`), a corrected-CHSH limits bullet (B1 numbers, D3 range, L2,
    L3); python/README.md:8 (quickstart prints raw S and a corrected S). L1-L7 where named; CHANGELOG entries for
    Fixes 1-10 and the B3 sign. docs/phase3-plan.md smaller-items table at M10: per-detector efficiency
    normalisation (B4); D3 union interval and a chsh-dead-time regime with identified dead times or a gated
    range (B2, G; owner's agreement, next to rotating seeds). Beside the D8 row: common-vs-confined misfit test
    for CHSH (B5).
12. **Propagation gate re-run**: `python -m dadl.analysis.propagation all --replicates 1000 --jobs N` (CI
    seed), all five regimes; record coverage, control, misfit and flag shares in the PR. Expect: coverage
    unchanged within noise (no interval changes in the gate's paths; the CHSH point at a bound is unchanged),
    flag shares in clean regimes <= 1%, hom-crosstalk still >= 90%.
