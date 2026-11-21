> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 6 (Bayesian computation, journal referee, metrology)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: R = journal referee (r6-ref.md rows 1-8 as R1-R8), T = metrologist (r6-metro.md rows 1-8 as
T1-T8), B = Bayesian computation (r6-bayes.md rows 1-8 as B1-B8). Council probes:
`scratchpad/council/r6-council/` (`ident.py` dark dead-time identification by regime, outputs `ident_*.txt`;
`tick2.py`, `tick3.py` tau_np coverage against a coarse tick, on and off the grid, output `tick2.txt` and the
transcript figures below).

## Checked numerically

All runs on HEAD 7bcc32b with the repo venv, through the package's own functions, one BLAS thread.

### R1: a falsely identified dark dead time reaches the corrected value (confirmed, wrong)

`ident.py` (dark runs `simulate_dark`, rng 10,000 + r, then `calibrate_dark`, as the referee's probe):

| setup, channel, dark run | identified | floor <= 0 (of identified) | identified intervals covering the truth |
|---|---|---|---|
| hom-snspd A, 20 s at 100 cps, tau 30 ns | 7/200 | 7 | **0/7** |
| hom-snspd A, 200 s at 100 cps | 7/100 | 7 | **0/7** |
| chsh-spad A+, 2 s at 500 cps, tau 22 ns | 24/200 | 22 | 20/22 (floor > 0: 1/2) |

Reading. When the true dead time is far below the dark gap spacing 1/(n f0) (every SNSPD), "identified" (lower
end > 0) is exactly the event that the interval misses: its nominal 2.5-5% tail, and a longer dark run does not
help (7% at 200 s). The data cannot tell that event from a real loose identification (chsh-spad 2 s: 20 of 22
loose identifications cover), so the identification rule in `dark.py` is not the place to fix it; the use-time
check is. Code read confirms the referee's chain: `forward.py:309` sets no floor when it is <= 0, so the
downward test is off; `_refuse_loose_bound` (`forward.py:580`) and the inline refusal (`:423`) apply to bounds
only; an identified dead time gets only `DEAD_FRACTION_HIGH` (`:435`, `:729`); `cli.py:613` warns only below
the declared dead time. With the lower 95% end as the floor when the 3-width floor is <= 0, the 7/7 SNSPD
artifacts are contradicted by the analysed gaps (85 kcps: about 1 - e^{-0.4} of gaps lie below 4.7 us); on the
22 loose chsh-spad artifacts the test can fire only on the 2 whose lower end exceeds the truth, which are real
misses (the refusal is overridable). The analysed data, not the dark run, carry the power, so the test is
strong exactly when the dead fraction matters. No test or gate depends on `DEAD_FRACTION_HIGH` (grep).

### T1: the tau_np interval against a coarse tick, and what the coverage gate then shows

Recorded gap d = floor(t2/T) - floor(t1/T) ticks, so the true gap lies in (d - 1, d + 1) ticks. With the true
tau_np off the tick grid, tau_np < m + T and tau_np > m - T - E (E the continuous edge excess): both ends need
one tick. `tick2.py`/`tick3.py`: continuous-time hard edge, 0.01 s at 2e6 and 5e6 cps, tags floored to T,
truth 22,000 ps (on the grid) or 22,000 + U(0, T) (off):

| T, truth | current [m - k_hi T, m - k_lo T] | widened [low - T, high + T] | upper end only + T |
|---|---|---|---|
| 80 ps, off grid | **0/200** | 200/200 | 135/200 |
| 25 ps, off grid | **0/400** | 395/400 (98.75%) | 192/400 |
| 4 ps, off grid | 325/400 (81%) | 397/400 (99.25%) | 378/400 |
| 25 ps, on grid | 386/400 (96.5%) | 400/400 | 386/400 |
| 4 ps, on grid | 194/200 | 200/200 | - |

The reviewer's case reproduces (`tick.py 21960 80`: 3/3 at [21920, 21920], 40 ps below the truth). Widening
only the upper end is not enough (48% at 25 ps): both ends must move.

**What the coverage gate then shows.** The gate simulates on the grid (integer tau_np, 1 ps tags), where the
current law is right and the widening is pure conservatism. Under the law P(D <= k) = 1 - e^{-lambda(k + 1/2)}
with the edge densities that reproduce the recorded figures (D = 0 in 79% and 45% of replicates), moving the
lower end one tick takes coverage from 99.2% to about 99.96% (series-spad) and from 97.3% to about 98.5%
(passive-quench); the on-grid probe rows agree (96.5% -> 100%). The at-least-95% gate still passes, but at
series-spad it then has no power at all: it can fail only if the law is grossly wrong. The check that matters
moves to an off-grid test (below). The gate is not widened or narrowed: it stays "at least 95%", and the
CHANGELOG records that the on-grid figures are conservative by construction.

### Other checks

- **T2: confirmed** by code read: `deadtime.py:575-596` checks only the dark artifact's kind, channel and status;
  `condition_mismatches` (`propagate.py:123`) is never called there. Note for the fix: that function treats keys
  starting `tagger_` or `bias`, `*_min_*`/`*_max_*` brackets and `temperature_K`; copied dark conditions must not
  use those prefixes, or later use-time checks would demand them of every dataset.
- **T3, T4, T5, T8: confirmed** by code read (`deadtime.py:563-564, 650, 662-663` "likelihood ratio" for
  2 ln LR; `DEFINITIONS["dead_time_ps"]` "shortest recorded gap"; `afterpulse_tau_sd_held` in ps without
  `_ps`; `artifact.py:230-240` hard-coded `uncertainty.kind`).
- **T6 key names:** `rate_min_cps`/`rate_max_cps` would match the `*_min_*` bracket rule and make every
  analysis that loads a folder with a series artifact demand `rate_cps` of its dataset. Use names outside the
  convention until M5 wires the check (ruling below).
- **R2, R3: confirmed** by code read: `cli.py:63-65` hides three labels only; `cli.py:909` adds the
  "(secondary ...)" suffix for CHSH only.
- **B1-B8:** not re-run (the reviewer's rank checks are extensive and the findings are about M5 or latent); the
  code references were spot-checked (`deadtime.py:699` draw cap, `propagate.py:383-386` independent
  permutations, `forward.py:18` closed form).

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| R1 | Falsely identified dark dead time (SNSPD: 7/200, 0/7 cover) reaches corrected HOM unrefused | Y/Y/Y | ADOPT-FIX | Fix 1 (lower-end floor, identified dead-fraction refusal, symmetric CLI warning) |
| R2 | Corrected-value labels printed on the raw line, one flat `labels` list | Y/Y/Y | ADOPT-FIX | Fix 5 |
| R3 | Forward-model HOM/fringe not marked secondary | Y/Y/Y | ADOPT-FIX | Fix 6 |
| R4 | "Calibration uncertainty" does not say what it excludes; prints 0.0000 | Y/Y/Y | ADOPT-FIX + ADOPT-LIMIT | Fix 7; wording L2 |
| R5 | README insufficient for a methods section | Y/Y/Y | ADOPT-FIX (docs) | Fix 9: methods section |
| R6 | `--heralding` writes A+/B+ only; four pairs needed for CHSH capture | N/N/Y (option now); Y/Y/Y (wording) | ADOPT-LIMIT + DEFER (`--pair`, M10) | Wording L4; M10 smaller item: `simulate --heralding --pair` |
| R7 | HOM refusals reuse CHSH wording; Python-only remedy | Y/Y/Y | ADOPT-FIX (text) | Fix 8 |
| R8 | HOM prints no selection caveat | Y/Y/Y | ADOPT-FIX | Fix 8 |
| T1 | tau_np interval zero width; misses an off-grid truth every time | Y/Y/Y | ADOPT-FIX | Fix 2 (both ends one tick, draws, `tag_tick_ps`, off-grid test); dark dead time: wording L3 (no `dark.py` change) |
| T2 | Dark artifact at other bias/temperature used silently by `calibrate deadtime` | Y/Y/Y | ADOPT-FIX | Fix 3 |
| T3 | No readable uncertainty budget for tau_p | Y/Y/Y | ADOPT-FIX | Fix 4(a) |
| T4 | "likelihood ratio" printed for 2 ln LR | Y/Y/Y | ADOPT-FIX | Fix 4(b), key renamed in the same bump |
| T5 | tau_np definition mixes quantity and estimator | Y/Y/Y | ADOPT-FIX (text) | Fix 4(c) |
| T6 | Validity range not machine-readable | Y/Y/N (record now, check at M5) | ADOPT-FIX (record) + DEFER (check, M5) | Fix 4(d); M5 entry item |
| T7 | Certificate lines: no level, too many digits, limit level | Y/Y/Y (a, b); Y/N/Y (c wording) | ADOPT-FIX (a, b) + ADOPT-LIMIT (c) | Fix 10; wording L5; M10 item: an unidentified tau_p coverage regime |
| T8 | "input rate ... cps" is the avalanche-candidate rate; `_sd_held` unit | Y/Y/Y | ADOPT-FIX | Fix 4(e) |
| B1 | M5's named NUTS users have no likelihood | Y/Y/Y | ADOPT-FIX (plan) | Fix 11(a); scope choice is owner item 2 |
| B2 | Multimodality invisible to R-hat from MLE starts | Y/Y/Y | ADOPT-FIX (plan) | Fix 11(b) |
| B3 | Hard boundaries / support parameters vs "any divergence" | Y/Y/Y | ADOPT-FIX (plan) | Fix 11(c) |
| B4 | Stored draws imply flat priors in fit coordinates (improper in logit p) | Y/Y/Y (limit now); Y/Y/Y (priors at M5) | ADOPT-LIMIT + DEFER (M5) | Wording L6; M5 entry item: declare proper priors per calibration, record them, bump where draws move |
| B5 | Series tau_p draws lose correlation (about 0.3) with dark afterpulse draws | Y/N/N (now) | DEFER (M5 entry) | Latent (`SERIES_UNUSED`); needs a shared row index in `CalibrationSet.draws`, a structural change best made with the M5 reader. M5 entry item with the reviewer's test |
| B6 | Union intervals drawn as one split normal (18% of draws in the gap) | Y/Y/N | ADOPT-FIX | Fix 4(f): cheap inside the deadtime bump already required |
| B7 | "Stored draws gated" (M1) not done; no draw-level test | Y/N/N (fast-tier test now) ; Y/Y/Y (plan) | DEFER (M5) + plan fix | Fix 11(d): move the plan row from M1 to M5, review item 7.17 stays open with that target |
| B8 | Benchmark speed summary selective | Y/Y/Y | ADOPT-FIX (docs) | Fix 11(e) |

Votes on R1. P: an "identified" label with a 2.5-5% false rate in the SNSPD regime, and a printed V with
sd 6, is wrong; the hard edge makes any analysed gap below the lower end a contradiction. M: local to
`forward.py` (not pinned), no gate depends on the label being replaced; condition: re-run the propagation
gate and published results to show no new refusals. U: a 660x-too-long dead time must not pass with exit 0.

Votes on T1. P: real dead times are analogue, so off-grid is the realistic case; zero width is a claim of no
uncertainty. M: yes, but only at the series call site (the bump is due anyway); the dark interval spans many
ticks, so a `dark.py` change (pinned "7", gates) buys nothing measurable. U: yes; the CLI's
`22000 [22000, 22000]` misleads.

## Limit wordings

- **L1** (README "Calibration artifacts", near the dead-time text): "An identified dark dead time can be wrong
  when the detector's dead time is far below the spacing of dark counts (an SNSPD at 100 cps: 3-7% of dark
  runs reports a microsecond-scale dead time as identified, and a longer dark run does not help). The analysis
  tests it against the analysed gaps and refuses it when they contradict it; calibrate an SNSPD's recovery from
  an illuminated inter-arrival histogram."
- **L2** (README "Raw values ... Corrected values" paragraph, and the docstring of `corrected_visibility` /
  `corrected_chsh`): "The calibration part is conditional on the models: it includes the calibration artifacts'
  statistical uncertainty only. It excludes each artifact's `uncertainty.not_included` (copied into the record
  as `calibration_excludes`), the IRF shape beyond the Gaussian the delay check accepted (HOM and fringe capture),
  efficiencies, afterpulsing, multi-pair emission, the splitter ratio, drift between calibration and run below
  the drift threshold, and the choice among several valid artifacts."
- **L3** (dark.py is pinned: README dead-time section and the `deadtime.py` docstring "Fit"): "Dead times are
  read on the tagger's tick. The series tau_np interval allows one tick either side (a recorded gap is within one
  tick of the true gap). The dark dead-time interval does not; it spans many ticks in any dark run, so the
  difference is below its resolution."
- **L4** (README "For a CHSH setup", and the `calibrate pairs` / quickstart printed hint): "Window capture for
  CHSH needs a delay-jitter artifact for each of the four channel pairs, so four heralding runs (re-route the
  fibres between them). `dadl simulate --heralding` writes the A+/B+ pair only."
- **L5** (`deadtime.py` docstring "Intervals" and the unidentified note): "An upper limit on tau_p uses the
  two-sided 95% threshold (2 ln LR = 3.84); on the tau_p >= 0 boundary that is a one-sided limit of about 97.5%.
  No coverage regime is unidentified, so this level is not gated."
- **L6** (README "Calibration draws", and `docs/phase3-plan.md` M5): "Stored draws are posteriors under flat
  priors in the fit's coordinates (log rate, logit p, log tau_a) inside the fit's box; rank checks confirm this
  (round 6). Near p = 0 that prior is improper without the box; M5 declares proper priors."

## Implementation order (adopted fixes)

1. **Fix 1 (R1, forward.py, cli.py).** (a) `check_dead_time_gaps`: when the 3-width floor is <= 0 and the dead
   time is identified, use the 95% lower end `tau.low` as the floor (counts record it with a new
   `gap_check_floor_kind` per channel: "3 widths" or "lower end"; None only for bounds). Same count rule
   max(`GAP_COUNT`, `GAP_SHARE` x gaps), same `StaleCalibration`/override. Docstring: the hard edge admits no
   analysed gap below the true dead time, so gaps below the lower end mean the interval missed; L1. (b) Refuse an
   identified dead time whose dead fraction at the analysed rate (point x recorded singles rate) exceeds
   `MAX_DEAD_FRACTION`, as for bounds (`GuardrailRefusal`, message: the first-order live fraction does not hold
   there; a series dead-time calibration is the route, M5); drop the identified `DEAD_FRACTION_HIGH` label
   paths (`:435`, `:729`) or keep the constant for records only. (c) `calibrate dark`: also warn when the lower
   end exceeds 3x the declared dead time ("... far above the declared ...; with a dead time far below the dark
   gap spacing, an 'identified' value is usually a chance miss: see L1"). Tests: the 7 SNSPD seeds
   (`calibrate_dark(simulate_dark(det, 20, default_rng(10_000 + r)), ..., seed=r)`, r in the `ident_hom20`
   set) are refused against the referee's hom tags (one seed suffices in the fast tier, all 7 slow-marked);
   the 20 covering chsh-spad 2 s artifacts are not refused against a chsh-spad analysed run; an identified
   22 ns at 1e7 cps is refused by (b); CLI warning text on the referee's channel A artifact.
2. **Fix 2 (T1, deadtime.py).** At the `_dead_time` call (`:694`): `est_np = Interval(point, max(0, low - tick),
   high + tick)` (point unchanged, clipped into the interval); draws `np_draws + tick * (u1 - u2)`, u ~ U(0, 1)
   (the triangular rounding of a difference of two floored times), then the existing tau_p cap. Record
   `conditions["tag_tick_ps"]` (not `tagger_*`: that prefix is a use-time check). `DEFINITIONS["dead_time_ps"]`
   and the docstring: "to within one tick". Tests: (a) the reviewer's case (`tick.py 21960 80`, 3 seeds) covers;
   (b) off-grid, T = 25 ps, 100 seeds of `tick2.py`: coverage >= 95% (council: 98.75% on 400); (c) the interval
   width is never 0. L3.
3. **Fix 3 (T2, deadtime.py, cli.py).** `calibrate_deadtime_series`: `mism = condition_mismatches(afterpulse.conditions,
   conditions)`; non-empty raises `StaleCalibration` naming them, unless `allow_condition_mismatch=True` (CLI
   `--allow-condition-mismatch REASON`), which adds a note with the reason. Copy the dark artifact's bias and
   temperature as `dark_bias_V`, `dark_temperature_K` (names outside the use-time prefixes). Test: the
   reviewer's 25 V/233 K dark artifact with a 30 V/300 K scan is refused, and passes with a note under the
   override.
4. **Fix 4 (T3-T6, T8, B6, deadtime.py, artifact.py).** (a) Record `dead_time_p_stat_low_ps`,
   `dead_time_p_stat_high_ps`, `dead_time_p_u_afterpulse_ps` (1 sd), `dtau_p_dp_ps`, `dtau_p_dtau_a`, and
   `afterpulse_p_tau_corr_used`; `uncertainty.kind` per kind (`deadtime-series`: "statistical, with afterpulse
   parameters propagated from the cited dark artifact by linear sensitivity"). (b) Notes and docstring
   "2 ln(likelihood ratio) = ..."; rename `lr_hard_edge` -> `two_ln_lr_hard_edge` (nothing reads it before
   M5). (c) T5 definition as suggested. (d) Record `scan_rate_low_cps`, `scan_rate_high_cps` (recorded count
   rates) and `scan_duration_low_s`, `scan_duration_high_s` (names outside the `_min_`/`_max_` convention until
   M5 wires the check). (e) Note text "avalanche-candidate rate (efficiency x photon rate + dark), s^-1";
   `afterpulse_tau_sd_held` -> `afterpulse_tau_ps_sd_held`. (f) Union intervals: draws from a two-component
   mixture, one split normal per branch, equal weights (conservative). Tests: (a) on the reviewer's `pq.py`
   replicate the stat ends are [28421.5, 41770.3] within 1 ps and the combined interval unchanged; (f) on
   `r6-bayes/union.py` replicate 158 under 2% of draws fall between the branches.
5. **Fix 5 (R2, result.py, analysis/*, forward.py, cli.py).** Record `labels` (raw and shared) and
   `corrected_labels` separately; the CLI prints the second under the corrected line; threshold labels are
   evaluated per value; calibration labels only for artifacts in `artifact_ids`. Test: the referee's HOM run
   (finding 1's artifact) has the threshold label in `corrected_labels` only; CHSH raw line has no
   "no validity range" label for unused artifacts. Update `test_published_results` expectations if their label
   lists move (values unchanged).
6. **Fix 6 (R3).** `corrected_visibility` adds the label "secondary; the headline is raw V" (HOM and fringe);
   the CLI prints "(secondary; not an indistinguishability claim)" on the HOM forward-model line. Test on both.
7. **Fix 7 (R4).** `counts["calibration_excludes"]` = union of the used artifacts' `uncertainty.not_included`;
   the CLI prints "calibration part excludes: ..." on one line; print the stat and calibration parts with 2
   significant figures (not "0.0000"). L2.
8. **Fix 8 (R7, R8).** Analysis-specific refusal text in `guardrails.py` (HOM: "corrected values are secondary;
   the headline V uses raw counts"); the bound refusal's remedy adds "on the CLI, move the dark artifacts out of
   the calibrations folder". HOM analyze prints "the dip and reference delays and the window must be declared
   before acquisition; choosing them from the data biases V upward".
9. **Fix 9 (R5, docs).** A "Methods (corrected values)" section in the README (or docs/methods.md, linked):
   per-outcome model s = (C - A)/(T L capture), L as in `forward.py`; capture from the fixed-grid formula with
   the fitted delay and jitter; Jeffreys Gamma(c + 1/2) count draws; 4,000 draws; the point is the plug-in at
   the artifacts' points, not the Monte Carlo centre (the interval need not be centred on it); 2.5/97.5%
   quantiles; stat and calibration parts as sds with the other held fixed (they need not add in quadrature);
   cite seed and `n_draws` from `monte_carlo`; Monte Carlo noise of the ends about 0.1 stat sd. L2 here too.
10. **Fix 10 (T7 a, b; base.py, cli.py).** Print "95% interval [..]" / "95% upper limit (..)" in
    `Interval.__str__`/the calibration `__str__`; round the printed uncertainty to 2 significant figures and the
    value to match (display only; stored values unchanged). L5.
11. **Fix 11 (plan and benchmark docs).** `docs/phase3-plan.md` M5: (a) entry item: name the model NUTS samples
    (likelihood, parameters, priors) or state "calibration fits only" (owner item 2), and re-run `bench/ppl` on
    it; (b) chains from overdispersed starts (Laplace covariance x 4-9, plus one per known branch or mode), branch
    enumeration and unions kept outside the sampler, the start rule recorded in `checks`; (c) dead time stays
    conditioned on m with exact draws; tau_p = tau_np sigmoid(z) and smooth transforms for the boxes, with proper
    priors; divergences recorded but the "any divergence" rule applied only after reparametrization is tested on
    passive-quench and series-spad; (d) the "stored draws gated" row moves from M1 to M5 (a rank-statistic test
    per propagated calibration, as `r6-bayes/sbc_dark.py`, run with `--sampler nuts`), 7.17 open with that
    target; entry items for B4 (priors), B5 (shared row index and the correlation test), T6 (rate bracket check,
    renaming the keys to the bracket convention), and the r5 V4 table (still listed). (e) Benchmark record: the
    reviewer's clause on steady-state speed. M10 smaller items: `simulate --heralding --pair` (R6); an
    unidentified tau_p coverage regime (T7 c).
12. **Method version and gates.** deadtime `METHOD_VERSION = "6"` (Fixes 2-4 together; one bump); regenerate
    `python/tests/method_versions.json`. `dark.py` is not touched (stays "7"). Re-run:
    - `python -m dadl.calibration.coverage deadtime --replicates 1000` on the CI seed, both regimes. tau_p: Fixes
      2-4 do not move tau_p's interval (2 adds one tick to tau_np's, which caps tau_p's upper end; expect the
      round-5 figures within a replicate or two). tau_np: expect about 99.9-100% (series-spad) and 98.5%
      (passive-quench), gate "at least 95%" unchanged; record in the CHANGELOG that on the simulator's grid the
      figure is conservative by construction and that the off-grid test (Fix 2b) checks the tick term.
    - The propagation gate and `tests/test_published_results.py`: Fix 1 adds refusals and Fixes 5-7 change
      records; the gate must show no new refusals (a dark artifact whose lower end exceeds its truth would now be
      refused: report any to the council rather than loosening the rule) and unchanged values.
    - The unit suite, including the new tests above.

## Owner items

1. **Behaviour change (R1 b).** An identified dead time above a 20% dead fraction was labelled and is now
   refused with no override; users at high rates (a 22 ns SPAD above about 9 Mcps) lose the corrected value
   until M5's series route. The owner should accept this in the CHANGELOG's release notes.
2. **M5 scope (B1).** Whether M5's NUTS targets a new joint likelihood (dip/fringe histogram with the IRF, or
   calibration and counts together) or calibration fits only is a plan decision; the council's plan text offers
   both and names the benchmark re-run either way.
3. **Carried from round 5:** the dev version (r5 owner item 1) is still 0.1.1 while method versions advance
   (deadtime 5 -> 6 here); the schema `$id` and D4 text (r5 owner items 2-3).

## Addendum: R1's lower-end floor breaks chsh-dead-time (ruling on the implementer's report, HEAD 0ec2a0b)

Probes in `r6-council/`: `regime_dark.py` (the regime's dark runs, output `regime_dark.json`), `snspd_pts.py`
(output `snspd_pts.txt`), `gate_revert.py` + `run_gate.py` (the propagation regime at HEAD 0ec2a0b with only
the lower-end floor removed, outputs `gate_*.json`).

### (a) Root cause: the council's rule, not dark.py

`regime_dark.py`: 800 dark artifacts of the chsh-dead-time detector (20 kcps, tau 50 ns, no afterpulsing,
0.1 s: about 2,000 gaps, dark gap spacing 1/(n f0) about 25 ns):

| identified | afterpulsing detected | identified with 3-width floor <= 0 | lower end above the truth | P(at least one of 4 channels) |
|---|---|---|---|---|
| 275/800 (34%) | 0 | 275/275 | **32/800 (4.0%)** | 15.1% |

The misses look like the implementer's (lower ends 51-113 ns, points 128-190 ns). `calibrate_dark` behaves as
designed. With afterpulsing not detected (`afterpulsing_ruled_out=False`) the interval is [m - k(0.95), m], all
5% in the lower tail, and its lower end sits above the truth in 4.0% of artifacts against a nominal 5%. The
"21%" is per replicate: a CHSH replicate has four dark artifacts, so 1 - 0.96^4 = 15%, and 42/200 is 21 +- 3%.
The dark coverage gate measures exactly this interval. This regime adds nothing new to it: it is the
documented short-run, undetected-afterpulsing case.

The defect is in the council's R1 ruling. Round 4 put the floor three widths out so that a working artifact is
refused about 1 in 10^4 times. R1(a) moved the floor to the 95% lower end whenever three widths reach zero.
That refuses the interval's nominal 5% tail per channel, about 15-20% of four-channel analyses with short dark
runs. The council's own probe showed this (2 of 22 loose chsh-spad artifacts), but the round-6 ruling
overlooked that it compounds over channels. A floor at or below zero means the dark run cannot exclude a
zero dead time at the 10^-4 level. No downward test with round 4's false-refusal rate exists for such an
artifact. The SNSPD false identifications (R1) are the same 5% tail event, so the data cannot tell them apart
from honest misses. What separates them is the consequence, which R1(b) already tests.

`snspd_pts.py` (the 7 SNSPD false identifications): points 15.2-26.3 us. R1(b) (an identified dead time above
a 20% dead fraction is refused) refuses all 7 at any analysed rate above 7.6-13.1 kcps per channel, including
the referee's 85 kcps run. Below those rates the live-fraction error is under 20% but can still matter (for
example 9% at 5 kcps). That residue needs a label, not a refusal.

### (b) Ruling

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| R1-a' | The lower-end floor refuses nominal 5% misses per channel (15-21% of chsh-dead-time replicates) | Y/Y/Y | ADOPT-FIX: withdraw R1(a)'s refusal | Fix 1' below. dark.py is unchanged and stays at "7" |
| R1-est | Change the dark estimator (do not identify here, or widen) | N/N/N | REJECT | The interval covers as stated (4.0% lower misses for a nominal 5%). Widening breaks its 95% claim and its gate. Reporting loose identifications as bounds would correct them at the bound (D3), where a miss keeps its m and so its bias, and the identified channels' information is lost |
| R1-regime | Keep the refusal and change the regime (longer dark runs) | N/N/N | REJECT | It would still refuse about 15% of real four-channel analyses with short dark runs. The regimes are declared (plan row at line 360, owner-only). The rule is wrong, not the regime |

Votes. P: 10^-4 is the standard for a downward refusal, and no downward test reaches it for a floor at or
below zero. The gross SNSPD case is a consequence problem, and R1(b) catches it. M: refusals that fire on 15%
of honest analyses fail the round-4 principle, and the gate already shows that the joint Monte Carlo absorbs
the nominal tail. U: an informative label is better than refusing one CHSH analysis in five. A user below the
R1(b) rate still learns that the analysed gaps contradict the dark artifact.

**Fix 1' (forward.py, README, CHANGELOG; replaces R1(a)).**
1. In `check_dead_time_gaps`, when three widths reach at or below zero, test the 95% lower end with the same
   count rule but **label, do not refuse**. The new label `DEAD_TIME_BELOW_DARK` reads: "the analysed gaps put
   the dead time of {name} below the dark artifact's 95% lower end (a 1-in-20 chance per channel for a dark
   run this short): its dead-time correction is likely too large; a longer dark run, or an illuminated
   calibration for an SNSPD, resolves it". It goes in `corrected_labels`.
2. Counts: `gap_check_floor_kind` is "3 widths" (refusing), "lower end (label only)" or None.
3. Keep R1(b) (identified dead fraction above 20% refused) and R1(c) (the CLI warning).
4. Docstring: restore round 4's sentence that a floor at or below zero is not refused. Add: "the 95% lower end
   is then tested for a label only: it misses in 5% of channels by design".
5. L1 in the README and the docstring changes from "refuses it when they contradict it" to "refuses it when it
   implies a dead fraction above 20% at the analysed rate, and labels it when the analysed gaps fall below its
   95% lower end".
6. Tests:
   - The 7 SNSPD seeds (rng 10,000 + r, r in 58, 64, 73, 99, 113, 157, 185) are refused against the referee's
     85 kcps hom tags by R1(b). One seed runs in the fast tier, all 7 are slow-marked.
   - SNSPD seed 113 against tags thinned to 5 kcps passes with `DEAD_TIME_BELOW_DARK`.
   - The 20 covering chsh-spad 2 s artifacts are neither refused nor labelled.
   - The chsh-dead-time replicates the implementer listed run to completion (no refusal).
7. The CHANGELOG round-6 entry replaces the "does not complete" text with the figures below.

### (c) Expected gate outcomes (measured)

The propagation regime chsh-dead-time ran at HEAD 0ec2a0b with only the lower-end floor removed (equivalent
to Fix 1' for gating: labels do not change coverage). The first 200 replicates reproduce the implementer's
forced-through run exactly.

| replicates | corrected | control | misfit artifacts | forward misfit | band (lower end: control must fall below it) |
|---|---|---|---|---|---|
| 0-199 | 191 (95.5%) | 170 (85.0%) | 2 | 0 | 89.5-99.0% |
| 0-399 | **382 (95.5%)** | **350 (87.5%)** | 3 (0.8%) | 0 | 91.0-98.25% |

Expected outcomes:
- **Propagation, all regimes.** The other four regimes stay as 0ec2a0b recorded (R1(a) never fired there).
  CI's 1,000-replicate chsh-dead-time is expected to pass, with corrected about 95.5% (band 92.6-97.1%) and
  control about 87% (must be below 92.6%).
- **R1(b).** It refused none of the 400 replicates. The regime's misses imply dead fractions below 20%.
- **Coverage gates.** Unaffected. dark.py is unchanged ("7"), deadtime stays "6", and no coverage regime
  uses `check_dead_time_gaps`. No re-run beyond the round-6 figures already recorded.
- **Published results and unit suite.** Must pass. The published results gain no label unless their analysed
  gaps contradict a loose dark artifact; any such label is reported, not suppressed.

No band changes.

### Implementation list (addendum)

1. Fix 1' in `forward.py`: lower end becomes label-only, add the `DEAD_TIME_BELOW_DARK` constant, the floor
   kinds and the docstring.
2. README L1 and the CHANGELOG round-6 entry: new wording, and the measured chsh-dead-time figures replace
   "does not complete".
3. Tests as in Fix 1' item 6. Remove or convert any test that expected the lower-end refusal.
4. Re-run: propagation `chsh --regimes chsh-dead-time` at 1,000 replicates (CI's eight shards), plus the
   other regimes at 200 as a check, `test_published_results` and the unit suite.
5. No method-version change.
