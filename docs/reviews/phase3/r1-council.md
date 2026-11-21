> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 1 (M4 series dead time and gates)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: S = point-process statistician (r1-stat.md), D = detector physicist (r1-phys.md), G = gate
statistician (r1-gate.md). Council probes: `scratchpad/council/r1-council/`.

## Checked numerically

All runs on branch m4/series-deadtime at 0fd3f20 with the repo venv, through the package's own functions.

- **S1 (interruption check spans runs): confirmed.** `r1-stat/p9_offset.py` rerun: same origin ppc 0.49 ok;
  runs 30 ms apart on one clock ppc 0 misfit; durations 20/20/15 ms ppc 0 misfit. Code: `deadtime.py` calls
  `with_interruption_check(..., *tags)` -> `stationarity.channels_pvalue`, which bounds all streams by the
  earliest and latest tag of all. Council fix tested (`r1-council/int_fix.py`): each run tested on its own
  bounds, Bonferroni over runs. Same/offset/unequal all give ppc 0.49 ok; a 20% hole cut inside one run still
  gives ppc 0, misfit.
- **S2 (survival tail near n tau_p = 1): confirmed.** `r1-stat/p1_survival.py` rerun against an 80-digit
  series: relative error of S at w = 20 tau_p is -2.7% (x = 0.8), -17% (0.9), -59% (1.0), -18% (1.1),
  -6.9e-4 (1.5). The second real root of lambda e^{-lambda tau} = a merges with the first at x = 1, so it has
  not decayed at 8 tau_p. Council fix tested (`r1-council/surv_fix.py`): keep the exact series to K tau_p,
  continue with both real modes (W0 and W_{-1} branches) matched to S(u0) and S'(u0) = -a S(u0 - tau); the
  divided difference written with `expm1` so it is stable as the roots merge, and lambda1 = lambda2 = 1/tau at
  the branch point (scipy's `lambertw(-1/e)` returns nan). Relative error <= 4.4e-8 for every x in
  0.2-3.0 including 0.99, 1.0, 1.01, up to w = 40 tau. Coverage effect not rerun by the council; the
  reviewer's K = 40 runs (79 -> 91/100, 83 -> 97/100) stand as the evidence; the fix is more accurate than K = 40.
- **S3 (jagged profile with afterpulsing): confirmed.** `r1-stat/p11_prof2.py` rerun: sd of the second
  difference of the profile NLL 0.002 (no afterpulsing), 1.381 (p = 0.05, tau_a = 20 ns, h = 1 ns), 0.009
  (h / 10). Council fix tested (`r1-council/prof_fix.py`): replace the point kernel q(lag) by its cell average
  over [lag - h/2, lag + h/2], which is exact and smooth in tau because int a S(t - tau) 1{t >= tau} dt over a
  cell is S(lo) - S(hi) (S = 1 below tau). Same data: sd of second difference 0.007, profile smooth, minimum at
  19.3 ns, matching the h/10 run.
- **S4 (pending afterpulses unmodelled): confirmed.** `r1-stat/p6_apfit.py` rerun with 20 replicates:
  p = 0.05, tau_a = 200 ns: covered 14/20, mean 18936 ps (bias -7.8 se), ppc < 0.01 in 0/20.
  p = 0.02, tau_a = 30 ns: 19/20, bias -0.8 se. Pending fraction p n tau_a / (1 + n tau_a) at 1e7 cps is 3.3%
  and 0.46%; the CI passive-quench regime peaks at 0.31%, so a 0.5% guard does not fire in CI.
- **D1 (intensity modulation): confirmed.** `r1-phys/mod.py` rerun: tau_p truth 20 ns, 100 Hz +-20% gives
  17974 [16608, 19305] ps, ppc 0.106 (passes, truth outside); +-40% is flagged (ppc 2.8e-91).
- **D4 (forward model ignores tau_p): confirmed.** Inverting `closed_form.series_measured_rate`, the forward
  live fraction 1 - m tau_np overestimates m/n by 0.10%/0.52% (20/50 ns at dead fraction 0.1/0.2),
  0.51%/2.82% (40/45 ns), 0.65%/3.67% (22/22 ns).
- **G1 (tau_np gate text): confirmed by arithmetic.** At 1e7 cps with tau_p 20 ns the edge density is
  n e^{-n tau_p} = 8.2e-6 per ps over about 1.3e5 gaps, about 1 gap per tick, summing to about 1.5 over the
  three runs: not "tens of thousands". The text sits in `coverage.deadtime_coverage`'s docstring and
  CHANGELOG.md:112-113 (not in dark.py, which need not change).
- **Method versions.** `tests/test_method_versions.py` hashes each calibration module's whole file, so any
  edit to `calibration/deadtime.py` (code or docstring) needs `METHOD_VERSION = "2"` and a new pin in
  `python/tests/method_versions.json`. No adopted fix touches `dark.py` or any other pinned module.
- Not rerun (accepted on the reviewers' evidence and code reading): S5, S6/D8, S7, D2, D3, D5, D6, D7, G2-G4.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| S1 | Interruption check treats separate scan runs as one acquisition: unequal durations or one running clock give status misfit on good data | Y/Y/Y | ADOPT-FIX | Fix 1. Method version deadtime 1 -> 2; coverage re-run (series-spad, passive-quench) |
| S2 | `survival` drops the second real mode past 8 tau_p; tail 20-80% low near n tau_p = 1; tau_p undercovers | Y/Y/Y | ADOPT-FIX | Fix 2. Same version bump and coverage re-run (CI regimes sit at n tau_p <= 0.24, so expected unchanged) |
| S3 | Afterpulse kernel jump snapped to the grid: profile NLL jagged by units, intervals partly grid noise | Y/Y/Y | ADOPT-FIX | Fix 3. Same version bump; coverage re-run matters for passive-quench (afterpulsing) |
| S4 | Afterpulses pending at the next gap not modelled: tau_p biased low (-1.3 sd at p 5%, tau_a 200 ns, 1e7 cps), status ok | Y/Y/Y (guard); Y/N/N (model now) | ADOPT-FIX (guard) + DEFER (model) | Fix 4: per-run pending fraction, note and condition above 0.5%. Modelling the pending candidate deferred to M10's smaller-items table (no milestone owns it; ships in 0.2.x if not done): it is a new gap law, a new regime and a gate re-run |
| S5 | Pearson check on thousands of fine groups has little power against localized misfits | Y/Y/Y (limit); N/N/Y (coarse test now) | ADOPT-LIMIT + DEFER | Wording L1 in the deadtime.py docstring. A coarse grouped test Bonferroni-combined with the fine one deferred to M10 smaller items: it changes the misfit share and needs the coverage gate re-run |
| S6 + D8 | tau_p point (and upper limit) can exceed tau_np | Y/Y/Y | ADOPT-FIX | Fix 5. Same version bump |
| S7 | With tau_p MLE = 0 the afterpulse uncertainty is not added to the upper limit | Y/N/Y | ADOPT-FIX | Fix 6. Same version bump |
| D1 | Slow intensity modulation biases tau_p low and passes the check; tolerance unstated | Y/Y/Y (limit); Y/N/N (Fano check now) | ADOPT-LIMIT + DEFER | Wording L2 (README, deadtime.py docstring, `calibrate deadtime` help). Count-variance stationarity check deferred to M10 smaller items (needs the model's renewal Fano factor and a gate re-run) |
| D2 | Series order (P then NP) is not an actively quenched SPAD's (NP hold-off then paralyzable reset window); regime comment says it is | Y/Y/Y (limit); Y/N/N (NP-then-P model) | ADOPT-LIMIT | Wording L3 (README) and Fix 8 (regime comments). The reverse-order model is not scheduled: listed under "Out of phase 3" |
| D3 | README says gradual passive-quench recovery "shows as a misfit"; often it is absorbed as tau_p close to tau_np | Y/Y/Y | ADOPT-LIMIT | Wording L3 (README) |
| D4 | Forward model ignores an identified tau_p: live fraction high by up to 2.8-3.7% per channel at the 20% gate | Y/N/Y (fix now); Y/Y/Y (state size) | ADOPT-LIMIT + DEFER | Wording L4 (README, CHANGELOG). Using live = m/n from the series closed form when the set holds an identified deadtime-series artifact deferred to M5 (joint fits rebuild the forward model's dead-time inputs); needs the propagation gate re-run |
| D5 | tau_np interval excludes the truth on a real tagger (tagger jitter after the dead time) | N/Y/Y | ADOPT-LIMIT | Wording L5 (README, `DEFINITIONS["dead_time_ps"]` in deadtime.py) |
| D6 | Default scan rates count photons, not avalanches; with efficiency 0.5 the scan only reaches n tau = 0.25 and tau_p goes unidentified | Y/N/Y | ADOPT-FIX | Fix 7 (scale defaults by efficiency; docstring and help say "photon rate at the detector") |
| D7 | Unrecorded avalanches restart tau_p but never afterpulse or emit crosstalk | Y/Y/Y | ADOPT-LIMIT | Wording L6 (sim.py docstring, README) |
| G1 | tau_np at-least gate: stated reason wrong (about 1.5 gaps per tick, not tens of thousands); simulator's integer phase over-covers; gate has little power | Y/Y/Y (text); Y/N/N (gate change now) | ADOPT-FIX (text) + DEFER (gate) | Fix 9. Gating tau_np two-sided around its exact expected coverage, or a goodness-of-fit gate on the D = m - tau_np histogram, plus a random sub-tick dead-time phase in the regime: deferred to M10 next to "Rotating seeds" (a gate-rule change needs the owner's agreement before running) |
| G2 | Plan's gate rules drifted from code (raw-interval check, 2% misfit cap, tau_np one-sided) | Y/Y/Y | ADOPT-FIX | Fix 10 (docs/phase3-plan.md only) |
| G3 | `coverage.merge` takes `gate`/`at_least` from the first shard unchecked | N/Y/Y | ADOPT-FIX | Fix 11 |
| G4 | One-sided rule reuses the two-sided lower edge: false-fail 1.6e-4, not 5e-4 | Y/Y/N | ADOPT-LIMIT | Wording L7 (coverage.py) |

### Limit wordings

- **L1** (deadtime.py module docstring, end of "Predictive check"): "The check groups thousands of fine bins on
  (groups - 1) degrees of freedom, so it sees gross shape misfits (soft recovery, reset-window paralysis,
  strong modulation) but has little power against a discrepancy confined to the tail or the edge: a few-%
  misfit there can pass while biasing tau_p."
- **L2** (README "Hard-edge dead time" bullet; deadtime.py docstring after "under steady light"; `calibrate
  deadtime` help): "The scan needs steady light: an attenuated CW laser with its intensity stable to about
  1-2% over each run, not pulsed light or a mains-powered lamp. Slow intensity noise looks like less
  paralysis: a +-20% 100 Hz modulation biased tau_p 10% low with the predictive check passing; +-40% is
  flagged."
- **L3** (README, replacing "Soft (continuous) recovery ... shows as a misfit."): "The series chain is a
  phenomenological filter with hard edges. It describes a diode or front-end recovery (paralyzable) followed
  by a longer electronics or tagger dead time (non-paralyzable). It does not describe the reset-window
  paralysis of an actively quenched module (photons in hold-off do nothing; one just after it re-quenches
  silently: non-paralyzable then paralyzable), twilight counts, soft recovery of SNSPDs or the gradual
  efficiency recovery of passive quenching. With enough gaps such a detector shows as a misfit. With few gaps,
  or when the recovery is close to fully paralyzable, it can pass with an effective tau_p near tau_np that
  describes the scanned rates but not rates above them."
- **L4** (README, replacing "The analyses' forward model uses the dark calibration's non-paralyzable dead time
  only."; same sentence in CHANGELOG "Stays out"): "The analyses' forward model uses the dark calibration's
  non-paralyzable dead time only. With paralysis the live fraction is overestimated: for tau_p close to
  tau_np by about 0.5% at a 10% dead fraction and 3-4% at the 20% limit, per channel (about twice that for a
  pair); ratios such as S and V cancel much of it."
- **L5** (README; `DEFINITIONS["dead_time_ps"]` gains "...; the shortest recorded gap, so it includes the
  tagger's timing jitter"): "tau_np is the lower edge of the dead time as recorded. A real tagger adds its own
  jitter after the dead time, so the shortest gap sits about 4 sqrt(2) sigma_tagger below the detector's dead
  time (about 0.1 ns for 20 ps rms); the interval does not cover that offset."
- **L6** (sim.py module docstring, README): "Avalanche candidates blocked by the dead time restart the
  paralyzable stage but neither afterpulse nor emit crosstalk."
- **L7** (coverage.py, beside the at-least rule in `CoverageResult._ok` or `deadtime_coverage`): "The at-least
  rule keeps the lower edge of the two-sided band, so its false-fail rate is about 1.6e-4 at 1,000 replicates,
  below the 5e-4 budget: conservative."

## Implementation order (adopted fixes)

Fixes 1-6 edit `python/dadl/calibration/deadtime.py`: do them in one PR with `METHOD_VERSION = "2"`, the new
pin in `python/tests/method_versions.json`, the CHANGELOG entry, and the coverage gate re-run for
`series-spad` and `passive-quench` (1,000 replicates each, CI seeds). The propagation gate is not affected
(it does not use deadtime-series artifacts).

1. **Per-run interruption check (S1).** `calibration/base.py`: add `with_interruption_check_runs(ppc, notes,
   runs)` (or a keyword `separate=True` on `with_interruption_check`) that computes, for each run,
   `stationarity.channels_pvalue(t)` on that run alone and combines `p_int = min(1, R * min(p_i))`; same note
   and `min(1, 2 min(ppc, p_int))` as now. `deadtime.calibrate_deadtime_series` calls it instead of
   `with_interruption_check(..., *tags)`. Tests (`tests/test_calibration_deadtime.py`): runs of 20/20/15 ms,
   and the same runs offset by 30 ms each on one clock, give status ok (ppc as for the aligned runs); a run
   with a 20% hole cut out still gives misfit.
2. **Two-mode survival tail (S2).** `deadtime.survival`: past u0 = K_TERMS tau, use
   `S = s0 e^{-l1 t} + (ds0 + l1 s0) e^{-l1 t} (-expm1(-(l2 - l1) t)) / (l2 - l1)` with t = w - u0,
   ds0 = -a series(u0 - tau), l1 = -W0(-a tau)/tau, l2 = -W_{-1}(-a tau)/tau (argument clamped to >= -1/e;
   l1 = l2 = 1/tau when lambertw is not finite or |l2 - l1| t < 1e-9, where the term is t e^{-l1 t}). Keep
   `_lambda` for other callers or drop it. Fix the docstring sentence on the decayed modes. Test: S against an
   80-digit `decimal` series at n tau in {0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2} and w up to 40 tau, relative
   error < 1e-6 (fails now at 0.9-1.1). Existing continuity test must still pass.
3. **Cell-averaged afterpulse kernel (S3).** `deadtime._Run._afterpulse`: replace
   `q = np.where(lags >= tau, a * survival(lags - tau), 0)` by
   `q = (survival(max(lags - h/2, 0)) - survival(lags + h/2)) / h`, and the same for `qt` in the tail branch
   (on t = u - w). Test: with p = 0.05, tau_a = 20 ns, tau_p 20 ns, rates 2/5/10 Mcps, 20 ms, seed fixed, the
   profile NLL on 18-22 ns in 100 ps steps has second-difference sd < 0.05 (now 1.38).
4. **Pending-afterpulse guard (S4).** `calibrate_deadtime_series`: when p > 0, compute
   `pending = max(p n_i tau_a / (1 + n_i tau_a))` over the fitted rates; record
   `conditions["afterpulse_pending_max"]`; above `PENDING_MAX = 0.005` add the note "about {pending:.1%} of
   gaps start with an afterpulse still pending, which the gap law does not model: tau_p is biased low
   (by about 1 sd at 3%); scan at lower rates or with a detector whose afterpulsing is faster". Also state it
   in the README bullet. Test: p = 0.05, tau_a = 200 ns at 1e7 cps carries the note and condition; the
   passive-quench regime does not.
5. **tau_p capped at tau_np (S6, D8).** After `est_np` is computed: `value = min(tau, est_np.value)`,
   `high = min(high, est_np.high)` (limit and identified cases), `tau_draws = np.minimum(tau_draws, np_draws)`
   in place of raising np_draws; when the profile interval reaches tau_hi add the note "tau_p reaches tau_np:
   the data do not separate the stages (consistent with a fully paralyzable dead time)". Test: the
   `r1-phys/low.py` seed-2 case (0.5/1/2 Mcps, tau_p 10 ns, tau_np 22 ns) and tau_p = tau_np = 50 ns give
   value and high <= dead_time_ps value and high, and every draw row has tau_p <= tau_np.
6. **Afterpulse uncertainty on the upper limit (S7).** In the sensitivity block, when the fit is unidentified
   (whatever the MLE), take the sensitivity of the upper limit: refit with `_fit_tau` and `_profile_interval`
   on `_Fit(binned, p_k, ta_k)` and use `(high_k - high) / step`; identified fits keep the tau_p sensitivity.
   Remove the `tau > 0` condition. Test: a hard-edge detector with an afterpulse artifact whose p has a
   large sd gives the "afterpulse uncertainty adds" note and an upper limit above the profile's.
7. **Default scan in avalanches (D6).** `simsetup.deadtime_scan_rates`: divide each rate by the smallest
   detector efficiency (`f / (longest * 1e-12 * eff)`), so the avalanche rate x dead time spans 0.05-0.5;
   docstring and `cli.py` `--rates` help say "photon rate at each detector (before efficiency)". Test: for the
   chsh-spad example, `max(rates) * eff * tau == 0.5`; `test_cli_scan_and_calibrate` still passes. Not a
   method version (simsetup is not pinned).
8. **Regime and simulator comments (D2, D7).** `calibration/coverage.py:414` becomes "A diode or front-end
   recovery that paralyzes (tau_p 20 ns) before a longer electronics dead time (50 ns)"; the passive-quench
   comment adds "a stand-in: it checks the fit is self-consistent, not that it matches a real RC recharge";
   sim.py docstring gets L6. Write wordings L1-L7 in the places named above.
9. **tau_np gate justification (G1).** `coverage.deadtime_coverage` docstring and CHANGELOG.md:112-113: replace
   "tens of thousands of gaps fall per tick ... almost always tau_np itself" with "the edge density is about
   1 gap per tick (series-spad 1.5, passive-quench 0.7), so the shortest gap equals tau_np in about 79% and 45%
   of replicates; the simulator's integer dead time makes the tick-quantile interval over-cover (design
   97.7% and 95.8%, observed 99.2% and 97.3%), so it is gated at least 95%; this gate has little power against
   edge-density errors (M10 item)". Add the deferred items (G1 gate, S4 model, S5 coarse test, D1 Fano check)
   to the plan's "Smaller items" table at M10, D4 under M5, and D2's reverse-order model under "Out of phase 3".
10. **Plan drift (G2).** docs/phase3-plan.md: gate item 2 drops "the raw interval is checked against its own
    expectation" (or says the uncertainty-off corrected control is the only extra series); M3 acceptance says
    "misfit artifacts at most 2% and forward flags at most 1% of clean runs"; M4 tests say "tau_p band as
    phase 2; tau_np at least 95% (one-sided, see coverage.deadtime_coverage)".
11. **Merge checks the gate (G3).** `coverage.merge`: per result name, raise `ValueError(f"{name}: shards
    disagree on the gate or its rule")` when any shard's `tuple(gate)` or `bool(at_least)` differs from the
    first. Test beside `test_coverage_shards_merge_to_the_whole`: a shard record with an edited `gate` is
    refused.
