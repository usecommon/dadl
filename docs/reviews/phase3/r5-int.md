> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Interval estimation with nuisance parameters: the passive-quench tau_p coverage drift

Probes: `scratchpad/council/r5-int/`. `probe.py` repeats `coverage.deadtime_coverage` replicate by
replicate (same `_children` streams, same dark and series calls). For each replicate it records the gate
interval (A), the held-branch profile interval before afterpulse widening (C), and a refit with
afterpulsing held at the true (p, tau_a) and no widening (B). Output files: `ci.jsonl` (CI seed 20261111, 0-999) and
`s4242.jsonl` (seed 4242, 0-999). `ana.py` gives coverage split by side, bias and the checks on the dark
artifact. `margin.py` measures the effect of the branch margin. `big.py` runs the same detector with
0.5 s runs, or a varied p and tau_a, with afterpulsing held at truth, and `big_*.jsonl` holds its
results. Everything runs on HEAD ef8bc0c through the package's own functions with one BLAS thread.

## Verdict

This is not replicate noise: the interval is slightly anti-conservative. Pooled over 2,000 replicates the gate interval covers **93.8%** (z = -2.46). Misses are lopsided: in **91** the interval lies below the truth, in **33** above it (about 50 each expected). `probe.py` reproduces the CI seed's 93.2% exactly. Seed 4242 gives 94.4%.

The deficit has two causes, and the afterpulse propagation is not one of them:

1. **False branch switches (about 0.45 points).** In 9 of 2,000 null replicates (7 of them on the CI seed) the
   both-branch fit beats the held fit by more than `BRANCH_MARGIN` = 3. The fit then switches to tau_p of
   about 30 ns, and the interval misses 40 ns in every one of these replicates. The held profile interval covers the truth in all 9.
2. **Unmodelled pending afterpulses (about 0.5 points at gate size, and growing with data).** The point is biased by
   -170 ps, a bias that holds its size with 25 times more data. It is -0.14 sd at the gate's 0.02 s runs and
   -0.67 sd at 0.5 s runs, where coverage drops to **88.1%** (141 of 160, all 19 misses low) and no note is given.

The earlier readings were 96.9% and 95.0%. Both were propped up by widening artifacts that later rounds
removed: the jagged snapped kernel (round 1) and the joined interval (round 3, p90 12.9 ns). The
underlying method has been about 94% throughout.

## Decomposition (pooled 2,000; CI seed in brackets)

| interval | coverage | interval above / below truth |
|---|---|---|
| A: gate (estimated afterpulse, widened) | **93.80%** [93.2] | 33 / 91 |
| C: held profile, estimated afterpulse, no widening | 92.40% [92.3] | 46 / 106 |
| B: true afterpulse held, no widening | 93.90% [93.4] | 32 / 90 |
| A with `BRANCH_MARGIN` = 5 (from records: a union in place of each switch) | **94.25%** [93.9] | |
| No afterpulsing in the detector (p = 0), same rates, 0.02 s (`big_noap_0.02`, 1,000) | **95.6%** | 20 / 24 |

- **The afterpulse propagation is honest.** The dark estimates are unbiased, and their draw sds are calibrated (z sd 1.02 for both
  p and tau_a; their intervals cover 95.3% and 94.8%). Leaving out the 69 switched or union replicates,
  the realised shift of the point between estimated and true afterpulse parameters has sd **414 ps**; the propagated
  linear shift has rms **397 ps**. The linear regression on (dp, dtau_a) leaves a 60 ps residual. A matches B
  (93.8 against 93.9). Without the widening (C), coverage falls to 92.4%. So the widening is needed, and it is
  the right size.
- **The chi-square calibration of the profile is fine at these sample sizes.** With afterpulsing off at the same rates and durations, coverage is
  95.6% with balanced misses, a bias of +61 ± 37 ps, and an empirical sd of 1,174 ps against an interval sd of 1,191 ps.
- **The cell-averaged kernel is not the cause.** With a 4 times finer grid (`BINS_PER_TAU_A` 80, `BINS_PER_MEAN_GAP` 800)
  the point moves by +4 ± 10 ps (16 paired replicates at 0.5 s).
- **Width is adequate.** Leaving out the switches, the empirical point sd is 1,217 ps against a nominal 1,352 ps (CI seed). The shortfall comes
  from location: the bias and the switches.

## Bias scaling (true afterpulse held, 0.5 s runs, 40-160 replicates each, `big.py`)

| detector | pending fraction at 6 Mcps | tau_p bias | coverage |
|---|---|---|---|
| p = 0 | 0 | -6 ± 20 ps (160) | 93.1% |
| p 2%, tau_a 10 ns | 0.11% | +13 ± 40 | 92.5% |
| p 2%, tau_a 30 ns (the regime) | 0.31% | **-163 ± 19** (160) | **88.1%** |
| p 4%, tau_a 30 ns | 0.61% | -263 ± 36 | 72.5% |
| p 2%, tau_a 90 ns | 0.70% | -808 ± 38 | 5% |

The bias grows with the pending fraction and with tau_a. At fixed data it does not shrink, so its size
relative to the sd grows as sqrt(gaps). At the CI regime's 0.02 s it is -0.14 sd. At the CLI default of 0.05 s
(`cli.py:338`) it is about -0.22 sd. At 0.5 s it is -0.67 sd.

## Findings

| # | Severity | Finding | Evidence | Where | Suggested fix |
|---|---|---|---|---|---|
| 1 | wrong | The branch margin of 3 NLL lets about 0.5% of below-maximum scans switch to the high branch (tau_p of about 30 ns, every one a miss). The docstring says the high branch "wins by at most 2.1 on 256 simulated scans". Its gain exceeds 3 in 9 of 2,000 replicates and 4 in 4 of 2,000, with a maximum of 4.97. This alone costs the gate 0.7 points on the CI seed (93.2 against 93.9). | `ci.jsonl`, `s4242.jsonl`, `margin.py`: switched replicates 165, 271, 425, 477, 494, 820 and 900 (CI seed) give tau_p 28.6-32.7 ns, intervals excluding 40, with gains 3.19-4.97, and the held profile covering in all of them. Below the maximum the gain between these two non-nested fits is about N(-D, 2D) (Cox), with D the Asimov separation. For D = 5.6 this predicts P(gain > 3) = 0.51%; observed 0.45%. | deadtime.py `BRANCH_MARGIN` (:115), and the module docstring "Two rates per gap law" | Set `BRANCH_MARGIN = 5`. Under the normal approximation the false-switch rate is Φ(-(M+D)/√(2D)). Its maximum over D, at D = M, is Φ(-√(2M)): **7.8e-4** for M = 5, against 7.2e-3 for M = 3. This holds whatever the scan's separation, so the margin does not have to be tuned per scan. Measured effect, from the records (exact, because a union covers when either fit covers): 0 of 2,000 switched, coverage 93.8 → **94.25%** (CI seed 93.2 → 93.9). Round-4 past-maximum scans (1/3/60 Mcps): 28 of 48 switch, against 36. The other 8 become unions, which cover (union median width 6.0 ns against 3.5 ns switched). Docstring: "wins by more than 3 in about 0.5% of scans; with margin 5 the false-switch rate is at most about 8e-4 for any separation (Cox's normal approximation for non-nested fits)". This is a method-version bump and a gate re-run. |
| 2 | wrong | Unmodelled pending afterpulses bias tau_p low by an amount that does not shrink with data, and the 0.5% note threshold is a fixed fraction. At the regime's own 0.31% pending, with 0.5 s runs, the interval covers 88.1% with all misses low and no note. The note text ("biased low (by about 1 sd at 3%)") understates it about tenfold at moderate data. At gate size this costs about 0.5 points and explains the miss asymmetry (90 low against 32 high in B). | Bias table above: -163 ± 19 ps at 0.31% pending, unchanged from -174 at the gate size. The bias is 0 with p = 0 and 0 with tau_a 10 ns, and unchanged with a 4 times finer grid. | deadtime.py `PENDING_MAX` (:113) and its note; module docstring "Afterpulsing" | Modelling the pending candidate is already deferred to M10 (round 1, S4), and it is the real fix. This new evidence argues for doing it before 0.2.0. Until then, make the guard scale with the data: fire when pending × √(gaps) > about 2. In the regime the bias is about 0.11 sd per unit of pending × √(gaps) (0.0031 × √(1.7e5) ≈ 1.3 gives -0.14 sd; at 0.5 s, 6.3 gives -0.67 sd). Note the tau_a dependence: at 90 ns the bias is twice what the pending fraction alone predicts. Reword the note to "tau_p is biased low by an amount that does not shrink with more data (about 0.7 sd at 0.3% pending with 1e6 gaps)". Add a 0.5 s passive-quench probe or regime in M10, because the 0.02 s gate cannot see this. |
| 3 | minor | The coverage history (96.9 → 93.8 → 95.0 → 93.2) reads as drift, but every number is on the same replicates. The higher readings came from widening artifacts (the snapped kernel, the joined interval), not from correct calibration. | The two 1,000-replicate seeds: 93.2 and 94.4. In 250-replicate blocks of the CI seed: 95.2, 93.6, 91.6, 92.4. | CHANGELOG round entries | When the gate is re-run after fix 1, record the low/high miss split and the switch count beside the coverage. |

## Checked and found sound

- The dark artifact's (p, tau_a) are unbiased (p mean 0.01998, tau_a 30,018 ps), and their draw sds match the empirical sds (z sd 1.02).
- The linear afterpulse sensitivity and the side-by-side `hypot` combination: the realised shift sd is 414 ps against 397 propagated. A profiled Gaussian constraint would not change the width measurably.
- The profile chi-square(1) interval at 2e4-1.2e5 gaps per run without afterpulsing: 95.6%, balanced.
- The cell-averaged convolution kernel: grid refinement changes the point by less than 10 ps.
- The union branch (within the margin): 103 of 2,000 replicates, 99 covered (CI seed 57 of 58).
- The CI-seed gate number reproduces exactly through the probe (93.2%), so the decomposition applies to the gate itself.
