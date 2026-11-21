> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 9 (red team 2, hardware twin, performance)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: Rd = r9-red.md rows 1-5, Tw = r9-twin.md rows 1-8, Pf = r9-perf.md rows 1-8. Council probes and
outputs are in `scratchpad/council/r9-council/`:
- `offpeak_seq.py` (+ `seqlib.py`), output `offpeak_seq.txt`: sequential Monte Carlo for the off-peak dispersion
  test (null size, clean chsh-spad, power, cost)
- `null_seed.py`, output `null_seed.txt`: the same test with a fixed seed against a seed derived from the counts
- `pulsed_chsh.py`, output `pulsed_chsh.txt`: power against pulsed accidentals at chsh-spad sparsity, 16 sides
- `gapup.py`, output `gapup.txt`: the upward gap test (shipped, red team's proposal, council form) on clean
  exponential edges, the round-4 power case, and the red team's echo data (`r9-red/m_echo100`, `m_echo_few`,
  `m_clean100`)
- `tw/probes/res80.py` (copy of the twin's), output `t3straddle.txt`: T3 delay calibration, fresh seed 4242

The refusal budget of round 8 (1e-3 per analysis on clean data; off-peak 2e-4, gap tests 1e-4 + 1e-4) is kept.
No ruling below adds a statistical refusal or raises a share. New refusals are deterministic consequence rules
(Rd3) and a split of an existing override (Rd2).

## Checked numerically

All runs used HEAD 2e15eba with the repo venv, one BLAS thread, and the package's own functions.

### Pf1 (wrong, confirmed): the off-peak dispersion half cannot refuse in CHSH

Arithmetic: `offpeak_p` = min(1, 2 min(dispersion, trend)). The exact dispersion p at n_sim 20,000 is at least
1/20,001, so the doubled p is at least 1.0e-4. CHSH flags at 2e-4 / 16 = 1.25e-5. **The dispersion half can never
refuse in CHSH**; only the trend half (analytic binomial) can. HOM and fringe (2 sides, 1e-4) can reach the
threshold only when 0 of 20,000 draws exceed. Round 8's "n_sim >= 20,000" was sized for HOM and not re-derived for
16 sides (r8-council's 0 / 800 on clean chsh-spad was true but uninformative about power).

**Council design: Besag-Clifford sequential Monte Carlo, seeded from the counts.**
- Draw multinomial null tables in growing chunks (64, doubling, capped at 16,384). Stop at the h-th exceedance,
  h = 20, with p = h / L (L draws so far). Otherwise stop at n_max = ceil(4 x 2 x n_sides / alpha) draws, with
  p = (g + 1) / (n_max + 1). n_max is 640,000 for CHSH (smallest p 1.6e-6, a fifth of the per-side dispersion
  threshold 6.25e-6) and 80,000 for HOM and fringe.
- The p is exactly valid for fresh draws: P(L >= m) = h / m when the true exceedance chance is uniform. Refusal
  happens only through the n_max branch (h / L < 6.25e-6 would need L > 3.2e6), that is, when at most 2 of 640,000
  draws exceed, with chance 3 / 640,001 = 4.7e-6 per clean side, below the 6.25e-6 share.
- **Seed from the counts**, not seed 0: `SeedSequence(counts as uint32)`. The same data still give the same p.
  `null_seed.txt` shows why: with seed 0 every call reuses one correlated stream (multinomial draws at nearby
  totals share it), so the stopping point's Monte Carlo error (about 1 / sqrt(h), 22%) does not average out across
  sides. 30,000 iid Poisson sides per mean, n_sides 16, P(p <= 0.05 / 0.01 / 1e-3 / 1e-4):

  | mean per window | seed 0 | seed from counts |
  |---|---|---|
  | 0.3 | 0.033 / 0.0065 / 0.0006 / 0.00007 | 0.041 / 0.0077 / 0.0008 / 0.00003 |
  | 1.3 | **0.0595** / 0.0090 / 0.0007 / 0.00013 | 0.0475 / 0.0095 / 0.0008 / 0.00013 |
  | 3.0 | **0.0556** / 0.0079 / 0.0008 / 0.00013 | 0.0476 / 0.0095 / 0.0010 / 0.00010 |

  Seeded from the counts, every level is at or below nominal within binomial error. With seed 0 the 5% level runs
  at 5.6-6.0% (about 6 sd high at 30,000 sides).
- Clean chsh-spad (`offpeak_seq.txt`; 100 simulated analyses of the gate regime, afterpulsing on, 16 sides each):
  0 refused. Adjusted p < 0.05 in 3.0% and < 0.01 in 1.0% of analyses, the same as the shipped test. Cost per
  analysis: **20.9 ms, against 1,852 ms shipped** (88x less).
- Power, one structured side among 2 (HOM) or 16 (CHSH), 60 trials each. Each cell is shipped / sequential:

  | structured side | 2 sides | 16 sides |
  |---|---|---|
  | comb, every 4th window, mean 1.0 (the r9-perf case) | 60 / 60 | **0 / 60** |
  | comb, every 8th window, mean 0.5 | 60 / 60 | **0 / 60** |
  | 25% of windows x3, mean 1.3 | 26 / 32 | **0 / 16** |
  | 25% of windows x2, mean 1.3 | 0 / 0 | 0 / 0 |
  | 50% of windows x1.6, mean 1.3 | 0 / 0 | 0 / 0 |

  Draws per flagged side run to n_max: 2.5-3.3 s per side.
- Pulsed accidentals at chsh-spad sparsity (`pulsed_chsh.txt`). One CHSH side holds a share f of pulsed
  accidentals (window weights from the comb geometry, 50 ps jitter), and 15 sides are flat, mean 1.3. Refused out
  of 30, shipped / sequential:

  | period | f = 0.1 | 0.2 | 0.3 | 0.5 | 1.0 |
  |---|---|---|---|---|---|
  | 13157.9 ps | 0 / 0 | 0 / 0 | 0 / 8 | 0 / 30 | 0 / 30 |
  | 12345.6 ps | 0 / 0 | 0 / 2 | 0 / 17 | 0 / 30 | 0 / 30 |

  The shipped test detects none of them in CHSH. The sequential test detects a pulsed share of half or more
  reliably, and a share of 20-30% sometimes. That residual is limit L4.

Cost: under the null the test stops after a mean of 137-307 draws (0.6-2.9 ms per side), against 20,000 draws
today. Pf1(b)'s cache is then unnecessary for the analysis. The calibration path (`calibrate_accidentals`,
n_sim 2,000) is unchanged.

Where: a new `sequential_dispersion_pvalue` in `analysis/_twopoint.py`, called by `offpeak_p` with
`n_sides = len(sides)`. `calibration/accidentals.py` is not touched, so **accidentals stays "3" with an unchanged
sha256**. `offpeak_p` gains the n_sides argument (an internal function).

### Rd1 (wrong, confirmed; the proposed fix is wrong as written)

`gapup.txt`.
- **The shipped test is defeated by one echo.** On `m_echo_few`, every channel's shortest analysed gap is exactly
  50 ns, equal to the dark bound, so the excess is 0. On `m_echo100` (more than 100 echoes per channel) the
  spacing g[100] - g[0] is 0 and the statistic is NaN, so the comparison is False.
- **The red team's fix as written refuses every clean analysis.** It reads g[k] with k = max(10, 1e-3 N) but keeps
  the g[0] threshold (11.2). Under a clean edge g[k] - tau is about k mean spacings, so the excess is about k.
  Clean exponential edges with the upper end at tau: refused 1,599 / 4,000 at N = 2,000, 4,000 / 4,000 at
  N = 50,000 and 400 / 400 at N = 1e6. Its "null rate is no worse, since g[k] >= g[0]" has the inequality the
  wrong way round.
- **Council form.** j = max(`GAP_COUNT`, `GAP_SHARE` x N) (the downward count rule), M = max(100, j + 1).
  excess_j = (g[j] - upper) x M / (g[j+M] - g[j]). Refuse when excess_j > (j + 1) F.isf(alpha / K, 2(j + 1), 2M),
  with alpha = `GAP_UP_ALPHA` 1e-4. Under a memoryless edge (gaps = tau + Exp) the ratio
  (g[j] - tau) / (g[j+M] - g[j]) is Gamma(j+1) / Gamma(M), so the threshold is exact. At j = 0, M = 100 it reduces
  to today's 100((K / alpha)^(1/100) - 1).
  - Null, upper end = tau (the worst matched case): 0 / 4,000, 0 / 4,000 and 0 / 400 refused. The p-value of the
    F law is uniform: P(p < 0.05) = 0.037 / 0.053 / 0.055 and P(p < 0.01) = 0.007 / 0.013 / 0.013.
  - With afterpulse-like excess near the edge and the upper end above tau: 0 / 2,000 and 0 / 300.
  - Red team's data, K = 4: `m_echo100` and `m_echo_few` are **refused on all four channels** (excess 1.5e5
    against thresholds near 1,250). `m_clean100` passes (excess -5,600 to -11,500).
  - Power, the round-4 case (22 ns artifact, 50 ns detector, 1 Mcps): refused at N = 2,000, 50,000 and 1e6, as
    today. At 10 kcps and 50,000 gaps today's test refuses narrowly (15 against 11.2) and the council form passes
    (80 against 100). This loss is the price of robustness to j echoes. It is stated in the docstring.
- Echoes above `GAP_SHARE` (0.1%) of the analysed gaps still pass. `_check_time_units` tolerates up to 0.2%. This
  is covered by the setup cross-check below and stated as limit L1.

### Tw1 (wrong, confirmed)

`t3straddle.txt`, seed 4242 (the twin used 8080 and two more), 80 ps dtime, delay 4372 ps, jitter 40 ps:

| case | delay [95%] | status |
|---|---|---|
| T3 P = 12500.3 | **4365.10 [4364.60, 4365.60]** | ok |
| T3 P = 50040 | 4368.45 [4367.95, 4368.95] | ok |
| T3 P = 49999.37 | 4367.69 (misfit p 5e-51) | misfit |
| T3 P = 50000, and T2 | 4372.03 [4371.53, 4372.52] | ok |

The bias is -(delay / P) x (P mod res). At 12500.3 that is -0.35 x 20.3 = -7.1 ps, and -6.9 ps was measured. With
"ok" status and an interval that excludes the truth, this is wrong. It is small against a 40 ps jitter, but the
delay feeds window capture and the drift check. Round 8 tested only small delays.

### Other checks

- Rd2: `forward.py:438` sets `labels = (OUTSIDE_VALIDITY,)`, and `problems` is used only in the exception text.
  Confirmed.
- Rd3: `guardrails.CW_SOURCE_TYPES = ("cw", "spdc-cw")` and `config.SOURCE_TYPES = ("spdc-cw", "cw")`. Every
  declarable source type counts as CW, so a T3 declaration is never missing for a setup that names a type.
  Confirmed.
- Rd4: `condition_mismatches({'detector_temperature_K': 233}, {...: 298})` and the same for `temperature_C` return
  `[]`; `temperature_K` 233 / 298 is refused. Confirmed.
- Rd5: `grep forward_check_note dadl/cli.py` gives no hits. Confirmed.
- Tw2 and Tw3: `tagfile.TAGGER_SETTING_PREFIXES` records `MeasDesc_GlobalResolution` for both modes, and
  `SETTING_RTOL` compares it at 1e-4. Confirmed by code read. The twin's refusal texts match.
- Tw4: 12 sites print `Policy.allow(` hints; only `guardrails.py:257` adds the CLI form. README:563 says the
  serials are "optional" (Tw7). Confirmed.
- Pf2, Pf3: the patch was not re-timed. The council compared the reviewer's result pickles: `o_head_*` and
  `o_pat_*` are byte-equal at python 1e5 / 1e6 and rust 1e5 / 1e6 / 1e7, and patched-python 1e7 equals HEAD-rust
  1e7 (`cmp`). The performance reviewer's byte-identity evidence (pickles at 1e5 / 1e6 / 1e7 on
  both backends, replicate hashes) was read and is sufficient for code that only removes sorts and replaces
  `np.unique` on monotone input. `as_sorted` returns the caller's array instead of a copy. No analysis path
  writes into a tag array (grep of in-place writes in chsh, _twopoint, hom, forward, forward_check, coincidences,
  histogram), but a test must pin that.
- Pf5: `ci.yml` has `timeout-minutes` only on the coverage and propagation steps; the unit jobs `python`,
  `python-slow` and `python-minimum` have none. Confirmed.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| Rd1 | One echo turns off the upward gap test; S 2.029 for 1.980, unrefused | Y/Y/Y | ADOPT-FIX (council form, not the red team's as written) + ADOPT-LIMIT | Fix 2; L1 |
| Rd1b | No analysis-time cross-check of the dark dead time against the setup's declared `dead_time_ps`; D3 says "the truth lies between" | Y/Y/Y | ADOPT-FIX (label) | Fix 2 |
| Rd2 | `calibration_outside_validity` masks the gap refusals; reason dropped | Y/N/Y | ADOPT-FIX (own override name, problems recorded) | Fix 3 |
| Rd3 | `spdc-cw` counts as the T3 clock declaration; templates bypass G3 | Y/Y/Y | ADOPT-FIX + ADOPT-LIMIT (README sentence) | Fix 5; L2 |
| Rd4 | Unrecognised condition names recorded but never compared | Y/Y/Y | ADOPT-FIX (label, help text) | Fix 8 |
| Rd5 | Forward check skipped on T3, note not printed | Y/Y/Y | ADOPT-FIX | Fix 8 |
| Tw1 | T3 delay biased by pairs straddling a sync edge, status ok | Y/Y/Y | ADOPT-FIX (delay "9" -> "10") | Fix 4 |
| Tw2 | T2 dark artifacts refused on T3 analyses (shared key) | Y/Y/Y | ADOPT-FIX | Fix 6 |
| Tw3 | A CW T3 clock change refuses every artifact | Y/Y/Y | ADOPT-FIX (label under the explicit `cw` declaration) | Fix 6 |
| Tw4 | Refusals point CLI users to `Policy.allow`; override names and the 0.5 K rule undocumented | Y/Y/Y | ADOPT-FIX | Fix 7 |
| Tw5 | SNSPD refusal suggests a remedy no command delivers | -/Y/Y | ADOPT-FIX (wording) | Fix 7 |
| Tw6 | One dark PTU file cannot be calibrated as given | -/Y/Y | ADOPT-FIX | Fix 8 |
| Tw7 | README calls serials optional | -/Y/Y | ADOPT-FIX | Fix 7 |
| Tw8 | First-run path not signposted; NO_VALIDITY repeated 8x | -/Y/Y | ADOPT-FIX | Fix 8 |
| Pf1 | Off-peak dispersion half cannot refuse in CHSH; doubled cost | Y/Y/Y | ADOPT-FIX (sequential MC, counts-seeded) + ADOPT-LIMIT (residual power) | Fix 1; L4 |
| Pf1(b) | Null-statistic cache in `accidentals.py` | -/N/- | REJECT (superseded by Fix 1; it would re-pin accidentals for no gain) | none |
| Pf2 | `occupied_slots` uses `np.unique` on sorted input | -/Y/Y | ADOPT-FIX (from `perf.patch`) | Fix 9 |
| Pf3 | Repeated sorts and copies; 3x peak memory | -/Y/Y | ADOPT-FIX (from `perf.patch`) | Fix 9 |
| Pf4 | Real datasets do not fit in memory; not stated | Y/Y/Y | ADOPT-LIMIT + DEFER (M9, chunked analyses) | L3; O16 |
| Pf5 | chsh-spad shards near timeout; unit jobs without timeouts | -/Y/Y | ADOPT-FIX (timeouts; Fix 1 and Fix 9 cut the cost) | Fix 10 |
| Pf6 | Series dead-time constants recomputed per call | -/Y/Y | ADOPT-FIX (re-pin deadtime sha256, version stays "7", output-unchanged note) | Fix 9 |
| Pf7 | Peak fit uses finite-difference gradients | Y/N/- | DEFER (M10 performance pass; not bit-identical, needs the stated acceptance test) | none now |
| Pf8 | Repeated full-array passes | -/Y/- | ADOPT-FIX (single `np.quantile`, memoised grid counts: bit-identical) + DEFER (difference extraction to M9) | Fix 9 |

**Votes on Rd2.**
- P: the gap tests are physics checks of the dead time, not a validity window. Lumping them with a date check
  means the README's benign advice disables a physics refusal. Own name.
- M: votes N on a new override name. Each correction name is API surface, and a label carrying the problem text
  would be enough. Accepts recording the problems in `counts` in either case.
- U: the README tells users to pass `calibration_outside_validity` for a date edge case, and the user then cannot
  see that it also accepted a contradicted dead time. A separate name is the only way the user's reason refers to
  the right thing.

**Votes on Pf7.** P: the analytic gradient is right, and the fit is a gate hot spot. M: it changes the optimiser
path. Without the 200-replicate acceptance run (forward-check p within optimiser tolerance, flagged counts equal)
it is a behaviour change with no review. Defer to the M10 performance pass, with that run as its entry test.
U abstains.

**Seats on items marked "-".** P abstains on CI, packaging and README-only items (process).

## Limit wordings

- **L1** (README dead-time limits; `forward.check_dead_time_gaps` docstring): "The analysis tests the dark dead
  time against the analysed gaps in both directions. Each test tolerates a few unusual gaps, up to 0.1% of a
  channel's gaps or 10, whichever is more. A rarer double trigger (a cable reflection, a ringing discriminator)
  cannot hide the edge. A more frequent one makes the edge look shorter than the detector's dead time. If you
  declare the detector's dead time in the setup (`detectors.<name>.dead_time_ps`), DADL labels a dark dead time
  far below it. To be sure, look at the dark artifact's notes, which say when a single short gap set its bound."
- **L2** (README T3 paragraph, replacing the sentence at README:94-98, extending round 8's L2): "On T3 input the
  accidentals are read off-peak only when you declare `source.type: cw`. The template's `spdc-cw` does not count,
  because it is the default. With pulsed light, a repetition period that divides the off-peak step is not detected
  by the flatness check. The off-peak estimate is then right only if the first off-peak shift also falls on the
  comb. Otherwise V and |E| read low (conservative for the threshold), with no label."
- **L3** (README Known limits; CHANGELOG): "Analyses hold every channel's tags in memory: about 8 bytes per tag,
  plus up to the same again for working copies. A 10-minute CHSH run at 1 Mcps per channel (6e8 tags per channel,
  16 channel-settings) needs about 80-150 GB, more than a laptop has. Analyse shorter acquisitions or a subset
  until chunked analysis lands (milestone M9). Series dead-time calibration is not limited this way."

- **L4** (README, after round 8's L2; `_twopoint.offpeak_check` docstring): "The off-peak flatness check
  finds pulsed accidentals from how unevenly they fill the off-peak windows. With about one off-peak coincidence
  per window, as in a typical CHSH run, it reliably finds them only when they are about half or more of one
  acquisition's off-peak counts. In council round 9's test it refused 30 of 30 runs at one half, 8-17 of 30 at
  30%, and 0-2 of 30 at 20%. A smaller pulsed share passes unlabelled. The off-peak windows then miss part of the
  same-pulse accidentals, so V and |E| read low (conservative for the threshold; council round 9 red team). Longer acquisitions raise the counts per window and the power."

## Implementation order (adopted fixes)

1. **Fix 1 (Pf1; analysis/_twopoint.py, analysis/chsh.py, hom.py, fringe.py, README).**
   - `sequential_dispersion_pvalue(counts, n_sides, alpha=OFFPEAK_ALPHA, h=20, n_max_factor=4)` as in the
     council's `seqlib.seq_dispersion_pvalue`. Seed it with `np.random.SeedSequence(counts.astype(np.uint32))`.
     Return the p and the draws used. m >= `EXACT_BELOW` keeps the chi-square path.
   - `offpeak_p(side, offsets, n_sides)` uses it. `offpeak_check` passes `len(sides)`. Record the total draws in
     `counts["offpeak_draws"]`.
   - Bound the cost of a pulsed dataset. A flagged side runs to n_max (about 3 s). Once one side is below its
     threshold, evaluate the remaining sides only far enough to know whether they are also below it (stop at
     g = 3 exceedances). Report the flagged side's p. A refused 16-side analysis then costs about 3 s, not 50.
   - Remove `OFFPEAK_N_SIM`, and correct its comment's rationale ("reachable at 2 sides only").
   - `accidentals.py` is unchanged ("3", same sha256).
   - Tests:
     - The perf probe's periodic side (every 4th window) among 16 sides is refused (today it is not).
     - The `pulsed_chsh.py` case at 13157.9 ps with pulsed share 0.5 among 16 sides is refused in at least 28
       of 30 seeds (slow; council 30 / 30).
     - 20,000 iid Poisson(1.3) sides, n_sides 16: P(p <= 0.01) <= 0.012 and P(p <= 1e-3) <= 0.0015 (slow).
     - The same counts give the same p twice.
     - The 20 covering chsh-spad replicates are not flagged.
2. **Fix 2 (Rd1, Rd1b; forward.py, cli.py, README).**
   - Upward test: j = max(`GAP_COUNT`, `GAP_SHARE` x N), M = max(`GAP_EDGE_GAPS`, j + 1). excess_j as above.
     `gap_edge_excess(K, j, M)` = (j + 1) `scipy.stats.f.isf(GAP_UP_ALPHA / K, 2(j + 1), 2M)`. A zero spacing
     with g[j] > upper is infinite excess, and a NaN never passes (fail closed). Docstring: the exact law, and the
     power trade-off at low rates.
   - Above-gaps label and downward test unchanged.
   - Rd1b, label in `cmd_analyze`: when the setup declares `dead_time_ps` for a channel and the dark artifact's
     dead time (point, or bound) is below 0.5 x the declared value, add the label "the dark dead time of {name}
     ({x}) is below half the declared {y}: a short double trigger may have set it; check the dark artifact's notes".
     This uses the same rule as `cmd_calibrate_dark`'s warning, in one helper.
   - `DEAD_TIME_BOUND`: replace "(the truth lies between)" with "(the truth lies between if the bound holds; a
     bound set by a few short dark gaps may not)".
   - Tests:
     - The red team's `m_echo_few` and `m_echo100` data (or a synthetic equivalent: 1 us edge, 50 ns bound, 2 and
       500 echoes in 1e6 gaps) are refused.
     - 4,000 clean exponential channels with the upper end at tau are not refused (slow; one fast seed).
     - The round-4 case (22 ns on 50 ns, 1 Mcps, 50,000 gaps) is still refused.
     - An all-equal gap run (zero spacing) fails closed.
3. **Fix 3 (Rd2; forward.py, guardrails.py, propagate.py, README, CHANGELOG).**
   - New correction name `dead_time_contradicted`. The gap tests' refusal names it ("on the CLI,
     --allow dead_time_contradicted=REASON"). When it is allowed, the label is "dead time contradicted by the
     analysed gaps, accepted: {problems}".
   - `calibration_outside_validity` no longer covers the gap tests (behaviour change, O7).
   - `counts["dead_time_gap_problems"]` records the problem strings, whether refused or accepted.
   - The "no validity range" and "outside its validity" labels: see Fix 8 (one collapsed label).
   - Test: the red team's `m_darkonly` case with `--allow calibration_outside_validity` is still refused. With
     `--allow dead_time_contradicted` it carries the new label and the problem text in counts.
4. **Fix 4 (Tw1; calibration/delay.py, forward_check.py, adapters/ptu.py).**
   - For T3 input, compute each pair difference's sync-count difference dN (from the tags and `t3_period_ps`, or
     carried from the reader). Bin dt - dN x r on the dtime lattice, with r = P - res x round(P / res). Fit one
     shared (delay, jitter) with each dN group's model mean shifted by -dN x r. The groups' weights are the
     observed counts.
   - The predictive check sums the groups. `forward_check.pair_pvalue` and `pooled_peak` do the same when
     `PairData` carries the T3 period. The analysis path supplies it from the TagSet metadata.
   - `METHOD_VERSION` "9" -> "10". Regenerate `tests/method_versions.json`.
   - Tests:
     - 80 ps res, 12500.3 ps period, 4.4 ns delay: the delay interval contains the truth in at least 19 of 20
       seeds, with no misfit (slow; one fast seed).
     - 49999.37 ps: not a misfit.
     - 50000 ps and T2: unchanged to the last digit.
5. **Fix 5 (Rd3; guardrails.py, cli.py, README, CHANGELOG).**
   - On T3 input only an explicit `source.type: cw` counts as the clock declaration. `spdc-cw` gets
     `T3_SOURCE_REFUSAL` with the added sentence "the template's spdc-cw is the default and does not count:
     declare source.type: cw if the sync is a clock only".
   - L2 in the README.
   - Behaviour change (O7).
   - Test: the red team's `t3cw.py` case with `source_type="spdc-cw"` is refused, and with "cw" it reaches the
     off-peak check.
6. **Fix 6 (Tw2, Tw3; tagfile.py, propagate.py).**
   - Record `tagger_mode` ("T2" / "T3") with the tagger settings.
   - `condition_mismatches` compares `tagger_MeasDesc_GlobalResolution` only when both sides record the same
     mode. Artifacts without `tagger_mode` keep today's comparison.
   - When the dataset is T3 under an explicit `cw` declaration, a period mismatch on a dark or pairs artifact is a
     label, not a refusal: "T3 sync clock changed since calibration ({a} -> {b}); with CW light these calibrations
     do not depend on it".
   - README: "the T3 sync clock need not match for CW; record dark runs in T2 or T3".
   - Tests:
     - The twin's `hom/art_t2dark` with the T3 HOM run is accepted.
     - `hom_sync10` under `cw` gives the label. Under `spdc-cw` it is refused by Fix 5 first.
7. **Fix 7 (Tw4, Tw5, Tw7; one helper `guardrails.override_hint(analysis, correction)` used at every site
   listed in Tw4; README).**
   - Every hint reads "Override with Policy.allow(...); on the CLI, --allow NAME=REASON".
   - README "Using your own data": list the override names (`calibration_outside_validity`,
     `dead_time_contradicted`, `misfit_forward_model`, `misfit_calibration`), the 0.5 K temperature rule, and
     that recording `temperature_K` at calibration makes it required at analysis.
   - Tw5 wording as the twin suggests, naming `.json` and `.posterior.npy`.
   - README:563: "required for calibrations of recorded data".
   - Test: a parametrised check that every `DadlError` whose text contains "Policy.allow(" also contains
     "--allow".
8. **Fix 8 (Rd4, Rd5, Tw6, Tw8; propagate.py, cli.py, tagfile.py, examples).**
   - Rd4: a recorded condition outside the compared forms gets the label "condition {k} recorded, not compared".
     A lone `_min_` or `_max_` is compared one-sided. `--condition` help lists the compared forms.
   - Rd5: the CLI prints `forward_check_note` under the corrected value.
   - Tw6: `calibrate dark` / `pairs` / `deadtime` read only the settings they use. A single `.ptu` path is read
     as setting 0.
   - Tw8:
     - `init` prints both next steps.
     - The examples gain a commented `hardware:` block.
     - `NO_VALIDITY` collapses into one label listing the IDs.
9. **Fix 9 (Pf2, Pf3, Pf6, Pf8; from `r9-perf/perf.patch` minus the accidentals hunk).**
   - `coincidences.as_sorted` and `_dedupe_sorted`. Use `as_sorted` at the chsh, _twopoint, forward and hom sites.
   - The deadtime `_series` / `_survival_consts` hunk.
   - One `np.quantile` call in `_dead_time_estimate`.
   - Grid counts memoised per (pair, width, delays) within one `chsh()` call.
   - **Re-pin deadtime**: regenerate its sha256 in `tests/method_versions.json`, keep version "7", and add a
     `notes` field "code changed, outputs byte-identical (council r9 Pf6)". If the guard test cannot hold a note,
     put it in CHANGELOG.
   - Tests:
     - A user's tag arrays are not modified by `chsh`, `hom_from_tags` or `check_dead_time_gaps` (compare a copy).
     - `occupied_slots` equals `np.unique` on sorted and unsorted input.
     - The r9-perf result hashes (3 replicates per regime) equal HEAD's.
10. **Fix 10 (Pf5; ci.yml).**
    - `timeout-minutes: 75` on `python`, `python-slow` and `python-minimum` (the slowest was 41.8 min).
    - Keep 4 chsh-spad shards. If the first post-fix run shows any chsh-spad shard above 30 min, go to 6.
11. **Method versions and gates.**
    - **Pinned methods.**
      - delay "9" -> "10" (Fix 4).
      - deadtime: re-pin the sha256 only, version "7" (Fix 9, outputs unchanged).
      - accidentals "3", dark "7", dcr "2" and heralding "6" are untouched, with sha256 unchanged.
    - **Re-run `coverage delay`** (1,000 replicates, CI seed). The gate simulates T2 or plain tags, so the figures
      must equal round 8's exactly. Any change is a bug in Fix 4's T2 path.
    - **Re-run `coverage deadtime`.** Figures must equal round 8's exactly (Fix 9).
    - **Re-run the propagation gate**, all five regimes at 1,000 replicates. Expected:
      - corrected and control coverage identical to the record (no value changes);
      - zero refusals;
      - the off-peak p distribution on clean chsh-spad and hom-snspd re-measured (round 8 table rows: p<0.01 at
        most 2%);
      - hom-crosstalk flagged at least 90%.
      Record the wall time per shard (Fix 10).
    - **Also re-run:** `test_published_results.py` (report any new label or refusal from Fixes 1-3 and 5); the
      unit suite; the r7-test `mutate.py` harness, adding mutants for the F threshold (j and M swapped, alpha /K
      dropped) and for the sequential stop (h / L replaced with h / n_max).

## Owner items (O1-O14 as in r8-council.md; changes and additions)

| # | Item | Source | Status at 2e15eba | Council recommendation |
|---|---|---|---|---|
| O7 (extended) | Accept the round-9 behaviour changes in release notes: on T3 only `source.type: cw` declares the clock (templates' `spdc-cw` now refused); the gap tests have their own override `dead_time_contradicted` (`calibration_outside_validity` no longer covers them); the off-peak check can now refuse in CHSH; the upward gap test tolerates a few short gaps | r9 Rd1-3, Pf1 | New | Accept |
| O15 | Whether the example and `dadl init` setups should keep `source.type: spdc-cw`, or ship `cw` with a comment for T3 users. Changing a template changes what new users declare | r9 Rd3 | New | Keep `spdc-cw`; the T3 refusal text explains |
| O16 | Extend M9 ("chunked reading") to chunked analyses: coincidences, shifted counts, cross-correlation histograms and gap statistics over time blocks, plus the per-pair difference extraction (Pf8). Plan text and milestone scope | r9 Pf4, Pf8 | New | Adopt into M9; L3 until then |
| O17 | Pf7 analytic gradient in the M10 performance pass, with its acceptance run as the entry test | r9 Pf7 | New | Schedule in M10 |
| O18 | CI budget: the gate jobs are 477-494 of about 700-790 job-minutes, and the 20-job concurrency limit sets the wall time. Whether to buy concurrency, or run the propagation gate nightly and on release rather than on every push | r9 Pf5 | New | Owner's choice; the council does not decide it |
