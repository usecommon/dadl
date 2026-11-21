> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Bayesian computation (M0 benchmark, posterior draws for joint Monte Carlo, M5 readiness)

Reviewer: Bayesian computation (round 6). HEAD 7bcc32b, repo venv, one BLAS thread. Probes in
`scratchpad/council/r6-bayes/`:
- `sbc_dark.py` + `ana.py`: rank checks of `calibrate_dark` draws in five regimes, outputs `sbc_*.jsonl`.
- `sbc_series.py` + `ana_ser.py`: rank checks of `calibrate_deadtime_series` draws, outputs `ser_*.jsonl`.
- `sbc_delay.py`: rank checks of `calibrate_delay` draws, output `sbc_delay.json`.
- `prior_box.py`, `prior_box2.py`: the flat prior in the fit coordinates near the detection threshold.
- `union.py`: draws of union (two-branch) intervals. `cap.py`: draws at the tau_p <= tau_np cap.

## Simulation-based calibration (rank statistics) of the stored draws

Each replicate draws a truth from a prior, simulates data with the package's simulator, runs the
calibration with the default 2,000 draws, and records the rank of the truth among the stored draws
(ties broken at random). If the draws are a calibrated posterior, the ranks are uniform. Selecting
replicates on the data (for example "afterpulsing detected") keeps them uniform, so the selection is
valid. `joint_d2` is the rank of the truth's Mahalanobis distance in (log r, logit p, log tau), a
joint check.

| calibration, regime (truth prior) | reps | parameter: mean rank, central 95% of draws covers, KS p |
|---|---|---|
| dark, 3,000 gaps, p 3-10% (SIR path) | 400 | rate 0.520, 93.8%, 0.05; p 0.525, 95.5%, 0.03; tau_a 0.484, 93.8%, 0.34; dead time 0.512, 97.5%, 0.18; joint 0.523, 95.5%, 0.29 |
| dark, 30,000 gaps, p 0.6-1.4% (SIR path, n p < 500) | 400 | rate 0.527, 91.7%, 0.14; p 0.538, 95.0%, 0.03; tau_a 0.520, 95.5%, 0.16; dead time 0.493, 93.0%, 0.45; joint 0.493, 94.0%, 0.65 |
| dark, 20 s, p 10-25% (Laplace path) | 400 | rate 0.487, 96.5%, 0.32; p 0.521, 95.0%, 0.21; tau_a 0.506, 94.0%, 0.63; dead time 0.500, 94.0%, 0.78; joint 0.489, 96.5%, 0.21 |
| dark, 600 gaps, p 8-20%, **p uniform** | 1,200 | p 0.543, 94.6%, **2e-6** |
| dark, 600 gaps, p 8-20%, **logit p uniform** | 1,200 | p 0.516, 95.7%, 0.14 |
| delay-jitter, ~400 coincidences, delay 1-5 ns, jitters 100-400 ps | 600 | delay 0.501, 95.8%, 0.89; jitter 0.527, 93.8%, 0.12 |
| deadtime-series, series-spad, tau_p U(8, 36) ns | 400 | tau_p 0.500, 95.0%, 0.89 (all identified, no unions) |
| deadtime-series, passive-quench, tau_p U(30, 44) ns, dark artifact from a 1 s run | 363 | tau_p 0.524, 94.2%, 0.08 (9 unions) |

Reading: the draws that M1-M4 actually push through the joint Monte Carlo (the dark `dead_time_ps`
tick-law draws and the delay/jitter Laplace draws, `forward._channel_dead_times`, `forward._capture`)
are calibrated, and so are the dark SIR draws and the series tau_p draws. The one skew (p, mean rank
0.543) goes away when the truth prior is uniform in logit p: the stored draws are posteriors under a
flat prior in the fit coordinates, exactly as `dark.py:264` says. The importance-sampling ESS was
3,800-3,990 of 4,000 in every regime above, and 2,100-3,100 at the detection threshold
(`prior_box2.py`, slow afterpulsing).

## Findings

| # | Severity | Finding | Evidence | Where | Suggested fix |
|---|---|---|---|---|---|
| 1 | questionable | **M5's named first users have no likelihood to sample.** The plan sends NUTS to "the M2/M3 joint fits at low counts", but M2/M3 has no fit: the corrected value is a closed-form transform of Jeffreys Gamma count draws times calibration draws, which is already exact Monte Carlo of a conjugate posterior. NUTS adds nothing there unless M5 first defines a joint likelihood (for example the whole dip or fringe histogram with the IRF, or calibration data and counts in one likelihood). The M0 record's own limit, "the M2 HOM joint model should be benchmarked again before M5 lands", is tracked nowhere (not in M5 or in the smaller-items table). | `forward.py:18` "the model is exactly identified and the fit is closed form"; no `optimize`/`minimize` in `forward.py`, `analysis/chsh.py`, `analysis/hom.py`; `2026-11-ppl-benchmark.md` Limits; `phase3-plan.md:195-211` | `docs/phase3-plan.md` M5 | M5 entry item: name the model NUTS will sample (likelihood, parameters, priors) and rerun `bench/ppl` on it, or say that M5 covers calibration fits only. Add the benchmark rerun to the M5 entry items. |
| 2 | questionable | **Multimodality cannot be detected by the planned convergence check.** M5 starts every chain at the MLE (benchmark finding 3) and calls R-hat > 1.01 a misfit. R-hat from chains that all start in one mode cannot see a second mode. DADL has two known bimodal posteriors: the delay model's broad-peak mode (R-hat 1.55 with dispersed starts in M0) and the series dead time's two branches per rate. In the passive-quench rank check 9 of 363 scans (2.5%) were unions, and in replicate 158 the branch intervals were disjoint (26.5-31.8 and 35.1-39.6 ns). A NUTS fit started at the low-branch MLE would report only that branch, while the Laplace path reports the union. | `bench/ppl/common.py:16-17`, benchmark finding 3; `phase3-plan.md:202-203`; `union.py` output | M5 `nuts.py` design | Start chains from overdispersed points (for example the Laplace covariance times 4 to 9, plus one chain per known branch or mode), and keep the explicit branch enumeration and union outside the sampler. Record the starting rule in `checks`. |
| 3 | questionable | **Hard boundaries and support-defining parameters will mark low-count fits misfit under "any divergence".** (a) The dead time sets the support (tau <= shortest gap m). The likelihood is discontinuous in it, so HMC cannot sample it; today it is plugged in and drawn from the exact tick law, and that should stay. (b) tau_p <= tau_np: at full paralysis (passive quench, tau_p = tau_np = 45 ns) a median of 48% (max 51%) of the stored tau_p draws sit exactly at tau_np (`cap.py`, 16 reps). This is a clipped point mass, not a truncated posterior, and a sampler on a box constraint there diverges. (c) The fit's boxes (`LOGIT_P_RANGE`, the tau_a search range) act as hard walls. | `deadtime.py:699` `np.minimum(tau_draws, np_draws)`; `dark.py` `_dead_time`; `phase3-plan.md:203` | M5 plan and `nuts.py` | Write into M5: the dead time stays conditioned on m with exact draws (no HMC); reparametrize tau_p = tau_np * sigmoid(z) and other boxes as smooth transforms with proper priors; record divergences but do not make the "any divergence" rule apply before reparametrization is tested on the passive-quench and series-spad regimes. |
| 4 | questionable | **The stored draws imply a prior that M5 must state, and a flat prior in logit p is improper.** The rank checks show the dark draws are posteriors under flat priors in (log r, logit p, log tau). With a uniform-in-p truth prior the p ranks skew (0.543, KS p 2e-6, 1,200 reps); with a logit-uniform one they do not (0.516, p 0.14). As p -> 0 the likelihood tends to a constant (the no-afterpulse model) for every tau_a, so a flat prior on logit p and log tau_a is improper. The Laplace/SIR draws hide this because the Gaussian proposal never reaches the plateau. Bounded by the fit's box, the plateau holds 0-1% of the mass at 2 ln LR 19-28, and the box posterior's 2.5% quantile of p sits below the stored draws' (0.012 against 0.020 at 2 ln LR 18.9). NUTS with "flat priors" would sample whatever box it is given. The benchmark used Normal(MLE, 10): data-centred, and not DADL's prior. | `ana.py` on `sbc_tiny1200.jsonl`, `sbc_tinylogit.jsonl`; `prior_box2.py` output; `run_numpyro.py:24-26` | M5 plan; `dark.py:264` docstring | Declare proper priors per calibration before M5 (for example Beta on p, log-uniform tau_a over the search range), state them in the artifact (`conditions` or a `prior` note), and bump the method version where they move the draws: about 0.1 sd in p at 10-60 afterpulses. |
| 5 | questionable | **Series tau_p draws lose their correlation with the dark artifact's afterpulse draws.** The afterpulse uncertainty enters tau_p by a linear sensitivity, as an independent normal term added to a split normal (`deadtime.py:663-669`), not as a function of the dark artifact's draw rows. `CalibrationSet.draws` then draws the dark and series artifacts with independent permutations (`propagate.py:383-386`). In passive-quench the afterpulse share of tau_p's sd is 0.26-0.41 (median 0.33, 363 reps), so any M5 model that uses both the dark p and tau_a and the series tau_p drops a correlation of about 0.3. Latent today (`SERIES_UNUSED`). | `ser_passive-quench.jsonl` (`ap_shift` notes against the interval width); `deadtime.py:636-669`; `propagate.py:376-399` | `deadtime.py`, `propagate.py` (M5 entry) | Build the series draws from the dark artifact's rows: tau_p_i = stat_i + sens . (ap_i - mean). In `draws()`, share one row index between an artifact and the artifact it derives from (the series artifact already carries the dark run's raw-tags inputs and, from the CLI, a `previous-artifact` input). Add a test: corr(series tau_p, dark p) across joint draws matches sens . cov / sd. |
| 6 | minor | **Union intervals are drawn as one split normal, which puts mass between disjoint branches.** In replicate 158 (branches 26.5-31.8 and 35.1-39.6 ns), 18% of the stored draws fall in the gap where neither fit is plausible. The draws still cover (ranks are fine), but they are not a posterior. Latent until M5 reads series draws. | `union.py` output; `deadtime.py:663-669` | `deadtime.py` | For a union, draw from a two-component mixture (one split normal per branch, weights from exp(-NLL) of each branch fit, or equal weights as the conservative choice). |
| 7 | minor | **"Stored draws gated" is assigned to M1 but was not done, and M5 adds no draw-level test.** The coverage harness checks intervals only (`n_draws=10`). Review item 7.17 is still "Open". M5 tests compare intervals (Laplace against NUTS), but a sampler can match 95% intervals and get the shape or the correlations wrong, and the joint Monte Carlo uses the whole distribution. | `phase3-plan.md:336`; `2026-11-review-cycles.md:184`; `coverage.py:256, 461` | plan, M5 tests | Add a rank-statistic test (as `sbc_dark.py`: a truth prior, 400 replicates, uniformity at 10 bins, Bonferroni over parameters) to the coverage fast tier for each calibration whose draws are propagated, and run it with `--sampler nuts` in M5. Close 7.17 with it, or move the plan row from M1 to M5. |
| 8 | minor | **The benchmark record's speed summary is selective.** "Speed where it matters" cites first-run time and the dark model. On the delay model in steady state, Stan was 1.6-3x faster than NumPyro on every OS, and NumPyro had the lowest ESS/s on macOS for both models (ranges overlap). The decision rests on install and first-run time, which the data support. | `2026-11-ppl-benchmark.md` Results table, "Why, in one line each" | benchmark record | Add one clause: "in steady state Stan was fastest on the delay model and the libraries overlap on macOS; NumPyro was chosen on install and first run". |

## Checked and found sound

- Benchmark fairness: same models and priors, all chains started at the MLE, 4 parallel chains in
  every library (`chain_method="parallel"` with 4 host devices, `cores=CHAINS`, `parallel_chains`),
  default target acceptance 0.8 everywhere, fresh binary-only venvs, second-run timing. The CmdStan
  draw reshape (`run_stan.py:46-47`) is chain-major and correct. The decision (NumPyro as an optional
  extra, Laplace default) follows from the install data.
- The quantile-shift tie (worst 0.07 half-widths) is what Monte Carlo noise gives at tail ESS of about
  1,000 (about 0.04 per quantile, max of about 0.1 over all runs). The 0.2 threshold cannot rank the
  libraries, and the record does not claim it does.
- The dark draws (Laplace and SIR), the dark dead-time tick-law draws, the delay/jitter Laplace draws
  and the series tau_p split-normal draws pass rank checks in all eight regimes above, including
  joint (Mahalanobis) ranks. The posterior correlations among (log r, logit p, log tau) are near 0
  (|mean| < 0.04) in these regimes, so the joint check has little power against lost correlation
  there.
- `CalibrationSet.draws`: one row index per artifact keeps within-artifact correlation; seeding per
  artifact ID gives independence across artifacts; it samples with replacement when n > 2,000
  (`N_MC` = 4,000), which is valid. Stand-in columns are returned as bounds and never drawn. The
  forward model uses bounds at the bound (D3).
- tau_p against tau_np correlation in the series draws: |corr| < 0.08. The tau_np interval is about
  2 ps against a tau_p interval of about 3 ns, so treating them as nearly independent is harmless.
- The series profile intervals cover 95.0% (series-spad) and 93.4% (passive-quench, 363 reps) under
  a prior over tau_p instead of a single truth.
