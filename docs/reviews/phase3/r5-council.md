> Phase 3 review record, copied verbatim from the review's working directory. Edited only to remove home and scratch paths (now `<scratch>/`, or "a home or scratch path" where the round-10 ruling names them), the repository owner's account (`<owner>`) and the at-sign (one matrix product, now `sens^T cov sens`). Probe directories, scripts and data files cited below were not copied. Summary: [../2026-11-phase3-council.md](../2026-11-phase3-council.md).

# Council ruling, review round 5 (interval coverage, provenance, claims)

Seats: P = physicist, M = maintainer, U = user advocate. A finding is adopted with two of three votes.
Finding ids: I = interval statistician (r5-int.md, rows 1-3 as I1-I3), V = provenance archivist (r5-prov.md,
P1-P7 as V1-V7), C = claims auditor (r5-claims.md, rows 1-13 as C1-C13). Council probes:
`scratchpad/council/r5-council/` (`pend.py` the pending-afterpulse law and an instrumented copy of the
simulator's event loop; `check_law.py`, `ks_control.py` the law checks, outputs `law_*.txt`, `ksc_*.txt`;
`cov.py` coverage with the law patched into `_Run.probs`, outputs `cov_big.jsonl`, `cov_gate.jsonl`;
`ana.py` the paired comparison with the reviewer's records).

## Checked numerically

All runs on HEAD ef8bc0c with the repo venv, through the package's own functions, one BLAS thread.

### The pending-afterpulse gap law (I2): derived, checked, adopted

**Derivation.** In the simulator an afterpulse candidate is scheduled at t + tau_np + E, E ~ Exp(tau_a). If
the next record comes first (u < E), then by memorylessness the candidate is still pending at the new record,
r ~ Exp(tau_a) after it, independent of the gap that just ended. So the state at the start of a gap is the
number m of pending candidates, each Exp(tau_a) from the record. A pending one sits at v = r - tau_np in
u = gap - tau_np, density phi(v) = e^{-(v + tau_np)/tau_a} / tau_a on v > -tau_np:

- v <= -tau_p: blocked by the dead time and out of the paralyzable window at u = 0: no effect;
- -tau_p < v < 0: blocked, but it restarts the paralyzable stage: the next record follows q(u - v),
  q(t) = a S(t - tau_p) (the existing renewal-from-a-candidate law);
- v > 0: exactly the own-afterpulse case with E replaced by v (recorded with density phi(v) e^{-n tau_p} S(v),
  or a restart followed by q).

To first order (m <= 1, the own-afterpulse x pending product dropped):

    g(u) = g_ap(u) + pi * Delta1(u)
    Delta1(u) = K Delta_own(u) - K (e^{tau_p/tau_a} - 1) g0(u) + (K / tau_a) int_0^{tau_p} e^{s/tau_a} q(u + s) ds
    K = e^{-tau_np/tau_a},  Delta_own = the existing bracket multiplying p (module docstring)

Delta1 integrates to 0 (the two middle terms cancel in mass). The pending fraction per gap:

    pi_own = p int f_E(e) (1 - S(e + tau_p)) de      (own afterpulse not reached when the gap ends)
    pi     = pi_own / (1 - K pi_own / p)             (a pending one persists through the next gap: memoryless)

**Pending fraction against the instrumented simulator** (count of pending candidates at each record):

| case | gaps | pi, simulated | pi, formula | current guard formula p n tau_a/(1 + n tau_a) | m >= 2 |
|---|---|---|---|---|---|
| passive-quench, 6e6 cps (tau_p 40, tau_np 45 ns, p 2%, tau_a 30 ns, dark 50 kcps), 5 s | 2.35e7 | 0.00284 | **0.00283** | 0.00307 | 7e-7 |
| same, tau_a 90 ns, 3 s | 1.41e7 | 0.00847 | **0.00842** | 0.00705 | 2e-5 |
| stress, 1e7 cps (tau_p 20, tau_np 50 ns, p 5%, tau_a 200 ns), 2 s | 1.34e7 | 0.0671 | **0.0679** | **0.0333** | 1.9e-3 |

The current guard understates pending twofold at long tau_a (it misses the persistence through a gap).

**Gap law against simulation** (Pearson chi-square, groups of >= 20 expected; KS on the binned cdf):

| case | all gaps, current law | all gaps, corrected law | gaps starting with one pending: vs g_ap + Delta1 | vs current law |
|---|---|---|---|---|
| passive-quench 6e6, 5 s | chi2 1476/1442, p 0.26 (KS 0.32) | 1475/1442, p 0.27 (KS 0.76) | 790/822, p 0.78 (66,526 gaps) | 2506/848, p 3e-163 |
| tau_a 90 ns, 3 s | 1859/1767, p 0.06 (KS **6e-28**) | 1710/1766, p 0.83 (KS 0.19) | 836/815, p 0.30 | 20252/964, p 0 |
| stress 1e7, 2 s | **6006/1794, p 0** (KS 0) | 1885/1783, p 0.046 (KS 0.24) | 1178/1128, p 0.15 | 58472/1382, p 0 |

The conditional KS is small (p 0.004, 0.018, 0.001 at D 0.007, 0.004, 0.002); the KS control on m = 0 gaps
against g_ap passes (p 0.87, 0.10), so the residue is the dropped own-afterpulse x pending term (order p of
Delta1), not the derivation. The stress case's marginal p 0.046 at 1.3e7 gaps reflects the second-order part
(pi p + P(m >= 2), about 0.5% there). Note the regime itself: the current law passes the check with 2.3e7
gaps, so the bias it causes cannot be seen by the predictive check: it has to be modelled or guarded.

**Coverage with the law patched in** (`cov.py`; the harness reproduces the reviewer's records exactly with
the term off: replicate 0, 40344.9 ps):

| scan | current law | pending modelled |
|---|---|---|
| 0.5 s runs, afterpulse held at truth, 160 reps (paired with `r5-int/big_ap_0.5.jsonl`) | 88.1%, misses below/above 19/0, bias -163 +- 19 ps | **92.5%, misses 5/7, bias +18 +- 19 ps**; paired shift +180 ps (sd 17) |
| CI-seed gate replicates 0-399 (dark artifact, widening, margin 3; paired with `r5-int/ci.jsonl` A) | 94.75%, misses 17/4 | **96.0%, misses 9/7**; paired shift +181 ps (sd 44); 11 replicates change (8 gain, 3 lose) |

The bias is gone and the misses are balanced. The 0.5 s residual (92.5%, 12 of 160) matches the reviewer's
no-afterpulse control at 0.5 s (93.1%), so it is not afterpulsing: within 1.5 sd of 95%, and an M10 probe
item. The two switched gate replicates (165, 271; gains 3.64 and 3.48 with the law) become unions under
margin 5 and their held profiles cover (r5-int), so the expected gate with both fixes is about 96.5% on these
400: inside the band 92.4-97.2%.

Cost: the probe's term is naive (a Python loop of survival calls per NLL) and makes a fit 5x slower; the
implementation must reuse what `_afterpulse` already computes (below).

### Other checks

- **I1 (branch margin): confirmed.** False switches in `r5-int/ci.jsonl` and `s4242.jsonl`: 9 of 2,000, gains
  3.09-4.97, all at tau_p 28.6-32.7 ns. The Cox bound: max over D of Phi(-(M + D)/sqrt(2D)) is at D = M,
  Phi(-sqrt(2M)): 7.2e-3 at M = 3, **7.8e-4 at M = 5**. Observed 0.45% against 0.51% predicted at D = 5.6.
- **V1: confirmed** by code read: `deadtime.py:672` inputs are the runs' raw-tags only; `previous-artifact`
  appears only in the schema; `Calibration` carries no artifact ID; dcr copies its inputs' inputs
  (`dcr.py:138`), the pattern to follow.
- **V2, V3, V5: confirmed** by code read: `__init__.py:12` "0.1.1"; `forward.py:362` seed and n_draws are
  arguments not recorded by `as_record`; `propagate.py:216` compares only `cal.conditions`.
- **C1, C4, C5, C6, C7, C8: confirmed** (`CHANGELOG.md:115, :140`; `design.md:621`; `coverage.py:505-512`;
  `ci.yml:11`; `CONTRIBUTING.md:49-50` against `MODULES`; `simsetup.py:73`). C6 is right on GitHub's
  semantics: one running and one pending run per group; a newer pending run replaces the older.

## Rulings

| id | finding | votes (P/M/U) | ruling | action |
|---|---|---|---|---|
| I1 | `BRANCH_MARGIN` 3 lets about 0.5% of below-maximum scans switch to tau_p about 30 ns, every one a miss | Y/Y/Y | ADOPT-FIX | Fix 1: margin 5, docstring with the Cox bound; test on CI-seed replicate 477 |
| I2 | Pending afterpulses unmodelled: tau_p bias -170 ps that does not shrink (88% at 0.5 s); guard fixed-fraction and its note understates tenfold | Y/Y/Y (model now); Y/Y/Y (residual guard) | ADOPT-FIX (model) | Fix 2: the first-order law above (derived, checked against simulation at both requested parameter sets, coverage measured: bias +18 +- 19 ps, misses balanced). Round 1's M10 deferral is withdrawn: the evidence shows a bias that grows relative to the sd with data and is invisible to the predictive check. Fix 3: a data-scaled guard on the unmodelled remainder |
| I3 | Coverage history reads as drift; record the miss split and switch count | Y/Y/Y | ADOPT-FIX (docs) | Fix 7: the re-run's entry records coverage, below/above misses, switches and unions |
| V1 | `deadtime-series` does not cite its dark artifact | Y/Y/Y | ADOPT-FIX | Fix 4 |
| V2 | Package still says 0.1.1; `valid_for.from` rule changed under the same version | Y/Y/Y (limit); version: owner | ADOPT-LIMIT + owner item | Wording L1; the dev version is an owner release decision (owner item 1) |
| V3 | Corrected records lack the Monte Carlo seed and draw count | Y/Y/Y | ADOPT-FIX | Fix 5 |
| V4 | Method versions never checked at use | Y/Y/N (defer); N/N/Y (now) | DEFER (M5) | Nothing reads `deadtime-series` until M5 (`SERIES_UNUSED`); add the minimum-version table (kind, method) -> version with a test as an M5 entry item, before M5 reads series artifacts |
| V5 | `valid_for.conditions` validates but is ignored | Y/Y/Y | ADOPT-FIX | Fix 6 (enforce: the package never writes it, so no existing artifact changes behaviour) |
| V6 | "readers from 0.1.x refuse" is false | Y/Y/Y | ADOPT-LIMIT | Wording L2 (schema `kind.description`; the plan's D4 text is the owner's: owner item 3) |
| V7 | Widened schema keeps the v0.1.0 `$id` | Y/Y/Y (freeze sentence); `$id`: owner | ADOPT-LIMIT + owner item | Wording L3 now; the `$id` target names a release (owner item 2) |
| C1 | Current passive-quench coverage figure recorded nowhere | Y/Y/Y | ADOPT-FIX (docs) | Fix 7: record the method-5 re-run figures; mark the M4 line "(before the review rounds)". Do not record the round-4 93.2% as current: it is superseded by this round |
| C2 | Gap-check floor quoted only for a favourable case | Y/Y/Y | ADOPT-LIMIT | Wording L4 (README:382-387, forward.py:253-258) |
| C3 | README "What works today (0.1.0)" lists unreleased features | Y/Y/Y | ADOPT-LIMIT | Wording L5 |
| C4 | design.md phase 3 "Not started" | Y/Y/Y | ADOPT-LIMIT | Wording L6 |
| C5 | Coverage `--jobs` help stale; `run()` families omit deadtime | Y/Y/Y | ADOPT-FIX (docs) | As suggested (coverage.py:630, :508) |
| C6 | CI comment "pushes to main all run" is false | Y/Y/Y | ADOPT-FIX (docs) | As suggested (ci.yml:11) |
| C7 | CONTRIBUTING calls every calibration module pinned | Y/Y/Y | ADOPT-FIX (docs) | As suggested (CONTRIBUTING:49-50) |
| C8 | `--rates` default holds for four settings only | Y/Y/Y | ADOPT-FIX (docs) | As suggested (cli.py:1016, CHANGELOG:137-138) |
| C9 | Counts-key table wrong about raw/corrected keys | Y/Y/Y | ADOPT-FIX (docs) | As suggested (README:389-399); add `monte_carlo` (Fix 5) to the corrected list |
| C10 | README says fringe has a `--calibrations` path | Y/Y/Y | ADOPT-FIX (docs) | As suggested (README:60-61) |
| C11 | Plan and site stale (a-d) | Y/Y/Y (a-c); (d) at M10 | ADOPT-FIX (a-c) + DEFER (d, M10) | (d) the site describes 0.1.1 as released; update it with the release |
| C12 | Docstrings/help miss refusals and limits | Y/Y/Y | ADOPT-FIX (docs) | As suggested (errors.py:24-26, deadtime.py docstring, `calibrate deadtime` epilog, cli.py:6) |
| C13 | design.md says afterpulsing is fitted jointly | Y/Y/Y | ADOPT-FIX (docs) | As suggested, plus "pending afterpulses modelled to first order" |

Votes on I2 "model now": P: the law is exact to first order and validated; a guard alone leaves a bias that
grows as sqrt(gaps). M: the change is local to `_Run` (one term), the method version bumps anyway for I1,
and the gate re-run is needed anyway; condition: the overhead target below. U: users with long runs get an
interval that misses low without a note today; a model beats a note.

### Limit wordings

- **L1** (CHANGELOG Unreleased, README "Calibration artifacts"): "Artifacts written by 0.1.x start their
  validity at creation (`valid_for.from` = created); this branch starts it at acquisition. For data taken
  between an old artifact's acquisition and its creation, pass the validity override or rewrite the artifact
  with `--supersedes`. Until the next release the package still reports version 0.1.1."
- **L2** (schema `kind.description`, both copies): "0.1.x readers load a `deadtime-series` artifact but nothing
  in them uses it; validators against the 0.1.0 schema refuse it."
- **L3** (schema top-level description, both copies): "Frozen: later versions may add optional fields (and
  the kind `deadtime-series`, decision D4)."
- **L4** (README:382-387 and forward.py:253-258, same words): after the 7% figure, "the floor widens as the
  dark run has fewer gaps: the bundled chsh-spad dark run (20 s at 500 cps) puts it 40-50% below the point,
  and the downward test then catches only gross errors".
- **L5** (README): heading "What works today (0.1.1, and unreleased on main)"; tag the deadtime,
  forward-model and `--calibrations` items "(unreleased)"; status "phase 1 done (unreleased), phase 2 done,
  phase 3 in progress (M0-M4, unreleased)".
- **L6** (design.md:621): "In progress (unreleased): M0-M4 (calibration sets, forward model, corrected
  HOM/fringe/CHSH, `propagation` gate, series dead time); M5-M10 planned in phase3-plan.md".
- **L7** (deadtime.py docstring "Afterpulsing", README dead-time section): "An afterpulse still pending when
  the next gap starts (from the previous avalanche, Exp(tau_a) after the record by memorylessness, a fraction
  pi of gaps recorded as `afterpulse_pending_max`) is modelled to first order: one pending candidate, without
  its product with the record's own afterpulse. The remainder, about pi (p + pi) of gaps, biases tau_p low by
  an amount that does not shrink with data (about 0.1 sd per unit of remainder x sqrt(gaps)); above 2 a note
  says so. Checked against the simulator at 6e6 cps (p 2%, tau_a 30 and 90 ns) and 1e7 cps (p 5%,
  tau_a 200 ns; pi 6.8%, where the remainder shows at 1.3e7 gaps, chi-square p 0.05)."

## Implementation order (adopted fixes)

1. **Fix 1 (I1, deadtime.py).** `BRANCH_MARGIN = 5.0`. Docstring "Two rates per gap law": replace "wins by at
   most 2.1 on 256 simulated scans" with "wins by more than 3 in about 0.5% of scans (9 of 2,000, round 5);
   with margin 5 the false-switch rate is at most about 8e-4 for any separation (Cox's normal approximation
   for non-nested fits: Phi(-sqrt(2M)))"; past the maximum "the margin is passed in 28 of 48 and the union
   covers the rest". Keep the comment "(2 ln LR = 10)". Test (slow-marked if > 30 s): CI-seed passive-quench
   replicate 477 (`_children(20261111, ...)`, gain 4.17 at margin 3) gives a union note and an interval that
   contains 40 ns.
2. **Fix 2 (I2, deadtime.py).** Add the pending term to `_Run.probs` when p > 0:
   `prob += pi(n, tau) * Delta1 * width`, with Delta1 and pi as in "Checked numerically"; tau_np enters as the
   fit's shortest gap m (pass it to `_Run` or `_Fit`). Efficiency (required; the naive probe is 5x): in
   `_afterpulse`, pi_own/p is sum(f) / h-weighted from the `survival(w + tau)` it already evaluates (extend
   `w` to 40 tau_a); the J term is a correlation of the already-computed `q_cell(lags)` with the weights
   e^{s/tau_a} on ceil(tau_p/h) cells (one `fftconvolve`, lags extended by tau_p), and the same matrix
   product as now beyond the uniform grid. Target: fit time at most 1.5x HEAD's on the passive-quench
   replicate (measure with `r4-council/prof.py`). Record `conditions["afterpulse_pending_max"]` as the new pi
   (max over runs). Tests: (a) law: simulate tau_p 20 ns, tau_np 50 ns, 1e7 cps, p 5%, tau_a 200 ns for 0.2 s;
   the binned Pearson chi-square (groups >= 20) against the corrected `probs` at the true n and tau_p has
   p > 1e-3, and against the law without the term p < 1e-6; (b) pi at those parameters is within 3% of 0.0671
   (the simulated value, `law_stress.txt`) and at the passive-quench 6e6 cps parameters within 3% of 0.00284.
3. **Fix 3 (I2 guard).** Replace `PENDING_MAX` by `UNMODELLED_Z = 2.0`: remainder r = max over runs of
   pi (p + pi); `conditions["afterpulse_unmodelled"] = r`; when r * sqrt(total gaps) > 2 the note "about
   {r:.2%} of gaps carry a second pending afterpulse or one beside the record's own, which the gap law does
   not model: tau_p is biased low by an amount that does not shrink with more data (about {0.1 r sqrt(gaps):.1f}
   sd here, twice that for tau_a beyond twice tau_p); scan at lower rates or with faster afterpulsing". Test:
   the note fires at the stress parameters with 1e6 gaps and not at the passive-quench regime.
4. **Fix 4 (V1).** `calibrate_deadtime_series`: inputs += the dark artifact's inputs (`afterpulse.inputs`,
   deduplicated, as `dcr.py:138`); `conditions` add `afterpulse_prob_sd_held`, `afterpulse_tau_sd_held` (the
   draws' sds). CLI `calibrate deadtime`: also add `(dark_artifact_id, "previous-artifact", str(channel))`
   from the loaded manifest (`artifact_id`). Test: two dark artifacts from different dark runs give
   different `inputs`; the CLI artifact lists the dark ID with role `previous-artifact`.
5. **Fix 5 (V3, forward.py / result.py).** `counts["monte_carlo"] = {"seed": seed, "n_draws": n_draws}` in
   `corrected_visibility` and `corrected_chsh`. Test: the record has both, and seeds 0 and 1 differ there.
6. **Fix 6 (V5, propagate.py).** Run `condition_mismatches` on `valid_for.conditions` too (merged with the
   top-level conditions; on a key in both, both must match). Test: the reviewer's re-minted dark artifact
   (`valid_for.conditions` 300 K, 25 V) is refused for a dataset at 250 K, 5 V, and labelled under the
   override.
7. **Method version and gates.** deadtime `METHOD_VERSION = "5"` (Fixes 1-4 together; one bump), regenerate
   `python/tests/method_versions.json`. Re-run `python -m dadl.calibration.coverage deadtime --replicates
   1000` on the CI seed for **passive-quench** (expected about 96%, band 92.4-97.2%) and **series-spad** (p = 0:
   Fixes 1-3 do not act; expect 95.5% unchanged, a confirmation). Record in the CHANGELOG round-5 entry:
   coverage, below/above misses, switched and union counts, misfit share, tau_np coverage (I3, C1); mark the
   M4 line's figures "(before the review rounds)". Propagation and published-results gates: no re-run
   (`deadtime-series` is not used by HOM/fringe/CHSH; Fixes 5-6 change records and refusals only, and no
   gate artifact has `valid_for.conditions`); the unit suite must pass. If passive-quench lands above 97.2%,
   do not widen the band: report to the council (the afterpulse widening may then be conservative).
8. **Docs.** L1-L7; C5-C13 as ruled; plan "Smaller items" at M10: remove "model the pending candidate", add
   "a 0.5 s passive-quench probe (0.5 s coverage 92.5%, matching the no-afterpulse control's 93.1%: check the
   profile's calibration at 1e6 gaps)" and "the second-order pending term if a regime needs it"; M5 entry:
   the minimum-method-version table (V4).

## Owner items

1. **Version (V2).** Moving `__version__` (and the core's) to a dev version such as 0.2.0.dev0 is a release
   decision; the council recommends it before more artifacts are written on this branch, since method
   versions alone cannot tell 0.1.1 from branch artifacts on the `valid_for.from` rule.
2. **Schema `$id` (V7).** Point `$id` at the release that widens the enum; it names a release, so the owner
   picks it with item 1.
3. **Decision D4 text (V6).** The plan's D4 says 0.1.x readers refuse `deadtime-series`; the owner's decision
   text should read as L2.
4. **Withdrawn deferral.** Round 1 deferred pending-afterpulse modelling to M10 with the owner's plan; this
   round adopts it in M4 (method 5). The owner should know the plan table changes.
