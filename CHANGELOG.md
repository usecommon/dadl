# Changelog

All notable changes. DADL is pre-alpha: in 0.x anything may change between
releases. The `calibration/1` artifact schema is frozen from 0.1.0 (new
optional fields may be added; changes to existing ones need `calibration/2`)
and is versioned separately from the packages.

## Unreleased

Artifacts written by 0.1.x start their validity at creation (`valid_for.from` = created); this
branch starts it at acquisition. For data taken between an old artifact's acquisition and its
creation, pass the validity override or rewrite the artifact with `--supersedes`. Until the next
release the package still reports version 0.1.1.

### Phase 3, M0: inference library benchmark

- `bench/ppl/` and the `PPL benchmark` workflow measure NumPyro, PyMC and Stan on Windows,
  Linux and macOS (binary-only install, size, compile and sampling speed, agreement with DADL's
  intervals). NumPyro is chosen as the optional MCMC backend (decision D2); see
  `docs/reviews/2026-11-ppl-benchmark.md`. No package code changed.

### Phase 3, M1: calibration sets

- `dadl.propagate.CalibrationSet` loads verified artifacts for an analysis and draws their
  posteriors jointly: columns of one artifact share a sample index; artifacts are drawn
  independently, keyed on their IDs. Stand-in columns are returned as bounds, never as draws
  (decision D3).
- Refusals before any calibration is used: misfit artifacts (override `misfit_calibration`);
  artifacts outside their `valid_for` range or whose recorded conditions do not cover the
  dataset, with `tagger_*` settings now compared at use (new `StaleCalibration`, override
  `calibration_outside_validity`); and artifacts calibrated on the analysed data (no override).
  Overrides label the result "policy modified".

### Phase 3, M2: forward model and corrected HOM visibility

- `dadl.forward.corrected_visibility`: the corrected two-point visibility fitted under a model of
  what the detectors recorded (live fraction of the channel pair from propagated dead times;
  accidentals from each acquisition's own off-peak bins), with joint Monte Carlo over counts and
  calibration draws. The interval is the 95% quantiles; the statistical and calibration parts
  are reported separately (`Estimate.stderr_statistical`, `stderr_calibration`).
- `hom_from_tags(..., calibrations=, channels=)` uses it; the subtraction path stays when no
  calibrations are given, and asking for both is refused. Results list `corrections` in the
  order applied and carry the artifact IDs. A dead time known only as an upper limit is
  corrected at the bound, with the value at zero alongside (decision D3).
- Dead-time correction is not in the default HOM rule (decision D10): it needs
  `Policy.allow("hom_visibility", "dead_time", reason)` and is labelled "policy modified".
- `dadl.analysis.fringe_from_tags`: Bell-type fringe visibility from a fringe-maximum and a
  fringe-minimum acquisition, raw under the default policy and corrected by the same forward model
  when the policy allows it. Acquisitions of unequal length (more than 1%) are refused.
- A dead time known only as an upper limit so loose that, at the recorded rate, it would be a
  dead fraction above 20% is refused: the range "between zero and the bound" carries no
  information (100 cps of SNSPD dark counts bound a 30 ns dead time only below about 10 us).
- Window capture: the fixed grid splits more pairs at a HOM reference delay (photons +-d apart)
  than at the dip, so V read about 1 sigma low and the corrected intervals undercovered (90.5%
  and 91.5% at 200 replicates). With a delay-jitter artifact for the channel pair and
  `hom_from_tags(..., reference_delay_ps=)`, the forward model now includes each acquisition's
  capture (closed form, `dadl.forward.grid_capture`), with delay and jitter drawn from the
  artifact. It is a new correction, `window_capture`, not in the default HOM rule (like dead
  time); without the artifact the corrected value is labelled.
- `dadl analyze SETUP --tags F --calibrations DIR [--allow dead_time=REASON] [--acquired T]
  [--condition K=V]`: HOM setups (reference delay from setup.yaml) print the corrected value from the forward model with its
  corrections in order and the artifact IDs. CHSH setups refuse `--calibrations` until M3.
- `dadl quickstart hom-snspd` calibrates simulated dark and pair runs and prints the corrected line from
  the forward model, saying why the SNSPD's dead time is not corrected.
- `python -m dadl.analysis.propagation hom`: the HOM regimes of the `propagation` gate
  (hom-snspd, bright-background SPADs), with the calibration-uncertainty-off control reported. A
  new CI job runs them on the fast tier (200 replicates).

### Phase 3, M3: corrected CHSH and the `propagation` gate

- `chsh(..., calibrations=, channels=)` adds a corrected S from the forward model
  (`dadl.forward.corrected_chsh`) as a labelled secondary value; the headline stays raw. Each
  setting pair's four outcome combinations get their own accidentals (off-peak bins), live
  fraction and grid capture. The live fraction's conditional term models dead-time selection:
  at 3e6 pairs/s and 50 ns, raw S read 1.6 sigma high on average over 12 runs and the corrected
  S 0.26 sigma low (+-0.29). The corrected point uses each dark artifact's reported dead time, the
  shortest gap: an upper edge, high by about one mean gap divided by the number of gaps. The
  interval comes from the posterior draws, so the corrected point sits low in it (about 0.4 sd for
  S in chsh-dead-time). Efficiencies are not corrected (fair sampling stays assumed).
- The default CHSH rule allows no correction: each of `dead_time`, `accidental_subtraction` and
  `window_capture` needs `Policy.allow("chsh", ...)`. Delay-jitter artifacts must cover all four
  channel pairs, and dark artifacts every channel, or the corrected value is refused rather than
  half-modelled.
- `dadl analyze` takes `--calibrations` for CHSH setups. `dadl quickstart` prints a corrected S
  from simulated pair runs under a separately reported, loosened policy; the raw headline keeps
  the default policy. Quickstart artifacts have a fixed creation time, so their IDs repeat.
- `python -m dadl.analysis.propagation {hom,chsh,all} --jobs N`: regimes hom-snspd, hom-bright,
  chsh-spad and chsh-dead-time, the last also the control regime (calibration uncertainty off
  must undercover). Misfit artifacts are accepted, labelled and counted (at most 2% of
  replicates).
  `--shard k/n --out FILE` runs one slice of the replicates and `merge FILES` gates the sum;
  CI runs chsh-dead-time in eight shards. Progress is printed per finished chunk.
- The forward model's predictive check (`dadl.forward_check`): wherever a delay-jitter artifact
  gives the peak shape (window capture), each analysed coincidence histogram (HOM, fringe, all 16
  of CHSH) is refitted with that shape (a Gaussian, or the mixture at the declared offsets) on a
  stepped background over |b - a| < W, and a Bonferroni-adjusted grouped Pearson p below 0.001
  raises the new `Misspecified` error unless the policy allows `misfit_forward_model` (then
  labelled). The p-value is in `counts["forward_check_p"]`. A reference run declared without its
  path delay is caught the same way.
- Misspecified regime `hom-crosstalk`: hom-bright with 0.1% crosstalk each way at 400 ps. It
  passes when the check flags at least 90% of replicates; clean regimes may be flagged in at
  most 1%.

### Phase 3, M4: series dead-time model

- Simulator: `DetectorParams.dead_time_p_ps` (setup.yaml `dead_time_p_ps`, default 0), a
  paralyzable stage before the non-paralyzable `dead_time_ps`. Every avalanche candidate
  (photon, dark count, afterpulse, crosstalk) restarts it, in the filter path and both event
  loops. Two `series` cases join the shared Rust/Python test vectors.
- `dadl.calibration.deadtime.calibrate_deadtime_series` (kind `deadtime-series`, decision D4;
  method `series-gap-mle` 1): tau_p and tau_np from two or more illuminated acquisitions at
  different rates. The gap law of the series chain is exact (a delay equation solved as a finite
  series, then its dominant mode); afterpulsing is held at the dark artifact's estimates with its
  uncertainty propagated by a linear sensitivity. tau_p's interval is the profile likelihood's;
  when the hard-edge model is not rejected, tau_p is an upper limit (unidentified, stand-in
  draws). The likelihood ratio against the hard-edge model is recorded.
- CLI: `dadl simulate SETUP --deadtime-scan [--rates R1,R2,...]` (one rate per setting) and
  `dadl calibrate deadtime SETUP --tags scan.npz [--dark DIR] --out DIR`.
- `coverage` regimes `series-spad` and `passive-quench` (1,000 replicates each in CI; locally,
  before the review rounds, tau_p covered 95.5% and 96.9%, tau_np 99.2% and 97.3%). tau_np is gated at least 95%: the
  edge density is about 1 gap per tick (series-spad 1.5, passive-quench 0.7), so the shortest gap
  equals tau_np in about 79% and 45% of replicates; the simulator's integer dead time makes the
  tick-quantile interval over-cover (design 97.7% and 95.8%, observed 99.2% and 97.3%), so it is
  gated at least 95%; this gate has little power against edge-density errors (M10 item).
- `python -m dadl.calibration.coverage` gains `--shard k/n --out FILE` and `merge FILES`, as the
  propagation gate has; CI runs passive-quench in twelve shards.
- Stays out: soft recovery of SNSPDs and passive-quench efficiency recovery. The analyses'
  forward model uses the dark calibration's non-paralyzable dead time only. With paralysis the
  live fraction is overestimated: for tau_p close to tau_np by about 0.5% at a 10% dead fraction
  and 3-4% at the 20% limit, per channel (about twice that for a pair); ratios such as S and V
  cancel much of it.
- Review council round 1 (method `series-gap-mle` 2):
  - the interruption check tests each scan run on its own bounds, Bonferroni over runs (runs of
    unequal length, or taken one after another on a running clock, were flagged misfit);
  - the survival function continues past 8 tau_p with both real modes (the tail was 20-60% low
    near n tau_p = 1);
  - the afterpulse kernel is averaged over each grid cell, so the profile of tau_p is smooth;
  - `afterpulse_pending_max` is recorded, with a note above 0.5% (pending afterpulses bias tau_p
    low and are not modelled);
  - tau_p's point, upper end and draws are capped at tau_np's;
  - the afterpulse uncertainty is added to tau_p's upper limit too;
  - default scan rates are photon rates scaled by the smallest efficiency, so the avalanche rate
    x dead time spans 0.05-0.5 over four settings (0.1 and 0.4 for two);
  - `coverage.merge` refuses shards that disagree on the gate or its rule.
  - coverage re-run locally at 1,000 replicates: tau_p covered 95.5% (series-spad) and 93.8%
    (passive-quench), tau_np 99.2% and 97.3%, misfits 0.3% and 0%.

- Review council round 2 (corrected values; no method version changes):
  - a dead time known only as an upper limit gives the corrected value and its interval at the
    bound (the two-point value was at zero, so it could fall outside its own interval); CHSH
    records `corrected_at_zero_dead_time` and `corrected_at_bound_dead_time` in its counts as HOM
    and fringe do. "Threshold crossed only after correction" also counts the range's ends, and a
    range wider than the corrected value's error is labelled;
  - a stand-in jitter (a coarse tagger's upper limit, status ok) is used at its bound and labelled
    when the capture depends on it (every corrected analysis raised KeyError);
  - raw CHSH fits each outcome pair's peak over the four settings and records the offsets and S
    reweighted by the window captures; delays that miss the peaks are labelled "unequal window
    capture" (A+ and B+ 400 ps late, declared 0: S 2.20 for 1.98, unlabelled before);
  - corrected CHSH refuses (`Misspecified`, override `misfit_forward_model`) calibrated delays that
    no longer match the analysed peaks: a 200 ps drift since calibration read corrected S 2.03-2.05
    for 1.98 and passed the predictive check;
  - the predictive check adds a start with the background from the outer bins (the floor start
    stalled on a clean histogram whose peak fills the window: p about 1e-206);
  - `hom_from_tags(calibrations=...)` requires `reference_delay_ps` (without it the reference was
    treated as unsplit and V read 10 sigma low, unlabelled);
  - fringe results carry a two-point label, record each channel's singles ratio between the max and
    min acquisitions, and are labelled when it differs beyond 4 sigma (+3% pump drift on the max run
    read 0.707 for 0.699);
  - a forward-model corrected visibility above 1 is labelled, not clipped;
  - corrected values cite only the artifacts the model read; `deadtime-series` artifacts are
    skipped (and labelled) for HOM, fringe and CHSH; the set-level "parameter known only as an upper
    limit" label is dropped (the forward model says which parameter); afterpulsing reported by a
    dark artifact is labelled as not modelled; a set loaded without a `DatasetContext` is labelled;
    the CHSH accidentals label says efficiencies are not corrected; the HOM two-point label names a
    reference delay inside the coherence time, and the CLI warns when setup.yaml's
    `coherence_time_ps` makes it so; the CLI's "V reads low" warning says the corrected value models
    the capture when it does;
  - an artifact's `valid_for` starts at `acquired` when recorded, not at `created`.

- Review council round 3 (method `series-gap-mle` 3):
  - CHSH capture balance: each outcome pair's pooled peak is refitted on a region 5 sigma wider than
    the window with bins of sigma/3 (narrow peaks sat in one or two bins; wide ones traded their
    offset against the background), with standard errors; S reweighted by the captures records its
    standard deviation (`capture_sd_data`). "Unequal window capture" needs the shift to exceed half
    a standard error and 3.29 of those deviations, or a low capture with the peak off the declared
    delay; a low capture with the peak on it is labelled "window narrow against the pair jitter"
    (matched delays at 495 ps jitter and a 0.8 ns window: 9/20 labelled before, 0/20 now);
  - the delay-drift refusal adds a noise term (the artifacts' capture draws and the peak fits) and
    records `delay_drift_diff` and `delay_drift_z`: false refusals on current artifacts were 2-65%,
    now 0/264 in the council's runs; 200-400 ps drifts are still refused;
  - the predictive check polishes a histogram's fit when its p is below 0.01, with a second start
    at the highest bin (L-BFGS-B stalled with the start delay tens of ps off: p 2e-6 on clean data);
  - the corrected CHSH label says what was modelled and what was not (it said dead time and window
    capture were modelled when neither was); a corrected S without dark artifacts on a raw value
    with dead-time selection is labelled so;
  - an artifact recorded for another detector serial than the setup's is outside its validity
    (`CalibrationSet.load(..., channels=)`, passed by the CLI); a dark dead time undercut by more than
    max(10, 0.1%) of the analysed gaps is refused (`forward.check_dead_time_gaps`), and the shortest
    gaps are recorded (`min_gap_ps`); both are overridable with `calibration_outside_validity`;
  - a recorded temperature or bracket the dataset does not give is a mismatch, as a bias is; an
    artifact with no validity range is labelled; corrected results record the dataset context they
    were checked against (`counts["dataset"]`); two artifacts for one channel pair in either order are
    ambiguous;
  - the CHSH CLI prints the dead-time range when a dead time is a bound;
  - docs: validity, tagger settings and bias are enforced at use, not only recorded;
  - series dead time (method `series-gap-mle` 3): each run's rate is fitted on both branches of
    n exp(-n tau_p) and the better kept (a 6e7 cps run past the count-rate maximum was fitted at
    7e6 and biased tau_p to tau_np, 4/4 intervals excluding the truth); the profile of tau_p is
    searched for every minimum and the interval spans those inside it. The branches differ only
    through the afterpulses, and below the count-rate maximum the high branch at a lower tau_p can
    fit better by up to about 4 NLL (a scan to 1e7 cps at tau_p 30 ns read 25.9 ns, interval
    excluding 30), so the interval is joined with that of the fit held to the low branch, with a
    note when this widens it. Without afterpulsing a run past the maximum notes both rates;
    `survival` and `gap_cdf` accept an integer tau_p.

- Review council round 4 (method `series-gap-mle` 4):
  - deadtime method version 4: the round-3 joined interval is removed; intervals on
    below-maximum scans are narrower (pq p90 5.1 ns from 12.9; s30 2.4 from 15.2), tau_p and the
    fitted rates no longer jump to the high branch (about 5% of scans in round 3). Each run is
    taken below the count-rate maximum unless the fit with runs past it is better by more than
    `BRANCH_MARGIN` = 3 in NLL; within the margin the interval spans both fits and a note names
    the runs (`branch_nll_gain` in the conditions). Round 3's reason for the join was wrong:
    below the maximum the high branch fits worse in expectation (Asimov: by 2.5-5.6 NLL), not
    better by a few. The afterpulse-sensitivity refits hold each run's branch (they jumped to the
    other minimum: a shift sd up to 1.5e4 ps). The calibration refuses duplicate timestamps, the
    same run twice, and an afterpulse artifact of another kind or channel, and notes a misfit one;
  - the dead-time gap check drops its 0.9 factor (a dead time up to about 16% too long passed)
    and adds an upward test: a shortest analysed gap more than 9 mean edge spacings above the
    dead time's upper end is refused (22 ns calibrated on a 50 ns detector passed), for upper
    limits too. The floor per channel is in the counts (`gap_check_floor_ps`);
  - the detector-serial check runs on the API paths (`chsh`, `hom_from_tags`, `fringe_from_tags`
    with calibrations), not only through `CalibrationSet.load(channels=)`;
  - corrected CHSH runs the predictive check only with delay-jitter artifacts, as HOM and fringe
    do (an accidentals-only corrected S was refused as a misfit);
  - a dark artifact for only some channels gives the named error in HOM and fringe (a bare
    `KeyError`, a CLI traceback, before);
  - `calibrate pairs` prints the window capture from the forward model's closed form (its own
    integral read about 0 for a narrow peak in a wide window); HOM's forward path records
    `subtracted` 1.0; `valid_for` without `from` starts at acquisition or creation, and without
    `to` has no end (a `KeyError` before); the schema says validity is enforced; the CLI help names
    CHSH; the README lists the counts keys;
  - CI: a `gates` job fails unless every coverage and propagation job passed (a failed shard
    left its merge job skipped, which counts as passing); chsh-spad propagation runs in four
    shards (it timed out at 45 minutes); propagation chunks keep every worker busy; the coverage
    gate prints a line per finished chunk; a newer push to a pull request cancels its running
    checks.

- Review council round 5 (method `series-gap-mle` 5):
  - an afterpulse still pending when a gap starts is modelled to first order (one pending
    candidate, Exp(tau_a) after the record by memorylessness; `afterpulse_pending_max` is now its
    fraction pi). Checked against the simulator at 6e6 and 1e7 cps; it removes a tau_p bias of
    about -170 ps that did not shrink with data (0.5 s runs: coverage 88.1% to 92.5%, bias
    +18 +- 19 ps). The remainder, about pi (p + pi) of gaps, is recorded as
    `afterpulse_unmodelled`, with a note when it times sqrt(gaps) exceeds 2 (the 0.5% note is gone);
  - `BRANCH_MARGIN` 5 (was 3): about 0.5% of below-maximum scans switched to tau_p near 30 ns,
    every one a miss; the false-switch rate is now at most about 8e-4;
  - `deadtime-series` artifacts cite the dark artifact's inputs, and its ID (role
    `previous-artifact`) from the CLI; the held sds of p and tau_a are in the conditions;
  - corrected HOM, fringe and CHSH records carry `counts["monte_carlo"]` (seed and draws);
  - `valid_for.conditions` is enforced at use, as the top-level conditions are;
  - coverage re-run locally at 1,000 replicates (CI seed): tau_p covered 95.5% (series-spad) and
    94.6% (passive-quench), misses below/above 19/26 and 28/26; no branch switches, 0 and 54
    unions; tau_np 99.2% and 97.3%; misfits 0.3% and 0.0%. These replace the earlier figures
    (round 4, not recorded then: 93.2% passive-quench, misses mostly low);
  - docs: dev-branch validity rule and version, schema wording on `deadtime-series`, gap-check
    floor at low dark counts, unreleased features marked, counts keys, CLI help and refusals.

- Review council round 6 (method `series-gap-mle` 6; `dark` unchanged at 7):
  - **Behaviour change (owner to accept in the release notes):** an identified dark dead time
    whose dead fraction at the analysed rate exceeds 20% is refused, as a bound was (it was only
    labelled; there is no override). A 22 ns SPAD above about 9 Mcps loses its corrected value
    until M5's series route. To correct the rest, move the dark artifacts out of the calibrations
    folder;
  - an SNSPD's 30 ns read as an identified 20 us in 7 of 200 default dark runs, none covering,
    and gave V = 0.55 +- 6.2 with exit 0; the dead-fraction refusal above refuses all 7 at the
    referee's 85 kcps. When three widths reach below zero the gap check tests the 95% lower end for
    a label only (`DEAD_TIME_BELOW_DARK`, in `corrected_labels`; `gap_check_floor_kind` "3 widths"
    (refusing) or "lower end (label only)"): the lower end misses in 5% of channels by design, and
    refusing it (as first implemented) refused 42 of 200 chsh-dead-time replicates (council round 6
    addendum). The 20 covering loose chsh-spad identifications are neither refused nor labelled.
    `calibrate dark` warns when the lower end is more than 3x the declared dead time;
  - results keep the corrected value's labels apart (`corrected_labels`, a new record column;
    `all_labels` is both): threshold labels are evaluated per value, and calibration labels come
    only from the artifacts the value used. The CLI prints them under the corrected line;
  - forward-model HOM and fringe values carry "secondary; the headline is raw V" and the CLI
    marks the HOM line "(secondary; not an indistinguishability claim)";
  - corrected records carry `calibration_excludes` (the used artifacts' `uncertainty.not_included`),
    printed on one line; stat and calibration parts print with 2 significant figures;
  - refusal text names the analysis (HOM: "corrected values are secondary; the headline V uses raw
    counts") and the CLI remedy; HOM analyze prints the declaration caveat; quickstart and
    `calibrate pairs` say that CHSH window capture needs four heralding runs;
  - series tau_np: the interval allows one tick either side and the draws add the rounding of two
    floored times (an off-grid truth was missed every time: 0 of 400 at 25 ps; now 98.75%);
    `tag_tick_ps` recorded. On the simulator's grid the widening is conservative by
    construction, so the coverage gate's tau_np figure says little; the off-grid test checks the
    tick term;
  - series: a dark artifact recorded under other conditions (bias, temperature, tagger settings)
    is refused unless `allow_condition_mismatch` (CLI `--allow-condition-mismatch REASON`) gives
    a reason, which is noted; `dark_bias_V` and `dark_temperature_K` are recorded;
  - series record: the tau_p budget (`dead_time_p_stat_low_ps`, `dead_time_p_stat_high_ps`,
    `dead_time_p_u_afterpulse_ps`, `dtau_p_dp_ps`, `dtau_p_dtau_a`, `afterpulse_p_tau_corr_used`),
    the scan range (`scan_rate_low_cps`, `scan_rate_high_cps`, `scan_duration_low_s`,
    `scan_duration_high_s`), `uncertainty.kind` per kind; `lr_hard_edge` renamed
    `two_ln_lr_hard_edge` and `afterpulse_tau_sd_held` `afterpulse_tau_ps_sd_held`; notes say
    "2 ln(likelihood ratio)" and "avalanche-candidate rate ... s^-1"; tau_np's definition
    separates the quantity from its estimator;
  - union intervals are drawn as an equal-weight mixture of the two branch fits' split normals
    (one split normal put 18% of draws in the gap between disjoint branches; now 2.8%, the
    components' tails);
  - printed intervals state their level ("95% interval", "95% upper limit") and round to two
    significant figures of the smaller half-width (display only);
  - docs: README methods section for corrected values and the limits found (identified dark dead
    time, model-conditional calibration part, tick, four heralding runs, flat priors of stored
    draws); phase 3 plan: M5 entry items (what NUTS samples, starts, boundaries, priors, shared
    row index, stored draws gated, scan range at use) and M10 items; benchmark record: steady-
    state speed;
  - coverage re-run locally at 1,000 replicates (CI seed): tau_p covered 95.5% (series-spad) and 94.6%
    (passive-quench), as in round 5; tau_np 99.9% and 99.0% (was 99.2% and 97.3%; conservative by
    construction on the grid, gate "at least 95%" unchanged); misfits 0.3% and 0.0%;
  - propagation gate at 200 replicates: hom-snspd, hom-bright, hom-crosstalk and chsh-spad
    unchanged (94.5/94.5, 95.0/95.5, 30.5/32.5 with 100% flagged, 93.5/93.5); chsh-dead-time at 400
    replicates: corrected 95.5%, control 87.5% (band 91.0-98.2%), misfit artifacts 0.8%, no flags,
    no refusals (the 1,000-replicate run is CI's).

- Review council round 7 (method `series-gap-mle` 7; `dark`, `delay`, `accidentals` unchanged):
  - T3 tick: a recorded T3 gap is within one dtime resolution of the true gap whatever the sync
    period, but a measured (non-whole) period left the tags a 1 ps tick, so tau_np widened by
    1 ps on 25 ps data (10 of 20 covered). `PtuFile.tick_ps` and `tag_tick_ps(folder)` give the
    dtime resolution in T3 (plus 1 ps when tags are rounded); `calibrate_deadtime_series(...,
    tick_ps=)` takes the larger of it and the tags' own, recorded as `tag_tick_ps` with
    `tag_tick_declared`; `calibrate deadtime` and `calibrate pairs` pass it. tau_np then covers 20 of 20;
  - the T3 sync period (`tagger_MeasDesc_GlobalResolution`) is compared to 1e-4, within a folder
    and at use (two files 10 Hz apart at 80 MHz were refused; 80 against 76 MHz still is);
    mismatches print 10 significant digits and both values;
  - T3 input from a sync source: analyses label it ("T3 input: DADL models CW light") and refuse a
    corrected value unless the source is declared CW (`source.type: cw`, `source_type="cw"`);
    `read_tags` carries the sync period in each acquisition's metadata (`RunTags`). A pulsed
    76 MHz HOM read V = 0.962 +- 0.001 for a true 1, unlabelled. Behaviour change: T3 users lose
    corrected values until they declare the source;
  - T3 files whose dtime range is shorter than the sync period get a `t3-range` event and note,
    and `read_tags` refuses them; T3 with `HWSync_Divider` > 1 is refused by `calibrate dark` and
    `calibrate deadtime` and labelled by analyses (convention unverified; M9 entry item);
  - delay-jitter and deadtime-series artifacts from another tagger (model or serial in
    `instrument`) are outside their validity; other kinds are labelled;
  - sharded gate merges recompute the gates from the code, refuse an empty merge, and refuse a
    result set other than the one the recorded family (or analysis) and regimes produce;
  - a posterior `.npy` header is checked against the manifest before the table is read (a
    192-byte file declaring 2^52 rows raised `MemoryError`);
  - docs and hygiene: T3 memory figure (about 170 bytes per record), the T3 limit in the READMEs,
    `persist-credentials: false` in `ppl-bench.yml`, the CmdStan binaries removed from the tree;
  - tests: gate verdicts on synthetic results, coverage merge tiling, channel order and declared
    delays on HOM, fringe and CHSH, wired labels observed through the analyses, behavioural oracles
    for the dead-time intervals, the predictive check's statistics; two seed-fragile tests reworked
    (CHSH selection at 2 s: 20 of 20 seeds; HOM with the capture modelled: 20 of 20 within 3 se);
    slow tests run in one Linux CI job;
  - known: the predictive check flags clean T3 tags with a measured sync period (p about 1e-100 at
    25 ps, 12500.3 ps); reported to the council, not loosened (xfail test);
  - coverage re-run locally at 1,000 replicates (CI seed): tau_p 95.5% (series-spad) and 94.6%
    (passive-quench), tau_np 99.9% and 99.0%, misfits 0.3% and 0.0%: identical to round 6 (the
    simulator's tags are on the 1 ps grid and no tick is declared);
  - propagation gate at 200 replicates, unchanged: hom-snspd 94.5/94.5, hom-bright 95.0/95.5,
    hom-crosstalk 30.5/32.5 with 100% flagged, chsh-spad 93.5/93.5, chsh-dead-time 95.5/85.0;
    both sharded merges run end to end; the 44 mutants the ruling names are all killed.
- Review council round 8 (method `xcorr-gaussian-mle` 9; `dark`, `deadtime`, `accidentals`,
  `dcr`, `heralding` unchanged):
  - T3 lattice: `PtuFile.lattice_ps` and `tag_lattice_ps(folder)` give the lattice pair delays sit
    on (the whole-ps dtime resolution in T3, the tick in T2); `tick_ps` keeps its round-7 meaning
    and serves dead-time intervals only. `calibrate pairs` passed the 26 ps gap bound as the tick
    on 25 ps data at a measured period: delay -12 ps, jitter covered 0 of 10, every artifact a
    misfit. It now passes the lattice and the sync period; the delay fit's off-grid note no longer
    fires on T3 cross-sync residues. `read_tags` carries `lattice_ps`; the predictive check and the
    CHSH peak fits bin on it, with edges half a lattice step off the delays' modal residue. T3
    with a resolution that is not whole ps skips the check (`forward_check_note`). The round-7
    xfail passes;
  - delay fit: the predictive check uses (groups - 1) degrees of freedom, as the forward check
    (groups - 5 gave p < 0.05 in 9.5% of clean SNSPD pair runs);
  - one refusal budget: on clean data an analysis is refused at analysis time with chance at most
    1e-3, split as the forward check 4e-4 (`MISFIT_ALPHA`), the CHSH drift 2e-4 (`DRIFT_Z` 3.72;
    the capture label keeps 3.29 as `CAPTURE_Z`), the off-peak check 2e-4 and the gap tests 2e-4
    (README, "Refusals and their rates"). Behaviour change;
  - dead-time gap tests: the downward floor is point - ln(K / 1e-4) / (N f) from the dark run's
    gap count and dark rate (`gap_check_floor_kind` "exponential law"; the 3-width floor refused
    clean channels at 2.5e-4 to 3.7e-4); the upward threshold is 100((K / 1e-4)^(1/100) - 1)
    edge spacings (11.2 at four channels). New label `DEAD_TIME_ABOVE_GAPS` when the dark dead
    time is longer than many analysed gaps; the correction keeps the dark value (limit stated;
    the fix is M5a, with the second-order pair term and a new control regime);
  - pulsed light on T2 or plain tags read V 10-40 sd off, unlabelled. Analyses that subtract
    accidentals now check each acquisition's off-peak windows for flatness (exact dispersion at
    20,000 simulations, and trend; `offpeak_p`) and refuse a value whose windows are not flat
    (`Misspecified`, override `misfit_forward_model`, then `OFFPEAK_ACCEPTED`). Behaviour change;
  - the CHSH drift override has its own label, `DELAY_DRIFT_ACCEPTED`. The propagation gate counts
    the predictive check, the drift and the off-peak check under their own keys (the 1% cap
    applies to their sum) and records refused replicates by type ("refused: <type>", not covered)
    instead of aborting; any refusal still fails a regime (tolerance: owner item O12);
  - wording: the 20% dead-fraction refusal gives its real reasons; user-facing text no longer
    names milestone M5;
  - plan: M5 split into M5a (series dead time at use) and M5b (NUTS); M2b, M3b and M4b take most
    of M10's smaller items; an owner-items section (O1-O14); D3 as-built note; dependency fixes;
  - CI: the slow-test job is in the `gates` job's needs;
  - delay coverage at 1,000 replicates (CI seed) identical in all six pairs regimes; delay misfits
    fell in sparse-snspd (0.3% to 0.0%) and broad-irf (0.1% to 0.0%). Deadtime coverage identical
    to round 7 (95.5/99.9 series-spad, 94.6/99.0 passive-quench);
  - propagation gate at 200 replicates (CI seed), no refusals and no flags in any clean regime:
    hom-snspd 94.5/94.5, hom-bright 95.0/95.5, chsh-spad 93.5/93.5, chsh-dead-time 95.5/85.0,
    hom-crosstalk 30.5/32.5 with 100% flagged (misfit artifacts 1.0% to 0.5%). Seed 4242:
    chsh-spad 93.0/93.0, chsh-dead-time 97.0/85.0, no refusals.

### Published results as executable tests

- `python/tests/test_published_results.py` rebuilds three published results
  from simulated time tags (design doc, the open high-risk review item):
  - accidental subtraction (Thompson, quant-ph/9903066, on Tittel et al.
    1997): a noisy-link source with accidentals about 30% of the fringe
    maximum reads raw V = 0.44 (exact interval below 1/sqrt 2); DADL's own
    off-peak estimate would lift it to about 0.8, and the default policy
    refuses that; a loosened policy labels it;
  - the coincidence-time loophole (Christensen et al., arXiv:1503.07573;
    Larsson et al., arXiv:1309.0712): a local model with setting-dependent
    detection times reads S = 2.77 (tens of Gaussian sigma) with a window
    centred on Alice's detections, 1.59 on DADL's planned 10 ns grid (both
    match the model's closed forms), and a data-chosen window is refused. A
    grid only a little wider than the attacker's stagger still reads S > 2:
    the fixed grid does not close the fair-coincidence assumption for CHSH,
    which the result's label already says;
  - the NIST calibration comparison (Polyakov and Migdall, SPIE 6583,
    65830A): the heralding calibration at the 0.18% (k = 1) precision regime
    recovers the injected efficiency, and agrees with a substitution
    calibration against a reference detector within combined uncertainty.
- Fact-check of the figures these tests use: 45% to 82% (with about 30%
  accidentals), the 71% limit, "over 2700-sigma" and the 0.14(14)%, 0.18%,
  0.17% and 0.5% (k = 1) NIST figures match their sources. Two wordings were
  corrected in the design doc and the site: the NIST methods *differ by*
  0.14(14)%, within the comparison's uncertainty (not "agree to within"), and
  the 2700-sigma violation was of the CH inequality with event-centred
  windows, with no significant violation (not "none") with a predefined one.
- No package code changed; no bugs found.

### PicoQuant PTU T3 mode

- `dadl.adapters.ptu` reads T3 files from the same devices as T2: PicoHarp
  300 (0x00010303), HydraHarp v1 (0x00010304), HydraHarp v2 (0x01010304),
  TimeHarp 260 N/P (0x00010305, 0x00010306) and generic T3 (MultiHarp,
  PicoHarp 330: 0x00010307). Record types, bit layouts and overflow rules
  follow PicoQuant's demo code (PicoQuant-Time-Tagged-File-Format-Demos,
  PTU/C/ptudemo.cc and PTU/Python/Read_PTU.py). Tags are nsync x sync period
  + dtime x resolution in int64 ps, exact when both are whole ps, rounded
  with a note otherwise (the period's fraction is kept exactly, so late
  syncs do not drift); metadata records the time model and its assumption
  of a constant sync period. Input numbers match T2 (HydraHarp family:
  record channel k is input k + 1); input 0 holds the syncs that started a
  recorded measurement, not every sync, and a `t3-sync` event says so.
  Overflows, markers (at their sync) and unknown special records are
  accounted for as in T2; PicoHarp T3 channels outside 1-4 and tags out of
  time order (a dtime beyond the header's sync period) are refused. Before,
  T3 files were refused.
- PTU folders (`--tags DIR`) accept T3 files; input 0 is read only when
  mapped. `tag_tick_ps` is the grid of both period and resolution, and the
  T3 dtime resolution is recorded with the tagger settings.
- `dadl.adapters.ptu_write.write_ptu(..., sync_period_ps=...)` writes
  simulated T3 files (`encode_t3_records`), byte-identical to the test
  suite's record-by-record writer.

### Crosstalk in the simulator (phase 1)

- `dadl.sim.Crosstalk(source, target, probability, delay_ps, jitter_sigma_ps)`:
  opt-in crosstalk links for `detect_channels`, `simulate_pairs`,
  `simulate_hom`, `simulate_chsh_setting` and `simulate_chsh`. Every
  recorded avalanche on the source can trigger one on the target at a
  jittered offset (a Si SPAD breakdown flash after 2 n_g L / c, or zero for
  electrical crosstalk); crosstalk avalanches start dead time, afterpulse
  and trigger crosstalk in turn, all in one time-domain event loop. Off by
  default: without links every channel is simulated as before and results
  for a given seed are unchanged.
- Closed-form tests: a thinned shifted copy on a dark channel, a
  cross-correlation excess of p x singles at the offset (none at the
  mirror offset), geometric cascades both ways, crosstalk lost inside the
  target's dead time, recorded crosstalk starting a dead time, crosstalk
  avalanches afterpulsing.
- Not yet: crosstalk in `setup.yaml`, a crosstalk calibration, or a label in
  the HOM analysis (in the two-point method an in-window crosstalk excess
  looks like coincidences). On `hom-snspd`, 1% zero-offset crosstalk A to B
  lowers V from 0.919 to 0.904 (README, Known limits).

## 0.1.1 (2026-11-14)

Release fix, no code changes to the packages. The 0.1.0 upload of
`dadl-quantum-core`, a new PyPI project, was refused: its pending trusted
publisher had the same repository, workflow and environment as
`dadl-quantum`'s existing publisher, which PyPI matched instead. The core now
publishes from its own environment, `pypi-core`. 0.1.0 reached crates.io only;
0.1.1 is the first release of `dadl-quantum-core` and the first phase-2
release of `dadl-quantum` on PyPI.

## 0.1.0 (2026-11-13)

Phase 2: calibrations with simulation-checked 95% intervals, content-addressed
artifacts, PicoQuant PTU input, a hardware twin and the optional compiled core,
after an outside-expert review of the assumptions and ten further review
cycles. Pre-alpha: read the README's "Known limits".

### Phase 2, milestone 1: data model and the dark calibration

- `dadl calibrate dark` and `dadl.calibration.calibrate_dark`: dark rate,
  afterpulse probability and time constant, and dead time, fitted jointly from
  one dark acquisition's inter-arrival gaps (maximum likelihood, Laplace
  posterior; an upper limit on afterpulsing when none is detected), with a
  posterior predictive check.
- `dadl.artifact`: content-addressed calibration artifacts (RFC 8785 canonical
  manifest, SHA-256 ID, `.npy` posterior table), verified on read. JCS test
  vectors in `tests/vectors/jcs.json`.
- `dadl.data` (`ChannelId`, `TagSet`, `Dataset`), `dadl.models`
  (`DetectorModel`), `dadl.stats` (Wilson, Clopper–Pearson, Beta and Gamma
  posteriors, Laplace).
- `schema/setup-1.schema.json` for `setup.yaml`. `calibration/1` (still a
  draft) gains `summary`, `notes` and `unidentified_params`, and the `npy`
  posterior format.
- `dadl simulate --dark` writes dark acquisitions for calibration.
- Simulator: afterpulsing is self-exciting (afterpulses cascade and every
  avalanche starts dead time), matching the closed form
  `closed_form.dark_afterpulse_rate`. Simulated tags for a given seed differ
  from 0.0.3; the quickstart prints S = 2.6941 instead of 2.6890.
- CI: a coverage job runs the 1,000-replicate gate on every pull request.
- New dependency: SciPy (>= 1.13).

### Phase 2, milestone 2: delay, jitter, accidentals, heralding, DCR(T)

- `dadl.histogram.cross_correlation` and `dadl::histogram::cross_correlation`
  (Rust): histogram of delays b - a within a range, held to new shared vectors.
- `calibrate_delay` (kind `delay-jitter`): relative delay and combined jitter
  from a Gaussian peak plus flat background fitted to the cross-correlation
  (Poisson maximum likelihood, Laplace posterior, grouped chi-square check).
- `calibrate_accidentals` (kind `accidentals`, new in `calibration/1`):
  expected accidentals in the declared fixed-grid window from 200 shifted
  windows beyond a 200 ns guard (CW) or at the side peaks of a pulsed pump;
  Gamma posterior and an off-peak dispersion check. R1·R2·Δt is a
  cross-check only.
- `calibrate_heralding` (kind `heralding-efficiency`): (C - A) / S with the
  accidentals' Gamma-Poisson uncertainty and a Beta posterior; an exact,
  deterministic interval.
- `calibrate_dcr_temperature` (kind `dcr-temperature`, new in
  `calibration/1`): Arrhenius fit of dark calibrations taken at several
  temperatures (`conditions["temperature_K"]`).
- `dadl calibrate pairs SETUP --tags F --herald A+ --signal B+ --out DIR`
  writes the delay/jitter, accidentals and heralding artifacts, using the
  setup's declared window and delays, and warns when the fitted delay
  disagrees with the declared one (it never moves the window).
- `simulate_pairs` and closed forms `arrhenius_dark_rate`,
  `accidentals_expected` and `heralding_ratio_expected`.
- Coverage gate: `python -m dadl.calibration.coverage {dark,pairs,dcr,all}`;
  the CI job now runs all six calibrations. The chance that a correct
  calibration fails its gate is lowered from 1% to 0.1% (band 92.6–97.1% for
  one parameter, 92.3–97.3% for four): the seed is fixed, so a chance failure
  blocks CI permanently. The delay/jitter gate read 96.9% for jitter on seed
  20261111, above the old 96.8% limit; three other seeds gave 95.3%, 94.7%
  and 94.0%. The accidentals gate first read 91.7% because mirrored off-peak
  windows are correlated; the windows no longer mirror each other (94.0%).

### Phase 2, milestone 3: PTU adapter and the hom-snspd example

- `dadl.adapters`: the `Adapter` protocol (Python subset of the Rust trait),
  `DeviceInfo`, `Capabilities` and typed `Event`s; `adapters/TEMPLATE.md`
  for contributors.
- `dadl.adapters.ptu` (`PtuFile`, `read_ptu`): DADL's own reader for
  PicoQuant PTU files in T2 mode (PicoHarp 300, HydraHarp v1/v2, TimeHarp
  260 N/P, MultiHarp). Exact integer ps conversion when the tick is whole
  ps; rollovers handled (including multi-rollover records); markers,
  rollovers and unknown records reported as events; every record accounted
  for; truncated files, trailing data, out-of-order records and T3 files are
  refused. Tests write fixtures with `python/tests/ptu_writer.py` and, when
  `ptufile` is installed, check that both readers give the same times.
- `analysis.type: hom` in `setup.yaml` (channels A and B; acquisitions "dip"
  and "reference" at `analysis.reference_delay_ps`; `source.coherence_time_ps`),
  `simulate_hom`, `dadl.analysis.hom.hom_from_tags` and the `hom-snspd`
  example: `dadl quickstart hom-snspd`, `dadl simulate` and `dadl analyze`
  work for HOM setups. The MVP tagger is now PicoQuant PTU instead of Swabian.
- Tag files take named settings (`read_tags(..., settings=...)`): HOM uses
  0 = dip, 1 = reference.
- Fix: the delay/jitter fit (method `xcorr-gaussian-mle` v2) treats each
  histogram bin as the whole-ps delays it holds. Version 1 was biased by
  +0.5 ps, visible only when the jitter is a few bin widths.
- `AdapterError` for malformed or unsupported vendor data.

### Phase 2, milestone 4: Python bindings to the Rust core

- `bindings/python`: a PyO3 crate (its own Cargo workspace) built with
  maturin into the optional package `dadl-quantum-core` (module
  `dadl_core`, abi3 wheels for Python 3.10+). It exposes dead time,
  occupied slots, fixed-grid coincidences, inter-arrival and
  cross-correlation histograms over NumPy int64 arrays, reading contiguous
  arrays in place. The GIL is held while the core runs, so no other Python
  thread can write to an array the core is reading.
- `dadl.backend()` and `DADL_BACKEND` (`python` / `rust`): `dadl` uses the
  compiled core when it is installed, else the pure-Python reference.
  Argument types are checked in Python for both; the shared vectors and
  random-stream comparisons run on both backends.
- `bindings/python/bench.py`: on 2,000,000 tags per channel (Apple M-series)
  the core is 61x faster for the non-paralyzable filter, 47x for fixed-grid
  coincidences, 48x for inter-arrival histograms and 5x for
  cross-correlation (already vectorized in NumPy). Recorded, not gated.
- CI: a bindings job builds the core on Windows, Linux and macOS and runs
  every test with `DADL_BACKEND=rust`. Release: abi3 wheels (manylinux2014
  x86-64 and aarch64, macOS universal2, Windows x64) and an sdist for
  `dadl-quantum-core`, published with Trusted Publishing.

### Hardware twin

- `dadl simulate SETUP --out DIR --format ptu [--device hydraharp|multiharp|picoharp|timeharp260p]`
  writes each simulated acquisition as a T2 PTU file (`setting0.ptu`, ...),
  with the device's tick, timer rollovers and header. `--tags DIR` reads such
  a folder in `analyze` and `calibrate`. The CHSH quickstart gives the same S
  through PTU files as through a tag file.
- `dadl.adapters.ptu_write.write_ptu`: a vectorized PTU writer, byte-for-byte
  equal to the independent writer in the test suite. Files it writes carry
  `CreatorSW_Name = DADL` and `DADL_Simulated = True`; the reader marks them
  `simulated`, and calibrations from them are recorded as `backend: sim`.
- `dadl.adapters.SimulatedTagger`: a live tagger emulated from a setup,
  streaming batches and injecting buffer overflows, dropped batches and
  configuration changes as typed events; `check_accounting` verifies that
  tags delivered plus tags reported lost equal tags generated.
- `dadl.simsetup` (`simulate_setup`, `simulate_setup_dark`), moved out of the CLI.
- Fix found by running the existing calibrations on twin data: the delay fit
  (method `xcorr-gaussian-mle` v3) finds the tagger tick from the data, rounds
  the bin width up to whole ticks, and models each bin accordingly. Version 2
  assumed 1 ps ticks and read 5 ps MultiHarp delays 2 ps high.
- Coverage analysis ([review record](docs/reviews/2026-11-phase2-review.md#follow-up-test-coverage-and-untested-regimes-2026-11-12)):
  the dark calibration was probed at six detector regimes. Its dead-time
  interval (method `interarrival-mle` v2) now uses exact quantiles of the
  shortest gap in whole ticks; v1 covered 89% with fast afterpulsing, v2 96%.
  A sparse, sharp coincidence peak no longer gives a NaN predictive check
  (which blocked saving the artifact); HOM says when accidentals are not
  subtracted.
- Coverage gates run at more than one regime: the dark calibration at a Si
  SPAD and a CMOS SPAD (20% afterpulsing), the pair calibrations at a CW SPAD
  source, a sparse low-rate SNSPD source and a broad (about 1.8 ns) IRF:
  12 gates, run in parallel (`--jobs`) in the existing CI check. In the sparse
  regime about 0.2 off-peak counts are expected in all, so the accidentals
  gate there requires coverage of at least the band's lower edge, 92.6% (declared
  before the first run).
- Fix found by the broad-IRF regime: the delay fit crashed on broad peaks
  (overflow, or a singular Hessian) because its starting width came from a
  fixed +-50 bins. The start now widens until it holds +-4 sigma, the fit is
  bounded, and when the background fits to zero the delay and jitter are
  reported conditional on it, with a note, instead of failing.
- Calibration probe campaign (three rounds; [review record](docs/reviews/2026-11-phase2-review.md#calibration-probe-campaign-2026-11-12)):
  - dark calibration v4: the dead-time interval averages the shortest-gap law
    over the uncertainty of f(0) and runs up to m when afterpulsing is not
    detected (it covered as little as 0% with a few hundred gaps); afterpulse
    intervals come from the profile likelihood below 5,000 gaps; the minimum
    is 20 gaps instead of 200, so cold DCR(T) runs calibrate; at low dark
    rates the dead time is reported as the bound it is, not as 0 ps;
  - delay: narrow peaks are refitted with finer bins (92% to 94-95%);
  - accidentals: an exact dispersion check for sparse counts (false alarms
    from 4.5% to 0.1%).
  - the coverage gates split each regime's replicates across processes
    (results identical to a serial run) and the profile intervals compute to
    the precision they need: all 12 gates in 95 s on 4 local cores, after the
    low-count fixes had pushed the CI job to 36 minutes.

### Twenty-round review of phase 2 ([docs/reviews/2026-11-phase2-review.md](docs/reviews/2026-11-phase2-review.md))

- Bindings hold the GIL while reading NumPy buffers (releasing it allowed a
  data race with other Python threads).
- Heralding efficiency counts heralds per slot (method `herald-ratio` v2):
  two pairs in one slot no longer halve it.
- Artifacts are written atomically; a manifest without a valid posterior
  entry raises `ArtifactError`.
- `dadl-quantum` uses `dadl-quantum-core` only when the versions match.
- CI builds the release wheels of `dadl-quantum-core` on every pull request.
- Notes: the dark rate includes undetected afterpulses (bounded by the upper
  limit); PTU files are read whole.

- Release: the crate is now published to crates.io by the release workflow with
  Trusted Publishing, in a `crates` environment, instead of by hand with a token.

### Outside-expert review of the assumptions ([docs/reviews/2026-11-expert-assumptions.md](docs/reviews/2026-11-expert-assumptions.md))

Five reviewers (detector physicist, statistician, metrologist, time-tagging
engineer, quantum optics experimentalist) rated 13 assumptions wrong and 14
questionable. Changes:

- **Heralding efficiency** (method `herald-ratio` v3): per-slot estimator
  η = (q − p)/(1 − p) with the accidental probability p propagated exactly;
  unbiased up to 20% accidentals (the old ratio read η(1 − p)).
  `calibrate pairs` refuses CHSH and HOM setups unless `--no-analysers`
  confirms a heralding run, and records delays, duration, singles rates and
  analysers in `conditions`. `dadl simulate --heralding` writes such a run.
- **Accidentals** (v2): the off-peak guard is at least 10 µs and ten afterpulse
  time constants past the dead time (afterpulses after true pairs had biased
  it +9% to +165%); a trend test joins the dispersion check.
- **Misfit status**: a predictive check below p = 0.001 sets
  `checks.status: misfit`; `dadl calibrate` refuses to write it unless
  `--accept-misfit` is given.
- **PTU data loss**: non-complete `TTResult_StopReason` codes become
  `data-loss` events and `incomplete` metadata, refused by the CLI; the
  duration comes from the data; PicoHarp channel values above 4 (sync 0, router inputs 1-4) are refused; the
  sync divider is reported.
- **Channel identity**: artifacts record the tagger serial and input, with the
  role and detector serial alongside; hardware setups need `hardware.tagger`.
- `source.type` is validated (`spdc-cw` only).
- **Exact visibility intervals** (Clopper–Pearson) with `NEAR_THRESHOLD_EXACT`
  and `LOW_COUNTS` labels; the delta method covered 63–93% near V = 1.
- **Dark calibration** (v5): "effective channel dead time"; up to three
  short gaps far below all others are dropped as double triggers (a single
  spurious gap near the dead time cannot be told apart and is noted, not
  dropped); profile intervals keyed on the afterpulse count;
  posterior draws reweighted to agree with the intervals; a stationarity check
  (split-half rate, longest gap) flags drift and latches as misfit.
- **DCR(T)**: refuses runs at different bias; records the temperature range;
  notes the SPAD-only validity and few-gap runs.
- **Delay** (v4): the tick comes from the adapter (the tag GCD is a
  cross-check); a two-level background models the dead-time step (bright
  background coverage 52% → 98%).
- `check_accounting` compares tag times, not counts; twin presets enforce the
  real input counts; a MultiHarp 150 N (80 ps) preset.
- Bell results carry "fair sampling and fair coincidence assumed", static
  settings and loophole labels, and coincidences per singles; HOM results are
  labelled as a two-point visibility at the declared delay.
- Artifacts gain an `uncertainty` block (statistical, 95%, components not
  included).
- Coverage gates: new regimes with afterpulsing SPADs, a bright background and
  an ideal heralding check against the true η; accidentals are judged against
  each replicate's own expectation.
- Docs: README uncertainty wording and "Known limits"; design.md notes on the
  dead-time model and misfit refusal; "signed" becomes "content-addressed,
  signing planned".

### Ten review cycles ([docs/reviews/2026-11-review-cycles.md](docs/reviews/2026-11-review-cycles.md))

Ten further cycles of AI-agent role-play review (three reviewers each, 30 in
all) found 47 problems rated wrong, all fixed or stated as limits; the
register lists every finding and its status. Method versions in this release:
dark `interarrival-mle` 7, delay `xcorr-gaussian-mle` 8, heralding
`herald-ratio` 6, accidentals `offpeak-grid` 3, DCR(T) `arrhenius-wls` 2; a
test now fails if a method's code changes without a new version.

**Breaking changes**

- `Window(width_ps, source, ...)`: `source` is required (no default "plan");
  `delays_ps` is read-only, and delays for unknown channels are refused.
- `write_tags(..., simulated=...)` is required.
- `write_artifact` refuses misfit calibrations unless `accept_misfit=True`.
- CHSH and HOM refuse duplicate timestamps and need one channel at delay 0;
  `bell_visibility` and `hom_visibility` need whole counts.
- `hom_from_tags` compares rates across acquisitions of different length;
  `hom_visibility(..., exposure_ratio=)`.
- Tag files must be in time order per channel (a clock reset is refused, not
  sorted).
- Corrected visibility above 1 is left out and labelled instead of raising.
- Artifact manifests: `summary` entries may have `value: null` with `bound`,
  plus `definition`; `posterior.stand_in`; named `inputs`; `acquired` and
  `instrument`; software versions; UCUM units. Non-converged fits are misfits.
- `as_record()` writes JSON cells and adds `dadl_version` and `simulated`.
- `dadl simulate` refuses to overwrite without `--force`; `--device` needs
  `--format ptu`; `--duration-s` applies to `--dark` and `--heralding` only.
- SciPy 1.13 or later.

**New**

- Labels: interrupted (latch, gaps at the start, middle or end of a run),
  dead-time selection, output imbalance, simulated input, accidentals
  subtracted, correction not computed, window not given, single basis.
- `calibrate pairs`: per-artifact refusal, window capture, coincidences to
  accidentals and pair rate (secondary), peak-placement check, declared dead
  time and afterpulse time constant passed; `calibrate dark`: per-channel
  handling, refusal of source-on and unit-error data.
- Record options: `--condition`, `--acquired`, `--created`, `--valid-until`,
  `--supersedes`; `load_calibration`; `dadl-quantum[core]` extra.
- Delay fit: exact tick model, two-level background, multi-start and Newton
  polishing, jitter upper limit below the tick; dark fit: time constants up to
  1 ms (free-running InGaAs); heralding: closed form with the conditional live
  fraction.
- Coverage gates: 10 regimes and 42 checks plus a misfit-share limit, split
  across CI jobs; `--regimes`.
- Input safety: YAML alias expansion and npz decompression limits, PTU header
  checks, histogram size limit in both backends, GIL declared for
  free-threaded Python.

## 0.0.3 (2026-11-11)

Applies the third review ([docs/reviews/2026-11-council.md](docs/reviews/2026-11-council.md)).

### Breaking

- `two_photon_visibility(..., analysis=...)` is replaced by
  `bell_visibility()` and `hom_visibility()`; the analysis type can no longer
  be chosen by the caller. HOM visibility now uses the dip convention
  (C_max − C_min) / C_max.
- Simulator and closed-form parameters are unit-suffixed to match
  `setup.yaml`: `dark_rate_cps`, `jitter_sigma_ps`, `pair_rate_per_s`.
- `load_setup` returns a typed `Setup` and refuses unknown or duplicate keys,
  wrong types (including `1.0e5`, which YAML reads as a string) and
  non-default policies.
- Rust: every public function returns `Result<_, DadlError>` instead of
  panicking; parameters are named with their units.
- `closed_form.series_hybrid_measured_rate` is replaced by
  `series_measured_rate` (the exact series model); the old formula is kept as
  `patil_usman_hybrid_rate`, documented as empirical.

### Fixed

- Corrected visibility error now propagates the accidentals' variance and
  uses the correct derivatives (the old error was about 10–30% too small).
- Corrected visibility above 1, negative accidentals and non-finite counts
  are refused.
- Unsorted input to dead-time filters and histograms is refused in every
  build, in Rust and Python (release builds used to return wrong answers).
- Overflow near the int64 limits is an error instead of a saturated or
  wrapped value (Rust used to report a false coincidence).
- Float timestamps are refused with a hint to convert to int64 picoseconds.
- The default policy fingerprint is pinned and covers the thresholds; rule
  tables are read-only; looser rules must go through `Policy.allow`, which
  validates the correction name.
- Licence texts ship in the wheel and sdist; the sdist test suite passes.

### Added

- `dadl quickstart [example|setup.yaml]` reads the YAML; `dadl init`,
  `dadl simulate`, `dadl analyze` (own tags in `.npz`/`.csv`); `--version`;
  `python -m dadl`; errors exit 1 with a message.
- Results record `policy_fingerprint`, `overrides` (with reasons) and
  `window`; threshold labels for values within 1.96σ of a threshold and for
  zero-σ results at a boundary; CHSH gets threshold labels.
- Shared Rust/Python test vectors (`tests/vectors/core.json`), a 200-seed
  pull study of the CHSH error, an exact series dead-time test.
- CI: release-mode tests, clippy, MSRV 1.74, sdist/wheel install test,
  Dependabot, a release check that the tag matches both versions and is on
  main, and separate build and publish jobs.
- CHANGELOG, CITATION.cff, a validated example calibration artifact.

## 0.0.2

First published release: `dadl-quantum` on PyPI and `dadl` on crates.io.
