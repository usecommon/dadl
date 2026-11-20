# DADL

**Detector-Aware Data Layer for quantum optics experiments.**

DADL is an open-source starter kit for single-photon experiments. It is
designed to keep a versioned model of every detector channel (efficiency, dark
counts, dead time, afterpulsing, crosstalk, delay, jitter) and carry that
model's uncertainty into Bell S, visibility and fidelity, always reporting raw
and corrected values side by side and refusing corrections that would make a
result look better than the data supports.

**Read the overview:** [usecommon.github.io/dadl](https://usecommon.github.io/dadl/),
a short technical report on the problem DADL solves, how it is used in an
experiment, its guardrails and its roadmap.

> **Status: pre-alpha (0.1.1): phase 1 done (unreleased), phase 2 done, phase 3 in progress (M0-M4,
> unreleased).** Anything may change in 0.x. The
> canonical phase plan is in [docs/design.md](docs/design.md#mvp-scope-and-milestones).

## What works today (0.1.1, and unreleased on main)

- A seeded time-tag simulator with injected efficiency, dark counts,
  non-paralyzable dead time (optionally after a paralyzable stage), jitter,
  delay and cascading afterpulsing.
- **Calibration (phase 2):** dark rate, afterpulse
  probability and time constant, and dead time fitted from a dark acquisition,
  with 95% intervals whose coverage is checked in CI at a few simulated
  regimes (see "Known limits"), written as
  content-addressed artifacts:
  `dadl simulate chsh-spad --dark --out dark.npz`, then
  `dadl calibrate dark chsh-spad --tags dark.npz --out artifacts/`. From a
  pair-source acquisition, `dadl calibrate pairs` fits the relative delay and
  combined jitter, off-peak accidentals in the declared window and the
  heralding efficiency; `calibrate_dcr_temperature` fits DCR(T) (Arrhenius)
  from dark calibrations at several temperatures. `dadl calibrate deadtime`
  (unreleased) fits a series dead time (paralyzable, then non-paralyzable) from
  acquisitions under steady light at two or more rates
  (`dadl simulate chsh-spad --deadtime-scan --out scan.npz`), with
  afterpulsing held at the dark artifact (`--dark artifacts/`).
- `setup.yaml` drives everything: copy it with `dadl init`, edit it, run it.
- Raw CHSH S from fixed-grid coincidences with a planned window, from
  simulated tags or **your own tags** (`dadl analyze`); with `--calibrations`
  (unreleased), a corrected value from the forward model too.
- Bell fringe and HOM dip visibility, with accidental subtraction only where
  the guardrail policy allows it. `dadl quickstart hom-snspd` simulates a HOM
  dip on two SNSPDs and analyses it from the tags.
- **PicoQuant PTU files** (T2 and T3: PicoHarp, HydraHarp, TimeHarp 260,
  MultiHarp) read by DADL's own adapter (phase 2; T3 unreleased). T3 files
  come from a sync (usually pulsed) source. DADL models CW light: in a T3 file
  the sync is used as a clock only. Corrected values are refused on T3 input
  unless the setup declares `source.type: cw` (the templates' `spdc-cw` does
  not count). The T3 sync clock need not match between calibration and run
  for CW; record dark runs in T2 or T3. Pulsed sources are not
  supported on any input (T2, T3 or plain tags): pulsed accidentals sit in
  side peaks at multiples of the repetition period. DADL detects them
  through the off-peak windows (not flat) and refuses the
  accidental-subtracted values, with limited power (see Known limits).
- Guardrails that refuse invalid corrections, and record the policy
  fingerprint, overrides and window on every result.
- **Hardware twin** (phase 2): `dadl simulate SETUP --out DIR
  --format ptu` writes the simulation as PicoQuant PTU files (marked
  simulated), and `SimulatedTagger` emulates a live tagger with buffer
  overflows, dropped batches and config changes. See `adapters/README.md`.
- **Optional compiled core** (phase 2): with `dadl-quantum-core`
  installed, dead time, coincidences and histograms run in Rust (about 40-60x
  faster; about 5x for cross-correlation); `dadl.backend()` says which is in use.

Not yet: signing, the explain engine, MCP.

Raw values (S, visibility) report a statistical Gaussian 1σ error;
visibilities also always carry an exact 95% interval. Corrected values from
the forward model (unreleased; HOM and CHSH with `--calibrations`; fringe through
`fringe_from_tags(calibrations=)`) include the
calibration uncertainty, by joint Monte Carlo over counts and calibration
draws ([Methods](#methods-corrected-values)). They are secondary values: the headline is
always raw. The calibration part is conditional on the models: it includes the calibration
artifacts' statistical uncertainty only. It excludes each artifact's `uncertainty.not_included`
(copied into the record as `calibration_excludes`), the IRF shape beyond the Gaussian the delay
check accepted (HOM and fringe capture), efficiencies, afterpulsing, multi-pair emission, the
splitter ratio, drift between calibration and run below the drift threshold, and the choice
among several valid artifacts. A result keeps the labels of the raw value (`labels`) apart from
those of the corrected value (`corrected_labels`).
Calibrations report 95% intervals that need not be symmetric, validated by
simulation; each artifact's `uncertainty` block lists the components not
included (tagger timebase, temperature sensor, connectors, declared delays).

### Known limits

The calibrations went through AI-agent role-play reviews by simulated outside
experts; there has been no independent human review yet
([assumptions register](docs/reviews/2026-11-expert-assumptions.md),
[review cycles](docs/reviews/2026-11-review-cycles.md)). Their limits, stated
plainly:

- **CW sources only.** `source.type` must be `spdc-cw` (on T3 input, `cw`).
  Pulsed sources and pulsed or thermal pair statistics are not modelled, on
  any input; where the off-peak windows show pulsed light, accidental-subtracted
  values are refused. On T3 input the accidentals are read off-peak only when
  you declare `source.type: cw`. The template's `spdc-cw` does not count,
  because it is the default. With pulsed light, a repetition period that
  divides the off-peak step is not detected by the flatness check. The
  off-peak estimate is then right only if the first off-peak shift also falls
  on the comb. Otherwise V and |E| read low (conservative for the threshold),
  with no label. The off-peak flatness check finds pulsed accidentals from how
  unevenly they fill the off-peak windows. With about one off-peak
  coincidence per window, as in a typical CHSH run, it reliably finds them
  only when they are about half or more of one acquisition's off-peak counts.
  In council round 9's test it refused 30 of 30 runs at one half, 8-17 of 30
  at 30%, and 0-2 of 30 at 20%. A smaller pulsed share passes unlabelled. The
  off-peak windows then miss part of the same-pulse accidentals, so V and |E|
  read low (conservative for the threshold; council round 9 red team). Longer
  acquisitions raise the counts per window and the power. Independent CW
  multi-pair coincidences are handled, leaving a bias of about
  −μη_a(1 − η_b)/2 (−1.3% at μ = 0.1 pairs per slot).
- **Single-exponential afterpulsing.** Detectors with a multi-component or
  power-law afterpulse kernel get an effective single-exponential value
  between the fast component and the total, with an interval that covers
  neither (an SPCM-like three-component kernel: 0.38–0.41% for a 0.30% fast
  component and a 0.50% total), and the check does not catch it at typical
  lengths. The afterpulse probability belongs to the hold-off setting it was
  measured at. Free-running InGaAs/InP detectors (µs hold-off, time constants
  of tens of µs) are fitted when the time constant is below a tenth of the
  mean dark gap; when dark rate × time constant exceeds about 0.01 (a hot
  InGaAs detector), the afterpulse probability reads high and its interval
  covers less than 95% (72–75% at 0.04), and the calibration says so.
- **Crosstalk is simulated, not calibrated.** The simulator
  takes opt-in crosstalk links between channels (`dadl.sim.Crosstalk`:
  probability, offset, jitter; off by default, and not yet a `setup.yaml`
  field), but no calibration fits a crosstalk matrix. A corrected value is
  checked for it only with a delay-jitter artifact (the forward model's
  predictive check refuses crosstalk in the window); raw values and corrected
  values without one are not flagged. On the `hom-snspd` example, 1% zero-offset crosstalk from A to B
  lowers V from 0.919 to 0.904 (0.889 with 1% each way), unlabelled; at an
  offset outside the coincidence window (a breakdown flash several ns late)
  it leaves V unchanged. Twilight counts are not modelled.
- **Jitter is a Gaussian sigma.** Datasheets usually quote the FWHM
  (2.355 sigma): convert before entering `jitter_sigma_ps`.
- **Hard-edge dead time.** The dark calibration reports one non-paralyzable
  effective channel dead time (the longer of detector and tagger dead time;
  afterpulses inside it are not seen), valid for actively quenched SPAD
  modules. `calibrate deadtime` adds the series model, a paralyzable stage
  before the non-paralyzable one, from illuminated runs; both stages are
  still hard edges. The series chain is a phenomenological filter with hard
  edges. It describes a diode or front-end recovery (paralyzable) followed by
  a longer electronics or tagger dead time (non-paralyzable). It does not
  describe the reset-window paralysis of an actively quenched module (photons
  in hold-off do nothing; one just after it re-quenches silently:
  non-paralyzable then paralyzable), twilight counts, soft recovery of SNSPDs
  or the gradual efficiency recovery of passive quenching. With enough gaps
  such a detector shows as a misfit. With few gaps, or when the recovery is
  close to fully paralyzable, it can pass with an effective tau_p near tau_np
  that describes the scanned rates but not rates above them.
  The scan needs steady light: an attenuated CW laser with its intensity
  stable to about 1-2% over each run, not pulsed light or a mains-powered
  lamp. Slow intensity noise looks like less paralysis: a +-20% 100 Hz
  modulation biased tau_p 10% low with the predictive check passing; +-40%
  is flagged.
  tau_np is the lower edge of the dead time as recorded. A real tagger adds
  its own jitter after the dead time, so the shortest gap sits about
  4 sqrt(2) sigma_tagger below the detector's dead time (about 0.1 ns for
  20 ps rms); the interval does not cover that offset.
  Dead times are read on the tagger's tick. The series tau_np interval allows
  one tick either side (a recorded gap is within one tick of the true gap).
  The dark dead-time interval does not; it spans many ticks in any dark run,
  so the difference is below its resolution.
  An afterpulse still pending when the next gap starts (from the previous
  avalanche, Exp(tau_a) after the record by memorylessness, a fraction pi of
  gaps recorded as `afterpulse_pending_max`) is modelled to first order: one
  pending candidate, without its product with the record's own afterpulse.
  The remainder, about pi (p + pi) of gaps, biases tau_p low by an amount
  that does not shrink with data (about 0.1 sd per unit of remainder x
  sqrt(gaps)); above 2 a note says so. Checked against the simulator at
  6e6 cps (p 2%, tau_a 30 and 90 ns) and 1e7 cps (p 5%, tau_a 200 ns; pi
  6.8%, where the remainder shows at 1.3e7 gaps, chi-square p 0.05). In the simulator, avalanche candidates blocked by
  the dead time restart the paralyzable stage but neither afterpulse nor
  emit crosstalk.
  Each scan run is taken below the count-rate maximum (n tau_p < 1, the
  designed scan) unless the fit with each run on its better branch is better
  by more than 5 in NLL (a likelihood ratio of about 150); then those runs
  are fitted past it. Within that margin the interval spans both fits and a
  note names the runs. Below the maximum the high branch fits worse in
  expectation and wins by more than 3 in about 0.5% of scans (9 of 2,000);
  with margin 5 the false-switch rate is at most about 8e-4 for any
  separation (Cox's normal approximation for non-nested fits). Past it
  (1/3/60 Mcps) the margin is passed in 28 of 48 and the union covers the
  rest. Past the maximum, coverage is checked by probes only, not by a CI
  regime.
  The analyses' forward model uses the dark calibration's non-paralyzable
  dead time only. With paralysis the live fraction is overestimated: for
  tau_p close to tau_np by about 0.5% at a 10% dead fraction and 3-4% at the
  20% limit, per channel (about twice that for a pair); ratios such as S and
  V cancel much of it.
  The forward model takes each channel's non-paralyzable dead time from the
  dark artifact. At high analysed rates the analysed gaps measure it more
  precisely, and often show it to be shorter than the dark value (label
  `DEAD_TIME_ABOVE_GAPS`). The pair live fraction is first order: the term
  r_c·tau·(fa + fb) is not included. In the chsh-dead-time regime the two
  effects cancel: the corrected S covers at 95.5% with mean z −0.43
  statistical se. Neither effect is corrected alone, because each alone
  biases the value more (council round 8). Reading the dead time from the
  analysed gaps, together with the second-order term, is milestone M5a.
- **Afterpulse time constant near the detection threshold** covers somewhat
  under 95% conditional on detection.
- **DCR(T)** is Arrhenius only, valid inside the calibrated temperature range
  at one bias, for SPADs.
- **Rate dependence.** Delay, jitter and efficiency are valid at the singles
  rates recorded in the artifact only.
- **Heralding efficiency** needs a heralding run with no analysers in the
  path (`--no-analysers`); on Bell data it measures the analyser projection.
- **Bell results** assume fair sampling and fair coincidence, with static
  settings: they do not close detection, locality or freedom-of-choice
  loopholes (nor the memory loophole, where a local model uses earlier
  trials). A planned, fixed window does not by itself close the
  coincidence-time loophole for CHSH: a local model whose timing shifts are
  a little shorter than the window still reads S > 2 (2.34 in
  `test_published_results.py`); only windows much wider than any plausible
  setting-dependent delay make that model fail. "Coincidences per singles" is printed for each arm: the fraction
  of that arm's clicks that post-selection kept, about the other arm's
  efficiency times the window capture. Compare the smaller of the two with
  the detection-loophole thresholds, about 0.83 (CHSH) or 2/3 (Eberhard). CHSH needs both outputs of each analyser (A+, A−, B+, B−); the
  two-detector, 16-measurement scheme with rotating polarisers is not
  supported.
- **Unequal output efficiencies can raise S.** Raw counts carry each
  side's output imbalance m = (η₊ − η₋)/(η₊ + η₋) into E, shifting S by up to
  2·m_a·m_b, upward when both sides lean the same way: 0.6/0.4 on both sides
  read S = 2.03 from a source at 1.98. Results where the singles show this
  carry an "output imbalance" label; swap outputs or normalise per detector.
  Corrected S does not correct detector efficiencies: an output imbalance
  biases it as it does raw S (0.6/0.4 on both sides: corrected S 2.033 for a
  1.980 source, 13 sigma high). The "output imbalance" label applies to both
  values.
- **Delays that miss the pair peaks can raise S.** A delay error on some
  channels changes the window capture of some outcome pairs and not others:
  A+ and B+ 400 ps late, declared 0, read S = 2.20 for a 1.98 source (2.05 at
  200 ps). Each outcome pair's peak is fitted over the four settings; S
  reweighted by those captures is in the counts, and a result whose raw S
  differs from it by more than half a standard error and by more than 3.29
  standard deviations of the peak fits carries an "unequal window capture"
  label. Declare the delays from a pair calibration. A window narrow against
  the jitter, with the peaks on the declared delays, is labelled separately
  ("window narrow against the pair jitter").
- **Corrected CHSH** is secondary, never the Bell claim. Calibrated delays
  that no longer match the analysed peaks (a drift since calibration) are
  refused as `Misspecified`: the predictive check refits the delay, so it
  would not see the drift. The drift check refuses only a drift that moves S
  by more than half its statistical error and more than 3.72 standard
  deviations of the check's own noise (two-sided 2e-4) (the delay artifacts' capture draws and
  the analysed peak fit). A smaller drift is not refused (100 ps at 350 ps
  jitter and a 2 ns window: 0.41 se, 0/20 refused); one within the artifacts'
  noise is covered by the calibration part of the interval. Longer delay
  calibrations make the check sharper. A dead time known only as an upper limit (an SNSPD,
  or a short dark run) is used at the bound for the value and its interval;
  the counts give the value at zero dead time as well (the truth lies
  between if the bound holds; a bound set by a few short dark gaps may not),
  and a range wider than the error is labelled: quote the range.
  The predictive check fits a Gaussian peak. A timing tail of a few percent
  (SPAD diffusion, SNSPD) fails it on a Bell run's 16 histograms even when the
  delay artifacts pass (3% of clicks with a 300 ps tail: 6 of 8 runs refused),
  although a tail common to all channels does not bias S. Until the IRF tail
  (D8) lands, real-detector Bell runs will usually need
  `--allow misfit_forward_model=...`, which also stops the check flagging
  crosstalk.
- **Dead time can raise S.** A recorded click makes the partner arm's dead
  state depend on the outcome, so |E| grows with η²Rτ: a V = 0.7 source
  (S = 1.98) read S = 2.14 with a 1 µs dead time at 3·10⁶ pairs/s and
  η = 0.5 (2.35 at η = 0.8). Results
  where coincidence rate × dead time exceeds 1% carry a "dead-time selection"
  label; tags whose gaps are impossible for the declared dead time (often a
  ns export read as ps) are refused.
- **Coverage gates** run 1,000 simulated replicates per calibration at a few
  regimes and catch gross failures (coverage below about 93%), not small ones.
- **Gaussian timing response.** The delay fit models a Gaussian IRF (with the
  tagger tick's quantization modelled exactly) on a two-level background. A detector with an exponential timing tail (most SNSPDs,
  many SPADs) is marked `misfit` and the delay artifact refused; the other
  pair calibrations are checked on their own.
- **SNSPDs.** `family` is a label only. Recovery, timing tails, bias-current
  dependence and latching signatures are not modelled; dark data cannot
  identify an SNSPD's dead time (reported as an upper limit, marked
  unidentified); DCR(T) refuses runs whose bias condition is recorded as a
  current (a `bias*` key not in volts). A latch or other
  interruption in any acquisition is detected (longest gap against the count
  rate) and marks calibrations `misfit` and analyses "interrupted". Many short
  dropouts (microseconds each, as from a self-resetting latch) are not caught
  by this single-gap test: 14% of a signal lost in 5–20 µs pieces read
  heralding 0.627 for 0.729, status ok. A source
  modulated by ±50% or more can trip the same test, as can free-running
  InGaAs detectors whose afterpulse time constant is comparable to the
  median gap (about 5x the nominal false-alarm rate at 1 µs dead time and
  τa = 10 µs). The predictive-check p-values are approximate (the delay
  check's runs high with sparse peaks); use them for trends, not as exact
  probabilities. Dark-rate drift within
  ±2% between the halves of a run is tolerated and recorded.
- **Heralding efficiency is grid-specific.** It counts pairs in the same fixed
  slot, so it includes the window capture (about 1 − E|Δt|/W) and the signal
  arm's live fraction given a recorded herald, and it is diluted by the
  herald channel's dark counts and afterpulses (20 kcps of herald darks read
  0.095 for a coupling of 0.2). `calibrate pairs` prints the capture implied
  by the fitted jitter when the delay fit passes its check.
- **Off-peak checks have limited power.** A side-window excess smaller than
  about twice the accidentals interval can pass the check unseen. Analyses
  that subtract accidentals check each acquisition's off-peak windows for
  flatness (`offpeak_p` in the counts) and refuse values whose windows are
  not flat (pulsed light, a periodic background), overridable as
  `misfit_forward_model`.
- **Corrected visibility** by subtraction (no calibrations) is refused when
  the accidentals exceed the minimum count (or the corrected value would
  exceed 1), which removes the datasets nearest V = 1; report the raw value
  with its exact interval there. The forward model keeps the unclipped value
  and labels it.
- **HOM** is two-fold and two-point: the visibility at the declared dip delay
  against a reference delay, on a 50:50 beamsplitter, with detectors that do
  not resolve photon number. Four-fold heralded HOM, delay-scan fits and
  splitter imbalance (which caps V at 2RT/(R² + T²)) are not handled. The
  reference delay must be small against the window. Acquisitions of different
  length are compared as rates; singles cannot normalise drift (they fall at
  the dip by design).
- **Provenance in the Python API.** The CLI labels simulated input and records
  input hashes; `read_tags` followed by `chsh` or `hom_from_tags` in Python
  does not (it returns arrays, not their provenance).
- **Hardware twin.** PicoHarp 300 and TimeHarp 260 P presets have two inputs,
  so they can run the two-channel HOM example but not chsh-spad. The twin
  puts the first detector on input 0, the sync input of HydraHarp and
  MultiHarp devices; recordings usually start detectors at input 1.
- **Whole-file PTU reading.** A PTU file is decoded in memory at about 100
  bytes per record in T2 (215 with many markers): a 4 GB file needs about
  100 GB. T3 needs up to about 170 bytes per record (it adds a sync tag per
  occupied sync period).
  Chunked reading is planned. Unrecognised special records and a divided sync
  input mapped as a detector are refused.
- **Analyses hold every channel's tags in memory**: about 8 bytes per tag,
  plus up to the same again for working copies. A 10-minute CHSH run at
  1 Mcps per channel (6e8 tags per channel, 16 channel-settings) needs about
  80-150 GB, more than a laptop has. Analyse shorter acquisitions or a subset
  until chunked analysis lands (milestone M9). Series dead-time calibration is
  not limited this way.
- **Tagger settings** (input offsets, trigger levels, dead-time settings, the
  tick) are recorded in each calibration's `conditions` from PTU headers.
  Tagger settings (`tagger_*`) and bias are compared when the calibration is
  used: a mismatch, or a setting recorded but not given for the dataset, is
  refused (overridable). The T2 tick and the T3 sync period share one header
  key; it is compared only between files of the same mode (`tagger_mode`).
  On T3 input under `source.type: cw`, a changed sync clock is labelled for
  dark and pair artifacts instead of refused.
- **Double clicks** in CHSH (both outputs of one side in a slot) count in two
  outcome combinations, which pulls E toward 0.
- **Reproducibility.** The same tags, setup and software stack give the same
  artifact ID when `created` is fixed (`--created`); posterior draws are
  seeded from the setup's `seed` mixed with the input hashes. Across
  platforms or NumPy builds, floating-point summation order can change the
  last digits, and so the ID (seen once in CI on Linux with NumPy 1.24).
- **Records, not certificates.** An artifact names the tagger, input and
  detector serial, the input file and setup hashes, the acquisition time
  and the conditions you give (`--condition temperature_K=233`), but it is
  not a calibration record in the ISO/IEC 17025 sense: uncertainty is
  statistical only, and traceability is recorded, not enforced. Validity is
  enforced (from acquisition, or creation without one, to `--valid-until`;
  `valid_for.conditions` too, as the top-level conditions are);
  an artifact with no validity range is labelled; supersession
  (`--supersedes`) is recorded, not enforced. The ID gives
  integrity, not authenticity: signing is planned. Parameters reported only
  as a bound have `value: null` and `bound: "upper"`, and draw columns that
  are stand-ins rather than a posterior are listed in `posterior.stand_in`.
  Artifacts written by 0.1.x start their validity at creation
  (`valid_for.from` = created); this branch starts it at acquisition. For
  data taken between an old artifact's acquisition and its creation, pass
  the validity override or rewrite the artifact with `--supersedes`. Until
  the next release the package still reports version 0.1.1.
- **Interactive use, not unattended operation yet.** Setups have exactly the
  four CHSH roles (or two HOM ones), so an 8-input tagger is calibrated as two
  setups; output is text, with no `--json`; nothing compares a calibration
  with earlier ones, finds the latest artifact for a detector, or skips
  rewriting an identical record (pass `--created` to reproduce one). Long dark
  runs with more than ±2% drift are refused as misfits and need splitting.
  CSV tag files take about 500 bytes per row in memory; a low-count dark fit
  peaks at about 0.5–0.75 GB. For trends, use the `p_shape`, `p_drift` and
  `p_gap` conditions rather than the combined `ppc_pvalue`. Tags must be in
  time order per channel in the file: a tagger clock reset is refused, not
  sorted.
- A failed predictive check marks the calibration `misfit`; the CLI refuses to
  write it unless `--accept-misfit` is given.

### Refusals and their rates

On clean data within the tested regimes, DADL refuses an analysis at analysis time with
probability at most 1 in 1,000, all checks together: the forward-model check 4e-4, the CHSH
delay-drift check 2e-4, the off-peak flatness check 2e-4, and the dead-time gap tests 2e-4. Each
check corrects for its own number of channels or acquisitions. Separately, a calibration whose
model check fails (p < 0.001 per artifact) is not written; a CHSH calibration set of eight
artifacts meets this in under 1% of sets, and the remedy is to re-take that calibration run.
Labels are not refusals: some fire on a few percent of clean analyses, and their text says how
often.

## What DADL will not do

- Make low-cost detectors "research-grade". Lost efficiency and jitter cannot
  be recovered in software.
- Close the detection loophole in software. Efficiency correction is
  equivalent to assuming fair sampling.
- Correct data in a QKD security path, or detect attacks on detectors. DADL
  has no blinding, after-gate, time-shift or efficiency-mismatch detection
  and does not support gated detectors; nothing it produces belongs in a QKD
  security argument.
- Verify, in phase 1, that a window declared as planned really was fixed
  before acquisition, or that the tags you pass in are unfiltered. It records
  what you declare; signed analysis plans arrive in phase 3.

## Quickstart

Requires Python 3.10 or later. (On macOS the system `python3` is 3.9, and pip
then reports "No matching distribution found for dadl-quantum"; install a
newer Python first.) The PyPI package is `dadl-quantum`; the import name and
the command are both `dadl`.

```bash
pip install dadl-quantum
dadl quickstart
```

This loads the bundled `chsh-spad` setup: a polarization-entangled SPDC
source on four detectors with known imperfections. It prints raw S with a
fixed, pre-declared window next to the closed-form expectation.

To change the experiment, copy the setup, edit it and run your copy:

```bash
dadl init my-experiment
dadl quickstart my-experiment/setup.yaml
```

## Using your own data

Export your time tags as a table with three integer columns, in `.npz`
(arrays named `setting`, `channel`, `timestamp_ps`) or `.csv` (header
`setting,channel,timestamp_ps`):

| Column | Meaning |
| --- | --- |
| `setting` | CHSH setting pair, 0–3, in the order (a,b), (a,b′), (a′,b), (a′,b′) |
| `channel` | Detector input number, mapped to A+, A−, B+, B− by `channels:` in `setup.yaml` |
| `timestamp_ps` | int64 time tag in picoseconds |

Set the window width and per-channel delays in `setup.yaml` *before* looking
at the data, then:

```bash
dadl analyze my-experiment/setup.yaml --tags tags.npz
```

`dadl simulate my-experiment/setup.yaml --out tags.npz` writes simulated tags
in the same format, so you can check the round trip first.

For a HOM setup, `--calibrations artifacts/` (artifacts written by `dadl calibrate`)
adds a corrected visibility from the forward model, with calibration uncertainty and the
artifact IDs. Dead-time and window-capture corrections (the latter from a
delay-jitter artifact, `dadl calibrate pairs`) are explicit policy changes for HOM:
`--allow dead_time="calibrated from separate dark runs"`, `--allow window_capture="..."`,
and the result says so.
With a delay-jitter artifact the forward model's peak shape is also checked against the
analysed coincidences; something in the window it lacks (crosstalk, a misdeclared reference
delay) is refused as `Misspecified` rather than corrected with an interval that does not know it.

For a CHSH setup, `--calibrations` adds a corrected S as a secondary value; the headline stays
raw. The default CHSH rule allows no correction, so each is an explicit policy change:
`--allow dead_time=... --allow accidental_subtraction=... --allow window_capture=...`. Dark artifacts must cover
every channel and delay-jitter artifacts all four channel pairs. Window capture for CHSH needs a
delay-jitter artifact for each of the four channel pairs, so four heralding runs (re-route the
fibres between them). `dadl simulate --heralding` writes the A+/B+ pair only.

One artifact per kind and channel set (in any order): with calibrations before and after the run,
pass the one to use; the drift between them is not propagated. An artifact recorded for another
detector serial than the setup names (also through the Python API) is refused as outside its
validity, and a dark dead time the analysed gaps contradict is refused on its own (both
overridable).

Refusals name their override; on the CLI each is `--allow NAME=REASON`, and the reason is
recorded. The names:

| Override | What it accepts |
| --- | --- |
| `calibration_outside_validity` | an artifact outside its validity range, conditions, detector or tagger |
| `dead_time_contradicted` | a dark dead time the analysed gaps contradict (labelled with the problems, recorded as `dead_time_gap_problems`) |
| `misfit_forward_model` | a failed predictive check, a CHSH delay drift, or off-peak windows that are not flat |
| `misfit_calibration` | an artifact whose own check failed |

Conditions: an artifact's `temperature_K` must be within 0.5 K of the dataset's, and recording
`temperature_K` at calibration (`--condition temperature_K=...`) makes it required at analysis.
`bias*` and `tagger_*` values must be equal, and a `<name>_<unit>` must lie within a recorded
`<name>_min_<unit>` / `<name>_max_<unit>` (one-sided when only one is recorded). Other recorded
names are labelled "recorded, not compared".

The dead-time check refuses a dark dead time that the analysed gaps contradict in either
direction. Downward: gaps shorter than its floor, point − ln(K / 1e-4) / (N f), with N the dark
run's gaps, f the lower end of its dark rate and K the analysed channels (the shortest of N dark
gaps sits above the dead time by Exp(N f)). That is about 4 ns below the point for a 1 s dark
run at 50 kcps and 265 ns for 0.1 s at 20 kcps, where the dark data cannot exclude anything;
artifacts without a recorded gap count keep the older floor, three lower half-widths below the
point. Upward: with j = max(10, 0.1% of the gaps) and M = max(100, j + 1), the gap j + 1 from the
shortest more than (j + 1) F⁻¹(1 − 1e-4/K; 2(j + 1), 2M) mean edge spacings above its upper end
(the exact law of an exponential edge; 100((K / 1e-4)^(1/100) − 1) at j = 0). Each direction
refuses a clean analysis with chance at most 1e-4 (far less upward: the upper end is the shortest
dark gap). A dead time between the floor and the truth passes; at low analysed rates (about
1,000 gaps) the upward test has no power, and at 10 kcps and 50,000 gaps it misses a 22 ns
artifact on a 50 ns detector.

The analysis tests the dark dead time against the analysed gaps in both directions. Each test
tolerates a few unusual gaps, up to 0.1% of a channel's gaps or 10, whichever is more. A rarer
double trigger (a cable reflection, a ringing discriminator) cannot hide the edge. A more frequent
one makes the edge look shorter than the detector's dead time. If you declare the detector's dead
time in the setup (`detectors.<name>.dead_time_ps`), DADL labels a dark dead time far below it. To
be sure, look at the dark artifact's notes, which say when a single short gap set its bound.
When the floor reaches below zero, nothing is refused downward; the 95% lower end is then tested
for a label only (`gap_check_floor_kind` "lower end (label only)"): it misses in 5% of channels
by design. A dark dead time longer than many analysed gaps is labelled
(`DEAD_TIME_ABOVE_GAPS`; see Known limits).

An identified dark dead time can be wrong when the detector's dead time is far below the
spacing of dark counts (an SNSPD at 100 cps: 3-7% of dark runs reports a microsecond-scale dead
time as identified, and a longer dark run does not help). The analysis refuses it when it
implies a dead fraction above 20% at the analysed rate, and labels it when the analysed gaps
fall below its 95% lower end ("the analysed gaps put the dead time ... below the dark artifact's
95% lower end"); calibrate an SNSPD's recovery from an illuminated inter-arrival histogram. A
dead time, identified or a bound, whose dead fraction at the analysed rate exceeds 20% is
refused: there the pair term and the unmodelled paralyzable stage are not trusted, and a dark
dead time this long is often a false identification. Correcting with a series dead-time
calibration is not yet supported. To correct the rest, move the dark artifacts out
of the calibrations folder.

CHSH results, raw and corrected, add these keys to `counts`:

| Key | What |
| --- | --- |
| `capture_offsets_ps` | each outcome pair's peak offset from its declared delay |
| `capture_reweighted_s` | S reweighted by the window captures from the peak fits |
| `capture_sd_data` | standard deviation of that S, from the peak fits |

Corrected results add these:

| Key | What |
| --- | --- |
| `delay_drift_diff`, `delay_drift_z` | analysed minus calibrated capture-reweighted S, and over its noise (CHSH) |
| `dataset` | the dataset context the calibrations were checked against |
| `forward_check_p` | the predictive check's p-value; absent without delay-jitter artifacts |
| `forward_check_note` | why the predictive check was skipped (T3 input with no integer lattice); the CLI prints it |
| `offpeak_p` | the off-peak flatness check's p-value, Bonferroni over the acquisitions (also on HOM and fringe values with accidentals subtracted) |
| `offpeak_draws` | the Monte Carlo draws the sequential off-peak dispersion test used (seeded from the counts) |
| `dead_time_gap_problems` | the dead-time gap tests' problems, when accepted with `dead_time_contradicted` |
| `monte_carlo` | the Monte Carlo `seed` and `n_draws` the corrected interval depends on |
| `min_gap_ps` | with dark artifacts: the shortest analysed gap per channel |
| `gap_check_floor_ps` | with dark artifacts: the dead-time floor per channel; None when not checked downward |
| `gap_check_floor_kind` | "exponential law" or "3 widths" (refusing), or "lower end (label only)", per channel; None when not checked downward |
| `calibration_excludes` | the used artifacts' `uncertainty.not_included`: what the calibration part leaves out |

### Methods (corrected values)

What a methods section needs for a corrected HOM, fringe or CHSH value (details in
`dadl/forward.py`):

- **Model.** Per acquisition and outcome pair, the signal rate is s = (C - A) / (T L capture):
  C the coincidences in the declared fixed-grid window, A the accidentals from the same
  acquisition's off-peak bins, T the duration, L = (1 - n_a tau_a)(1 - n_b tau_b)(1 + r_c
  min(tau_a, tau_b)) the live fraction of the pair (n the recorded singles rates, tau the dark
  artifacts' dead times, r_c the recorded coincidence rate), and capture the fixed grid's
  capture E[max(0, 1 - |D| / W)] for a Gaussian arrival difference D at the fitted delay and
  jitter of the delay-jitter artifact. HOM V = 1 - s_dip / s_ref; fringe V = (s_max - s_min) /
  (s_max + s_min); CHSH E is the normalised signed sum of the four outcome signals and S the
  signed sum of the E.
- **Draws.** 4,000 joint Monte Carlo draws: counts and off-peak sums from Jeffreys posteriors
  Gamma(c + 1/2); calibration parameters from the artifacts' stored draws (joint within an
  artifact, independent across artifacts). A parameter known only as an upper limit is used at
  the bound, and the value at zero is reported too.
- **Point and interval.** The point is the plug-in value at the counts and the artifacts' point
  values, not the Monte Carlo centre; the interval is the 2.5% and 97.5% quantiles of the joint
  draws, so it need not be centred on the point (for dead time the point sits about 0.4 sd low).
- **Parts.** The statistical and calibration parts are standard deviations with the other held
  at its point value; they need not add in quadrature to the total.
- **Citing.** Quote the Monte Carlo `seed` and `n_draws` from `counts["monte_carlo"]` and the
  artifact IDs. The Monte Carlo noise of the interval ends is about 0.1 statistical sd.
- **Calibration draws.** Stored draws are posteriors under flat priors in the fit's coordinates
  (log rate, logit p, log tau_a) inside the fit's box; rank checks confirm this (round 6). Near
  p = 0 that prior is improper without the box; M5 declares proper priors.
- **What it leaves out.** The calibration part is conditional on the models, as above
  (`calibration_excludes` lists the artifacts' own exclusions).

From Python:

```python
import numpy as np
import dadl

setup = dadl.load_setup("my-experiment/setup.yaml")
runs = dadl.read_tags("tags.npz", setup.channels, setup.setting_keys)   # or build the dict yourself:
# runs[("a", "b")] = {"A+": np.int64 array of ps, "A-": ..., "B+": ..., "B-": ...}
result = dadl.chsh(runs, setup.window, settings_deg=setup.settings_deg)  # the planned window and delays
print(result)               # raw S, labels
print(result.as_record())   # flat record with labels, policy fingerprint, window

print(dadl.bell_visibility(1450, 550))            # raw fringe visibility
print(dadl.hom_visibility(1000, 120, accidentals=40))  # dip visibility, corrected value labelled
```

PicoQuant `.ptu` files (T2 or T3 mode) read directly: `from dadl.adapters import read_ptu`, then
`read_ptu("run.ptu")` returns int64 ps per channel, with rollovers and markers reported as events
(see `adapters/README.md`). From the command line, pass a folder of `setting0.ptu`,
`setting1.ptu`, ... as `--tags`. On HydraHarp, MultiHarp and TimeHarp 260 files, input 0 is the
sync input and detector inputs are numbered from 1, as in PicoQuant's software; map those
numbers under `channels:`. T3 times are rebuilt as sync count x sync period + start-stop time,
assuming the header's sync period held for the whole run, and in T3 input 0 holds only the syncs
that started a recorded measurement. Other formats: export to the generic tag file.

Before calibrating recorded data, edit the copied `setup.yaml`:

```yaml
backend: file                 # not `sim`: calibrations of recorded data are refused under `sim`
hardware:
  tagger: HydraHarp-1234567   # the tagger's serial, recorded as the channel's device
  detectors: {A+: SPCM-0042}  # detector serials: required for calibrations of recorded data
```

The fixed coincidence grid starts at t = 0 of the tag file: slot k holds
times t with k·W ≤ t − delay < (k + 1)·W, and the record says so
(`grid_origin_ps`). One channel must have delay 0, so a common offset cannot
move the grid. Shifting all tags (for example re-exporting from the first
tag) still moves pairs across slot edges: each shift changes S within its
statistical error, but choosing the best of ten shifts doubles the chance of
a spurious S > 2 + 1.96σ. Keep the export that was planned.

## Developing DADL

From a clone of this repository:

```bash
pip install -e "./python[test]"
pytest python/tests
cargo test --workspace
pip install ./bindings/python    # optional: the compiled core (needs Rust)
```

The Rust core is pinned to the toolchain in `rust-toolchain.toml`; its
minimum supported Rust version is 1.74. The pure-Python package is the
reference implementation; `bindings/python` (PyO3, built with maturin) puts
the Rust core behind the same functions as the optional package
`dadl-quantum-core`. Both backends are held to the shared cases in
`tests/vectors/core.json`, and `python bindings/python/bench.py` compares
their throughput (about 40-60x for dead time, coincidences and inter-arrival
histograms; about 5x for the cross-correlation the calibrations use).

## Repository layout

| Path | Contents |
| --- | --- |
| `core/` | Rust core (`dadl` crate): dead time, fixed-grid coincidences, histograms |
| `bindings/python/` | PyO3 bindings: the optional `dadl-quantum-core` package (module `dadl_core`) |
| `python/` | Python package `dadl`: simulator, analyses, guardrails, config, tag files, CLI, tests |
| `tests/vectors/` | Test cases shared by the Rust and Python suites |
| `schema/` | JSON Schemas for calibration artifacts (`calibration/1`, frozen at 0.1.0) and `setup.yaml` (`setup-1`), with a validated example |
| `examples/` | Example setups, starting with `chsh-spad` (also bundled in the package) |
| `adapters/` | Vendor adapters (phase 2 onward) |
| `mcp/` | MCP agent surface (v0.2+, seeking contributors) |
| `docs/` | Design doc, review records and the GitHub Pages site (`docs/index.html`) |

## Citing

See [CITATION.cff](CITATION.cff). Please cite the version you used and keep
DADL's labels with any number you publish.

## Licence

Dual-licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT licence ([LICENSE-MIT](LICENSE-MIT))

at your option. Unless you state otherwise, any contribution you submit for
inclusion is dual-licensed as above, without additional terms.
