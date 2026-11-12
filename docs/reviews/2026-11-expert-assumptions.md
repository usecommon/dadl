# Fifth review: five outside-expert reviewers on the assumptions (2026-11-12)

Before 0.1.0, five reviewers each examined one field's assumptions, read-only and independently.
Each started with no knowledge of the project or of the earlier reviews' reasoning. They are AI
agents playing outside experts, not people; a human domain review is still the top open item.

| Reviewer | Scope |
| --- | --- |
| Detector physicist (SPAD, SNSPD) | Dead time, afterpulsing, dark counts, timing response, DCR(T) |
| Statistician | Intervals, posteriors, predictive checks, the coverage gates |
| Metrologist (NMI detector calibration) | Measurands, units, uncertainty, traceability, records |
| Time-tagging engineer | PTU reader and writer, ticks, data loss, the hardware twin |
| Quantum optics experimentalist | Sources, HOM and CHSH models, heralding, guardrails |

Each assumption was rated **sound**, **acceptable with stated limits**, **questionable** (needs
evidence) or **wrong** (wrong numbers or misleading in realistic use). The reviewers backed their
ratings with their own simulations through the package's functions, about 40 scripts in all. The
register below merges findings that several reviewers made, keeping the strictest rating. Line
references point to commit d587f4b.

## Verdicts

| Reviewer | Release 0.1.0? |
| --- | --- |
| Detector physicist | Yes as pre-alpha, after making predictive-check failures loud and correcting design.md:189, with six public caveats |
| Statistician | Yes as pre-alpha, with four public caveats; fix the accidentals guard first if cheap |
| Metrologist | Not with heralding efficiency in its current form; fix channel identity and record delays |
| Time-tagging engineer | Yes, after refusing or flagging a non-zero PTU stop reason and taking duration from the data |
| Quantum optics experimentalist | Yes as a CW-only tool, after three blocking fixes and public caveats |

No reviewer would release it unchanged. All five agree on releasing as a pre-alpha after a short
list of fixes, with the limits stated publicly.

## Rated wrong (fix before release)

| # | Assumption | Reviewers | Evidence | Fix |
| --- | --- | --- | --- | --- |
| W1 | Running `calibrate pairs` on Bell (CHSH) data gives the signal arm's system efficiency; the quickstart suggests exactly this | Metrologist, quantum optics | 0.06–0.33 across the four settings for detectors of true efficiency 0.5, tight intervals, label "detector x coupling" | Drop the suggestion; refuse CHSH and HOM setups (or require an explicit no-analyser acknowledgement); record analysers in `conditions` |
| W2 | (C − λ)/S is unbiased for heralding efficiency | Metrologist (algebra confirmed) | Coincidences are counted per slot, so (C − λ)/S → η(1 − p). Truth outside the interval at p ≥ 1%. The gate's truth was the estimator's own expectation at p ≈ 2e-5 | η = (C/S − p)/(1 − p) with p = λ/S, propagated through the posterior; a high-accidental gate against the true η |
| W3 | The 200 ns guard keeps off-peak windows free of afterpulse structure | Statistician, detector physicist | Afterpulses following true pairs fill the side windows: +9% at the package's own Si-SPAD parameters, +61% to +165% with slower or stronger afterpulsing; the dispersion check misses it. Over-subtraction inflates the corrected HOM visibility | Tie the guard to the channel's afterpulse time constant, or place side windows tens of µs out; add afterpulsing to a pair-gate regime |
| W4 | A failing predictive check needs no action | Detector physicist, quantum optics | Artifacts are written and printed as normal at p = 0 (latch, IRF tail, wrong DCR(T) model, pulsed data with CW offsets), contrary to design.md:436 ("accepted only after posterior predictive checks pass") | CLI warning, a `misfit` status in the artifact, and refusal where the design says so |
| W5 | PTU files are complete unless truncated; duration is `MeasDesc_AcquisitionTime` | Time-tagging engineer | `TTResult_StopReason` codes for FIFO overrun, transfer failure and lost clocks are passed on as a neutral event and then dropped; the requested acquisition time gave a rate 0.0028 cps instead of 1000 cps | Map loss codes to typed events, mark `incomplete`, warn or refuse in the CLI; take the duration from the data |
| W6 | `source.type` selects the source model | Quantum optics | Any value, including "pulsed", silently gives the CW simulation | Validate it (CW only in 0.1.0) |
| W7 | The artifact's channel `device` identifies the physical detector | Metrologist | It holds the role label ("A+"), contradicting the schema; records cannot follow hardware through a swap | Record device and tagger serials separately from role labels; require them for hardware artifacts |
| W8 | Delta-method Gaussian errors suit Bell and HOM visibility at any count | Statistician | Coverage 63–93% near V = 1 or at low counts; threshold "crossings" without the near-threshold label 2–4x more often than 2.5% | Exact (Clopper–Pearson-based) intervals and a minimum count |
| W9 | The shortest dark gap is "detector plus TDC" dead time; the afterpulse probability is the detector's | Time-tagging engineer, detector physicist | It measures the longer of the two: about 80–95 ns on HydraHarp and PicoHarp, hiding about 65–70% of afterpulses with no warning; one spurious double trigger collapses it to 2 ns | Name it "effective channel dead time", record the tagger, state that p is conditional on it; make the estimate robust to a single outlier |
| W10 | The stored posterior draws match the reported intervals | Statistician | Below 5,000 gaps the Laplace draws for p and τa disagree with the profile intervals (draw coverage as low as 73%); latent until phase 3 reads the draws | Draws consistent with the intervals, and a test that their quantiles agree |
| W11 | `check_accounting` proves nothing was lost | Time-tagging engineer | Passes with lost-plus-duplicated, shifted, swapped or extra data; needs simulator truth no real device has | Compare times, not counts; restate the adapter rule for real devices as record accounting plus vendor loss flags |
| W12 | Twin presets write files a real device could produce | Time-tagging engineer | PicoHarp 300 and TimeHarp presets write four inputs (the PicoHarp 300 has two); "5 ps for a MultiHarp" holds only for some models (150 N: 80 ps) | Enforce input counts per preset; add an 80 ps preset; fix the docs |
| W13 | README: "uncertainties are statistical only (Gaussian, 1σ)" | Metrologist | Calibrations report non-Gaussian 95% intervals | Separate wording for analyses (1σ) and calibrations (95%, statistical, simulation-validated) |

## Rated questionable (needs evidence or a clear statement)

| # | Assumption | Reviewers | Evidence | Needed |
| --- | --- | --- | --- | --- |
| Q1 | Single-exponential afterpulse kernel | Physicist, statistician | SPCM-like two-component kernel: p 0.47% against 0.70%, 0% coverage, check silent at typical lengths; CMOS-like: dark rate +6% | Sum of exponentials, or report p as "fast component"; realistic misspecification test |
| Q2 | Hard-edge non-paralyzable dead time for SNSPDs and passive-quench SPADs | Physicist | Soft recovery: estimate is neither onset nor effective dead time, narrow interval; rate corrections −1.4% to −4.1% at Mcps | Restrict to actively quenched modules or fit recovery |
| Q3 | Multi-pair emission and pulsed sources | Quantum optics, metrologist | η 0.470 at μ = 0.1, 0.414 at μ = 0.3 (true 0.5), tight intervals; grid not locked to the laser: η 0.25 with pulses on slot edges | Pulsed source simulator, sync-locked slots, pump-power extrapolation |
| Q4 | Fixed-grid Bell analysis closes the coincidence-time loophole | Quantum optics | A local model with every photon recorded reaches S = 2.825 under the default policy, labelled only "fair sampling assumed"; delays can be fitted from the analysed data | Report singles and C/singles with S; "fair coincidence assumed"; refuse delays calibrated on the analysed tags; state static settings and open locality loopholes |
| Q5 | Two-point HOM visibility at the declared delay | Quantum optics | Dip off by 0.3 ps: V 0.882 against 0.919; ±5% pump drift: ±3σ | Scan fit, or label as visibility at the declared delay; normalise by singles |
| Q6 | Tick from the GCD of tags | Time-tagging engineer | A non-tick offset collapses it (−5 ps delay bias); 80 ps ticks fail the delay fit's own check | Take the tick from the adapter; GCD as cross-check only |
| Q7 | τa interval and detection switch | Statistician | τa covers 84–92% conditional on detection near the threshold, even at 20,000 gaps | Key the profile switch on afterpulse count; unified interval or stated conditional coverage |
| Q8 | DCR(T): exact temperatures, constant bias, Arrhenius only, any detector family | Physicist, metrologist, statistician | Narrow scans extrapolate about 3x low at 243 K; 0.5 K sensor offset ≈ 5%; inputs near the 20-gap minimum cover only 70–73% | Record T range, bias and sensor; refuse extrapolation and SNSPDs; warn near the minimum |
| Q9 | Intervals cover statistics only; Type B components omitted silently | Metrologist | Tagger timebase, sensor, connector, declared delays absent; no field says so | An `uncertainty` field listing included and omitted components |
| Q10 | Records support reuse and validity | Metrologist | Delays, acquisition time, singles rates, conditions not recorded; `valid_for` and `supersedes` unused and unenforced; site says "signed" | Record them; implement or remove the refusal claim; "content-addressed, signing planned" |
| Q11 | The coverage gates support "coverage checked in CI" | Statistician | Gates catch only gross failures (≤ about 91%) and sit in favourable regimes; the rule was widened once after a result; an independent seed passes all 12 | Qualify the README claim; gate weak-afterpulse and afterpulsing-pair regimes; freeze the rule, pool rotating seeds |
| Q12 | Stationary dark counts, uninterrupted acquisitions | Physicist | ±20% drift: check flags 20%, false afterpulse detection 5%; a 6 s latch: rate 91.5 against 100 | Split-half rate check; refuse long gaps |
| Q13 | Timing and efficiency independent of count rate | Physicist | Count-rate timing shifts (SPCM) and SNSPD reset effects | Record singles rates; valid at the calibration rate only |
| Q14 | Simulated-file marker survives | Time-tagging engineer | Lost in mixed folders, if one header entry is removed, in npz/csv exports | Also use `CreatorSW_Name`; refuse mixed folders; carry the flag |

## Rated sound

The register of things the reviewers found sound is short but substantive:
- PTU record types, wrap values, channel numbering and tick use match PicoQuant's demo code.
- The dead-time interval's tick law and the conservative predictive checks in the dark fit.
- The off-peak window definition matches the analysis grid.
- Content addressing gives integrity.
- Bunched photons on non-number-resolving detectors are modelled correctly.
- The guardrails do refuse accidental subtraction in Bell claims.
- Per-gate false-failure rates are controlled (family-wise 0.89% across 12 gates).
- All 12 gates pass on the statistician's independent seed.

## Accepted with stated limits

Accepted on condition the limit is stated publicly (several are not yet):
- Gaussian IRF.
- Werner-state CHSH model with the same visibility in every basis.
- Non-paralyzable dead time for actively quenched SPAD modules.
- Whole-file PTU reading.
- Jeffreys intervals (with coverage dips to about 90% at some means).
- Two-stage DCR(T) fit at normal counts.
- No crosstalk or twilight counts.
- 50:50 beamsplitter.
- Dark-count dilution of heralds.
- The design's series dead-time model is not yet implemented (design.md:189 must say so).
