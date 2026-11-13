# Review cycles after the outside-expert review (2026-11-12)

After the fixes from the [outside-expert review of the assumptions](2026-11-expert-assumptions.md),
ten further review cycles ran before 0.1.0. Each cycle used a fresh panel of reviewers, with no
knowledge of the project or of earlier cycles, each through a different expert lens. Reviewers
worked read-only and backed findings with their own scripts through the package's functions.
Findings were fixed, or stated as limits, before the next cycle started.

Reviewers are AI agents playing outside experts, not people; a human domain review is still the
top open item.

Severity: **wrong** (wrong numbers or misleading in realistic use), **questionable** (needs
evidence or a stated limit), **minor** (documentation, polish).

## Cycle 1: SNSPD specialist, frequentist statistician, first-time user

| # | Severity | Finding | Status |
| --- | --- | --- | --- |
| 1.1 | Wrong | A latch during a pair or HOM acquisition went undetected; a 20% latch in the HOM dip run raised V by 13σ; a 10% latch on a heralding signal read η 0.657 against 0.730, status ok | Fixed: interruption test (longest gap against the exponential tail of gaps above the median, robust to dead time and afterpulsing) in delay, accidentals and heralding (misfit) and in CHSH and HOM ("interrupted" label). False alarms 0–0.5% at p < 0.001 across the gate regimes; a 2% latch always caught |
| 1.2 | Wrong | `calibrate pairs` stopped at the first misfit; any SNSPD timing tail made the delay fit a misfit, so no accidentals or heralding artifact was written | Fixed: each artifact is checked and refused on its own. Gaussian IRF stated as a limit; a tail model is open |
| 1.3 | Questionable | Dark-only dead time of an SNSPD read 23x the true value, not marked unidentified | Fixed: reported as an upper limit only and marked unidentified |
| 1.4 | Questionable | DCR(T) accepted current-biased (SNSPD) runs; only `bias_V` was compared | Fixed: every `bias*` key compared; bias currents refused |
| 1.5 | Minor | `family` selects nothing | Documented as a label (schema, examples, README) |
| 1.6 | Minor | design.md oversold the hom-snspd example | Corrected |
| 1.7 | Wrong | The closed-form heralding truth used the unconditional live fraction: 1.3% low at 1 Mpairs/s, 20 ns, η 0.8 (0/30 inside the interval); the docstring invited correcting with the average dead fraction | Fixed: conditional factor 1 + η_a η_b R min(τ_a, τ_b), within 0.06% at five rates and dead-time pairs; docstring and a note corrected; high-rate SNSPD gate regime added |
| 1.8 | Questionable | The off-peak check cannot see biases of about 2σ of the accidentals interval | Stated as a limit (README) |
| 1.9 | Questionable | Corrected visibility refused when A exceeded the minimum count, censoring 26–35% of true high-visibility datasets; raw V < 0 refused | Fixed: raw value kept with its exact interval, corrected value left out and labelled; V < 0 allowed; LOW_COUNTS keyed on n_min − A |
| 1.10 | Questionable | Gates have little power (true 93% passes 72% of the time) | Open (Q11), stated in README |
| 1.11 | Minor | Jeffreys heralding interval dips to 88% at 100 heralds near η = 0 or 1; docstring said "exact" | Docstring corrected ("posterior mixture, computed exactly"); dip stated |
| 1.12 | Minor | Shifted windows lose \|o\|/T of the record (coverage 37% at 10 Mcps over 3 ms) | Fixed: exposure Σ(1 − \|o\|/T); 93.5% |
| 1.13 | Minor | Accidentals gate truth counted tags, the estimator counts occupied slots | Fixed: slot-occupancy truth |
| 1.14 | Minor | HOM exact interval assumes equal exposure and integer counts | Stated in the docstring |
| 1.15 | Wrong | A copied setup (`backend: sim`) labelled recorded data as simulated, device "sim" | Fixed: refused unless the tags are marked simulated; `backend: file` and `hardware.tagger` documented; tagger serial checked against PTU `HW_SerialNo` |
| 1.16 | Wrong | `calibrate dark` on source-on tags wrote a confident dark rate (25,400 cps against 500) | Fixed: refused when channels share coincidences; warning above 10x the declared dark rate |
| 1.17 | Questionable | Heralding efficiency depends on the window (0.30–0.47 for 1–8 ns) | `calibrate pairs` prints the implied window capture and η / capture; README limit |
| 1.18 | Questionable | The fixed grid's origin (t = 0 of the file) is unrecorded; shifting tags moves S by up to 1.7σ | Documented in README; recording a declared grid phase is open |
| 1.19 | Questionable | PicoHarp 300 router channel 4 refused | Fixed: channels 0–4, as in PicoQuant's demo code |
| 1.20 | Minor | HydraHarp numbering (sync = 0) undocumented; unhelpful error | Documented; error message explains |
| 1.21 | Minor | PicoHarp twin traceback; README snippet `dadl.adapters` failed; version heading; quickstart help | Fixed |

## Cycle 2: Bell-test experimentalist, time-tagger firmware engineer, research software engineer

| # | Severity | Finding | Status |
| --- | --- | --- | --- |
| 2.1 | Wrong | `bell_visibility` reported V = 0.80 past 1/√2 with no label; its planned-window rule was never enforced | Fixed: fair-sampling and single-basis labels always; optional `window` checked against the policy and recorded, else a "window not given" label |
| 2.2 | Wrong | The quickstart's closed-form S ignored `settings_deg` (2.69 printed beside a measured 0.01); the fixed sign vector was not recorded | Fixed: expectation from Σ sign·V cos 2(α − β); warning when the angles cannot exceed 2 with the fixed signs; signs and angles recorded in the result |
| 2.3 | Questionable | A common offset on every delay moves the grid origin unseen; best of ten offsets doubled false S > 2 + 1.96σ | Fixed: CHSH and HOM need a reference channel at delay 0; `grid_origin_ps` recorded; README states the selection effect |
| 2.4 | Minor | Threshold labels ignored S < −2 | Fixed: \|S\| compared |
| 2.5 | Minor | Double clicks count in two combinations | Documented (docstring, README) |
| 2.6 | Minor | README grid wording; example dropped declared delays; non-integer counts accepted; "(default)" policy text | Fixed |
| 2.7 | Questionable | Coarse-tick delay intervals covered 0–12% (80 ps tick, 21 ps jitter), status ok | Fixed: exact Gaussian ⊗ triangle bin probabilities (method v5); 93–98% coverage, bias < 0.15 ps |
| 2.8 | Questionable | Tagger input offsets and trigger levels never recorded, though design.md promised them | Recorded in `conditions` (`tagger_*`) from PTU headers; files that disagree are refused; comparison at use is open |
| 2.9 | Questionable | Unknown special records and a divided sync input were dropped silently by the analysis path | Fixed: both refused |
| 2.10 | Questionable | About 100 bytes per record in memory; metadata checks read whole files 2–5 times | Header-only reads for metadata; memory stated in README; chunked reading open |
| 2.11 | Minor | Doc drift (flat background, PicoHarp `describe` channels) | Fixed |
| 2.12 | Questionable | `calibrate pairs` crashed (capture 0) after writing; with a 15 ns delay error it wrote η = 0, status ok | Fixed: guard; a delay mismatch notes the accidentals and refuses the heralding artifact |
| 2.13 | Questionable | Artifacts with values ≥ 2^53 written but unreadable | Fixed: read back as floats |
| 2.14 | Questionable | Verification checked parsed content, not bytes; duplicate keys accepted | Fixed: bytes must be canonical; duplicate keys refused |
| 2.15 | Questionable | Same data gave a new artifact ID every run (`created` now, unseeded draws) | Fixed: seeds from input hashes; `SOURCE_DATE_EPOCH`; contract in README |
| 2.16 | Minor | Bools accepted as integers | Fixed |
| 2.17 | Minor | Manifests not checked (0 draws, string conditions, valid_until < created) | Fixed in `build_manifest` |
| 2.18 | Minor | CSV out-of-range integer traceback; BOM refused | Fixed |
| 2.19 | Minor | Corrupt `HW_InpChannels` exhausted memory | Fixed: capped at 256 |
| 2.20 | Minor | Library ValueErrors reached CLI users as tracebacks; huge counts leaked a SciPy TypeError | Fixed |
| 2.21 | Minor | Schemas not in the wheel; version still 0.0.3 | Schemas shipped (`dadl/schemas`, guarded copy); version bump at release |

## Cycle 3: Bayesian computation specialist, calibration-lab quality manager, HOM experimentalist

| # | Severity | Finding | Status |
| --- | --- | --- | --- |
| 3.1 | Wrong | With the setup's shared seed, independent channels' draws were correlated at 0.93–1.00 (spread of a difference 180x too small) | Fixed: the input hashes always enter the draw seed; test of independence |
| 3.2 | Wrong | Jitter pinned at a tick/4 bound (0% coverage below 20 ps on an 80 ps tick), Laplace at the bound, status ok | Fixed: bound tick/100, multi-start in σ (a second optimum at σ → 0), a jitter not resolved against the tick is a profile upper limit (unidentified, stand-in draws) with the delay conditional; delay 93–97%, jitter 95–97% for 10–35 ps |
| 3.3 | Questionable | Stand-in draws (uniform below an upper limit) not machine-readable | Fixed: `posterior.stand_in`; limit-only summaries have `value: null`, `bound: "upper"` |
| 3.4 | Minor | Heralding stored no draws of the accidental probability, correlated with η | Fixed: joint `accidental_probability` column |
| 3.5 | Minor | No importance-sampling diagnostic | Fixed: effective sample size in `conditions`, note below 20% |
| 3.6 | Wrong | Unidentified dead time stored as a point value | Fixed (3.3) |
| 3.7 | Wrong | `SOURCE_DATE_EPOCH` silently backdated records; no acquisition date | Fixed: the environment is not read; `--created`; `acquired` from the PTU header or `--acquired` |
| 3.8 | Wrong | The calibrated detector could be left unidentified (only the tagger required) | Fixed: detector serial required for hardware calibrations |
| 3.9 | Questionable | Dark records had no duration, setting, temperature or bias; no instrument identity | Fixed: duration and setting recorded, `--condition key=value`, `instrument` block from PTU headers |
| 3.10 | Questionable | Input hashes not traceable to files; no setup hash | Fixed: file hashes (`raw-file`) and setup hash (`config`) recorded |
| 3.11 | Questionable | Validity and supersession on paper only; invalid `supersedes` accepted | `supersedes` validated and must exist in `--out`; `--valid-until`; docs say recorded, not enforced |
| 3.12 | Questionable | Integrity only; summary not cross-checked against the posterior | Stated: integrity, not authenticity (signing planned) |
| 3.13 | Questionable | `prior_dominated_params` always empty; software versions incomplete | Removed; numpy, scipy and Python versions recorded |
| 3.14 | Minor | Units were prose | Fixed: UCUM units plus a `definition` per parameter |
| 3.15 | Minor | Site "Signed calibration records", supersession wording; design.md example fields not written | Corrected |
| 3.16 | Wrong | Unequal dip and reference exposure biased V (0.960 for 0.92, 50σ) unlabelled | Fixed: exposure ratio from the tag spans; V compares rates, exact interval included |
| 3.17 | Wrong (register) | Q5's "normalise by singles" is wrong for HOM with non-PNR detectors | Corrected in the register and hom.py |
| 3.18 | Questionable | `hom_from_tags` dropped the not-computed and low-count labels | Fixed |
| 3.19 | Questionable | V reads low when the reference delay is not small against the window | Stated; CLI warning above 0.1 W |
| 3.20 | Minor | HOM limits missing from README; `read_tags` example; `hom_from_tags` not exported; subtracted V unqualified | Fixed |

## Cycle 4: QKD security engineer, numerical analyst, teaching-lab instructor

| # | Severity | Finding | Status |
| --- | --- | --- | --- |
| 4.1 | Wrong | A latch at the start or end of a run left no internal gap: heralding 0.457 for 0.653 and HOM +8σ, status ok (1.1 only partly fixed) | Fixed: each channel's edge gaps are tested against its acquisition's bounds; leading and trailing latches caught, no false alarms in 1,200 clean runs |
| 4.2 | Wrong (docs) | The site described a health monitor, blinding detection and signed QKD records as present | Corrected: marked planned; README, SECURITY.md and the site say plainly that DADL has no side-channel detection and nothing it produces belongs in a QKD security argument |
| 4.3 | Questionable | ETSI GS QKD 011 claim without gated detectors or efficiency mismatch | Marked planned; non-support stated |
| 4.4 | Minor | Records with contradicting fields (status vs p, value outside its interval) verified | Fixed: checked on read; "verify" documented as self-consistency, not authenticity |
| 4.5 | Minor | Crafted PTU overflow records wrapped int64 silently | Fixed: refused |
| 4.6 | Minor | A file without a stop reason was treated as complete | An event says completeness is unknown |
| 4.7 | Wrong | The Hessian step on the delay scaled with \|d\| (100 ps at 1 µs): intervals 2–40x too wide or "not identified" at fibre delays | Fixed: absolute step scaled to the jitter; widths identical from 1 ns to 5 µs |
| 4.8 | Questionable | Importance reweighting built a 4000 × n array (2.6 GB at 40,000 gaps) | Fixed: computed in blocks |
| 4.9 | Minor | The jitter upper limit used Δnll = 1.92 (97.5%) | Fixed: 1.35, a one-sided 95% limit |
| 4.10 | Wrong | Dead time raised raw S above 2 with no label (S = 2.14 ± 0.004 from a V = 0.7 source at 1 µs, 3·10⁶ pairs/s) | Labelled when coincidence rate × shortest gap exceeds 1%; README limit; quickstart text |
| 4.11 | Wrong | Tags in ns read as ps gave S = 2.23 ± 0.01 | Fixed: refused when gaps are impossible for the declared dead time, naming the likely unit error |
| 4.12 | Questionable | Corrected HOM values carried no label in the record | Fixed: "accidentals subtracted ... secondary" label |
| 4.13 | Minor | Stale site listing (S = 2.689, old label); teaching claim covered dead time and afterpulsing | Regenerated; claim narrowed |
| 4.14 | Minor | Coincidences per singles and the memory loophole unexplained; narrow-window error message; 16-measurement CHSH unsupported | Explained in README; message names the window; scope stated |
| 4.15 | Minor | design.md: one afterpulse generation; afterpulse closed-form test listed as outstanding | Corrected |

CI: the coverage gate exceeded its 30-minute limit on GitHub runners after cycle 1 added regimes. The
bright-background regime (45% of the run time) now acquires 0.2 s instead of 0.5 s and the limit is
50 minutes; all 34 checks pass locally.
