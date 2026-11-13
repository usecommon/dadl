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
