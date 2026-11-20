# dadl-quantum

Detector-Aware Data Layer for quantum optics experiments: an open-source
starter kit for single-photon experiments. Import name and command: `dadl`.

```bash
pip install dadl-quantum
dadl quickstart                 # simulate the bundled chsh-spad setup, print raw S and a corrected S
dadl init my-experiment         # copy its setup.yaml to edit
dadl analyze my-experiment/setup.yaml --tags tags.npz   # raw S from your own tags
```

Pre-alpha (0.x): a time-tag simulator with injected detector imperfections;
fixed-grid coincidence counting; CHSH and Bell/HOM visibility analyses with
guardrails that refuse invalid corrections; calibrations of dark rate,
afterpulsing, dead time, DCR(T), delay and jitter, accidentals and heralding
efficiency from your own tags, written as content-addressed artifacts
(`dadl calibrate dark`, `dadl calibrate pairs`); PicoQuant PTU (T2 and T3) input; and
an optional compiled core (`pip install "dadl-quantum[core]"`). Raw values carry
statistical uncertainties; corrected values include calibration uncertainty.
Read the known limits in the repository README before relying on a number.

T3 files come from a sync (usually pulsed) source. DADL models CW light: in a
T3 file the sync is used as a clock only. Corrected values are refused on T3
input unless the setup declares `source.type: cw`. Pulsed sources are not
supported on any input (T2, T3 or plain tags): pulsed accidentals sit in side
peaks at multiples of the repetition period. On T2 and plain tags DADL detects
them through the off-peak windows (not flat) and refuses the
accidental-subtracted values when pulsed accidentals are a large share of the
off-peak counts (about half or more in a typical CHSH run; smaller shares pass
unlabelled, with V and |E| reading low; see Known limits); a period that is a
whole multiple of the window grid is not detected, but then the off-peak
estimate is also right.

- Overview: https://usecommon.github.io/dadl/
- Source, design doc and "using your own data": https://github.com/usecommon/dadl
