# dadl-quantum

Detector-Aware Data Layer for quantum optics experiments: an open-source
starter kit for single-photon experiments. Import name and command: `dadl`.

```bash
pip install dadl-quantum
dadl quickstart                 # simulate the bundled chsh-spad setup, print raw S
dadl init my-experiment         # copy its setup.yaml to edit
dadl analyze my-experiment/setup.yaml --tags tags.npz   # raw S from your own tags
```

Pre-alpha (0.x): a time-tag simulator with injected detector imperfections;
fixed-grid coincidence counting; CHSH and Bell/HOM visibility analyses with
guardrails that refuse invalid corrections; calibrations of dark rate,
afterpulsing, dead time, DCR(T), delay and jitter, accidentals and heralding
efficiency from your own tags, written as content-addressed artifacts
(`dadl calibrate dark`, `dadl calibrate pairs`); PicoQuant PTU (T2) input; and
an optional compiled core (`pip install "dadl-quantum[core]"`). Uncertainties
are statistical only; read the known limits in the repository README before
relying on a number.

- Overview: https://usecommon.github.io/dadl/
- Source, design doc and "using your own data": https://github.com/usecommon/dadl
