# dadl-quantum

Detector-Aware Data Layer for quantum optics experiments: an open-source
starter kit for single-photon experiments. Import name and command: `dadl`.

```bash
pip install dadl-quantum
dadl quickstart                 # simulate the bundled chsh-spad setup, print raw S
dadl init my-experiment         # copy its setup.yaml to edit
dadl analyze my-experiment/setup.yaml --tags tags.npz   # raw S from your own tags
```

Phase 1 (pre-alpha, 0.0.x): a time-tag simulator with injected detector
imperfections, fixed-grid coincidence counting, CHSH and Bell/HOM visibility
analyses, guardrails that refuse invalid corrections, and a labelled result
type that records the policy and window used. Uncertainties are statistical
only; calibration fits arrive in phase 2.

- Overview: https://usecommon.github.io/dadl/
- Source, design doc and "using your own data": https://github.com/usecommon/dadl
