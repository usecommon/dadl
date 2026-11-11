# dadl (Python package)

Detector-Aware Data Layer for quantum optics experiments.

Phase 1 is pure NumPy: a time-tag simulator with injected detector
imperfections, fixed-grid coincidence counting, CHSH and two-photon
visibility analyses, the guardrail policy, and the labelled result type.

```bash
pip install -e ".[test]"
dadl quickstart chsh-spad
pytest
```

See the repository README and `docs/design.md` for the full design.
