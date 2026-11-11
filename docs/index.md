# DADL

**Detector-Aware Data Layer for quantum optics experiments.**

DADL keeps a versioned model of every single-photon detector channel and
carries its uncertainty into Bell S, visibility and fidelity. It always
reports raw and corrected values side by side, and refuses corrections that
are invalid for the analysis.

- [Design doc](design.md)
- Quickstart: `pip install -e "./python[test]"`, then `dadl quickstart chsh-spad`

How-to guides for each extension point (adapters, detector models, analyses,
guardrail policies, telemetry sources) are an MVP deliverable and will be
added here as each extension point lands.
