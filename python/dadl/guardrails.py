"""Guardrail policy.

Guardrails are refusals in code, not warnings. Each analysis type declares
which corrections may appear (always as a labelled secondary value, never as
the headline) and where its coincidence window may come from.

Users may loosen a rule, but never silently: an override needs a written
reason, and any policy whose content differs from the shipped default is
labelled "policy modified" on every result. In later phases an override also
needs a human approval token, and agents cannot write policy files.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field

from dadl.errors import GuardrailRefusal

POLICY_MODIFIED = "policy modified"
THRESHOLD_AFTER_CORRECTION = "threshold crossed only after correction"
UNVALIDATED_ANALYSIS = "unvalidated analysis"

#: Thresholds that carry a physical claim, per analysis type.
THRESHOLDS: dict[str, float] = {
    "chsh": 2.0,
    "bell_visibility": 1.0 / math.sqrt(2.0),
    "hom_visibility": 0.5,
}


@dataclass(frozen=True)
class AnalysisRule:
    allowed_corrections: frozenset[str]
    window_sources: frozenset[str]

    def to_json(self) -> dict[str, list[str]]:
        return {
            "allowed_corrections": sorted(self.allowed_corrections),
            "window_sources": sorted(self.window_sources),
        }


#: Applied to any analysis type the policy does not know, e.g. plugin analyses.
STRICTEST = AnalysisRule(allowed_corrections=frozenset(), window_sources=frozenset({"plan"}))

DEFAULT_RULES: dict[str, AnalysisRule] = {
    # Bell-type: raw counts only, fixed windows declared in the analysis plan.
    "chsh": STRICTEST,
    "bell_visibility": STRICTEST,
    # HOM and tomography: corrections allowed as labelled secondary values.
    "hom_visibility": AnalysisRule(
        allowed_corrections=frozenset({"accidental_subtraction", "afterpulse"}),
        window_sources=frozenset({"plan"}),
    ),
    "tomography": AnalysisRule(
        allowed_corrections=frozenset({"accidental_subtraction"}),
        window_sources=frozenset({"plan"}),
    ),
}


@dataclass(frozen=True)
class Override:
    analysis: str
    correction: str
    reason: str


@dataclass(frozen=True)
class Policy:
    rules: dict[str, AnalysisRule] = field(default_factory=lambda: dict(DEFAULT_RULES))
    overrides: tuple[Override, ...] = ()

    @classmethod
    def default(cls) -> Policy:
        return cls()

    # -- content hashing ---------------------------------------------------

    def fingerprint(self) -> str:
        """SHA-256 of the canonical rule content (overrides applied)."""
        content = {name: rule.to_json() for name, rule in sorted(self._effective_rules().items())}
        blob = json.dumps(content, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(blob).hexdigest()

    @property
    def is_modified(self) -> bool:
        return self.fingerprint() != _DEFAULT_FINGERPRINT

    # -- loosening ---------------------------------------------------------

    def allow(self, analysis: str, correction: str, reason: str) -> Policy:
        """Return a new policy that allows `correction` for `analysis`.

        A written reason is mandatory. Results under the new policy carry the
        "policy modified" label.
        """
        if not reason or not reason.strip():
            raise ValueError("A policy override needs a written reason.")
        return Policy(rules=dict(self.rules), overrides=self.overrides + (Override(analysis, correction, reason.strip()),))

    # -- checks ------------------------------------------------------------

    def rule_for(self, analysis: str) -> AnalysisRule:
        return self._effective_rules().get(analysis, STRICTEST)

    def check_correction(self, analysis: str, correction: str) -> None:
        if correction not in self.rule_for(analysis).allowed_corrections:
            raise GuardrailRefusal(
                f"Correction '{correction}' is not allowed for '{analysis}'. "
                "Bell-type headline numbers use raw counts only; see docs/design.md, "
                "'Scientific guardrails'."
            )

    def check_window(self, analysis: str, window_source: str) -> None:
        if window_source not in self.rule_for(analysis).window_sources:
            raise GuardrailRefusal(
                f"A coincidence window from '{window_source}' is not allowed for '{analysis}'. "
                "Windows must be fixed and declared in the analysis plan before acquisition; "
                "windows chosen from the data open the coincidence-time loophole."
            )

    def labels_for(self, analysis: str) -> tuple[str, ...]:
        labels: list[str] = []
        if self.is_modified:
            labels.append(POLICY_MODIFIED)
        if analysis not in DEFAULT_RULES:
            labels.append(UNVALIDATED_ANALYSIS)
        return tuple(labels)

    def _effective_rules(self) -> dict[str, AnalysisRule]:
        rules = dict(self.rules)
        for o in self.overrides:
            base = rules.get(o.analysis, STRICTEST)
            rules[o.analysis] = AnalysisRule(
                allowed_corrections=base.allowed_corrections | {o.correction},
                window_sources=base.window_sources,
            )
        return rules


_DEFAULT_FINGERPRINT = Policy().fingerprint()


def threshold_labels(analysis: str, raw: float, corrected: float | None) -> tuple[str, ...]:
    """Label a result whose corrected value crosses a physical threshold that raw does not."""
    threshold = THRESHOLDS.get(analysis)
    if threshold is None or corrected is None:
        return ()
    if raw < threshold <= corrected:
        return (THRESHOLD_AFTER_CORRECTION,)
    return ()
