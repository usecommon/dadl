"""Guardrail policy.

Guardrails are refusals in code, not warnings. Each analysis type declares
which corrections may appear (always as a labelled secondary value, never as
the headline) and where its coincidence window may come from.

Users may loosen a rule, but never silently: loosening goes through
`Policy.allow`, which needs a written reason; every result records the
policy fingerprint and the overrides with their reasons, and any policy
that differs from the shipped default is labelled "policy modified".

What this does not do in phase 1: it cannot verify that a window declared as
coming from the plan really was fixed before acquisition, and it trusts the
counts and tags it is given. Signed analysis plans (phase 3) close the first
gap; dataset IDs computed from raw tags in the core (phase 2) the second.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

from dadl.errors import GuardrailRefusal
from dadl.result import Estimate

POLICY_MODIFIED = "policy modified"
THRESHOLD_AFTER_CORRECTION = "threshold crossed only after correction"
NEAR_THRESHOLD = "within 1.96 sigma of the threshold"
GAUSSIAN_INVALID = "Gaussian error invalid at boundary (stderr is zero)"
NEAR_THRESHOLD_EXACT = "the exact 95% interval contains the threshold"
LOW_COUNTS = "low counts: the Gaussian error is unreliable; use the exact interval"
UNVALIDATED_ANALYSIS = "unvalidated analysis"
SINGLE_BASIS = (
    "single-basis fringe: V > 1/sqrt(2) implies S > 2 only for a Werner state with this visibility in two "
    "mutually unbiased bases"
)
WINDOW_UNDECLARED = "coincidence window not given: the counts are trusted as given"
ACCIDENTALS_SUBTRACTED = (
    "corrected value: accidentals subtracted (off-peak estimate); secondary, not an indistinguishability or "
    "Bell claim"
)
DEAD_TIME_SELECTION = (
    "dead-time selection: coincidence rate x shortest gap exceeds 1%; a recorded click makes the partner arm's "
    "dead state depend on the outcome, which can raise |E| and S"
)
SIMULATED_INPUT = "simulated data: the input tags were written by DADL's simulator, not recorded"
CORRECTION_NOT_COMPUTED = (
    "corrected value not computed: the accidentals exceed the minimum count; report the raw value and its "
    "exact interval"
)
INTERRUPTED = (
    "interrupted: a channel has a gap far longer than its count rate allows (latch, shutdown or lost data); "
    "the result is not reliable"
)

#: Corrections the library knows how to apply. `Policy.allow` refuses others.
KNOWN_CORRECTIONS = frozenset({"accidental_subtraction", "afterpulse", "dead_time", "efficiency"})

#: Thresholds that carry a physical claim, per analysis type.
#: bell_visibility: V > 1/sqrt(2) implies S > 2 only for a Werner state with
#: that visibility in both the H/V and the diagonal basis.
#: hom_visibility: dip visibility (C_max - C_min) / C_max above 0.5 is non-classical.
THRESHOLDS: Mapping[str, float] = MappingProxyType(
    {
        "chsh": 2.0,
        "bell_visibility": 1.0 / math.sqrt(2.0),
        "hom_visibility": 0.5,
    }
)


@dataclass(frozen=True)
class AnalysisRule:
    allowed_corrections: frozenset[str]
    window_sources: frozenset[str]

    def to_json(self) -> dict[str, list[str]]:
        return {
            "allowed_corrections": sorted(self.allowed_corrections),
            "window_sources": sorted(self.window_sources),
        }

    def is_at_least_as_strict_as(self, other: AnalysisRule) -> bool:
        return self.allowed_corrections <= other.allowed_corrections and self.window_sources <= other.window_sources


#: Applied to any analysis type the policy does not know, e.g. plugin analyses.
STRICTEST = AnalysisRule(allowed_corrections=frozenset(), window_sources=frozenset({"plan"}))

DEFAULT_RULES: Mapping[str, AnalysisRule] = MappingProxyType(
    {
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
)

#: SHA-256 of the shipped default policy. Pinned here (not computed at import)
#: so that editing the rule tables in place cannot redefine "default";
#: tests/test_guardrails.py checks that it matches.
DEFAULT_POLICY_FINGERPRINT = "0b6cd03bc7f74cc7ae8d63c714764d4d12c3b11c19a7e1a8be3bb066e705b4ec"


@dataclass(frozen=True)
class Override:
    analysis: str
    correction: str
    reason: str

    def describe(self) -> str:
        return f"{self.analysis}: {self.correction} ({self.reason})"


@dataclass(frozen=True, eq=False)
class Policy:
    rules: Mapping[str, AnalysisRule] = field(default_factory=lambda: DEFAULT_RULES)
    overrides: tuple[Override, ...] = ()

    def __post_init__(self) -> None:
        # Rules are rebuilt from their raw fields: a subclass overriding to_json or the strictness
        # test made a loosened policy fingerprint as the default (red-team review, cycle 5).
        clean: dict[str, AnalysisRule] = {}
        for name, rule in dict(self.rules).items():
            if type(rule) is not AnalysisRule:
                raise TypeError(f"policy rule for {name!r} must be an AnalysisRule, not {type(rule).__name__}")
            clean[str(name)] = AnalysisRule(
                frozenset(str(c) for c in rule.allowed_corrections), frozenset(str(w) for w in rule.window_sources)
            )
        rules = MappingProxyType(clean)
        for name, rule in rules.items():
            baseline = DEFAULT_RULES.get(name, STRICTEST)
            if not (
                rule.allowed_corrections <= baseline.allowed_corrections
                and rule.window_sources <= baseline.window_sources
            ):
                raise ValueError(
                    f"Policy rules for '{name}' are looser than the default. "
                    "Loosen a rule with Policy.allow(analysis, correction, reason), so the reason is recorded."
                )
        object.__setattr__(self, "rules", rules)

    @classmethod
    def default(cls) -> Policy:
        return cls()

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Policy) and (self.fingerprint(), self.overrides) == (
            other.fingerprint(),
            other.overrides,
        )

    def __hash__(self) -> int:
        return hash((self.fingerprint(), self.overrides))

    # -- content hashing ---------------------------------------------------

    def fingerprint(self) -> str:
        """SHA-256 of the canonical effective rules (overrides applied) and the thresholds."""
        content = {
            "rules": {name: rule.to_json() for name, rule in sorted(self._effective_rules().items())},
            "thresholds": {name: repr(value) for name, value in sorted(THRESHOLDS.items())},
        }
        blob = json.dumps(content, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(blob).hexdigest()

    @property
    def is_modified(self) -> bool:
        return self.fingerprint() != DEFAULT_POLICY_FINGERPRINT

    # -- loosening ---------------------------------------------------------

    def allow(self, analysis: str, correction: str, reason: str) -> Policy:
        """Return a new policy that allows `correction` for `analysis`.

        A written reason is mandatory. Results under the new policy carry the
        "policy modified" label and record the override with its reason.
        """
        if correction not in KNOWN_CORRECTIONS:
            raise ValueError(f"unknown correction {correction!r}; known: {sorted(KNOWN_CORRECTIONS)}")
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("A policy override needs a written reason.")
        return Policy(rules=self.rules, overrides=self.overrides + (Override(analysis, correction, reason.strip()),))

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

    def provenance(self) -> dict[str, object]:
        """Fields every result records: fingerprint and overrides with reasons."""
        return {
            "policy_fingerprint": self.fingerprint(),
            "overrides": tuple(o.describe() for o in self.overrides),
        }

    def _effective_rules(self) -> dict[str, AnalysisRule]:
        rules = dict(self.rules)
        for o in self.overrides:
            base = rules.get(o.analysis, STRICTEST)
            rules[o.analysis] = AnalysisRule(
                allowed_corrections=base.allowed_corrections | {o.correction},
                window_sources=base.window_sources,
            )
        return rules


def _mirror(est: Estimate) -> Estimate:
    low = None if est.high is None else -est.high
    high = None if est.low is None else -est.low
    return Estimate(-est.value, est.stderr, low, high)


def threshold_labels(analysis: str, raw: Estimate, corrected: Estimate | None = None) -> tuple[str, ...]:
    """Labels that qualify a value against the analysis's physical threshold.

    - "threshold crossed only after correction": raw is below, corrected at or above.
    - "within 1.96 sigma of the threshold": the 95 % Gaussian interval of the
      headline (raw) value, or of the corrected value, contains the threshold.
    - "Gaussian error invalid at boundary": a stderr of zero (e.g. a correlator
      of exactly +-1 or n_min = 0), where the Gaussian interval means nothing.
    """
    labels: list[str] = []
    if analysis == "chsh" and raw.value < 0:
        # |S| > 2 violates the inequality whatever the sign convention (Bell review, cycle 2:
        # S = -2.01 +- 0.02 carried no label).
        raw = _mirror(raw)
        corrected = None if corrected is None else _mirror(corrected)
    if raw.stderr == 0 or (corrected is not None and corrected.stderr == 0):
        labels.append(GAUSSIAN_INVALID)
    threshold = THRESHOLDS.get(analysis)
    if threshold is None:
        return tuple(labels)
    if corrected is not None and raw.value < threshold <= corrected.value:
        labels.append(THRESHOLD_AFTER_CORRECTION)
    for est in (raw, corrected):
        if est is None:
            continue
        if est.low is not None and est.high is not None:
            if est.low <= threshold <= est.high:
                labels.append(NEAR_THRESHOLD_EXACT)
                break
        elif abs(est.value - threshold) < 1.96 * est.stderr:
            labels.append(NEAR_THRESHOLD)
            break
    return tuple(labels)
