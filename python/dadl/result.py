"""The labelled result type.

A DADL result always carries its raw value, an optional corrected value, the
labels that qualify it, and the provenance needed to judge it: the policy
fingerprint, any policy overrides with their reasons, and the coincidence
window.

It is not implicitly convertible to a float, so a corrected number is not
quoted without its raw value and labels by accident (see docs/design.md,
"Labels travel with numbers"). This is a guard against accidents, not against
forgery: code running in the same process can always build or edit a result.
"""

from __future__ import annotations

from dataclasses import dataclass, field

#: What `Estimate.stderr` means in phase 1. Calibration uncertainty is not yet propagated.
UNCERTAINTY_KIND = "statistical only (Gaussian, 1 sigma)"


@dataclass(frozen=True)
class Estimate:
    """A value with its one-sigma standard error (statistical only in phase 1)."""

    value: float
    stderr: float

    def __str__(self) -> str:
        return f"{self.value:.4f} ± {self.stderr:.4f}"


@dataclass(frozen=True)
class LabelledValue:
    """Raw value first, optional corrected value, labels and provenance."""

    analysis: str
    raw: Estimate
    corrected: Estimate | None = None
    labels: tuple[str, ...] = field(default_factory=tuple)
    artifact_ids: tuple[str, ...] = field(default_factory=tuple)
    #: SHA-256 of the effective guardrail policy (rules and thresholds).
    policy_fingerprint: str | None = None
    #: Policy overrides in force, as "analysis: correction (reason)".
    overrides: tuple[str, ...] = field(default_factory=tuple)
    #: Coincidence window used: width, delays and declared source.
    window: dict[str, object] | None = None

    def __float__(self) -> float:
        raise TypeError(
            "A DADL result is not a bare number: use .raw.value or .corrected.value "
            "explicitly, and keep .labels with it."
        )

    def __str__(self) -> str:
        parts = [f"{self.analysis}: raw {self.raw}"]
        if self.corrected is not None:
            parts.append(f"corrected {self.corrected}")
        if self.labels:
            parts.append("[" + "; ".join(self.labels) + "]")
        return "  ".join(parts)

    def as_record(self) -> dict[str, object]:
        """Flat record for CSV or DataFrame export. Labels are a mandatory column."""
        window = self.window or {}
        return {
            "analysis": self.analysis,
            "raw_value": self.raw.value,
            "raw_stderr": self.raw.stderr,
            "corrected_value": None if self.corrected is None else self.corrected.value,
            "corrected_stderr": None if self.corrected is None else self.corrected.stderr,
            "uncertainty": UNCERTAINTY_KIND,
            "labels": "; ".join(self.labels),
            "artifact_ids": "; ".join(self.artifact_ids),
            "policy_fingerprint": self.policy_fingerprint,
            "overrides": "; ".join(self.overrides),
            "window_width_ps": window.get("width_ps"),
            "window_delays_ps": window.get("delays_ps"),
            "window_source": window.get("source"),
        }
