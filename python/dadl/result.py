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
    """A value with its one-sigma standard error (statistical only in phase 1).

    `low` and `high`, when present, are an exact 95% interval (e.g. Clopper-Pearson for
    visibilities); threshold labels use them in preference to value +- 1.96 stderr, which is
    unreliable at low counts.
    """

    value: float
    stderr: float
    low: float | None = None
    high: float | None = None

    def __str__(self) -> str:
        text = f"{self.value:.4f} ± {self.stderr:.4f}"
        if self.low is not None and self.high is not None:
            text += f" (exact 95%: [{self.low:.4f}, {self.high:.4f}])"
        return text

    def __format__(self, spec: str) -> str:
        """Format value and stderr with `spec` (e.g. f"{r.raw:.3f}"); an empty spec is str()."""
        if not spec:
            return str(self)
        text = f"{format(self.value, spec)} ± {format(self.stderr, spec)}"
        if self.low is not None and self.high is not None:
            text += f" (exact 95%: [{format(self.low, spec)}, {format(self.high, spec)}])"
        return text


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
    #: The counts behind the value (e.g. singles, coincidences and their ratio per setting), so a
    #: reader can judge what post-selection kept.
    counts: dict[str, object] | None = None
    #: SHA-256 of the inputs (setup file, tag files), so a result can be tied to what it read.
    inputs: dict[str, str] | None = None

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
        """Flat record for CSV or DataFrame export. Labels are a mandatory column.

        List and dict cells are JSON (labels contain "; " themselves, and Python repr did not
        survive a CSV round trip: archivist review, cycle 8).
        """
        import json

        from dadl import __version__
        from dadl.guardrails import SIMULATED_INPUT

        def js(v: object) -> str | None:
            return None if v is None else json.dumps(v, sort_keys=True, default=str)

        window = self.window or {}
        return {
            "dadl_version": __version__,
            "analysis": self.analysis,
            "simulated": SIMULATED_INPUT in self.labels,
            "raw_value": self.raw.value,
            "raw_stderr": self.raw.stderr,
            "raw_low": self.raw.low,
            "raw_high": self.raw.high,
            "corrected_value": None if self.corrected is None else self.corrected.value,
            "corrected_stderr": None if self.corrected is None else self.corrected.stderr,
            "corrected_low": None if self.corrected is None else self.corrected.low,
            "corrected_high": None if self.corrected is None else self.corrected.high,
            "uncertainty": UNCERTAINTY_KIND,
            "labels": js(list(self.labels)),
            "artifact_ids": js(list(self.artifact_ids)),
            "policy_fingerprint": self.policy_fingerprint,
            "overrides": js(list(self.overrides)),
            "window_width_ps": window.get("width_ps"),
            "window_delays_ps": js(dict(window["delays_ps"])) if window.get("delays_ps") is not None else None,
            "window_source": window.get("source"),
            "window_grid_origin_ps": window.get("grid_origin_ps"),
            "counts": js(self.counts),
            "inputs": js(self.inputs),
        }
