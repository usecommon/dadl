"""The labelled result type.

A DADL result always carries its raw value, an optional corrected value and
the labels that qualify it. It deliberately refuses to become a bare float, so
a corrected number cannot be quoted without its raw value and labels by
accident (see docs/design.md, "Labels travel with numbers").
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Estimate:
    """A value with its one-sigma standard error."""

    value: float
    stderr: float

    def __str__(self) -> str:
        return f"{self.value:.4f} ± {self.stderr:.4f}"


@dataclass(frozen=True)
class LabelledValue:
    """Raw value first, optional corrected value, and labels."""

    analysis: str
    raw: Estimate
    corrected: Estimate | None = None
    labels: tuple[str, ...] = field(default_factory=tuple)
    artifact_ids: tuple[str, ...] = field(default_factory=tuple)

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
        return {
            "analysis": self.analysis,
            "raw_value": self.raw.value,
            "raw_stderr": self.raw.stderr,
            "corrected_value": None if self.corrected is None else self.corrected.value,
            "corrected_stderr": None if self.corrected is None else self.corrected.stderr,
            "labels": "; ".join(self.labels),
            "artifact_ids": "; ".join(self.artifact_ids),
        }
