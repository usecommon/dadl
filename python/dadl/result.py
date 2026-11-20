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

#: What `Estimate.stderr` means for raw values (and corrected ones without propagated calibrations).
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
    #: Corrected values from the forward model (phase 3): the statistical and calibration parts
    #: of `stderr`, each with the other part held at its point value. None for raw values.
    stderr_statistical: float | None = None
    stderr_calibration: float | None = None

    @property
    def uncertainty_kind(self) -> str:
        if self.stderr_calibration is None:
            return UNCERTAINTY_KIND
        return "statistical and calibration (joint Monte Carlo, 1 sigma; interval is 95% quantiles)"

    def __str__(self) -> str:
        text = f"{self.value:.4f} ± {self.stderr:.4f}"
        if self.stderr_calibration is not None:
            # A corrected value's "±" is one sd; its interval is 95% (council round 10, Ux11).
            text += " (1 sd)"
        if self.low is not None and self.high is not None:
            kind = "95%" if self.stderr_calibration is not None else "exact 95%"
            text += f" ({kind}: [{self.low:.4f}, {self.high:.4f}])"
        if self.stderr_calibration is not None:
            # Two significant figures: "0.0000" read as exactly known (council round 6, R4). A part
            # below 1e-12 of the value is floating-point noise, printed as 0 (council round 10, Rg5).
            cal = self.stderr_calibration
            shown = "0" if cal == 0 or abs(cal) < 1e-12 * abs(self.value) else f"{cal:.2g}"
            text += (f" [stat {self.stderr_statistical:.2g}, calibration part {shown} (statistical part of the "
                     "artifacts only)]")
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
    """Raw value first, optional corrected value, labels and provenance.

    `labels` qualify the raw value or both values (provenance, simulated input, the analysis's
    caveats); `corrected_labels` qualify the corrected value only (council round 6, R2: a corrected
    value's threshold label was printed beside the raw one). `all_labels` is both.
    """

    analysis: str
    raw: Estimate
    corrected: Estimate | None = None
    labels: tuple[str, ...] = field(default_factory=tuple)
    #: Labels of the corrected value only.
    corrected_labels: tuple[str, ...] = field(default_factory=tuple)
    artifact_ids: tuple[str, ...] = field(default_factory=tuple)
    #: Corrections behind `corrected`, in the order applied (design, "Minimum metadata").
    corrections: tuple[str, ...] = field(default_factory=tuple)
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

    @property
    def all_labels(self) -> tuple[str, ...]:
        """`labels` then `corrected_labels`, without repeats."""
        return tuple(dict.fromkeys(self.labels + self.corrected_labels))

    def __float__(self) -> float:
        raise TypeError(
            "A DADL result is not a bare number: use .raw.value or .corrected.value "
            "explicitly, and keep .labels with it."
        )

    def __str__(self) -> str:
        parts = [f"{self.analysis}: raw {self.raw}"]
        if self.labels:
            parts.append("[" + "; ".join(self.labels) + "]")
        if self.corrected is not None:
            parts.append(f"corrected {self.corrected}")
            if self.corrected_labels:
                parts.append("[" + "; ".join(self.corrected_labels) + "]")
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
            "corrected_uncertainty": None if self.corrected is None else self.corrected.uncertainty_kind,
            "corrected_stderr_statistical": None if self.corrected is None else self.corrected.stderr_statistical,
            "corrected_stderr_calibration": None if self.corrected is None else self.corrected.stderr_calibration,
            "corrections": js(list(self.corrections)),
            "labels": js(list(self.labels)),
            "corrected_labels": js(list(self.corrected_labels)),
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
