"""Two-photon fringe visibility, with optional accidental subtraction.

For Bell-type visibility claims (analysis "bell_visibility") accidental
subtraction is refused: it can move a non-violating dataset past the 71%
threshold. For HOM it is allowed, as a labelled secondary value.
"""

from __future__ import annotations

import math

from dadl.guardrails import Policy, threshold_labels
from dadl.result import Estimate, LabelledValue


def _visibility(n_max: float, n_min: float) -> tuple[float, float]:
    total = n_max + n_min
    if total <= 0:
        raise ValueError("counts must be positive")
    v = (n_max - n_min) / total
    # Poisson counts: var(V) = 4 n_max n_min / (n_max + n_min)^3
    return v, math.sqrt(4.0 * n_max * n_min / total**3)


def two_photon_visibility(
    n_max: float,
    n_min: float,
    *,
    analysis: str,
    accidentals: float | None = None,
    policy: Policy | None = None,
) -> LabelledValue:
    """Raw visibility, and the accidental-subtracted value only where the policy allows it."""
    policy = policy or Policy.default()
    raw_v, raw_err = _visibility(n_max, n_min)
    corrected: Estimate | None = None
    if accidentals is not None:
        policy.check_correction(analysis, "accidental_subtraction")
        denom = n_max + n_min - 2.0 * accidentals
        if denom <= 0:
            raise ValueError("accidentals exceed the measured counts")
        cor_v = (n_max - n_min) / denom
        corrected = Estimate(cor_v, raw_err * (n_max + n_min) / denom)
    labels = policy.labels_for(analysis) + threshold_labels(
        analysis, raw_v, None if corrected is None else corrected.value
    )
    return LabelledValue(analysis=analysis, raw=Estimate(raw_v, raw_err), corrected=corrected, labels=labels)
