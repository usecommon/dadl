"""Two-photon visibility: Bell-type fringe visibility and HOM dip visibility.

The analysis type is fixed by the function you call, not by an argument, so a
Bell dataset cannot be filed under the looser HOM rule.

- `bell_visibility`: fringe visibility V = (N_max - N_min) / (N_max + N_min).
  Accidental subtraction is refused under the default policy: it can move a
  non-violating dataset past the 1/sqrt(2) threshold. (V > 1/sqrt(2) implies
  S > 2 only for a Werner state with that visibility in both the H/V and the
  diagonal basis.)
- `hom_visibility`: dip visibility V = (C_max - C_min) / C_max, where C_max is
  the coincidence rate far from the dip. Accidental subtraction is allowed as
  a labelled secondary value.

Uncertainties are first-order (delta method) with Poisson counts. The
accidentals estimate has its own variance: by default Poisson (sigma_A^2 = A),
or `accidentals_stderr` if you know it better (e.g. from many off-peak bins).
"""

from __future__ import annotations

import math

from dadl import _checks
from dadl.guardrails import Policy, threshold_labels
from dadl.result import Estimate, LabelledValue


def _counts(n_max: float, n_min: float) -> tuple[float, float]:
    n_max = _checks.finite(n_max, "n_max", minimum=0.0)
    n_min = _checks.finite(n_min, "n_min", minimum=0.0)
    if n_max + n_min <= 0:
        raise ValueError("counts must not both be zero")
    if n_min > n_max:
        raise ValueError("n_min must not exceed n_max")
    return n_max, n_min


def _accidentals(accidentals: float, accidentals_stderr: float | None, n_min: float) -> tuple[float, float]:
    a = _checks.finite(accidentals, "accidentals", minimum=0.0)
    if a > n_min:
        raise ValueError(
            f"accidentals ({a:g}) exceed the minimum count ({n_min:g}); the corrected visibility would exceed 1"
        )
    if accidentals_stderr is None:
        var_a = a  # Poisson
    else:
        var_a = _checks.finite(accidentals_stderr, "accidentals_stderr", minimum=0.0) ** 2
    return a, var_a


def fringe_visibility_corrected(m: float, n: float, a: float, var_a: float) -> Estimate:
    """(M - m) / (M + m - 2A) with first-order error, M and m Poisson.

    With D = M + m - 2A: dV/dM = 2(m - A)/D^2, dV/dm = -2(M - A)/D^2, dV/dA = 2(M - m)/D^2.
    """
    big, small = m, n
    d = big + small - 2.0 * a
    if d <= 0:
        raise ValueError("accidentals leave no signal")
    v = (big - small) / d
    dv_dbig = 2.0 * (small - a) / d**2
    dv_dsmall = -2.0 * (big - a) / d**2
    dv_da = 2.0 * (big - small) / d**2
    var = dv_dbig**2 * big + dv_dsmall**2 * small + dv_da**2 * var_a
    return Estimate(v, math.sqrt(var))


def hom_visibility_corrected(c_max: float, c_min: float, a: float, var_a: float) -> Estimate:
    """(C_max - C_min) / (C_max - A) with first-order error.

    dV/dC_max = (C_min - A)/(C_max - A)^2, dV/dC_min = -1/(C_max - A), dV/dA = (C_max - C_min)/(C_max - A)^2.
    """
    d = c_max - a
    if d <= 0:
        raise ValueError("accidentals leave no signal")
    v = (c_max - c_min) / d
    var = ((c_min - a) / d**2) ** 2 * c_max + (1.0 / d) ** 2 * c_min + ((c_max - c_min) / d**2) ** 2 * var_a
    return Estimate(v, math.sqrt(var))


def _finish(
    analysis: str,
    raw: Estimate,
    corrected: Estimate | None,
    policy: Policy,
) -> LabelledValue:
    labels = policy.labels_for(analysis) + threshold_labels(analysis, raw, corrected)
    return LabelledValue(analysis=analysis, raw=raw, corrected=corrected, labels=labels, **policy.provenance())


def bell_visibility(
    n_max: float,
    n_min: float,
    *,
    accidentals: float | None = None,
    accidentals_stderr: float | None = None,
    policy: Policy | None = None,
) -> LabelledValue:
    """Raw fringe visibility for a Bell-type claim.

    Passing `accidentals` raises GuardrailRefusal under the default policy. Only
    a policy loosened with `Policy.allow(..., reason=...)` returns a corrected
    value, and that result is labelled "policy modified".
    """
    policy = policy or Policy.default()
    n_max, n_min = _counts(n_max, n_min)
    total = n_max + n_min
    raw = Estimate((n_max - n_min) / total, math.sqrt(4.0 * n_max * n_min / total**3))
    corrected = None
    if accidentals is not None:
        policy.check_correction("bell_visibility", "accidental_subtraction")
        a, var_a = _accidentals(accidentals, accidentals_stderr, n_min)
        corrected = fringe_visibility_corrected(n_max, n_min, a, var_a)
    return _finish("bell_visibility", raw, corrected, policy)


def hom_visibility(
    c_max: float,
    c_min: float,
    *,
    accidentals: float | None = None,
    accidentals_stderr: float | None = None,
    policy: Policy | None = None,
) -> LabelledValue:
    """HOM dip visibility V = (C_max - C_min) / C_max, optionally with accidentals subtracted."""
    policy = policy or Policy.default()
    c_max, c_min = _counts(c_max, c_min)
    if c_max <= 0:
        raise ValueError("c_max must be positive")
    # V = 1 - m/M: var = (m/M^2)^2 M + (1/M)^2 m
    raw = Estimate((c_max - c_min) / c_max, math.sqrt(c_min**2 / c_max**3 + c_min / c_max**2))
    corrected = None
    if accidentals is not None:
        policy.check_correction("hom_visibility", "accidental_subtraction")
        a, var_a = _accidentals(accidentals, accidentals_stderr, c_min)
        corrected = hom_visibility_corrected(c_max, c_min, a, var_a)
    return _finish("hom_visibility", raw, corrected, policy)
