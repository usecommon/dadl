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

import dataclasses
import math
from typing import TYPE_CHECKING

from dadl import _checks
from dadl.stats import clopper_pearson
from dadl.guardrails import (
    ACCIDENTALS_SUBTRACTED,
    CORRECTION_NOT_COMPUTED,
    LOW_COUNTS,
    SINGLE_BASIS,
    WINDOW_UNDECLARED,
    Policy,
    threshold_labels,
)
from dadl.result import Estimate, LabelledValue

if TYPE_CHECKING:
    from dadl.analysis.chsh import Window


def _counts(n_max: float, n_min: float) -> tuple[float, float]:
    n_max = _checks.finite(n_max, "n_max", minimum=0.0)
    n_min = _checks.finite(n_min, "n_min", minimum=0.0)
    if n_max != round(n_max) or n_min != round(n_min):
        # Bell review, cycle 2: rates or averaged counts mis-scale the Poisson errors and the
        # exact interval silently rounded them.
        raise ValueError("counts must be whole numbers of coincidences (not rates or averages)")
    if max(n_max, n_min) > 2**53:
        raise ValueError("counts above 2^53 cannot be represented exactly")
    if n_max + n_min <= 0:
        raise ValueError("counts must not both be zero")
    # n_min > n_max is allowed: a fringe near zero visibility or a HOM peak gives V < 0, and
    # refusing it (12% of datasets at V = 0.03) would censor what gets reported.
    return n_max, n_min


def _accidentals(accidentals: float, accidentals_stderr: float | None, n_min: float) -> tuple[float, float]:
    a = _checks.finite(accidentals, "accidentals", minimum=0.0)
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


def hom_visibility_corrected(
    c_max: float, c_min: float, a: float, var_a: float, var_c_min: float | None = None
) -> Estimate:
    """(C_max - C_min) / (C_max - A) with first-order error.

    dV/dC_max = (C_min - A)/(C_max - A)^2, dV/dC_min = -1/(C_max - A), dV/dA = (C_max - C_min)/(C_max - A)^2.
    `var_c_min` defaults to C_min (Poisson); it differs when C_min was rescaled to the reference exposure.
    """
    d = c_max - a
    if d <= 0:
        raise ValueError("accidentals leave no signal")
    v = (c_max - c_min) / d
    var_m = c_min if var_c_min is None else var_c_min
    var = ((c_min - a) / d**2) ** 2 * c_max + (1.0 / d) ** 2 * var_m + ((c_max - c_min) / d**2) ** 2 * var_a
    return Estimate(v, math.sqrt(var))


#: Below this smaller count the delta-method error is unreliable (statistician review: coverage
#: 63-93% near V = 1 or at a handful of counts), so results carry the LOW_COUNTS label.
LOW_COUNT_MIN = 10


def _exact_bell(n_max: float, n_min: float) -> tuple[float, float]:
    """Exact 95% interval for V = 1 - 2 m / n: Clopper-Pearson on m out of n = M + m (counts rounded)."""
    n, m = int(round(n_max + n_min)), int(round(n_min))
    ci = clopper_pearson(m, n)
    return 1.0 - 2.0 * ci.high, 1.0 - 2.0 * ci.low


def _exact_hom(c_max: float, c_min: float, exposure_ratio: float = 1.0) -> tuple[float, float]:
    """Exact 95% interval for V = 1 - rho, rho = rate_min / rate_max.

    With exposures T_dip = r T_ref, given m + M the dip count m is binomial in r rho / (1 + r rho).
    C_max and C_min are integer counts.
    """
    n, m = int(round(c_max + c_min)), int(round(c_min))
    ci = clopper_pearson(m, n)

    def rho(p: float) -> float:
        return math.inf if p >= 1.0 else p / ((1.0 - p) * exposure_ratio)

    return 1.0 - rho(ci.high), 1.0 - rho(ci.low)


def _finish(
    analysis: str,
    raw: Estimate,
    corrected: Estimate | None,
    policy: Policy,
    min_count: float | None = None,
    not_computed: bool = False,
) -> LabelledValue:
    labels = policy.labels_for(analysis) + threshold_labels(analysis, raw, corrected)
    if min_count is not None and min_count < LOW_COUNT_MIN:
        labels += (LOW_COUNTS,)
    if not_computed:
        labels += (CORRECTION_NOT_COMPUTED,)
    if corrected is not None:
        # Teaching-lab review, cycle 4: the record said nothing about what the corrected value was.
        labels += (ACCIDENTALS_SUBTRACTED,)
    return LabelledValue(analysis=analysis, raw=raw, corrected=corrected, labels=labels, **policy.provenance())


def bell_visibility(
    n_max: float,
    n_min: float,
    *,
    accidentals: float | None = None,
    accidentals_stderr: float | None = None,
    policy: Policy | None = None,
    window: "Window | None" = None,
) -> LabelledValue:
    """Raw fringe visibility for a Bell-type claim.

    Always labelled with the fair-sampling and single-basis caveats (Bell review, cycle 2: a
    threshold crossing was reported with no label). Pass the coincidence `window`: it is checked
    against the policy (a planned window only) and recorded; without one, the result says so.

    Passing `accidentals` raises GuardrailRefusal under the default policy. Only
    a policy loosened with `Policy.allow(..., reason=...)` returns a corrected
    value, and that result is labelled "policy modified".
    """
    policy = policy or Policy.default()
    n_max, n_min = _counts(n_max, n_min)
    total = n_max + n_min
    low, high = _exact_bell(n_max, n_min)
    raw = Estimate((n_max - n_min) / total, math.sqrt(4.0 * n_max * n_min / total**3), low, high)
    corrected = None
    not_computed = False
    min_count = n_min
    if accidentals is not None:
        policy.check_correction("bell_visibility", "accidental_subtraction")
        a, var_a = _accidentals(accidentals, accidentals_stderr, n_min)
        # Statistician review (cycle 1): refusing when A > the minimum count dropped 26-35% of
        # true high-visibility datasets, the best ones. The raw value (exact interval) is kept and
        # the corrected value is left out, labelled; the low-count label keys on the signal.
        if a > n_min or a >= n_max:
            not_computed = True
        else:
            corrected = fringe_visibility_corrected(n_max, n_min, a, var_a)
        min_count = n_min - a
    result = _finish("bell_visibility", raw, corrected, policy, min_count, not_computed)
    from dadl.analysis.chsh import FAIR_SAMPLING

    extra = (FAIR_SAMPLING, SINGLE_BASIS)
    if window is None:
        extra += (WINDOW_UNDECLARED,)
    else:
        policy.check_window("bell_visibility", window.source)
    return dataclasses.replace(
        result, labels=extra + result.labels, window=None if window is None else window.as_record()
    )


def hom_visibility(
    c_max: float,
    c_min: float,
    *,
    accidentals: float | None = None,
    accidentals_stderr: float | None = None,
    policy: Policy | None = None,
    exposure_ratio: float = 1.0,
) -> LabelledValue:
    """HOM dip visibility V = 1 - (C_min / T_dip) / (C_max / T_ref), optionally with accidentals subtracted.

    `exposure_ratio` is T_dip / T_ref (1 for equal acquisitions). The HOM review (cycle 3) found a
    dip run half as long as the reference read V = 0.960 for a true 0.92 (50 sigma) with no label.
    `accidentals` are counted in the dip acquisition.
    """
    policy = policy or Policy.default()
    c_max, c_min = _counts(c_max, c_min)
    if c_max <= 0:
        raise ValueError("c_max must be positive")
    r = _checks.finite(exposure_ratio, "exposure_ratio")
    if not r > 0:
        raise ValueError("exposure_ratio must be positive")
    # In reference-exposure units: m' = m / r, var(m') = m / r^2. V = 1 - m'/M.
    m_ref = c_min / r
    low, high = _exact_hom(c_max, c_min, r)
    raw = Estimate(1.0 - m_ref / c_max, math.sqrt(m_ref**2 / c_max**3 + c_min / (r * r * c_max**2)), low, high)
    corrected = None
    not_computed = False
    min_count = c_min
    if accidentals is not None:
        policy.check_correction("hom_visibility", "accidental_subtraction")
        a, var_a = _accidentals(accidentals, accidentals_stderr, c_min)
        # Statistician review (cycle 1): refusing when A > the minimum count dropped 26-35% of
        # true high-visibility datasets, the best ones. The raw value (exact interval) is kept and
        # the corrected value is left out, labelled; the low-count label keys on the signal.
        if a > c_min or a / r >= c_max:
            not_computed = True
        else:
            corrected = hom_visibility_corrected(c_max, m_ref, a / r, var_a / (r * r), var_c_min=c_min / (r * r))
        min_count = c_min - a
    return _finish("hom_visibility", raw, corrected, policy, min_count, not_computed)
