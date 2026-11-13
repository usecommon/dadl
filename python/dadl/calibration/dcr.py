"""Temperature dependence of the dark count rate, DCR(T).

Method `arrhenius-wls`: thermal generation follows

    DCR(T) = DCR_ref * exp(-(Ea / k) (1/T - 1/T_ref)),

so log DCR is linear in 1/T. Each input is a dark calibration
(`dark-afterpulse-deadtime`) of the same channel at a recorded temperature
(`conditions["temperature_K"]`); its dark-rate estimate and interval give
log DCR and its standard error, and a weighted least-squares line gives
DCR_ref and the activation energy Ea with their joint Gaussian posterior.
Using the dark calibrations, not raw counts, keeps afterpulses and dead
time out of the rates.

T_ref defaults to the inverse of the weighted mean of 1/T (weights 1/se^2),
where the two parameters are uncorrelated. The predictive check is the chi-square of the residuals
with n - 2 degrees of freedom, so at least three temperatures are needed.

Phase 2 limits: generation only; the field-dependent tunnelling term of
DCR(T, V) and bias-voltage scans are not fitted.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from scipy import stats as _st

from dadl.calibration.base import N_DRAWS, Calibration, draws_rng
from dadl.calibration.dark import KIND as DARK_KIND
from dadl.closed_form import K_B_EV
from dadl.stats import LEVEL, Laplace

KIND = "dcr-temperature"
METHOD = "arrhenius-wls"
METHOD_VERSION = "1"
MIN_TEMPERATURES = 3

UNITS = {"dcr_ref_cps": "s-1", "activation_energy_ev": "eV"}
DEFINITIONS = {
    "dcr_ref_cps": "dark count rate at the reference temperature (Arrhenius fit)",
    "activation_energy_ev": "Arrhenius activation energy of the dark count rate",
}


def calibrate_dcr_temperature(
    darks: Sequence[Calibration],
    *,
    reference_k: float | None = None,
    seed: int | None = None,
    n_draws: int = N_DRAWS,
    accept_misfit_inputs: bool = False,
) -> Calibration:
    """Fit DCR_ref (at `reference_k`) and the activation energy from dark calibrations at several temperatures.

    Misfit dark calibrations are refused unless `accept_misfit_inputs` (applications review, cycle
    6: two misfit cold runs pulled Ea to 0.280 [0.278, 0.282] for a true 0.30).
    """
    darks = list(darks)
    bad = [f"{d.conditions.get('temperature_K', '?')} K" for d in darks if d.status == "misfit"]
    if bad and not accept_misfit_inputs:
        raise ValueError(
            f"the dark calibrations at {', '.join(bad)} are misfits; their dark rates are not reliable. "
            "Pass accept_misfit_inputs=True to fit them anyway."
        )
    if len(darks) < MIN_TEMPERATURES:
        raise ValueError(f"need dark calibrations at {MIN_TEMPERATURES} or more temperatures, got {len(darks)}")
    for d in darks:
        if d.kind != DARK_KIND:
            raise ValueError(f"expected {DARK_KIND} calibrations, got {d.kind!r}")
        if "temperature_K" not in d.conditions:
            raise ValueError("every dark calibration needs conditions['temperature_K']")
    channels = {d.channels for d in darks}
    if len(channels) != 1:
        raise ValueError("all dark calibrations must be of the same channel")
    backends = {d.backend for d in darks}
    if len(backends) != 1:
        raise ValueError("cannot mix simulated and hardware calibrations")

    bias_keys = sorted({k for d in darks for k in d.conditions if k.startswith("bias")})
    current = [k for k in bias_keys if not k.endswith(("_V", "_mV", "_kV"))]
    if current:
        # SNSPD specialist review: SNSPDs are current-biased and their dark counts are not
        # thermally activated; an Arrhenius fit gave Ea = 0.0007 eV with a tight interval.
        raise ValueError(
            f"the dark calibrations record a bias current ({', '.join(current)}): DCR(T) is an Arrhenius model "
            "for voltage-biased SPADs, not for current-biased detectors such as SNSPDs"
        )
    for key in bias_keys:
        values = {d.conditions.get(key) for d in darks}
        if len(values) > 1:
            raise ValueError(
                f"the dark calibrations were taken at different {key} {sorted(str(v) for v in values)}: at fixed "
                "temperature the dark rate depends on bias, so an Arrhenius fit across them is not meaningful"
            )
    temps = np.array([float(d.conditions["temperature_K"]) for d in darks])
    if np.any(temps <= 0):
        raise ValueError("temperature_K must be positive (kelvin)")
    if np.unique(temps).size < MIN_TEMPERATURES:
        raise ValueError(f"need {MIN_TEMPERATURES} or more distinct temperatures")
    z = float(_st.norm.ppf(0.5 + LEVEL / 2))
    y, se = [], []
    for d in darks:
        iv = d.estimates["dark_rate_cps"]
        if not (iv.value > 0 and iv.low > 0):
            raise ValueError("every dark rate must be positive to take its logarithm")
        y.append(math.log(iv.value))
        se.append((math.log(iv.high) - math.log(iv.low)) / (2 * z))
    y_arr, w = np.array(y), 1.0 / np.square(se)
    t_ref = float(np.sum(w) / np.sum(w / temps)) if reference_k is None else float(reference_k)
    x = (1.0 / temps - 1.0 / t_ref) / K_B_EV  # log DCR = log DCR_ref - Ea * x
    design = np.column_stack([np.ones_like(x), -x])
    info = design.T @ (w[:, None] * design)
    beta = np.linalg.solve(info, design.T @ (w * y_arr))
    laplace = Laplace.from_hessian(beta, info)
    resid = y_arr - design @ beta
    chi2 = float(np.sum(w * resid**2))
    dof = len(darks) - 2

    notes = [
        f"{len(darks)} dark calibrations from {temps.min():g} K to {temps.max():g} K; chi-square {chi2:.3g} on {dof} dof",
        f"valid from {temps.min():g} K to {temps.max():g} K only: outside the scan other mechanisms (trap-assisted "
        "tunnelling) can dominate and an Arrhenius extrapolation can be several times off",
        "an effective activation energy at constant bias, for SPADs; not meaningful for SNSPDs",
    ]
    if bad:
        notes.append(f"misfit dark calibrations were accepted as inputs ({', '.join(bad)})")
    few = [d.conditions.get("gaps") for d in darks if d.conditions.get("gaps") is not None and d.conditions["gaps"] < 100]
    if few:
        notes.append(
            f"{len(few)} input(s) had fewer than 100 dark gaps: near the 20-gap minimum, runs that pass are biased "
            "upward (and re-acquiring until one passes biases more); prefer longer cold acquisitions"
        )
    rng = draws_rng(seed, *(sha for d in darks for sha, _ in d.inputs))
    draws = laplace.draws(rng, n_draws)
    inputs = tuple(dict.fromkeys(i for d in darks for i in d.inputs))
    return Calibration(
        kind=KIND,
        method=METHOD,
        method_version=METHOD_VERSION,
        channels=next(iter(channels)),
        estimates={
            "dcr_ref_cps": laplace.interval(0, math.exp),
            "activation_energy_ev": laplace.interval(1, lambda v: v),
        },
        units=dict(UNITS),
        definitions=dict(DEFINITIONS),
        samples=np.column_stack([np.exp(draws[:, 0]), draws[:, 1]]),
        sample_params=("dcr_ref_cps", "activation_energy_ev"),
        converged=True,
        ppc_pvalue=float(_st.chi2.sf(chi2, dof)),
        inputs=inputs,
        backend=next(iter(backends)),
        conditions={
            "reference_K": t_ref,
            "temperature_min_K": float(temps.min()),
            "temperature_max_K": float(temps.max()),
            **{k: float(darks[0].conditions[k]) for k in bias_keys if k in darks[0].conditions},
        },
        notes=tuple(notes),
    )


__all__ = ["KIND", "calibrate_dcr_temperature"]
