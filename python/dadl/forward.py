"""Forward model for two-point visibilities with propagated calibrations (phase 3, M2).

The corrected value is the visibility of the two-photon signal fitted under a model of what the
detectors recorded, never a chain of subtractions. For one acquisition (a run at a fixed setting
or path delay) of duration T with recorded singles n_a, n_b and C coincidences in the declared
fixed-grid window:

    E[C] = T * L * (s + a)

where s is the signal coincidence rate the visibility is defined on, a the accidental rate in the
window, and L the live fraction of the pair of channels: a photon pair is recorded only if both
detectors are live, (1 - n_a tau_a)(1 - n_b tau_b), raised by the conditional factor
1 + r_c min(tau_a, tau_b) for a recorded coincidence rate r_c (as in
`closed_form.heralding_ratio_expected`). The accidental count A of each acquisition is measured in
the same likelihood from its own off-peak bins (`calibration.accidentals.shifted_coincidences`),
not taken from another run, so calibration and analysis data stay disjoint.

With two acquisitions the model is exactly identified and the fit is closed form:

    HOM:     V = 1 - s_dip / s_ref
    fringe:  V = (s_max - s_min) / (s_max + s_min)
    s = (C - A) / (T L)

Uncertainty by joint Monte Carlo: counts from their Jeffreys posteriors (Gamma(c + 1/2)), the
off-peak sums likewise, and the dead times from the dark artifacts of a `CalibrationSet` (joint
within an artifact, independent across artifacts). The interval is the 2.5% and 97.5% quantiles
of the joint draws; the statistical and calibration parts are reported separately as standard
deviations with the other part held at its point value.

Corrections applied, in order: dead time (if dark artifacts are given), then accidentals. Each is
checked against the policy for the analysis; dead time is not in the default HOM rule, so it
needs `Policy.allow(analysis, "dead_time", reason)`.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from dadl.calibration.dark import KIND as DARK_KIND
from dadl.data import ChannelId
from dadl.guardrails import Policy
from dadl.propagate import CalibrationSet
from dadl.result import Estimate

#: Joint Monte Carlo draws per corrected value.
N_MC = 4000
#: Above this dead fraction per channel the first-order live fraction is not trusted.
MAX_DEAD_FRACTION = 0.2

FORWARD_MODEL = "corrected by the forward model (dead time, accidentals), calibration uncertainty included"
DEAD_FRACTION_HIGH = "dead fraction above 20%: the live-fraction model is first order"
DEAD_TIME_BOUND = "dead time known only as an upper limit: corrected at the bound; counts give the value at zero"


@dataclass(frozen=True)
class RunCounts:
    """What the forward model needs from one acquisition."""

    coincidences: int
    #: Sum of coincidences over the off-peak shifted windows, and how many shifts.
    offpeak_sum: int
    offpeak_shifts: int
    duration_s: float
    singles_a: int
    singles_b: int

    def __post_init__(self) -> None:
        if self.duration_s <= 0:
            raise ValueError("duration_s must be positive")
        if self.offpeak_shifts < 1:
            raise ValueError("offpeak_shifts must be at least 1")
        for name in ("coincidences", "offpeak_sum", "singles_a", "singles_b"):
            if getattr(self, name) < 0 or int(getattr(self, name)) != getattr(self, name):
                raise ValueError(f"{name} must be a non-negative whole number")

    @property
    def accidentals(self) -> float:
        return self.offpeak_sum / self.offpeak_shifts


def live_fraction(
    singles_a_cps: np.ndarray | float,
    singles_b_cps: np.ndarray | float,
    dead_a_s: np.ndarray | float,
    dead_b_s: np.ndarray | float,
    coincidence_cps: np.ndarray | float,
) -> np.ndarray:
    """Probability that a pair arriving in the window is recorded on both channels."""
    fa = np.asarray(singles_a_cps) * np.asarray(dead_a_s)
    fb = np.asarray(singles_b_cps) * np.asarray(dead_b_s)
    conditional = 1.0 + np.asarray(coincidence_cps) * np.minimum(dead_a_s, dead_b_s)
    return (1.0 - fa) * (1.0 - fb) * conditional


@dataclass(frozen=True)
class Corrected:
    estimate: Estimate
    corrections: tuple[str, ...]
    labels: tuple[str, ...]
    artifact_ids: tuple[str, ...]


def _signal(run: RunCounts, c, offsum, tau_a, tau_b) -> np.ndarray:
    a = offsum / run.offpeak_shifts
    t = run.duration_s
    live = live_fraction(run.singles_a / t, run.singles_b / t, tau_a, tau_b, run.coincidences / t)
    return (c - a) / (t * live)


def _dead_times(
    calibrations: CalibrationSet | None, channels: tuple[ChannelId, ChannelId], n: int, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Draws and point values of the two dead times (s), or None without dark artifacts."""
    if calibrations is None:
        return None
    darks = [m for m in calibrations.members if m.kind == DARK_KIND]
    if not darks:
        return None
    joint = calibrations.draws(n, seed)
    out = []
    for ch in channels:
        cal = calibrations.get(DARK_KIND, ch)
        if "dead_time_ps" in cal.stand_in or "dead_time_ps" in cal.unidentified:
            # Decision D3: a bound, not a posterior. The caller reports both ends.
            lo, hi = 0.0, float(cal.estimates["dead_time_ps"].high)
            out.append((np.full(n, hi * 1e-12), lo * 1e-12, hi * 1e-12))
        else:
            col = joint.column(DARK_KIND, "dead_time_ps", ch) * 1e-12
            out.append((col, float(cal.estimates["dead_time_ps"].value) * 1e-12, None))
    (da, pa, _), (db, pb, _) = out
    return da, db, np.array(pa), np.array(pb)


def corrected_visibility(
    analysis: str,
    high: RunCounts,
    low: RunCounts,
    channels: tuple[ChannelId, ChannelId],
    *,
    calibrations: CalibrationSet | None,
    policy: Policy | None = None,
    seed: int = 0,
    n_draws: int = N_MC,
) -> Corrected:
    """Corrected two-point visibility for `analysis` ("hom_visibility" or "bell_visibility").

    `high` is the acquisition with more signal (HOM reference delay, fringe maximum), `low` the
    other (HOM dip, fringe minimum).
    """
    if analysis not in ("hom_visibility", "bell_visibility"):
        raise ValueError(f"no forward model for {analysis!r}")
    policy = policy or Policy.default()
    corrections: list[str] = []
    labels: list[str] = [FORWARD_MODEL]
    dead = _dead_times(calibrations, channels, n_draws, seed)
    if dead is not None:
        policy.check_correction(analysis, "dead_time")
        corrections.append("dead_time")
    policy.check_correction(analysis, "accidental_subtraction")
    corrections.append("accidental_subtraction")
    if calibrations is not None:
        labels += [lab for lab in calibrations.labels if lab not in labels]

    rng = np.random.default_rng(np.random.SeedSequence([int(seed) & 0xFFFFFFFF, 0x46574D]))
    if dead is None:
        ta = tb = np.zeros(n_draws)
        pa = pb = np.array(0.0)
    else:
        ta, tb, pa, pb = dead
    if dead is not None and any(
        "dead_time_ps" in c.stand_in or "dead_time_ps" in c.unidentified
        for c in (calibrations.get(DARK_KIND, ch) for ch in channels)
    ):
        labels.append(DEAD_TIME_BOUND)
    for run in (high, low):
        for n_s, tau in ((run.singles_a, pa), (run.singles_b, pb)):
            if n_s / run.duration_s * float(tau) > MAX_DEAD_FRACTION and DEAD_FRACTION_HIGH not in labels:
                labels.append(DEAD_FRACTION_HIGH)

    def value(sh, sl):
        if analysis == "hom_visibility":
            return 1.0 - sl / sh
        return (sh - sl) / (sh + sl)

    def counts_draw(run):
        c = rng.gamma(run.coincidences + 0.5, size=n_draws)
        off = rng.gamma(run.offpeak_sum + 0.5, size=n_draws)
        return c, off

    ch, oh = counts_draw(high)
    cl, ol = counts_draw(low)
    joint = value(_signal(high, ch, oh, ta, tb), _signal(low, cl, ol, ta, tb))
    point = float(value(
        _signal(high, high.coincidences, high.offpeak_sum, pa, pb),
        _signal(low, low.coincidences, low.offpeak_sum, pa, pb),
    ))
    stat = value(_signal(high, ch, oh, pa, pb), _signal(low, cl, ol, pa, pb))
    calib = value(
        _signal(high, high.coincidences, high.offpeak_sum, ta, tb),
        _signal(low, low.coincidences, low.offpeak_sum, ta, tb),
    )
    joint = joint[np.isfinite(joint)]
    if joint.size < n_draws // 2:
        raise ValueError("accidentals leave no signal in most draws")
    lo, hi = (float(q) for q in np.quantile(joint, [0.025, 0.975]))
    est = Estimate(
        point,
        float(np.std(joint)),
        lo,
        hi,
        stderr_statistical=float(np.std(stat[np.isfinite(stat)])),
        stderr_calibration=float(np.std(calib)) if dead is not None else 0.0,
    )
    ids = calibrations.artifact_ids if calibrations is not None else ()
    return Corrected(est, tuple(corrections), tuple(labels), ids)


def bound_range(
    analysis: str,
    high: RunCounts,
    low: RunCounts,
    channels: tuple[ChannelId, ChannelId],
    *,
    calibrations: CalibrationSet,
    policy: Policy | None = None,
) -> Mapping[str, float] | None:
    """Decision D3: when a dead time is known only as an upper limit, the corrected value at zero
    and at the bound (point values), or None when every dead time is a posterior."""
    darks = [calibrations.get(DARK_KIND, ch) for ch in channels if any(m.kind == DARK_KIND for m in calibrations.members)]
    bounded = [c for c in darks if "dead_time_ps" in c.stand_in or "dead_time_ps" in c.unidentified]
    if not bounded:
        return None

    def at(scale: float) -> float:
        taus = []
        for c in darks:
            if c in bounded:
                taus.append(scale * c.estimates["dead_time_ps"].high * 1e-12)
            else:
                taus.append(c.estimates["dead_time_ps"].value * 1e-12)
        sh = _signal(high, high.coincidences, high.offpeak_sum, taus[0], taus[1])
        sl = _signal(low, low.coincidences, low.offpeak_sum, taus[0], taus[1])
        v = 1.0 - sl / sh if analysis == "hom_visibility" else (sh - sl) / (sh + sl)
        return float(v)

    return {"at_zero": at(0.0), "at_bound": at(1.0)}


__all__ = [
    "DEAD_FRACTION_HIGH",
    "DEAD_TIME_BOUND",
    "FORWARD_MODEL",
    "N_MC",
    "Corrected",
    "RunCounts",
    "bound_range",
    "corrected_visibility",
    "live_fraction",
]

