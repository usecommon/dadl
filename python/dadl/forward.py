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
from scipy import special

from dadl.calibration.dark import KIND as DARK_KIND
from dadl.calibration.delay import KIND as DELAY_KIND
from dadl.data import ChannelId
from dadl.errors import GuardrailRefusal
from dadl.guardrails import Policy
from dadl.propagate import CalibrationSet
from dadl.result import Estimate

#: Joint Monte Carlo draws per corrected value.
N_MC = 4000
#: Above this dead fraction per channel the first-order live fraction is not trusted.
MAX_DEAD_FRACTION = 0.2

FORWARD_MODEL = "corrected by the forward model (corrections listed in order), calibration uncertainty included"
DEAD_FRACTION_HIGH = "dead fraction above 20%: the live-fraction model is first order"
CAPTURE_UNMODELLED = (
    "window capture not modelled (no delay-jitter calibration for the channel pair): the fixed grid splits more "
    "pairs at the reference delay than at the dip, and V reads low"
)
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
    #: Arrival-time offsets (b - a, ps) of the signal pairs on top of the window's declared
    #: delays, as an equal-weight mixture: (0,) at a HOM dip or a fringe setting, (+d, -d) at a
    #: HOM path delay d (either photon may reach either output).
    pair_offsets_ps: tuple[float, ...] = (0.0,)

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


def _duration(ps: float) -> str:
    for unit, scale in (("ms", 1e9), ("us", 1e6), ("ns", 1e3)):
        if ps >= scale:
            return f"{ps / scale:.3g} {unit}"
    return f"{ps:.3g} ps"


def grid_capture(offset_ps, sigma_ps, width_ps: float) -> np.ndarray:
    """Probability that a pair with arrival difference D ~ N(offset, sigma) lands in one slot of a
    fixed grid of width W: E[max(0, 1 - |D| / W)], in closed form (vectorised)."""
    mu = np.asarray(offset_ps, dtype=float)
    sig = np.maximum(np.asarray(sigma_ps, dtype=float), 1e-9)
    w = float(width_ps)

    def part(lo, hi, sign):
        # integral over [lo, hi] of (1 - sign * d / W) N(d; mu, sig)
        a, b = (lo - mu) / sig, (hi - mu) / sig
        mass = special.ndtr(b) - special.ndtr(a)
        first = mu * mass - sig * (np.exp(-0.5 * b * b) - np.exp(-0.5 * a * a)) / np.sqrt(2 * np.pi)
        return mass - sign * first / w

    return part(0.0, w, 1.0) + part(-w, 0.0, -1.0)


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


def _signal(run: RunCounts, c, offsum, tau_a, tau_b, capture=1.0) -> np.ndarray:
    a = offsum / run.offpeak_shifts
    t = run.duration_s
    live = live_fraction(run.singles_a / t, run.singles_b / t, tau_a, tau_b, run.coincidences / t)
    return (c - a) / (t * live * capture)


def _capture(
    calibrations: CalibrationSet | None,
    channels: tuple[ChannelId, ChannelId],
    runs: tuple[RunCounts, RunCounts],
    window_ps: int | None,
    declared_delay_ps: float,
    n: int,
    seed: int,
):
    """Per-run grid capture (draws, point value) from the delay-jitter artifact of the channel
    pair, or None when there is none or no window was given."""
    if calibrations is None or window_ps is None:
        return None
    hits = [m for m in calibrations.members if m.kind == DELAY_KIND and set(m.channels) == set(channels)]
    if not hits:
        return None
    cal = hits[0].calibration
    sign = 1.0 if tuple(cal.channels) == tuple(channels) else -1.0
    joint = calibrations.draws(n, seed)
    delay = sign * joint.column(DELAY_KIND, "delay_ps", channels[0])
    jitter = joint.column(DELAY_KIND, "jitter_sigma_ps", channels[0])
    d0 = sign * float(cal.estimates["delay_ps"].value)
    j0 = float(cal.estimates["jitter_sigma_ps"].value)
    out = []
    for run in runs:
        offs = run.pair_offsets_ps
        draws = sum(grid_capture(o + delay - declared_delay_ps, jitter, window_ps) for o in offs) / len(offs)
        point = sum(float(grid_capture(o + d0 - declared_delay_ps, j0, window_ps)) for o in offs) / len(offs)
        out.append((draws, point))
    return out


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
    window_ps: int | None = None,
    declared_delay_ps: float = 0.0,
) -> Corrected:
    """Corrected two-point visibility for `analysis` ("hom_visibility" or "bell_visibility").

    `high` is the acquisition with more signal (HOM reference delay, fringe maximum), `low` the
    other (HOM dip, fringe minimum). With `window_ps` (the fixed-grid width) and a delay-jitter
    artifact for the channel pair, each acquisition's grid capture is modelled from the
    calibrated delay and jitter, `declared_delay_ps` (b - a, the window's) and the runs'
    `pair_offsets_ps`.
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
    capture = _capture(calibrations, channels, (high, low), window_ps, declared_delay_ps, n_draws, seed)
    if capture is not None:
        policy.check_correction(analysis, "window_capture")
        corrections.append("window_capture")
    elif any(any(o != 0 for o in run.pair_offsets_ps) for run in (high, low)):
        labels.append(CAPTURE_UNMODELLED)
    if calibrations is not None:
        labels += [lab for lab in calibrations.labels if lab not in labels]

    rng = np.random.default_rng(np.random.SeedSequence([int(seed) & 0xFFFFFFFF, 0x46574D]))
    if dead is None:
        ta = tb = np.zeros(n_draws)
        pa = pb = np.array(0.0)
    else:
        ta, tb, pa, pb = dead
    if dead is not None:
        bounded = [
            (ch, float(c.estimates["dead_time_ps"].high))
            for ch, c in ((ch, calibrations.get(DARK_KIND, ch)) for ch in channels)
            if "dead_time_ps" in c.stand_in or "dead_time_ps" in c.unidentified
        ]
        if bounded:
            labels.append(DEAD_TIME_BOUND)
        for ch, hi_ps in bounded:
            side = "singles_a" if ch == channels[0] else "singles_b"
            worst = max(getattr(run, side) / run.duration_s for run in (high, low)) * hi_ps * 1e-12
            if worst > MAX_DEAD_FRACTION:
                # A bound this loose spans "no correction" to a dead fraction the first-order live
                # fraction cannot describe (above 1 it is negative): the range carries no information.
                raise GuardrailRefusal(
                    f"the dead time of {ch} is known only below {_duration(hi_ps)} (the dark data did not identify "
                    f"it); at the recorded rate that bound is a dead fraction of {worst:.3g}, above "
                    f"{MAX_DEAD_FRACTION:.0%}, so a correction between zero and the bound says nothing. Calibrate the "
                    "dead time from data that identify it, or load the calibrations without the dark artifacts "
                    f"(CalibrationSet.load(..., kinds=...)) to correct the accidentals only."
                )
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

    (kh, kh0), (kl, kl0) = capture if capture is not None else ((1.0, 1.0), (1.0, 1.0))
    ch, oh = counts_draw(high)
    cl, ol = counts_draw(low)
    joint = value(_signal(high, ch, oh, ta, tb, kh), _signal(low, cl, ol, ta, tb, kl))
    point = float(value(
        _signal(high, high.coincidences, high.offpeak_sum, pa, pb, kh0),
        _signal(low, low.coincidences, low.offpeak_sum, pa, pb, kl0),
    ))
    stat = value(_signal(high, ch, oh, pa, pb, kh0), _signal(low, cl, ol, pa, pb, kl0))
    calib = value(
        _signal(high, high.coincidences, high.offpeak_sum, ta, tb, kh),
        _signal(low, low.coincidences, low.offpeak_sum, ta, tb, kl),
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
        stderr_calibration=float(np.std(calib)) if dead is not None or capture is not None else 0.0,
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
    window_ps: int | None = None,
    declared_delay_ps: float = 0.0,
) -> Mapping[str, float] | None:
    """Decision D3: when a dead time is known only as an upper limit, the corrected value at zero
    and at the bound (point values), or None when every dead time is a posterior."""
    darks = [calibrations.get(DARK_KIND, ch) for ch in channels if any(m.kind == DARK_KIND for m in calibrations.members)]
    bounded = [c for c in darks if "dead_time_ps" in c.stand_in or "dead_time_ps" in c.unidentified]
    if not bounded:
        return None

    cap = _capture(calibrations, channels, (high, low), window_ps, declared_delay_ps, 1, 0)
    kh, kl = (cap[0][1], cap[1][1]) if cap is not None else (1.0, 1.0)

    def at(scale: float) -> float:
        taus = []
        for c in darks:
            if c in bounded:
                taus.append(scale * c.estimates["dead_time_ps"].high * 1e-12)
            else:
                taus.append(c.estimates["dead_time_ps"].value * 1e-12)
        sh = _signal(high, high.coincidences, high.offpeak_sum, taus[0], taus[1], kh)
        sl = _signal(low, low.coincidences, low.offpeak_sum, taus[0], taus[1], kl)
        v = 1.0 - sl / sh if analysis == "hom_visibility" else (sh - sl) / (sh + sl)
        return float(v)

    return {"at_zero": at(0.0), "at_bound": at(1.0)}


__all__ = [
    "CAPTURE_UNMODELLED",
    "DEAD_FRACTION_HIGH",
    "DEAD_TIME_BOUND",
    "FORWARD_MODEL",
    "N_MC",
    "Corrected",
    "RunCounts",
    "bound_range",
    "corrected_visibility",
    "grid_capture",
    "live_fraction",
]

