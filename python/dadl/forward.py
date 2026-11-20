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

The calibration part is conditional on the models: it includes the calibration artifacts'
statistical uncertainty only. It excludes each artifact's `uncertainty.not_included` (copied into
the record as `calibration_excludes`), the IRF shape beyond the Gaussian the delay check accepted
(HOM and fringe capture), efficiencies, afterpulsing, multi-pair emission, the splitter ratio,
drift between calibration and run below the drift threshold, and the choice among several valid
artifacts.

Corrections applied, in order: dead time (if dark artifacts are given), then accidentals. Each is
checked against the policy for the analysis; dead time is not in the default HOM rule, so it
needs `Policy.allow(analysis, "dead_time", reason)`.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np
from scipy import special

from dadl.calibration.dark import KIND as DARK_KIND
from dadl.calibration.delay import KIND as DELAY_KIND
from dadl.coincidences import as_sorted as _as_sorted
from dadl.data import ChannelId
from dadl.errors import GuardrailRefusal
from dadl.guardrails import Policy
from dadl.propagate import CalibrationSet
from dadl.result import Estimate

#: Joint Monte Carlo draws per corrected value.
N_MC = 4000
#: Above this dead fraction per channel the pair term and the unmodelled paralyzable stage are not
#: trusted, and a dark dead time this long is often a false identification.
MAX_DEAD_FRACTION = 0.2

FORWARD_MODEL = "corrected by the forward model (corrections listed in order), calibration uncertainty included"
#: Kept for old records only: an identified dead time above this dead fraction is now refused (council round 6).
DEAD_FRACTION_HIGH = "dead fraction above 20%: the live-fraction model is first order"
#: The analysed gaps fall below a loose dark dead time's 95% lower end (council round 6 addendum):
#: labelled, not refused, since that lower end misses in 5% of channels by design.
DEAD_TIME_BELOW_DARK = (
    "the analysed gaps put the dead time of {name} below the dark artifact's 95% lower end (a 1-in-20 chance per "
    "channel for a dark run this short): its dead-time correction is likely too large; a longer dark run, or an "
    "illuminated calibration for an SNSPD, resolves it"
)
#: The dark dead time is longer than many analysed gaps (council round 8, Ph4). It fires on most
#: clean analyses, so it says how often (council round 10, Rg2), and it is left out for a channel that
#: already carries `DEAD_TIME_BELOW_DARK` or a gap problem.
DEAD_TIME_ABOVE_GAPS = (
    "the dark dead time of {name} ({x}) is longer than {k} of the analysed gaps (shortest {g}); this fires on most "
    "analyses whose dark run is short against the analysed one (council round 8, 400 of 400 chsh-dead-time "
    "replicates). The correction uses the dark value and its uncertainty (see Known limits)"
)
#: The same for a dark dead time known only as an upper limit (council round 10, Rg2).
DEAD_TIME_ABOVE_GAPS_BOUND = (
    "the dark run bounds the dead time of {name} below {x}; {k} analysed gaps are shorter, so the bound is loose; "
    "the correction uses the bound (see Known limits)"
)
#: HOM and fringe forward-model values (council round 6, R3), as CHSH's model label says for S.
SECONDARY_V = "secondary; the headline is raw V"
CAPTURE_UNMODELLED = (
    "window capture not modelled (no delay-jitter calibration for the channel pair): the fixed grid splits more "
    "pairs at the reference delay than at the dip, and V reads low"
)
DEAD_TIME_BOUND = (
    "dead time known only as an upper limit: value and interval at the bound; counts give the value at zero "
    "(the truth lies between if the bound holds; a bound set by a few short dark gaps may not)"
)
#: The override of the analysis-time gap tests (council round 9, Rd2), and its label.
DEAD_TIME_OVERRIDE = "dead_time_contradicted"
DEAD_TIME_CONTRADICTED = "dead time contradicted by the analysed gaps, accepted: {problems}"
#: A dark dead time below this fraction of the setup's declared one is labelled (council round 9, Rd1b).
DEAD_TIME_DECLARED_FRACTION = 0.5
DEAD_TIME_FAR_BELOW = (
    "the dark dead time of {name} ({x}) is below half the declared {y}: a short double trigger may have set it; "
    "check the dark artifact's notes"
)
DEAD_TIME_RANGE_WIDE = (
    "the range from zero to the dead-time bound moves the corrected value by more than its standard error: "
    "quote the range"
)
JITTER_BOUND = "delay jitter known only as an upper limit: window capture at the bound"
#: HOM and fringe without dark artifacts, in CHSH's words (`chsh_model_label`; council round 10, Ux12).
DEAD_TIME_UNMODELLED = "dead time (no dark artifacts) not modelled"
CORRECTED_ABOVE_ONE = "corrected value above 1: the accidentals fluctuate above the minimum count; quote the interval"
AFTERPULSING_UNMODELLED = (
    "afterpulsing reported by a dark calibration is not in the forward model (no kernel propagated)"
)
#: Relative change of a window capture between jitter 0 and its bound above which the bound is labelled.
JITTER_BOUND_RTOL = 1e-3


def _policy_labels(policy: Policy) -> list[str]:
    """`POLICY_MODIFIED` on a corrected value under a modified policy (council round 10, Ux2: the
    raw value says it is unchanged instead, `Policy.labels_for`)."""
    from dadl.guardrails import POLICY_MODIFIED

    return [POLICY_MODIFIED] if policy.is_modified else []


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
    """Probability that a pair arriving in the window is recorded on both channels.

    The forward model takes each channel's non-paralyzable dead time from the dark artifact. At
    high analysed rates the analysed gaps measure it more precisely, and often show it to be
    shorter than the dark value (label `DEAD_TIME_ABOVE_GAPS`). The pair live fraction is first
    order: the term r_c·tau·(fa + fb) is not included. In the chsh-dead-time regime the two
    effects cancel: the corrected S covers at 95.5% with mean z -0.43 statistical se. Neither
    effect is corrected alone, because each alone biases the value more (council round 8).
    Reading the dead time from the analysed gaps, together with the second-order term, is
    milestone M5a."""
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
    #: Added to the result's counts: {"monte_carlo": {"seed", "n_draws"}} (council round 5, V3) and
    #: `calibration_excludes`, the used artifacts' `uncertainty.not_included` (council round 6, R4).
    counts: Mapping[str, object] = field(default_factory=dict)


def _mc_counts(seed: int, n_draws: int, calibrations: CalibrationSet | None = None, ids=()) -> dict[str, object]:
    """The Monte Carlo seed and draw count (the corrected interval depends on both), and what the
    used artifacts say their uncertainty leaves out (`calibration_excludes`)."""
    out: dict[str, object] = {"monte_carlo": {"seed": int(seed), "n_draws": int(n_draws)}}
    if calibrations is not None:
        excluded: list[str] = []
        for m in calibrations.members:
            if m.artifact_id in ids:
                excluded += list(((m.manifest or {}).get("uncertainty") or {}).get("not_included") or ())
        out["calibration_excludes"] = list(dict.fromkeys(excluded))
    return out


def _calibration_labels(calibrations: CalibrationSet | None, ids) -> list[str]:
    """The set's labels that concern the artifacts the value used (council round 6, R2)."""
    if calibrations is None:
        return []
    return list(calibrations.labels_for(ids))


def _signal(run: RunCounts, c, offsum, tau_a, tau_b, capture=1.0) -> np.ndarray:
    a = offsum / run.offpeak_shifts
    t = run.duration_s
    live = live_fraction(run.singles_a / t, run.singles_b / t, tau_a, tau_b, run.coincidences / t)
    return (c - a) / (t * live * capture)


def _delay_jitter(cal, joint, n: int):
    """(delay draws, jitter draws, delay point, jitter point, jitter is a bound) of a delay-jitter
    artifact, in its own channel order. A stand-in jitter (a coarse tagger's upper limit) is not
    drawn: it is used at the bound."""
    key = tuple(cal.channels)
    delay = joint.columns[(DELAY_KIND, "delay_ps", key)]
    if "jitter_sigma_ps" in cal.stand_in or "jitter_sigma_ps" in cal.unidentified:
        hi = float(cal.estimates["jitter_sigma_ps"].high)
        return delay, np.full(n, hi), float(cal.estimates["delay_ps"].value), hi, True
    jitter = joint.columns[(DELAY_KIND, "jitter_sigma_ps", key)]
    return delay, jitter, float(cal.estimates["delay_ps"].value), float(cal.estimates["jitter_sigma_ps"].value), False


def _bound_matters(at_bound: float, at_zero: float) -> bool:
    return abs(at_bound - at_zero) > JITTER_BOUND_RTOL * max(abs(at_bound), 1e-12)


def _capture(
    calibrations: CalibrationSet | None,
    channels: tuple[ChannelId, ChannelId],
    runs: tuple[RunCounts, RunCounts],
    window_ps: int | None,
    declared_delay_ps: float,
    n: int,
    seed: int,
):
    """Per-run grid capture (draws, point value, jitter bound matters) from the delay-jitter
    artifact of the channel pair, or None when there is none or no window was given."""
    if calibrations is None or window_ps is None:
        return None
    hits = [m for m in calibrations.members if m.kind == DELAY_KIND and set(m.channels) == set(channels)]
    if not hits:
        return None
    cal = hits[0].calibration
    sign = 1.0 if tuple(cal.channels) == tuple(channels) else -1.0
    delay, jitter, d0, j0, bounded = _delay_jitter(cal, calibrations.draws(n, seed), n)
    delay, d0 = sign * delay, sign * d0
    out = []
    for run in runs:
        offs = run.pair_offsets_ps
        draws = sum(grid_capture(o + delay - declared_delay_ps, jitter, window_ps) for o in offs) / len(offs)
        point = sum(float(grid_capture(o + d0 - declared_delay_ps, j0, window_ps)) for o in offs) / len(offs)
        wide = False
        if bounded:
            zero = sum(float(grid_capture(o + d0 - declared_delay_ps, 0.0, window_ps)) for o in offs) / len(offs)
            wide = _bound_matters(point, zero)
        out.append((draws, point, wide))
    return out


def _dead_times(
    calibrations: CalibrationSet | None, channels: tuple[ChannelId, ChannelId], n: int, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Draws and point values of the two dead times (s), or None without dark artifacts.

    The point is each dark artifact's reported dead time, the shortest gap: an upper edge, high by
    about one mean gap divided by the number of gaps. The interval comes from the posterior draws,
    so the corrected point sits low in it (about 0.4 sd for S in chsh-dead-time). A dead time known
    only as an upper limit is used at the bound, for the point and every draw (decision D3)."""
    # One helper for every path: a partial dark set raises the named ValueError (council round 4).
    dead = _channel_dead_times(calibrations, {"A": channels[0], "B": channels[1]}, n, seed)
    if dead is None:
        return None
    (da, pa, _), (db, pb, _) = dead["A"], dead["B"]
    return da, db, np.array(pa), np.array(pb)


#: Without a recorded gap count the floor sits this many times the interval's lower half-width
#: below the point (artifacts written before the count was recorded).
GAP_FLOOR_WIDTHS = 3.0
#: Per-analysis false-refusal rates of the downward and upward gap tests (council round 8: their
#: share of the 1e-3 per-analysis refusal budget, each Bonferroni over the analysed channels).
GAP_DOWN_ALPHA = 1e-4
GAP_UP_ALPHA = 1e-4
#: More gaps below the floor than max(GAP_COUNT, GAP_SHARE x gaps) refuse the dead time.
GAP_COUNT = 10
GAP_SHARE = 1e-3
#: The edge spacing is read from this many shortest gaps (the test needs one more).
GAP_EDGE_GAPS = 100


def gap_floor_ps(point_ps: float, gaps: float, dark_rate_low_cps: float, channels: int) -> float:
    """The downward floor of a dark dead time (council round 8, St1): point - ln(K / alpha) / (N f),
    with N the dark run's gaps, f the lower end of its dark rate and K the analysed channels. The
    shortest of N dark gaps sits above tau by Exp(N f), and afterpulsing only raises the edge
    density, so the floor misses the truth with chance at most alpha / K per channel."""
    if gaps <= 0 or dark_rate_low_cps <= 0:
        return -math.inf
    return point_ps - math.log(channels / GAP_DOWN_ALPHA) / (gaps * dark_rate_low_cps) * 1e12


def gap_edge_excess(channels: int, j: int = 0, m: int = GAP_EDGE_GAPS) -> float:
    """The upward test's threshold (council round 9, Rd1): (j + 1) F.isf(alpha / K, 2 (j + 1), 2 m).

    Under a memoryless edge (gaps = tau + Exp) the excess (g[j] - tau) m / (g[j + m] - g[j]) is
    m Gamma(j + 1) / Gamma(m), that is (j + 1) times an F(2 (j + 1), 2 m) variable, so a matched
    upper end exceeds the threshold with chance alpha / K per channel. At j = 0 and m = 100 it is
    100 ((K / alpha)^(1/100) - 1) (10.4 at K = 2, 11.2 at K = 4)."""
    from scipy import stats as _st

    return (j + 1) * float(_st.f.isf(GAP_UP_ALPHA / channels, 2 * (j + 1), 2 * m))


def upward_gap_excess(gaps: np.ndarray, upper: float) -> tuple[float, int, int] | None:
    """The upward gap test's statistic on sorted gaps (council round 9, Rd1): (excess, j, M) with
    j = max(`GAP_COUNT`, `GAP_SHARE` x gaps), M = max(`GAP_EDGE_GAPS`, j + 1) and excess =
    (g[j] - upper) M / (g[j + M] - g[j]), infinite at a zero spacing; None when there are too few
    gaps or g[j] is not above `upper`. Compare it with `gap_edge_excess` (K, j, M)."""
    j = int(max(GAP_COUNT, GAP_SHARE * gaps.size))
    m = max(GAP_EDGE_GAPS, j + 1)
    if gaps.size <= j + m or not gaps[j] > upper:
        return None
    spread = float(gaps[j + m] - gaps[j])
    return ((float(gaps[j]) - upper) * m / spread if spread > 0 else math.inf), j, m


def dead_time_far_below(cal, declared_ps: float) -> float | None:
    """The dark artifact's dead time (point, or bound when not identified) when it is below
    `DEAD_TIME_DECLARED_FRACTION` of the setup's declared dead time, else None. One rule for the
    CLI's calibrate warning and analyze label (council round 9, Rd1b)."""
    if not declared_ps or declared_ps <= 0:
        return None
    est = cal.estimates["dead_time_ps"]
    bound = "dead_time_ps" in cal.stand_in or "dead_time_ps" in cal.unidentified
    used = float(est.high) if bound else float(est.value)
    return used if used < DEAD_TIME_DECLARED_FRACTION * float(declared_ps) else None


def check_dead_time_gaps(
    calibrations: CalibrationSet | None,
    channels: Mapping[str, ChannelId],
    runs,
    *,
    analysis: str,
    policy: Policy | None = None,
) -> tuple[tuple[str, ...], dict[str, object]]:
    """Compare each dark dead time with the analysed tags: (labels, counts).

    A non-paralyzable dead time tau leaves no gap shorter than tau on its channel, and the shortest
    gaps sit just above it. The check refuses a dark dead time that the analysed gaps contradict in
    either direction, at most `GAP_DOWN_ALPHA` + `GAP_UP_ALPHA` per analysis on clean data (the
    README's refusal budget).

    Downward: more than max(`GAP_COUNT`, `GAP_SHARE` x gaps) of the analysed gaps (each
    acquisition in `runs`, a dict of channel name -> tags) below the floor `gap_floor_ps` (another
    detector's dark artifact gave a 1 us dead time to a run with 30 ns gaps; council round 3). The
    floor comes from the dark run's own gap count and dark rate (`conditions["gaps"]`), not from
    the interval's width: the 3-width floor missed 2.5e-4 to 3.7e-4 per channel when the dark fit
    falsely found afterpulsing (council round 8, St1). It is about 4 ns below the point for a 1 s
    dark run at 50 kcps and 265 ns below it for 0.1 s at 20 kcps, where the dark data cannot
    exclude anything. An artifact without the gap count keeps the 3-width floor. A dead time known
    only as an upper limit has no floor; a floor at or below zero is not refused. The 95% lower end
    is then tested for a label only (`DEAD_TIME_BELOW_DARK`, same count rule): it misses in 5% of
    channels by design. An identified dark dead time can be wrong when the detector's dead time is
    far below the spacing of dark counts (an SNSPD at 100 cps: 3-7% of dark runs reports a
    microsecond-scale dead time as identified). The analysis refuses it when it implies a dead
    fraction above 20% at the analysed rate (`corrected_visibility`, `corrected_chsh`), and labels
    it when the analysed gaps fall below its 95% lower end. A lower end at zero is not checked.

    Upward: with g the sorted analysed gaps, j = max(`GAP_COUNT`, `GAP_SHARE` x gaps) (the
    downward count rule) and M = max(`GAP_EDGE_GAPS`, j + 1), the excess is (g[j] - upper end) M /
    (g[j + M] - g[j]); the dead time is refused when it exceeds `gap_edge_excess` (K, j, M) (a 22 ns
    artifact on a 50 ns detector; council round 4). This applies to upper limits too. A zero spacing
    is an infinite excess, and a NaN never passes. On clean data the upper end is the shortest dark
    gap, never below the truth, so the test fires far less often than its bound. Reading g[j], not
    g[0], makes it robust to j short gaps: one echo (a cable reflection) set g[0] to the dark bound
    and turned the test off (council round 9, Rd1). The price is power at low rates: at 10 kcps and
    50,000 gaps the round-4 case passes (excess 80 against 100; g[0] gave 15 against 11.2).

    The test tolerates a few unusual gaps, up to 0.1% of a channel's gaps or 10, whichever is more.
    A rarer double trigger (a cable reflection, a ringing discriminator) cannot hide the edge. A
    more frequent one makes the edge look shorter than the detector's dead time. If the setup
    declares the detector's dead time (`detectors.<name>.dead_time_ps`), the CLI labels a dark dead
    time far below it (`dead_time_far_below`). The dark artifact's notes say when a single short gap
    set its bound.

    Above the gaps: a dark dead time (point, or bound) longer than more than max(`GAP_COUNT`,
    `GAP_SHARE` x gaps) of the analysed gaps is labelled `DEAD_TIME_ABOVE_GAPS` (council round 8,
    Ph4), or `DEAD_TIME_ABOVE_GAPS_BOUND` for a bound; the correction keeps the dark value (README,
    Known limits). The label is left out on a channel that already carries `DEAD_TIME_BELOW_DARK` or a
    gap problem (council round 10, Rg2).

    A refused dead time raises `StaleCalibration`, unless the policy allows
    `dead_time_contradicted` (then `DEAD_TIME_CONTRADICTED`, with the problems; council round 9,
    Rd2: `calibration_outside_validity` no longer covers it). The problems are also recorded in
    `dead_time_gap_problems` when accepted. The counts carry the shortest
    analysed gap (`min_gap_ps`), the floor (`gap_check_floor_ps`) and which floor it is
    (`gap_check_floor_kind`: "exponential law" or "3 widths" (refusing), "lower end (label only)",
    or None when not checked downward) per channel."""
    from dadl.errors import StaleCalibration
    from dadl.guardrails import override_hint

    if calibrations is None or not any(m.kind == DARK_KIND for m in calibrations.members):
        return (), {}
    policy = policy or Policy.default()
    darks = {}
    for name, cid in channels.items():
        try:
            darks[name] = calibrations.get(DARK_KIND, cid)
        except KeyError:
            continue
    k_channels = max(len(darks), 1)
    min_gap: dict[str, float] = {}
    floors: dict[str, float | None] = {}
    kinds: dict[str, str | None] = {}
    problems: list[str] = []
    below: list[str] = []
    for name, cal in darks.items():
        cid = channels[name]
        gaps = np.concatenate([np.diff(_as_sorted(np.asarray(run[name], dtype=np.int64))) for run in runs] or [np.zeros(0)])
        if gaps.size == 0:
            continue
        gaps = np.sort(gaps)
        min_gap[name] = float(gaps[0])
        tau = cal.estimates["dead_time_ps"]
        bound = "dead_time_ps" in cal.stand_in or "dead_time_ps" in cal.unidentified
        count_rule = max(GAP_COUNT, GAP_SHARE * gaps.size)
        floor, kind = None, None
        if not bound:
            n_dark = (cal.conditions or {}).get("gaps")
            rate = cal.estimates.get("dark_rate_cps")
            if n_dark is not None and rate is not None:
                floor = gap_floor_ps(float(tau.value), float(n_dark), float(rate.low), k_channels)
                kind = "exponential law"
            else:
                floor, kind = float(tau.value) - GAP_FLOOR_WIDTHS * (float(tau.value) - float(tau.low)), "3 widths"
            if floor <= 0:
                floor, kind = float(tau.low), "lower end (label only)"
            if floor <= 0:
                floor, kind = None, None
        floors[name], kinds[name] = floor, kind
        flagged = False  # this channel already carries DEAD_TIME_BELOW_DARK or a gap problem
        if floor is not None:
            k = int(np.searchsorted(gaps, floor))
            if k > count_rule:
                flagged = True
                if kind != "lower end (label only)":
                    problems.append(
                        f"dead time {tau.value / 1000:.4g} ns of {name} ({cid}): {k} analysed gaps are shorter, "
                        f"down to {gaps[0] / 1000:.4g} ns"
                    )
                else:
                    below.append(DEAD_TIME_BELOW_DARK.format(name=f"{name} ({cid})"))
        upper = float(tau.high)
        used = upper if bound else float(tau.value)
        above = int(np.searchsorted(gaps, used))
        edge = upward_gap_excess(gaps, upper)
        if edge is not None:
            excess, j, m = edge
            if not excess <= gap_edge_excess(k_channels, j, m):
                what = f"dead time {tau.value / 1000:.4g}" if not bound else f"dead time below {upper / 1000:.4g}"
                problems.append(
                    f"{what} ns of {name} ({cid}): the analysed gap {j + 1} from the shortest is "
                    f"{gaps[j] / 1000:.4g} ns, longer than its upper end by {excess:.0f} mean edge spacings"
                )
                flagged = True
        if above > count_rule and not flagged:
            text = DEAD_TIME_ABOVE_GAPS_BOUND if bound else DEAD_TIME_ABOVE_GAPS
            below.append(text.format(name=f"{name} ({cid})", x=_duration(used), k=above, g=_duration(float(gaps[0]))))
    labels: tuple[str, ...] = ()
    counts: dict[str, object] = {}
    if problems:
        if DEAD_TIME_OVERRIDE not in policy.rule_for(analysis).allowed_corrections:
            raise StaleCalibration(
                f"the dark artifacts do not describe the analysed detectors: {'; '.join(problems)}. "
                f"{override_hint(analysis, DEAD_TIME_OVERRIDE)}"
            )
        labels = (DEAD_TIME_CONTRADICTED.format(problems="; ".join(problems)),)
        counts["dead_time_gap_problems"] = list(problems)
    labels += tuple(below)
    if min_gap:
        counts["min_gap_ps"] = min_gap
        counts["gap_check_floor_ps"] = floors
        counts["gap_check_floor_kind"] = kinds
    return labels, counts


def _used(
    calibrations: CalibrationSet | None, dark_channels, pairs
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """The artifact IDs the forward model read (dark artifacts of `dark_channels`, delay-jitter
    artifacts of the channel `pairs`) and the labels they carry."""
    if calibrations is None:
        return (), ()
    ids: list[str] = []
    labels: list[str] = []
    for m in calibrations.members:
        if m.kind == DARK_KIND and set(m.channels) & set(dark_channels):
            ids.append(m.artifact_id)
            cal = m.calibration
            if "afterpulse_prob" in cal.estimates and not (
                "afterpulse_prob" in cal.stand_in or "afterpulse_prob" in cal.unidentified
            ):
                labels.append(AFTERPULSING_UNMODELLED)
        elif m.kind == DELAY_KIND and any(set(m.channels) == set(p) for p in pairs):
            ids.append(m.artifact_id)
    return tuple(ids), tuple(dict.fromkeys(labels))


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
    labels: list[str] = [FORWARD_MODEL, SECONDARY_V] + _policy_labels(policy)
    dead = _dead_times(calibrations, channels, n_draws, seed)
    if dead is not None:
        policy.check_correction(analysis, "dead_time")
        corrections.append("dead_time")
    else:
        labels.append(DEAD_TIME_UNMODELLED)
    policy.check_correction(analysis, "accidental_subtraction")
    corrections.append("accidental_subtraction")
    capture = _capture(calibrations, channels, (high, low), window_ps, declared_delay_ps, n_draws, seed)
    if capture is not None:
        policy.check_correction(analysis, "window_capture")
        corrections.append("window_capture")
        if any(wide for _, _, wide in capture):
            labels.append(JITTER_BOUND)
    elif any(any(o != 0 for o in run.pair_offsets_ps) for run in (high, low)):
        labels.append(CAPTURE_UNMODELLED)

    rng = np.random.default_rng(np.random.SeedSequence([int(seed) & 0xFFFFFFFF, 0x46574D]))
    if dead is None:
        ta = tb = np.zeros(n_draws)
        pa = pb = np.array(0.0)
    else:
        ta, tb, pa, pb = dead
    if dead is not None:
        any_bound = False
        for ch, side in ((channels[0], "singles_a"), (channels[1], "singles_b")):
            c = calibrations.get(DARK_KIND, ch)
            rate = max(getattr(run, side) / run.duration_s for run in (high, low))
            if "dead_time_ps" in c.stand_in or "dead_time_ps" in c.unidentified:
                any_bound = True
                _refuse_loose_bound(str(ch), None, float(c.estimates["dead_time_ps"].high), rate)
            else:
                _refuse_dead_fraction(str(ch), None, float(c.estimates["dead_time_ps"].value), rate)
        if any_bound:
            labels.append(DEAD_TIME_BOUND)

    def value(sh, sl):
        if analysis == "hom_visibility":
            return 1.0 - sl / sh
        return (sh - sl) / (sh + sl)

    def counts_draw(run):
        c = rng.gamma(run.coincidences + 0.5, size=n_draws)
        off = rng.gamma(run.offpeak_sum + 0.5, size=n_draws)
        return c, off

    (kh, kh0, _), (kl, kl0, _) = capture if capture is not None else ((1.0, 1.0, False), (1.0, 1.0, False))
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
    ids, used_labels = _used(calibrations, channels if dead is not None else (), [channels] if capture is not None else [])
    labels += [lab for lab in _calibration_labels(calibrations, ids) + list(used_labels) if lab not in labels]
    lo, hi = (float(q) for q in np.quantile(joint, [0.025, 0.975]))
    if point > 1.0:
        # Kept unclipped: clipping would hide how far the accidentals fluctuated.
        labels.append(CORRECTED_ABOVE_ONE)
    est = Estimate(
        point,
        float(np.std(joint)),
        lo,
        hi,
        stderr_statistical=float(np.std(stat[np.isfinite(stat)])),
        stderr_calibration=float(np.std(calib)) if dead is not None or capture is not None else 0.0,
    )
    return Corrected(est, tuple(corrections), tuple(labels), ids, _mc_counts(seed, n_draws, calibrations, ids))


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
    dead = _channel_dead_times(calibrations, {"A": channels[0], "B": channels[1]}, 1, 0)
    if dead is None:
        return None
    darks = [calibrations.get(DARK_KIND, ch) for ch in channels]
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


# --- CHSH ----------------------------------------------------------------------------------------

#: Signs of the four outcome combinations in E = (N++ + N-- - N+- - N-+) / N.
OUTCOME_SIGNS: Mapping[tuple[str, str], int] = {("+", "+"): 1, ("-", "-"): 1, ("+", "-"): -1, ("-", "+"): -1}
#: What a corrected S models, built from its corrections (`chsh_model_label`; council round 3: the
#: text said dead time and window capture were modelled when neither was applied).
ACCIDENTALS_MODELLED = (
    "{modelled} modelled; {unmodelled}detector efficiencies not (an output imbalance biases corrected S as it does "
    "raw S); secondary, not a Bell claim (the headline is raw S)"
)
_CHSH_PARTS = (
    ("accidental_subtraction", "accidentals", None),
    ("dead_time", "dead time", "dead time (no dark artifacts)"),
    ("window_capture", "window capture", "window capture (no delay-jitter artifacts)"),
)


def _and(words: list[str]) -> str:
    return ", ".join(words[:-1]) + " and " + words[-1] if len(words) > 1 else "".join(words)


def chsh_model_label(corrections) -> str:
    """`ACCIDENTALS_MODELLED` for a corrected S with these corrections."""
    done = [name for key, name, _ in _CHSH_PARTS if key in corrections]
    missing = [why for key, _, why in _CHSH_PARTS if why and key not in corrections]
    return ACCIDENTALS_MODELLED.format(
        modelled=_and(done), unmodelled=f"{_and(missing)} not modelled; " if missing else ""
    )


def _channel_dead_times(
    calibrations: CalibrationSet | None, channels: Mapping[str, ChannelId], n: int, seed: int
) -> dict[str, tuple[np.ndarray, float, float | None]] | None:
    """Per channel name: (draws in s, point in s, upper limit in ps or None), or None without
    dark artifacts. With any dark artifact, every channel needs one.

    The point is each dark artifact's reported dead time, the shortest gap: an upper edge, high by
    about one mean gap divided by the number of gaps. The interval comes from the posterior draws,
    so the corrected point sits low in it (about 0.4 sd for S in chsh-dead-time). A dead time known
    only as an upper limit is used at the bound (decision D3).

    The dead time is the dark artifact's, even where the analysed gaps show it to be shorter
    (`DEAD_TIME_ABOVE_GAPS`); see `live_fraction` for why it is not truncated at the shortest
    analysed gap (council round 8, Ph4 and Ph5; milestone M5a)."""
    if calibrations is None or not any(m.kind == DARK_KIND for m in calibrations.members):
        return None
    joint = calibrations.draws(n, seed)
    out = {}
    for name, cid in channels.items():
        try:
            cal = calibrations.get(DARK_KIND, cid)
        except KeyError as exc:
            raise ValueError(f"dark artifacts are given but none for {name} ({cid}): every channel needs one") from exc
        if "dead_time_ps" in cal.stand_in or "dead_time_ps" in cal.unidentified:
            hi = float(cal.estimates["dead_time_ps"].high)
            out[name] = (np.full(n, hi * 1e-12), hi * 1e-12, hi)
        else:
            out[name] = (
                joint.column(DARK_KIND, "dead_time_ps", cid) * 1e-12,
                float(cal.estimates["dead_time_ps"].value) * 1e-12,
                None,
            )
    return out


#: How to correct the rest when the dark dead time cannot be used (council round 6, R7).
WITHOUT_DARK = (
    "load the calibrations without the dark artifacts (CalibrationSet.load(..., kinds=...); on the CLI, move the "
    "dark artifacts out of the calibrations folder, each .json with its .posterior.npy) to correct the rest"
)


def _who(name: str, cid: ChannelId | None) -> str:
    return name if cid is None else f"{name} ({cid})"


def _refuse_loose_bound(name: str, cid: ChannelId | None, hi_ps: float, rate_cps: float) -> None:
    """A bound this loose spans "no correction" to a dead fraction the first-order live fraction
    cannot describe (above 1 it is negative): the range carries no information."""
    worst = rate_cps * hi_ps * 1e-12
    if worst > MAX_DEAD_FRACTION:
        raise GuardrailRefusal(
            f"the dead time of {_who(name, cid)} is known only below {_duration(hi_ps)} (the dark data did not "
            f"identify it); at the recorded rate that bound gives rate x dead time = {worst:.3g} (limit "
            f"{MAX_DEAD_FRACTION:g}), so a correction between zero and the bound says nothing. No DADL calibration "
            f"can supply it yet (series dead-time artifacts are not used by analyses): {WITHOUT_DARK}."
        )


def _refuse_dead_fraction(name: str, cid: ChannelId | None, point_ps: float, rate_cps: float) -> None:
    """An identified dead time is refused as a bound is, above `MAX_DEAD_FRACTION` at the recorded
    rate (council round 6, R1: a falsely identified 20 us dark dead time gave V = 0.55 +- 6.2 with
    exit 0, labelled only)."""
    worst = rate_cps * point_ps * 1e-12
    if worst > MAX_DEAD_FRACTION:
        raise GuardrailRefusal(
            f"the dead time of {_who(name, cid)} is {_duration(point_ps)}; at the recorded rate, rate x dead time = "
            f"{worst:.3g} (limit {MAX_DEAD_FRACTION:g}): there the pair term and the unmodelled "
            f"paralyzable stage are not trusted, and a dark dead time this long is often a false identification. "
            f"Correcting with a series dead-time calibration is not yet supported (see Known limits); "
            f"{WITHOUT_DARK}."
        )


def _pair_capture(
    calibrations: CalibrationSet | None,
    a: ChannelId,
    b: ChannelId,
    window_ps: int | None,
    declared_delay_ps: float,
    n: int,
    seed: int,
):
    """(draws, point, jitter bound matters) of the grid capture of pairs on channels a and b at
    zero offset, or None."""
    if calibrations is None or window_ps is None:
        return None
    hits = [m for m in calibrations.members if m.kind == DELAY_KIND and set(m.channels) == {a, b}]
    if not hits:
        return None
    cal = hits[0].calibration
    sign = 1.0 if tuple(cal.channels) == (a, b) else -1.0
    delay, jitter, d0, j0, bounded = _delay_jitter(cal, calibrations.draws(n, seed), n)
    draws = grid_capture(sign * delay - declared_delay_ps, jitter, window_ps)
    point = float(grid_capture(sign * d0 - declared_delay_ps, j0, window_ps))
    wide = bounded and _bound_matters(point, float(grid_capture(sign * d0 - declared_delay_ps, 0.0, window_ps)))
    return draws, point, wide


def _chsh_point(counts, signs, taus: Mapping[str, float], caps: Mapping[tuple[str, str], float]) -> float:
    """Point S from point counts, dead times (s) per channel name and captures per outcome pair."""
    s = 0.0
    for key, sign in signs.items():
        num = den = 0.0
        for (x, y), o_sign in OUTCOME_SIGNS.items():
            run = counts[key][(x, y)]
            v = float(_signal(run, run.coincidences, run.offpeak_sum, taus["A" + x], taus["B" + y], caps[(x, y)]))
            num += o_sign * v
            den += v
        s += sign * num / den
    return s


def chsh_bound_range(
    counts: Mapping[tuple[str, str], Mapping[tuple[str, str], RunCounts]],
    channels: Mapping[str, ChannelId],
    signs: Mapping[tuple[str, str], int],
    *,
    calibrations: CalibrationSet | None,
    window_ps: int | None = None,
    delays_ps: Mapping[str, float] | None = None,
) -> Mapping[str, float] | None:
    """Decision D3 for CHSH: the point S with every bounded dead time at zero and at its bound
    (the others at their point values), or None when no dead time is a bound."""
    dead = _channel_dead_times(calibrations, channels, 1, 0)
    if dead is None:
        return None
    bounded = {name for name, (_, _, hi) in dead.items() if hi is not None}
    if not bounded:
        return None
    delays = dict(delays_ps or {})
    caps = {}
    for x, y in OUTCOME_SIGNS:
        cap = _pair_capture(
            calibrations, channels["A" + x], channels["B" + y], window_ps,
            delays.get("B" + y, 0.0) - delays.get("A" + x, 0.0), 1, 0,
        )
        caps[(x, y)] = cap[1] if cap is not None else 1.0

    def at(scale: float) -> float:
        taus = {name: (scale * point if name in bounded else point) for name, (_, point, _) in dead.items()}
        return _chsh_point(counts, signs, taus, caps)

    return {"at_zero": at(0.0), "at_bound": at(1.0)}


def corrected_chsh(
    counts: Mapping[tuple[str, str], Mapping[tuple[str, str], RunCounts]],
    channels: Mapping[str, ChannelId],
    signs: Mapping[tuple[str, str], int],
    *,
    calibrations: CalibrationSet | None,
    policy: Policy | None = None,
    seed: int = 0,
    n_draws: int = N_MC,
    window_ps: int | None = None,
    delays_ps: Mapping[str, float] | None = None,
) -> Corrected:
    """Corrected CHSH S from the forward model, per setting pair and outcome combination.

    `counts[setting][(x, y)]` holds the coincidences of channels A{x} and B{y} (with their
    singles as `singles_a`, `singles_b`), `channels` maps "A+", "A-", "B+", "B-" to their
    ChannelIds. Each combination's signal is s = (C - A) / (T L capture): accidentals from its
    own off-peak bins, L the live fraction of that channel pair from the propagated dead times
    (its conditional term is the dead-time selection: detectors that fire together are dead
    together), and the grid capture from a delay-jitter artifact of the pair. E is the normalised
    signed sum of the four signals; S the signed sum of the E. Efficiencies are not corrected:
    E normalises them out only if they are equal (fair sampling).

    Corrected S does not correct detector efficiencies: an output imbalance biases it as it does
    raw S (0.6/0.4 on both sides: corrected S 2.033 for a 1.980 source, 13 sigma high). The
    "output imbalance" label applies to both values.
    """
    policy = policy or Policy.default()
    delays = dict(delays_ps or {})
    corrections: list[str] = []
    labels: list[str] = [FORWARD_MODEL] + _policy_labels(policy)
    dead = _channel_dead_times(calibrations, channels, n_draws, seed)
    if dead is not None:
        policy.check_correction("chsh", "dead_time")
        corrections.append("dead_time")
        bounded = {name: hi for name, (_, _, hi) in dead.items() if hi is not None}
        if bounded:
            labels.append(DEAD_TIME_BOUND)
        for name, (_, point_s, hi) in dead.items():
            side = "singles_a" if name.startswith("A") else "singles_b"
            rate = max(getattr(run, side) / run.duration_s for setting in counts.values() for (x, y), run in setting.items() if (("A" + x) if side == "singles_a" else ("B" + y)) == name)
            if hi is not None:
                _refuse_loose_bound(name, channels[name], hi, rate)
            else:
                _refuse_dead_fraction(name, channels[name], point_s * 1e12, rate)
    policy.check_correction("chsh", "accidental_subtraction")
    corrections.append("accidental_subtraction")
    captures = {
        (x, y): _pair_capture(
            calibrations, channels["A" + x], channels["B" + y], window_ps,
            delays.get("B" + y, 0.0) - delays.get("A" + x, 0.0), n_draws, seed,
        )
        for (x, y) in OUTCOME_SIGNS
    }
    have = [k for k, v in captures.items() if v is not None]
    if have:
        if len(have) != len(OUTCOME_SIGNS):
            missing = [f"A{x}/B{y}" for (x, y), v in captures.items() if v is None]
            raise ValueError(
                f"delay-jitter artifacts are given but not for {missing}: every channel pair needs one; take a "
                "heralding run for each missing pair, or move the delay-jitter artifacts (.json and .posterior.npy) "
                "out of the folder to correct without window capture"
            )
        policy.check_correction("chsh", "window_capture")
        corrections.append("window_capture")
        if any(v[2] for v in captures.values()):
            labels.append(JITTER_BOUND)

    rng = np.random.default_rng(np.random.SeedSequence([int(seed) & 0xFFFFFFFF, 0x434853]))
    zeros = (np.zeros(n_draws), 0.0, None)

    def tau(name):
        return dead[name] if dead is not None else zeros

    joint = np.zeros(n_draws)
    stat = np.zeros(n_draws)
    calib = np.zeros(n_draws)
    point = 0.0
    for key, sign in signs.items():
        if key not in counts:
            raise ValueError(f"missing counts for setting pair {key}")
        num = {"j": 0.0, "s": 0.0, "c": 0.0, "p": 0.0}
        den = {"j": 0.0, "s": 0.0, "c": 0.0, "p": 0.0}
        for (x, y), o_sign in OUTCOME_SIGNS.items():
            run = counts[key][(x, y)]
            (ta, pa, _), (tb, pb, _) = tau("A" + x), tau("B" + y)
            cap = captures[(x, y)] or (1.0, 1.0)
            c = rng.gamma(run.coincidences + 0.5, size=n_draws)
            off = rng.gamma(run.offpeak_sum + 0.5, size=n_draws)
            sig = {
                "j": _signal(run, c, off, ta, tb, cap[0]),
                "s": _signal(run, c, off, pa, pb, cap[1]),
                "c": _signal(run, run.coincidences, run.offpeak_sum, ta, tb, cap[0]),
                "p": _signal(run, run.coincidences, run.offpeak_sum, pa, pb, cap[1]),
            }
            for k, v in sig.items():
                num[k] = num[k] + o_sign * v
                den[k] = den[k] + v
        e = {k: num[k] / den[k] for k in num}
        joint = joint + sign * e["j"]
        stat = stat + sign * e["s"]
        calib = calib + sign * e["c"]
        point += sign * float(e["p"])
    finite = joint[np.isfinite(joint)]
    if finite.size < n_draws // 2:
        raise ValueError("accidentals leave no signal in most draws")
    lo, hi = (float(q) for q in np.quantile(finite, [0.025, 0.975]))
    calibrated = dead is not None or bool(have)
    est = Estimate(
        point,
        float(np.std(finite)),
        lo,
        hi,
        stderr_statistical=float(np.std(stat[np.isfinite(stat)])),
        stderr_calibration=float(np.std(calib)) if calibrated else 0.0,
    )
    pairs = [(channels["A" + x], channels["B" + y]) for (x, y) in OUTCOME_SIGNS] if have else []
    ids, used_labels = _used(calibrations, list(channels.values()) if dead is not None else (), pairs)
    labels += [lab for lab in _calibration_labels(calibrations, ids) + list(used_labels) if lab not in labels]
    labels.insert(1, chsh_model_label(corrections))
    return Corrected(est, tuple(corrections), tuple(labels), ids, _mc_counts(seed, n_draws, calibrations, ids))


__all__ = [
    "ACCIDENTALS_MODELLED",
    "AFTERPULSING_UNMODELLED",
    "CORRECTED_ABOVE_ONE",
    "DEAD_TIME_RANGE_WIDE",
    "JITTER_BOUND",
    "chsh_bound_range",
    "check_dead_time_gaps",
    "chsh_model_label",
    "OUTCOME_SIGNS",
    "corrected_chsh",
    "CAPTURE_UNMODELLED",
    "DEAD_FRACTION_HIGH",
    "DEAD_TIME_ABOVE_GAPS",
    "DEAD_TIME_ABOVE_GAPS_BOUND",
    "DEAD_TIME_UNMODELLED",
    "DEAD_TIME_BELOW_DARK",
    "DEAD_TIME_BOUND",
    "FORWARD_MODEL",
    "SECONDARY_V",
    "WITHOUT_DARK",
    "N_MC",
    "Corrected",
    "RunCounts",
    "bound_range",
    "corrected_visibility",
    "grid_capture",
    "live_fraction",
]

