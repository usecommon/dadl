"""CHSH from fixed-grid coincidences.

Headline S uses raw counts only. Every coincidence-normalized S is labelled
"fair sampling assumed", because normalizing by coincidences post-selects on
detection.

Phase 1 limits, stated so nobody over-reads the result:
- the window's `source` is declared by the caller, not verified; signed
  analysis plans (phase 3) will bind it to a plan fixed before acquisition;
- the tags are trusted as given (no dataset ID yet);
- the raw error is statistical only and Gaussian (the corrected value's includes calibration).

A delay error on some channels changes the window capture of some outcome pairs and not others
and moves S. Each outcome pair's peak is fitted over the four settings (`capture_balance`); S
reweighted by those captures is in the counts, with its standard deviation from the peak fits
(`capture_sd_data`). `CAPTURE_IMBALANCE` is added when the reweighted S differs from raw S by more
than half a standard error and by more than `CAPTURE_Z` = 3.29 of those standard deviations (two-sided 0.001),
or when a capture is below `MIN_CAPTURE` with its peak off the declared delay by more than half
its jitter. A capture below `MIN_CAPTURE` with the peak on the declared delay is a window narrow
against the jitter, not a delay error: `CAPTURE_LOW`.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np

from dadl import _checks
from dadl.coincidences import grid_coincidences
from dadl.errors import InsufficientCounts
from dadl.guardrails import (
    DEAD_TIME_SELECTION,
    EFFICIENCY_IMBALANCE,
    INTERRUPTED,
    Policy,
    corrected_threshold_labels,
    t3_input,
    threshold_labels,
)
from dadl.stationarity import interrupted
from dadl.result import Estimate, LabelledValue

if TYPE_CHECKING:
    from dadl.data import ChannelId
    from dadl.propagate import CalibrationSet

FAIR_SAMPLING = "fair sampling and fair coincidence assumed (post-selected on coincidences)"
#: A fixed grid closes the event-centred window loophole but not the fair-coincidence one: a
#: local model that records every photon and shifts timing with the setting reached S = 2.8
#: (expert review). Static settings leave locality and freedom of choice open.
STATIC_SETTINGS = "static settings: locality, freedom-of-choice and memory loopholes open; not a loophole-free test"

#: A delay error on some channels changes the window capture of some outcome pairs and not others,
#: which moves S (A+ and B+ 400 ps late, declared 0: S = 2.20 for a 1.98 source; council round 2).
CAPTURE_IMBALANCE = (
    "unequal window capture across outcome pairs: the declared delays miss the pair peaks (offsets in counts), "
    "shifting S by about {dS:+.3f}; declare the delays from a pair calibration"
)
CAPTURE_LOW = (
    "window narrow against the pair jitter: capture {c:.2f} on {pair} with the peak on the declared delay (offsets "
    "in counts); S is post-selected more strongly on the timing"
)
#: On the corrected S when the raw value carries `DEAD_TIME_SELECTION` and no dark artifact is used.
DEAD_TIME_UNCORRECTED = "dead time not modelled: the dead-time selection bias of the raw value stays in the corrected value"
#: Smallest window capture of an outcome pair before it is labelled.
MIN_CAPTURE = 0.5
#: Two-sided 0.001 normal quantile: the capture-balance label's noise threshold.
CAPTURE_Z = 3.29
#: Two-sided 2e-4 normal quantile: the delay-drift refusal's noise threshold, its share of the 1e-3
#: per-analysis refusal budget (council round 8; it was 3.29, about 5e-4 per clean analysis).
DRIFT_Z = 3.72
#: The delay drift accepted under `misfit_forward_model` (council round 8, St7: it shared the
#: predictive check's label, so the two could not be counted apart).
DELAY_DRIFT_ACCEPTED = (
    "delay drift accepted: the calibrated delays do not match the analysed pair peaks (policy modified)"
)
#: Capture draws of the delay-jitter artifacts behind the delay-drift check's noise.
DRIFT_DRAWS = 400

#: Below this many coincidences per setting the Gaussian error is not trusted.
MIN_COUNTS_PER_SETTING = 100

#: Channel names a CHSH run must contain.
CHANNELS: tuple[str, ...] = ("A+", "A-", "B+", "B-")

#: Setting pairs and their signs in S = E(a,b) - E(a,b') + E(a',b) + E(a',b').
SIGNS: dict[tuple[str, str], int] = {("a", "b"): 1, ("a", "b'"): -1, ("a'", "b"): 1, ("a'", "b'"): 1}


@dataclass(frozen=True)
class Window:
    """A fixed coincidence window. `source` must be "plan" for Bell-type analyses.

    `source` has no default: saying "plan" is the user's declaration, not the library's (API
    review, cycle 7). `delays_ps` maps channel name to its declared delay in integer ps and is
    read-only after construction.
    """

    width_ps: int
    source: str
    delays_ps: Mapping[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "width_ps", _checks.integer(self.width_ps, "window width_ps", minimum=1))
        if not isinstance(self.source, str):
            raise TypeError("window source must be a string")
        delays = {str(ch): _checks.integer(d, f"delay for {ch}") for ch, d in dict(self.delays_ps).items()}
        object.__setattr__(self, "delays_ps", MappingProxyType(delays))

    def as_record(self) -> dict[str, object]:
        return {
            "width_ps": self.width_ps,
            "delays_ps": dict(sorted(self.delays_ps.items())),
            "source": self.source,
            # Slot k holds tags with k W <= t - delay < (k + 1) W, t in the tag file's own clock.
            "grid_origin_ps": 0,
        }


#: Gap quantiles behind the dead-time estimate; the largest resulting lower bound is used.
DEAD_TIME_QUANTILES = (0.05, 0.15, 0.30)


def _dead_time_estimate(gaps: np.ndarray) -> float:
    """Dead time from the gap distribution: the q-quantile of tau + Exp(R) gaps is
    tau - ln(1 - q) / R, so the exponential part is subtracted (the bare low quantile read about
    0.01 / R at low counts and labelled every short run; audit, cycle 6)."""
    if gaps.size < 20:
        return 0.0
    # The largest of the lower bounds from several quantiles: a cluster of echoes or double triggers
    # below the dead time sets the low quantiles but not the higher ones (second red-team review,
    # cycle 9: 5.5% ringing echoes removed the label at the 5% quantile).
    best = 0.0
    for p in DEAD_TIME_QUANTILES:
        q = float(np.quantile(gaps, p))
        excess = float(gaps.mean()) - q
        if excess <= 0:
            continue
        # For gaps tau + Exp(R): mean - Q_p = (1 + ln(1 - p)) / R.
        rate = (1.0 + math.log1p(-p)) / excess
        estimate = q + math.log1p(-p) / rate
        # A lower bound, 3 standard errors of the quantile below: the label then needs the dead
        # time to be established, not a fluctuation of a short run's quantile.
        # Quantile noise plus the mean gap's noise, which enters through the rate.
        lever = math.log1p(-p) / (1.0 + math.log1p(-p))
        se = math.sqrt(p * (1 - p) / gaps.size / ((1 - p) * rate) ** 2 + lever**2 * float(gaps.var()) / gaps.size)
        best = max(best, estimate - 4.0 * se)  # 4 sigma: the maximum is taken over many bounds
    return best


def dead_time_selection(runs: dict, coincidences: dict | None = None) -> float:
    """Largest coincidence rate x dead time over the settings, each channel's dead time estimated
    from its gap distribution and the longest taken (the bare minimum let one duplicated tag remove
    the label: red-team review, cycle 5; the minimum over channels missed mismatched dead times).
    The shift in |E| is about 0.05-0.08 eta^2 R tau (audit, cycle 6).

    Teaching-lab review, cycle 4: with a 1 us dead time at 3e6 pairs/s a V = 0.7 source (S = 1.98)
    read S = 2.14 +- 0.004; the shift grows with eta^2 R tau (about 0.05-0.08 eta^2 R tau in |E|).
    """
    worst = 0.0
    for key, tags in runs.items():
        arrays = [np.asarray(tags[ch]) for ch in CHANNELS if np.asarray(tags[ch]).size > 1]
        if not arrays:
            continue
        # Each channel's dead time as a low gap quantile (robust to a few duplicates or double
        # triggers), and the longest of them: the bias comes from the slow channel, and a mismatch
        # between channels hid it when the minimum was used (point-process review, cycle 5).
        tau = max(_dead_time_estimate(np.diff(a)) for a in arrays)
        span = max(float(a[-1]) for a in arrays) - min(float(a[0]) for a in arrays)
        if span <= 0:
            continue
        n_c = (coincidences or {}).get(key)
        if n_c is None:
            continue
        worst = max(worst, n_c / span * tau)
    return worst


def check_no_duplicates(runs: dict, channels: tuple[str, ...]) -> None:
    """Refuse duplicate timestamps on one channel: a detector cannot click twice at once, and
    duplicates can hide the dead-time-selection label (in the library, not only the CLI: audit,
    cycle 6)."""
    for key, tags in runs.items():
        for ch in channels:
            arr = np.asarray(tags.get(ch, ()))
            if arr.size > 1 and np.any(np.diff(np.sort(arr)) == 0):
                raise ValueError(f"{ch} in acquisition {key}: duplicate timestamps on one channel")


def check_delay_channels(window: Window, channels: tuple[str, ...]) -> None:
    """Refuse delays declared for channels the analysis does not have (a misspelt name was ignored
    and its channel ran at delay 0: API review, cycle 7)."""
    unknown = sorted(set(window.delays_ps) - set(channels))
    if unknown:
        raise ValueError(f"delays declared for unknown channel(s) {unknown}; the channels are {list(channels)}")


def check_reference_channel(window: Window, channels: tuple[str, ...]) -> None:
    """Refuse delays with no channel at 0: a common offset only moves the fixed grid's origin.

    Bell review, cycle 2: adding one offset to every delay leaves the relative delays (and the
    delay check) unchanged but moves pairs across slot edges; picking the best of ten such
    offsets doubled the rate of S > 2 + 1.96 sigma at a true S of 2.00.
    """
    delays = [window.delays_ps.get(ch, 0) for ch in channels]
    if min(delays) != 0:
        raise ValueError(
            f"declared delays {dict(zip(channels, delays))} have no reference channel at 0 ps; subtract "
            f"{min(delays)} ps from every delay (a common offset only moves the coincidence grid's origin)"
        )


def correlator(tags: dict[str, np.ndarray], window: Window) -> tuple[float, float, int]:
    """E = (N++ + N-- - N+- - N-+) / N and its variance (1 - E^2) / N.

    Double clicks: a slot with both outputs of one side and a click on the other counts in two
    outcome combinations, so N can exceed the coincident slots (0.06% of slots at 10^6 pairs/s
    on the example). This pulls E toward 0 (conservative for CHSH) and is not the
    fixed-random-assignment rule a CH/Eberhard analysis needs; that rule is chosen with them.
    """
    missing = [ch for ch in CHANNELS if ch not in tags]
    if missing:
        raise ValueError(f"run is missing channel(s) {missing}; a CHSH run needs {list(CHANNELS)}")
    d = window.delays_ps

    def n(a: str, b: str) -> int:
        return grid_coincidences(tags[a], tags[b], window.width_ps, d.get(a, 0), d.get(b, 0))

    npp, nmm, npm, nmp = n("A+", "B+"), n("A-", "B-"), n("A+", "B-"), n("A-", "B+")
    total = npp + nmm + npm + nmp
    if total < MIN_COUNTS_PER_SETTING:
        raise InsufficientCounts(
            f"{total} coincidences in a setting; need at least {MIN_COUNTS_PER_SETTING}. "
            "Check that timestamps are int64 picoseconds, that the declared delays align the channels, and that "
            "the window is not narrow against the timing jitter."
        )
    e = (npp + nmm - npm - nmp) / total
    return e, (1.0 - e * e) / total, total


def setting_counts(tags: dict[str, np.ndarray], window: Window) -> dict[str, float]:
    """Singles, coincidences, coincidences per singles on each arm, and each side's output imbalance.

    The detection-loophole thresholds need the condition on both arms, so the smaller ratio is the
    one to compare (senior-physicist review, cycle 10: the larger one was reported, 0.85 for an arm
    at 0.28).
    """
    _, _, total = correlator(tags, window)
    singles_a = int(tags["A+"].size + tags["A-"].size)
    singles_b = int(tags["B+"].size + tags["B-"].size)
    per_a = total / max(1.0, singles_a)
    per_b = total / max(1.0, singles_b)
    return {
        "coincidences": float(total),
        "singles_a": float(singles_a),
        "singles_b": float(singles_b),
        "coincidences_per_singles_a": per_a,
        "coincidences_per_singles_b": per_b,
        "coincidences_per_singles": min(per_a, per_b),
        "imbalance_a": (tags["A+"].size - tags["A-"].size) / max(1.0, singles_a),
        "imbalance_b": (tags["B+"].size - tags["B-"].size) / max(1.0, singles_b),
    }


def _pair_name(x: str, y: str) -> str:
    return f"A{x}/B{y}"


def _setting_outcome_counts(runs, window: Window) -> dict:
    """Grid coincidences per setting pair and outcome combination."""
    from dadl.forward import OUTCOME_SIGNS

    d = window.delays_ps
    return {
        key: {
            (x, y): grid_coincidences(runs[key]["A" + x], runs[key]["B" + y], window.width_ps, d.get("A" + x, 0), d.get("B" + y, 0))
            for x, y in OUTCOME_SIGNS
        }
        for key in SIGNS
    }


def reweighted_s(counts, captures: Mapping[tuple[str, str], float]):
    """S with each outcome combination's count (`_setting_outcome_counts`) divided by its window
    capture. Captures may be arrays of draws; then so is S."""
    from dadl.forward import OUTCOME_SIGNS

    s = 0.0
    for key, sign in SIGNS.items():
        n = {k: counts[key][k] / captures[k] for k in OUTCOME_SIGNS}
        s = s + sign * sum(o * n[k] for k, o in OUTCOME_SIGNS.items()) / sum(n.values())
    return s


def capture_balance(runs, window: Window, counts=None, lattice_ps: int | None = None):
    """S reweighted by each outcome pair's window capture, per pair (offset, jitter, capture,
    offset se, jitter se), and the reweighted S's standard deviation from the peak fits.

    Each outcome pair's peak is fitted once (`forward_check.pooled_peak`) to its histogram summed
    over the four settings, and its capture is `forward.grid_capture(offset, jitter, W)`. The
    standard deviation is the delta method's: each pair's offset, then its jitter, moved by its
    standard error, S recomputed, the squared changes summed. A pair below
    `forward_check.MIN_COUNTS` is skipped; then the reweighted S and its deviation are None.
    `lattice_ps` is the runs' lattice of pair delays (`dadl.data.runs_lattice`; None or 0: the
    tags' own step)."""
    from dadl.forward import OUTCOME_SIGNS, grid_capture
    from dadl.forward_check import PairData, pooled_peak

    d = window.delays_ps
    out: dict[tuple[str, str], tuple[float, float, float, float, float]] = {}
    for x, y in OUTCOME_SIGNS:
        a, b = "A" + x, "B" + y
        declared = float(d.get(b, 0) - d.get(a, 0))
        fit = pooled_peak([
            PairData(f"{k[0]},{k[1]} {a}/{b}", runs[k][a], runs[k][b], declared, window.width_ps, lattice_ps=lattice_ps or None)
            for k in SIGNS
        ])
        if fit is not None:
            cap = float(grid_capture(fit.offset_ps, fit.sigma_ps, window.width_ps))
            out[(x, y)] = (fit.offset_ps, fit.sigma_ps, cap, fit.offset_se_ps, fit.sigma_se_ps)
    if len(out) != len(OUTCOME_SIGNS):
        return None, out, None
    counts = counts if counts is not None else _setting_outcome_counts(runs, window)
    caps = {k: v[2] for k, v in out.items()}
    base = float(reweighted_s(counts, caps))
    var = 0.0
    for k, (offset, sigma, _, se_offset, se_sigma) in out.items():
        for do, ds in ((se_offset, 0.0), (0.0, se_sigma)):
            moved = dict(caps)
            moved[k] = float(grid_capture(offset + do, sigma + ds, window.width_ps))
            var += (float(reweighted_s(counts, moved)) - base) ** 2
    return base, out, math.sqrt(var)


def capture_labels(raw: Estimate, s_balanced, balance, sd_data) -> tuple[str, ...]:
    """`CAPTURE_IMBALANCE`, `CAPTURE_LOW` or nothing (module docstring)."""
    shift = raw.value - s_balanced if s_balanced is not None else 0.0
    moved = s_balanced is not None and abs(shift) > max(0.5 * raw.stderr, CAPTURE_Z * (sd_data or 0.0))
    low = {k: v for k, v in balance.items() if v[2] < MIN_CAPTURE}
    if moved or any(abs(v[0]) > v[1] / 2 for v in low.values()):
        return (CAPTURE_IMBALANCE.format(dS=shift),)
    if low:
        k = min(low, key=lambda k: low[k][2])
        return (CAPTURE_LOW.format(c=low[k][2], pair=_pair_name(*k)),)
    return ()


def chsh(
    runs: dict[tuple[str, str], dict[str, np.ndarray]],
    window: Window,
    policy: Policy | None = None,
    *,
    settings_deg: dict[str, float] | None = None,
    calibrations: "CalibrationSet | None" = None,
    channels: "dict[str, ChannelId] | None" = None,
    seed: int = 0,
    source_type: str | None = None,
) -> LabelledValue:
    """Raw CHSH S from one run per setting pair, with a fixed window from the analysis plan.

    The signs are fixed (`SIGNS`) and recorded with the result, with the declared analyser angles
    when given, so the sign assignment is part of the record rather than a later choice.

    `runs` maps each setting pair ("a", "b"), ("a", "b'"), ("a'", "b"),
    ("a'", "b'") to a dict of channel name -> int64 ps tags for A+, A-, B+, B-.

    With `calibrations` (a `dadl.propagate.CalibrationSet`) and `channels` (each channel name's
    ChannelId), a corrected S from the forward model (`dadl.forward.corrected_chsh`) is added as a
    labelled secondary value; the headline stays raw. The default policy allows no correction
    for CHSH, so each one needs `Policy.allow("chsh", correction, reason)`. On PTU T3 input the
    corrected S is refused unless `source_type` declares CW light ("cw"), and the result is
    labelled (`dadl.guardrails.t3_input`).
    """
    t3_labels, _ = t3_input(runs, source_type, calibrations is not None)
    from dadl.data import runs_lattice

    lattice = runs_lattice(runs)
    policy = policy or Policy.default()
    policy.check_window("chsh", window.source)
    check_delay_channels(window, CHANNELS)
    check_reference_channel(window, CHANNELS)
    check_no_duplicates(runs, CHANNELS)
    # Counting does not need sorted tags, but the dead-time and interruption guards do (unsorted
    # input dropped both labels: API review, cycle 7).
    runs = {k: {ch: np.sort(np.asarray(v)) for ch, v in r.items()} for k, r in runs.items()}
    s = 0.0
    var = 0.0
    for key, sign in SIGNS.items():
        if key not in runs:
            raise ValueError(f"missing run for setting pair {key}; need {list(SIGNS)}")
        e, v, _ = correlator(runs[key], window)
        s += sign * e
        var += v
    raw = Estimate(s, math.sqrt(var))
    labels = (FAIR_SAMPLING, STATIC_SETTINGS) + policy.labels_for("chsh") + threshold_labels("chsh", raw) + t3_labels
    if interrupted(*([runs[k][ch] for ch in CHANNELS] for k in SIGNS)):
        labels += (INTERRUPTED,)
    counts: dict[str, object] = {f"{k[0]},{k[1]}": setting_counts(runs[k], window) for k in SIGNS}
    counts["signs"] = {f"{k[0]},{k[1]}": s for k, s in SIGNS.items()}
    if dead_time_selection(runs, {k: counts[f"{k[0]},{k[1]}"]["coincidences"] for k in SIGNS}) > 0.01:
        labels += (DEAD_TIME_SELECTION,)
    # Unequal output efficiencies within a side shift raw E by about m_a m_b (1 - E^2) and S by up
    # to 2 m_a m_b, upward when the two sides' imbalances share a sign (0.6/0.4 on both sides: S =
    # 2.033 from a source at 1.98). m is read from the singles, uniform for a maximally entangled state.
    m_a = float(np.mean([counts[f"{k[0]},{k[1]}"]["imbalance_a"] for k in SIGNS]))
    m_b = float(np.mean([counts[f"{k[0]},{k[1]}"]["imbalance_b"] for k in SIGNS]))
    if 2.0 * abs(m_a * m_b) > 0.5 * raw.stderr:
        labels += (EFFICIENCY_IMBALANCE,)
    outcome_counts = _setting_outcome_counts(runs, window)
    s_balanced, balance, sd_data = capture_balance(runs, window, outcome_counts, lattice)
    counts["capture_offsets_ps"] = {_pair_name(*k): v[0] for k, v in balance.items()}
    counts["capture_reweighted_s"] = s_balanced
    counts["capture_sd_data"] = sd_data
    labels += capture_labels(raw, s_balanced, balance, sd_data)
    if settings_deg is not None:
        counts["settings_deg"] = {k: float(settings_deg[k]) for k in ("a", "a'", "b", "b'")}
    corrected = None
    extra: dict[str, object] = {}
    corrected_labels: tuple[str, ...] = ()
    if calibrations is not None:
        corrected, extra, more, more_counts = _corrected(
            runs, window, policy, calibrations, channels, seed, raw, (outcome_counts, s_balanced, balance, sd_data),
            lattice,
        )
        counts.update(more_counts)
        if DEAD_TIME_SELECTION in labels and "dead_time" not in extra["corrections"]:
            more += (DEAD_TIME_UNCORRECTED,)
        # Council round 6 (R2): the corrected value's labels are kept apart from the raw value's.
        corrected_labels = tuple(lab for lab in more if lab not in labels)
    return LabelledValue(
        analysis="chsh",
        raw=raw,
        corrected=corrected,
        labels=labels,
        corrected_labels=corrected_labels,
        window=window.as_record(),
        counts=counts,
        **extra,
        **policy.provenance(),
    )


def _corrected(
    runs, window: Window, policy: Policy, calibrations, channels, seed: int, raw: Estimate, capture, lattice=None
):
    """The forward-model S, its result fields, its labels (threshold crossing included) and the
    counts it adds. `capture` is (outcome counts, reweighted S, balance, its sd) from `chsh`;
    `lattice` the runs' lattice of pair delays (`dadl.data.runs_lattice`)."""
    from dadl.analysis._twopoint import offpeak_check, run_counts, span_ps
    from dadl.calibration.accidentals import shifted_coincidences, side_offsets
    from dadl.forward import DEAD_TIME_RANGE_WIDE, OUTCOME_SIGNS, check_dead_time_gaps, chsh_bound_range, corrected_chsh

    if channels is None or set(channels) != set(CHANNELS):
        raise ValueError(f"with calibrations, pass channels mapping each of {list(CHANNELS)} to its ChannelId")
    from dadl.propagate import check_detectors

    d = window.delays_ps
    offsets = side_offsets(window.width_ps)
    detector_labels = check_detectors(calibrations, channels.values(), analysis="chsh", policy=policy)
    check_p, check_labels = _forward_check(runs, window, policy, calibrations, channels, lattice)
    gap_labels, gap_counts = check_dead_time_gaps(calibrations, channels, [runs[k] for k in SIGNS], analysis="chsh", policy=policy)
    reduced = {}
    sides = {}
    for key in SIGNS:
        tags = runs[key]
        span = max(span_ps(tags, (a, b)) for a, b in (("A+", "B+"), ("A-", "B-")))
        per = {}
        for x, y in OUTCOME_SIGNS:
            a, b = "A" + x, "B" + y
            c = grid_coincidences(tags[a], tags[b], window.width_ps, d.get(a, 0), d.get(b, 0))
            side = shifted_coincidences(tags[a], tags[b], window.width_ps, d.get(a, 0), d.get(b, 0), offsets)
            sides[f"{key[0]},{key[1]} {a}/{b}"] = side
            per[(x, y)] = run_counts(tags, c, side, span, channels=(a, b))
        reduced[key] = per
    offpeak_p, offpeak_labels = offpeak_check("chsh", sides, offsets, policy)
    grid = {"window_ps": window.width_ps, "delays_ps": {k: float(v) for k, v in d.items()}}
    fwd = corrected_chsh(reduced, channels, SIGNS, calibrations=calibrations, policy=policy, seed=seed, **grid)
    counts: dict[str, object] = {**({"forward_check_p": check_p} if check_p is not None else {}), **gap_counts,
                                 "offpeak_p": offpeak_p, **fwd.counts}
    if lattice == 0:
        from dadl.forward_check import NO_LATTICE

        counts["forward_check_note"] = NO_LATTICE
    if calibrations.dataset_record() is not None:
        counts["dataset"] = calibrations.dataset_record()
    labels = fwd.labels + detector_labels + check_labels + offpeak_labels + gap_labels
    drift_labels, drift_counts = _delay_drift(window, policy, calibrations, channels, seed, capture, fwd.estimate)
    labels += drift_labels
    counts.update(drift_counts)
    ends: tuple[float, ...] = ()
    rng = chsh_bound_range(reduced, channels, SIGNS, calibrations=calibrations, **grid)
    if rng is not None:
        counts["corrected_at_zero_dead_time"] = rng["at_zero"]
        counts["corrected_at_bound_dead_time"] = rng["at_bound"]
        ends = (rng["at_zero"], rng["at_bound"])
        if abs(rng["at_bound"] - rng["at_zero"]) > fwd.estimate.stderr:
            labels += (DEAD_TIME_RANGE_WIDE,)
    labels += corrected_threshold_labels("chsh", raw, fwd.estimate, ends)
    labels = tuple(dict.fromkeys(labels))
    return fwd.estimate, {"corrections": fwd.corrections, "artifact_ids": fwd.artifact_ids}, labels, counts


def _delay_drift(window: Window, policy: Policy, calibrations, channels, seed: int, capture, corrected: Estimate):
    """Compare S reweighted by the analysed peaks' captures with S reweighted by the delay-jitter
    artifacts' point captures: (labels, counts). A drift since calibration changes the window
    capture, and the predictive check refits the delay, so it does not see it.

    The difference's noise has two parts: the artifacts' (the sd of S reweighted by each of
    `DRIFT_DRAWS` capture draws) and the analysed peak fits' (`capture_balance`). Refused
    (`Misspecified`) when the difference exceeds half the corrected value's statistical error and
    `DRIFT_Z` times that noise, unless the policy allows `misfit_forward_model`; then labelled
    `DELAY_DRIFT_ACCEPTED`.

    The drift check refuses only a drift that moves S by more than half its statistical error and
    more than 3.72 standard deviations of the check's own noise (two-sided 2e-4) (the delay artifacts' capture draws
    and the analysed peak fit). A smaller drift is not refused (100 ps at 350 ps jitter and a 2 ns
    window: 0.41 se, 0/20 refused); one within the artifacts' noise is covered by the calibration
    part of the interval. Longer delay calibrations make the check sharper.

    On clean data the statistic has a small positive mean (+0.0005, about 0.1 statistical se, from
    the peak fit at offsets of +-400 ps); it is well inside the threshold."""
    from dadl.errors import Misspecified
    from dadl.forward import OUTCOME_SIGNS, _pair_capture
    from dadl.forward_check import OVERRIDE, delay_start

    counts, s_balanced, balance, sd_data = capture
    if s_balanced is None:
        return (), {}
    d = window.delays_ps
    artifact, draws = {}, {}
    for x, y in OUTCOME_SIGNS:
        a, b = "A" + x, "B" + y
        cap = _pair_capture(
            calibrations, channels[a], channels[b], window.width_ps, float(d.get(b, 0) - d.get(a, 0)), DRIFT_DRAWS, seed
        )
        if cap is None:
            return (), {}
        draws[(x, y)], artifact[(x, y)] = cap[0], cap[1]
    s_artifact = float(reweighted_s(counts, artifact))
    sd = math.hypot(float(np.std(reweighted_s(counts, draws))), sd_data)
    diff = s_balanced - s_artifact
    z = diff / sd if sd > 0 else 0.0
    out = {"delay_drift_diff": diff, "delay_drift_z": z}
    if abs(diff) <= max(0.5 * corrected.stderr_statistical, DRIFT_Z * sd):
        return (), out
    if OVERRIDE in policy.rule_for("chsh").allowed_corrections:
        return (DELAY_DRIFT_ACCEPTED,), out
    worst = max(OUTCOME_SIGNS, key=lambda k: abs(artifact[k] - balance[k][2]))
    a, b = "A" + worst[0], "B" + worst[1]
    declared = float(d.get(b, 0) - d.get(a, 0))
    calibrated = delay_start(calibrations, channels[a], channels[b])[0]
    raise Misspecified(
        f"the calibrated delays do not match the analysed pair peaks ({a}/{b}: calibrated {calibrated:.0f}, analysed "
        f"{declared + balance[worst][0]:.0f} ps; S {s_artifact:.4f} with the calibrated captures, {s_balanced:.4f} "
        f"with the analysed, z = {z:.1f}): a drift since calibration changes the window capture. Recalibrate the "
        f"delays, or override with Policy.allow('chsh', {OVERRIDE!r}, reason)."
    )


def _forward_check(
    runs, window: Window, policy: Policy, calibrations, channels, lattice=None
) -> tuple[float | None, tuple[str, ...]]:
    """The forward model's predictive check over all 16 histograms (4 settings x 4 outcome pairs),
    or (None, ()) without delay-jitter artifacts: the shape check needs the calibrated jitter, as in
    HOM and fringe (council round 4; an accidentals-only corrected S was refused). Skipped (None,
    ()) when `lattice` is 0: T3 input with no integer lattice."""
    if lattice == 0:
        return None, ()
    from dadl.forward import OUTCOME_SIGNS
    from dadl.forward_check import PairData, delay_start, enforce, run_check

    d = window.delays_ps
    pairs = []
    for key in SIGNS:
        tags = runs[key]
        for x, y in OUTCOME_SIGNS:
            a, b = "A" + x, "B" + y
            declared = float(d.get(b, 0) - d.get(a, 0))
            d_cal, jitter = delay_start(calibrations, channels[a], channels[b])
            if jitter is None:
                return None, ()
            pairs.append(PairData(
                f"{key[0]},{key[1]} {a}/{b}", np.sort(tags[a]), np.sort(tags[b]), declared, window.width_ps,
                d_cal - declared, jitter, lattice_ps=lattice,
            ))
    result = run_check(pairs)
    return result.pvalue, enforce(result, "chsh", policy)
