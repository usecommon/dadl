"""CHSH from fixed-grid coincidences.

Headline S uses raw counts only. Every coincidence-normalized S is labelled
"fair sampling assumed", because normalizing by coincidences post-selects on
detection.

Phase 1 limits, stated so nobody over-reads the result:
- the window's `source` is declared by the caller, not verified; signed
  analysis plans (phase 3) will bind it to a plan fixed before acquisition;
- the tags are trusted as given (no dataset ID yet);
- the error is statistical only and Gaussian.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np

from dadl import _checks
from dadl.coincidences import grid_coincidences
from dadl.errors import InsufficientCounts
from dadl.guardrails import DEAD_TIME_SELECTION, EFFICIENCY_IMBALANCE, INTERRUPTED, Policy, threshold_labels
from dadl.stationarity import interrupted
from dadl.result import Estimate, LabelledValue

FAIR_SAMPLING = "fair sampling and fair coincidence assumed (post-selected on coincidences)"
#: A fixed grid closes the event-centred window loophole but not the fair-coincidence one: a
#: local model that records every photon and shifts timing with the setting reached S = 2.8
#: (expert review). Static settings leave locality and freedom of choice open.
STATIC_SETTINGS = "static settings: locality, freedom-of-choice and memory loopholes open; not a loophole-free test"

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


def chsh(
    runs: dict[tuple[str, str], dict[str, np.ndarray]],
    window: Window,
    policy: Policy | None = None,
    *,
    settings_deg: dict[str, float] | None = None,
) -> LabelledValue:
    """Raw CHSH S from one run per setting pair, with a fixed window from the analysis plan.

    The signs are fixed (`SIGNS`) and recorded with the result, with the declared analyser angles
    when given, so the sign assignment is part of the record rather than a later choice.

    `runs` maps each setting pair ("a", "b"), ("a", "b'"), ("a'", "b"),
    ("a'", "b'") to a dict of channel name -> int64 ps tags for A+, A-, B+, B-.
    """
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
    labels = (FAIR_SAMPLING, STATIC_SETTINGS) + policy.labels_for("chsh") + threshold_labels("chsh", raw)
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
    if settings_deg is not None:
        counts["settings_deg"] = {k: float(settings_deg[k]) for k in ("a", "a'", "b", "b'")}
    return LabelledValue(
        analysis="chsh",
        raw=raw,
        labels=labels,
        window=window.as_record(),
        counts=counts,
        **policy.provenance(),
    )
