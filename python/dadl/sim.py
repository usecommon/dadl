"""Time-tag simulator with known, injected detector imperfections.

Every parameter that the calibration engine will later estimate is set
explicitly here, so tests can compare estimates against the truth. Parameter
names match setup.yaml (unit-suffixed).

Units: time tags are int64 picoseconds; rates are per second.

Detection is self-exciting, as in docs/design.md ("Detection process"):
every recorded avalanche (photon, dark count or afterpulse) starts a
non-paralyzable dead time and, with probability `afterpulse_prob`, schedules
an afterpulse at the end of its dead time plus an exponential delay, so
afterpulses cascade. An afterpulse that falls inside a later avalanche's dead
time is lost.

Series dead time (opt-in, `dead_time_p_ps` > 0): a paralyzable stage comes
first. Every avalanche candidate (photon, dark count, afterpulse or crosstalk)
restarts it, and a candidate passes only when the previous candidate, passed or
not, is at least `dead_time_p_ps` earlier; what passes then meets the
non-paralyzable `dead_time_ps` (design, "Detection process": a paralyzable stage
followed by a non-paralyzable one, tau_np > tau_p).

Crosstalk (opt-in, `Crosstalk`): with probability `probability`, every
recorded avalanche on the source channel triggers a crosstalk avalanche
candidate on the target channel `delay_ps` later (Gaussian jitter
`jitter_sigma_ps`, never earlier than its cause). A nonzero delay models a Si
SPAD breakdown flash reaching the other detector after 2 n_g L / c; zero
models electrical or pixel crosstalk. The candidate is lost inside the
target's dead time; recorded, it is an avalanche like any other: it starts a
dead time, can afterpulse, and can trigger crosstalk in turn, so links in both
directions cascade. All channels then run in one event loop in time order.
Without crosstalk (the default) every channel is simulated on its own, as
before, and results for a given seed are unchanged.

Not simulated yet: twilight counts and a separate TDC dead time. Arrival
times are drawn as float64 ps, exact to the picosecond only up to 2^53 ps
(about 2.5 hours of simulated time per run).
"""

from __future__ import annotations

import heapq
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from dadl import _checks
from dadl.analysis.chsh import CHANNELS, SIGNS
from dadl.deadtime import apply_nonparalyzable, apply_paralyzable

PS_PER_S = 10**12

#: Analyzer angles in degrees for the standard CHSH settings on a Phi+ state.
CHSH_SETTINGS: dict[str, float] = {"a": 0.0, "a'": 45.0, "b": 22.5, "b'": 67.5}
CHSH_PAIRS: tuple[tuple[str, str], ...] = tuple(SIGNS)

__all__ = [
    "CHANNELS",
    "CHSH_PAIRS",
    "CHSH_SETTINGS",
    "Crosstalk",
    "DetectorParams",
    "PairSource",
    "detect",
    "detect_channels",
    "poisson_times",
    "simulate_chsh",
    "simulate_chsh_setting",
    "simulate_dark",
    "simulate_hom",
    "simulate_pairs",
]


@dataclass(frozen=True)
class DetectorParams:
    """Ground truth for one detector channel."""

    efficiency: float = 0.5  # system efficiency (detector x coupling)
    dark_rate_cps: float = 100.0  # counts per second
    dead_time_ps: int = 22_000  # non-paralyzable
    jitter_sigma_ps: float = 350.0  # Gaussian sigma
    afterpulse_prob: float = 0.0  # per recorded avalanche
    afterpulse_tau_ps: float = 50_000.0  # exponential delay after the dead time
    delay_ps: int = 0  # fixed channel delay
    dead_time_p_ps: int = 0  # paralyzable stage before the non-paralyzable one (series model)

    def __post_init__(self) -> None:
        _checks.finite(self.efficiency, "efficiency", minimum=0.0, maximum=1.0)
        _checks.finite(self.dark_rate_cps, "dark_rate_cps", minimum=0.0)
        object.__setattr__(self, "dead_time_ps", _checks.integer(self.dead_time_ps, "dead_time_ps", minimum=0))
        _checks.finite(self.jitter_sigma_ps, "jitter_sigma_ps", minimum=0.0)
        p = _checks.finite(self.afterpulse_prob, "afterpulse_prob", minimum=0.0)
        if p >= 1.0:
            raise ValueError("afterpulse_prob must be in [0, 1)")
        _checks.finite(self.afterpulse_tau_ps, "afterpulse_tau_ps", minimum=0.0)
        object.__setattr__(self, "delay_ps", _checks.integer(self.delay_ps, "delay_ps"))
        tau_p = _checks.integer(self.dead_time_p_ps, "dead_time_p_ps", minimum=0)
        if tau_p > self.dead_time_ps:
            raise ValueError("dead_time_p_ps must not exceed dead_time_ps (the series model has tau_np >= tau_p)")
        object.__setattr__(self, "dead_time_p_ps", tau_p)


@dataclass(frozen=True)
class PairSource:
    """A CW polarization-entangled SPDC source producing a Werner-like Phi+ state."""

    pair_rate_per_s: float = 1e5  # pairs per second at the crystal
    #: CHSH: E(alpha, beta) = V cos 2(alpha - beta). HOM: the dip depth (indistinguishability).
    visibility: float = 0.95
    #: HOM only: the dip is exp(-(delay / coherence_time)^2).
    coherence_time_ps: float = 1.0

    def __post_init__(self) -> None:
        _checks.finite(self.pair_rate_per_s, "pair_rate_per_s", minimum=0.0)
        _checks.finite(self.visibility, "visibility", minimum=0.0, maximum=1.0)
        c = _checks.finite(self.coherence_time_ps, "coherence_time_ps", minimum=0.0)
        if c <= 0:
            raise ValueError("coherence_time_ps must be positive")


@dataclass(frozen=True)
class Crosstalk:
    """Crosstalk from one channel to another (ground truth for one direction of one pair).

    Every recorded avalanche on `source` triggers, with probability
    `probability`, a crosstalk avalanche candidate on `target` at
    `delay_ps` + N(0, `jitter_sigma_ps`) after it, in recorded time (the
    offset of the peak it makes in the target - source cross-correlation),
    floored at 0: crosstalk never precedes its cause. For both directions
    give two links.
    """

    source: str
    target: str
    probability: float
    delay_ps: float = 0.0
    jitter_sigma_ps: float = 0.0

    def __post_init__(self) -> None:
        if self.source == self.target:
            raise ValueError("crosstalk needs two different channels (use afterpulse_prob for a channel's own)")
        p = _checks.finite(self.probability, "probability", minimum=0.0)
        if p >= 1.0:
            raise ValueError("crosstalk probability must be in [0, 1)")
        _checks.finite(self.delay_ps, "delay_ps", minimum=0.0)
        _checks.finite(self.jitter_sigma_ps, "jitter_sigma_ps", minimum=0.0)


def poisson_times(rate: float, duration_s: float, rng: np.random.Generator) -> np.ndarray:
    """Sorted float64 arrival times (ps) of a homogeneous Poisson process."""
    _checks.finite(rate, "rate", minimum=0.0)
    _checks.finite(duration_s, "duration_s", minimum=0.0)
    n = rng.poisson(rate * duration_s)
    return np.sort(rng.uniform(0.0, duration_s * PS_PER_S, n))


def detect(
    photon_times_ps: np.ndarray,
    det: DetectorParams,
    duration_s: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Turn photon arrival times into recorded int64 time tags for one detector."""
    arrivals = _primaries(photon_times_ps, det, duration_s, rng)
    if det.afterpulse_prob == 0:
        if det.dead_time_p_ps > 0:
            arrivals = apply_paralyzable(arrivals, det.dead_time_p_ps)
        return apply_nonparalyzable(arrivals, det.dead_time_ps)
    return _self_exciting(
        arrivals, det.dead_time_ps, det.afterpulse_prob, det.afterpulse_tau_ps, rng, det.dead_time_p_ps
    )


def _primaries(
    photon_times_ps: np.ndarray, det: DetectorParams, duration_s: float, rng: np.random.Generator
) -> np.ndarray:
    """Sorted int64 primary avalanche candidates: detected photons (jittered, delayed) and dark counts."""
    photons = np.asarray(photon_times_ps, dtype=np.float64)
    kept = photons[rng.random(photons.size) < det.efficiency]
    if det.jitter_sigma_ps > 0:
        kept = kept + rng.normal(0.0, det.jitter_sigma_ps, kept.size)
    kept = kept + det.delay_ps
    darks = poisson_times(det.dark_rate_cps, duration_s, rng)
    arrivals = np.sort(np.rint(np.concatenate([kept, darks]))).astype(np.int64)
    if det.dead_time_ps <= 0:
        # Even an ideal detector records one tag per picosecond: two arrivals rounded to the same
        # ps are one click (analyses refuse duplicate timestamps).
        arrivals = np.unique(arrivals)
    return arrivals


def detect_channels(
    photons: Mapping[str, np.ndarray],
    detectors: Mapping[str, DetectorParams],
    duration_s: float,
    rng: np.random.Generator,
    crosstalk: Sequence[Crosstalk] = (),
) -> dict[str, np.ndarray]:
    """Recorded tags of several detectors, coupled by `crosstalk` links (see the module docstring).

    Without crosstalk this is `detect` per channel, in the order of `photons`.
    """
    missing = [ch for ch in photons if ch not in detectors]
    if missing:
        raise ValueError(f"missing detector parameters for {missing}")
    unknown = sorted({x for link in crosstalk for x in (link.source, link.target)} - set(photons))
    if unknown:
        raise ValueError(f"crosstalk names channel(s) {unknown}, not among {list(photons)}")
    if not crosstalk:
        return {ch: detect(photons[ch], detectors[ch], duration_s, rng) for ch in photons}
    names = list(photons)
    primaries = [_primaries(photons[ch], detectors[ch], duration_s, rng) for ch in names]
    out = _coupled(names, primaries, [detectors[ch] for ch in names], crosstalk, rng)
    return dict(zip(names, out))


class _Draws:
    """Random numbers drawn in blocks, one at a time (the event loop needs scalars)."""

    def __init__(self, rng: np.random.Generator, block: int = 4096) -> None:
        self._rng, self._block = rng, block
        self._u = self._n = self._e = np.empty(0)
        self._iu = self._in = self._ie = 0

    def uniform(self) -> float:
        if self._iu == self._u.size:
            self._u, self._iu = self._rng.random(self._block), 0
        self._iu += 1
        return float(self._u[self._iu - 1])

    def normal(self) -> float:
        if self._in == self._n.size:
            self._n, self._in = self._rng.standard_normal(self._block), 0
        self._in += 1
        return float(self._n[self._in - 1])

    def exponential(self) -> float:
        if self._ie == self._e.size:
            self._e, self._ie = self._rng.standard_exponential(self._block), 0
        self._ie += 1
        return float(self._e[self._ie - 1])


def _coupled(
    names: list[str],
    primaries: list[np.ndarray],
    dets: list[DetectorParams],
    crosstalk: Sequence[Crosstalk],
    rng: np.random.Generator,
) -> list[np.ndarray]:
    """Dead time, cascading afterpulses and cascading crosstalk across channels, event by event.

    Primaries of all channels are merged in time order (ties: channel order);
    triggered candidates (afterpulses, crosstalk) wait in a min-heap and go
    first at equal times, as in `_self_exciting`.
    """
    index = {name: i for i, name in enumerate(names)}
    links: list[list[tuple[int, float, float, float]]] = [[] for _ in names]
    for link in crosstalk:
        links[index[link.source]].append(
            (index[link.target], float(link.probability), float(link.delay_ps), float(link.jitter_sigma_ps))
        )
    times = np.concatenate(primaries) if primaries else np.empty(0, dtype=np.int64)
    chans = np.concatenate([np.full(p.size, i, dtype=np.int64) for i, p in enumerate(primaries)])
    order = np.argsort(times, kind="stable")
    prim_t, prim_c = times[order].tolist(), chans[order].tolist()
    # A zero dead time still records at most one tag per ps per channel (as `detect` does).
    dead = [max(int(d.dead_time_ps), 1) for d in dets]
    para = [int(d.dead_time_p_ps) for d in dets]
    last_input = [-(1 << 62)] * len(names)
    ap = [(float(d.afterpulse_prob), float(d.afterpulse_tau_ps)) for d in dets]
    draws = _Draws(rng)
    out: list[list[int]] = [[] for _ in names]
    next_allowed = [-(1 << 62)] * len(names)
    pending: list[tuple[int, int]] = []  # (time, channel) of triggered candidates
    i, n = 0, len(prim_t)
    while i < n or pending:
        if pending and (i >= n or pending[0][0] <= prim_t[i]):
            t, c = heapq.heappop(pending)
        else:
            t, c = prim_t[i], prim_c[i]
            i += 1
        if para[c]:
            previous, last_input[c] = last_input[c], t
            if t - previous < para[c]:
                continue
        if t < next_allowed[c]:
            continue
        out[c].append(t)
        next_allowed[c] = t + dead[c]
        p_ap, tau = ap[c]
        if p_ap > 0 and draws.uniform() < p_ap:
            heapq.heappush(pending, (t + dets[c].dead_time_ps + round(tau * draws.exponential()), c))
        for target, p, delay, sigma in links[c]:
            if draws.uniform() < p:
                offset = delay + sigma * draws.normal() if sigma > 0 else delay
                heapq.heappush(pending, (t + max(0, round(offset)), target))
    return [np.array(o, dtype=np.int64) for o in out]


def _self_exciting(
    arrivals: np.ndarray, dead_time_ps: int, p: float, tau_ps: float, rng: np.random.Generator,
    dead_time_p_ps: int = 0,
) -> np.ndarray:
    """Dead time plus cascading afterpulses, event by event.

    `arrivals` are sorted primary avalanche candidates (photons and dark
    counts). Each recorded avalanche blocks the next `dead_time_ps` and, with
    probability `p`, adds an afterpulse candidate at its dead-time end plus an
    exponential delay with mean `tau_ps`. With `dead_time_p_ps` every candidate
    first meets the paralyzable stage (see the module docstring).
    """
    last_input: int | None = None
    out: list[int] = []
    pending: list[int] = []  # afterpulse candidates (min-heap)
    next_allowed: int | None = None
    block = 4096
    u = rng.random(block)
    e = rng.exponential(tau_ps, block)
    k = 0
    primaries = arrivals.tolist()
    i, n = 0, len(primaries)
    while i < n or pending:
        if pending and (i >= n or pending[0] <= primaries[i]):
            t = heapq.heappop(pending)
        else:
            t = primaries[i]
            i += 1
        if dead_time_p_ps:
            previous, last_input = last_input, t
            if previous is not None and t - previous < dead_time_p_ps:
                continue
        if next_allowed is not None and t < next_allowed:
            continue
        out.append(t)
        next_allowed = t + dead_time_ps
        if k == block:
            u = rng.random(block)
            e = rng.exponential(tau_ps, block)
            k = 0
        if u[k] < p:
            heapq.heappush(pending, next_allowed + round(float(e[k])))
        k += 1
    return np.array(out, dtype=np.int64)


def simulate_dark(det: DetectorParams, duration_s: float, rng: np.random.Generator) -> np.ndarray:
    """Recorded tags of one detector with no light: dark counts, dead time and afterpulsing."""
    _checks.finite(duration_s, "duration_s", minimum=0.0)
    if duration_s <= 0:
        raise ValueError("duration_s must be positive")
    return detect(np.empty(0), det, duration_s, rng)


def simulate_pairs(
    pair_rate_per_s: float,
    det_a: DetectorParams,
    det_b: DetectorParams,
    duration_s: float,
    rng: np.random.Generator,
    crosstalk: Sequence[Crosstalk] = (),
) -> tuple[np.ndarray, np.ndarray]:
    """Recorded tags of two detectors, each receiving one photon of every pair (a CW SPDC source).

    Channel b's photon arrives at the same time as channel a's; the relative
    delay of the recorded tags is det_b.delay_ps - det_a.delay_ps and the
    width of the coincidence peak comes from both detectors' jitter.
    `crosstalk` links name the channels "a" and "b".
    """
    _checks.finite(duration_s, "duration_s", minimum=0.0)
    if duration_s <= 0:
        raise ValueError("duration_s must be positive")
    t = poisson_times(pair_rate_per_s, duration_s, rng)
    out = detect_channels({"a": t, "b": t}, {"a": det_a, "b": det_b}, duration_s, rng, crosstalk)
    return out["a"], out["b"]


def simulate_hom(
    source: PairSource,
    det_a: DetectorParams,
    det_b: DetectorParams,
    path_delay_ps: float,
    duration_s: float,
    rng: np.random.Generator,
    crosstalk: Sequence[Crosstalk] = (),
) -> dict[str, np.ndarray]:
    """Recorded tags on the two outputs ("A", "B") of a Hong-Ou-Mandel beamsplitter.

    One photon of each pair is delayed by `path_delay_ps`. The photons leave
    by different ports with probability (1 - V exp(-(delay / tau_c)^2)) / 2,
    V = source.visibility and tau_c = source.coherence_time_ps; otherwise
    both leave by the same port, where the detector (not photon-number
    resolving) records at most one of them. `crosstalk` links name the
    channels "A" and "B".
    """
    _checks.finite(path_delay_ps, "path_delay_ps")
    _checks.finite(duration_s, "duration_s", minimum=0.0)
    if duration_s <= 0:
        raise ValueError("duration_s must be positive")
    t1 = poisson_times(source.pair_rate_per_s, duration_s, rng)
    t2 = t1 + path_delay_ps
    g = math.exp(-((path_delay_ps / source.coherence_time_ps) ** 2))
    split = rng.random(t1.size) < 0.5 * (1.0 - source.visibility * g)
    first_to_a = rng.random(t1.size) < 0.5
    to_a = np.concatenate([
        t1[split & first_to_a], t2[split & ~first_to_a],  # one photon each way
        t1[~split & first_to_a], t2[~split & first_to_a],  # both to A
    ])
    to_b = np.concatenate([
        t2[split & first_to_a], t1[split & ~first_to_a],
        t1[~split & ~first_to_a], t2[~split & ~first_to_a],  # both to B
    ])
    photons = {"A": np.sort(to_a), "B": np.sort(to_b)}
    return detect_channels(photons, {"A": det_a, "B": det_b}, duration_s, rng, crosstalk)


def simulate_chsh_setting(
    source: PairSource,
    alpha_deg: float,
    beta_deg: float,
    detectors: dict[str, DetectorParams],
    duration_s: float,
    rng: np.random.Generator,
    crosstalk: Sequence[Crosstalk] = (),
) -> dict[str, np.ndarray]:
    """Recorded tags on channels A+, A-, B+, B- for one pair of analyzer angles (`crosstalk` among them)."""
    missing = set(CHANNELS) - set(detectors)
    if missing:
        raise ValueError(f"missing detector parameters for {sorted(missing)}")
    _checks.finite(duration_s, "duration_s", minimum=0.0)
    if duration_s <= 0:
        raise ValueError("duration_s must be positive")
    t = poisson_times(source.pair_rate_per_s, duration_s, rng)
    n = t.size
    a_plus = rng.random(n) < 0.5
    p_same = 0.5 * (1.0 + source.visibility * math.cos(2.0 * math.radians(alpha_deg - beta_deg)))
    same = rng.random(n) < p_same
    b_plus = np.where(same, a_plus, ~a_plus)
    photons = {"A+": t[a_plus], "A-": t[~a_plus], "B+": t[b_plus], "B-": t[~b_plus]}
    return detect_channels({ch: photons[ch] for ch in CHANNELS}, detectors, duration_s, rng, crosstalk)


def simulate_chsh(
    source: PairSource,
    detectors: dict[str, DetectorParams],
    duration_s: float,
    rng: np.random.Generator,
    settings: dict[str, float] | None = None,
    crosstalk: Sequence[Crosstalk] = (),
) -> dict[tuple[str, str], dict[str, np.ndarray]]:
    """One acquisition per CHSH setting pair, keyed like ("a", "b'")."""
    angles = CHSH_SETTINGS if settings is None else settings
    missing = [s for s in CHSH_SETTINGS if s not in angles]
    if missing:
        raise ValueError(f"settings are missing angle(s) for {missing}")
    return {
        (sa, sb): simulate_chsh_setting(source, angles[sa], angles[sb], detectors, duration_s, rng, crosstalk)
        for sa, sb in CHSH_PAIRS
    }
