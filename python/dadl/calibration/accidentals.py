"""Off-peak accidental coincidences for a fixed-grid coincidence window.

Method `offpeak-grid`: the same fixed-grid coincidence count the analyses
use (`dadl.coincidences.grid_coincidences`, declared width and delays) is
repeated with channel b shifted by many window widths, where no pair can
coincide. Each shifted count is an independent estimate of the accidentals
in the on-peak window (docs/design.md, "Accidentals"):

- CW pumps: shifts of guard + j * step for j = 0 .. n_side - 1, and
  -(guard + j * step) for j = n_side .. 2 n_side - 1, with step about
  guard / n_side in whole windows, so the windows span one to three guards;
- pulsed pumps: the side peaks, shifts of k * rep_period for the n_side
  multiples beyond the guard, and -k * rep_period for the next n_side.

The two sides never mirror each other. Two pairs a time o apart give one
accidental at +o and another at -o, so mirrored windows are positively
correlated and would make the interval too narrow (seen as 91.7% coverage in
the first run of the phase 2 gate).

The guard keeps the shifted windows clear of dead-time and afterpulse
structure that follows true coincidences: afterpulses after a pair detection
land in windows a few afterpulse time constants away and bias lambda high (by
9% with a 50 ns, 3% afterpulse at a 200 ns guard; method v1). The default is
10 us; given the channel's afterpulse time constant and dead time (from its
dark calibration) it is at least dead time + 10 tau_a. With M shifted counts summing to S,
the expected accidental count lambda has the Jeffreys Gamma(S + 1/2, M)
posterior. R1 * R2 * width * duration is reported only as a cross-check in
the notes, since it assumes stationary, uncorrelated streams.

The predictive check combines two tests of the M shifted counts (the smaller
p-value, doubled): dispersion (structure such as drift or a periodic
background missed by a CW setting) and a trend test, nearer windows against
farther ones, which catches a smooth afterpulse tail that dispersion misses.
"""

from __future__ import annotations

import numpy as np
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.base import N_DRAWS, Calibration, with_interruption_check, draws_rng
from dadl.coincidences import occupied_slots
from dadl.data import ChannelId, tags_sha256
from dadl.histogram import cross_correlation
from dadl.stats import gamma_posterior

KIND = "accidentals"
METHOD = "offpeak-grid"
METHOD_VERSION = "2"
GUARD_PS = 10_000_000
#: Afterpulse time constants the guard keeps clear of.
GUARD_TAUS = 10
N_SIDE = 100

UNITS = {"accidental_coincidences": "{count}"}
DEFINITIONS = {
    "accidental_coincidences": "expected accidental coincidences in the declared fixed-grid window over this acquisition"
}


def side_offsets(width_ps: int, *, rep_period_ps: int | None = None, n_side: int = N_SIDE, guard_ps: int = GUARD_PS) -> np.ndarray:
    """Shifts of channel b (ps) at which off-peak coincidences are counted: n_side on each side, never mirrored."""
    width_ps = _checks.integer(width_ps, "width_ps", minimum=1)
    n_side = _checks.integer(n_side, "n_side", minimum=1)
    guard_ps = _checks.integer(guard_ps, "guard_ps", minimum=0)
    if rep_period_ps is None:
        start = -(-max(guard_ps, width_ps) // width_ps) * width_ps
        step = max(1, -(-guard_ps // (n_side * width_ps))) * width_ps
        steps = start + step * np.arange(2 * n_side, dtype=np.int64)
    else:
        period = _checks.integer(rep_period_ps, "rep_period_ps", minimum=1)
        first = max(1, -(-guard_ps // period))
        steps = period * np.arange(first, first + 2 * n_side, dtype=np.int64)
    return np.concatenate([-steps[n_side:][::-1], steps[:n_side]])


def shifted_coincidences(
    a: np.ndarray, b: np.ndarray, width_ps: int, delay_a_ps: int, delay_b_ps: int, offsets: np.ndarray
) -> np.ndarray:
    """grid_coincidences(a, b, width, delay_a, delay_b + offset) for each offset."""
    offsets = np.asarray(offsets, dtype=np.int64)
    slots_a = occupied_slots(a, width_ps, delay_a_ps)
    base_b = occupied_slots(b, width_ps, delay_b_ps)
    whole = offsets % width_ps == 0
    out = np.empty(offsets.size, dtype=np.int64)
    if whole.any():
        # Shifting b by a whole number of windows moves its slots by -offset / width, so it
        # coincides with a where base_b - slot_a = offset / width: one histogram gives every count.
        shift = offsets[whole] // width_ps
        reach = int(np.abs(shift).max()) + 1
        hist = cross_correlation(slots_a, base_b, 1, reach)
        out[whole] = hist[(shift + reach).astype(np.intp)].astype(np.int64)
    for i in np.flatnonzero(~whole):
        slots_b = occupied_slots(b, width_ps, delay_b_ps + int(offsets[i]))
        out[i] = np.intersect1d(slots_a, slots_b, assume_unique=True).size
    return out


#: Below this mean count per window the chi-square approximation to the dispersion test fails
#: (it flagged 4.5% of good sparse runs at p < 0.01), so the exact conditional test is used.
EXACT_BELOW = 5.0


def guard_for(afterpulse_tau_ps: float | None = None, dead_time_ps: float | None = None, guard_ps: int = GUARD_PS) -> int:
    """The guard: at least dead time + GUARD_TAUS afterpulse time constants, and never below `guard_ps`."""
    if afterpulse_tau_ps is None:
        return int(guard_ps)
    return int(max(guard_ps, (dead_time_ps or 0.0) + GUARD_TAUS * float(afterpulse_tau_ps)))


def trend_pvalue(counts: np.ndarray, offsets: np.ndarray) -> float:
    """One-sided p-value that the nearer half of the windows holds more counts than the farther half."""
    order = np.argsort(np.abs(np.asarray(offsets)), kind="stable")
    half = order.size // 2
    near, far = int(np.sum(np.asarray(counts)[order[:half]])), int(np.sum(np.asarray(counts)[order[half : 2 * half]]))
    if near + far == 0:
        return 1.0
    return float(_st.binom.sf(near - 1, near + far, 0.5))


def offpeak_pvalue(counts: np.ndarray, offsets: np.ndarray) -> float:
    """Predictive check of the shifted counts: min(dispersion, trend) p-value, doubled (Bonferroni)."""
    return min(1.0, 2.0 * min(dispersion_pvalue(counts), trend_pvalue(counts, offsets)))


def dispersion_pvalue(counts: np.ndarray, n_sim: int = 2000) -> float:
    """Upper-tail p-value of the Poisson dispersion statistic (small values flag structure).

    With few counts per window the p-value is exact: given the total, Poisson counts are
    multinomial with equal cells, which is simulated with a fixed seed.
    """
    counts = np.asarray(counts)
    m = float(np.mean(counts))
    if counts.size < 2 or m <= 0:
        return 1.0
    stat = float(np.sum((counts - m) ** 2) / m)
    if m >= EXACT_BELOW:
        return float(_st.chi2.sf(stat, counts.size - 1))
    total = int(counts.sum())
    sims = np.random.default_rng(0).multinomial(total, np.full(counts.size, 1.0 / counts.size), size=n_sim)
    sim_stat = np.sum((sims - m) ** 2, axis=1) / m
    return float((1 + np.sum(sim_stat >= stat - 1e-9)) / (n_sim + 1))


def _span_s(*streams: np.ndarray) -> float:
    lo = min(int(s[0]) for s in streams if s.size)
    hi = max(int(s[-1]) for s in streams if s.size)
    return (hi - lo) / 1e12


def side_exposure(offsets: np.ndarray, *streams: np.ndarray) -> float:
    """Exposure of the shifted windows in on-peak window units: sum of (1 - |offset| / T).

    A window shifted by o loses the |o| / T of the record that falls past its end, which the
    on-peak window does not (statistician review, cycle 1: -0.65% at 10 Mcps over 3 ms,
    coverage 37%). Negligible for acquisitions of a second or more.
    """
    span_ps = _span_s(*streams) * 1e12 if any(s.size for s in streams) else 0.0
    offsets = np.abs(np.asarray(offsets, dtype=np.float64))
    if span_ps <= 0:
        return float(offsets.size)
    return float(np.clip(1.0 - offsets / span_ps, 0.0, 1.0).sum())


def calibrate_accidentals(
    a: np.ndarray,
    b: np.ndarray,
    channels: tuple[ChannelId, ChannelId],
    *,
    width_ps: int,
    delay_a_ps: int = 0,
    delay_b_ps: int = 0,
    rep_period_ps: int | None = None,
    n_side: int = N_SIDE,
    guard_ps: int = GUARD_PS,
    afterpulse_tau_ps: float | None = None,
    dead_time_ps: float | None = None,
    backend: str = "hardware",
    seed: int | None = None,
    n_draws: int = N_DRAWS,
    conditions: dict[str, float] | None = None,
) -> Calibration:
    """Expected accidental coincidences in the declared on-peak window over this acquisition."""
    if len(channels) != 2:
        raise ValueError("channels must be (channel a, channel b)")
    a = _checks.timestamps(a, "a")
    b = _checks.timestamps(b, "b")
    _checks.check_sorted(a, "a")
    _checks.check_sorted(b, "b")
    guard = guard_for(afterpulse_tau_ps, dead_time_ps, guard_ps)
    offsets = side_offsets(width_ps, rep_period_ps=rep_period_ps, n_side=n_side, guard_ps=guard)
    counts = shifted_coincidences(a, b, width_ps, delay_a_ps, delay_b_ps, offsets)
    total, m = int(counts.sum()), counts.size
    rng = draws_rng(seed, a, b)
    est, draws = gamma_posterior(total, side_exposure(offsets, a, b), rng, n_draws)
    notes = [
        f"{total} coincidences in {m} shifted windows ({'pulsed' if rep_period_ps else 'CW'} offsets) "
        f"from {int(np.abs(offsets).min())} to {int(np.abs(offsets).max())} ps"
    ]
    if afterpulse_tau_ps is None:
        notes.append(
            f"guard {guard} ps without a measured afterpulse time constant; afterpulsing slower than about "
            f"{guard / GUARD_TAUS:.3g} ps would bias the estimate high (the trend check looks for it)"
        )
    if a.size and b.size:
        t = _span_s(a, b)
        if t > 0:
            r1, r2 = a.size / t, b.size / t
            notes.append(f"cross-check only: R1 R2 width duration = {r1 * r2 * width_ps * 1e-12 * t:.4g}")
    return Calibration(
        kind=KIND,
        method=METHOD,
        method_version=METHOD_VERSION,
        channels=tuple(channels),
        estimates={"accidental_coincidences": est},
        units=dict(UNITS),
        definitions=dict(DEFINITIONS),
        samples=draws[:, None],
        sample_params=("accidental_coincidences",),
        converged=True,
        ppc_pvalue=with_interruption_check(offpeak_pvalue(counts, offsets), notes, a, b),
        inputs=((tags_sha256(a), "raw-tags"), (tags_sha256(b), "raw-tags")),
        backend=backend,
        conditions={"width_ps": float(width_ps), "guard_ps": float(guard), **(conditions or {})},
        notes=tuple(notes),
    )


__all__ = ["KIND", "calibrate_accidentals", "shifted_coincidences", "side_offsets"]
