"""Off-peak accidental coincidences for a fixed-grid coincidence window.

Method `offpeak-grid`: the same fixed-grid coincidence count the analyses
use (`dadl.coincidences.grid_coincidences`, declared width and delays) is
repeated with channel b shifted by many window widths, where no pair can
coincide. Each shifted count is an independent estimate of the accidentals
in the on-peak window (docs/design.md, "Accidentals"):

- CW pumps: shifts of guard + j * width for j = 0 .. n_side - 1, and
  -(guard + j * width) for j = n_side .. 2 n_side - 1;
- pulsed pumps: the side peaks, shifts of k * rep_period for the n_side
  multiples beyond the guard, and -k * rep_period for the next n_side.

The two sides never mirror each other. Two pairs a time o apart give one
accidental at +o and another at -o, so mirrored windows are positively
correlated and would make the interval too narrow (seen as 91.7% coverage in
the first run of the phase 2 gate).

The guard (default 200 ns) keeps the shifted windows clear of dead-time and
afterpulse structure next to the peak. With M shifted counts summing to S,
the expected accidental count lambda has the Jeffreys Gamma(S + 1/2, M)
posterior. R1 * R2 * width * duration is reported only as a cross-check in
the notes, since it assumes stationary, uncorrelated streams.

The predictive check is a dispersion test of the M shifted counts (chi-square
with M - 1 degrees of freedom): structure in the off-peak region, such as
drift or a periodic background missed by a CW setting, shows as
over-dispersion.
"""

from __future__ import annotations

import numpy as np
from scipy import stats as _st

from dadl import _checks
from dadl.calibration.base import N_DRAWS, Calibration
from dadl.coincidences import occupied_slots
from dadl.data import ChannelId, tags_sha256
from dadl.histogram import cross_correlation
from dadl.stats import gamma_posterior

KIND = "accidentals"
METHOD = "offpeak-grid"
METHOD_VERSION = "1"
GUARD_PS = 200_000
N_SIDE = 100

UNITS = {"accidental_coincidences": "counts per acquisition"}


def side_offsets(width_ps: int, *, rep_period_ps: int | None = None, n_side: int = N_SIDE, guard_ps: int = GUARD_PS) -> np.ndarray:
    """Shifts of channel b (ps) at which off-peak coincidences are counted: n_side on each side, never mirrored."""
    width_ps = _checks.integer(width_ps, "width_ps", minimum=1)
    n_side = _checks.integer(n_side, "n_side", minimum=1)
    guard_ps = _checks.integer(guard_ps, "guard_ps", minimum=0)
    if rep_period_ps is None:
        start = -(-max(guard_ps, width_ps) // width_ps) * width_ps
        steps = start + width_ps * np.arange(2 * n_side, dtype=np.int64)
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
    offsets = side_offsets(width_ps, rep_period_ps=rep_period_ps, n_side=n_side, guard_ps=guard_ps)
    counts = shifted_coincidences(a, b, width_ps, delay_a_ps, delay_b_ps, offsets)
    total, m = int(counts.sum()), counts.size
    rng = np.random.default_rng(seed)
    est, draws = gamma_posterior(total, float(m), rng, n_draws)
    notes = [f"{total} coincidences in {m} shifted windows ({'pulsed' if rep_period_ps else 'CW'} offsets)"]
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
        samples=draws[:, None],
        sample_params=("accidental_coincidences",),
        converged=True,
        ppc_pvalue=dispersion_pvalue(counts),
        inputs=((tags_sha256(a), "raw-tags"), (tags_sha256(b), "raw-tags")),
        backend=backend,
        conditions={"width_ps": float(width_ps), **(conditions or {})},
        notes=tuple(notes),
    )


__all__ = ["KIND", "calibrate_accidentals", "shifted_coincidences", "side_offsets"]
