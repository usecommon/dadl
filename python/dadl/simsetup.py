"""Simulate the acquisitions a setup file describes (the physics half of the hardware twin).

`simulate_setup` returns the `runs` structure `dadl.tagfile` reads and
writes: {setting key: {channel name: int64 ps tags}}. `dadl simulate` writes
it as a generic tag file, or as one PicoQuant PTU file per acquisition
(`dadl.adapters.ptu_write`), so simulated data can enter DADL exactly as a
recording would.
"""

from __future__ import annotations

import numpy as np

from dadl.config import Setup
from dadl.errors import DadlError
from dadl.sim import simulate_chsh, simulate_dark, simulate_hom


def _require_sim(setup: Setup) -> None:
    if setup.backend != "sim":
        raise DadlError(f"{setup.name}: backend is {setup.backend!r}; only 'sim' setups can be simulated")


def simulate_setup(setup: Setup) -> dict:
    """Every acquisition of the setup's analysis (CHSH: four setting pairs; HOM: dip and reference)."""
    _require_sim(setup)
    rng = np.random.default_rng(setup.seed)
    if setup.analysis_type == "hom":
        delays = {"dip": 0.0, "reference": setup.reference_delay_ps}
        return {
            key: simulate_hom(
                setup.source, setup.detectors["A"], setup.detectors["B"], delays[key],
                setup.duration_s_per_setting, rng,
            )
            for key in setup.setting_keys
        }
    return simulate_chsh(setup.source, setup.detectors, setup.duration_s_per_setting, rng, settings=setup.settings_deg)


def simulate_setup_dark(setup: Setup, duration_s: float) -> dict:
    """One dark acquisition per channel (source off), as the first setting; the others are empty."""
    _require_sim(setup)
    if not duration_s > 0:
        raise DadlError("duration_s must be positive")
    streams = np.random.SeedSequence(setup.seed).spawn(len(setup.channels))
    dark = {
        name: simulate_dark(setup.detectors[name], duration_s, np.random.default_rng(stream))
        for (_num, name), stream in zip(sorted(setup.channels.items()), streams)
    }
    empty = {name: np.empty(0, dtype=np.int64) for name in dark}
    return {key: (dark if i == 0 else empty) for i, key in enumerate(setup.setting_keys)}


#: Default input rates of a dead-time scan, as fractions of 1 / dead_time_ps of the slowest channel.
DEADTIME_SCAN_FRACTIONS = (0.05, 0.15, 0.3, 0.5)


def deadtime_scan_rates(setup: Setup) -> tuple[float, ...]:
    """Default photon rates (per s, at each detector) of a dead-time scan: one per setting, up to
    four, spread so that rate x dead time runs from 0.05 to 0.5 (paralysis grows with it)."""
    longest = max((d.dead_time_ps for d in setup.detectors.values()), default=0)
    if longest <= 0:
        raise DadlError("the setup's detectors have no dead time to scan")
    k = min(len(setup.setting_keys), len(DEADTIME_SCAN_FRACTIONS))
    fractions = DEADTIME_SCAN_FRACTIONS[:k] if k > 2 else (0.1, 0.4)
    return tuple(f / (longest * 1e-12) for f in fractions)


def simulate_setup_deadtime_scan(setup: Setup, rates: tuple[float, ...], duration_s: float) -> dict:
    """Steady light straight onto every detector, one rate per setting (for `calibrate deadtime`).

    `rates` are photon rates at each detector (per s, before its efficiency); settings past the
    last rate are empty.
    """
    from dadl.sim import detect, poisson_times

    _require_sim(setup)
    if not duration_s > 0:
        raise DadlError("duration_s must be positive")
    if not 2 <= len(rates) <= len(setup.setting_keys):
        raise DadlError(f"give 2 to {len(setup.setting_keys)} rates (one per setting of this setup), got {len(rates)}")
    if any(not r > 0 for r in rates):
        raise DadlError("rates must be positive")
    names = [name for _num, name in sorted(setup.channels.items())]
    streams = np.random.SeedSequence(setup.seed).spawn(len(rates) * len(names))
    out = {}
    for i, key in enumerate(setup.setting_keys):
        if i >= len(rates):
            out[key] = {name: np.empty(0, dtype=np.int64) for name in names}
            continue
        run = {}
        for j, name in enumerate(names):
            rng = np.random.default_rng(streams[i * len(names) + j])
            run[name] = detect(poisson_times(rates[i], duration_s, rng), setup.detectors[name], duration_s, rng)
        out[key] = run
    return out


def simulate_setup_heralding(setup: Setup, duration_s: float) -> dict:
    """A heralding run: the pair source straight onto one detector per arm, no analysers or beamsplitter.

    CHSH setups use A+ (herald) and B+ (signal); HOM setups use A and B. The
    run is the first setting; the others are empty.
    """
    from dadl.sim import simulate_pairs

    _require_sim(setup)
    if not duration_s > 0:
        raise DadlError("duration_s must be positive")
    herald, signal = ("A", "B") if setup.analysis_type == "hom" else ("A+", "B+")
    rng = np.random.default_rng(setup.seed)
    a, b = simulate_pairs(setup.source.pair_rate_per_s, setup.detectors[herald], setup.detectors[signal], duration_s, rng)
    run = {name: np.empty(0, dtype=np.int64) for name in setup.channels.values()}
    run[herald], run[signal] = a, b
    empty = {name: np.empty(0, dtype=np.int64) for name in run}
    return {key: (run if i == 0 else empty) for i, key in enumerate(setup.setting_keys)}
