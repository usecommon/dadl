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
