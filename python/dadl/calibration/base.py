"""The result type shared by all calibrations."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np

from dadl.data import ChannelId
from dadl.stats import Interval

#: Posterior draws stored per artifact (docs/design.md, "Inference").
N_DRAWS = 2000
#: A predictive-check p-value below this marks the calibration a misfit: the model does not
#: describe the data, and its intervals are not to be trusted (docs/design.md: artifacts are
#: accepted only after the predictive check passes). 0.1% keeps false alarms rare.
MISFIT_ALPHA = 0.001


def draws_rng(seed: int | None, *inputs: object) -> np.random.Generator:
    """The generator for posterior draws: `seed` mixed with a hash of the inputs.

    Software review (cycle 2): with seed=None the same data gave a different artifact ID on every
    run, so the inputs always enter. Bayesian review (cycle 3): with only the shared setup seed,
    two channels' artifacts had draws correlated at 0.93-1.00 (the spread of a difference 180x too
    small), so the inputs enter even when a seed is given.
    """
    import hashlib

    from dadl.data import tags_sha256

    h = hashlib.sha256()
    for x in inputs:
        h.update((tags_sha256(x) if isinstance(x, np.ndarray) else repr(x)).encode())
    words = [int.from_bytes(h.digest()[i : i + 4], "little") for i in range(0, 16, 4)]
    entropy = words if seed is None else [int(seed) & 0xFFFFFFFF, int(seed) >> 32 & 0xFFFFFFFF, *words]
    return np.random.default_rng(np.random.SeedSequence(entropy))


def with_interruption_check(ppc: float, notes: list[str], *streams: np.ndarray) -> float:
    """Combine a model check with the per-channel interruption test (smaller p doubled).

    Appends a note when a channel's longest gap points to a latch or interruption.
    """
    from dadl.stationarity import channels_pvalue

    p_int, longest = channels_pvalue(*streams)
    if 2.0 * p_int < MISFIT_ALPHA:
        # Same threshold as the status below (point-process review, cycle 5: the note fired at
        # p < 0.001 while the status needed p < 0.0005). A strongly modulated source (+-50%) can
        # also trip this test.
        notes.append(
            f"a channel has a {longest / 1e12:.3g} s gap, far longer than its count rate allows (p = {p_int:.2g}): "
            "the rate was not stationary or the acquisition was interrupted (latch, overload shutdown, lost "
            "data); the rates and ratios include the dead stretch"
        )
    return min(1.0, 2.0 * min(ppc, p_int))


@dataclass(frozen=True)
class Calibration:
    """Estimates, intervals and posterior draws for one calibration run."""

    kind: str
    method: str
    method_version: str
    channels: tuple[ChannelId, ...]
    #: Point estimate and 95% interval per parameter, in the units named by `units`.
    estimates: Mapping[str, Interval]
    units: Mapping[str, str]
    #: Posterior draws, one row per draw, columns in `sample_params` order (natural units).
    samples: np.ndarray
    sample_params: tuple[str, ...]
    converged: bool
    #: Posterior predictive p-value of the fit (small values flag misfit).
    ppc_pvalue: float
    #: (sha256, role) or (sha256, role, name) of every input.
    inputs: tuple[tuple[str, ...], ...]
    backend: str = "hardware"
    conditions: Mapping[str, float] = field(default_factory=dict)
    #: Parameters the data could not identify (reported as limits or not at all).
    unidentified: tuple[str, ...] = ()
    #: Draw columns that are conservative stand-ins (e.g. uniform below an upper limit), not a
    #: posterior: phase 3 must not propagate them as one (Bayesian review, cycle 3).
    stand_in: tuple[str, ...] = ()
    #: What each parameter is, beside its unit (calibration-lab review, cycle 3: units were prose).
    definitions: Mapping[str, str] = field(default_factory=dict)
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        # A calibration of simulated channels is simulated, whatever the default backend said (API
        # review, cycle 7: device "sim" was written with backend "hardware").
        if self.backend == "hardware" and self.channels and all(ch.device == "sim" for ch in self.channels):
            object.__setattr__(self, "backend", "sim")
        samples = np.asarray(self.samples, dtype=float)
        if samples.ndim != 2 or samples.shape[1] != len(self.sample_params):
            raise ValueError("samples must be a 2-D array with one column per sample parameter")
        samples = samples.copy()
        samples.setflags(write=False)
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "estimates", MappingProxyType(dict(self.estimates)))
        object.__setattr__(self, "units", MappingProxyType(dict(self.units)))
        object.__setattr__(self, "conditions", MappingProxyType(dict(self.conditions)))
        if self.backend not in ("sim", "hardware"):
            raise ValueError(f"backend must be 'sim' or 'hardware', got {self.backend!r}")

    @property
    def status(self) -> str:
        """"ok", or "misfit" when the predictive check fails (p < MISFIT_ALPHA) or the fit did not
        converge (second red-team review, cycle 9: a non-converged delay fit was written as ok)."""
        return "misfit" if self.ppc_pvalue < MISFIT_ALPHA or not self.converged else "ok"

    def __str__(self) -> str:
        lines = [f"{self.kind} ({self.method} {self.method_version}) for {', '.join(map(str, self.channels))}"]
        for name, iv in self.estimates.items():
            shown = f"<= {iv.high:.6g} (upper limit only)" if name in self.unidentified and iv.low == 0 else str(iv)
            lines.append(f"  {name:<20} {shown}  {self.units.get(name, '')}".rstrip())
        lines.append(f"  converged: {self.converged}; predictive check p = {self.ppc_pvalue:.3g}")
        if self.status == "misfit":
            lines.append(
                "  MISFIT: the model does not describe these data (predictive check p < "
                f"{MISFIT_ALPHA:g}); do not use these intervals"
            )
        lines.extend(f"  note: {n}" for n in self.notes)
        return "\n".join(lines)
