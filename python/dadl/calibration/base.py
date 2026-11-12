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
    #: (sha256, role) of every input.
    inputs: tuple[tuple[str, str], ...]
    backend: str = "hardware"
    conditions: Mapping[str, float] = field(default_factory=dict)
    #: Parameters the data could not identify (reported as limits or not at all).
    unidentified: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
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

    def __str__(self) -> str:
        lines = [f"{self.kind} ({self.method} {self.method_version}) for {', '.join(map(str, self.channels))}"]
        for name, iv in self.estimates.items():
            lines.append(f"  {name:<20} {iv}  {self.units.get(name, '')}".rstrip())
        lines.append(f"  converged: {self.converged}; predictive check p = {self.ppc_pvalue:.3f}")
        lines.extend(f"  note: {n}" for n in self.notes)
        return "\n".join(lines)
