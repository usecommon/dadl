"""Analyses. Each returns a LabelledValue with the raw value first."""

from dadl.analysis.chsh import Window, chsh
from dadl.analysis.visibility import two_photon_visibility

__all__ = ["Window", "chsh", "two_photon_visibility"]
