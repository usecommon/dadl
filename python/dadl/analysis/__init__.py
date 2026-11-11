"""Analyses. Each returns a LabelledValue with the raw value first.

The analysis type is fixed by the function, never passed in by the caller.
"""

from dadl.analysis.chsh import Window, chsh
from dadl.analysis.visibility import bell_visibility, hom_visibility

__all__ = ["Window", "bell_visibility", "chsh", "hom_visibility"]
