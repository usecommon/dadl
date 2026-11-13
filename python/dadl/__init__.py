"""DADL: Detector-Aware Data Layer for quantum optics experiments.

Calibrations live in `dadl.calibration` and artifacts in `dadl.artifact`;
they are not imported here. (SciPy is loaded for the exact visibility
intervals.)

`dadl.backend()` says which implementation runs the per-tag functions:
"rust" when the optional compiled core (`pip install dadl-quantum-core`) is
installed, else "python". Set DADL_BACKEND=python or =rust to choose.
"""

__version__ = "0.0.3"

from dadl._backend import name as backend  # noqa: E402
from dadl.analysis.chsh import Window, chsh  # noqa: E402
from dadl.analysis.hom import hom_from_tags  # noqa: E402
from dadl.analysis.visibility import bell_visibility, hom_visibility  # noqa: E402
from dadl.config import Setup, load_setup  # noqa: E402
from dadl.data import ChannelId, Dataset, TagSet  # noqa: E402
from dadl.errors import ConfigError, DadlError, GuardrailRefusal, InsufficientCounts  # noqa: E402
from dadl.guardrails import Policy  # noqa: E402
from dadl.result import Estimate, LabelledValue  # noqa: E402
from dadl.tagfile import read_tags, write_tags  # noqa: E402

__all__ = [
    "ChannelId",
    "ConfigError",
    "Dataset",
    "DadlError",
    "Estimate",
    "GuardrailRefusal",
    "InsufficientCounts",
    "LabelledValue",
    "Policy",
    "Setup",
    "TagSet",
    "Window",
    "__version__",
    "backend",
    "bell_visibility",
    "chsh",
    "hom_from_tags",
    "hom_visibility",
    "load_setup",
    "read_tags",
    "write_tags",
]
