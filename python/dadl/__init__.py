"""DADL: Detector-Aware Data Layer for quantum optics experiments."""

__version__ = "0.0.3"

from dadl.analysis.chsh import Window, chsh  # noqa: E402
from dadl.analysis.visibility import bell_visibility, hom_visibility  # noqa: E402
from dadl.config import Setup, load_setup  # noqa: E402
from dadl.errors import ConfigError, DadlError, GuardrailRefusal, InsufficientCounts  # noqa: E402
from dadl.guardrails import Policy  # noqa: E402
from dadl.result import Estimate, LabelledValue  # noqa: E402
from dadl.tagfile import read_tags, write_tags  # noqa: E402

__all__ = [
    "ConfigError",
    "DadlError",
    "Estimate",
    "GuardrailRefusal",
    "InsufficientCounts",
    "LabelledValue",
    "Policy",
    "Setup",
    "Window",
    "__version__",
    "bell_visibility",
    "chsh",
    "hom_visibility",
    "load_setup",
    "read_tags",
    "write_tags",
]
