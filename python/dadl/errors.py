"""Typed errors. The same types are raised by the library and, later, over MCP."""


class DadlError(Exception):
    """Base class for all DADL errors."""


class GuardrailRefusal(DadlError):
    """A correction, window or analysis choice is invalid for the analysis type."""


class InsufficientCounts(DadlError):
    """Counts are below the estimator's validity threshold."""


class ConfigError(DadlError):
    """A setup or policy file is missing, malformed or of an unsupported version."""


class AdapterError(DadlError):
    """A vendor file or device stream is malformed, truncated or of an unsupported format."""


class StaleCalibration(DadlError):
    """A calibration artifact does not cover the dataset: outside its validity time range, or
    recorded under conditions (tagger settings, temperature, bias) that do not match it."""
