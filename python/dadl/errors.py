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
