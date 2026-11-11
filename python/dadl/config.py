"""Loading setup.yaml files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from dadl.errors import ConfigError

SUPPORTED_SCHEMA_VERSIONS = {1}
REQUIRED_KEYS = ("schema_version", "name", "source", "detectors", "analysis")


def load_setup(path: str | Path) -> dict[str, Any]:
    """Load and minimally validate a setup file.

    Full JSON Schema validation (dadl check-config) and migrations
    (dadl migrate) arrive in a later phase.
    """
    path = Path(path)
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ConfigError(f"setup file not found: {path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"setup file is not valid YAML: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError("setup file must be a mapping at the top level")
    missing = [k for k in REQUIRED_KEYS if k not in data]
    if missing:
        raise ConfigError(f"setup file is missing keys: {', '.join(missing)}")
    if data["schema_version"] not in SUPPORTED_SCHEMA_VERSIONS:
        raise ConfigError(
            f"schema_version {data['schema_version']!r} is not supported; "
            f"supported: {sorted(SUPPORTED_SCHEMA_VERSIONS)}"
        )
    window = data["analysis"].get("window", {})
    if window.get("source") != "plan":
        raise ConfigError("analysis.window.source must be 'plan': windows are declared before acquisition")
    return data
