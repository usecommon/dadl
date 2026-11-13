"""A calibration method's code may not change without a new method version.

Archivist review, cycle 8: methods changed across several commits under one version label, so
artifacts with different statistics were indistinguishable. Each calibration module's source hash
is pinned here with its METHOD_VERSION. After changing a module, bump its METHOD_VERSION and run
`python python/tests/test_method_versions.py` to update the table.
"""

import hashlib
import importlib
import json
from pathlib import Path

TABLE = Path(__file__).with_name("method_versions.json")
MODULES = ("accidentals", "dark", "dcr", "delay", "heralding")


def current() -> dict[str, dict[str, str]]:
    out = {}
    for name in MODULES:
        mod = importlib.import_module(f"dadl.calibration.{name}")
        text = Path(mod.__file__).read_text(encoding="utf-8").replace("\r\n", "\n")
        out[name] = {"version": mod.METHOD_VERSION, "sha256": hashlib.sha256(text.encode()).hexdigest()}
    return out


def test_method_code_changes_come_with_a_new_version():
    pinned = json.loads(TABLE.read_text(encoding="utf-8"))
    for name, now in current().items():
        was = pinned[name]
        if now["sha256"] != was["sha256"]:
            assert now["version"] != was["version"], (
                f"dadl.calibration.{name} changed but METHOD_VERSION is still {now['version']}: bump it, "
                "then regenerate python/tests/method_versions.json"
            )
            raise AssertionError(f"regenerate python/tests/method_versions.json for {name} (version {now['version']})")
        assert now["version"] == was["version"], f"{name}: version changed; regenerate the table"


if __name__ == "__main__":
    TABLE.write_text(json.dumps(current(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {TABLE}")
