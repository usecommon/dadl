"""Which implementation runs the per-tag functions: the Rust core or the Python reference.

`dadl-quantum-core` (import name `dadl_core`) is an optional compiled
package. When it is installed it is used; otherwise the pure-Python
reference runs. Both are held to the shared vectors in tests/vectors/core.json.

    DADL_BACKEND=python   always use the reference
    DADL_BACKEND=rust     require the compiled core (ImportError if missing)
    (unset)               the core if installed, else the reference

Argument types are checked in Python for both backends; sorting, range and
overflow checks run in whichever backend does the work, and raise the same
exception types.
"""

from __future__ import annotations

import contextlib
import os
from collections.abc import Iterator
from types import ModuleType

_core: ModuleType | None = None
_forced = os.environ.get("DADL_BACKEND", "").strip().lower()
if _forced not in ("", "python", "rust"):
    raise ImportError(f"DADL_BACKEND must be 'python' or 'rust', got {_forced!r}")
if _forced != "python":
    try:
        import dadl_core as _loaded  # type: ignore[import-not-found]
    except ImportError:
        if _forced == "rust":
            raise ImportError(
                "DADL_BACKEND=rust but the compiled core is not installed: pip install dadl-quantum-core"
            ) from None
        _loaded = None
    _core = _loaded
_installed = _core


def core() -> ModuleType | None:
    """The compiled core module in use, or None for the Python reference."""
    return _core


def name() -> str:
    """"rust" or "python"."""
    return "rust" if _core is not None else "python"


def available() -> tuple[str, ...]:
    return ("python", "rust") if _installed is not None else ("python",)


@contextlib.contextmanager
def use(backend: str) -> Iterator[None]:
    """Temporarily switch backend (for tests and benchmarks)."""
    global _core
    if backend not in ("python", "rust"):
        raise ValueError("backend must be 'python' or 'rust'")
    if backend == "rust" and _installed is None:
        raise ImportError("the compiled core (dadl-quantum-core) is not installed")
    previous = _core
    _core = _installed if backend == "rust" else None
    try:
        yield
    finally:
        _core = previous
