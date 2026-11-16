"""Shared settings and I/O for the sampler runners (NumPy and the standard library only)."""

from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np

CHAINS = 4
WARMUP = 1000
DRAWS = 1000

# Unconstrained parameter names; every library starts all chains at the MLE (meta prior_centre),
# as DADL would, since the delay model has a second, broad-peak mode far from it.
DARK_SITES = ("log_r", "logit_p", "log_tau")
DELAY_SITES = ("log_n", "mu", "log_sigma", "log_b")


def load(data_dir: str):
    d = Path(data_dir)
    meta = json.loads((d / "meta.json").read_text())
    return meta, dict(np.load(d / "dark.npz")), dict(np.load(d / "delay.npz"))


def save(out_dir: str, library: str, version: str, import_s: float, results: dict, extra: dict | None = None) -> None:
    """results: model -> (first run s, second run s, {param: array (chains, draws)})."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    record = {
        "library": library,
        "version": version,
        "python": platform.python_version(),
        "os": sys.platform,
        "import_s": import_s,
        "chains": CHAINS,
        "warmup": WARMUP,
        "draws": DRAWS,
        "models": {},
        **(extra or {}),
    }
    for model, (first, second, draws) in results.items():
        np.savez(out / f"{model}_draws.npz", **{k: np.asarray(v, dtype=float) for k, v in draws.items()})
        record["models"][model] = {"first_run_s": first, "second_run_s": second}
    (out / "timings.json").write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))
