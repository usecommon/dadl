"""Diagnostics for one runner's output. Usage: python analyze.py DATA_DIR RUN_DIR

Runs in its own environment (NumPy and ArviZ), so the measured library's install stays clean.
Adds to RUN_DIR/result.json: per model the worst split R-hat, the smallest bulk and tail ESS,
bulk ESS per second of the second run, and how far the posterior 2.5%, 50% and 97.5% quantiles
sit from the reference interval, in units of the reference half-width. Install figures, if the
workflow wrote RUN_DIR/install.json, are merged in.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import arviz as az
import numpy as np


def main(data_dir: str, run_dir: str) -> None:
    run = Path(run_dir)
    meta = json.loads((Path(data_dir) / "meta.json").read_text())
    record = json.loads((run / "timings.json").read_text())
    install = run / "install.json"
    if install.exists():
        record["install"] = json.loads(install.read_text())
    for model, info in record["models"].items():
        draws = dict(np.load(run / f"{model}_draws.npz"))
        idata = az.from_dict(posterior=draws)
        rhat = az.rhat(idata)
        bulk = az.ess(idata, method="bulk")
        tail = az.ess(idata, method="tail")
        ref = meta[model]["reference"]
        shift = {}
        for name, (lo, mid, hi) in ref.items():
            q = np.quantile(draws[name].ravel(), [0.025, 0.5, 0.975])
            half = 0.5 * (hi - lo)
            shift[name] = float(max(abs(q[0] - lo), abs(q[1] - mid), abs(q[2] - hi)) / half)
        info["rhat_max"] = float(max(float(rhat[v]) for v in draws))
        info["ess_bulk_min"] = float(min(float(bulk[v]) for v in draws))
        info["ess_tail_min"] = float(min(float(tail[v]) for v in draws))
        info["ess_bulk_per_s"] = info["ess_bulk_min"] / info["second_run_s"]
        info["quantile_shift_max"] = max(shift.values())
        info["quantile_shift"] = shift
    (run / "result.json").write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
