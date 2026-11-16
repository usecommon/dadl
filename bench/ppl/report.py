"""Collect every run's result.json (or install.json when the run failed) into one Markdown table.

Usage: python report.py RESULTS_DIR > table.md
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def fmt(v, nd=1):
    return "" if v is None else f"{v:.{nd}f}"


def main(root: str) -> None:
    rows = []
    for d in sorted(Path(root).glob("*/")):
        res, inst = d / "result.json", d / "install.json"
        rec = json.loads(res.read_text()) if res.exists() else {"install": json.loads(inst.read_text()) if inst.exists() else {}}
        label = json.loads((d / "label.json").read_text()) if (d / "label.json").exists() else {"name": d.name}
        rows.append((label, rec))
    print("| Run | Install (s) | Size (MB) | Import (s) | Model | First run (s) | ESS/s | R-hat | Quantile shift | Status |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for label, rec in rows:
        inst = rec.get("install", {})
        install_s = inst.get("pip_s")
        size = inst.get("installed_mb")
        if "cmdstan_build_s" in inst:
            install_s = (install_s or 0) + inst["cmdstan_build_s"]
            size = (size or 0) + inst.get("cmdstan_mb", 0)
        status = "ok"
        if inst and not inst.get("pip_ok", True):
            status = "pip install failed (no wheel)"
        elif inst and inst.get("cmdstan_ok") is False:
            status = "CmdStan build failed"
        elif "models" not in rec:
            status = "run failed"
        models = rec.get("models") or {"-": {}}
        for model, m in models.items():
            bad = m and (m.get("rhat_max", 1) > 1.01 or m.get("quantile_shift_max", 0) > 0.2)
            print(
                f"| {label['name']} | {fmt(install_s)} | {fmt(size, 0)} | {fmt(rec.get('import_s'))} | {model} "
                f"| {fmt(m.get('first_run_s'))} | {fmt(m.get('ess_bulk_per_s'), 0)} | {fmt(m.get('rhat_max'), 3)} "
                f"| {fmt(m.get('quantile_shift_max'), 2)} | {'check' if bad else status} |"
            )


if __name__ == "__main__":
    main(sys.argv[1])
