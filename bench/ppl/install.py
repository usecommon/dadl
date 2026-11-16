"""Install one library into a fresh venv and record what it took.

Usage: python install.py LIBRARY VENV_DIR RUN_DIR

Installs with `pip install --only-binary=:all:`, so a library that would need a compiler to
build from source fails here instead of quietly compiling: that is the "clean install on a lab
PC" test. For `cmdstanpy`, CmdStan itself is then built with `install_cmdstan` (on Windows with
`--compiler`, which also downloads a toolchain), and that is timed separately. Writes
RUN_DIR/install.json; exits 0 even on failure so the record is kept.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import venv
from pathlib import Path

# Pinned so the install does not query the GitHub API for the latest release (rate-limited on CI).
CMDSTAN_VERSION = "2.40.0"
PACKAGES = {"numpyro": ["numpyro"], "pymc": ["pymc"], "cmdstanpy": ["cmdstanpy"]}


def size_mb(path: Path) -> float:
    total = 0
    for root, _dirs, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total / 1e6


def main(library: str, venv_dir: str, run_dir: str) -> None:
    env_dir = Path(venv_dir)
    out = Path(run_dir)
    out.mkdir(parents=True, exist_ok=True)
    venv.create(env_dir, with_pip=True, clear=True)
    py = env_dir / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    subprocess.run([str(py), "-m", "pip", "install", "-q", "--upgrade", "pip"], check=False)
    base = size_mb(env_dir)
    rec: dict = {"library": library, "pip_only_binary": True}
    t = time.perf_counter()
    proc = subprocess.run([str(py), "-m", "pip", "install", "--only-binary=:all:", *PACKAGES[library]], capture_output=True, text=True)
    rec["pip_s"] = time.perf_counter() - t
    rec["pip_ok"] = proc.returncode == 0
    rec["installed_mb"] = size_mb(env_dir) - base
    if not rec["pip_ok"]:
        rec["pip_error"] = (proc.stderr or proc.stdout)[-2000:]
    if library == "cmdstanpy" and rec["pip_ok"]:
        home = (out / "cmdstan-home").resolve()
        home.mkdir(exist_ok=True)
        args = [str(py), "-m", "cmdstanpy.install_cmdstan", "--dir", str(home), "--cores", "4", "--version", CMDSTAN_VERSION]
        if os.name == "nt":
            args.append("--compiler")
        t = time.perf_counter()
        proc = subprocess.run(args, capture_output=True, text=True)
        rec["cmdstan_build_s"] = time.perf_counter() - t
        rec["cmdstan_ok"] = proc.returncode == 0
        rec["cmdstan_mb"] = size_mb(home)
        if not rec["cmdstan_ok"]:
            rec["cmdstan_error"] = (proc.stderr or proc.stdout)[-2000:]
        else:
            rec["cmdstan_dir"] = str(next(home.glob("cmdstan-*")))
    rec["python"] = sys.version.split()[0]
    (out / "install.json").write_text(json.dumps(rec, indent=2))
    print(json.dumps({k: v for k, v in rec.items() if not k.endswith("_error")}, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:4])
