"""CmdStanPy runner. Usage: python run_stan.py DATA_DIR OUT_DIR

Needs CmdStan installed (`python -m cmdstanpy.install_cmdstan`, which needs make and a C++
compiler; on Windows `--compiler` also installs a toolchain). The model compile time is recorded
separately from sampling.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

T0 = time.perf_counter()
import cmdstanpy  # noqa: E402

IMPORT_S = time.perf_counter() - T0

from common import CHAINS, DARK_SITES, DELAY_SITES, DRAWS, WARMUP, load, save  # noqa: E402

HERE = Path(__file__).parent / "models"


def run(stan_file, data, names, seed, init):
    t = time.perf_counter()
    model = cmdstanpy.CmdStanModel(stan_file=str(stan_file), force_compile=True)
    compile_s = time.perf_counter() - t

    def once(sd):
        t = time.perf_counter()
        fit = model.sample(data=data, chains=CHAINS, parallel_chains=CHAINS, iter_warmup=WARMUP, iter_sampling=DRAWS, seed=sd, inits=init, show_progress=False)
        return fit, time.perf_counter() - t

    _, first = once(seed)
    fit, second = once(seed + 1)
    draws = {}
    for n in names:
        arr = fit.stan_variable(n)  # (chains * draws,)
        draws[n] = arr.reshape(CHAINS, DRAWS)
    return draws, compile_s + first, second, compile_s


def main(data_dir: str, out_dir: str) -> None:
    meta, dark, delay = load(data_dir)
    results, compile_times = {}, {}
    data = {"N": int(dark["s_us"].size), "s": dark["s_us"].tolist(), "centre": meta["dark"]["prior_centre"]}
    d, a, b, c = run(HERE / "dark.stan", data, ("r_per_us", "p", "tau_us"), 1, dict(zip(DARK_SITES, meta["dark"]["prior_centre"])))
    results["dark"], compile_times["dark"] = (a, b, d), c
    data = {
        "K": int(delay["counts"].size),
        "centres": delay["centres"].tolist(),
        "width": float(delay["width"]),
        "counts": [int(v) for v in delay["counts"]],
        "centre": meta["delay"]["prior_centre"],
    }
    d, a, b, c = run(HERE / "delay.stan", data, ("n", "mu", "sigma", "b"), 2, dict(zip(DELAY_SITES, meta["delay"]["prior_centre"])))
    results["delay"], compile_times["delay"] = (a, b, d), c
    save(out_dir, "cmdstanpy", f"{cmdstanpy.__version__} (CmdStan {cmdstanpy.cmdstan_version()})", IMPORT_S, results, {"compile_s": compile_times})


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
