"""PyMC runner. Usage: python run_pymc.py DATA_DIR OUT_DIR

Set PYTENSOR_FLAGS=cxx= to measure PyMC without a C compiler (pure Python and NumPy backend).
"""

from __future__ import annotations

import os
import sys
import time

T0 = time.perf_counter()
import pymc as pm  # noqa: E402
import pytensor  # noqa: E402
import pytensor.tensor as pt  # noqa: E402

IMPORT_S = time.perf_counter() - T0

from common import CHAINS, DARK_SITES, DELAY_SITES, DRAWS, WARMUP, load, save  # noqa: E402


def dark_model(s, centre):
    with pm.Model() as m:
        x = pm.Normal("log_r", centre[0], 10.0)
        y = pm.Normal("logit_p", centre[1], 10.0)
        z = pm.Normal("log_tau", centre[2], 10.0)
        r, p, ta = pt.exp(x), pm.math.sigmoid(y), pt.exp(z)
        d = r * (1 - p) + p * (r + 1 / ta) * pt.exp(-s / ta)
        pm.Potential("ll", pt.sum(pt.log(d)) - r * s.sum())
        pm.Deterministic("r_per_us", r)
        pm.Deterministic("p", p)
        pm.Deterministic("tau_us", ta)
    return m


def delay_model(centres, width, counts, centre):
    with pm.Model() as m:
        log_n = pm.Normal("log_n", centre[0], 10.0)
        mu = pm.Normal("mu", centre[1], 10.0)
        log_sigma = pm.Normal("log_sigma", centre[2], 10.0)
        log_b = pm.Normal("log_b", centre[3], 10.0)
        sigma = pt.exp(log_sigma)
        lam = pt.exp(log_b) + pt.exp(log_n) * width * pt.exp(-0.5 * ((centres - mu) / sigma) ** 2) / (sigma * pt.sqrt(2 * pt.pi))
        pm.Poisson("obs", lam, observed=counts)
        pm.Deterministic("n", pt.exp(log_n))
        pm.Deterministic("sigma", sigma)
        pm.Deterministic("b", pt.exp(log_b))
    return m


def run(model, names, seed, init):
    def once(sd):
        t = time.perf_counter()
        with model:
            idata = pm.sample(DRAWS, tune=WARMUP, chains=CHAINS, cores=CHAINS, random_seed=sd, initvals=init, progressbar=False, compute_convergence_checks=False)
        return idata, time.perf_counter() - t

    _, first = once(seed)
    idata, second = once(seed + 1)
    return {n: idata.posterior[n].values for n in names}, first, second


def main(data_dir: str, out_dir: str) -> None:
    meta, dark, delay = load(data_dir)
    results = {}
    results["dark"] = run(dark_model(dark["s_us"], meta["dark"]["prior_centre"]), ("r_per_us", "p", "tau_us"), 1, dict(zip(DARK_SITES, meta["dark"]["prior_centre"])))
    m = delay_model(delay["centres"], float(delay["width"]), delay["counts"], meta["delay"]["prior_centre"])
    results["delay"] = run(m, ("n", "mu", "sigma", "b"), 2, dict(zip(DELAY_SITES, meta["delay"]["prior_centre"])))
    results = {k: (a, b, d) for k, (d, a, b) in results.items()}
    cxx = pytensor.config.cxx
    name = "pymc" if cxx else "pymc-nocompiler"
    save(out_dir, name, pm.__version__, IMPORT_S, results, {"pytensor_cxx": cxx or "", "pytensor_flags": os.environ.get("PYTENSOR_FLAGS", "")})


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
