"""Write the benchmark data and the reference intervals the samplers are compared with.

Run in an environment with dadl installed (`pip install -e python`):

    python bench/ppl/prepare.py OUT_DIR

Two models, the same for every library:

- `dark`: DADL's dark inter-arrival likelihood (`dadl.models.DarkAfterpulseModel`) on about
  20,000 gaps from `dadl.sim.simulate_dark`, past the dead time that `calibrate_dark` finds.
  Parameters log r, logit p, log tau_a, with times in microseconds. The reference is DADL's own
  `calibrate_dark` (profile-likelihood intervals).
- `delay`: a Poisson histogram of a Gaussian peak on a flat background, the shape of the delay
  calibration's cross-correlation. It stands in for the phase 3 HOM joint model, which does not
  exist yet. Parameters log N, mu, log sigma, log b (ns). The reference is the maximum-likelihood
  fit with a Laplace interval.

Priors in every library: Normal(MLE, 10) on each parameter in these coordinates, wide enough
that the data dominate; the MLE only centres them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import optimize, stats

from dadl.calibration.dark import calibrate_dark
from dadl.data import ChannelId
from dadl.sim import DetectorParams, simulate_dark

SEED = 20261116
PS_PER_US = 1e6

DARK_TRUTH = DetectorParams(dark_rate_cps=1000.0, dead_time_ps=22_000, afterpulse_prob=0.05, afterpulse_tau_ps=50_000.0)
DARK_SECONDS = 20.0

DELAY_TRUTH = {"n": 5000.0, "mu": 0.4, "sigma": 0.3, "b": 20.0}  # counts, ns, ns, counts per bin
DELAY_EDGES = np.linspace(-10.0, 10.0, 401)  # 50 ps bins


def dark(out: Path, rng: np.random.Generator) -> dict:
    tags = simulate_dark(DARK_TRUTH, DARK_SECONDS, rng)
    cal = calibrate_dark(tags, ChannelId("bench", 0), backend="sim", seed=SEED)
    dead = cal.estimates["dead_time_ps"].value
    s_us = (np.diff(tags) - dead) / PS_PER_US
    s_us = s_us[s_us >= 0]
    est = cal.estimates
    r_us = est["dark_rate_cps"].value / 1e6
    p = est["afterpulse_prob"].value
    tau_us = est["afterpulse_tau_ps"].value / PS_PER_US
    np.savez(out / "dark.npz", s_us=s_us)
    ref = {
        # natural parameters, as the runners report them
        "r_per_us": [est["dark_rate_cps"].low / 1e6, r_us, est["dark_rate_cps"].high / 1e6],
        "p": [est["afterpulse_prob"].low, p, est["afterpulse_prob"].high],
        "tau_us": [est["afterpulse_tau_ps"].low / PS_PER_US, tau_us, est["afterpulse_tau_ps"].high / PS_PER_US],
    }
    centre = [float(np.log(r_us)), float(np.log(p / (1 - p))), float(np.log(tau_us))]
    return {"n_gaps": int(s_us.size), "prior_centre": centre, "reference": ref, "reference_method": "calibrate_dark profile intervals"}


def _delay_nll(theta: np.ndarray, centres: np.ndarray, width: float, counts: np.ndarray) -> float:
    log_n, mu, log_sigma, log_b = theta
    sigma = np.exp(log_sigma)
    lam = np.exp(log_b) + np.exp(log_n) * width * stats.norm.pdf(centres, mu, sigma)
    return float(np.sum(lam - counts * np.log(lam)))


def delay(out: Path, rng: np.random.Generator) -> dict:
    centres = 0.5 * (DELAY_EDGES[1:] + DELAY_EDGES[:-1])
    width = float(DELAY_EDGES[1] - DELAY_EDGES[0])
    t = DELAY_TRUTH
    lam = t["b"] + t["n"] * width * stats.norm.pdf(centres, t["mu"], t["sigma"])
    counts = rng.poisson(lam)
    np.savez(out / "delay.npz", centres=centres, width=np.array(width), counts=counts)
    x0 = np.array([np.log(t["n"]) + 0.2, 0.0, np.log(0.5), np.log(max(counts.min(), 1))])
    fit = optimize.minimize(_delay_nll, x0, args=(centres, width, counts), method="BFGS")
    # Hessian by central differences, Laplace interval on each coordinate
    h = 1e-4
    k = fit.x.size
    hess = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            def f(di, dj):
                x = fit.x.copy()
                x[i] += di
                x[j] += dj
                return _delay_nll(x, centres, width, counts)
            hess[i, j] = (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)
    se = np.sqrt(np.diag(np.linalg.inv(hess)))
    z = stats.norm.ppf(0.975)
    names = ("n", "mu", "sigma", "b")
    natural = (np.exp, lambda v: v, np.exp, np.exp)
    ref = {nm: [float(fn(fit.x[i] - z * se[i])), float(fn(fit.x[i])), float(fn(fit.x[i] + z * se[i]))] for i, (nm, fn) in enumerate(zip(names, natural))}
    return {"n_bins": int(counts.size), "prior_centre": [float(v) for v in fit.x], "reference": ref, "reference_method": "MLE with Laplace interval"}


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    meta = {"seed": SEED, "dark": dark(out, rng), "delay": delay(out, rng)}
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps({k: v for k, v in meta.items() if k != "seed"}, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
