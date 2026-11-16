"""NumPyro (JAX, CPU) runner. Usage: python run_numpyro.py DATA_DIR OUT_DIR"""

from __future__ import annotations

import sys
import time

T0 = time.perf_counter()
import numpyro  # noqa: E402

numpyro.set_host_device_count(4)
numpyro.enable_x64()
import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import numpyro.distributions as dist  # noqa: E402
from numpyro.infer import MCMC, NUTS, init_to_value  # noqa: E402

IMPORT_S = time.perf_counter() - T0

from common import CHAINS, DARK_SITES, DELAY_SITES, DRAWS, WARMUP, load, save  # noqa: E402


def dark_model(s, centre):
    x = numpyro.sample("log_r", dist.Normal(centre[0], 10.0))
    y = numpyro.sample("logit_p", dist.Normal(centre[1], 10.0))
    z = numpyro.sample("log_tau", dist.Normal(centre[2], 10.0))
    r, p, ta = jnp.exp(x), jax.nn.sigmoid(y), jnp.exp(z)
    d = r * (1 - p) + p * (r + 1 / ta) * jnp.exp(-s / ta)
    numpyro.factor("ll", jnp.sum(jnp.log(d)) - r * jnp.sum(s))
    numpyro.deterministic("r_per_us", r)
    numpyro.deterministic("p", p)
    numpyro.deterministic("tau_us", ta)


def delay_model(centres, width, counts, centre):
    log_n = numpyro.sample("log_n", dist.Normal(centre[0], 10.0))
    mu = numpyro.sample("mu", dist.Normal(centre[1], 10.0))
    log_sigma = numpyro.sample("log_sigma", dist.Normal(centre[2], 10.0))
    log_b = numpyro.sample("log_b", dist.Normal(centre[3], 10.0))
    sigma = jnp.exp(log_sigma)
    lam = jnp.exp(log_b) + jnp.exp(log_n) * width * jnp.exp(-0.5 * ((centres - mu) / sigma) ** 2) / (sigma * jnp.sqrt(2 * jnp.pi))
    numpyro.sample("obs", dist.Poisson(lam), obs=counts)
    numpyro.deterministic("n", jnp.exp(log_n))
    numpyro.deterministic("sigma", sigma)
    numpyro.deterministic("b", jnp.exp(log_b))


def run(model, kwargs, names, seed, init):
    mcmc = MCMC(NUTS(model, init_strategy=init_to_value(values=init)), num_warmup=WARMUP, num_samples=DRAWS, num_chains=CHAINS, chain_method="parallel", progress_bar=False)
    t = time.perf_counter()
    mcmc.run(jax.random.PRNGKey(seed), **kwargs)
    first = time.perf_counter() - t
    t = time.perf_counter()
    mcmc.run(jax.random.PRNGKey(seed + 1), **kwargs)
    second = time.perf_counter() - t
    samples = mcmc.get_samples(group_by_chain=True)
    return {n: jax.device_get(samples[n]) for n in names}, first, second


def main(data_dir: str, out_dir: str) -> None:
    meta, dark, delay = load(data_dir)
    s = jnp.asarray(dark["s_us"])
    timings = {}
    draws, a, b = run(dark_model, {"s": s, "centre": meta["dark"]["prior_centre"]}, ("r_per_us", "p", "tau_us"), 1, dict(zip(DARK_SITES, meta["dark"]["prior_centre"])))
    timings["dark"] = (a, b, draws)
    kw = {"centres": jnp.asarray(delay["centres"]), "width": float(delay["width"]), "counts": jnp.asarray(delay["counts"]), "centre": meta["delay"]["prior_centre"]}
    draws, a, b = run(delay_model, kw, ("n", "mu", "sigma", "b"), 2, dict(zip(DELAY_SITES, meta["delay"]["prior_centre"])))
    timings["delay"] = (a, b, draws)
    save(out_dir, "numpyro", numpyro.__version__, IMPORT_S, timings)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
