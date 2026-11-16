# Inference library benchmark (phase 3, M0)

Measures NumPyro, PyMC and Stan (CmdStanPy) on the same two models, to choose the optional MCMC
backend for phase 3 (decision D2 in [docs/phase3-plan.md](../../docs/phase3-plan.md)). Not
shipped in any package. Results: [docs/reviews/2026-11-ppl-benchmark.md](../../docs/reviews/2026-11-ppl-benchmark.md).

## Models

- `dark`: DADL's dark inter-arrival likelihood (`dadl.models.DarkAfterpulseModel`), about
  20,000 gaps from `dadl.sim.simulate_dark` (1,000 cps, afterpulse p = 0.05, tau = 50 ns).
  Reference: DADL's own `calibrate_dark` profile intervals.
- `delay`: Poisson counts in 400 bins of a Gaussian peak on a flat background, the delay
  calibration's shape. It stands in for the phase 3 HOM joint model, which does not exist yet.
  Reference: maximum likelihood with a Laplace interval.

Priors are Normal(MLE, 10) in the sampled coordinates. Every chain starts at the MLE, as DADL
would: the delay model has a second, broad-peak mode, and the libraries' default random
initialisation put some NumPyro chains there (R-hat 1.55 in a first try).

Each library runs 4 chains of 1,000 warmup and 1,000 draws, twice. The first run includes
compilation; ESS per second uses the second.

## Files

| File | Runs in | Does |
| --- | --- | --- |
| `prepare.py` | an environment with dadl | writes `data/` (tags, histogram, references) |
| `install.py` | any Python | fresh venv, `pip install --only-binary=:all:` (no building from source), size; builds CmdStan for Stan |
| `run_numpyro.py`, `run_pymc.py`, `run_stan.py` | the clean library venv | samples, writes draws and timings |
| `analyze.py` | the dadl environment plus `arviz<1` | R-hat, bulk and tail ESS, quantile shift against the reference |
| `report.py` | any Python | one Markdown table of all runs |

Quantile shift: the largest distance of the posterior 2.5%, 50% and 97.5% quantiles from the
reference interval's ends and centre, in units of its half-width. Under 0.2 is agreement.

## Running

CI: the `PPL benchmark` workflow (manual, and on pull requests that touch `bench/ppl/`) runs
every library on Windows, Linux and macOS with Python 3.10 and 3.13, plus PyMC with its C
backend off (`PYTENSOR_FLAGS=cxx=`), and posts the table as the job summary.

Locally:

```bash
python bench/ppl/prepare.py /tmp/ppl/data
python bench/ppl/install.py numpyro /tmp/ppl/venv /tmp/ppl/run
cd bench/ppl && /tmp/ppl/venv/bin/python run_numpyro.py /tmp/ppl/data /tmp/ppl/run
```
