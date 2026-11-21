# Inference library benchmark (phase 3, M0)

Nov 16, 2026 · decision D2 in [phase3-plan.md](../phase3-plan.md) · harness in
[bench/ppl/](../../bench/ppl/README.md) · CI run 43965078466 (`PPL benchmark` workflow)

## Decision

**NumPyro**, as the optional extra `dadl-quantum[mcmc]`, confirming the preference taken going in.
Laplace stays the default and the base install keeps NumPy and SciPy only. PyMC is the fallback
if JAX stops shipping Windows wheels. Stan is not adopted.

Why, in one line each:

- **Clean install on a lab PC:** NumPyro (with JAX) and PyMC install from binary wheels on
  Windows, Linux and macOS, with no compiler, in under 35 s. Stan needs CmdStan built from source:
  2–3 minutes, 0.5–2.6 GB, and a C++ toolchain (on Windows, one it downloads itself).
- **Speed where it matters:** NumPyro's first run, compilation included, is 1–3 s on every
  platform. PyMC's is 4–38 s, slowest on Windows; Stan adds 5–12 s of model compilation per model
  and was the slowest on the 20,000-gap dark likelihood (150–300 ESS/s against 770–1,850 for
  NumPyro). In steady state Stan was fastest on the delay model (1.6–3x NumPyro on every OS) and
  the libraries overlap on macOS; NumPyro was chosen on install and first run (council round 6).
- **Correctness is a tie:** all 42 runs converged (R-hat ≤ 1.004) and every posterior's 2.5%,
  50% and 97.5% quantiles were within 0.07 reference half-widths of DADL's own intervals.

## Setup

Two models, the same for every library (details in the harness README):

- `dark`: DADL's dark inter-arrival likelihood on 20,881 simulated gaps (1,000 cps, afterpulse
  p = 0.05, tau = 50 ns), compared with `calibrate_dark`'s profile intervals.
- `delay`: Poisson counts in 400 bins of a Gaussian peak on a flat background, compared with the
  maximum-likelihood Laplace interval. It stands in for the M2 HOM joint model, which does not
  exist yet; M2 should rerun this benchmark on the real model.

4 chains × (1,000 warmup + 1,000 draws), every chain started at the MLE, priors Normal(MLE, 10).
Each library installed into a fresh venv with `pip install --only-binary=:all:`. GitHub-hosted
runners (windows-latest, ubuntu-latest, macos-latest), Python 3.10 and 3.13, plus PyMC with its C
backend off (`PYTENSOR_FLAGS=cxx=`) on 3.13. ESS/s is the smallest bulk ESS over the parameters,
per second of a second run (compilation cached).

## Results

Ranges over the two Python versions. "First run" includes compilation (for Stan, model
compilation is in addition).

| | NumPyro | PyMC | PyMC, no compiler | Stan (CmdStanPy) |
| --- | --- | --- | --- | --- |
| Versions | 0.22.0 (0.19.0 on 3.10) | 6.4.0 (5.25.1 on 3.10) | 6.4.0 | 1.3.2, CmdStan 2.40.0 |
| Binary-only install | yes, all | yes, all | yes, all | pip yes; CmdStan built from source |
| Install time, s (Win / Linux / macOS) | 15–17 / 10–12 / 9–16 | 25–28 / 19 / 13–34 | 33 / 22 / 24 | 187–206 / 109–126 / 140–150 |
| Size, MB (Win / Linux / macOS) | 426–474 / 559–585 / 387–415 | 353–477 / 411–559 / 293–431 | 477 / 559 / 431 | 2,754–2,772 / 1,291–1,304 / 578–602 |
| Import, s | 0.9–2.3 | 1.0–26 (first import on 3.10 compiles) | 1.1–1.3 | 0.2–1.0 |
| Model compile, s | (in first run) | (in first run) | (in first run) | 5–12 per model |
| First run, dark, s (Win / Linux / macOS) | 1.4–1.6 / 1.1–1.6 / 1.9–2.3 | 21–38 / 9–21 / 4–24 | 20 / 11 / 8 | 33 / 19–24 / 23–26 |
| ESS/s, dark (Win / Linux / macOS) | 832–1,275 / 868–1,846 / 772–1,425 | 308–519 / 666–1,140 / 834–3,045 | 376 / 610 / 1,562 | 148–183 / 261–302 / 198–203 |
| ESS/s, delay (Win / Linux / macOS) | 2,033–2,382 / 2,454–2,616 / 1,569–3,002 | 602–962 / 2,036–2,885 / 2,810–9,041 | 695 / 2,264 / 3,349 | 3,758–3,897 / 6,311–8,406 / 4,685–5,766 |
| Worst R-hat | 1.003 | 1.004 | 1.004 | 1.003 |
| Worst quantile shift (half-widths) | 0.07 | 0.06 | 0.06 | 0.04 |
| Licence | Apache-2.0 (JAX Apache-2.0) | Apache-2.0 (PyTensor BSD-3-Clause) | same | BSD-3-Clause |

All licences are compatible with DADL's MIT OR Apache-2.0.

## Findings

1. **Stan on Windows breaks with a second MinGW on PATH.** CmdStan built with the toolchain
   `install_cmdstan --compiler` downloads, but models were then compiled with the runner's own
   GCC 15 (`C:\mingw64`, first on PATH) and failed to link (undefined `std::` references). It
   worked once `cmdstanpy.utils.cxx_toolchain_path` put CmdStan's toolchain first. Lab PCs often
   carry another MinGW (Git for Windows, Strawberry Perl, vendor SDKs), so this would recur.
2. **CmdStan's version lookup needs the GitHub API.** On macOS runners `install_cmdstan` failed
   with "Cannot connect to CmdStan github repo" until the version was pinned.
3. **Default initialisation found a second mode.** With the libraries' default random starting
   points, some NumPyro chains settled in a broad-peak mode of the delay model (R-hat 1.55). Every
   library was then started at the MLE, as DADL will do. The phase 3 sampler must start from the
   maximum-likelihood fit and check R-hat across chains (M5 already plans this).
4. **PyMC without a compiler is usable.** PyMC 6 with the C backend off ran at 37–78% of its
   compiled speed (ESS/s, Python 3.13), so a compiler is not a hard requirement for PyMC either.
5. **Python 3.10 gets older releases.** NumPyro 0.22 and PyMC 6 need Python 3.11+, so 3.10
   resolves to NumPyro 0.19 and PyMC 5.25. The `[mcmc]` extra should set a floor that works on
   3.10 (NumPyro ≥ 0.19) or drop 3.10 for it.
6. **JAX is large.** NumPyro with JAX is 390–590 MB, against about 130 MB for DADL with SciPy.
   That is acceptable for an optional extra, not for the base install.

## Limits

- GitHub-hosted runners, not lab PCs; one run per configuration, so the speeds are indicative
  (macOS ESS/s varied about 2× between reruns).
- Two small models. The M2 HOM joint model should be benchmarked again before M5 lands.
- CPU only; GPU JAX was not measured and is not needed for these model sizes.
- Stan's dark-model speed includes CmdStanPy writing the 20,000 gaps to a JSON file on each run,
  which is how it would be used.
