# dadl-quantum-core

The compiled per-tag core of [DADL](https://github.com/usecommon/dadl): dead-time
filters, fixed-grid coincidences, inter-arrival and cross-correlation
histograms, from the Rust crate `dadl`, for NumPy int64 arrays.

You do not import it yourself. Install it next to `dadl-quantum`:

```bash
pip install dadl-quantum dadl-quantum-core
```

and `dadl` uses it automatically (`dadl.backend()` returns `"rust"`). Set
`DADL_BACKEND=python` to force the pure-Python reference, or
`DADL_BACKEND=rust` to fail loudly when the core is missing. Both backends are
held to the same shared test vectors.

Dual-licensed MIT OR Apache-2.0.
