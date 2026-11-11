# Adapters

Vendor adapters arrive in phase 2, starting with one tagger (Swabian) and an
adapter template. PicoQuant, qutools and IDQ follow as "v0.2+ / seeking
contributors".

**Until then:** export your tags to the generic format that `dadl analyze`
reads: `.npz` or `.csv` with integer columns `setting`, `channel`,
`timestamp_ps`. See "Using your own data" in the repository README and
`python/dadl/tagfile.py`.

Planned design:

- In-tree adapters implement the Rust `Adapter` trait.
- Python-only vendor SDKs implement the equivalent Python Protocol and run
  through the `PyAdapter` shim as a lower-throughput tier.

Proprietary vendor SDKs are never committed to this repository. Every adapter
must pass the shared conformance suite on recorded replay files, which are
downloaded from release assets with hash checks rather than stored in git.

See `docs/design.md`, "Adapter contract".
