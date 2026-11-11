# Adapters

Vendor adapters arrive in phase 2 (Swabian and PicoQuant first).

- In-tree adapters implement the Rust `Adapter` trait.
- Python-only vendor SDKs implement the equivalent Python Protocol and run
  through the `PyAdapter` shim as a lower-throughput tier.
- Out-of-tree Rust adapters ship as their own wheels and pass data through
  the Arrow C Stream interface.

Proprietary vendor SDKs are never committed to this repository. Every adapter
must pass the shared conformance suite on recorded replay files, which are
downloaded from release assets with hash checks rather than stored in git.

See `docs/design.md`, "Adapter contract".
