//! Per-tag processing core for DADL, the Detector-Aware Data Layer.
//!
//! Every loop that touches individual time tags lives here. Python reaches it
//! through the optional `dadl-quantum-core` package (`bindings/python`); the
//! pure-Python package is the reference implementation, and both are held to
//! the shared test vectors in `tests/vectors/core.json`.
//!
//! Contract, identical in Rust and Python:
//! - timestamps are `i64` picoseconds on a declared clock;
//! - dead-time and histogram inputs must be sorted ascending, and this is
//!   checked in every build profile; coincidence inputs need not be sorted;
//! - arithmetic that would leave the `i64` range returns
//!   [`DadlError::Overflow`] instead of wrapping or saturating;
//! - invalid input returns a [`DadlError`]; no public function panics on input.
//!
//! Phase 1 contents:
//! - [`deadtime`]: non-paralyzable and paralyzable dead-time filters.
//! - [`coincidences`]: fixed-grid (clock-locked) coincidence counting, the
//!   only coincidence method allowed for Bell-type analyses.
//! - [`histogram`]: inter-arrival histograms for dead-time and afterpulse fits,
//!   and cross-correlation (delay) histograms for delay, jitter and heralding.

pub mod coincidences;
pub mod deadtime;
pub mod error;
pub mod histogram;

pub use error::DadlError;

/// A time tag in picoseconds.
pub type Timestamp = i64;
