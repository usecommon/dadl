//! Per-tag processing core for DADL.
//!
//! Every loop that touches individual time tags lives here, so that no
//! per-tag loop ever runs in Python. Timestamps are `i64` picoseconds on a
//! declared clock, and every input slice is expected to be sorted ascending.
//!
//! Phase 1 contents:
//! - [`deadtime`]: non-paralyzable and paralyzable dead-time filters.
//! - [`coincidences`]: fixed-grid (clock-locked) coincidence counting, the
//!   only coincidence method allowed for Bell-type analyses.
//! - [`histogram`]: inter-arrival histograms for dead-time and afterpulse fits.

pub mod coincidences;
pub mod deadtime;
pub mod histogram;

/// A time tag in picoseconds.
pub type Timestamp = i64;
