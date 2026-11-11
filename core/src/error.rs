//! Errors returned by the per-tag functions.
//!
//! Every check runs in every build profile: a release build must refuse the
//! same inputs as a debug build. The Python reference raises the matching
//! exceptions (`ValueError` / `OverflowError`), pinned by the shared vectors
//! in `tests/vectors/core.json`.

use std::fmt;

/// Why a per-tag function refused its input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DadlError {
    /// Input that must be sorted ascending is not; `index` is the first
    /// element smaller than its predecessor.
    Unsorted {
        /// Index of the first out-of-order element.
        index: usize,
    },
    /// A dead time was negative.
    NegativeDeadTime,
    /// A window or bin width was zero or negative.
    NonPositiveWidth,
    /// Arithmetic on timestamps would leave the `i64` range.
    Overflow,
}

impl fmt::Display for DadlError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DadlError::Unsorted { index } => {
                write!(f, "times must be sorted ascending (first violation at index {index})")
            }
            DadlError::NegativeDeadTime => write!(f, "dead time must be non-negative"),
            DadlError::NonPositiveWidth => write!(f, "width must be positive"),
            DadlError::Overflow => write!(f, "timestamp arithmetic would leave the i64 range"),
        }
    }
}

impl std::error::Error for DadlError {}

pub(crate) fn check_sorted(times: &[i64]) -> Result<(), DadlError> {
    match times.windows(2).position(|w| w[1] < w[0]) {
        Some(i) => Err(DadlError::Unsorted { index: i + 1 }),
        None => Ok(()),
    }
}

/// For sorted input: the difference between the last and first tag must fit in `i64`.
pub(crate) fn check_span(times: &[i64]) -> Result<(), DadlError> {
    if let (Some(&first), Some(&last)) = (times.first(), times.last()) {
        last.checked_sub(first).ok_or(DadlError::Overflow)?;
    }
    Ok(())
}
