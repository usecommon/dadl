//! Dead-time filters.
//!
//! These apply an ideal dead time to a sorted stream of arrivals. They are
//! used by the simulator and by tests of the closed-form rate models
//! (non-paralyzable m = n / (1 + n tau), paralyzable m = n exp(-n tau)).
//! Times and dead times are in picoseconds.

use crate::error::{check_sorted, check_span, DadlError};
use crate::Timestamp;

/// Non-paralyzable dead time: after each recorded event, arrivals within
/// `tau_ps` are lost and do not extend the dead time.
///
/// ```
/// let kept = dadl::deadtime::apply_nonparalyzable(&[0, 5, 10, 25, 30], 10).unwrap();
/// assert_eq!(kept, vec![0, 10, 25]);
/// ```
///
/// # Errors
///
/// [`DadlError::NegativeDeadTime`] if `tau_ps < 0`;
/// [`DadlError::Unsorted`] if `times` is not sorted ascending.
#[must_use = "the filtered tags are returned, not applied in place"]
pub fn apply_nonparalyzable(times: &[Timestamp], tau_ps: Timestamp) -> Result<Vec<Timestamp>, DadlError> {
    if tau_ps < 0 {
        return Err(DadlError::NegativeDeadTime);
    }
    check_sorted(times)?;
    let mut out = Vec::with_capacity(times.len());
    // i128 so that t + tau never wraps or saturates near i64::MAX.
    let mut next_allowed: Option<i128> = None;
    for &t in times {
        match next_allowed {
            Some(n) if i128::from(t) < n => continue,
            _ => {
                out.push(t);
                next_allowed = Some(i128::from(t) + i128::from(tau_ps));
            }
        }
    }
    Ok(out)
}

/// Paralyzable dead time: every arrival, recorded or not, restarts the dead
/// time, so an event is recorded only if the previous arrival is at least
/// `tau_ps` earlier.
///
/// ```
/// let kept = dadl::deadtime::apply_paralyzable(&[0, 5, 10, 25, 30], 10).unwrap();
/// assert_eq!(kept, vec![0, 25]);
/// ```
///
/// # Errors
///
/// [`DadlError::NegativeDeadTime`] if `tau_ps < 0`;
/// [`DadlError::Unsorted`] if `times` is not sorted ascending;
/// [`DadlError::Overflow`] if the time span of `times` does not fit in `i64`.
#[must_use = "the filtered tags are returned, not applied in place"]
pub fn apply_paralyzable(times: &[Timestamp], tau_ps: Timestamp) -> Result<Vec<Timestamp>, DadlError> {
    if tau_ps < 0 {
        return Err(DadlError::NegativeDeadTime);
    }
    check_sorted(times)?;
    check_span(times)?;
    let mut out = Vec::with_capacity(times.len());
    let mut last: Option<Timestamp> = None;
    for &t in times {
        // Sorted and span-checked, so t - l cannot overflow.
        let recorded = last.map_or(true, |l| t - l >= tau_ps);
        if recorded {
            out.push(t);
        }
        last = Some(t);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonparalyzable_drops_events_inside_dead_time() {
        assert_eq!(apply_nonparalyzable(&[0, 5, 10, 25, 30], 10).unwrap(), vec![0, 10, 25]);
    }

    #[test]
    fn paralyzable_extends_dead_time_on_every_arrival() {
        assert_eq!(apply_paralyzable(&[0, 5, 10, 25, 30], 10).unwrap(), vec![0, 25]);
    }

    #[test]
    fn unsorted_input_is_refused_in_every_build() {
        assert_eq!(apply_nonparalyzable(&[30, 0, 5, 10], 10), Err(DadlError::Unsorted { index: 1 }));
        assert_eq!(apply_paralyzable(&[30, 0, 5, 10], 10), Err(DadlError::Unsorted { index: 1 }));
    }

    #[test]
    fn negative_dead_time_is_refused() {
        assert_eq!(apply_nonparalyzable(&[0, 1], -1), Err(DadlError::NegativeDeadTime));
        assert_eq!(apply_paralyzable(&[0, 1], -1), Err(DadlError::NegativeDeadTime));
    }
}
