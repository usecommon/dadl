//! Dead-time filters.
//!
//! These apply an ideal dead time to a sorted stream of arrivals. They are
//! used by the simulator and by tests of the closed-form rate models
//! (non-paralyzable m = n / (1 + n tau), paralyzable m = n exp(-n tau)).

use crate::Timestamp;

/// Non-paralyzable dead time: after each recorded event, arrivals within
/// `tau` are lost and do not extend the dead time.
///
/// `times` must be sorted ascending. `tau` must be non-negative.
pub fn apply_nonparalyzable(times: &[Timestamp], tau: Timestamp) -> Vec<Timestamp> {
    assert!(tau >= 0, "dead time must be non-negative");
    debug_assert!(is_sorted(times), "times must be sorted ascending");
    let mut out = Vec::with_capacity(times.len());
    let mut next_allowed: Option<Timestamp> = None;
    for &t in times {
        match next_allowed {
            Some(n) if t < n => continue,
            _ => {
                out.push(t);
                next_allowed = Some(t.saturating_add(tau));
            }
        }
    }
    out
}

/// Paralyzable dead time: every arrival, recorded or not, restarts the dead
/// time, so an event is recorded only if the previous arrival is at least
/// `tau` earlier.
///
/// `times` must be sorted ascending. `tau` must be non-negative.
pub fn apply_paralyzable(times: &[Timestamp], tau: Timestamp) -> Vec<Timestamp> {
    assert!(tau >= 0, "dead time must be non-negative");
    debug_assert!(is_sorted(times), "times must be sorted ascending");
    let mut out = Vec::with_capacity(times.len());
    let mut last: Option<Timestamp> = None;
    for &t in times {
        let recorded = match last {
            None => true,
            Some(l) => t.saturating_sub(l) >= tau,
        };
        if recorded {
            out.push(t);
        }
        last = Some(t);
    }
    out
}

fn is_sorted(times: &[Timestamp]) -> bool {
    times.windows(2).all(|w| w[0] <= w[1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nonparalyzable_drops_events_inside_dead_time() {
        let times = [0, 5, 10, 25, 30];
        assert_eq!(apply_nonparalyzable(&times, 10), vec![0, 10, 25]);
    }

    #[test]
    fn paralyzable_extends_dead_time_on_every_arrival() {
        let times = [0, 5, 10, 25, 30];
        assert_eq!(apply_paralyzable(&times, 10), vec![0, 25]);
    }

    #[test]
    fn zero_dead_time_keeps_everything() {
        let times = [0, 0, 1, 2];
        assert_eq!(apply_nonparalyzable(&times, 0), times.to_vec());
        assert_eq!(apply_paralyzable(&times, 0), times.to_vec());
    }

    #[test]
    fn empty_input() {
        assert!(apply_nonparalyzable(&[], 10).is_empty());
        assert!(apply_paralyzable(&[], 10).is_empty());
    }
}
