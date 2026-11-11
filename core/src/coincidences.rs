//! Fixed-grid (clock-locked) coincidence counting.
//!
//! Time is divided into slots of a fixed `width` that is declared before the
//! data is seen. A coincidence is a slot occupied by both channels after each
//! channel's declared delay is removed. Windows centred on detection events
//! open the coincidence-time loophole, so Bell-type analyses must use this
//! method (see docs/design.md, "Scientific guardrails").

use crate::Timestamp;
use std::cmp::Ordering;

/// Count slots of width `width` that contain at least one event from each
/// channel, after subtracting each channel's declared delay.
///
/// `width` must be positive. Inputs need not be sorted.
pub fn grid_coincidences(
    a: &[Timestamp],
    b: &[Timestamp],
    width: Timestamp,
    delay_a: Timestamp,
    delay_b: Timestamp,
) -> u64 {
    assert!(width > 0, "coincidence window width must be positive");
    let slots_a = occupied_slots(a, width, delay_a);
    let slots_b = occupied_slots(b, width, delay_b);
    let (mut i, mut j, mut count) = (0usize, 0usize, 0u64);
    while i < slots_a.len() && j < slots_b.len() {
        match slots_a[i].cmp(&slots_b[j]) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
        }
    }
    count
}

/// Sorted, de-duplicated slot indices occupied by `times` after removing `delay`.
pub fn occupied_slots(times: &[Timestamp], width: Timestamp, delay: Timestamp) -> Vec<i64> {
    assert!(width > 0, "slot width must be positive");
    let mut slots: Vec<i64> = times
        .iter()
        .map(|&t| t.saturating_sub(delay).div_euclid(width))
        .collect();
    slots.sort_unstable();
    slots.dedup();
    slots
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_shared_slots() {
        let a = [0, 1000, 2500];
        let b = [100, 1900, 9000];
        assert_eq!(grid_coincidences(&a, &b, 1000, 0, 0), 2);
    }

    #[test]
    fn declared_delay_shifts_slots() {
        // b arrives 900 ps late; with its delay declared, both events share slot 0.
        let a = [50];
        let b = [950];
        assert_eq!(grid_coincidences(&a, &b, 1000, 0, 0), 1);
        let b_late = [1050];
        assert_eq!(grid_coincidences(&a, &b_late, 1000, 0, 0), 0);
        assert_eq!(grid_coincidences(&a, &b_late, 1000, 0, 1000), 1);
    }

    #[test]
    fn several_events_in_one_slot_count_once() {
        let a = [10, 20, 30];
        let b = [40];
        assert_eq!(grid_coincidences(&a, &b, 1000, 0, 0), 1);
    }

    #[test]
    fn negative_times_use_floor_division() {
        let a = [-1];
        let b = [-999];
        assert_eq!(grid_coincidences(&a, &b, 1000, 0, 0), 1);
        let c = [0];
        assert_eq!(grid_coincidences(&a, &c, 1000, 0, 0), 0);
    }
}
