//! Fixed-grid (clock-locked) coincidence counting.
//!
//! Time is divided into slots of a fixed width that is declared before the
//! data is seen. A coincidence is a slot occupied by both channels after each
//! channel's declared delay is removed. Windows centred on detection events
//! open the coincidence-time loophole, so Bell-type analyses must use this
//! method (see docs/design.md, "Scientific guardrails").
//!
//! Inputs need not be sorted (this module is the exception to the crate-wide
//! sorted-input rule). A streaming O(n) version for sorted input is planned.

use crate::error::DadlError;
use crate::Timestamp;
use std::cmp::Ordering;

/// A slot index: `floor((t - delay) / width)`.
pub type Slot = i64;

/// Count slots of width `width_ps` that contain at least one event from each
/// channel, after subtracting each channel's declared delay.
///
/// ```
/// let n = dadl::coincidences::grid_coincidences(&[0, 1000, 2500], &[100, 1900, 9000], 1000, 0, 0).unwrap();
/// assert_eq!(n, 2);
/// ```
///
/// # Errors
///
/// [`DadlError::NonPositiveWidth`] if `width_ps <= 0`;
/// [`DadlError::Overflow`] if removing a delay moves a tag outside `i64`.
pub fn grid_coincidences(
    a: &[Timestamp],
    b: &[Timestamp],
    width_ps: Timestamp,
    delay_a_ps: Timestamp,
    delay_b_ps: Timestamp,
) -> Result<u64, DadlError> {
    let slots_a = occupied_slots(a, width_ps, delay_a_ps)?;
    let slots_b = occupied_slots(b, width_ps, delay_b_ps)?;
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
    Ok(count)
}

/// Sorted, de-duplicated slot indices occupied by `times` after removing `delay_ps`.
///
/// # Errors
///
/// [`DadlError::NonPositiveWidth`] if `width_ps <= 0`;
/// [`DadlError::Overflow`] if `t - delay_ps` leaves the `i64` range for any tag.
pub fn occupied_slots(times: &[Timestamp], width_ps: Timestamp, delay_ps: Timestamp) -> Result<Vec<Slot>, DadlError> {
    if width_ps <= 0 {
        return Err(DadlError::NonPositiveWidth);
    }
    let mut slots = times
        .iter()
        .map(|&t| t.checked_sub(delay_ps).map(|x| x.div_euclid(width_ps)).ok_or(DadlError::Overflow))
        .collect::<Result<Vec<Slot>, DadlError>>()?;
    slots.sort_unstable();
    slots.dedup();
    Ok(slots)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn declared_delay_shifts_slots() {
        // b arrives 900 ps late; with its delay declared, both events share slot 0.
        assert_eq!(grid_coincidences(&[50], &[950], 1000, 0, 0), Ok(1));
        assert_eq!(grid_coincidences(&[50], &[1050], 1000, 0, 0), Ok(0));
        assert_eq!(grid_coincidences(&[50], &[1050], 1000, 0, 1000), Ok(1));
    }

    #[test]
    fn several_events_in_one_slot_count_once() {
        assert_eq!(grid_coincidences(&[10, 20, 30], &[40], 1000, 0, 0), Ok(1));
    }

    #[test]
    fn saturation_no_longer_makes_false_coincidences() {
        // 0.0.2 saturated here and reported 1 coincidence; the correct answer is an error.
        let r = grid_coincidences(&[i64::MIN + 5], &[i64::MIN + 500], 1000, 1000, 0);
        assert_eq!(r, Err(DadlError::Overflow));
    }

    #[test]
    fn zero_width_is_refused() {
        assert_eq!(grid_coincidences(&[0], &[0], 0, 0, 0), Err(DadlError::NonPositiveWidth));
    }
}
