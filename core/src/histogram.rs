//! Histograms used by calibration fits.

use crate::error::{check_sorted, check_span, DadlError};
use crate::Timestamp;

/// Histogram of gaps between consecutive events on one channel.
///
/// Bin `k` counts gaps in `[k * bin_width_ps, (k + 1) * bin_width_ps)`. Gaps
/// beyond the last bin are dropped. The dead-time region shows as empty
/// leading bins and afterpulsing as an excess over the exponential tail.
///
/// ```
/// let h = dadl::histogram::interarrival(&[0, 10, 30, 60], 10, 4).unwrap();
/// assert_eq!(h, vec![0, 1, 1, 1]);
/// ```
///
/// # Errors
///
/// [`DadlError::NonPositiveWidth`] if `bin_width_ps <= 0`;
/// [`DadlError::Unsorted`] if `times` is not sorted ascending;
/// [`DadlError::Overflow`] if the time span of `times` does not fit in `i64`.
#[must_use = "the histogram is returned"]
pub fn interarrival(times: &[Timestamp], bin_width_ps: Timestamp, n_bins: usize) -> Result<Vec<u64>, DadlError> {
    if bin_width_ps <= 0 {
        return Err(DadlError::NonPositiveWidth);
    }
    check_sorted(times)?;
    check_span(times)?;
    let mut hist = vec![0u64; n_bins];
    for pair in times.windows(2) {
        // Sorted and span-checked: the gap is non-negative and fits in i64.
        let k = (pair[1] - pair[0]) / bin_width_ps;
        if let Ok(k) = usize::try_from(k) {
            if k < n_bins {
                hist[k] += 1;
            }
        }
    }
    Ok(hist)
}

/// Histogram of delays `b[j] - a[i]` over all pairs of tags within `range_ps`.
///
/// Delays `d` with `-range_ps <= d < range_ps` are counted; bin `k` covers
/// `[-range_ps + k * bin_width_ps, -range_ps + (k + 1) * bin_width_ps)`, and
/// the last of the `ceil(2 * range_ps / bin_width_ps)` bins is cut at
/// `range_ps`. A pair source shows as a peak at the relative delay of the two
/// channels, with the jitter of both as its width, on a flat background of
/// accidental pairs. Both inputs must be sorted; the cost is linear in the
/// number of tags plus the number of counted pairs.
///
/// ```
/// let h = dadl::histogram::cross_correlation(&[0, 100], &[30, 120], 20, 40).unwrap();
/// // delays 30 and 20 fall in [20, 40); 120 and -70 are out of range
/// assert_eq!(h, vec![0, 0, 0, 2]);
/// ```
///
/// # Errors
///
/// [`DadlError::NonPositiveWidth`] if `bin_width_ps <= 0` or `range_ps <= 0`;
/// [`DadlError::Unsorted`] if `a` or `b` is not sorted ascending;
/// [`DadlError::Overflow`] if `2 * range_ps`, or any tag of `b` plus or minus
/// `range_ps`, leaves the `i64` range.
#[must_use = "the histogram is returned"]
pub fn cross_correlation(
    a: &[Timestamp],
    b: &[Timestamp],
    bin_width_ps: Timestamp,
    range_ps: Timestamp,
) -> Result<Vec<u64>, DadlError> {
    if bin_width_ps <= 0 || range_ps <= 0 {
        return Err(DadlError::NonPositiveWidth);
    }
    check_sorted(a)?;
    check_sorted(b)?;
    let span = range_ps.checked_mul(2).ok_or(DadlError::Overflow)?;
    let n_bins = usize::try_from((span - 1) / bin_width_ps + 1).map_err(|_| DadlError::Overflow)?;
    let mut hist = vec![0u64; n_bins];
    let mut start = 0usize;
    for &tb in b {
        let lo = tb.checked_sub(range_ps).ok_or(DadlError::Overflow)?;
        let hi = tb.checked_add(range_ps).ok_or(DadlError::Overflow)?;
        // First a with a > tb - range_ps, i.e. delay tb - a < range_ps.
        while start < a.len() && a[start] <= lo {
            start += 1;
        }
        for &ta in &a[start..] {
            if ta > hi {
                break;
            }
            // ta in (lo, hi]: tb - ta + range_ps lies in [0, span), so it fits in i64.
            let k = (tb - ta + range_ps) / bin_width_ps;
            hist[usize::try_from(k).map_err(|_| DadlError::Overflow)?] += 1;
        }
    }
    Ok(hist)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_inputs() {
        assert_eq!(interarrival(&[], 10, 2), Ok(vec![0, 0]));
        assert_eq!(interarrival(&[5], 10, 2), Ok(vec![0, 0]));
    }

    #[test]
    fn cross_correlation_counts_every_pair_in_range() {
        // delays b - a: 5-0=5, 5-10=-5, 12-0=12 (out), 12-10=2
        assert_eq!(
            cross_correlation(&[0, 10], &[5, 12], 5, 10),
            Ok(vec![0, 1, 1, 1])
        );
        // a partial last bin: range 10, width 6 -> bins [-10,-4), [-4,2), [2,8), [8,10)
        assert_eq!(
            cross_correlation(&[0], &[9, 2], 6, 10).unwrap_err(),
            DadlError::Unsorted { index: 1 }
        );
        assert_eq!(
            cross_correlation(&[0], &[2, 9], 6, 10),
            Ok(vec![0, 0, 1, 1])
        );
        assert_eq!(
            cross_correlation(&[0], &[0], 6, 0),
            Err(DadlError::NonPositiveWidth)
        );
        assert_eq!(
            cross_correlation(&[0], &[i64::MAX], 6, 10),
            Err(DadlError::Overflow)
        );
    }

    #[test]
    fn extreme_span_is_refused() {
        assert_eq!(interarrival(&[i64::MIN, i64::MAX], 10, 4), Err(DadlError::Overflow));
    }
}
