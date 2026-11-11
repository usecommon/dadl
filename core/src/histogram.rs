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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn short_inputs() {
        assert_eq!(interarrival(&[], 10, 2), Ok(vec![0, 0]));
        assert_eq!(interarrival(&[5], 10, 2), Ok(vec![0, 0]));
    }

    #[test]
    fn extreme_span_is_refused() {
        assert_eq!(interarrival(&[i64::MIN, i64::MAX], 10, 4), Err(DadlError::Overflow));
    }
}
