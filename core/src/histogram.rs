//! Histograms used by calibration fits.

use crate::Timestamp;

/// Histogram of gaps between consecutive events on one channel.
///
/// Bin `k` counts gaps in `[k * bin_width, (k + 1) * bin_width)`. Gaps beyond
/// the last bin are dropped. The dead-time region shows as empty leading bins
/// and afterpulsing as an excess over the exponential tail.
///
/// `times` must be sorted ascending; `bin_width` must be positive.
pub fn interarrival(times: &[Timestamp], bin_width: Timestamp, n_bins: usize) -> Vec<u64> {
    assert!(bin_width > 0, "bin width must be positive");
    let mut hist = vec![0u64; n_bins];
    for pair in times.windows(2) {
        let gap = pair[1] - pair[0];
        if gap < 0 {
            continue;
        }
        let k = (gap / bin_width) as u64;
        if k < n_bins as u64 {
            hist[k as usize] += 1;
        }
    }
    hist
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_gaps() {
        let times = [0, 10, 30, 60];
        assert_eq!(interarrival(&times, 10, 4), vec![0, 1, 1, 1]);
    }

    #[test]
    fn drops_gaps_beyond_range() {
        let times = [0, 100];
        assert_eq!(interarrival(&times, 10, 4), vec![0, 0, 0, 0]);
    }

    #[test]
    fn short_inputs() {
        assert_eq!(interarrival(&[], 10, 2), vec![0, 0]);
        assert_eq!(interarrival(&[5], 10, 2), vec![0, 0]);
    }
}
