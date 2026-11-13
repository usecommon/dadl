//! `dadl_core`: Python bindings to the DADL per-tag core.
//!
//! Inputs are 1-D NumPy int64 arrays. C-contiguous arrays are read in place
//! (no copy); others are copied once. The GIL is held while the core runs:
//! releasing it while borrowing a NumPy buffer would let another Python thread
//! write to the array during the read, which is a data race.
//! Errors map to the exceptions the pure-Python reference raises, so both
//! backends pass the shared vectors in tests/vectors/core.json.

use dadl::{coincidences, deadtime, histogram, DadlError};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::exceptions::{PyOverflowError, PyValueError};
use pyo3::prelude::*;
use std::borrow::Cow;

fn to_py(e: DadlError) -> PyErr {
    match e {
        DadlError::Overflow => PyOverflowError::new_err(e.to_string()),
        _ => PyValueError::new_err(e.to_string()),
    }
}

/// A slice view of the array, copied only if it is not contiguous.
fn view<'a>(arr: &'a PyReadonlyArray1<'_, i64>) -> Cow<'a, [i64]> {
    match arr.as_slice() {
        Ok(s) => Cow::Borrowed(s),
        Err(_) => Cow::Owned(arr.as_array().to_vec()),
    }
}

/// Non-paralyzable dead time: drop arrivals within `tau_ps` of the last recorded one.
#[pyfunction]
fn apply_nonparalyzable<'py>(py: Python<'py>, times: PyReadonlyArray1<'py, i64>, tau_ps: i64) -> PyResult<Bound<'py, PyArray1<i64>>> {
    let t = view(&times);
    let out = deadtime::apply_nonparalyzable(&t, tau_ps).map_err(to_py)?;
    Ok(out.into_pyarray(py))
}

/// Paralyzable dead time: keep an arrival only if the previous arrival is at least `tau_ps` earlier.
#[pyfunction]
fn apply_paralyzable<'py>(py: Python<'py>, times: PyReadonlyArray1<'py, i64>, tau_ps: i64) -> PyResult<Bound<'py, PyArray1<i64>>> {
    let t = view(&times);
    let out = deadtime::apply_paralyzable(&t, tau_ps).map_err(to_py)?;
    Ok(out.into_pyarray(py))
}

/// Fixed-grid coincidences: slots of width `width_ps` occupied by both channels after their delays.
#[pyfunction]
fn grid_coincidences(
    a: PyReadonlyArray1<'_, i64>,
    b: PyReadonlyArray1<'_, i64>,
    width_ps: i64,
    delay_a_ps: i64,
    delay_b_ps: i64,
) -> PyResult<u64> {
    let (a, b) = (view(&a), view(&b));
    coincidences::grid_coincidences(&a, &b, width_ps, delay_a_ps, delay_b_ps).map_err(to_py)
}

/// Sorted, de-duplicated slot indices occupied by `times` after removing `delay_ps`.
#[pyfunction]
fn occupied_slots<'py>(py: Python<'py>, times: PyReadonlyArray1<'py, i64>, width_ps: i64, delay_ps: i64) -> PyResult<Bound<'py, PyArray1<i64>>> {
    let t = view(&times);
    let out = coincidences::occupied_slots(&t, width_ps, delay_ps).map_err(to_py)?;
    Ok(out.into_pyarray(py))
}

/// Histogram of gaps between consecutive tags.
#[pyfunction]
fn interarrival<'py>(py: Python<'py>, times: PyReadonlyArray1<'py, i64>, bin_width_ps: i64, n_bins: usize) -> PyResult<Bound<'py, PyArray1<u64>>> {
    let t = view(&times);
    let out = histogram::interarrival(&t, bin_width_ps, n_bins).map_err(to_py)?;
    Ok(out.into_pyarray(py))
}

/// Histogram of delays b - a within +-range_ps.
#[pyfunction]
fn cross_correlation<'py>(
    py: Python<'py>,
    a: PyReadonlyArray1<'py, i64>,
    b: PyReadonlyArray1<'py, i64>,
    bin_width_ps: i64,
    range_ps: i64,
) -> PyResult<Bound<'py, PyArray1<u64>>> {
    let (a, b) = (view(&a), view(&b));
    let out = histogram::cross_correlation(&a, &b, bin_width_ps, range_ps).map_err(to_py)?;
    Ok(out.into_pyarray(py))
}

// The functions read NumPy buffers borrowed as slices while the GIL is held; on a free-threaded
// interpreter that protection exists only if the module says it needs the GIL (Rust review, cycle 6).
#[pymodule(gil_used = true)]
fn dadl_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(apply_nonparalyzable, m)?)?;
    m.add_function(wrap_pyfunction!(apply_paralyzable, m)?)?;
    m.add_function(wrap_pyfunction!(grid_coincidences, m)?)?;
    m.add_function(wrap_pyfunction!(occupied_slots, m)?)?;
    m.add_function(wrap_pyfunction!(interarrival, m)?)?;
    m.add_function(wrap_pyfunction!(cross_correlation, m)?)?;
    Ok(())
}
