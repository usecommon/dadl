//! Runs the shared vectors in tests/vectors/core.json; pytest runs the same file.

use dadl::{coincidences, deadtime, histogram, DadlError};
use serde_json::Value;

const VECTORS: &str = include_str!("../../tests/vectors/core.json");

fn ints(v: &Value) -> Vec<i64> {
    v.as_array().expect("array").iter().map(|x| x.as_i64().expect("i64")).collect()
}

fn error_name(e: DadlError) -> &'static str {
    match e {
        DadlError::Unsorted { .. } => "Unsorted",
        DadlError::NegativeDeadTime => "NegativeDeadTime",
        DadlError::NonPositiveWidth => "NonPositiveWidth",
        DadlError::Overflow => "Overflow",
        DadlError::TooManyBins { .. } => "TooManyBins",
    }
}

fn run(case: &Value) -> Result<Value, DadlError> {
    let i = |k: &str| case[k].as_i64().expect(k);
    Ok(match case["op"].as_str().expect("op") {
        "nonparalyzable" => Value::from(deadtime::apply_nonparalyzable(&ints(&case["times"]), i("tau"))?),
        "paralyzable" => Value::from(deadtime::apply_paralyzable(&ints(&case["times"]), i("tau"))?),
        "grid_coincidences" => Value::from(coincidences::grid_coincidences(
            &ints(&case["a"]),
            &ints(&case["b"]),
            i("width"),
            i("delay_a"),
            i("delay_b"),
        )?),
        "interarrival" => Value::from(histogram::interarrival(
            &ints(&case["times"]),
            i("bin_width"),
            usize::try_from(i("n_bins")).expect("n_bins"),
        )?),
        "cross_correlation" => Value::from(histogram::cross_correlation(
            &ints(&case["a"]),
            &ints(&case["b"]),
            i("bin_width"),
            i("range"),
        )?),
        op => panic!("unknown op {op}"),
    })
}

#[test]
fn shared_vectors() {
    let doc: Value = serde_json::from_str(VECTORS).expect("valid JSON");
    let cases = doc["cases"].as_array().expect("cases");
    assert!(!cases.is_empty());
    for case in cases {
        let name = case["name"].as_str().unwrap_or("?");
        let expect = &case["expect"];
        match (run(case), expect.get("error")) {
            (Ok(got), None) => assert_eq!(&got, expect, "case {name}"),
            (Err(e), Some(want)) => assert_eq!(error_name(e), want.as_str().unwrap(), "case {name}"),
            (got, _) => panic!("case {name}: expected {expect}, got {got:?}"),
        }
    }
}
