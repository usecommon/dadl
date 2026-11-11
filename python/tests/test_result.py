"""The labelled result type keeps raw values, labels and provenance attached."""

import pytest

from dadl.result import UNCERTAINTY_KIND, Estimate, LabelledValue


def make():
    return LabelledValue(
        analysis="hom_visibility",
        raw=Estimate(0.45, 0.02),
        corrected=Estimate(0.82, 0.04),
        labels=("policy modified",),
        policy_fingerprint="ab" * 32,
        overrides=("hom_visibility: afterpulse (test)",),
        window={"width_ps": 1000, "delays_ps": {}, "source": "plan"},
    )


def test_refuses_bare_float():
    with pytest.raises(TypeError):
        float(make())


def test_string_puts_raw_first_and_keeps_labels():
    text = str(make())
    assert text.index("raw") < text.index("corrected")
    assert "policy modified" in text


def test_export_record_has_mandatory_label_and_provenance_columns():
    record = make().as_record()
    assert record["raw_value"] == 0.45
    assert record["labels"] == "policy modified"
    assert record["policy_fingerprint"] == "ab" * 32
    assert record["overrides"] == "hom_visibility: afterpulse (test)"
    assert record["window_width_ps"] == 1000
    assert record["uncertainty"] == UNCERTAINTY_KIND
