"""The calibration schema, its example artifact and the design doc agree."""

import json
from pathlib import Path

import pytest

jsonschema = pytest.importorskip("jsonschema")

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "schema" / "calibration-1.schema.json"
EXAMPLE = ROOT / "schema" / "examples" / "calibration-dark-afterpulse-deadtime.json"

pytestmark = pytest.mark.skipif(not SCHEMA.exists(), reason="repository schema/ not present")


def load(p):
    return json.loads(p.read_text(encoding="utf-8"))


def test_schema_is_valid_draft_2020_12():
    jsonschema.Draft202012Validator.check_schema(load(SCHEMA))


def test_example_artifact_validates():
    validator = jsonschema.Draft202012Validator(load(SCHEMA), format_checker=jsonschema.FormatChecker())
    validator.validate(load(EXAMPLE))


def test_design_doc_embeds_the_example_verbatim():
    design = (ROOT / "docs" / "design.md").read_text(encoding="utf-8")
    assert json.dumps(load(EXAMPLE), indent=2) in design


def test_timestamps_must_be_utc():
    schema = load(SCHEMA)
    bad = {**load(EXAMPLE), "created": "last tuesday"}
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema, format_checker=jsonschema.FormatChecker())
