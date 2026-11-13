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


SETUP_SCHEMA = ROOT / "schema" / "setup-1.schema.json"


def test_setup_schema_accepts_the_examples():
    yaml = pytest.importorskip("yaml")
    schema = load(SETUP_SCHEMA)
    jsonschema.Draft202012Validator.check_schema(schema)
    paths = [ROOT / d / name / "setup.yaml" for d in ("examples", "python/dadl/examples") for name in ("chsh-spad", "hom-snspd")]
    for path in paths:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        data["channels"] = {str(k): v for k, v in data["channels"].items()}  # YAML int keys -> JSON strings
        jsonschema.validate(data, schema)


@pytest.mark.parametrize(
    "patch",
    [
        {"schema_version": 2},
        {"policy": "loose"},
        {"extra": 1},
        {"analysis": {"type": "chsh", "window": {"width_ps": 1000, "source": "data"}}},
    ],
)
def test_setup_schema_rejects_bad_setups(patch):
    yaml = pytest.importorskip("yaml")
    data = yaml.safe_load((ROOT / "examples" / "chsh-spad" / "setup.yaml").read_text(encoding="utf-8"))
    data["channels"] = {str(k): v for k, v in data["channels"].items()}
    data.update(patch)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(data, load(SETUP_SCHEMA))


def test_packaged_schemas_match_the_repository_copies():
    # Software review, cycle 2: the schemas ship in the wheel so artifacts can be validated offline.
    from importlib import resources

    for name in ("calibration-1.schema.json", "setup-1.schema.json"):
        packaged = resources.files("dadl").joinpath("schemas", name).read_bytes()
        assert packaged == (ROOT / "schema" / name).read_bytes()
