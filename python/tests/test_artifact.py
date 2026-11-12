"""Content-addressed calibration artifacts and RFC 8785 canonical JSON."""

import datetime as dt
import json
from pathlib import Path

import numpy as np
import pytest

from dadl.artifact import (
    ArtifactError,
    artifact_id,
    build_manifest,
    canonical_json,
    read_artifact,
    write_artifact,
)
from dadl.calibration import calibrate_dark
from dadl.data import ChannelId
from dadl.sim import DetectorParams, simulate_dark

ROOT = Path(__file__).resolve().parents[2]
VECTORS = ROOT / "tests" / "vectors" / "jcs.json"
SCHEMA = ROOT / "schema" / "calibration-1.schema.json"
CREATED = dt.datetime(2026, 11, 11, 12, 0, 0, tzinfo=dt.timezone.utc)


@pytest.fixture(scope="module")
def cal():
    det = DetectorParams(dark_rate_cps=1000.0, dead_time_ps=22_000, afterpulse_prob=0.03, afterpulse_tau_ps=50_000.0)
    tags = simulate_dark(det, 10.0, np.random.default_rng(1))
    return calibrate_dark(tags, ChannelId("SPCM-0001", 3), backend="sim", seed=1, n_draws=200)


@pytest.mark.skipif(not VECTORS.exists(), reason="repository tests/vectors not present")
def test_jcs_vectors():
    v = json.loads(VECTORS.read_text(encoding="utf-8"))
    for text, expected in v["numbers"]:
        value = json.loads(text)
        if any(c in text for c in ".e"):
            value = float(value)
        assert canonical_json(value).decode() == expected, text
    for text, expected in v["objects"]:
        assert canonical_json(json.loads(text)).decode("utf-8") == expected


def test_jcs_refuses_what_json_cannot_hold():
    for bad in (float("nan"), float("inf"), 2**60, {1: "x"}, object()):
        with pytest.raises(ArtifactError):
            canonical_json(bad)


def test_id_is_stable_and_ignores_signature(cal):
    m1, _ = build_manifest(cal, created=CREATED)
    m2, _ = build_manifest(cal, created=CREATED)
    assert artifact_id(m1) == artifact_id(m2)
    signed = {**m1, "signature": {"alg": "ed25519", "key_id": "k", "sig": "s"}}
    assert artifact_id(signed) == artifact_id(m1)
    later, _ = build_manifest(cal, created=CREATED + dt.timedelta(seconds=1))
    assert artifact_id(later) != artifact_id(m1)


def test_round_trip(cal, tmp_path):
    path = write_artifact(cal, tmp_path, created=CREATED)
    assert path.name == f"{artifact_id(build_manifest(cal, created=CREATED)[0])}.json"
    manifest, samples = read_artifact(path)
    np.testing.assert_array_equal(samples, cal.samples)
    assert manifest["summary"]["dark_rate_cps"]["unit"] == "1/s"
    assert manifest["channels"] == [{"device": "SPCM-0001", "input": 3}]
    # Writing the same artifact again is a no-op, never an overwrite.
    assert write_artifact(cal, tmp_path, created=CREATED) == path


def test_tampering_is_detected(cal, tmp_path):
    path = write_artifact(cal, tmp_path, created=CREATED)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["summary"]["dark_rate_cps"]["value"] = 1.0
    path.write_bytes(canonical_json(manifest))
    with pytest.raises(ArtifactError, match="does not match the file name"):
        read_artifact(path)

    path = write_artifact(cal, tmp_path / "b", created=CREATED)
    post = path.with_name(path.stem + ".posterior.npy")
    s = np.load(post)
    s[0, 0] += 1
    np.save(post, s)
    with pytest.raises(ArtifactError, match="posterior hash"):
        read_artifact(path)


def test_naive_timestamps_are_refused(cal):
    with pytest.raises(ArtifactError):
        build_manifest(cal, created=dt.datetime(2026, 11, 11))


@pytest.mark.skipif(not SCHEMA.exists(), reason="repository schema/ not present")
def test_manifest_validates_against_the_schema(cal):
    jsonschema = pytest.importorskip("jsonschema")
    manifest, _ = build_manifest(cal, created=CREATED, valid_until=CREATED + dt.timedelta(days=7))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker()).validate(manifest)


def test_write_leaves_no_temporary_files_and_refuses_bad_manifests(cal, tmp_path):
    path = write_artifact(cal, tmp_path, created=CREATED)
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted([path.name, path.stem + ".posterior.npy"])
    manifest = json.loads(path.read_text(encoding="utf-8"))
    del manifest["posterior"]
    bad = tmp_path / "bad"
    bad.mkdir()
    target = bad / f"{artifact_id(manifest)}.json"
    target.write_bytes(canonical_json(manifest))
    with pytest.raises(ArtifactError, match="posterior"):
        read_artifact(target)
