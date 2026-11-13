"""Calibration artifacts: immutable, content-addressed manifests plus a posterior table.

An artifact is two files in one directory:

- `<id>.json`: the manifest (schema `calibration/1`, see
  schema/calibration-1.schema.json), written as RFC 8785 (JCS) canonical
  JSON. The artifact ID is the SHA-256 of these bytes, so the file name, the
  ID and the content cannot drift apart.
- `<id>.posterior.npy`: the posterior draws, one row per draw and one column
  per name in `posterior.params`. The manifest records the SHA-256 of the
  canonical layout (little-endian float64, column after column in `params`
  order), not of the file bytes, so the hash survives a change of container.

Signing is planned (v0.2+); the `signature` field is excluded from the hash
when it arrives. Reading an artifact verifies both hashes.
"""

from __future__ import annotations

import contextlib
import datetime as _dt
import hashlib
import io
import json
import math
import os
import tempfile
from decimal import Decimal
from pathlib import Path
from typing import Any

import numpy as np

from dadl import __version__
from dadl.calibration.base import Calibration
from dadl.errors import DadlError

SCHEMA = "calibration/1"


class ArtifactError(DadlError):
    """An artifact is malformed or fails its hash check."""


# --- RFC 8785 (JSON Canonicalization Scheme) -------------------------------------------------


def _jcs_number(x: float | int) -> str:
    """ECMAScript Number-to-String, as JCS requires."""
    if isinstance(x, bool):
        raise TypeError("bool is not a number")
    if isinstance(x, int):
        if abs(x) > 2**53:
            raise ArtifactError(f"integer {x} is outside the IEEE double safe range")
        return str(x)
    if not math.isfinite(x):
        raise ArtifactError(f"{x} cannot be represented in JSON")
    if x == 0:
        return "0"
    sign = "-" if x < 0 else ""
    t = Decimal(repr(abs(x))).as_tuple()
    digits = "".join(map(str, t.digits)).rstrip("0") or "0"
    # value = 0.<digits> x 10^n
    n = len("".join(map(str, t.digits))) + t.exponent
    k = len(digits)
    if k <= n <= 21:
        s = digits + "0" * (n - k)
    elif 0 < n <= 21:
        s = digits[:n] + "." + digits[n:]
    elif -6 < n <= 0:
        s = "0." + "0" * (-n) + digits
    else:
        e = n - 1
        mant = digits if k == 1 else digits[0] + "." + digits[1:]
        s = f"{mant}e{'+' if e > 0 else '-'}{abs(e)}"
    return sign + s


def canonical_json(obj: Any) -> bytes:
    """Serialize `obj` as RFC 8785 canonical JSON (UTF-8 bytes)."""

    def ser(v: Any) -> str:
        if v is None:
            return "null"
        if v is True:
            return "true"
        if v is False:
            return "false"
        if isinstance(v, (int, float)):
            return _jcs_number(v)
        if isinstance(v, (np.integer, np.floating)):
            return _jcs_number(v.item())
        if isinstance(v, str):
            return json.dumps(v, ensure_ascii=False)
        if isinstance(v, (list, tuple)):
            return "[" + ",".join(ser(i) for i in v) + "]"
        if isinstance(v, dict):
            keys = list(v)
            if not all(isinstance(k, str) for k in keys):
                raise ArtifactError("JSON object keys must be strings")
            keys.sort(key=lambda k: k.encode("utf-16-be"))
            return "{" + ",".join(json.dumps(k, ensure_ascii=False) + ":" + ser(v[k]) for k in keys) + "}"
        raise ArtifactError(f"cannot serialize {type(v).__name__} in an artifact")

    return ser(obj).encode("utf-8")


def artifact_id(manifest: dict[str, Any]) -> str:
    """SHA-256 of the canonical manifest, excluding any signature."""
    body = {k: v for k, v in manifest.items() if k != "signature"}
    return hashlib.sha256(canonical_json(body)).hexdigest()


# --- building, writing, reading ---------------------------------------------------------------


def _utc(ts: _dt.datetime) -> str:
    if ts.tzinfo is None:
        raise ArtifactError("timestamps must be timezone-aware")
    return ts.astimezone(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def posterior_bytes(samples: np.ndarray) -> bytes:
    """The posterior table as .npy bytes (little-endian float64, C order)."""
    buf = io.BytesIO()
    np.save(buf, np.ascontiguousarray(samples, dtype="<f8"), allow_pickle=False)
    return buf.getvalue()


def posterior_sha256(samples: np.ndarray) -> str:
    """SHA-256 of the canonical posterior layout: little-endian float64 columns in order."""
    cols = np.ascontiguousarray(np.asarray(samples, dtype=float).T, dtype="<f8")
    return hashlib.sha256(cols.tobytes()).hexdigest()


def build_manifest(
    cal: Calibration,
    *,
    created: _dt.datetime | None = None,
    valid_until: _dt.datetime | None = None,
    supersedes: str | None = None,
) -> tuple[dict[str, Any], bytes]:
    """Manifest and posterior bytes for a calibration.

    `created` defaults to $SOURCE_DATE_EPOCH when set (reproducible builds), else now (UTC). With a
    fixed `created` the same data and code give the same artifact ID: seeds default to a hash of
    the inputs.
    """
    if created is None:
        epoch = os.environ.get("SOURCE_DATE_EPOCH")
        if epoch:
            try:
                created = _dt.datetime.fromtimestamp(int(epoch), _dt.timezone.utc)
            except ValueError as exc:
                raise ArtifactError(f"SOURCE_DATE_EPOCH must be an integer, got {epoch!r}") from exc
        else:
            created = _dt.datetime.now(_dt.timezone.utc)
    if cal.samples.ndim != 2 or cal.samples.shape[0] < 1:
        raise ArtifactError("an artifact needs at least one posterior draw")
    bad = sorted(k for k, v in cal.conditions.items() if isinstance(v, bool) or not isinstance(v, (int, float)))
    if bad:
        raise ArtifactError(f"conditions must be numbers (schema calibration/1); not numeric: {bad}")
    if valid_until is not None and valid_until < created:
        raise ArtifactError("valid_until is earlier than created")
    post = posterior_bytes(cal.samples)
    manifest: dict[str, Any] = {
        "schema": SCHEMA,
        "kind": cal.kind,
        "channels": [ch.as_record() for ch in cal.channels],
        "method": {"name": cal.method, "version": cal.method_version},
        "created": _utc(created),
        "backend": cal.backend,
        "inputs": [{"sha256": sha, "role": role} for sha, role in cal.inputs],
        "posterior": {
            "format": "npy",
            "sha256": posterior_sha256(cal.samples),
            "n_samples": int(cal.samples.shape[0]),
            "params": list(cal.sample_params),
        },
        "summary": {
            name: {"value": iv.value, "low": iv.low, "high": iv.high, "level": iv.level, "unit": cal.units.get(name, "")}
            for name, iv in cal.estimates.items()
        },
        "checks": {
            "converged": cal.converged,
            "ppc_pvalue": cal.ppc_pvalue,
            "status": cal.status,
            "prior_dominated_params": [],
        },
        "uncertainty": {
            "kind": "statistical (counting) only",
            "level": 0.95,
            "validation": "frequentist coverage checked on DADL's own simulator at the gated regimes",
            "not_included": [
                "tagger timebase and non-linearity",
                "temperature and bias sensor uncertainty",
                "connector reseat and coupling drift",
                "declared delays and window placement",
                "model misspecification beyond the predictive check",
            ],
        },
        "unidentified_params": list(cal.unidentified),
        "notes": list(cal.notes),
        "supersedes": supersedes,
        "software": {"dadl": __version__},
    }
    if cal.conditions:
        manifest["conditions"] = dict(cal.conditions)
    if valid_until is not None:
        manifest["valid_for"] = {"from": manifest["created"], "to": _utc(valid_until)}
    return manifest, post


def write_artifact(cal: Calibration, directory: str | Path, **kwargs: Any) -> Path:
    """Write `<id>.json` and `<id>.posterior.npy` into `directory`; return the manifest path.

    Existing artifacts are never overwritten: the same ID means the same bytes.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    manifest, post = build_manifest(cal, **kwargs)
    aid = artifact_id(manifest)
    mpath = directory / f"{aid}.json"
    ppath = directory / f"{aid}.posterior.npy"
    for path, data in ((ppath, post), (mpath, canonical_json(manifest))):
        if path.exists():
            if path.read_bytes() != data:
                raise ArtifactError(f"{path} exists with different content")
            continue
        # Write to a temporary file and rename, so a crash never leaves a partial artifact.
        fd, tmp = tempfile.mkstemp(dir=directory, prefix=f".{path.name}.", suffix=".tmp")
        try:
            with os.fdopen(fd, "wb") as fh:
                fh.write(data)
            os.replace(tmp, path)
        except BaseException:
            with contextlib.suppress(OSError):
                os.unlink(tmp)
            raise
    return mpath


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k, v in pairs:
        if k in out:
            raise ValueError(f"duplicate key {k!r}")
        out[k] = v
    return out


def _jcs_int(text: str) -> int | float:
    """JCS writes whole doubles of 2^53 and above as integer literals; read them back as floats."""
    value = int(text)
    return value if abs(value) <= 2**53 else float(value)


def read_artifact(path: str | Path) -> tuple[dict[str, Any], np.ndarray]:
    """Read and verify an artifact manifest and its posterior table."""
    path = Path(path)
    try:
        raw = path.read_bytes()
        manifest = json.loads(raw.decode("utf-8"), object_pairs_hook=_no_duplicates, parse_int=_jcs_int)
    except (OSError, ValueError) as exc:
        raise ArtifactError(f"cannot read artifact {path}: {exc}") from exc
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise ArtifactError(f"{path}: not a {SCHEMA} manifest")
    post_meta = manifest.get("posterior")
    if not isinstance(post_meta, dict) or not {"sha256", "n_samples", "params"} <= post_meta.keys():
        raise ArtifactError(f"{path}: the manifest has no valid posterior entry")
    aid = artifact_id(manifest)
    if path.name != f"{aid}.json":
        raise ArtifactError(f"{path}: content hash {aid} does not match the file name")
    # The ID is the hash of the file's bytes, so the bytes must be the canonical form itself: a
    # re-indented file, or one with a duplicate key that parsers resolve differently, is refused
    # (software review, cycle 2).
    if raw != canonical_json(manifest):
        raise ArtifactError(f"{path}: the file is not in canonical (RFC 8785) form")
    ppath = path.with_name(f"{aid}.posterior.npy")
    try:
        post = ppath.read_bytes()
    except OSError as exc:
        raise ArtifactError(f"posterior table missing: {ppath}") from exc
    try:
        samples = np.load(io.BytesIO(post), allow_pickle=False)
    except ValueError as exc:
        raise ArtifactError(f"{ppath}: not a .npy table: {exc}") from exc
    if samples.ndim != 2 or posterior_sha256(samples) != manifest["posterior"]["sha256"]:
        raise ArtifactError(f"{ppath}: posterior hash does not match the manifest")
    if samples.shape != (manifest["posterior"]["n_samples"], len(manifest["posterior"]["params"])):
        raise ArtifactError(f"{ppath}: table shape {samples.shape} does not match the manifest")
    return manifest, samples
