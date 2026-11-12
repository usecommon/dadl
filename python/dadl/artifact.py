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

import datetime as _dt
import hashlib
import io
import json
import math
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
    """Manifest and posterior bytes for a calibration. `created` defaults to now (UTC)."""
    created = created or _dt.datetime.now(_dt.timezone.utc)
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
        "checks": {"converged": cal.converged, "ppc_pvalue": cal.ppc_pvalue, "prior_dominated_params": []},
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
        path.write_bytes(data)
    return mpath


def read_artifact(path: str | Path) -> tuple[dict[str, Any], np.ndarray]:
    """Read and verify an artifact manifest and its posterior table."""
    path = Path(path)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ArtifactError(f"cannot read artifact {path}: {exc}") from exc
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise ArtifactError(f"{path}: not a {SCHEMA} manifest")
    aid = artifact_id(manifest)
    if path.name != f"{aid}.json":
        raise ArtifactError(f"{path}: content hash {aid} does not match the file name")
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
