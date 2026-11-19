"""Calibration sets: the bridge from recorded artifacts to analyses (phase 3, M1).

A `CalibrationSet` loads verified artifacts for one analysis and refuses those that must not be
used on its dataset:

- **Misfit** artifacts (the predictive check failed or the fit did not converge): refused unless
  the policy allows `misfit_calibration` for the analysis (then labelled "policy modified").
- **Outside validity** (design, "Calibration choice"): an artifact whose `valid_for` time range
  does not contain the dataset's acquisition time, or whose recorded conditions do not match the
  dataset's, raises `StaleCalibration` unless the policy allows `calibration_outside_validity`.
  Conditions compared (`condition_mismatches`): every `tagger_*` setting must be equal (cycle 2.8:
  recorded but never compared at use); `bias*` values must be equal; a `<name>_min_<unit>` /
  `<name>_max_<unit>` pair must bracket the dataset's `<name>_<unit>`; a single `temperature_K`
  must be within `TEMPERATURE_TOLERANCE_K`. A setting, bias, bracket or temperature the artifact
  records but the dataset does not give is a mismatch (council round 3: leaving the condition out
  passed the check). Keys the artifact does not record are not compared. An artifact with no
  validity range is labelled `NO_VALIDITY`.
- **Another detector** (council round 3): when the artifact's channel record and the setup's
  channel (`channels=`) both name a detector and they differ, the artifact is outside its
  validity as above. `ChannelId` equality is the tagger input; the detector is checked here.
- **Calibrated on the analysed data** (design, "Joint posteriors and data reuse"; expert review
  Q4): an artifact whose raw input hashes overlap the dataset's raises `GuardrailRefusal`. There
  is no override: reusing the data narrows the intervals, and the fix is a joint fit.

Disjointness compares whole tag-array and file hashes: a calibration made on a slice or a copy
with any change of the analysed run is not caught. Without a `DatasetContext` none of the
use-time checks runs (disjointness, validity, conditions); the set is labelled `NOT_CHECKED`.

One artifact per kind and channel set (in any order): with calibrations before and after the run,
pass the one to use; the drift between them is not propagated.

Draws (`CalibrationSet.draws`): the columns of one artifact share a sample index, so parameters
fitted together keep their correlation; different artifacts are drawn independently, each from
a generator keyed on the seed and its own artifact ID (the rule of `draws_rng`, cycle 3.1).
Stand-in columns (uniform below an upper limit, not a posterior) are never drawn: they are
returned as bounds, and phase 3 reports a corrected value at the bound and at zero, as a
labelled range (decision D3). The forward model labels which parameter is a bound
(`forward.DEAD_TIME_BOUND`, `forward.JITTER_BOUND`); the set does not.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from dadl.artifact import artifact_id, load_calibration, read_artifact
from dadl.calibration.base import Calibration
from dadl.data import ChannelId
from dadl.errors import ConfigError, GuardrailRefusal, StaleCalibration
from dadl.guardrails import POLICY_MODIFIED, Policy

#: Largest difference between a calibration's recorded temperature and the dataset's. A SPAD's
#: dark rate changes by about 5-10% per kelvin near room temperature, so 0.5 K keeps the change
#: within a few percent.
TEMPERATURE_TOLERANCE_K = 0.5
#: Relative tolerance for settings that must be equal (floats from headers and the command line).
EQUAL_RTOL = 1e-9

MISFIT_ACCEPTED = "misfit calibration accepted"
OUTSIDE_VALIDITY = "calibration outside its validity"
NO_VALIDITY = "calibration {aid} has no validity range: not checked against the acquisition time"
#: Kept for compatibility; no longer added by `CalibrationSet.load` (the forward model labels the
#: parameter that is a bound).
BOUND_ONLY = "parameter known only as an upper limit"
NOT_CHECKED = (
    "calibrations not checked against the dataset (no DatasetContext: disjointness, validity and conditions not tested)"
)
SERIES_UNUSED = (
    "deadtime-series artifacts present but not used: the forward model takes the dark dead time until M5"
)
#: Analyses whose forward model does not read deadtime-series artifacts (skipped when `kinds` is None).
_NO_SERIES = ("hom_visibility", "bell_visibility", "chsh")
_SERIES_KIND = "deadtime-series"

_RAW_ROLES = ("raw-tags", "raw-file")


@dataclass(frozen=True)
class DatasetContext:
    """What a calibration is checked against: when and under what conditions the analysed data
    were taken, and the hashes of its raw inputs (`dadl.data.tags_sha256` of each channel's tags,
    and SHA-256 of the files)."""

    acquired: _dt.datetime | None = None
    conditions: Mapping[str, float] = field(default_factory=dict)
    input_sha256: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        if self.acquired is not None and self.acquired.tzinfo is None:
            raise ValueError("acquired must be timezone-aware (UTC)")
        object.__setattr__(self, "conditions", MappingProxyType({str(k): float(v) for k, v in self.conditions.items()}))
        object.__setattr__(self, "input_sha256", frozenset(str(h) for h in self.input_sha256))


@dataclass(frozen=True)
class Member:
    """One loaded artifact."""

    artifact_id: str
    calibration: Calibration
    manifest: Mapping[str, Any]

    @property
    def kind(self) -> str:
        return self.calibration.kind

    @property
    def channels(self) -> tuple[ChannelId, ...]:
        return self.calibration.channels


def _parse_utc(text: str) -> _dt.datetime:
    return _dt.datetime.fromisoformat(text.replace("Z", "+00:00"))


def condition_mismatches(recorded: Mapping[str, float], dataset: Mapping[str, float]) -> list[str]:
    """Why a calibration's recorded conditions do not cover the dataset's (empty if they do)."""
    out: list[str] = []
    for key, value in sorted(recorded.items()):
        if key.startswith("tagger_") or key.startswith("bias"):
            if key not in dataset:
                out.append(f"{key} = {value:g} recorded, not given for the dataset")
            elif not np.isclose(dataset[key], value, rtol=EQUAL_RTOL, atol=0.0):
                out.append(f"{key} = {value:g} recorded, {dataset[key]:g} for the dataset")
        elif "_min_" in key:
            name, unit = key.split("_min_", 1)
            high_key, target = f"{name}_max_{unit}", f"{name}_{unit}"
            if high_key not in recorded:
                continue
            if target not in dataset:
                out.append(f"{key} = {value:g} recorded, not given for the dataset")
            elif not recorded[key] <= dataset[target] <= recorded[high_key]:
                out.append(f"{target} = {dataset[target]:g} outside the calibrated {recorded[key]:g} to {recorded[high_key]:g}")
        elif key == "temperature_K":
            if key not in dataset:
                out.append(f"{key} = {value:g} recorded, not given for the dataset")
            elif abs(dataset[key] - value) > TEMPERATURE_TOLERANCE_K:
                out.append(f"temperature_K = {value:g} recorded, {dataset[key]:g} for the dataset (tolerance {TEMPERATURE_TOLERANCE_K:g} K)")
    return out


def detector_mismatches(member: Member, setup: Mapping[ChannelId, ChannelId]) -> list[str]:
    """Channels whose artifact record and setup both name a detector, and differ."""
    out = []
    for ch in member.channels:
        used = setup.get(ch)
        if ch.detector and used is not None and used.detector and used.detector != ch.detector:
            out.append(f"{ch}: calibrated on detector {ch.detector}, analysing {used.detector}")
    return out


def check_detectors(
    calibrations: "CalibrationSet", channels: Iterable[ChannelId], *, analysis: str, policy: Policy | None = None
) -> tuple[str, ...]:
    """The detector check of `CalibrationSet.load(channels=)` on the analysis channels, for sets
    loaded without them (council round 4: the API paths skipped it). Refuses an artifact
    calibrated on another detector (`StaleCalibration`), or labels it under the
    `calibration_outside_validity` override."""
    setup = {c: c for c in channels}
    rule = (policy or Policy.default()).rule_for(analysis)
    labels: list[str] = []
    for m in calibrations.members:
        if set(m.channels) & set(setup):
            labels += _outside(detector_mismatches(m, setup), m.calibration.kind, m.artifact_id[:12], analysis, rule)
    return tuple(dict.fromkeys(labels))


def _check(
    member: Member, analysis: str, policy: Policy, dataset: DatasetContext | None, setup: Mapping[ChannelId, ChannelId]
) -> list[str]:
    """Refuse or label one artifact for this analysis; return its labels."""
    labels: list[str] = []
    cal, aid = member.calibration, member.artifact_id[:12]
    rule = policy.rule_for(analysis)
    if cal.status == "misfit":
        if "misfit_calibration" not in rule.allowed_corrections:
            raise GuardrailRefusal(
                f"{cal.kind} artifact {aid} is a misfit (predictive check p = {cal.ppc_pvalue:.3g}); its intervals are "
                f"not to be trusted. Override with Policy.allow({analysis!r}, 'misfit_calibration', reason)."
            )
        labels.append(MISFIT_ACCEPTED)
    raw = {i[0] for i in cal.inputs if len(i) > 1 and i[1] in _RAW_ROLES}
    problems: list[str] = detector_mismatches(member, setup)
    if dataset is None:
        return labels + _outside(problems, cal.kind, aid, analysis, rule)
    shared = raw & dataset.input_sha256
    if shared:
        raise GuardrailRefusal(
            f"{cal.kind} artifact {aid} was calibrated on the data being analysed ({len(shared)} shared input(s)): "
            "calibration data must be disjoint from analysis data, or the intervals are too narrow. Calibrate "
            "on a separate acquisition, or fit the nuisance parameters and the metric in one likelihood."
        )
    valid = member.manifest.get("valid_for")
    if not valid:
        labels.append(NO_VALIDITY.format(aid=aid))
    else:
        if dataset.acquired is None:
            problems.append("the artifact has a validity range but the dataset's acquisition time was not given")
        else:
            # `from` and `to` are optional in the schema: validity starts at acquisition, or creation,
            # without `from`, and has no end without `to` (council round 4: a bare KeyError).
            start = valid.get("from") or member.manifest.get("acquired") or member.manifest.get("created")
            end = valid.get("to")
            if (start is not None and dataset.acquired < _parse_utc(start)) or (
                end is not None and dataset.acquired > _parse_utc(end)
            ):
                problems.append(f"acquired {dataset.acquired.isoformat()} is outside {start or 'any time'} to "
                                f"{end or 'any time'}")
    problems += condition_mismatches(cal.conditions, dataset.conditions)
    return labels + _outside(problems, cal.kind, aid, analysis, rule)


def _outside(problems: list[str], kind: str, aid: str, analysis: str, rule) -> list[str]:
    """Refuse an artifact that does not cover the dataset, or label it under the override."""
    if not problems:
        return []
    if "calibration_outside_validity" not in rule.allowed_corrections:
        raise StaleCalibration(
            f"{kind} artifact {aid} does not cover this dataset: {'; '.join(problems)}. Override with "
            f"Policy.allow({analysis!r}, 'calibration_outside_validity', reason)."
        )
    return [OUTSIDE_VALIDITY]


def _artifact_files(sources: Iterable[str | Path] | str | Path) -> list[Path]:
    if isinstance(sources, (str, Path)):
        sources = [sources]
    files: list[Path] = []
    for src in sources:
        p = Path(src)
        if p.is_dir():
            files += sorted(f for f in p.glob("*.json") if not f.name.startswith("."))
        else:
            files.append(p)
    return files


@dataclass(frozen=True)
class JointDraws:
    """Draws of every posterior column of a set, `n` rows, keyed by (kind, param, channels)."""

    n: int
    columns: Mapping[tuple[str, str, tuple[ChannelId, ...]], np.ndarray]
    #: (kind, param, channels) -> (0, upper limit) for stand-in parameters, which are not drawn.
    bounds: Mapping[tuple[str, str, tuple[ChannelId, ...]], tuple[float, float]]

    def column(self, kind: str, param: str, channel: ChannelId | None = None) -> np.ndarray:
        return self.columns[self._key(self.columns, kind, param, channel)]

    def bound(self, kind: str, param: str, channel: ChannelId | None = None) -> tuple[float, float]:
        return self.bounds[self._key(self.bounds, kind, param, channel)]

    @staticmethod
    def _key(table: Mapping, kind: str, param: str, channel: ChannelId | None) -> tuple:
        hits = [k for k in table if k[0] == kind and k[1] == param and (channel is None or channel in k[2])]
        if len(hits) != 1:
            where = f" for {channel}" if channel is not None else ""
            raise KeyError(f"{len(hits)} matches for {kind}.{param}{where}")
        return hits[0]


@dataclass(frozen=True)
class CalibrationSet:
    """Artifacts checked for one analysis and dataset."""

    analysis: str
    members: tuple[Member, ...]
    labels: tuple[str, ...]
    #: Whether the use-time checks ran against a `DatasetContext`.
    dataset_checked: bool = True
    #: The `DatasetContext` the set was checked against, if any.
    dataset: DatasetContext | None = None

    def dataset_record(self) -> dict[str, object] | None:
        """The dataset context the artifacts were checked against, for a result's counts."""
        if self.dataset is None:
            return None
        acquired = self.dataset.acquired
        return {
            "acquired": acquired.isoformat() if acquired is not None else None,
            "conditions": dict(sorted(self.dataset.conditions.items())),
        }

    @classmethod
    def load(
        cls,
        sources: Iterable[str | Path] | str | Path,
        *,
        analysis: str,
        policy: Policy | None = None,
        dataset: DatasetContext | None = None,
        channels: Iterable[ChannelId] | None = None,
        kinds: Iterable[str] | None = None,
    ) -> CalibrationSet:
        """Load and check artifacts from files or directories of `<id>.json` manifests.

        `channels` and `kinds` select which artifacts to keep (an artifact is kept if any of its
        channels is listed). Two kept artifacts of the same kind and channels (in any order) are
        ambiguous and refused. The detector of each listed channel, when given, is compared with the
        artifact's. Without `dataset`, only the misfit check runs and the set is labelled
        `NOT_CHECKED`; analyses must pass one. With `kinds` None, deadtime-series artifacts are
        skipped for the HOM, fringe and CHSH analyses (their forward model does not read them) and
        the set is labelled `SERIES_UNUSED`.
        """
        policy = policy or Policy.default()
        setup = {c: c for c in channels} if channels is not None else {}
        want_ch = set(setup) if channels is not None else None
        want_kind = set(kinds) if kinds is not None else None
        members: list[Member] = []
        skipped_series = False
        for path in _artifact_files(sources):
            manifest, _ = read_artifact(path)
            cal = load_calibration(path)
            if want_kind is not None and cal.kind not in want_kind:
                continue
            if want_kind is None and cal.kind == _SERIES_KIND and analysis in _NO_SERIES:
                skipped_series = True
                continue
            if want_ch is not None and not want_ch & set(cal.channels):
                continue
            members.append(Member(artifact_id(manifest), cal, MappingProxyType(manifest)))
        seen: dict[tuple[str, frozenset[ChannelId]], str] = {}
        for m in members:
            key = (m.kind, frozenset(m.channels))
            if key in seen and seen[key] != m.artifact_id:
                raise ConfigError(
                    f"two {m.kind} artifacts for {', '.join(map(str, m.channels))} ({seen[key][:12]}, "
                    f"{m.artifact_id[:12]}): pass the one to use"
                )
            seen[key] = m.artifact_id
        members = list({m.artifact_id: m for m in members}.values())
        labels: list[str] = []
        for m in members:
            labels += _check(m, analysis, policy, dataset, setup)
        if skipped_series:
            labels.append(SERIES_UNUSED)
        if dataset is None:
            labels.append(NOT_CHECKED)
        if policy.is_modified:
            labels.insert(0, POLICY_MODIFIED)
        return cls(analysis, tuple(members), tuple(dict.fromkeys(labels)), dataset is not None, dataset)

    @property
    def artifact_ids(self) -> tuple[str, ...]:
        return tuple(m.artifact_id for m in self.members)

    def get(self, kind: str, channel: ChannelId | None = None) -> Calibration:
        hits = [m.calibration for m in self.members if m.kind == kind and (channel is None or channel in m.channels)]
        if len(hits) != 1:
            raise KeyError(f"{len(hits)} {kind} calibrations{'' if channel is None else f' for {channel}'} in the set")
        return hits[0]

    def draws(self, n: int, seed: int) -> JointDraws:
        """`n` joint draws: one shared row index per artifact, independent across artifacts."""
        if n < 1:
            raise ValueError("n must be at least 1")
        columns: dict[tuple, np.ndarray] = {}
        bounds: dict[tuple, tuple[float, float]] = {}
        for m in self.members:
            cal = m.calibration
            words = [int.from_bytes(hashlib.sha256(m.artifact_id.encode()).digest()[i : i + 4], "little") for i in range(0, 16, 4)]
            rng = np.random.default_rng(np.random.SeedSequence([int(seed) & 0xFFFFFFFF, int(seed) >> 32 & 0xFFFFFFFF, *words]))
            rows = cal.samples.shape[0]
            idx = rng.permutation(rows)[:n] if n <= rows else rng.integers(0, rows, n)
            for j, name in enumerate(cal.sample_params):
                key = (cal.kind, name, cal.channels)
                if name in cal.stand_in:
                    bounds[key] = (0.0, float(cal.estimates[name].high))
                    continue
                col = cal.samples[idx, j].copy()
                col.setflags(write=False)
                columns[key] = col
        return JointDraws(n, MappingProxyType(columns), MappingProxyType(bounds))


__all__ = [
    "BOUND_ONLY",
    "MISFIT_ACCEPTED",
    "NO_VALIDITY",
    "NOT_CHECKED",
    "SERIES_UNUSED",
    "OUTSIDE_VALIDITY",
    "TEMPERATURE_TOLERANCE_K",
    "CalibrationSet",
    "DatasetContext",
    "JointDraws",
    "Member",
    "condition_mismatches",
    "check_detectors",
    "detector_mismatches",
]
