"""Calibration sets: the bridge from recorded artifacts to analyses (phase 3, M1).

A `CalibrationSet` loads verified artifacts for one analysis and refuses those that must not be
used on its dataset:

- **Misfit** artifacts (the predictive check failed or the fit did not converge): refused unless
  the policy allows `misfit_calibration` for the analysis (then labelled "policy modified").
- **Outside validity** (design, "Calibration choice"): an artifact whose `valid_for` time range
  does not contain the dataset's acquisition time, or whose recorded conditions do not match the
  dataset's, raises `StaleCalibration` unless the policy allows `calibration_outside_validity`.
  Conditions compared (`condition_mismatches`): every `tagger_*` setting must be equal (cycle 2.8:
  recorded but never compared at use; the T3 sync period to 1e-4, `SETTING_RTOL`); `bias*` values must be equal; a `<name>_min_<unit>` /
  `<name>_max_<unit>` pair must bracket the dataset's `<name>_<unit>`; a single `temperature_K`
  must be within `TEMPERATURE_TOLERANCE_K`. A setting, bias, bracket or temperature the artifact
  records but the dataset does not give is a mismatch (council round 3: leaving the condition out
  passed the check). Keys the artifact does not record are not compared. An artifact with no
  validity range is labelled `NO_VALIDITY`.
- **Another tagger** (council round 7, G6): when the artifact's `instrument` and the dataset's
  both give the tagger model or serial and they differ, a delay-jitter or deadtime-series
  artifact is outside its validity as above; other kinds are labelled `OTHER_TAGGER`.
- **Another detector** (council round 3): when the artifact's channel record and the setup's
  channel (`channels=`) both name a detector and they differ, the artifact is outside its
  validity as above. `ChannelId` equality is the tagger input; the detector is checked here.
- **Method version** (council round 5, V4): an artifact whose method version is below
  `METHOD_MINIMUM` for its (kind, method) is refused (`GuardrailRefusal`) unless the policy allows
  `outdated_calibration_method` (then labelled `OUTDATED_ACCEPTED`); one below `METHOD_LABELLED`'s
  version loads with a label saying what it lacks. Methods the tables do not list are not checked.
  Delay artifacts from 0.1.x (method 8) are labelled, not refused (council round 10, L3).
- **Scan range** (council round 6, T6): a deadtime-series artifact records the count rates it was
  scanned at (`rate_min_cps`, `rate_max_cps`); the analysed channel's recorded count rate
  (`DatasetContext.rates_cps`, or a `rate_cps` condition for every channel) must lie within them,
  as any `_min_`/`_max_` bracket: an effective tau_p describes the scanned rates only.
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
a generator keyed on the seed and its own artifact ID (the rule of `draws_rng`, cycle 3.1). The
exception is a deadtime-series artifact and the dark artifact it cites (council round 6, B5): its
row i was built from the dark artifact's row i, so both take the dark artifact's row index and
tau_p keeps its correlation with the afterpulse draws.
Stand-in columns (uniform below an upper limit, not a posterior) are never drawn: they are
returned as bounds, and phase 3 reports a corrected value at the bound and at zero, as a
labelled range (decision D3). The forward model labels which parameter is a bound
(`forward.DEAD_TIME_BOUND`, `forward.JITTER_BOUND`); the set does not.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import math
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
from dadl.guardrails import POLICY_MODIFIED, Policy, override_hint

#: Largest difference between a calibration's recorded temperature and the dataset's. A SPAD's
#: dark rate changes by about 5-10% per kelvin near room temperature, so 0.5 K keeps the change
#: within a few percent.
TEMPERATURE_TOLERANCE_K = 0.5
#: Relative tolerance for settings that must be equal (floats from headers and the command line).
EQUAL_RTOL = 1e-9
#: Per-setting tolerances (council round 7, G2). In T3, MeasDesc_GlobalResolution is the source's
#: repetition period, a laser property that drifts at the ppm level, not a setting that places
#: tags: 1e-4 admits the drift and refuses another repetition rate. In T2 it is the tick, whose
#: values are discrete, so 1e-4 is still exact there.
SETTING_RTOL = {"tagger_MeasDesc_GlobalResolution": 1e-4}
#: The T2 tick or T3 sync period (one header key for both), and the mode that says which it is
#: (council round 9, Tw2: a T2 dark artifact's tick was compared with a T3 dataset's period).
PERIOD_KEY = "tagger_MeasDesc_GlobalResolution"
MODE_KEY = "tagger_mode"
#: Kinds whose values do not depend on the T3 sync clock under CW light (council round 9, Tw3).
CLOCK_FREE_KINDS = ("dark-afterpulse-deadtime", "delay-jitter", "accidentals", "heralding-efficiency")
T3_CLOCK_CHANGED = (
    "T3 sync clock changed since calibration ({a} -> {b}); with CW light these calibrations do not depend on it"
)
#: A recorded condition the comparison does not understand (council round 9, Rd4).
NOT_COMPARED = "condition {key} recorded, not compared"
#: Oldest method version accepted at use, per (kind, method) (council round 5, V4). Set to the
#: oldest version a release or the main branch wrote: every artifact in use today loads. A fix that
#: makes older artifacts wrong raises the minimum (and says so in the CHANGELOG).
METHOD_MINIMUM: Mapping[tuple[str, str], int] = MappingProxyType({
    ("accidentals", "offpeak-grid"): 3,
    ("dark-afterpulse-deadtime", "interarrival-mle"): 7,
    ("dcr-temperature", "arrhenius-wls"): 2,
    ("delay-jitter", "xcorr-gaussian-mle"): 8,
    ("heralding-efficiency", "herald-ratio"): 6,
    ("deadtime-series", "series-gap-mle"): 8,
})
#: Versions at or above the minimum that load with a label: (kind, method) -> (the first version
#: that needs none, what the older ones lack).
METHOD_LABELLED: Mapping[tuple[str, str], tuple[int, str]] = MappingProxyType({
    # Council round 10 (Mg7, L3): 0.1.x delay artifacts still load.
    ("delay-jitter", "xcorr-gaussian-mle"): (
        10, "its misfit flag was computed before council round 8's (groups - 1) degrees of freedom; refit it if "
        "that matters"),
    ("deadtime-series", "series-gap-mle"): (
        9, "its tau_p draws do not share the dark artifact's rows and its scan rates are not checked at use; refit it"),
})
OLD_METHOD = "{kind} artifact {aid} from {method} version {version}: {why}"
OUTDATED_ACCEPTED = "outdated calibration method accepted"
#: The scanned count-rate bracket of a deadtime-series artifact, and the dataset key it brackets.
SCAN_RATE_KEYS = ("rate_min_cps", "rate_max_cps")
RATE_KEY = "rate_cps"
#: Kinds whose artifact belongs to the tagger that recorded it (its channel skew and own dead
#: time): calibrated on another tagger, they do not cover the dataset (council round 7, G6).
TAGGER_KINDS = ("delay-jitter", "deadtime-series")
OTHER_TAGGER = "calibration taken on another tagger"

MISFIT_ACCEPTED = "misfit calibration accepted"
OUTSIDE_VALIDITY = "calibration outside its validity"
NO_VALIDITY = "calibration {aid} has no validity range: not checked against the acquisition time"
#: `NO_VALIDITY` for several artifacts, in one label (council round 9, Tw8).
NO_VALIDITY_MANY = "calibrations {aids} have no validity range: not checked against the acquisition time"
#: Kept for compatibility; no longer added by `CalibrationSet.load` (the forward model labels the
#: parameter that is a bound).
BOUND_ONLY = "parameter known only as an upper limit"
NOT_CHECKED = (
    "calibrations not checked against the dataset (no DatasetContext: disjointness, validity and conditions not tested)"
)
SERIES_UNUSED = (
    "deadtime-series artifacts present but not used: the forward model takes the dark dead time (not yet supported; "
    "see Known limits)"
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
    #: The tagger that recorded the dataset (`model`, `serial`), as the artifact's `instrument`.
    instrument: Mapping[str, str] = field(default_factory=dict)
    #: T3 input under an explicit `source.type: cw`: the sync is a clock only, and a changed clock
    #: is labelled for `CLOCK_FREE_KINDS` instead of refused (council round 9, Tw3).
    t3_cw: bool = False
    #: Each analysed channel's recorded count rate (s^-1), checked against a deadtime-series
    #: artifact's scanned rates (council round 6, T6). Overrides a `rate_cps` condition.
    rates_cps: Mapping[ChannelId, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.acquired is not None and self.acquired.tzinfo is None:
            raise ValueError("acquired must be timezone-aware (UTC)")
        object.__setattr__(self, "conditions", MappingProxyType({str(k): float(v) for k, v in self.conditions.items()}))
        object.__setattr__(self, "instrument", MappingProxyType(
            {str(k): str(v) for k, v in self.instrument.items() if v not in (None, "")}))
        object.__setattr__(self, "input_sha256", frozenset(str(h) for h in self.input_sha256))
        object.__setattr__(self, "rates_cps", MappingProxyType({k: float(v) for k, v in self.rates_cps.items()}))


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


def setting_rtol(key: str) -> float:
    """The relative tolerance within which a recorded setting must equal the dataset's."""
    return SETTING_RTOL.get(key, EQUAL_RTOL)


def _modes_differ(recorded: Mapping[str, float], dataset: Mapping[str, float]) -> bool:
    return MODE_KEY in recorded and MODE_KEY in dataset and recorded[MODE_KEY] != dataset[MODE_KEY]


#: What to do about a recorded condition the dataset does not give (council round 10, Ux5).
CONDITION_NOT_GIVEN = "{key} = {value} recorded, not given for the dataset: pass --condition {target}=<value at acquisition>"


def condition_mismatches(
    recorded: Mapping[str, float], dataset: Mapping[str, float], *, hint: bool = False
) -> list[str]:
    """Why a calibration's recorded conditions do not cover the dataset's (empty if they do).

    `tagger_MeasDesc_GlobalResolution` is the T2 tick or the T3 sync period: it is compared only
    when both sides record the same `tagger_mode`, or either records none (council round 9, Tw2).
    A lone `<name>_min_<unit>` or `<name>_max_<unit>` is compared one-sided (council round 9, Rd4).
    Values are printed to 10 significant digits (council round 7, G2: `%g` showed two sync
    periods 1e-7 apart as equal). With `hint`, a condition not given for the dataset says which
    `--condition` to pass (`CONDITION_NOT_GIVEN`)."""
    out: list[str] = []

    def not_given(key: str, value: float, target: str) -> str:
        if hint:
            return CONDITION_NOT_GIVEN.format(key=key, value=f"{value:.10g}", target=target)
        return f"{key} = {value:.10g} recorded, not given for the dataset"

    for key, value in sorted(recorded.items()):
        if key == MODE_KEY or (key == PERIOD_KEY and _modes_differ(recorded, dataset)):
            continue
        if key.startswith("tagger_") or key.startswith("bias"):
            if key not in dataset:
                out.append(not_given(key, value, key))
            elif not np.isclose(dataset[key], value, rtol=setting_rtol(key), atol=0.0):
                out.append(f"{key} = {value:.10g} recorded, {dataset[key]:.10g} for the dataset")
        elif "_min_" in key or "_max_" in key:
            low = "_min_" in key
            name, unit = key.split("_min_" if low else "_max_", 1)
            low_key, high_key, target = f"{name}_min_{unit}", f"{name}_max_{unit}", f"{name}_{unit}"
            if not low and low_key in recorded:
                continue  # the pair is compared once, at its _min_ key
            lo, hi = recorded.get(low_key, -math.inf), recorded.get(high_key, math.inf)
            if target not in dataset:
                out.append(not_given(key, value, target))
            elif not lo <= dataset[target] <= hi:
                out.append(f"{target} = {dataset[target]:.10g} outside the calibrated {lo:.10g} to {hi:.10g}")
        elif key == "temperature_K":
            if key not in dataset:
                out.append(not_given(key, value, key))
            elif abs(dataset[key] - value) > TEMPERATURE_TOLERANCE_K:
                out.append(f"temperature_K = {value:.10g} recorded, {dataset[key]:.10g} for the dataset (tolerance "
                           f"{TEMPERATURE_TOLERANCE_K:g} K)")
    return out


def _compared(key: str) -> bool:
    return (key.startswith(("tagger_", "bias")) or "_min_" in key or "_max_" in key or key == "temperature_K")


def uncompared_conditions(recorded: Mapping[str, float], dataset: Mapping[str, float]) -> list[str]:
    """Recorded conditions outside the compared forms that look like a condition of use: the dataset
    gives the same key (`analyze --condition`), or the name says temperature (council round 9,
    Rd4: `detector_temperature_K` 233 against 298 passed silently). The fit's own diagnostics are
    not listed, nor `dark_*` keys: a deadtime-series artifact's copy of its dark artifact's
    conditions, compared with the scan's when it was calibrated (council round 10, Rg3)."""
    return sorted(k for k in recorded
                  if not _compared(k) and not k.startswith("dark_") and (k in dataset or "temperature" in k))


def clock_changes(recorded: Mapping[str, float], dataset: Mapping[str, float]) -> str | None:
    """The T3 sync period change between artifact and dataset, as (recorded -> dataset), or None."""
    if PERIOD_KEY not in recorded or PERIOD_KEY not in dataset or _modes_differ(recorded, dataset):
        return None
    a, b = recorded[PERIOD_KEY], dataset[PERIOD_KEY]
    if np.isclose(b, a, rtol=setting_rtol(PERIOD_KEY), atol=0.0):
        return None
    return T3_CLOCK_CHANGED.format(a=f"{a:.10g} s", b=f"{b:.10g} s")


def collapse_no_validity(labels: Iterable[str]) -> tuple[str, ...]:
    """The labels with every `NO_VALIDITY` label merged into one that lists the IDs (council round 9,
    Tw8: CHSH repeated it eight times), at the first one's place."""
    head, tail = NO_VALIDITY.split("{aid}")
    labels = list(dict.fromkeys(labels))
    ids = [lab[len(head) : len(lab) - len(tail)] for lab in labels if lab.startswith(head) and lab.endswith(tail)]
    if len(ids) < 2:
        return tuple(labels)
    merged = NO_VALIDITY_MANY.format(aids=", ".join(ids))
    out: list[str] = []
    for lab in labels:
        if lab.startswith(head) and lab.endswith(tail):
            if merged not in out:
                out.append(merged)
        else:
            out.append(lab)
    return tuple(out)


def instrument_mismatches(recorded: Mapping[str, Any] | None, dataset: Mapping[str, Any] | None) -> list[str]:
    """The tagger model and serial where the artifact's `instrument` record and the dataset's both
    give them and differ (council round 7, G6)."""
    out = []
    for key in ("model", "serial"):
        a, b = (recorded or {}).get(key), (dataset or {}).get(key)
        if a and b and str(a).strip() != str(b).strip():
            out.append(f"tagger {key} {a} recorded, {b} for the dataset")
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


def _method_version(member: Member, analysis: str, rule) -> list[str]:
    """Refuse an artifact below its method's minimum version, or label it (council round 5, V4)."""
    cal, aid = member.calibration, member.artifact_id[:12]
    key = (cal.kind, cal.method)
    if key not in METHOD_MINIMUM:
        return []
    text = str(cal.method_version).strip()
    version = int(text) if text.isdigit() else None
    minimum = METHOD_MINIMUM[key]
    if version is None or version < minimum:
        if "outdated_calibration_method" not in rule.allowed_corrections:
            raise GuardrailRefusal(
                f"{cal.kind} artifact {aid} was written by {cal.method} version {text}; the oldest version accepted "
                f"at use is {minimum}: refit it with this version of DADL. "
                f"{override_hint(analysis, 'outdated_calibration_method')}"
            )
        return [OUTDATED_ACCEPTED]
    if key in METHOD_LABELLED and version < METHOD_LABELLED[key][0]:
        return [OLD_METHOD.format(kind=cal.kind, aid=aid, method=cal.method, version=version,
                                  why=METHOD_LABELLED[key][1])]
    return []


def _dataset_conditions(member: Member, dataset: DatasetContext) -> Mapping[str, float]:
    """The dataset's conditions for one artifact: for a deadtime-series artifact that records its
    scanned rates, `rate_cps` is its channel's analysed rate when the dataset gives it (T6)."""
    if member.kind != _SERIES_KIND or not any(k in member.calibration.conditions for k in SCAN_RATE_KEYS):
        return dataset.conditions
    rates = [dataset.rates_cps[c] for c in member.channels if c in dataset.rates_cps]
    if not rates:
        return dataset.conditions
    # A series artifact has one channel.
    return {**dataset.conditions, RATE_KEY: rates[0]}


def _check(
    member: Member, analysis: str, policy: Policy, dataset: DatasetContext | None, setup: Mapping[ChannelId, ChannelId],
    accepted: list[str] | None = None,
) -> list[str]:
    """Refuse or label one artifact for this analysis; return its labels."""
    labels: list[str] = []
    cal, aid = member.calibration, member.artifact_id[:12]
    rule = policy.rule_for(analysis)
    labels += _method_version(member, analysis, rule)
    if cal.status == "misfit":
        if "misfit_calibration" not in rule.allowed_corrections:
            raise GuardrailRefusal(
                f"{cal.kind} artifact {aid} is a misfit (predictive check p = {cal.ppc_pvalue:.3g}); its intervals are "
                f"not to be trusted. {override_hint(analysis, 'misfit_calibration')}"
            )
        labels.append(MISFIT_ACCEPTED)
    raw = {i[0] for i in cal.inputs if len(i) > 1 and i[1] in _RAW_ROLES}
    problems: list[str] = detector_mismatches(member, setup)
    if dataset is None:
        return labels + _outside(problems, cal.kind, aid, analysis, rule, accepted)
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
    # valid_for.conditions are enforced too, merged with the top-level conditions; a key in both
    # must match in both (council round 5, V5).
    extra = (valid or {}).get("conditions") or {}
    merged = {**cal.conditions, **extra}
    given = _dataset_conditions(member, dataset)
    found = condition_mismatches(cal.conditions, given, hint=True)
    if extra:
        found += condition_mismatches(merged, given, hint=True)
    changed = clock_changes(merged, given) if dataset.t3_cw and cal.kind in CLOCK_FREE_KINDS else None
    if changed:
        found = [f for f in found if not f.startswith(f"{PERIOD_KEY} = ")]
        labels.append(changed)
    labels += [NOT_COMPARED.format(key=k) for k in uncompared_conditions(merged, given)]
    problems += list(dict.fromkeys(found))
    tagger = instrument_mismatches(member.manifest.get("instrument"), dataset.instrument)
    if tagger and cal.kind in TAGGER_KINDS:
        problems += tagger
    elif tagger:
        labels.append(OTHER_TAGGER)
    return labels + _outside(problems, cal.kind, aid, analysis, rule, accepted)


def _outside(
    problems: list[str], kind: str, aid: str, analysis: str, rule, accepted: list[str] | None = None
) -> list[str]:
    """Refuse an artifact that does not cover the dataset, or label it under the override (and add
    the problems to `accepted`, for `CalibrationSet.overridden`)."""
    if not problems:
        return []
    if "calibration_outside_validity" not in rule.allowed_corrections:
        unknown = any("pass --condition" in p for p in problems)
        raise StaleCalibration(
            f"{kind} artifact {aid} does not cover this dataset: {'; '.join(problems)}. "
            f"{'If a condition is unknown: ' if unknown else ''}{override_hint(analysis, 'calibration_outside_validity')}"
        )
    if accepted is not None:
        accepted.extend(f"{kind} artifact {aid}: {p}" for p in problems)
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
    #: Labels each artifact's checks added, by artifact ID (`labels_for`).
    member_labels: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    #: What `calibration_outside_validity` let through: "<kind> artifact <id>: <problem>", one per
    #: problem (printed by `dadl analyze`; council round 10, Ux1).
    overridden: tuple[str, ...] = ()

    def labels_for(self, artifact_ids: Iterable[str]) -> tuple[str, ...]:
        """The set's labels, without those that only artifacts outside `artifact_ids` added (council
        round 6, R2: a CHSH result carried "no validity range" for artifacts it did not use)."""
        ids = set(artifact_ids)
        other = {lab for aid, labs in self.member_labels.items() if aid not in ids for lab in labs}
        own = {lab for aid, labs in self.member_labels.items() if aid in ids for lab in labs}
        return collapse_no_validity(lab for lab in self.labels if lab not in other or lab in own)

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
        per: dict[str, tuple[str, ...]] = {}
        accepted: list[str] = []
        for m in members:
            per[m.artifact_id] = tuple(_check(m, analysis, policy, dataset, setup, accepted))
            labels += per[m.artifact_id]
        if skipped_series:
            labels.append(SERIES_UNUSED)
        if dataset is None:
            labels.append(NOT_CHECKED)
        if policy.is_modified:
            labels.insert(0, POLICY_MODIFIED)
        return cls(analysis, tuple(members), tuple(dict.fromkeys(labels)), dataset is not None, dataset,
                   MappingProxyType(per), tuple(accepted))

    @property
    def artifact_ids(self) -> tuple[str, ...]:
        return tuple(m.artifact_id for m in self.members)

    def get(self, kind: str, channel: ChannelId | None = None) -> Calibration:
        hits = [m.calibration for m in self.members if m.kind == kind and (channel is None or channel in m.channels)]
        if len(hits) != 1:
            raise KeyError(f"{len(hits)} {kind} calibrations{'' if channel is None else f' for {channel}'} in the set")
        return hits[0]

    def cited_dark(self, member: Member) -> Member | None:
        """The set's dark artifact a deadtime-series artifact was calibrated with, or None: the one
        its `previous-artifact` input names (`dadl calibrate deadtime --dark`), else the one on its
        channel whose raw inputs it cites (the Python API copies them)."""
        if member.kind != _SERIES_KIND:
            return None
        darks = [m for m in self.members if m.kind == "dark-afterpulse-deadtime" and set(m.channels) & set(member.channels)]
        named = {i[0] for i in member.calibration.inputs if len(i) > 1 and i[1] == "previous-artifact"}
        for d in darks:
            if d.artifact_id in named:
                return d
        cited = set(member.calibration.inputs)
        for d in darks:
            raw = {tuple(i) for i in d.calibration.inputs if len(i) > 1 and i[1] in _RAW_ROLES}
            if raw and raw <= cited:
                return d
        return None

    def draws(self, n: int, seed: int) -> JointDraws:
        """`n` joint draws: one shared row index per artifact, independent across artifacts, except
        that a deadtime-series artifact takes the row index of the dark artifact it cites (B5)."""
        if n < 1:
            raise ValueError("n must be at least 1")
        columns: dict[tuple, np.ndarray] = {}
        bounds: dict[tuple, tuple[float, float]] = {}
        index: dict[str, np.ndarray] = {}
        for m in self.members:
            rows = m.calibration.samples.shape[0]
            words = [int.from_bytes(hashlib.sha256(m.artifact_id.encode()).digest()[i : i + 4], "little") for i in range(0, 16, 4)]
            rng = np.random.default_rng(np.random.SeedSequence([int(seed) & 0xFFFFFFFF, int(seed) >> 32 & 0xFFFFFFFF, *words]))
            index[m.artifact_id] = rng.permutation(rows)[:n] if n <= rows else rng.integers(0, rows, n)
        for m in self.members:
            dark = self.cited_dark(m)
            # Series row i holds the dark artifact's row i mod its rows: shared when the dark rows fit.
            if dark is not None and dark.calibration.samples.shape[0] <= m.calibration.samples.shape[0]:
                index[m.artifact_id] = index[dark.artifact_id]
        for m in self.members:
            cal = m.calibration
            idx = index[m.artifact_id]
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
    "NO_VALIDITY_MANY",
    "NOT_CHECKED",
    "SERIES_UNUSED",
    "OTHER_TAGGER",
    "OUTSIDE_VALIDITY",
    "METHOD_LABELLED",
    "METHOD_MINIMUM",
    "OLD_METHOD",
    "OUTDATED_ACCEPTED",
    "TEMPERATURE_TOLERANCE_K",
    "CalibrationSet",
    "DatasetContext",
    "JointDraws",
    "Member",
    "condition_mismatches",
    "check_detectors",
    "detector_mismatches",
    "instrument_mismatches",
    "setting_rtol",
]
