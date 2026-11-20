"""Command-line interface.

    dadl quickstart [chsh-spad | path/to/setup.yaml]   simulate and analyse an example
    dadl init [DIR]                                     copy the example setup.yaml into DIR
    dadl simulate SETUP --out tags.npz [--dark]         write simulated tags in the tag-file format
    dadl analyze SETUP --tags tags.npz|tags.csv         raw CHSH S or HOM visibility from your own tags
                 [--calibrations DIR --allow C=REASON]  HOM and CHSH: corrected by the forward model
    dadl calibrate dark SETUP --tags dark.npz --out DIR fit dark rate, afterpulsing and dead time
    dadl calibrate deadtime SETUP --tags scan.npz --out DIR [--dark DIR]
                                                        series (paralyzable + non-paralyzable) dead time
    dadl calibrate pairs SETUP --tags F --herald A+ --signal B+ --out DIR
                                                        delay, jitter, accidentals, heralding efficiency
    dadl version
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as _datetime
import math
import os
import re
import shutil
import sys
from importlib import resources
from pathlib import Path

import numpy as np

from dadl import __version__
from dadl.analysis.chsh import chsh
from dadl.guardrails import ACCIDENTALS_SUBTRACTED, POLICY_MODIFIED, SIMULATED_INPUT
from dadl.closed_form import chsh_expected_grid, chsh_ideal, hom_visibility_expected
from dadl.config import Setup, load_setup
from dadl.errors import DadlError
from dadl.result import UNCERTAINTY_KIND, LabelledValue
from dadl.simsetup import simulate_setup, simulate_setup_dark, simulate_setup_heralding
from dadl.tagfile import ptu_file, read_tags, tag_tick_ps, tagger_settings, tags_are_simulated, write_ptu_dir, write_tags

SITE = "https://usecommon.github.io/dadl/"
EXAMPLES = ("chsh-spad", "hom-snspd")


def example_path(name: str) -> Path:
    """Path of a bundled example's setup.yaml."""
    if name not in EXAMPLES:
        raise DadlError(f"unknown example {name!r}; available: {', '.join(EXAMPLES)}")
    return Path(str(resources.files("dadl") / "examples" / name / "setup.yaml"))


def _resolve(target: str) -> Path:
    return example_path(target) if target in EXAMPLES else Path(target)


_simulate = simulate_setup


def _print_result(result: LabelledValue, symbol: str = "S") -> None:
    print(f"  Raw {symbol} = {result.raw}   ({UNCERTAINTY_KIND})")
    # Labels about the corrected value are printed with it, not on the raw line (senior-physicist
    # review, cycle 10).
    from dadl.forward import FORWARD_MODEL, chsh_model_label

    hidden = (ACCIDENTALS_SUBTRACTED, FORWARD_MODEL, chsh_model_label(result.corrections))
    raw_labels = [lab for lab in result.labels if lab not in hidden]
    print(f"  Labels: {'; '.join(raw_labels) or 'none'}")
    if result.analysis == "chsh" and result.counts:
        ratios = ", ".join(
            f"{k} {v['coincidences_per_singles_a']:.3f}/{v['coincidences_per_singles_b']:.3f}"
            for k, v in result.counts.items()
            if isinstance(v, dict) and "coincidences_per_singles_a" in v
        )
        print(
            f"  Coincidences per singles, arm A/arm B (what post-selection kept; compare the smaller with a "
            f"detection-loophole threshold): {ratios}"
        )
    w = result.window or {}
    print(f"  Window: {w.get('width_ps')} ps, delays {w.get('delays_ps') or 'none'}, source '{w.get('source')}' (declared, not verified)")
    print(f"  Policy: {result.policy_fingerprint[:12]}{' (modified)' if result.overrides or POLICY_MODIFIED in result.labels else ' (default)'}")


def _warn_angles(setup: Setup) -> None:
    ideal = chsh_ideal(1.0, setup.settings_deg)
    if ideal < 2.0:
        print(
            f"  warning: with the fixed signs S = E(a,b) - E(a,b') + E(a',b) + E(a',b'), the declared angles give at "
            f"most S = {ideal:.3f} for a perfect Phi+ state (2.828 is optimal). Relabelling runs after looking at "
            "the data is a choice outside the plan; set settings_deg to match the sign convention instead."
        )


def _input_hashes(setup_path: Path, tags: str) -> dict[str, str]:
    import hashlib

    out = {}
    if setup_path.is_file():
        out[setup_path.name] = hashlib.sha256(setup_path.read_bytes()).hexdigest()
    path = Path(tags)
    for f in sorted(path.glob("setting*.ptu")) if path.is_dir() else [path]:
        if f.is_file():
            out[f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


def _print_hom_extras(result: LabelledValue, counts: dict, setup: Setup) -> None:
    ratio = counts["duration_dip_s"] / counts["duration_reference_s"] if counts["duration_reference_s"] else 1.0
    if abs(ratio - 1) > 0.01:
        print(
            f"  Exposure: dip {counts['duration_dip_s']:.4g} s, reference {counts['duration_reference_s']:.4g} s; "
            "V compares rates."
        )
    if result.corrections:
        _print_forward(result, counts)
    elif result.corrected is not None:
        print(
            f"  Accidentals subtracted ({counts['accidentals']:.3g} +- {counts['accidentals_stderr']:.2g}, off-peak): "
            f"V = {result.corrected}  (secondary; not an indistinguishability claim)"
        )
    elif counts["accidentals"] > 0:
        print(f"  Accidentals not subtracted (estimate {counts['accidentals']:.3g}); see the labels.")
    ref = abs(setup.reference_delay_ps or 0)
    if ref > 0.1 * setup.window.width_ps:
        effect = (
            "the raw V reads low; the corrected value models the capture."
            if "window_capture" in result.corrections
            else "V reads low."
        )
        print(
            f"  warning: the reference delay ({ref:g} ps) is not small against the {setup.window.width_ps} ps window: "
            f"the fixed grid loses part of the split peak at +-delay, and {effect}"
        )
    coherence = (setup.data.get("source") or {}).get("coherence_time_ps") if isinstance(setup.data, dict) else None
    if coherence and ref < 3 * float(coherence):
        print(
            f"  warning: the reference delay ({ref:g} ps) is less than three coherence times ({float(coherence):g} ps "
            "each): the reference still interferes, and V reads low."
        )


#: Creation and acquisition time of the quickstart's simulated calibrations: fixed, so the artifact
#: IDs it prints are the same on every run.
QUICKSTART_CREATED = _datetime.datetime(2026, 1, 1, tzinfo=_datetime.timezone.utc)
#: Length of the simulated dark run per channel that the quickstart calibrates from.
QUICKSTART_DARK_S = 10.0
#: Length of the simulated pair run (no beamsplitter) that calibrates the delay and jitter.
QUICKSTART_PAIRS_S = 0.2


def _quickstart_hom_corrected(runs: dict, setup: Setup) -> tuple[LabelledValue, dict, str | None]:
    """HOM through the forward model with calibrations from simulated dark and pair runs.

    The dark runs (dead time) and a pair run (delay and jitter, for the window capture) are calibrated into a temporary folder and loaded as an analysis would load
    them. If the dark data only bound the dead time too loosely to use (an SNSPD's 100 cps), the
    refusal is reported and the accidentals are corrected alone.
    """
    import tempfile

    from dadl.analysis.hom import hom_from_tags
    from dadl.artifact import write_artifact
    from dadl.calibration import calibrate_dark
    from dadl.errors import GuardrailRefusal
    from dadl.guardrails import Policy
    from dadl.propagate import CalibrationSet, DatasetContext

    from dadl.calibration.delay import KIND as DELAY_KIND
    from dadl.calibration.delay import calibrate_delay

    dark = simulate_setup_dark(setup, QUICKSTART_DARK_S)[setup.setting_keys[0]]
    pairs = simulate_setup_heralding(setup, QUICKSTART_PAIRS_S)[setup.setting_keys[0]]
    channels = (setup.channel_id("A"), setup.channel_id("B"))
    acquired = QUICKSTART_CREATED
    dataset = DatasetContext(acquired=acquired)
    policy = (
        Policy.default()
        .allow("hom_visibility", "dead_time", "quickstart: dead time from simulated dark runs")
        .allow("hom_visibility", "window_capture", "quickstart: delay and jitter from a simulated pair run")
    )
    with tempfile.TemporaryDirectory() as tmp:
        for name, ch in zip(("A", "B"), channels):
            cal = calibrate_dark(dark[name], ch, backend="sim", seed=setup.seed)
            write_artifact(cal, tmp, accept_misfit=True, created=acquired)
        delay = calibrate_delay(pairs["A"], pairs["B"], channels, backend="sim", seed=setup.seed)
        write_artifact(delay, tmp, accept_misfit=True, created=acquired)
        common = {"channels": channels, "seed": setup.seed, "reference_delay_ps": setup.reference_delay_ps}
        cset = CalibrationSet.load(tmp, analysis="hom_visibility", policy=policy, dataset=dataset)
        try:
            result, counts = hom_from_tags(runs, setup.window, policy=policy, calibrations=cset, **common)
            return result, counts, None
        except GuardrailRefusal as exc:
            note = str(exc).split(", so a correction", 1)[0]
        without_dark = CalibrationSet.load(tmp, analysis="hom_visibility", policy=policy, dataset=dataset, kinds=(DELAY_KIND,))
        result, counts = hom_from_tags(runs, setup.window, policy=policy, calibrations=without_dark, **common)
    return result, counts, note


#: Council round 4: the old note ("a few hundred per second do not identify it") was fixed text
#: whatever the setup's dark rate; the quickstart simply does not calibrate it for CHSH.
QUICKSTART_CHSH_DARK_CPS_NOTE = "a dark calibration is not attempted in the quickstart"


def _quickstart_chsh_corrected(runs: dict, setup: Setup) -> LabelledValue:
    """CHSH with a corrected S from the forward model: each channel pair's delay and jitter from a
    simulated pair run (no analysers), accidentals from the off-peak bins. The policy is loosened
    for both corrections, as any corrected CHSH must be, and the result says so."""
    import tempfile

    from dadl.artifact import write_artifact
    from dadl.calibration.delay import calibrate_delay
    from dadl.guardrails import Policy
    from dadl.propagate import CalibrationSet, DatasetContext
    from dadl.sim import simulate_pairs

    ids = {name: setup.channel_id(name) for name in ("A+", "A-", "B+", "B-")}
    policy = (
        Policy.default()
        .allow("chsh", "accidental_subtraction", "quickstart: accidentals from the off-peak bins")
        .allow("chsh", "window_capture", "quickstart: delay and jitter from simulated pair runs")
    )
    rng = np.random.default_rng(np.random.SeedSequence(setup.seed, spawn_key=(0x434855,)))
    acquired = QUICKSTART_CREATED
    with tempfile.TemporaryDirectory() as tmp:
        for x in "+-":
            for y in "+-":
                a, b = "A" + x, "B" + y
                ta, tb = simulate_pairs(setup.source.pair_rate_per_s, setup.detectors[a], setup.detectors[b], QUICKSTART_PAIRS_S, rng)
                cal = calibrate_delay(ta, tb, (ids[a], ids[b]), backend="sim", seed=setup.seed)
                write_artifact(cal, tmp, accept_misfit=True, created=acquired)
        cset = CalibrationSet.load(tmp, analysis="chsh", policy=policy, dataset=DatasetContext(acquired=acquired))
    return chsh(
        runs, setup.window, policy, settings_deg=setup.settings_deg, calibrations=cset, channels=ids, seed=setup.seed
    )


def _quickstart_hom(path: Path, setup: Setup) -> int:
    runs = _simulate(setup)
    result, counts, dead_time_note = _quickstart_hom_corrected(runs, setup)
    det = setup.detectors["A"]
    src = setup.source
    print(f"DADL {__version__} quickstart: {setup.name} (simulated data, seed {setup.seed})")
    print(f"  Setup: {path}")
    print(
        f"  Source: {src.pair_rate_per_s:.3g} pairs/s, indistinguishability {src.visibility}, coherence time "
        f"{src.coherence_time_ps:g} ps; dip and reference ({setup.reference_delay_ps:g} ps) acquisitions of "
        f"{setup.duration_s_per_setting} s"
    )
    print(
        f"  Detector A{'' if setup.detectors['B'] == det else ' (B differs)'}: efficiency {det.efficiency}, "
        f"dark {det.dark_rate_cps:.0f} cps, dead time {det.dead_time_ps / 1000:g} ns, jitter {det.jitter_sigma_ps:.0f} ps"
    )
    print()
    print(f"  Coincidences: {counts['c_min']:.0f} at the dip, {counts['c_max']:.0f} at the reference")
    _print_result(result, "V")
    _print_hom_extras(result, counts, setup)
    if dead_time_note:
        print(f"  Dead time not corrected: {dead_time_note}")
    expected = hom_visibility_expected(src.visibility, setup.reference_delay_ps, src.coherence_time_ps)
    print(f"  Closed-form expectation (no accidentals, dead time or window losses): V = {expected:.4f}")
    print()
    print("  Edit a copy with `dadl init --example hom-snspd`, then `dadl quickstart setup.yaml`.")
    print(f"  PicoQuant .ptu files: dadl.adapters.read_ptu. Overview: {SITE}")
    return 0


def cmd_quickstart(args: argparse.Namespace) -> int:
    path = _resolve(args.target)
    setup = load_setup(path)
    if setup.analysis_type == "hom":
        return _quickstart_hom(path, setup)
    runs = _simulate(setup)
    result = chsh(runs, setup.window, settings_deg=setup.settings_deg)
    corrected = _quickstart_chsh_corrected(runs, setup)
    det = setup.detectors["A+"]
    same = all(d == det for d in setup.detectors.values())
    print(f"DADL {__version__} quickstart: {setup.name} (simulated data, seed {setup.seed})")
    print(f"  Setup: {path}")
    print(
        f"  Source: {setup.source.pair_rate_per_s:.3g} pairs/s, visibility {setup.source.visibility}, "
        f"{setup.duration_s_per_setting} s per setting"
    )
    print(
        f"  Detectors{'' if same else ' (A+ shown)'}: efficiency {det.efficiency}, dark {det.dark_rate_cps:.0f} cps, "
        f"dead time {det.dead_time_ps / 1000:g} ns, jitter {det.jitter_sigma_ps:.0f} ps, "
        f"afterpulsing {det.afterpulse_prob:.1%}"
    )
    _warn_angles(setup)
    print()
    _print_result(result)
    _print_forward(corrected, corrected.counts or {}, "S")
    print(f"    under a loosened policy {corrected.policy_fingerprint[:12]}: " + "; ".join(corrected.overrides))
    worst = max(d.dead_time_ps * 1e-12 * (setup.source.pair_rate_per_s * d.efficiency / 2 + d.dark_rate_cps) for d in setup.detectors.values())
    print(
        f"  Dead time not corrected: {QUICKSTART_CHSH_DARK_CPS_NOTE}; at these rates about {worst:.2%} of the time is dead."
    )
    if same:
        expected = chsh_expected_grid(
            pair_rate_per_s=setup.source.pair_rate_per_s,
            duration_s=setup.duration_s_per_setting,
            efficiency=det.efficiency,
            dark_rate_cps=det.dark_rate_cps,
            jitter_sigma_ps=det.jitter_sigma_ps,
            width_ps=setup.window.width_ps,
            visibility=setup.source.visibility,
            settings_deg=setup.settings_deg,
        )
        print(f"  Closed-form expectation (ignores dead time, which can raise S, and afterpulsing): S = {expected:.4f}")
    print()
    print("  Edit a copy of this setup with `dadl init`, then `dadl quickstart setup.yaml`.")
    print("  Calibrate from dark data: `dadl simulate chsh-spad --dark --out dark.npz`, then")
    print("  `dadl calibrate dark chsh-spad --tags dark.npz --out artifacts/`. From a heralding run (no")
    print("  analysers): `dadl simulate chsh-spad --heralding --duration-s 2 --out her.npz`, then")
    print("  `dadl calibrate pairs chsh-spad --tags her.npz --no-analysers --herald A+ --signal B+ --out artifacts/`.")
    print(f"  Overview: {SITE}")
    return 0


def cmd_init(args: argparse.Namespace) -> int:
    dest_dir = Path(args.dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / "setup.yaml"
    if dest.exists():
        raise DadlError(f"{dest} already exists; not overwriting it")
    shutil.copyfile(example_path(args.example), dest)
    print(f"Wrote {dest} (from the {args.example} example).")
    print(f"Edit it, then run: dadl quickstart {dest}")
    return 0


def cmd_simulate(args: argparse.Namespace) -> int:
    setup = load_setup(_resolve(args.setup))
    if args.duration_s is not None and not (args.dark or args.heralding or args.deadtime_scan):
        raise DadlError(
            "--duration-s applies to --dark, --heralding and --deadtime-scan; a CHSH or HOM run takes its duration "
            "from the setup"
        )
    if args.rates is not None and not args.deadtime_scan:
        raise DadlError("--rates applies to --deadtime-scan")
    default_s = 0.05 if args.deadtime_scan else 20.0
    args.duration_s = default_s if args.duration_s is None else args.duration_s
    if not args.duration_s > 0:
        raise DadlError("--duration-s must be positive")
    if sum(map(bool, (args.dark, args.heralding, args.deadtime_scan))) > 1:
        raise DadlError("choose one of --dark, --heralding and --deadtime-scan")
    if args.deadtime_scan:
        from dadl.simsetup import deadtime_scan_rates, simulate_setup_deadtime_scan

        rates = _rates(args.rates) if args.rates is not None else deadtime_scan_rates(setup)
        runs = simulate_setup_deadtime_scan(setup, rates, args.duration_s)
        shown = ", ".join(f"{r:.3g}" for r in rates)
        what = f"dead-time scan tags ({args.duration_s:g} s per rate; photon rates {shown} per s, settings 0-{len(rates) - 1})"
    elif args.dark:
        runs = simulate_setup_dark(setup, args.duration_s)
        what = f"dark tags ({args.duration_s:g} s per channel, setting 0)"
    elif args.heralding:
        runs = simulate_setup_heralding(setup, args.duration_s)
        what = f"heralding-run tags ({args.duration_s:g} s, no analysers or beamsplitter, setting 0)"
    else:
        runs = simulate_setup(setup)
        what = "HOM tags (setting 0 dip, 1 reference)" if setup.analysis_type == "hom" else "CHSH tags"
    # Security review, cycle 9: simulate replaced recorded PTU files and wrote through symlinks.
    out_path = Path(args.out)
    if args.format != "ptu" and args.device is not None:
        raise DadlError("--device applies only with --format ptu")
    args.device = args.device or "hydraharp"
    if args.format == "ptu" and out_path.suffix in (".npz", ".csv"):
        raise DadlError(f"--format ptu writes a folder of PTU files; {args.out} looks like a tag file name")
    targets = [ptu_file(out_path, k) for k in range(len(setup.setting_keys))] if args.format == "ptu" else [out_path]
    existing = [str(t) for t in targets if os.path.lexists(t)]
    if existing and not args.force:
        raise DadlError(f"{', '.join(existing[:3])} already exist(s); pass --force to overwrite them with simulated data")
    if args.force:
        for t in existing:
            if Path(t).is_symlink():
                Path(t).unlink()
    if args.format == "ptu":
        try:
            files = write_ptu_dir(args.out, runs, setup.channels, setup.setting_keys, device=args.device)
        except ValueError as exc:
            raise DadlError(str(exc)) from exc
        print(f"Wrote simulated {what} for {setup.name} as {len(files)} {args.device} PTU files in {args.out}")
        print("  The files are marked as simulated (CreatorSW_Name DADL, DADL_Simulated) and read back as such.")
        return 0
    out = write_tags(args.out, runs, setup.channels, setup.setting_keys, simulated=True)
    print(f"Wrote simulated {what} for {setup.name} to {out}")
    return 0


def _rates(text: str) -> tuple[float, ...]:
    try:
        return tuple(float(x) for x in text.split(","))
    except ValueError as exc:
        raise DadlError(f"--rates: comma-separated photon rates per s, got {text!r}") from exc


def _add_record_options(p: argparse.ArgumentParser) -> None:
    """Options that complete the artifact as a record (calibration-lab review, cycle 3)."""
    p.add_argument("--condition", action="append", default=[], metavar="KEY=VALUE",
                   help="a numeric condition to record, e.g. temperature_K=233.1 or bias_V=24.5 (repeatable)")
    p.add_argument("--acquired", help="when the data were taken (ISO 8601 with a UTC offset); a PTU file's local creation time is recorded separately")
    p.add_argument("--created", help="record creation time (ISO 8601), to reproduce an artifact; default now")
    p.add_argument(
        "--valid-until",
        help="end of the intended validity (ISO 8601); an analysis acquired outside it is refused unless "
        "calibration_outside_validity is allowed",
    )
    p.add_argument("--supersedes", help="ID of the artifact this one replaces (must exist in --out)")


def _when(text: str | None, name: str):
    import datetime as dt

    if text is None:
        return None
    try:
        t = dt.datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as exc:
        raise DadlError(f"--{name}: not an ISO 8601 time: {text!r}") from exc
    if t.tzinfo is None:
        # Archivist review, cycle 8: a local time stored as UTC breaks matching with other logs.
        raise DadlError(f"--{name}: give a UTC offset (e.g. {text}Z or {text}+01:00)")
    return t


def _record(args: argparse.Namespace, setup_path: Path) -> tuple[dict, dict[str, float], tuple]:
    """Write options, extra conditions and extra inputs for this run's artifacts."""
    import datetime as dt
    import hashlib

    conditions = _conditions(args.condition)
    conditions.update(tagger_settings(args.tags))
    tags_path = Path(args.tags)
    files = sorted(tags_path.glob("setting*.ptu")) if tags_path.is_dir() else [tags_path]
    inputs = tuple((hashlib.sha256(f.read_bytes()).hexdigest(), "raw-file", f.name) for f in files if f.is_file())
    if setup_path.is_file():
        inputs += ((hashlib.sha256(setup_path.read_bytes()).hexdigest(), "config", setup_path.name),)
    acquired, instrument = _when(args.acquired, "acquired"), {}
    if tags_path.is_dir() and files:
        from dadl.adapters.ptu import PtuFile

        header = PtuFile(files[0]).header
        instrument = {"model": header.get("HW_Type"), "serial": header.get("HW_SerialNo"),
                      "firmware": header.get("HW_Version")}
        stamp = header.get("File_CreatingTime")
        if isinstance(stamp, float) and stamp > 0:
            # The acquisition PC's local clock, zone unknown: recorded as such, not as UTC.
            local = dt.datetime(1899, 12, 30) + dt.timedelta(days=stamp)
            instrument["file_created_local"] = local.strftime("%Y-%m-%dT%H:%M:%S") + " (local time, zone unknown)"
    supersedes = args.supersedes
    if supersedes and not re.fullmatch(r"[0-9a-f]{64}", supersedes):
        raise DadlError(f"--supersedes must be an artifact ID (64 hex digits), got {supersedes!r}")
    if supersedes and not (Path(args.out) / f"{supersedes}.json").exists():
        raise DadlError(f"--supersedes: no artifact {supersedes} in {args.out}")
    write = {
        "created": _when(args.created, "created"),
        "acquired": acquired,
        "valid_until": _when(args.valid_until, "valid-until"),
        "supersedes": supersedes,
        "instrument": {k: str(v) for k, v in instrument.items() if v not in (None, "")},
    }
    return write, conditions, inputs


def _write_checked(cal, out: str, accept_misfit: bool, **write) -> str:
    """Write the artifact, refusing a misfit unless the user accepts it (it is then marked misfit)."""
    from dadl.artifact import write_artifact

    if cal.status == "misfit" and not accept_misfit:
        raise DadlError(
            f"{cal.kind} for {', '.join(map(str, cal.channels))}: the model does not describe these data "
            f"(predictive check p = {cal.ppc_pvalue:.3g}). No artifact written. Inspect the data (drift, latching, "
            "a timing tail, afterpulse structure, a pulsed source), or pass --accept-misfit to keep it marked misfit."
        )
    return str(write_artifact(cal, out, accept_misfit=accept_misfit, **write))


def _backend(setup, tags_path: str) -> str:
    """"sim" for tags DADL marked as simulated, else "hardware".

    A setup still saying `backend: sim` (as `dadl init` copies it) with tags that are not marked
    simulated is refused: the new-user review found recorded data written up as simulated, with
    the channel named "sim:0" instead of the tagger.
    """
    if tags_are_simulated(tags_path):
        return "sim"
    if setup.backend == "sim":
        raise DadlError(
            f"{tags_path} is not marked as simulated (DADL marks the .npz and PTU files it simulates), but the "
            "setup says `backend: sim`. For recorded data set `backend: file` and `hardware: {tagger: <serial>}` "
            "in setup.yaml. (Files written by `dadl simulate` are marked; a file re-saved by other software loses "
            "the mark.)"
        )
    path = Path(tags_path)
    if path.is_dir() and setup.tagger:
        from dadl.adapters.ptu import PtuFile

        serials = {str(PtuFile(f).header.get("HW_SerialNo", "")).strip() for f in sorted(path.glob("setting*.ptu"))}
        serials.discard("")
        # Whole-token comparison: "123" must not match "HydraHarp-1234567" (twin review, cycle 9).
        tagger = setup.tagger
        wrong = sorted(sn for sn in serials if not (sn == tagger or tagger.endswith(("-" + sn, " " + sn, "_" + sn))))
        if wrong:
            raise DadlError(
                f"the PTU files were recorded on tagger serial {', '.join(wrong)}, but setup.yaml says "
                f"hardware.tagger: {setup.tagger}"
            )
    return "hardware"


def _check_time_units(runs: dict, setup) -> None:
    """Refuse tags whose gaps are impossible in picoseconds for the declared dead times.

    Teaching-lab review, cycle 4: tags exported in ns and read as ps gave S = 2.23 +- 0.01 with a
    window really 1000x wider; 43% of the gaps were shorter than the declared dead time.
    """
    for key, tags in runs.items():
        for name, arr in tags.items():
            det = setup.detectors.get(name)
            arr = np.asarray(arr)
            if arr.size > 1 and np.any(np.diff(arr) == 0):
                # Red-team review, cycle 5: a duplicated row hid the dead-time-selection label.
                raise DadlError(
                    f"{name} in acquisition {key}: duplicate timestamps on one channel (a detector cannot click "
                    "twice at once); check how the file was spliced or exported"
                )
            if det is None or det.dead_time_ps <= 0 or arr.size < 100:
                continue
            short = float(np.mean(np.diff(arr) < det.dead_time_ps / 2))
            # Recorded ps data have essentially no gaps under half the dead time; 0.2% catches ns
            # exports of dark data too (second red-team review, cycle 9: 0.5% passed at 5%).
            if short > 0.002:
                raise DadlError(
                    f"{name} in acquisition {key}: {short:.1%} of the gaps between tags are shorter than half the "
                    f"declared dead time ({det.dead_time_ps} ps). Are the timestamps in picoseconds (not ns)? If "
                    "they are, the channel has double triggers or echoes (a cable reflection, a discriminator "
                    "ringing): fix the hardware, or remove them and say so; do not lower dead_time_ps to pass."
                )


SOURCE_ON_RANGE_PS = 50_000


def _refuse_source_on(tags: dict, setup) -> None:
    """Refuse a "dark" acquisition whose channels share coincidences: the source was on.

    The new-user review ran `calibrate dark` on source-on CHSH tags and got a confident dark rate
    of 25,400 cps (true 500), written as a normal artifact.
    """
    from dadl.histogram import cross_correlation

    # Searched over +-50 ns, not only at the declared delays: dark calibration usually comes before
    # the delays are declared (second red-team review, cycle 9: a 3 ns cable delay hid the peak).
    names = [n for n in sorted(tags) if np.asarray(tags[n]).size > 1]
    width = setup.window.width_ps
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            a, b = np.asarray(tags[x]), np.asarray(tags[y])
            span = float(max(a[-1], b[-1]) - min(a[0], b[0]))
            if span <= 0:
                continue
            hist = cross_correlation(a, b, width, SOURCE_ON_RANGE_PS).astype(float)
            expected = a.size * b.size * width / span
            peak = float(hist.max())
            n_bins = hist.size
            # 6 sigma plus a look-elsewhere margin over the bins searched.
            if peak - expected > 10 and (peak - expected) / np.sqrt(expected + 1.0) > 6 + math.log(n_bins):
                raise DadlError(
                    f"channels {x} and {y} share a coincidence peak of {peak:.0f} in one {width} ps bin where "
                    f"{expected:.3g} are expected by chance: the source was on. A dark calibration needs an "
                    "acquisition with the source off (simulated: dadl simulate SETUP --dark)."
                )


def cmd_calibrate_dark(args: argparse.Namespace) -> int:
    from dadl.calibration import calibrate_dark

    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels, setup.setting_keys)
    keys = setup.setting_keys
    if not 0 <= args.setting < len(keys):
        raise DadlError(f"--setting must be 0..{len(keys) - 1}")
    tags = runs[keys[args.setting]]
    number = {name: num for num, name in setup.channels.items()}
    wanted = args.channel or [setup.channels[n] for n in sorted(setup.channels)]
    unknown = [c for c in wanted if c not in number]
    if unknown:
        raise DadlError(f"unknown channel(s) {unknown}; the setup's channels are {sorted(number)}")
    backend = _backend(setup, args.tags)
    _check_time_units(runs, setup)
    _refuse_source_on(tags, setup)
    print(f"DADL {__version__} calibrate dark: {setup.name}, tags from {args.tags} (setting {args.setting})")
    write, extra, inputs = _record(args, _resolve(args.setup))
    # Each channel is calibrated and refused on its own, and a failure names its channel (operator
    # review, cycle 8: the first bad channel stopped the rest).
    failed: list[str] = []
    for name in wanted:
        arr = tags[name]
        run = {"setting": float(args.setting)}
        if arr.size > 1:
            run["duration_s"] = (int(arr[-1]) - int(arr[0])) / 1e12
        channel = setup.channel_id(name)  # a setup error stops the command
        try:
            cal = calibrate_dark(arr, channel, backend=backend, seed=setup.seed, conditions={**run, **extra})
        except (DadlError, ValueError) as exc:
            print()
            print(f"  {name}: not calibrated: {exc}")
            failed.append(name)
            continue
        cal = dataclasses.replace(cal, inputs=cal.inputs + inputs)
        print()
        print(cal)
        dead_declared = setup.detectors[name].dead_time_ps if name in setup.detectors else 0
        dead_high = cal.estimates["dead_time_ps"].high
        if dead_declared > 0 and dead_high < 0.8 * dead_declared:
            print(
                f"  warning: the dead time (at most {dead_high:.0f} ps) is well below the setup's declared "
                f"{dead_declared} ps. A spurious double trigger or reflection can set it; check the shortest gaps."
            )
        declared = setup.detectors[name].dark_rate_cps if name in setup.detectors else 0.0
        fitted = cal.estimates["dark_rate_cps"].low
        if declared > 0 and fitted > 100 * declared:
            # Second red-team review, cycle 9: dark tags exported in ns read 500,613 cps for 500, ok.
            print(
                f"  {name}: not written: the dark rate is more than 100x the setup's declared {declared:g} cps. "
                "Timestamps in ns rather than ps, light on the detector, or a wrong declared rate?"
            )
            failed.append(name)
            continue
        if declared > 0 and fitted > 10 * declared:
            print(
                f"  warning: the dark rate is more than 10x the setup's declared {declared:g} cps. Check that the "
                "source was off and the detector shielded from light."
            )
        try:
            path = _write_checked(cal, args.out, args.accept_misfit, **write)
        except DadlError as exc:
            print(f"  not written: {exc}")
            failed.append(name)
        else:
            print(f"  artifact: {path}")
    if failed:
        raise DadlError(f"{len(failed)} of {len(wanted)} channels not calibrated or not written: {', '.join(failed)}")
    return 0


def _dark_artifact_for(source: str, channel):
    """The channel's dark artifact from a folder or file, or a refusal naming what is there."""
    from dadl.artifact import load_calibration
    from dadl.calibration.dark import KIND as DARK_KIND

    path = Path(source)
    files = sorted(path.glob("*.json")) if path.is_dir() else [path]
    found = []
    for f in files:
        try:
            cal = load_calibration(f)
        except (DadlError, ValueError, OSError):
            continue
        if cal.kind == DARK_KIND and tuple(cal.channels) == (channel,):
            found.append((f, cal))
    if not found:
        raise DadlError(f"no dark artifact for {channel} in {source}")
    if len(found) > 1:
        raise DadlError(f"{len(found)} dark artifacts for {channel} in {source}; pass the one to use as a file")
    return found[0]


def cmd_calibrate_deadtime(args: argparse.Namespace) -> int:
    from dadl.calibration.deadtime import calibrate_deadtime_series

    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels, setup.setting_keys)
    keys = setup.setting_keys
    settings = args.settings if args.settings else [
        i for i, k in enumerate(keys) if any(np.asarray(t).size > 1 for t in runs[k].values())
    ]
    bad = [i for i in settings if not 0 <= i < len(keys)]
    if bad:
        raise DadlError(f"--settings must be in 0..{len(keys) - 1}")
    number = {name: num for num, name in setup.channels.items()}
    wanted = args.channel or [setup.channels[n] for n in sorted(setup.channels)]
    unknown = [c for c in wanted if c not in number]
    if unknown:
        raise DadlError(f"unknown channel(s) {unknown}; the setup's channels are {sorted(number)}")
    backend = _backend(setup, args.tags)
    _check_time_units(runs, setup)
    print(f"DADL {__version__} calibrate deadtime: {setup.name}, tags from {args.tags} (settings {settings})")
    write, extra, inputs = _record(args, _resolve(args.setup))
    failed: list[str] = []
    for name in wanted:
        channel = setup.channel_id(name)
        acquisitions = [runs[keys[i]][name] for i in settings]
        afterpulse = None
        dark_input: tuple = ()
        try:
            if args.dark:
                from dadl.artifact import artifact_id, read_artifact

                dark_path, afterpulse = _dark_artifact_for(args.dark, channel)
                dark_input = ((artifact_id(read_artifact(dark_path)[0]), "previous-artifact", str(channel)),)
                print(f"  {name}: afterpulsing from {dark_path}")
            cal = calibrate_deadtime_series(
                acquisitions, channel, afterpulse=afterpulse, backend=backend, seed=setup.seed,
                conditions={"settings": float(len(settings)), **extra},
            )
        except (DadlError, ValueError) as exc:
            print()
            print(f"  {name}: not calibrated: {exc}")
            failed.append(name)
            continue
        cal = dataclasses.replace(cal, inputs=cal.inputs + dark_input + inputs)
        print()
        print(cal)
        try:
            path = _write_checked(cal, args.out, args.accept_misfit, **write)
        except DadlError as exc:
            print(f"  not written: {exc}")
            failed.append(name)
        else:
            print(f"  artifact: {path}")
    if failed:
        raise DadlError(f"{len(failed)} of {len(wanted)} channels not calibrated or not written: {', '.join(failed)}")
    return 0


def cmd_calibrate_pairs(args: argparse.Namespace) -> int:
    from dadl.calibration import calibrate_accidentals, calibrate_delay, calibrate_heralding

    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels, setup.setting_keys)
    keys = setup.setting_keys
    if not 0 <= args.setting < len(keys):
        raise DadlError(f"--setting must be 0..{len(keys) - 1}")
    _check_time_units(runs, setup)
    tags = runs[keys[args.setting]]
    if setup.analysis_type in ("chsh", "hom") and not args.no_analysers:
        raise DadlError(
            f"{setup.name} is a {setup.analysis_type} setup: its acquisitions pass through analysers or a "
            "beamsplitter, which enter the heralding ratio directly (a Bell setting can read 0.06 for an "
            "efficiency of 0.5). Record a heralding run with the source straight onto the detectors and pass "
            "--no-analysers to confirm it (simulated: dadl simulate SETUP --heralding)."
        )
    number = {name: num for num, name in setup.channels.items()}
    for name in (args.herald, args.signal):
        if name not in number:
            raise DadlError(f"unknown channel {name!r}; the setup's channels are {sorted(number)}")
    if args.herald == args.signal:
        raise DadlError("--herald and --signal must be different channels")
    chs = (setup.channel_id(args.herald), setup.channel_id(args.signal))
    a, b = tags[args.herald], tags[args.signal]
    window = setup.window
    d_a, d_b = window.delays_ps.get(args.herald, 0), window.delays_ps.get(args.signal, 0)
    backend = _backend(setup, args.tags)
    span_s = (max(int(a[-1]), int(b[-1])) - min(int(a[0]), int(b[0]))) / 1e12 if a.size and b.size else 0.0
    run_conditions = {
        "delay_herald_ps": float(d_a),
        "delay_signal_ps": float(d_b),
        "duration_s": span_s,
        "rate_herald_cps": a.size / span_s if span_s else 0.0,
        "rate_signal_cps": b.size / span_s if span_s else 0.0,
        "analysers_in_path": 0.0,
        "setting": float(args.setting),
    }
    write, extra, inputs = _record(args, _resolve(args.setup))
    run_conditions.update(extra)
    common = {"backend": backend, "seed": setup.seed, "conditions": run_conditions}
    print(f"DADL {__version__} calibrate pairs: {setup.name}, {args.herald} -> {args.signal}, tags from {args.tags}")
    print(f"  Window {window.width_ps} ps and delays {d_a}/{d_b} ps as declared in the setup.")
    # The declared detectors' dead time and afterpulse time constant set the off-peak guard and the
    # dead-time notes (source-physics review, cycle 7: the CLI never passed them).
    dets = [setup.detectors[n] for n in (args.herald, args.signal) if n in setup.detectors]
    dead = max((float(d.dead_time_ps) for d in dets), default=0.0) or None
    tau_a = max((float(d.afterpulse_tau_ps) for d in dets if d.afterpulse_prob > 0), default=0.0) or None
    if dead or tau_a:
        print(f"  Off-peak guard and dead-time notes use the setup's declared dead time ({dead or 0:g} ps) and "
              f"afterpulse time constant ({tau_a or 0:g} ps).")
    detector_args = {"dead_time_ps": dead, "afterpulse_tau_ps": tau_a}
    cals = [
        calibrate_delay(
            a, b, chs, range_ps=args.range_ps, bin_width_ps=args.bin_width_ps, tick_ps=tag_tick_ps(args.tags), **common
        ),
        calibrate_accidentals(
            a, b, chs, width_ps=window.width_ps, delay_a_ps=d_a, delay_b_ps=d_b, **detector_args, **common
        ),
        calibrate_heralding(
            a, b, chs, width_ps=window.width_ps, delay_herald_ps=d_a, delay_signal_ps=d_b, **detector_args, **common
        ),
    ]
    # Each artifact is checked and refused on its own: a timing tail that makes the delay fit a
    # misfit says nothing about the accidentals or the heralding ratio (SNSPD specialist review).
    refused = []
    fitted = cals[0].estimates["delay_ps"].value
    off_by = fitted - (d_b - d_a)
    misplaced = abs(off_by) > window.width_ps / 4
    if misplaced:
        # Software review, cycle 2: a delay declared 15 ns off wrote heralding_efficiency 0 with
        # status ok. The heralding ratio needs the pairs in the declared slot.
        note = (
            f"the declared delay differs from the fitted one by {off_by:.0f} ps (more than a quarter window): "
            "true pairs are not in the declared slot"
        )
        cals[1] = dataclasses.replace(cals[1], notes=(*cals[1].notes, note))
        cals[2] = dataclasses.replace(cals[2], notes=(*cals[2].notes, note))
    for cal in cals:
        print()
        print(cal)
        if misplaced and cal is cals[2]:
            # Marked misfit, so --accept-misfit writes it as such rather than as ok (second
            # red-team review, cycle 9).
            cal = dataclasses.replace(cal, ppc_pvalue=0.0)
        try:
            path = _write_checked(dataclasses.replace(cal, inputs=cal.inputs + inputs), args.out, args.accept_misfit, **write)
        except DadlError as exc:
            refused.append(cal.kind)
            print(f"  not written: {exc}")
        else:
            print(f"  artifact: {path}")
    if misplaced:
        print()
        print(
            f"  The fitted delay ({fitted:.0f} ps) differs from the declared one ({d_b - d_a} ps) by more than a "
            "quarter window. Declare it in the plan before analysis; DADL never moves the window to fit the data."
        )
    from dadl.forward import grid_capture

    # The forward model's closed form (council round 4: a quad integral read about 0 for a narrow
    # peak in a wide window).
    capture = float(grid_capture(fitted - (d_b - d_a), cals[0].estimates["jitter_sigma_ps"].value, window.width_ps))
    eta = cals[2].estimates["heralding_efficiency"].value
    print()
    if cals[0].status == "misfit":
        # Source-physics review, cycle 7: a Gaussian capture from a rejected fit read 7% low.
        print("  Window capture not estimated: the delay fit's Gaussian timing model does not describe these data.")
    elif capture > 0.01:
        print(
            f"  Window capture: about {capture:.3f} of true pairs fall in the same {window.width_ps} ps slot, from the "
            f"fitted jitter and delay offset. heralding_efficiency / capture = {eta / capture:.3f}, a figure closer to "
            "a sliding-window or Klyshko efficiency (an estimate, not recorded in the artifact; dead-time losses "
            "and dilution by herald dark counts and afterpulses are still included)."
        )
    else:
        print(f"  Window capture: {capture:.3g}; the declared window and delays miss the pair peak.")
    from dadl.coincidences import grid_coincidences

    n_c = grid_coincidences(a, b, window.width_ps, d_a, d_b)
    acc = cals[1].estimates["accidental_coincidences"].value
    if n_c > 0 and span_s > 0:
        car = f"{n_c / acc:.4g}" if acc > 0 else "no accidentals seen"
        rate = a.size * b.size / n_c / span_s
        if cals[0].status != "misfit" and capture > 0.01:
            # Senior-physicist review, cycle 10: without the window capture this read 27% high.
            rate_text = f"{rate * capture:.4g} per s (singles x singles x capture / coincidences"
        else:
            rate_text = f"{rate:.4g} per s (singles x singles / coincidences, high by 1/capture"
        print(
            f"  Source figures (secondary, not recorded): coincidences-to-accidentals {car}; pair rate {rate_text}; "
            "dark counts in the singles also raise it)."
        )
    if refused:
        raise DadlError(f"{len(refused)} of {len(cals)} artifacts not written (misfit): {', '.join(refused)}")
    return 0


def _conditions(items: list[str]) -> dict[str, float]:
    conditions: dict[str, float] = {}
    for item in items:
        key, sep, value = item.partition("=")
        try:
            conditions[key.strip()] = float(value)
        except ValueError:
            sep = ""
        if not sep or not key.strip():
            raise DadlError(f"--condition needs KEY=NUMBER, got {item!r}")
    return conditions


def _policy(analysis: str, allow: list[str]):
    from dadl.guardrails import Policy

    policy = Policy.default()
    for item in allow:
        correction, sep, reason = item.partition("=")
        if not sep or not correction.strip() or not reason.strip():
            raise DadlError(f"--allow needs CORRECTION=REASON (a reason is recorded), got {item!r}")
        policy = policy.allow(analysis, correction.strip(), reason.strip())
    return policy


def _calibration_set(args: argparse.Namespace, runs: dict, setup_path: Path, analysis: str, policy, channels):
    """Artifacts from --calibrations, checked against this dataset (time, conditions, raw inputs) and
    the setup's channels (detector serials)."""
    import hashlib

    from dadl.data import tags_sha256
    from dadl.propagate import CalibrationSet, DatasetContext

    hashes = {tags_sha256(np.asarray(a)) for run in runs.values() for a in run.values()}
    hashes |= set(_input_hashes(setup_path, args.tags).values())
    tags_path = Path(args.tags)
    if tags_path.is_file():
        hashes.add(hashlib.sha256(tags_path.read_bytes()).hexdigest())
    conditions = {**tagger_settings(args.tags), **_conditions(args.condition)}
    dataset = DatasetContext(acquired=_when(args.acquired, "acquired"), conditions=conditions, input_sha256=frozenset(hashes))
    return CalibrationSet.load(args.calibrations, analysis=analysis, policy=policy, dataset=dataset, channels=channels)


def _print_forward(result: LabelledValue, counts: dict, symbol: str = "V") -> None:
    """The corrected line of a forward-model result: value, corrections in order, artifacts."""
    order = ", then ".join(c.replace("_", " ") for c in result.corrections)
    note = "  (secondary; not a Bell claim)" if result.analysis == "chsh" else ""
    print(f"  Corrected (forward model: {order}): {symbol} = {result.corrected}{note}")
    if "corrected_at_zero_dead_time" in counts:
        print(
            f"  Dead time known only as an upper limit: {symbol} = {counts['corrected_at_zero_dead_time']:.4f} at zero, "
            f"{counts['corrected_at_bound_dead_time']:.4f} at the bound (the corrected value and interval use the bound)"
        )
    if result.artifact_ids:
        print("  Calibrations: " + ", ".join(a[:12] for a in result.artifact_ids))


def cmd_analyze(args: argparse.Namespace) -> int:
    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels, setup.setting_keys)
    _check_time_units(runs, setup)
    print(f"DADL {__version__} analyze: {setup.name}, tags from {args.tags}")
    # Red-team review, cycle 5: simulated tags gave an unmarked headline, and nothing tied the result
    # to the setup and tags it read.
    extra = {"inputs": _input_hashes(_resolve(args.setup), args.tags)}
    simulated = tags_are_simulated(args.tags)
    analysis = "hom_visibility" if setup.analysis_type == "hom" else "chsh"
    policy = _policy(analysis, args.allow)
    names = ("A", "B") if setup.analysis_type == "hom" else ("A+", "A-", "B+", "B-")
    used = [setup.channel_id(name) for name in names]
    cset = _calibration_set(args, runs, _resolve(args.setup), analysis, policy, used) if args.calibrations else None
    if setup.analysis_type == "hom":
        from dadl.analysis.hom import hom_from_tags

        channels = (setup.channel_id("A"), setup.channel_id("B")) if cset is not None else None
        result, counts = hom_from_tags(
            runs, setup.window, policy=policy, calibrations=cset, channels=channels, seed=setup.seed,
            reference_delay_ps=setup.reference_delay_ps,
        )
        result = dataclasses.replace(
            result, labels=result.labels + ((SIMULATED_INPUT,) if simulated else ()), **extra
        )
        print("  The tags are trusted as given.")
        print()
        print(f"  Coincidences: {counts['c_min']:.0f} at the dip, {counts['c_max']:.0f} at the reference")
        _print_result(result, "V")
        _print_hom_extras(result, counts, setup)
        return 0
    ids = {name: setup.channel_id(name) for name in ("A+", "A-", "B+", "B-")} if cset is not None else None
    result = chsh(
        runs, setup.window, settings_deg=setup.settings_deg, policy=policy, calibrations=cset, channels=ids, seed=setup.seed
    )
    result = dataclasses.replace(result, labels=result.labels + ((SIMULATED_INPUT,) if simulated else ()), **extra)
    _warn_angles(setup)
    print("  The tags are trusted as given; the headline is raw S." if cset is not None else "  The tags are trusted as given; raw S only.")
    print()
    _print_result(result)
    if result.corrected is not None:
        _print_forward(result, result.counts or {}, "S")
    print("  Inputs: " + ", ".join(f"{k} {v[:12]}" for k, v in result.inputs.items()))
    print(
        "  Choosing the window, delays or run after seeing the data inflates false violations (best of 72 "
        "window and delay choices: 8% at a true S of 2); declare them before acquisition."
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dadl",
        description="Detector-Aware Data Layer for quantum optics experiments.",
        epilog=f"Overview and docs: {SITE}",
    )
    parser.add_argument("--version", action="version", version=f"dadl {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    qs = sub.add_parser(
        "quickstart",
        help="simulate an example (or your setup.yaml) and print its raw result (CHSH S or HOM visibility), and a "
        "corrected value from simulated calibrations",
        description="Simulate an experiment (CHSH or HOM) from a setup file and analyse it with a fixed, planned window.",
    )
    qs.add_argument(
        "target",
        nargs="?",
        default="chsh-spad",
        help=f"a bundled example ({', '.join(EXAMPLES)}) or a path to a setup.yaml (default: chsh-spad)",
    )
    qs.set_defaults(func=cmd_quickstart)

    init = sub.add_parser("init", help="copy an example setup.yaml into a folder to edit")
    init.add_argument("dir", nargs="?", default=".", help="destination folder (default: current folder)")
    init.add_argument("--example", default="chsh-spad", choices=EXAMPLES)
    init.set_defaults(func=cmd_init)

    sim = sub.add_parser("simulate", help="write simulated tags from a setup in the tag-file format")
    sim.add_argument("setup", help="setup.yaml or a bundled example name")
    sim.add_argument("--out", required=True, help="output file (.npz or .csv), or a folder with --format ptu")
    sim.add_argument(
        "--format", choices=["tags", "ptu"], default="tags",
        help="tags: one generic tag file (default); ptu: one PicoQuant T2 file per acquisition (hardware twin)",
    )
    from dadl.adapters.ptu_write import DEVICES

    sim.add_argument("--device", choices=sorted(DEVICES), default=None, help="PTU device to mimic (default hydraharp)")
    sim.add_argument("--force", action="store_true", help="overwrite existing output files")
    sim.add_argument("--dark", action="store_true", help="simulate dark acquisitions (source off) for calibration")
    sim.add_argument(
        "--heralding", action="store_true",
        help="simulate a heralding run: the pair source straight onto one detector per arm (for calibrate pairs)",
    )
    sim.add_argument(
        "--deadtime-scan", action="store_true",
        help="steady light on every detector at several rates, one per setting (for calibrate deadtime)",
    )
    sim.add_argument(
        "--rates", default=None, metavar="R1,R2,...",
        help="--deadtime-scan: photon rate at each detector (before efficiency) per s, one per setting (default: "
        "avalanche rate x dead time from 0.05 to 0.5 over four settings; 0.1 and 0.4 for two)",
    )
    sim.add_argument(
        "--duration-s", type=float, default=None,
        help="acquisition length for --dark (per channel) or --heralding (default 20), or per rate for --deadtime-scan "
        "(default 0.05); a CHSH or HOM run uses the setup",
    )
    sim.set_defaults(func=cmd_simulate)

    an = sub.add_parser("analyze", help="raw CHSH S or HOM visibility from your own tag file (corrected too with "
                        "--calibrations)")
    an.add_argument("setup", help="setup.yaml giving the window, delays and channel map")
    an.add_argument("--tags", required=True, help="tag file (.npz or .csv: columns setting, channel, timestamp_ps) or a folder of PTU files setting0.ptu, ...")
    an.add_argument(
        "--calibrations", action="append", metavar="DIR",
        help="calibration artifacts (folder or file, repeatable): HOM and CHSH gain a corrected value from the "
        "forward model",
    )
    an.add_argument("--acquired", help="when the tags were taken (ISO 8601 with a UTC offset), checked against each artifact's validity")
    an.add_argument("--condition", action="append", default=[], metavar="KEY=VALUE",
                    help="a numeric condition of this acquisition (repeatable), checked against each artifact's conditions")
    an.add_argument(
        "--allow", action="append", default=[], metavar="CORRECTION=REASON",
        help="loosen the default policy for this analysis (e.g. dead_time='calibrated from dark runs'); labelled 'policy modified'",
    )
    an.set_defaults(func=cmd_analyze)

    cal = sub.add_parser("calibrate", help="estimate detector parameters and write calibration artifacts")
    cal_sub = cal.add_subparsers(dest="calibration", required=True)
    dark = cal_sub.add_parser(
        "dark",
        help="dark rate, afterpulsing and dead time from a dark acquisition",
        description="Fit each channel's dark rate, afterpulse probability and time constant, and dead time "
        "from the gaps between dark counts. Writes one content-addressed artifact per channel.",
    )
    dark.add_argument("setup", help="setup.yaml (channel map) or a bundled example name")
    dark.add_argument("--tags", required=True, help="tag file, or a folder of PTU files, with the dark acquisition")
    dark.add_argument("--setting", type=int, default=0, help="which setting index holds the dark data (default 0)")
    dark.add_argument("--channel", action="append", help="channel name to calibrate (repeatable; default all)")
    dark.add_argument("--out", required=True, help="directory for the artifacts")
    dark.add_argument(
        "--accept-misfit", action="store_true",
        help="write artifacts whose predictive check fails (p < 0.001), marked misfit; refused otherwise",
    )
    _add_record_options(dark)
    dark.set_defaults(func=cmd_calibrate_dark)

    dt = cal_sub.add_parser(
        "deadtime",
        help="series dead time (paralyzable then non-paralyzable) from illuminated runs at several rates",
        description="Fit each channel's paralyzable (dead_time_p_ps) and non-paralyzable (dead_time_ps) dead times "
        "from the gaps of two or more acquisitions under steady light at different rates (simulated: dadl simulate "
        "SETUP --deadtime-scan). Afterpulsing is held at the channel's dark artifact (--dark). Writes one "
        "content-addressed artifact per channel, of kind deadtime-series. The scan needs steady light: an "
        "attenuated CW laser with its intensity stable to about 1-2% over each run, not pulsed light or a "
        "mains-powered lamp. Slow intensity noise looks like less paralysis: a +-20% 100 Hz modulation biased "
        "tau_p 10% low with the predictive check passing; +-40% is flagged.",
        epilog="HOM, fringe and CHSH do not use this artifact yet (they take the dark dead time; M5).",
    )
    dt.add_argument("setup", help="setup.yaml (channel map) or a bundled example name")
    dt.add_argument("--tags", required=True, help="tag file, or a folder of PTU files, one rate per setting")
    dt.add_argument("--settings", type=int, nargs="+", help="setting indices holding the runs (default: every non-empty one)")
    dt.add_argument("--channel", action="append", help="channel name to calibrate (repeatable; default all)")
    dt.add_argument("--dark", help="folder or file with the channels' dark artifacts (afterpulsing); default: none")
    dt.add_argument("--out", required=True, help="directory for the artifacts")
    dt.add_argument(
        "--accept-misfit", action="store_true",
        help="write artifacts whose predictive check fails (p < 0.001), marked misfit; refused otherwise",
    )
    _add_record_options(dt)
    dt.set_defaults(func=cmd_calibrate_deadtime)

    pairs = cal_sub.add_parser(
        "pairs",
        help="delay, jitter, accidentals and heralding efficiency from a pair-source acquisition",
        description="Fit the delay and combined jitter of SIGNAL relative to HERALD from their cross-correlation, "
        "estimate off-peak accidentals in the declared window, and the heralding efficiency of SIGNAL. "
        "Uses the window width and delays declared in the setup. The acquisition must be a heralding run "
        "with no analysers or beamsplitter between source and detectors; for CHSH and HOM setups confirm "
        "this with --no-analysers. Writes three content-addressed artifacts.",
    )
    pairs.add_argument("setup", help="setup.yaml (channel map, window) or a bundled example name")
    pairs.add_argument("--tags", required=True, help="tag file, or a folder of PTU files, with the pair-source acquisition")
    pairs.add_argument("--setting", type=int, default=0, help="which setting index holds the data (default 0)")
    pairs.add_argument("--herald", required=True, help="herald channel name")
    pairs.add_argument("--signal", required=True, help="signal channel name")
    pairs.add_argument("--range-ps", type=int, default=20_000, help="delay search range, +- ps (default 20000)")
    pairs.add_argument("--bin-width-ps", type=int, default=20, help="delay histogram bin width (default 20)")
    pairs.add_argument("--out", required=True, help="directory for the artifacts")
    pairs.add_argument(
        "--accept-misfit", action="store_true",
        help="write artifacts whose predictive check fails (p < 0.001), marked misfit; refused otherwise",
    )
    pairs.add_argument(
        "--no-analysers", action="store_true",
        help="confirm the acquisition had no analysers or beamsplitter in the path (required for CHSH and HOM setups)",
    )
    _add_record_options(pairs)
    pairs.set_defaults(func=cmd_calibrate_pairs)

    ver = sub.add_parser("version", help="print the version")
    ver.set_defaults(func=lambda _args: print(__version__) or 0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except (DadlError, ValueError, OverflowError, OSError) as exc:
        # Library input errors (ValueError, OverflowError) reach CLI users as a message, not a
        # traceback (software review, cycle 2).
        print(f"dadl: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
