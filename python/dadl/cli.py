"""Command-line interface.

    dadl quickstart [chsh-spad | path/to/setup.yaml]   simulate and analyse an example
    dadl init [DIR]                                     copy the example setup.yaml into DIR
    dadl simulate SETUP --out tags.npz [--dark]         write simulated tags in the tag-file format
    dadl analyze SETUP --tags tags.npz|tags.csv         raw CHSH S or HOM visibility from your own tags
    dadl calibrate dark SETUP --tags dark.npz --out DIR fit dark rate, afterpulsing and dead time
    dadl calibrate pairs SETUP --tags F --herald A+ --signal B+ --out DIR
                                                        delay, jitter, accidentals, heralding efficiency
    dadl version
"""

from __future__ import annotations

import argparse
import dataclasses
import shutil
import sys
from importlib import resources
from pathlib import Path

import numpy as np

from dadl import __version__
from dadl.analysis.chsh import chsh
from dadl.guardrails import POLICY_MODIFIED, SIMULATED_INPUT
from dadl.closed_form import chsh_expected_grid, chsh_ideal, hom_visibility_expected
from dadl.config import Setup, load_setup
from dadl.errors import DadlError
from dadl.result import UNCERTAINTY_KIND, LabelledValue
from dadl.simsetup import simulate_setup, simulate_setup_dark, simulate_setup_heralding
from dadl.tagfile import read_tags, tag_tick_ps, tagger_settings, tags_are_simulated, write_ptu_dir, write_tags

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
    print(f"  Labels: {'; '.join(result.labels) or 'none'}")
    if result.analysis == "chsh" and result.counts:
        ratios = ", ".join(
            f"{k} {v['coincidences_per_singles']:.3f}"
            for k, v in result.counts.items()
            if isinstance(v, dict) and "coincidences_per_singles" in v
        )
        print(f"  Coincidences per singles (what post-selection kept): {ratios}")
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
    if result.corrected is not None:
        print(
            f"  Accidentals subtracted ({counts['accidentals']:.3g} +- {counts['accidentals_stderr']:.2g}, off-peak): "
            f"V = {result.corrected}  (secondary; not an indistinguishability claim)"
        )
    elif counts["accidentals"] > 0:
        print(f"  Accidentals not subtracted (estimate {counts['accidentals']:.3g}); see the labels.")
    ref = abs(setup.reference_delay_ps or 0)
    if ref > 0.1 * setup.window.width_ps:
        print(
            f"  warning: the reference delay ({ref:g} ps) is not small against the {setup.window.width_ps} ps window: "
            "the fixed grid loses part of the split peak at +-delay, and V reads low."
        )


def _quickstart_hom(path: Path, setup: Setup) -> int:
    from dadl.analysis.hom import hom_from_tags

    runs = _simulate(setup)
    result, counts = hom_from_tags(runs, setup.window)
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
    print(f"  Corrected values arrive in phase 3. Overview: {SITE}")
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
    if not args.duration_s > 0:
        raise DadlError("--duration-s must be positive")
    if args.dark and args.heralding:
        raise DadlError("choose one of --dark and --heralding")
    if args.dark:
        runs = simulate_setup_dark(setup, args.duration_s)
        what = f"dark tags ({args.duration_s:g} s per channel, setting 0)"
    elif args.heralding:
        runs = simulate_setup_heralding(setup, args.duration_s)
        what = f"heralding-run tags ({args.duration_s:g} s, no analysers or beamsplitter, setting 0)"
    else:
        runs = simulate_setup(setup)
        what = "HOM tags (setting 0 dip, 1 reference)" if setup.analysis_type == "hom" else "CHSH tags"
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


def _add_record_options(p: argparse.ArgumentParser) -> None:
    """Options that complete the artifact as a record (calibration-lab review, cycle 3)."""
    p.add_argument("--condition", action="append", default=[], metavar="KEY=VALUE",
                   help="a numeric condition to record, e.g. temperature_K=233.1 or bias_V=24.5 (repeatable)")
    p.add_argument("--acquired", help="when the data were taken (ISO 8601); default: from the PTU header, if any")
    p.add_argument("--created", help="record creation time (ISO 8601), to reproduce an artifact; default now")
    p.add_argument("--valid-until", help="end of the intended validity (ISO 8601); recorded, not enforced")
    p.add_argument("--supersedes", help="ID of the artifact this one replaces (must exist in --out)")


def _when(text: str | None, name: str):
    import datetime as dt

    if text is None:
        return None
    try:
        t = dt.datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as exc:
        raise DadlError(f"--{name}: not an ISO 8601 time: {text!r}") from exc
    return t if t.tzinfo else t.replace(tzinfo=dt.timezone.utc)


def _record(args: argparse.Namespace, setup_path: Path) -> tuple[dict, dict[str, float], tuple]:
    """Write options, extra conditions and extra inputs for this run's artifacts."""
    import datetime as dt
    import hashlib

    conditions: dict[str, float] = {}
    for item in args.condition:
        key, sep, value = item.partition("=")
        try:
            conditions[key.strip()] = float(value)
        except ValueError:
            sep = ""
        if not sep or not key.strip():
            raise DadlError(f"--condition needs KEY=NUMBER, got {item!r}")
    conditions.update(tagger_settings(args.tags))
    tags_path = Path(args.tags)
    files = sorted(tags_path.glob("setting*.ptu")) if tags_path.is_dir() else [tags_path]
    inputs = tuple((hashlib.sha256(f.read_bytes()).hexdigest(), "raw-file") for f in files if f.is_file())
    if setup_path.is_file():
        inputs += ((hashlib.sha256(setup_path.read_bytes()).hexdigest(), "config"),)
    acquired, instrument = _when(args.acquired, "acquired"), {}
    if tags_path.is_dir() and files:
        from dadl.adapters.ptu import PtuFile

        header = PtuFile(files[0]).header
        instrument = {"model": header.get("HW_Type"), "serial": header.get("HW_SerialNo"),
                      "firmware": header.get("HW_Version")}
        stamp = header.get("File_CreatingTime")
        if acquired is None and isinstance(stamp, float) and stamp > 0:
            acquired = dt.datetime(1899, 12, 30, tzinfo=dt.timezone.utc) + dt.timedelta(days=stamp)
    supersedes = args.supersedes
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
            "in setup.yaml; for simulated data write .npz or PTU instead of .csv."
        )
    path = Path(tags_path)
    if path.is_dir() and setup.tagger:
        from dadl.adapters.ptu import PtuFile

        serials = {str(PtuFile(f).header.get("HW_SerialNo", "")).strip() for f in sorted(path.glob("setting*.ptu"))}
        serials.discard("")
        wrong = sorted(sn for sn in serials if sn not in setup.tagger)
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
            if short > 0.05:
                raise DadlError(
                    f"{name} in acquisition {key}: {short:.0%} of the gaps between tags are shorter than half the "
                    f"declared dead time ({det.dead_time_ps} ps). Are the timestamps in picoseconds (not ns)? "
                    "If they are, correct dead_time_ps in setup.yaml."
                )


def _refuse_source_on(tags: dict, setup) -> None:
    """Refuse a "dark" acquisition whose channels share coincidences: the source was on.

    The new-user review ran `calibrate dark` on source-on CHSH tags and got a confident dark rate
    of 25,400 cps (true 500), written as a normal artifact.
    """
    from dadl.coincidences import grid_coincidences

    names = [n for n in sorted(tags) if np.asarray(tags[n]).size > 1]
    width = setup.window.width_ps
    delays = setup.window.delays_ps
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            a, b = np.asarray(tags[x]), np.asarray(tags[y])
            span = float(max(a[-1], b[-1]) - min(a[0], b[0]))
            if span <= 0:
                continue
            c = grid_coincidences(a, b, width, delays.get(x, 0), delays.get(y, 0))
            expected = a.size * b.size * width / span
            if c - expected > 10 and (c - expected) / np.sqrt(expected + 1.0) > 6:
                raise DadlError(
                    f"channels {x} and {y} share {c} coincidences where {expected:.3g} are expected by chance: "
                    "the source was on. A dark calibration needs an acquisition with the source off "
                    "(simulated: dadl simulate SETUP --dark)."
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
    for name in wanted:
        arr = tags[name]
        run = {"setting": float(args.setting)}
        if arr.size > 1:
            run["duration_s"] = (int(arr[-1]) - int(arr[0])) / 1e12
        cal = calibrate_dark(arr, setup.channel_id(name), backend=backend, seed=setup.seed, conditions={**run, **extra})
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
        if declared > 0 and fitted > 10 * declared:
            print(
                f"  warning: the dark rate is more than 10x the setup's declared {declared:g} cps. Check that the "
                "source was off and the detector shielded from light."
            )
        path = _write_checked(cal, args.out, args.accept_misfit, **write)
        print(f"  artifact: {path}")
    return 0


def _grid_capture(offset_ps: float, sigma_ps: float, width_ps: int) -> float:
    """Probability that a pair lands in one fixed-grid slot: E[max(0, 1 - |D| / W)], D ~ N(offset, sigma)."""
    from scipy import integrate, stats

    if sigma_ps <= 0:
        return max(0.0, 1.0 - abs(offset_ps) / width_ps)
    pdf = stats.norm(offset_ps, sigma_ps).pdf
    value, _ = integrate.quad(lambda d: (1.0 - abs(d) / width_ps) * pdf(d), -width_ps, width_ps, points=[0.0], limit=200)
    return float(value)


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
        if misplaced and cal is cals[2] and not args.accept_misfit:
            refused.append(cal.kind)
            print("  not written: the declared delay does not put pairs in the declared slot; declare the fitted delay")
            continue
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
    capture = _grid_capture(fitted - (d_b - d_a), cals[0].estimates["jitter_sigma_ps"].value, window.width_ps)
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
        print(
            f"  Source figures (secondary, not recorded): coincidences-to-accidentals {car}; pair rate from "
            f"singles x singles / coincidences {a.size * b.size / n_c / span_s:.4g} per s (diluted and dead-time "
            "limited like the heralding figure)."
        )
    if refused:
        raise DadlError(f"{len(refused)} of {len(cals)} artifacts not written (misfit): {', '.join(refused)}")
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels, setup.setting_keys)
    _check_time_units(runs, setup)
    print(f"DADL {__version__} analyze: {setup.name}, tags from {args.tags}")
    # Red-team review, cycle 5: simulated tags gave an unmarked headline, and nothing tied the result
    # to the setup and tags it read.
    extra = {"inputs": _input_hashes(_resolve(args.setup), args.tags)}
    simulated = tags_are_simulated(args.tags)
    if setup.analysis_type == "hom":
        from dadl.analysis.hom import hom_from_tags

        result, counts = hom_from_tags(runs, setup.window)
        result = dataclasses.replace(
            result, labels=result.labels + ((SIMULATED_INPUT,) if simulated else ()), **extra
        )
        print("  The tags are trusted as given.")
        print()
        print(f"  Coincidences: {counts['c_min']:.0f} at the dip, {counts['c_max']:.0f} at the reference")
        _print_result(result, "V")
        _print_hom_extras(result, counts, setup)
        return 0
    result = chsh(runs, setup.window, settings_deg=setup.settings_deg)
    result = dataclasses.replace(result, labels=result.labels + ((SIMULATED_INPUT,) if simulated else ()), **extra)
    _warn_angles(setup)
    print("  The tags are trusted as given; raw S only.")
    print()
    _print_result(result)
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
        help="simulate an example (or your setup.yaml) and print its raw result (CHSH S or HOM visibility)",
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

    sim.add_argument("--device", choices=sorted(DEVICES), default="hydraharp", help="PTU device to mimic (default hydraharp)")
    sim.add_argument("--dark", action="store_true", help="simulate dark acquisitions (source off) for calibration")
    sim.add_argument(
        "--heralding", action="store_true",
        help="simulate a heralding run: the pair source straight onto one detector per arm (for calibrate pairs)",
    )
    sim.add_argument("--duration-s", type=float, default=20.0, help="dark acquisition length per channel (default 20)")
    sim.set_defaults(func=cmd_simulate)

    an = sub.add_parser("analyze", help="raw CHSH S or HOM visibility from your own tag file")
    an.add_argument("setup", help="setup.yaml giving the window, delays and channel map")
    an.add_argument("--tags", required=True, help="tag file (.npz or .csv: columns setting, channel, timestamp_ps) or a folder of PTU files setting0.ptu, ...")
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
    except (DadlError, ValueError, OverflowError) as exc:
        # Library input errors (ValueError, OverflowError) reach CLI users as a message, not a
        # traceback (software review, cycle 2).
        print(f"dadl: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
