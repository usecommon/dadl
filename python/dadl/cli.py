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
import shutil
import sys
from importlib import resources
from pathlib import Path


from dadl import __version__
from dadl.analysis.chsh import chsh
from dadl.closed_form import chsh_expected_grid, hom_visibility_expected
from dadl.config import Setup, load_setup
from dadl.errors import DadlError
from dadl.result import UNCERTAINTY_KIND, LabelledValue
from dadl.simsetup import simulate_setup, simulate_setup_dark
from dadl.tagfile import read_tags, tags_are_simulated, write_ptu_dir, write_tags

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
    w = result.window or {}
    print(f"  Window: {w.get('width_ps')} ps, delays {w.get('delays_ps') or 'none'}, source '{w.get('source')}' (declared, not verified)")
    print(f"  Policy: {result.policy_fingerprint[:12]}{' (modified)' if result.overrides else ' (default)'}")


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
    if result.corrected is not None:
        print(
            f"  Accidentals subtracted ({counts['accidentals']:.3g} +- {counts['accidentals_stderr']:.2g}, off-peak): "
            f"V = {result.corrected}"
        )
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
    result = chsh(runs, setup.window)
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
        )
        print(f"  Closed-form expectation (ignores dead time and afterpulsing): S = {expected:.4f}")
    print()
    print("  Edit a copy of this setup with `dadl init`, then `dadl quickstart setup.yaml`.")
    print("  Calibrate from dark data: `dadl simulate chsh-spad --dark --out dark.npz`, then")
    print("  `dadl calibrate dark chsh-spad --tags dark.npz --out artifacts/`; from pair data,")
    print("  `dadl calibrate pairs chsh-spad --tags tags.npz --herald A+ --signal B+ --out artifacts/`.")
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
    if args.dark:
        runs = simulate_setup_dark(setup, args.duration_s)
        what = f"dark tags ({args.duration_s:g} s per channel, setting 0)"
    else:
        runs = simulate_setup(setup)
        what = "HOM tags (setting 0 dip, 1 reference)" if setup.analysis_type == "hom" else "CHSH tags"
    if args.format == "ptu":
        files = write_ptu_dir(args.out, runs, setup.channels, setup.setting_keys, device=args.device)
        print(f"Wrote simulated {what} for {setup.name} as {len(files)} {args.device} PTU files in {args.out}")
        print("  The files are marked as simulated (CreatorSW_Name DADL, DADL_Simulated) and read back as such.")
        return 0
    out = write_tags(args.out, runs, setup.channels, setup.setting_keys)
    print(f"Wrote simulated {what} for {setup.name} to {out}")
    return 0


def cmd_calibrate_dark(args: argparse.Namespace) -> int:
    from dadl.artifact import write_artifact
    from dadl.calibration import calibrate_dark
    from dadl.data import ChannelId

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
    backend = "sim" if setup.backend == "sim" or tags_are_simulated(args.tags) else "hardware"
    print(f"DADL {__version__} calibrate dark: {setup.name}, tags from {args.tags} (setting {args.setting})")
    for name in wanted:
        cal = calibrate_dark(tags[name], ChannelId(name, number[name]), backend=backend, seed=setup.seed)
        path = write_artifact(cal, args.out)
        print()
        print(cal)
        print(f"  artifact: {path}")
    return 0


def cmd_calibrate_pairs(args: argparse.Namespace) -> int:
    from dadl.artifact import write_artifact
    from dadl.calibration import calibrate_accidentals, calibrate_delay, calibrate_heralding
    from dadl.data import ChannelId

    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels, setup.setting_keys)
    keys = setup.setting_keys
    if not 0 <= args.setting < len(keys):
        raise DadlError(f"--setting must be 0..{len(keys) - 1}")
    tags = runs[keys[args.setting]]
    number = {name: num for num, name in setup.channels.items()}
    for name in (args.herald, args.signal):
        if name not in number:
            raise DadlError(f"unknown channel {name!r}; the setup's channels are {sorted(number)}")
    if args.herald == args.signal:
        raise DadlError("--herald and --signal must be different channels")
    chs = (ChannelId(args.herald, number[args.herald]), ChannelId(args.signal, number[args.signal]))
    a, b = tags[args.herald], tags[args.signal]
    window = setup.window
    d_a, d_b = window.delays_ps.get(args.herald, 0), window.delays_ps.get(args.signal, 0)
    backend = "sim" if setup.backend == "sim" or tags_are_simulated(args.tags) else "hardware"
    common = {"backend": backend, "seed": setup.seed}
    print(f"DADL {__version__} calibrate pairs: {setup.name}, {args.herald} -> {args.signal}, tags from {args.tags}")
    print(f"  Window {window.width_ps} ps and delays {d_a}/{d_b} ps as declared in the setup.")
    cals = [
        calibrate_delay(a, b, chs, range_ps=args.range_ps, bin_width_ps=args.bin_width_ps, **common),
        calibrate_accidentals(a, b, chs, width_ps=window.width_ps, delay_a_ps=d_a, delay_b_ps=d_b, **common),
        calibrate_heralding(
            a, b, chs, width_ps=window.width_ps, delay_herald_ps=d_a, delay_signal_ps=d_b, **common
        ),
    ]
    for cal in cals:
        path = write_artifact(cal, args.out)
        print()
        print(cal)
        print(f"  artifact: {path}")
    fitted = cals[0].estimates["delay_ps"].value
    if abs(fitted - (d_b - d_a)) > window.width_ps / 4:
        print()
        print(
            f"  The fitted delay ({fitted:.0f} ps) differs from the declared one ({d_b - d_a} ps) by more than a "
            "quarter window. Declare it in the plan before analysis; DADL never moves the window to fit the data."
        )
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels, setup.setting_keys)
    print(f"DADL {__version__} analyze: {setup.name}, tags from {args.tags}")
    if setup.analysis_type == "hom":
        from dadl.analysis.hom import hom_from_tags

        result, counts = hom_from_tags(runs, setup.window)
        print("  The tags are trusted as given.")
        print()
        print(f"  Coincidences: {counts['c_min']:.0f} at the dip, {counts['c_max']:.0f} at the reference")
        _print_result(result, "V")
        if result.corrected is not None:
            print(f"  Accidentals subtracted ({counts['accidentals']:.3g}, off-peak): V = {result.corrected}")
        return 0
    result = chsh(runs, setup.window)
    print("  Phase 1 trusts the tags as given and reports raw S only.")
    print()
    _print_result(result)
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
        description="Simulate a CHSH experiment from a setup file and analyse it with a fixed, planned window.",
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
    dark.set_defaults(func=cmd_calibrate_dark)

    pairs = cal_sub.add_parser(
        "pairs",
        help="delay, jitter, accidentals and heralding efficiency from a pair-source acquisition",
        description="Fit the delay and combined jitter of SIGNAL relative to HERALD from their cross-correlation, "
        "estimate off-peak accidentals in the declared window, and the heralding efficiency of SIGNAL. "
        "Uses the window width and delays declared in the setup. For a heralding efficiency, record "
        "without polarisation analysers in the path. Writes three content-addressed artifacts.",
    )
    pairs.add_argument("setup", help="setup.yaml (channel map, window) or a bundled example name")
    pairs.add_argument("--tags", required=True, help="tag file, or a folder of PTU files, with the pair-source acquisition")
    pairs.add_argument("--setting", type=int, default=0, help="which setting index holds the data (default 0)")
    pairs.add_argument("--herald", required=True, help="herald channel name")
    pairs.add_argument("--signal", required=True, help="signal channel name")
    pairs.add_argument("--range-ps", type=int, default=20_000, help="delay search range, +- ps (default 20000)")
    pairs.add_argument("--bin-width-ps", type=int, default=20, help="delay histogram bin width (default 20)")
    pairs.add_argument("--out", required=True, help="directory for the artifacts")
    pairs.set_defaults(func=cmd_calibrate_pairs)

    ver = sub.add_parser("version", help="print the version")
    ver.set_defaults(func=lambda _args: print(__version__) or 0)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except DadlError as exc:
        print(f"dadl: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
