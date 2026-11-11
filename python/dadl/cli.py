"""Command-line interface.

    dadl quickstart [chsh-spad | path/to/setup.yaml]   simulate and analyse an example
    dadl init [DIR]                                     copy the example setup.yaml into DIR
    dadl simulate SETUP --out tags.npz                  write simulated tags in the tag-file format
    dadl analyze SETUP --tags tags.npz|tags.csv         raw CHSH S from your own tags
    dadl version
"""

from __future__ import annotations

import argparse
import shutil
import sys
from importlib import resources
from pathlib import Path

import numpy as np

from dadl import __version__
from dadl.analysis.chsh import chsh
from dadl.closed_form import chsh_expected_grid
from dadl.config import Setup, load_setup
from dadl.errors import DadlError
from dadl.result import UNCERTAINTY_KIND, LabelledValue
from dadl.sim import simulate_chsh
from dadl.tagfile import read_tags, write_tags

SITE = "https://usecommon.github.io/dadl/"
EXAMPLES = ("chsh-spad",)


def example_path(name: str) -> Path:
    """Path of a bundled example's setup.yaml."""
    if name not in EXAMPLES:
        raise DadlError(f"unknown example {name!r}; available: {', '.join(EXAMPLES)}")
    return Path(str(resources.files("dadl") / "examples" / name / "setup.yaml"))


def _resolve(target: str) -> Path:
    return example_path(target) if target in EXAMPLES else Path(target)


def _simulate(setup: Setup):
    if setup.backend != "sim":
        raise DadlError(f"{setup.name}: backend is {setup.backend!r}; only 'sim' setups can be simulated")
    rng = np.random.default_rng(setup.seed)
    return simulate_chsh(
        setup.source, setup.detectors, setup.duration_s_per_setting, rng, settings=setup.settings_deg
    )


def _print_result(result: LabelledValue) -> None:
    print(f"  Raw S = {result.raw}   ({UNCERTAINTY_KIND})")
    print(f"  Labels: {'; '.join(result.labels)}")
    w = result.window or {}
    print(f"  Window: {w.get('width_ps')} ps, delays {w.get('delays_ps') or 'none'}, source '{w.get('source')}' (declared, not verified)")
    print(f"  Policy: {result.policy_fingerprint[:12]}{' (modified)' if result.overrides else ' (default)'}")


def cmd_quickstart(args: argparse.Namespace) -> int:
    path = _resolve(args.target)
    setup = load_setup(path)
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
    print(f"  Calibration fits and corrected values arrive in phase 2. Overview: {SITE}")
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
    out = write_tags(args.out, _simulate(setup), setup.channels)
    print(f"Wrote simulated tags for {setup.name} to {out}")
    return 0


def cmd_analyze(args: argparse.Namespace) -> int:
    setup = load_setup(_resolve(args.setup))
    runs = read_tags(args.tags, setup.channels)
    result = chsh(runs, setup.window)
    print(f"DADL {__version__} analyze: {setup.name}, tags from {args.tags}")
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
        help="simulate an example (or your setup.yaml) and print raw CHSH S",
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
    sim.add_argument("--out", required=True, help="output file, .npz or .csv")
    sim.set_defaults(func=cmd_simulate)

    an = sub.add_parser("analyze", help="raw CHSH S from your own tag file")
    an.add_argument("setup", help="setup.yaml giving the window, delays and channel map")
    an.add_argument("--tags", required=True, help="tag file, .npz or .csv (columns setting, channel, timestamp_ps)")
    an.set_defaults(func=cmd_analyze)

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
