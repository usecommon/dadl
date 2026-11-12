"""Setup files, tag files and the CLI."""

from pathlib import Path

import numpy as np
import pytest

from dadl.cli import example_path, main
from dadl.config import load_setup
from dadl.errors import ConfigError
from dadl.tagfile import read_tags, write_tags

ROOT_EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "chsh-spad" / "setup.yaml"
EXAMPLE = example_path("chsh-spad")

GOOD = """\
schema_version: 1
name: x
seed: 5
source: {pair_rate_per_s: 100000, visibility: 0.9}
detectors:
  A+: &d {efficiency: 0.5, dark_rate_cps: 100, dead_time_ps: 0, jitter_sigma_ps: 50}
  A-: *d
  B+: *d
  B-: *d
analysis: {type: chsh, duration_s_per_setting: 0.05, window: {width_ps: 2000, source: plan}}
"""


def write(tmp_path, text):
    p = tmp_path / "setup.yaml"
    p.write_text(text, encoding="utf-8")
    return p


def test_bundled_example_loads_with_numeric_values():
    setup = load_setup(EXAMPLE)
    assert setup.source.pair_rate_per_s == 100000
    assert setup.window.source == "plan"
    assert setup.channels == {0: "A+", 1: "A-", 2: "B+", 3: "B-"}


@pytest.mark.skipif(not ROOT_EXAMPLE.exists(), reason="repository examples/ not present")
def test_repository_example_matches_bundled_copy():
    assert ROOT_EXAMPLE.read_bytes() == EXAMPLE.read_bytes()


@pytest.mark.parametrize(
    "replace, match",
    [
        (("source: plan", "source: data"), "plan"),
        (("schema_version: 1", "schema_version: 99"), "schema_version"),
        (("schema_version: 1", "schema_version: true"), "schema_version"),
        (("pair_rate_per_s: 100000", "pair_rate_per_s: 1.0e5"), "pair_rate_per_s"),
        (("analysis: {type", "analysis: x\nignored: {type"), "analysis"),
        (("width_ps: 2000", "width_ps: -5"), "width_ps"),
        (("efficiency: 0.5", "efficiency: 0.5, efficiency: 0.6"), "duplicate"),
        (("name: x", "name: x\nname: y"), "duplicate"),
        (("dead_time_ps: 0", "dead_time_ps: 0, colour: red"), "colour"),
        (("dark_rate_cps: 100", "dark_rate_cps: .nan"), "dark_rate_cps"),
    ],
)
def test_bad_setups_raise_config_error(tmp_path, replace, match):
    old, new = replace
    assert old in GOOD
    with pytest.raises(ConfigError, match=match):
        load_setup(write(tmp_path, GOOD.replace(old, new, 1)))


def test_directory_and_bad_encoding_raise_config_error(tmp_path):
    with pytest.raises(ConfigError):
        load_setup(tmp_path)
    p = tmp_path / "bad.yaml"
    p.write_bytes(b"\xff\xfe\x00name")
    with pytest.raises(ConfigError):
        load_setup(p)


def test_version_command_and_flag(capsys):
    assert main(["version"]) == 0
    assert capsys.readouterr().out.strip()
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0


def test_quickstart_runs(capsys):
    assert main(["quickstart"]) == 0
    out = capsys.readouterr().out
    assert "Raw S =" in out
    assert "fair sampling assumed" in out
    assert "usecommon.github.io/dadl" in out


def test_quickstart_reads_the_yaml(tmp_path, capsys):
    hi = write(tmp_path, GOOD)
    assert main(["quickstart", str(hi)]) == 0
    s_hi = capsys.readouterr().out
    lo = tmp_path / "lo.yaml"
    lo.write_text(GOOD.replace("visibility: 0.9", "visibility: 0.5"), encoding="utf-8")
    assert main(["quickstart", str(lo)]) == 0
    s_lo = capsys.readouterr().out
    value = lambda text: float(text.split("Raw S = ")[1].split()[0])  # noqa: E731
    assert value(s_hi) > value(s_lo) + 0.5


def test_errors_exit_1_with_message(tmp_path, capsys):
    assert main(["quickstart", str(tmp_path / "missing.yaml")]) == 1
    assert "dadl: error:" in capsys.readouterr().err
    with pytest.raises(SystemExit) as exc:
        main(["no-such-command"])
    assert exc.value.code == 2


def test_init_copies_example_and_refuses_overwrite(tmp_path, capsys):
    assert main(["init", str(tmp_path)]) == 0
    assert (tmp_path / "setup.yaml").read_bytes() == EXAMPLE.read_bytes()
    assert main(["init", str(tmp_path)]) == 1


@pytest.mark.parametrize("suffix", [".npz", ".csv"])
def test_simulate_then_analyze_round_trip(tmp_path, capsys, suffix):
    setup = write(tmp_path, GOOD)
    tags = tmp_path / f"tags{suffix}"
    assert main(["simulate", str(setup), "--out", str(tags)]) == 0
    capsys.readouterr()
    assert main(["quickstart", str(setup)]) == 0
    quick = capsys.readouterr().out.split("Raw S = ")[1].split()[0]
    assert main(["analyze", str(setup), "--tags", str(tags)]) == 0
    analyzed = capsys.readouterr().out.split("Raw S = ")[1].split()[0]
    assert quick == analyzed


def test_tag_file_round_trip_and_channel_map(tmp_path):
    runs = {k: {ch: np.array([1, 5, 9], dtype=np.int64) * (i + 1) for i, ch in enumerate(("A+", "A-", "B+", "B-"))}
            for k in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    cmap = {7: "A+", 3: "A-", 1: "B+", 0: "B-"}
    path = write_tags(tmp_path / "t.csv", runs, cmap)
    back = read_tags(path, cmap)
    for k in runs:
        for ch in runs[k]:
            assert back[k][ch].tolist() == runs[k][ch].tolist()


def test_tag_file_errors(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text("setting,channel,timestamp_ps\n0,9,100\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="channel"):
        read_tags(p)
    p.write_text("time,chan\n1,2\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="header"):
        read_tags(p)
    p.write_text("setting,channel,timestamp_ps\n0,0,1.5\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        read_tags(p)


def test_simulate_dark_then_calibrate(tmp_path, capsys):
    from dadl.artifact import read_artifact

    tags = tmp_path / "dark.npz"
    assert main(["simulate", "chsh-spad", "--dark", "--duration-s", "10", "--out", str(tags)]) == 0
    capsys.readouterr()
    runs = read_tags(tags)
    assert all(runs[("a", "b'")][ch].size == 0 for ch in ("A+", "B-"))
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(tags), "--channel", "B+", "--out", str(tmp_path / "art")]) == 0
    out = capsys.readouterr().out
    assert "dark_rate_cps" in out and "artifact:" in out
    (path,) = (tmp_path / "art").glob("*.json")
    manifest, samples = read_artifact(path)
    assert manifest["channels"] == [{"device": "B+", "input": 2}]
    assert manifest["backend"] == "sim"
    assert manifest["summary"]["dark_rate_cps"]["low"] < 500 < manifest["summary"]["dark_rate_cps"]["high"]


def test_calibrate_rejects_unknown_channel(tmp_path, capsys):
    tags = tmp_path / "dark.npz"
    main(["simulate", "chsh-spad", "--dark", "--duration-s", "1", "--out", str(tags)])
    assert main(["calibrate", "dark", "chsh-spad", "--tags", str(tags), "--channel", "C", "--out", str(tmp_path)]) == 1
    assert "unknown channel" in capsys.readouterr().err


def test_calibrate_pairs_writes_three_artifacts(tmp_path, capsys):
    from dadl.artifact import read_artifact

    tags = tmp_path / "tags.npz"
    assert main(["simulate", "chsh-spad", "--out", str(tags)]) == 0
    out_dir = tmp_path / "art"
    argv = ["calibrate", "pairs", "chsh-spad", "--tags", str(tags), "--herald", "A+", "--signal", "B+", "--out", str(out_dir)]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "delay_ps" in out and "accidental_coincidences" in out and "heralding_efficiency" in out
    kinds = sorted(read_artifact(p)[0]["kind"] for p in out_dir.glob("*.json"))
    assert kinds == ["accidentals", "delay-jitter", "heralding-efficiency"]
    assert main(argv[:7] + ["--signal", "A+", "--out", str(out_dir)]) == 1
    assert "different channels" in capsys.readouterr().err
