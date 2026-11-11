"""Example setup file and CLI smoke tests."""

from pathlib import Path

import pytest

from dadl.cli import main
from dadl.config import load_setup
from dadl.errors import ConfigError

EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "chsh-spad" / "setup.yaml"


def test_example_setup_loads():
    setup = load_setup(EXAMPLE)
    assert setup["schema_version"] == 1
    assert setup["analysis"]["window"]["source"] == "plan"


def test_window_from_data_rejected_in_config(tmp_path):
    bad = tmp_path / "setup.yaml"
    bad.write_text(
        "schema_version: 1\nname: x\nsource: {}\ndetectors: {}\n"
        "analysis: {type: chsh, window: {width_ps: 1000, source: data}}\n",
        encoding="utf-8",
    )
    with pytest.raises(ConfigError):
        load_setup(bad)


def test_unsupported_schema_version(tmp_path):
    bad = tmp_path / "setup.yaml"
    bad.write_text("schema_version: 99\nname: x\nsource: {}\ndetectors: {}\nanalysis: {}\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_setup(bad)


def test_version_command(capsys):
    assert main(["version"]) == 0
    assert capsys.readouterr().out.strip()


def test_quickstart_runs(capsys):
    assert main(["quickstart", "chsh-spad"]) == 0
    out = capsys.readouterr().out
    assert "Raw S =" in out
    assert "fair sampling assumed" in out
