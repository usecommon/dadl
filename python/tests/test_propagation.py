"""The propagation gate's machinery (phase 3, M2-M3). The replicate runs are in CI."""

from __future__ import annotations

import json

import pytest

from dadl.analysis.propagation import (
    CONTROL,
    CONTROL_REGIMES,
    CORRECTED,
    MISFITS,
    REGIMES,
    main,
    merge,
    propagation,
    run,
    shard_range,
)


@pytest.mark.parametrize("regime", sorted(REGIMES))
def test_replicates_run_and_parts_reproduce_the_whole(regime):
    whole = propagation(2, seed=7, regime=regime)
    assert set(whole.covered) == {CORRECTED, CONTROL, MISFITS} and whole.replicates == 2
    # Replicate i depends only on (seed, regime, i), so a run split into parts adds up exactly.
    parts = [propagation(2, seed=7, regime=regime, part=range(i, i + 1)) for i in (0, 1)]
    assert {k: sum(p.covered[k] for p in parts) for k in whole.covered} == whole.covered


def test_truths():
    for reg in REGIMES.values():
        if reg.analysis == "hom":
            assert reg.truth == pytest.approx(reg.source.visibility, abs=1e-6)
        else:
            assert reg.truth == pytest.approx(reg.source.visibility * 2 * 2**0.5)
    assert set(CONTROL_REGIMES) <= set(REGIMES)


def test_unknown_regime_and_cli(capsys):
    with pytest.raises(ValueError, match="unknown"):
        run("hom", 2, regimes=["chsh-spad"])
    main(["hom", "--replicates", "2", "--regimes", "hom-snspd"])
    assert CONTROL in capsys.readouterr().out


def test_shards_merge_to_the_whole(tmp_path, capsys):
    whole = run("hom", 3, seed=7, regimes=["hom-snspd"])[0]
    files = []
    for k in (1, 2):
        out = tmp_path / f"part{k}.json"
        assert main(["hom", "--replicates", "3", "--seed", "7", "--regimes", "hom-snspd", "--shard", f"{k}/2", "--out", str(out)]) == 0
        files.append(str(out))
    assert "shard 2/2" in capsys.readouterr().out
    records = [json.loads(open(f).read()) for f in files]
    assert merge(records)[0].covered == whole.covered
    code = main(["merge", *files])
    assert "hom-snspd" in capsys.readouterr().out and code in (0, 1)
    with pytest.raises(ValueError, match="tile"):
        merge(records[:1] + records[:1])
    with pytest.raises(ValueError, match="not all"):
        merge(records[:1])
    with pytest.raises(ValueError, match="k/n"):
        shard_range("2", 10)
    assert shard_range("4/4", 1000) == range(750, 1000)
