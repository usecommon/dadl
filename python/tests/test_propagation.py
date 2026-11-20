"""The propagation gate's machinery (phase 3, M2-M3). The replicate runs are in CI."""

from __future__ import annotations

import json

import pytest

from dadl.analysis.propagation import (
    CONTROL,
    CONTROL_REGIMES,
    KEYS,
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
    assert set(whole.covered) == set(KEYS) and whole.replicates == 2
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


def test_misspecified_regime_is_flagged():
    from dadl.analysis.propagation import FLAGGED

    assert REGIMES["hom-crosstalk"].misspecified and not REGIMES["hom-bright"].misspecified
    r = propagation(2, seed=7, regime="hom-crosstalk")
    assert r.covered[FLAGGED] == 2


def _result(regime, corrected=940, control=940, misfits=0, flagged=0):
    from dadl.analysis.propagation import CORRECTED, FLAGGED, MISFITS, result_name
    from dadl.calibration.coverage import CoverageResult, gate

    covered = {CORRECTED: corrected, CONTROL: control, MISFITS: misfits, FLAGGED: flagged}
    return CoverageResult(result_name(regime), 1000, covered, gate(1000, 1))


@pytest.mark.parametrize("regime, counts, verdict", [
    # Council round 7 (Q1): the verdicts on synthetic results, one rule at a time.
    ("hom-bright", {}, True),
    ("hom-bright", {"corrected": 900}, False),
    ("hom-bright", {"corrected": 990}, False),
    ("hom-bright", {"flagged": 20}, False),
    ("hom-bright", {"flagged": 10}, True),
    ("hom-bright", {"misfits": 30}, False),
    ("hom-bright", {"misfits": 20}, True),
    # The control regime must undercover with the calibration uncertainty off.
    ("chsh-dead-time", {"control": 940}, False),
    ("chsh-dead-time", {"control": 800}, True),
    # The misspecified regime is gated on its flags only.
    ("hom-crosstalk", {"flagged": 890, "corrected": 300}, False),
    ("hom-crosstalk", {"flagged": 950, "corrected": 300}, True),
])
def test_gate_verdicts(regime, counts, verdict):
    from dadl.analysis.propagation import passed

    assert passed(_result(regime, **counts)) is verdict


def test_the_control_drops_the_calibration_uncertainty():
    from dadl.analysis.propagation import covers
    from dadl.result import Estimate

    c = Estimate(1.0, 0.05, 0.902, 1.098, stderr_statistical=0.01, stderr_calibration=0.049)
    assert covers(c, 1.04) == (True, False)
    assert covers(c, 1.01) == (True, True)
    assert covers(c, 1.2) == (False, False)


def test_merge_refuses_empty_lax_and_foreign_shards():
    # Council round 7 (S1): an empty shard set exited 0, and a merge took its band from the file.
    from dadl.analysis.propagation import shard_record

    with pytest.raises(ValueError, match="family and regimes|analysis and regimes"):
        merge([{"seed": 20261111, "replicates": 1000, "range": [0, 1000], "results": []}])
    lax = {"seed": 1, "replicates": 1000, "range": [0, 1000], "results": [
        {"calibration": "hom_visibility [hom-snspd]", "covered": dict.fromkeys(KEYS, 10), "gate": [0, 1000]}]}
    with pytest.raises(ValueError, match="analysis and regimes"):
        merge([lax])
    good = shard_record([_result("hom-snspd")], 1, range(0, 1000), "hom", ["hom-snspd"])
    assert merge([good])[0].gate == _result("hom-snspd").gate
    with pytest.raises(ValueError, match="no results"):
        merge([{**good, "results": []}])
    with pytest.raises(ValueError, match="not the code's"):
        merge([{**good, "results": [{**good["results"][0], "gate": [0, 1000]}]}])
    with pytest.raises(ValueError, match="missing"):
        merge([{**good, "regimes": ["hom-bright", "hom-snspd"]}])
    with pytest.raises(ValueError, match="unexpected"):
        merge([{**good, "results": good["results"] + [
            {"calibration": "chsh [chsh-spad]", "covered": dict.fromkeys(KEYS, 0)}]}])


@pytest.mark.slow
@pytest.mark.parametrize("i", [19, 87, 131])
def test_chsh_dead_time_replicates_with_a_dark_miss_complete(i):
    # Council round 6 addendum: the lower-end gap test refused 42 of 200 chsh-dead-time replicates
    # (CI seed), each a dark dead time whose 95% lower end missed the 50 ns truth (by design in 5%
    # of channels). It labels them now; the replicate runs to completion.
    from dadl.analysis.propagation import REGIMES, _replicate
    from dadl.calibration.coverage import _children

    child = next(iter(_children(20261118, 200, range(i, i + 1), "chsh-dead-time")))
    assert len(_replicate(REGIMES["chsh-dead-time"], child)) == 4
