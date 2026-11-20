"""The propagation gate's machinery (phase 3, M2-M3). The replicate runs are in CI."""

from __future__ import annotations

import json
import math
from pathlib import Path

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


def _result(regime, corrected=940, control=940, misfits=0, flagged=0, drift=0, offpeak=0, refused=0):
    from dadl.analysis.propagation import CORRECTED, DRIFT, FLAGGED, MISFITS, OFFPEAK, REFUSED, result_name
    from dadl.calibration.coverage import CoverageResult, gate

    covered = {CORRECTED: corrected, CONTROL: control, MISFITS: misfits, FLAGGED: flagged, DRIFT: drift,
               OFFPEAK: offpeak}
    if refused:
        covered[REFUSED + "StaleCalibration"] = refused
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
    # Council round 8 (St2, St7): the flag cap covers the three checks together; any refusal fails.
    ("hom-bright", {"flagged": 4, "drift": 4, "offpeak": 4}, False),
    ("hom-bright", {"flagged": 4, "drift": 3, "offpeak": 3}, True),
    ("hom-bright", {"refused": 1}, False),
    ("hom-crosstalk", {"flagged": 950, "corrected": 300, "refused": 1}, False),
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


def test_a_refused_replicate_is_counted_and_fails_the_regime(monkeypatch, capsys):
    # Council round 8 (St2): a refusal aborted the whole gate with a traceback. It is counted by
    # type, as not covered, and fails the regime.
    import dadl.analysis.propagation as P
    from dadl.errors import StaleCalibration

    real, calls = P._analyse, []

    def planted(regime, child):
        calls.append(1)
        if len(calls) == 1:
            raise StaleCalibration("planted: the dark artifacts do not describe the analysed detectors")
        return real(regime, child)

    monkeypatch.setattr(P, "_analyse", planted)
    r = propagation(2, seed=7, regime="hom-snspd")
    assert r.covered[P.REFUSED + "StaleCalibration"] == 1 and r.covered[P.CORRECTED] <= 1
    assert P.refusals(r) == {"StaleCalibration": 1} and not P.passed(r)
    assert "refused: StaleCalibration" in P.describe(r)
    assert "planted" in capsys.readouterr().err


def test_merge_tiling_keeps_the_refusal_keys():
    from dadl.analysis.propagation import REFUSED, shard_record

    first = shard_record([_result("hom-snspd", 470, 470, refused=2)], 1, range(0, 500), "hom", ["hom-snspd"])
    second = shard_record([_result("hom-snspd", 470, 470)], 1, range(500, 1000), "hom", ["hom-snspd"])
    for rec in (first, second):
        for res in rec["results"]:
            res["covered"] = dict(res["covered"])
    (merged,) = merge([first, second])
    assert merged.covered[REFUSED + "StaleCalibration"] == 2 and merged.covered["corrected"] == 940
    with pytest.raises(ValueError, match="does not count"):
        bad = {**second, "results": [{**second["results"][0], "covered": {"corrected": 1}}]}
        merge([first, bad])


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
    assert not any(k.startswith("refused") for k in _replicate(REGIMES["chsh-dead-time"], child))


@pytest.mark.parametrize("i", [69, 360])
def test_r8_stat_gap_refusals_now_complete(i):
    # Council round 8 (St1, Fix 3): the 3-width floor refused these clean chsh-dead-time replicates
    # (seed 4242: "dead time 74.04 ns of B-", "dead time 128 ns"); the exponential-law floor does not.
    from dadl.analysis.propagation import REGIMES, _replicate
    from dadl.calibration.coverage import _children

    child = next(iter(_children(4242, 400, range(i, i + 1), "chsh-dead-time")))
    assert not any(k.startswith("refused") for k in _replicate(REGIMES["chsh-dead-time"], child))


def test_a_chsh_dead_time_replicate_carries_the_dead_time_above_gaps_label():
    # Council round 8 (Ph4, Fix 7): the dark dead times sit above the analysed shortest gap on
    # essentially every channel of this regime; the corrected value says so.
    from dadl.analysis.propagation import REGIMES, _analyse
    from dadl.calibration.coverage import _children
    # Council round 10 (Rg2): a bound has its own wording, which says the bound is loose.
    from dadl.forward import DEAD_TIME_ABOVE_GAPS, DEAD_TIME_ABOVE_GAPS_BOUND

    child = next(iter(_children(20261118, 200, range(0, 1), "chsh-dead-time")))
    result, _ = _analyse(REGIMES["chsh-dead-time"], child)
    heads = [t.split("{name}")[0] for t in (DEAD_TIME_ABOVE_GAPS, DEAD_TIME_ABOVE_GAPS_BOUND)]
    assert any(lab.startswith(h) for lab in result.corrected_labels for h in heads)


# --- Council round 9 ----------------------------------------------------------------------------

#: The flattened values of three gate replicates per regime, off-peak keys left out, pinned in
#: ``tests/data/r8_gate_replicates.json``: the round-8 code with only the delay method version raised to
#: 10 (artifact IDs seed the calibration draws) gives these, so the performance changes (as_sorted,
#: one-pass de-duplication, memoised grid counts, one quantile call, cached series constants) change no
#: value (council round 9, Pf2-Pf8). Round 9 pinned a SHA-256 of full-precision ``repr`` floats, which
#: differs on Linux libm/BLAS in the last bits (council round 10, Mg1); the fixture compares integers,
#: strings, booleans and None exactly and floats to a relative 1e-9.
#:
#: Generated on macOS, one BLAS thread, at e8c1aba after checking that the five round-9 hashes still
#: matched (``8d494923702d38aa``, ``d5d75abf34eaf9af``, ``5aa8e5cff09ddac9``, ``6dff3feb3e607303``,
#: ``20c150e54b341803``), with, from ``python/``:
#:
#:     OPENBLAS_NUM_THREADS=1 python -c "import sys; sys.path.insert(0, 'tests'); \
#:         import test_propagation as t; t.write_r8_fixture()"
R8_FIXTURE = Path(__file__).parent / "data" / "r8_gate_replicates.json"
R8_REGIMES = ("chsh-dead-time", "chsh-spad", "hom-bright", "hom-crosstalk", "hom-snspd")
R8_REL_TOL, R8_ABS_TOL = 1e-9, 1e-12


def _flatten(x, out, path=""):
    import numpy as np

    if isinstance(x, dict):
        for k in sorted(x):
            _flatten(x[k], out, f"{path}.{k}")
    elif isinstance(x, (list, tuple)):
        for i, v in enumerate(x):
            _flatten(v, out, f"{path}[{i}]")
    elif isinstance(x, np.generic):
        out.append([path, x.item()])
    elif hasattr(x, "__dict__") and not isinstance(x, (str, bytes)):
        _flatten(vars(x), out, path)
    else:
        out.append([path, x])


def _flatten_replicates(regime):
    from dadl.analysis import propagation as P

    out = []
    for i, child in enumerate(P._children(20261118, 1000, range(0, 3), regime)):
        r, mis = P._analyse(P.REGIMES[regime], child)
        counts = {k: v for k, v in (r.counts or {}).items() if k not in ("offpeak_p", "offpeak_draws")}
        _flatten((r.raw, r.corrected, counts, mis), out, f"r{i}")
    return out


def write_r8_fixture(path=R8_FIXTURE):
    path.parent.mkdir(exist_ok=True)
    blocks = []
    for r in R8_REGIMES:
        rows = ",\n".join(f"  {json.dumps(row)}" for row in _flatten_replicates(r))
        blocks.append(f"{json.dumps(r)}: [\n{rows}\n]")
    path.write_text("{\n" + ",\n".join(blocks) + "\n}\n")


def _same(a, b):
    """(equal, relative difference): floats to R8_REL_TOL, NaN equal to NaN, everything else exactly."""
    if type(a) is float and type(b) is float:
        if math.isnan(a) or math.isnan(b):
            return (True, 0.0) if math.isnan(a) and math.isnan(b) else (False, math.inf)
        rel = abs(a - b) / max(abs(a), abs(b)) if a != b else 0.0
        return math.isclose(a, b, rel_tol=R8_REL_TOL, abs_tol=R8_ABS_TOL), rel
    equal = type(a) is type(b) and a == b
    return equal, (0.0 if equal else math.inf)


def compare_to_fixture(got, pinned):
    """Problems between two flattened value lists: path differences, then the five worst values."""
    gp, pp = [p for p, _ in got], [p for p, _ in pinned]
    if gp != pp:
        missing, extra = sorted(set(pp) - set(gp)), sorted(set(gp) - set(pp))
        return [f"paths differ: missing {missing[:5]}, extra {extra[:5]}" + ("" if missing or extra else ", order")]
    bad = []
    for (path, a), (_, b) in zip(got, pinned):
        ok, rel = _same(a, b)
        if not ok:
            bad.append((rel, path, a, b))
    bad.sort(key=lambda t: -t[0])
    return [f"{path}: got {a!r}, pinned {b!r} (relative difference {rel:.3g})" for rel, path, a, b in bad[:5]]


@pytest.mark.slow
@pytest.mark.parametrize("regime", R8_REGIMES)
def test_the_performance_changes_leave_the_gate_replicates_unchanged(regime):
    pinned = json.loads(R8_FIXTURE.read_text())[regime]
    problems = compare_to_fixture(_flatten_replicates(regime), pinned)
    assert not problems, "\n".join(problems)


def test_the_replicate_fixture_comparison_catches_a_changed_count_and_a_changed_value():
    # Council round 10 (Mg1, Fix 1): the test of the test. A count changed by 1 and a float changed by
    # 1e-6 relative both fail; last-bit noise and NaN against NaN pass.
    pinned = json.loads(R8_FIXTURE.read_text())["hom-snspd"]
    assert compare_to_fixture(pinned, pinned) == []
    i_int = next(i for i, (_, v) in enumerate(pinned) if type(v) is int)
    i_float = next(i for i, (_, v) in enumerate(pinned) if isinstance(v, float) and v != 0 and math.isfinite(v))
    count = [list(x) for x in pinned]
    count[i_int][1] += 1
    assert len(compare_to_fixture(count, pinned)) == 1
    value = [list(x) for x in pinned]
    value[i_float][1] *= 1 + 1e-6
    assert len(compare_to_fixture(value, pinned)) == 1
    noise = [list(x) for x in pinned]
    noise[i_float][1] *= 1 + 1e-14
    assert compare_to_fixture(noise, pinned) == []
    assert compare_to_fixture([["a", math.nan]], [["a", math.nan]]) == []
    assert compare_to_fixture([["a", 1]], [["a", 1.0]]) != []
    assert compare_to_fixture([["b", 1]], [["a", 1]]) != []


# --- Council round 10 ---------------------------------------------------------------------------


def test_describe_says_the_control_rule_needs_about_400_replicates():
    # Council round 10 (Rg4, Fix 6): at 100 replicates chsh-dead-time's control (true 87.5%) lands
    # inside the band about one run in three; `passed` is unchanged (a gate rule, O12).
    from dadl.analysis import propagation as P

    def result(n, regime):
        covered = {k: 0 for k in KEYS}
        covered[P.CORRECTED] = int(0.95 * n)
        covered[CONTROL] = int(0.93 * n)
        return P.CoverageResult(P.result_name(regime), n, covered, P.gate(n, 1))

    note = "control rule not decidable below about 400 replicates"
    assert note in P.describe(result(100, "chsh-dead-time")) and not P.passed(result(100, "chsh-dead-time"))
    assert note not in P.describe(result(400, "chsh-dead-time"))
    assert note not in P.describe(result(100, "chsh-spad"))
