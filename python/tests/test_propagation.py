"""The propagation gate's machinery (phase 3, M2). The 200-replicate HOM tier runs in CI."""

from __future__ import annotations

import pytest

from dadl.analysis.propagation import CONTROL, CORRECTED, HOM_REGIMES, hom_propagation, main, run


@pytest.mark.parametrize("regime", sorted(HOM_REGIMES))
def test_replicates_run_and_parts_reproduce_the_whole(regime):
    whole = hom_propagation(2, seed=7, regime=regime)
    assert set(whole.covered) == {CORRECTED, CONTROL} and whole.replicates == 2
    # Replicate i depends only on (seed, regime, i), so a run split into parts adds up exactly.
    parts = [hom_propagation(2, seed=7, regime=regime, part=range(i, i + 1)) for i in (0, 1)]
    assert {k: sum(p.covered[k] for p in parts) for k in whole.covered} == whole.covered


def test_truth_is_the_two_point_visibility():
    for reg in HOM_REGIMES.values():
        assert reg.truth == pytest.approx(reg.source.visibility, abs=1e-6)


def test_unknown_regime_and_cli(capsys):
    with pytest.raises(ValueError, match="unknown"):
        run(2, regimes=["nope"])
    main(["hom", "--replicates", "2", "--regimes", "hom-snspd"])
    assert CONTROL in capsys.readouterr().out
