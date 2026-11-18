"""The forward model's predictive check (phase 3, M3). Its rates over 1,000 replicates are in the
`propagation` gate; these are the mechanics."""

from __future__ import annotations

import numpy as np
import pytest

from dadl.errors import Misspecified
from dadl.forward_check import FORWARD_MISFIT_ACCEPTED, OVERRIDE, PairData, enforce, run_check
from dadl.guardrails import Policy
from dadl.sim import Crosstalk, DetectorParams, PairSource, simulate_hom

SRC = PairSource(pair_rate_per_s=2.0e6, visibility=0.92, coherence_time_ps=2.0)
DET = DetectorParams(efficiency=0.5, dark_rate_cps=20_000.0, dead_time_ps=60_000, jitter_sigma_ps=10.0)
JITTER = 10.0 * 2**0.5
XT = (Crosstalk("A", "B", 0.001, 400.0, 50.0), Crosstalk("B", "A", 0.001, 400.0, 50.0))


def _pairs(seed: int, crosstalk=(), declared_offsets=None):
    rng = np.random.default_rng(seed)
    out = []
    for name, d in (("dip", 0.0), ("reference", 20.0)):
        run = simulate_hom(SRC, DET, DET, d, 0.2, rng, crosstalk)
        offsets = declared_offsets if declared_offsets is not None and d else ((d, -d) if d else (0.0,))
        out.append(PairData(name, run["A"], run["B"], 0.0, 1_000, 0.0, JITTER, offsets))
    return out


def test_clean_runs_pass_and_crosstalk_is_flagged():
    clean = run_check(_pairs(1))
    assert set(clean.pvalues) == {"dip", "reference"} and not clean.flagged
    bad = run_check(_pairs(1, XT))
    assert bad.flagged and bad.pvalue < 1e-6


def test_a_misdeclared_reference_delay_is_flagged():
    r = run_check(_pairs(2, declared_offsets=(0.0,)))
    assert r.flagged and r.worst() == "reference"


def test_enforce_refuses_unless_allowed():
    bad = run_check(_pairs(3, XT))
    with pytest.raises(Misspecified, match=OVERRIDE):
        enforce(bad, "hom_visibility", Policy.default())
    loose = Policy.default().allow("hom_visibility", OVERRIDE, "test")
    assert enforce(bad, "hom_visibility", loose) == (FORWARD_MISFIT_ACCEPTED,)
    assert enforce(run_check(_pairs(3)), "hom_visibility", Policy.default()) == ()


def test_sparse_histograms_are_skipped():
    empty = np.array([], dtype=np.int64)
    r = run_check([PairData("none", empty, np.arange(10, dtype=np.int64), 0.0, 1_000)])
    assert r.pvalues == {} and r.pvalue == 1.0 and not r.flagged
