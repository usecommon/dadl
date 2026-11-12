"""Visibility analyses: error propagation, input checks, and the split entry points."""

import math

import numpy as np
import pytest

from dadl.analysis.visibility import bell_visibility, hom_visibility
from dadl.guardrails import GAUSSIAN_INVALID, NEAR_THRESHOLD, Policy

THOMPSON = (145.0, 55.0, 45.12)

LOOSE = Policy.default().allow("bell_visibility", "accidental_subtraction", reason="reproducing Thompson's example")


def test_corrected_fringe_error_with_exact_accidentals():
    m, n, a = THOMPSON
    r = bell_visibility(m, n, accidentals=a, accidentals_stderr=0.0, policy=LOOSE)
    assert r.corrected.value == pytest.approx(0.8196, abs=1e-3)
    assert r.corrected.stderr == pytest.approx(0.1245, abs=1e-3)


def test_corrected_fringe_error_with_poisson_accidentals():
    m, n, a = THOMPSON
    r = bell_visibility(m, n, accidentals=a, policy=LOOSE)
    assert r.corrected.stderr == pytest.approx(0.160, abs=2e-3)


def test_corrected_fringe_error_matches_monte_carlo():
    m, n, a = THOMPSON
    rng = np.random.default_rng(3)
    big, small, acc = rng.poisson(m, 100_000), rng.poisson(n, 100_000), rng.poisson(a, 100_000)
    d = big + small - 2.0 * acc
    ok = d > 0
    v = (big[ok] - small[ok]) / d[ok]
    # Trim the heavy tail of a ratio; compare the central spread (robust sigma).
    q16, q84 = np.percentile(v, [16, 84])
    r = bell_visibility(m, n, accidentals=a, policy=LOOSE)
    assert (q84 - q16) / 2 == pytest.approx(r.corrected.stderr, rel=0.10)


def test_hom_uses_dip_convention():
    r = hom_visibility(60.0, 40.0)
    assert r.raw.value == pytest.approx(1 / 3)
    assert r.raw.stderr == pytest.approx(math.sqrt(40**2 / 60**3 + 40 / 60**2))
    r = hom_visibility(60.0, 40.0, accidentals=15.0, accidentals_stderr=0.0)
    assert r.corrected.value == pytest.approx(20 / 45)


@pytest.mark.parametrize("fn", [bell_visibility, hom_visibility])
def test_corrected_visibility_cannot_exceed_one(fn):
    kwargs = {"policy": LOOSE} if fn is bell_visibility else {}
    with pytest.raises(ValueError):
        fn(100.0, 10.0, accidentals=20.0, **kwargs)


@pytest.mark.parametrize("args", [(-1.0, 5.0), (float("nan"), 5.0), (5.0, float("inf")), (0.0, 0.0), (5.0, 10.0)])
def test_bad_counts_refused(args):
    with pytest.raises((TypeError, ValueError)):
        bell_visibility(*args)


def test_negative_accidentals_refused():
    with pytest.raises(ValueError):
        hom_visibility(100.0, 10.0, accidentals=-50.0)


def test_zero_minimum_is_labelled_as_invalid_gaussian():
    r = bell_visibility(100.0, 0.0)
    assert r.raw.stderr == 0
    assert GAUSSIAN_INVALID in r.labels


def test_near_threshold_is_labelled():
    # V = 0.70: the exact 95% interval for the raw value contains 1/sqrt(2).
    from dadl.guardrails import NEAR_THRESHOLD_EXACT

    r = bell_visibility(85.0, 15.0)
    assert NEAR_THRESHOLD_EXACT in r.labels and r.raw.low < 0.7071 < r.raw.high
    # The corrected value (delta method) still uses the 1.96 sigma rule.
    assert NEAR_THRESHOLD not in r.labels


def test_low_counts_get_an_exact_interval_and_a_label():
    # Statistician review W8: near V = 1 the delta method covered 63-93%; the exact interval >= 95%.
    from dadl.analysis.visibility import hom_visibility
    from dadl.guardrails import LOW_COUNTS

    r = hom_visibility(50.0, 1.0)
    assert LOW_COUNTS in r.labels and r.raw.low < 0.98 <= r.raw.high <= 1.0
    assert LOW_COUNTS not in hom_visibility(5000.0, 400.0).labels


def test_results_record_policy_provenance():
    m, n, a = THOMPSON
    r = bell_visibility(m, n, accidentals=a, policy=LOOSE)
    assert r.policy_fingerprint == LOOSE.fingerprint()
    assert r.overrides and "reproducing Thompson" in r.overrides[0]
    rec = r.as_record()
    assert rec["policy_fingerprint"] == LOOSE.fingerprint()
    assert "reproducing Thompson" in rec["overrides"]
