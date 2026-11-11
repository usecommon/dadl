"""Guardrail regression cases from published results.

These must keep passing. Changing them needs two maintainer approvals
(see CONTRIBUTING.md).
"""

import numpy as np
import pytest

from dadl.analysis.chsh import Window, chsh
from dadl.analysis.visibility import two_photon_visibility
from dadl.errors import GuardrailRefusal
from dadl.guardrails import (
    POLICY_MODIFIED,
    THRESHOLD_AFTER_CORRECTION,
    UNVALIDATED_ANALYSIS,
    Policy,
)

# Thompson, quant-ph/9903066: with ~30% accidentals, subtraction moved
# visibility from about 45% to 82%, past the 71% Bell threshold.
N_MAX, N_MIN = 145.0, 55.0  # raw V = 0.45
ACCIDENTALS = 45.12  # corrected V = 90 / (200 - 2 * 45.12) = 0.82


def test_raw_visibility_matches_published_case():
    result = two_photon_visibility(N_MAX, N_MIN, analysis="bell_visibility")
    assert result.raw.value == pytest.approx(0.45)
    assert result.corrected is None


def test_accidental_subtraction_refused_for_bell_visibility():
    with pytest.raises(GuardrailRefusal):
        two_photon_visibility(N_MAX, N_MIN, analysis="bell_visibility", accidentals=ACCIDENTALS)


def test_loosened_policy_is_never_silent():
    policy = Policy.default().allow(
        "bell_visibility", "accidental_subtraction", reason="reproducing Thompson's example"
    )
    result = two_photon_visibility(
        N_MAX, N_MIN, analysis="bell_visibility", accidentals=ACCIDENTALS, policy=policy
    )
    assert result.raw.value == pytest.approx(0.45)
    assert result.corrected.value == pytest.approx(0.82, abs=1e-3)
    assert POLICY_MODIFIED in result.labels
    assert THRESHOLD_AFTER_CORRECTION in result.labels


def test_override_needs_a_reason():
    with pytest.raises(ValueError):
        Policy.default().allow("bell_visibility", "accidental_subtraction", reason="  ")


def test_default_policy_is_not_modified():
    assert not Policy.default().is_modified
    assert Policy.default().fingerprint() == Policy().fingerprint()


def test_hom_allows_accidental_subtraction_as_labelled_secondary():
    result = two_photon_visibility(60.0, 40.0, analysis="hom_visibility", accidentals=15.0)
    assert result.raw.value == pytest.approx(0.2)
    assert result.corrected is not None
    assert POLICY_MODIFIED not in result.labels


def test_unknown_analysis_gets_strictest_policy():
    with pytest.raises(GuardrailRefusal):
        Policy.default().check_correction("my_plugin_analysis", "accidental_subtraction")
    assert UNVALIDATED_ANALYSIS in Policy.default().labels_for("my_plugin_analysis")


def test_window_chosen_from_data_is_refused():
    # Coincidence-time loophole (arXiv:1503.07573): windows chosen from the data
    # can fake a violation, so Bell analyses accept only planned windows.
    empty = {ch: np.array([], dtype=np.int64) for ch in ("A+", "A-", "B+", "B-")}
    runs = {key: empty for key in (("a", "b"), ("a", "b'"), ("a'", "b"), ("a'", "b'"))}
    with pytest.raises(GuardrailRefusal):
        chsh(runs, Window(width_ps=1000, source="data"))
