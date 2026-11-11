"""Guardrail regression cases from published results, and bypass attempts.

These must keep passing. Changing them needs maintainer review under the rule
in CONTRIBUTING.md.
"""

import dataclasses

import numpy as np
import pytest

from dadl.analysis.chsh import Window, chsh
from dadl.analysis.visibility import bell_visibility, hom_visibility
from dadl.errors import GuardrailRefusal
from dadl.guardrails import (
    DEFAULT_POLICY_FINGERPRINT,
    DEFAULT_RULES,
    NEAR_THRESHOLD,
    POLICY_MODIFIED,
    STRICTEST,
    THRESHOLD_AFTER_CORRECTION,
    THRESHOLDS,
    UNVALIDATED_ANALYSIS,
    AnalysisRule,
    Policy,
    threshold_labels,
)
from dadl.result import Estimate

# Thompson, quant-ph/9903066: with ~30% accidentals, subtraction moved
# visibility from about 45% to 82%, past the 71% Bell threshold.
N_MAX, N_MIN = 145.0, 55.0  # raw V = 0.45
ACCIDENTALS = 45.12  # corrected V = 90 / (200 - 2 * 45.12) = 0.82


def test_raw_visibility_matches_published_case():
    result = bell_visibility(N_MAX, N_MIN)
    assert result.raw.value == pytest.approx(0.45)
    assert result.corrected is None


def test_accidental_subtraction_refused_for_bell_visibility():
    with pytest.raises(GuardrailRefusal):
        bell_visibility(N_MAX, N_MIN, accidentals=ACCIDENTALS)


def test_bell_data_cannot_be_relabelled_as_hom():
    # 0.0.2 accepted analysis="hom_visibility" as an argument; the analysis type
    # is now fixed by the function, so the old bypass is a TypeError.
    with pytest.raises(TypeError):
        bell_visibility(1450, 550, analysis="hom_visibility", accidentals=450)  # type: ignore[call-arg]


def test_loosened_policy_is_never_silent():
    policy = Policy.default().allow(
        "bell_visibility", "accidental_subtraction", reason="reproducing Thompson's example"
    )
    result = bell_visibility(N_MAX, N_MIN, accidentals=ACCIDENTALS, policy=policy)
    assert result.raw.value == pytest.approx(0.45)
    assert result.corrected.value == pytest.approx(0.82, abs=1e-3)
    assert POLICY_MODIFIED in result.labels
    assert THRESHOLD_AFTER_CORRECTION in result.labels
    assert result.overrides == ("bell_visibility: accidental_subtraction (reproducing Thompson's example)",)


def test_override_needs_a_reason():
    with pytest.raises(ValueError):
        Policy.default().allow("bell_visibility", "accidental_subtraction", reason="  ")


def test_override_needs_a_known_correction():
    with pytest.raises(ValueError):
        Policy.default().allow("chsh", "window:data", reason="because")


def test_default_policy_fingerprint_is_pinned():
    assert Policy.default().fingerprint() == DEFAULT_POLICY_FINGERPRINT
    assert not Policy.default().is_modified


def test_rule_tables_are_read_only():
    with pytest.raises(TypeError):
        THRESHOLDS["hom_visibility"] = 0.0  # type: ignore[index]
    with pytest.raises(TypeError):
        DEFAULT_RULES["chsh"] = AnalysisRule(frozenset({"efficiency"}), frozenset({"data"}))  # type: ignore[index]
    with pytest.raises(TypeError):
        Policy.default().rules["chsh"] = STRICTEST  # type: ignore[index]


def test_looser_rules_must_go_through_allow():
    loose = AnalysisRule(allowed_corrections=frozenset(), window_sources=frozenset({"plan", "data"}))
    with pytest.raises(ValueError, match="allow"):
        Policy(rules={**DEFAULT_RULES, "chsh": loose})


def test_stricter_custom_rules_are_accepted_but_labelled():
    strict = Policy(rules={**DEFAULT_RULES, "hom_visibility": STRICTEST})
    assert strict.is_modified
    with pytest.raises(GuardrailRefusal):
        hom_visibility(60.0, 40.0, accidentals=5.0, policy=strict)


def test_policies_are_hashable_and_compare_by_content():
    assert Policy() == Policy.default()
    assert hash(Policy()) == hash(Policy.default())


def test_hom_allows_accidental_subtraction_as_labelled_secondary():
    result = hom_visibility(60.0, 40.0, accidentals=15.0)
    assert result.raw.value == pytest.approx(1 / 3)
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


def test_chsh_threshold_labels():
    assert NEAR_THRESHOLD in threshold_labels("chsh", Estimate(2.05, 0.04))
    assert NEAR_THRESHOLD not in threshold_labels("chsh", Estimate(2.6, 0.01))


def test_results_cannot_be_used_as_bare_numbers():
    r = bell_visibility(N_MAX, N_MIN)
    with pytest.raises(TypeError):
        float(r)
    # Not a forgery guard: explicit extraction and in-process edits remain possible,
    # which the docs state. This pins that the label survives an honest replace().
    assert dataclasses.replace(r).labels == r.labels
