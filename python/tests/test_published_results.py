"""Published results as executable tests (docs/design.md, "Review findings": the open high-risk item).

Each section rebuilds one published result from simulated time tags and checks what DADL does
with it. The figures quoted in the docstrings were checked against the sources on 2026-11-15;
docs/design.md cites the same sources in its References.

1. Accidental subtraction: Thompson, quant-ph/9903066, on Tittel et al. (1997).
2. The coincidence-time loophole: Christensen et al., arXiv:1503.07573 (NIST and the University
   of Illinois), and Larsson et al., arXiv:1309.0712.
3. The NIST calibration comparison: Polyakov and Migdall, Proc. SPIE 6583, 65830A (2007),
   "High accuracy calibration of photon-counting detectors on demand".

These are simulations of the published mechanisms, not the published data. Changing them needs
maintainer review under the guardrail rule in CONTRIBUTING.md.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from dadl.analysis.chsh import FAIR_SAMPLING, SIGNS, Window, chsh
from dadl.analysis.visibility import bell_visibility
from dadl.calibration import calibrate_accidentals, calibrate_heralding
from dadl.closed_form import heralding_ratio_expected, nonparalyzable_measured_rate
from dadl.coincidences import grid_coincidences
from dadl.data import ChannelId
from dadl.errors import GuardrailRefusal
from dadl.guardrails import POLICY_MODIFIED, THRESHOLD_AFTER_CORRECTION, Policy
from dadl.sim import DetectorParams, PairSource, detect, poisson_times, simulate_chsh_setting, simulate_pairs

BELL_VISIBILITY_THRESHOLD = 1.0 / math.sqrt(2.0)

# --- 1. Accidental subtraction (Thompson, quant-ph/9903066) ---------------------------------------
#
# Thompson reads Tittel et al.'s 1997 Geneva data: accidentals about 30% of the rate, and "the
# adjustment changes the visibility from about 45% to 82% ... sufficient to bring the value above
# the Bell limit of 71%". Here a source with pair visibility 0.82 sends its photons through a noisy
# link (dark/background counts on every detector) so that accidentals are about 31% of the fringe
# maximum and 45% of the mean count, as in the published case. Raw V = c V_pair / (c + A) with c the
# mean true count and A the accidentals per run: c / (c + A) = 0.45 / 0.82.

FRINGE_SOURCE = PairSource(pair_rate_per_s=1e5, visibility=0.82)
FRINGE_WIDTH_PS = 5_000
FRINGE_DURATION_S = 4.0
#: Singles r per channel with r^2 W = 0.822 x (R eta^2 / 4): A / c = 0.45 / 0.82 inverted.
FRINGE_DET = DetectorParams(efficiency=0.1, dark_rate_cps=1.94e5, jitter_sigma_ps=50.0, dead_time_ps=22_000)


@pytest.fixture(scope="module")
def fringe():
    """N_max (analysers parallel) and N_min (crossed) on A+ B+, with off-peak accidentals per run."""
    rng = np.random.default_rng(1999)
    dets = {ch: FRINGE_DET for ch in ("A+", "A-", "B+", "B-")}
    out = {}
    for name, beta in (("max", 0.0), ("min", 90.0)):
        tags = simulate_chsh_setting(FRINGE_SOURCE, 0.0, beta, dets, FRINGE_DURATION_S, rng)
        a, b = tags["A+"], tags["B+"]
        n = grid_coincidences(a, b, FRINGE_WIDTH_PS)
        cal = calibrate_accidentals(
            a, b, (ChannelId("sim", 0), ChannelId("sim", 2)), width_ps=FRINGE_WIDTH_PS, backend="sim", seed=1
        )
        out[name] = (n, cal.estimates["accidental_coincidences"].value)
    return out


def test_accidental_subtraction_from_simulated_tags_is_refused(fringe):
    """Thompson (quant-ph/9903066), Tittel 1997 data: 45% raw, 82% after subtraction, threshold 71%.

    From simulated time tags (not hand-set counts): the raw fringe visibility stays below 1/sqrt(2)
    with its exact interval, subtraction of DADL's own off-peak accidentals estimate would cross
    it, and the default policy refuses the subtracted value.
    """
    (n_max, a_max), (n_min, a_min) = fringe["max"], fringe["min"]
    accidentals = 0.5 * (a_max + a_min)
    # The published regime: accidentals about 30% of the fringe maximum.
    assert 0.25 < accidentals / n_max < 0.36

    raw = bell_visibility(n_max, n_min, window=Window(FRINGE_WIDTH_PS, "plan"))
    assert raw.raw.value == pytest.approx(0.45, abs=3 * raw.raw.stderr)
    assert raw.raw.high < BELL_VISIBILITY_THRESHOLD  # raw, with its exact 95% interval, does not cross

    with pytest.raises(GuardrailRefusal):
        bell_visibility(n_max, n_min, accidentals=accidentals, window=Window(FRINGE_WIDTH_PS, "plan"))


def test_subtracted_visibility_reproduces_the_published_82_percent_only_under_a_labelled_override(fringe):
    """Thompson (quant-ph/9903066): the subtraction reaches about 82%. DADL computes it only under a
    loosened policy, and then labels it "policy modified" and "threshold crossed only after
    correction", with the raw 45% kept first."""
    (n_max, a_max), (n_min, a_min) = fringe["max"], fringe["min"]
    accidentals = 0.5 * (a_max + a_min)
    policy = Policy.default().allow(
        "bell_visibility", "accidental_subtraction", reason="reproducing Thompson quant-ph/9903066"
    )
    r = bell_visibility(n_max, n_min, accidentals=accidentals, policy=policy, window=Window(FRINGE_WIDTH_PS, "plan"))
    assert r.raw.value < BELL_VISIBILITY_THRESHOLD < r.corrected.value
    # Subtraction recovers the injected pair visibility (it is the pairs' V, not the data's).
    assert r.corrected.value == pytest.approx(FRINGE_SOURCE.visibility, abs=3 * r.corrected.stderr + 0.01)
    assert POLICY_MODIFIED in r.labels
    assert THRESHOLD_AFTER_CORRECTION in r.corrected_labels


# --- 2. The coincidence-time loophole (arXiv:1503.07573, arXiv:1309.0712) ---------------------------
#
# Christensen et al. built a classical source whose detection *times* depend on the local setting
# (two lasers per side, each polarized orthogonal to one setting, fired in staggered groups), so
# that a window centred on Alice's detections drops exactly the events that would spoil the
# correlation. Their CH analysis with moving windows showed B_CH up to 0.49 ("over 2700-sigma,
# assuming Gaussian statistics") for window radii between 1 and 3 us, and no significant violation
# with a predefined window locked to a synchronisation signal. Larsson et al. prove that fixed,
# locally assigned time slots close the loophole for CH/Eberhard-type inequalities.
#
# The model below is the same mechanism adapted to CHSH (DADL's analysis). Each emission carries a
# hidden type k. For k = 0..3 the outcomes are deterministic and spoil exactly one CHSH term (any
# local deterministic assignment spoils an odd number), and the setting-dependent delays put that
# term's two detections 2h apart while every other setting pair is at most h apart. A fraction G of
# emissions are uncorrelated noise at equal times. Every emission is detected (efficiency 1), so the
# only loophole is the coincidence time. With all events counted, S = 2 (1 - G) = 1.5.
#
# Event-centred window of radius r (h < r < 2h): the spoiled term's events are dropped, and
#   S = 4 (3/4)(1 - G) / ((3/4)(1 - G) + G) = 2.77 for G = 1/4.
# Fixed grid of width W >> h, emission times uniform against the grid: a pair Delta apart shares a
# slot with probability 1 - Delta / W, and summing over types gives exactly
#   S = 2 (1 - G) / (1 - (1 - G) h / W).

CTL_H_PS = 750  # the attacker's stagger, tuned to a 1 ns event-centred radius
CTL_RADIUS_PS = 1_000
CTL_NOISE = 0.25  # G
CTL_RATE = 5e4  # emissions per second per setting run
CTL_DURATION_S = 0.4
CTL_DET = DetectorParams(efficiency=1.0, dark_rate_cps=0.0, jitter_sigma_ps=20.0, dead_time_ps=22_000)


def _ctl_types() -> list[tuple[tuple[str, str], dict[str, int]]]:
    """For each CHSH term, a deterministic local assignment that spoils that term only."""
    out = []
    for spoiled in SIGNS:
        for a, a2, b, b2 in itertools.product((1, -1), repeat=4):
            outcome = {"a": a, "a'": a2, "b": b, "b'": b2}
            spoils = [key for key, sign in SIGNS.items() if sign * outcome[key[0]] * outcome[key[1]] < 0]
            if spoils == [spoiled]:
                out.append((spoiled, outcome))
                break
    assert len(out) == 4
    return out


CTL_TYPES = _ctl_types()


def _ctl_delays(spoiled: tuple[str, str], sa: str, sb: str) -> tuple[int, int]:
    """Alice's and Bob's delays: the spoiled term's settings 2h apart, every other pair within h."""
    da = 0 if sa == spoiled[0] else CTL_H_PS
    db = 2 * CTL_H_PS if sb == spoiled[1] else CTL_H_PS
    return da, db


def simulate_ctl_run(sa: str, sb: str, rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Recorded tags on A+, A-, B+, B- for one setting pair of the local coincidence-time model."""
    t0 = poisson_times(CTL_RATE, CTL_DURATION_S, rng)
    n = t0.size
    kind = rng.integers(0, 4, n)
    noise = rng.random(n) < CTL_NOISE
    flip = np.where(rng.random(n) < 0.5, 1, -1)  # a global sign flip keeps every output balanced
    out_a = np.empty(n, dtype=np.int64)
    out_b = np.empty(n, dtype=np.int64)
    t_a = t0.copy()
    t_b = t0.copy()
    for k, (spoiled, outcome) in enumerate(CTL_TYPES):
        sel = (kind == k) & ~noise
        da, db = _ctl_delays(spoiled, sa, sb)
        out_a[sel] = outcome[sa] * flip[sel]
        out_b[sel] = outcome[sb] * flip[sel]
        t_a[sel] += da
        t_b[sel] += db
    out_a[noise] = np.where(rng.random(int(noise.sum())) < 0.5, 1, -1)
    out_b[noise] = np.where(rng.random(int(noise.sum())) < 0.5, 1, -1)
    return {
        "A+": detect(t_a[out_a > 0], CTL_DET, CTL_DURATION_S, rng),
        "A-": detect(t_a[out_a < 0], CTL_DET, CTL_DURATION_S, rng),
        "B+": detect(t_b[out_b > 0], CTL_DET, CTL_DURATION_S, rng),
        "B-": detect(t_b[out_b < 0], CTL_DET, CTL_DURATION_S, rng),
    }


@pytest.fixture(scope="module")
def ctl_runs():
    rng = np.random.default_rng(1503_07573)
    return {key: simulate_ctl_run(*key, rng) for key in SIGNS}


def event_centred_s(runs: dict, radius_ps: int) -> tuple[float, float]:
    """CHSH S with a window of +-radius centred on each of Alice's detections (the loophole).

    This is the data-dependent coincidence rule DADL does not offer; it is written out here only to
    show what it does to the model above.
    """
    s, var = 0.0, 0.0
    for key, sign in SIGNS.items():
        tags = runs[key]
        bob = np.concatenate([tags["B+"], tags["B-"]])
        bob_out = np.concatenate([np.ones(tags["B+"].size), -np.ones(tags["B-"].size)])
        order = np.argsort(bob, kind="stable")
        bob, bob_out = bob[order], bob_out[order]
        products = []
        for ch, sign_a in (("A+", 1), ("A-", -1)):
            t = tags[ch]
            i = np.clip(np.searchsorted(bob, t), 1, bob.size - 1)
            nearest = np.where(np.abs(bob[i] - t) < np.abs(bob[i - 1] - t), i, i - 1)
            hit = np.abs(bob[nearest] - t) <= radius_ps
            products.append(sign_a * bob_out[nearest[hit]])
        p = np.concatenate(products)
        e = float(p.mean())
        s += sign * e
        var += (1 - e * e) / p.size
    return s, math.sqrt(var)


def ctl_expected_event_centred(radius_ps: float) -> float:
    """Expected S of the model for an event-centred radius, jitter neglected."""
    g, h = CTL_NOISE, CTL_H_PS
    if radius_ps < h:  # only the type with no shift on this term, and the noise, coincide
        good, bad = (1 - g) / 4, 0.0
    elif radius_ps < 2 * h:  # the spoiled type is dropped
        good, bad = 3 * (1 - g) / 4, 0.0
    else:  # everything coincides
        good, bad = 3 * (1 - g) / 4, (1 - g) / 4
    return 4 * (good - bad) / (good + bad + g)


def ctl_expected_grid(width_ps: float) -> float:
    """Expected S of the model on DADL's fixed grid of width W (W >= 2h, jitter neglected)."""
    g, h = CTL_NOISE, CTL_H_PS
    return 2 * (1 - g) / (1 - (1 - g) * h / width_ps)


def test_event_centred_window_fakes_a_violation_with_a_local_model(ctl_runs):
    """arXiv:1503.07573 (Christensen et al., NIST and Illinois): a classical source with
    setting-dependent detection times shows a large apparent violation with windows centred on the
    detections ("over 2700-sigma, assuming Gaussian statistics" in their CH analysis), but only for
    radii in a band set by the source's stagger (1-3 us there).

    Here, adapted to CHSH: S = 2.77 for h < r < 2h, tens of Gaussian sigma above 2, from a model
    whose S over all events is 1.5; outside the band no violation.
    """
    s, err = event_centred_s(ctl_runs, CTL_RADIUS_PS)
    assert s == pytest.approx(ctl_expected_event_centred(CTL_RADIUS_PS), abs=4 * err)
    assert s > 2.7 and (s - 2) / err > 30  # a "violation" by tens of Gaussian sigma
    for radius in (CTL_H_PS // 2, 3 * CTL_H_PS):
        s_out, err_out = event_centred_s(ctl_runs, radius)
        assert s_out == pytest.approx(ctl_expected_event_centred(radius), abs=4 * err_out)
        assert s_out < 2


def test_planned_fixed_grid_shows_no_violation_for_the_same_tags(ctl_runs):
    """arXiv:1503.07573: with a predefined window locked to a synchronisation signal "we do not see a
    statistically significant Bell violation"; arXiv:1309.0712: fixed, locally assigned time slots
    avoid the coincidence-time loophole.

    DADL's clock-locked grid with a planned 10 ns width (W >> h) on the very tags that gave S = 2.77
    returns the model's grid expectation, 1.59, below 2 by many sigma.
    """
    width = 10_000
    result = chsh(ctl_runs, Window(width, "plan"))
    expected = ctl_expected_grid(width)
    assert result.raw.value == pytest.approx(expected, abs=4 * result.raw.stderr)
    assert result.raw.value + 3 * result.raw.stderr < 2.0


def test_data_dependent_window_is_refused_on_loophole_data(ctl_runs):
    """arXiv:1503.07573: the moving (event-centred) window is what creates the violation, so a
    window whose source is the data is refused before anything is counted."""
    with pytest.raises(GuardrailRefusal):
        chsh(ctl_runs, Window(2 * CTL_RADIUS_PS, "data"))


def test_fixed_grid_alone_does_not_close_the_loophole_for_chsh():
    """The limit, stated so nobody over-reads the fixed grid. Larsson et al. (arXiv:1309.0712)
    prove closure for CH/Eberhard inequalities normalised without post-selection; DADL's CHSH is
    normalised by coincidences (the fair-coincidence assumption). With a grid only a little wider
    than the stagger, the same model reads S > 2 on a planned grid, so the result must carry the
    fair-sampling and fair-coincidence label. (A source that knew the grid's phase could do more.)
    """
    rng = np.random.default_rng(1309_0712)
    runs = {key: simulate_ctl_run(*key, rng) for key in SIGNS}
    width = 2 * CTL_H_PS + 100
    result = chsh(runs, Window(width, "plan"))
    expected = ctl_expected_grid(width)
    assert expected > 2.0
    assert result.raw.value == pytest.approx(expected, abs=4 * result.raw.stderr)
    assert FAIR_SAMPLING in result.labels


# --- 3. NIST calibration comparison (Polyakov and Migdall, SPIE 6583, 65830A, 2007) -----------------
#
# "We have reported a relative standard uncertainty for the correlated-photon method of 0.18 % (k=1)
# and for the substitution method of 0.17 % (k=1) ... the two methods agree to within the
# uncertainty of the comparison, with the difference being 0.14(14) %." The same paper proposes a
# low-cost substitution procedure with 0.5 % (k=1).
#
# Here: a signal detector with injected efficiency 0.6 is calibrated twice. (1) Correlated photons
# (Klyshko): DADL's heralding calibration on simulated SPDC pairs, with statistics chosen to land in
# the 0.18 % (k=1) regime. Its measurand is "as detected" (window capture, dead time, herald dark
# counts); the known detector model converts it to the detector efficiency, as NIST corrects for
# those effects. (2) Substitution: the same detector and a reference detector of known efficiency
# (taken as exact here) take turns on one attenuated, Poisson laser beam; their dark- and
# dead-time-corrected rates give the efficiency ratio.

ETA_SIGNAL = 0.6
NIST_HERALD = DetectorParams(efficiency=0.8, dark_rate_cps=20.0, jitter_sigma_ps=20.0, dead_time_ps=22_000)
NIST_SIGNAL = DetectorParams(
    efficiency=ETA_SIGNAL, dark_rate_cps=20.0, jitter_sigma_ps=20.0, dead_time_ps=22_000, delay_ps=1_234
)
NIST_PAIR_RATE = 2e5
NIST_DURATION_S = 1.5
NIST_WIDTH_PS = 5_000
#: The k = 1 relative uncertainties NIST reports for the correlated-photon and substitution methods.
NIST_KLYSHKO_REL_K1 = 0.0018
NIST_SUBSTITUTION_REL_K1 = 0.0017


@pytest.fixture(scope="module")
def klyshko():
    herald, signal = simulate_pairs(
        NIST_PAIR_RATE, NIST_HERALD, NIST_SIGNAL, NIST_DURATION_S, np.random.default_rng(6583)
    )
    cal = calibrate_heralding(
        herald,
        signal,
        (ChannelId("sim", 0), ChannelId("sim", 1)),
        width_ps=NIST_WIDTH_PS,
        delay_signal_ps=NIST_SIGNAL.delay_ps,
        dead_time_ps=NIST_SIGNAL.dead_time_ps,
        backend="sim",
        seed=1,
    )
    as_detected = heralding_ratio_expected(
        NIST_PAIR_RATE,
        NIST_HERALD.efficiency,
        NIST_SIGNAL.efficiency,
        NIST_HERALD.dark_rate_cps,
        NIST_SIGNAL.dark_rate_cps,
        NIST_HERALD.dead_time_ps * 1e-12,
        NIST_SIGNAL.dead_time_ps * 1e-12,
        math.hypot(NIST_HERALD.jitter_sigma_ps, NIST_SIGNAL.jitter_sigma_ps),
        NIST_WIDTH_PS,
    )
    return cal.estimates["heralding_efficiency"], as_detected


def substitution(rng: np.random.Generator, eta_ref: float = 0.7, flux: float = 1e6, duration_s: float = 1.0):
    """Substitution against a reference detector of known efficiency: eta = eta_ref x (R_dut / R_ref),
    each rate dark-subtracted (from its own dark run) and corrected for its known non-paralyzable
    dead time. Returns the estimate and its relative standard uncertainty (counting only)."""
    reference = DetectorParams(efficiency=eta_ref, dark_rate_cps=20.0, jitter_sigma_ps=20.0, dead_time_ps=22_000)

    def incident(det: DetectorParams) -> tuple[float, float]:
        tags = detect(poisson_times(flux, duration_s, rng), det, duration_s, rng)
        dark = detect(np.empty(0), det, duration_s, rng)
        tau = det.dead_time_ps * 1e-12
        rate = tags.size / duration_s / (1 - tags.size / duration_s * tau)  # undo non-paralyzable loss
        dark_rate = dark.size / duration_s
        # Poisson counting error, scaled by the dead-time correction's derivative 1 / (1 - m tau)^2.
        m = tags.size / duration_s
        var = (m / duration_s) / (1 - m * tau) ** 4 + dark_rate / duration_s
        return rate - dark_rate, var

    r_ref, v_ref = incident(reference)
    r_dut, v_dut = incident(NIST_SIGNAL)
    eta = eta_ref * r_dut / r_ref
    rel = math.sqrt(v_ref / r_ref**2 + v_dut / r_dut**2)
    return eta, rel


def test_klyshko_calibration_recovers_the_injected_efficiency_at_nist_precision(klyshko):
    """Polyakov and Migdall (SPIE 6583, 65830A): correlated-photon calibration with a relative
    standard uncertainty of 0.18 % (k = 1).

    DADL's heralding calibration matches the as-detected truth, its 95% interval is in the same
    precision regime (a half-width of 1.96 x 0.18 % within a factor 1.5), and after the known
    detector-model factors it recovers the injected detector efficiency to within that interval.
    """
    iv, as_detected = klyshko
    # Within three standard uncertainties, not "inside the 95% interval": one seed is not a coverage
    # check (60 replicates of this setup: mean 0.59683(12) for a truth of 0.59690, 58/60 covered).
    sigma = 0.5 * (iv.high - iv.low) / 1.96
    assert abs(iv.value - as_detected) < 3 * sigma
    half_rel = 0.5 * (iv.high - iv.low) / iv.value
    target = 1.96 * NIST_KLYSHKO_REL_K1
    assert target / 1.5 < half_rel < target * 1.5
    factor = as_detected / ETA_SIGNAL  # window capture x live fraction x herald dilution
    assert abs(factor - 1) < 0.02  # small, as in a metrology setup
    eta = iv.value / factor
    assert abs(eta - ETA_SIGNAL) / ETA_SIGNAL < 3 * NIST_KLYSHKO_REL_K1


def test_klyshko_and_substitution_agree_like_the_nist_comparison(klyshko):
    """Polyakov and Migdall (SPIE 6583, 65830A): the two primary methods differ by 0.14(14) %,
    within the comparison's uncertainty, with 0.18 % and 0.17 % (k = 1) each.

    DADL's Klyshko-style result and a substitution against a reference detector, on the same
    simulated detector, differ by less than three combined standard uncertainties, and by under
    0.6 % (the NIST difference plus three of its uncertainties is 0.56 %).
    """
    iv, as_detected = klyshko
    eta_k = iv.value / (as_detected / ETA_SIGNAL)
    rel_k = 0.5 * (iv.high - iv.low) / iv.value / 1.96
    eta_s, rel_s = substitution(np.random.default_rng(65830))
    assert rel_s < 2 * NIST_SUBSTITUTION_REL_K1  # the substitution run is in the same regime
    diff = (eta_k - eta_s) / eta_s
    combined = math.hypot(rel_k, rel_s)
    assert abs(diff) < 3 * combined
    assert abs(diff) < 0.0014 + 3 * 0.0014
    assert abs(eta_s - ETA_SIGNAL) / ETA_SIGNAL < 3 * rel_s


def test_substitution_helper_matches_the_dead_time_closed_form():
    """Guard on the reference method itself: its dead-time model is the simulator's (non-paralyzable)."""
    tau = NIST_SIGNAL.dead_time_ps * 1e-12
    measured = nonparalyzable_measured_rate(6e5, tau)
    assert measured / (1 - measured * tau) == pytest.approx(6e5)
