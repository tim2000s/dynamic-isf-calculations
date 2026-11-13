import numpy as np
import pandas as pd

from inv019.outcome_test import (ARMS, LIMITS, MARGIN, RULES, blended_tdd, matched_factor,
                                 outcomes, paired, verdict)

STEPS_PER_DAY = 288


def test_the_blend_waits_for_seven_days_and_returns_a_steady_dose_unchanged():
    per_step = 48.0 / STEPS_PER_DAY
    cum = list(np.arange(7 * STEPS_PER_DAY + 1) * per_step)
    assert blended_tdd(cum[:-1]) is None
    assert abs(blended_tdd(cum) - 48.0) < 1e-9


def test_the_blend_weights_the_last_four_hours_most():
    quiet = np.full(8 * STEPS_PER_DAY, 48.0 / STEPS_PER_DAY)
    late = quiet.copy()
    late[-12:] += 0.5          # 6 U in the last hour
    early = quiet.copy()
    early[-12 * 20:-12 * 19] += 0.5     # the same 6 U twenty hours back
    tdd = lambda u: blended_tdd([0.0] + list(np.cumsum(u)))
    assert tdd(late) - 48.0 > 3 * (tdd(early) - 48.0) > 0


def test_each_rule_returns_the_profile_where_it_is_anchored():
    isf, tdd = 50.0, 40.0
    f = matched_factor(isf, tdd)
    at = lambda rule, g: float(RULES[rule](g, isf, tdd, f, 110.0))
    assert at("profile", 300.0) == isf
    assert abs(at("candidate", 140.0) - isf) < 1e-9
    assert abs(at("candidate_180", 180.0) - isf) < 1e-9
    assert abs(at("log_matched", 140.0) - isf) < 1e-9
    assert abs(at("sigmoid", 110.0) - isf) < 1e-6
    assert at("flat", 250.0) == isf and at("flat", 80.0) > isf


def test_the_candidate_falls_with_glucose_and_the_limits_hold_it():
    isf, at = 50.0, lambda rule, g: float(RULES[rule](g, 50.0, 40.0, 1.0, 110.0))
    assert at("candidate", 100.0) > at("candidate", 200.0) > at("candidate", 380.0) == at("candidate", 450.0)
    for g in (60.0, 100.0, 140.0, 250.0, 380.0):
        assert isf / LIMITS[1] - 1e-9 <= at("candidate_limits", g) <= isf / LIMITS[0] + 1e-9
        assert isf / LIMITS[1] - 1e-9 <= at("log_trio", g) <= isf / LIMITS[0] + 1e-9
    assert at("candidate_limits", 380.0) > at("candidate", 380.0)


def test_every_arm_names_a_rule_and_the_placebo_is_the_profile_moved_by_one_per_cent():
    assert all(rule in RULES for rule, _, _ in ARMS.values())
    assert ARMS["placebo"] == ("profile", 1.01, 1.0)


def test_outcomes_count_time_by_threshold():
    o = outcomes([50, 60, 100, 200, 300])
    assert (o["tbr54"], o["tbr70"], o["tir"], o["tar180"], o["tar250"]) == (20.0, 40.0, 20.0, 40.0, 20.0)


def test_the_verdict_weighs_size_as_well_as_sign():
    rng = np.random.default_rng(1)
    m = MARGIN["tbr54"]
    assert verdict(m * 4 + rng.normal(0, 0.05, 40), m) == "harm"
    assert verdict(rng.normal(0, 0.05, 40), m) == "cleared"
    assert verdict(rng.normal(0, 3.0, 40), m) == "unresolved"
    # a difference that is certain and smaller than the margin is not harm
    assert verdict(m / 5 + rng.normal(0, 0.01, 40), m) == "cleared"


def test_a_persons_difference_is_the_mean_over_the_seeds_where_both_runs_completed():
    row = lambda who, seed, arm, v, done=True: dict(subject=who, seed=seed, arm=arm, tbr54=v, completed=done)
    runs = pd.DataFrame([row("a", 1, "x", 1.0), row("a", 1, "ref", 0.0), row("a", 2, "x", 4.0), row("a", 2, "ref", 1.0),
                         row("b", 1, "x", 2.0), row("b", 1, "ref", 1.5), row("b", 2, "x", 9.0, False), row("b", 2, "ref", 0.0)])
    assert list(paired(runs, "x", "ref", "tbr54")) == [2.0, 0.5]
    assert list(paired(runs[runs.seed == 1], "x", "ref", "tbr54")) == [1.0, 0.5]
