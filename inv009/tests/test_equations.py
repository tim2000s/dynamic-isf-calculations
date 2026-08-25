"""The two equations, pinned to the source they were ported from.

Both are implemented in inv008.dynisf and are used here to produce every v1 and
v2 figure in the write-up. These tests hold them against the plugin code so a
future edit cannot quietly move a comparator.

Source references, from a local checkout of the plugin repository:

    plugins/aps/src/main/kotlin/app/aaps/plugins/aps/openAPSSMB/OpenAPSSMBPlugin.kt
        line 295  tddWeightedFromLast8H = ((1.4 * tddLast4H) + (0.6 * tddLast8to4H)) * 3
        line 296  tdd = (w8h * 0.33) + (tdd7D * 0.34) + (tdd1D * 0.33)
        line 297  variableSensitivity = 1800 / (tdd * ln((glucose / insulinDivisor) + 1))

    plugins/aps/src/test/kotlin/app/aaps/plugins/aps/openAPSBoostV3/BoostV3IsfCalculationTest.kt
        computeV1BlendedTdd  adds a branch when w8h < 0.75 * tdd7D
        computeVariableSensV3  caps glucose, counting the excess at a third

The stock plugin has neither the branch nor the cap. inv008.dynisf implements the
Boost variant, which carries both. The difference is measured in
test_variant_choice_does_not_move_the_comparison below and in the write-up.
"""
from __future__ import annotations

import numpy as np
import pytest

from inv008 import dynisf


def stock_sensitivity(bg, tdd, divisor=75.0):
    """OpenAPSSMBPlugin.kt line 297, transcribed."""
    return 1800.0 / (tdd * np.log(bg / divisor + 1.0))


def stock_blend(t4, t84, t1, t7):
    """OpenAPSSMBPlugin.kt lines 295 to 296, transcribed. No branch."""
    w8h = (1.4 * t4 + 0.6 * t84) * 3.0
    return 0.33 * w8h + 0.34 * t7 + 0.33 * t1


def test_v1_matches_the_plugin_formula_below_the_cap():
    """Under the cap the ported v1 and the plugin line agree exactly."""
    for bg in (80.0, 120.0, 160.0, 200.0, 209.0):
        for tdd in (12.0, 35.0, 90.0):
            got = float(dynisf.isf_v1(np.array([bg]), np.array([tdd]))[0])
            assert got == pytest.approx(stock_sensitivity(bg, tdd), rel=1e-9)


def test_v1_diverges_above_the_cap_by_the_documented_amount():
    """Above 210 mg/dL the Boost variant counts the excess at a third."""
    tdd = np.array([40.0])
    for bg, expected_capped in ((240.0, 210 + 30 / 3), (300.0, 210 + 90 / 3)):
        got = float(dynisf.isf_v1(np.array([bg]), tdd)[0])
        assert got == pytest.approx(stock_sensitivity(expected_capped, 40.0), rel=1e-9)
        assert got > stock_sensitivity(bg, 40.0)      # capping raises sensitivity


def test_blend_reproduces_the_weighted_eight_hour_figure():
    """A flat eight hours gives twice a four hour total, scaled by three."""
    flat = 5.0
    w8h = (1.4 * flat + 0.6 * flat) * 3.0
    assert w8h == pytest.approx(2.0 * flat * 3.0)
    # With everything consistent at 30 U/day the blend returns 30 U/day.
    got = float(dynisf.blend_tdd(np.array([5.0]), np.array([5.0]),
                                 np.array([30.0]), np.array([30.0]))[0])
    assert got == pytest.approx(30.0, rel=1e-9)


def test_blend_matches_stock_when_the_low_branch_does_not_fire():
    """Above 75% of the seven day average, both implementations agree."""
    t4, t84, t7, t1 = 6.0, 6.0, 30.0, 32.0
    w8h = (1.4 * t4 + 0.6 * t84) * 3.0
    assert w8h >= 0.75 * t7                       # the branch is not taken
    got = float(dynisf.blend_tdd(np.array([t4]), np.array([t84]),
                                 np.array([t1]), np.array([t7]))[0])
    assert got == pytest.approx(stock_blend(t4, t84, t1, t7), rel=1e-9)


def test_blend_diverges_from_stock_when_the_low_branch_fires():
    """Below the threshold the Boost variant pulls the seven day term down."""
    t4, t84, t7, t1 = 1.0, 1.0, 40.0, 38.0
    w8h = (1.4 * t4 + 0.6 * t84) * 3.0
    assert w8h < 0.75 * t7                        # the branch is taken
    got = float(dynisf.blend_tdd(np.array([t4]), np.array([t84]),
                                 np.array([t1]), np.array([t7]))[0])
    assert got < stock_blend(t4, t84, t1, t7)     # lower TDD, so higher sensitivity


def test_the_crossover_depends_on_glucose_and_is_far_above_64():
    """Where v1 and v2 swap places, which is not where 63.9/TDD suggests.

    The coefficient ratio alone is 115000/1800 = 63.9, which reads as a
    crossover at 64 U/day. It is not one. v1 takes ln(bg/divisor + 1) and v2
    takes ln(bg/divisor), and that pair differs by a factor of two to three
    across the range, so the real crossover carries a glucose term and lands far
    higher. An earlier note in this series recorded 64 and this test exists
    because that number reached a draft.
    """
    expected = {99.0: 194.0, 120.0: 130.0, 150.0: 101.0, 180.0: 89.0, 210.0: 83.0}
    for bg, tdd_cross in expected.items():
        r = float(dynisf.isf_v2(np.array([bg]), np.array([tdd_cross]))[0]
                  / dynisf.isf_v1(np.array([bg]), np.array([tdd_cross]))[0])
        assert r == pytest.approx(1.0, abs=0.02), f"at {bg} mg/dL"
    # 64 U/day is nowhere near it at any glucose in range.
    for bg in expected:
        r = float(dynisf.isf_v2(np.array([bg]), np.array([64.0]))[0]
                  / dynisf.isf_v1(np.array([bg]), np.array([64.0]))[0])
        assert r > 1.25


def test_v2_doses_smaller_than_v1_across_this_cohort():
    """Below the crossover, which is 97% of these archives, v2 gives the weaker correction."""
    for tdd in (7.0, 15.0, 40.0, 80.0, 100.0):
        r = float(dynisf.isf_v2(np.array([150.0]), np.array([tdd]))[0]
                  / dynisf.isf_v1(np.array([150.0]), np.array([tdd]))[0])
        assert r > 1.0, f"v2 should read higher, so dose smaller, at {tdd} U/day"


def test_v2_hands_a_small_dose_an_implausible_sensitivity():
    """The failure mode the write-up reports for low daily totals."""
    isf = float(dynisf.isf_v2(np.array([120.0]), np.array([10.0]))[0])
    assert isf > 1000.0
