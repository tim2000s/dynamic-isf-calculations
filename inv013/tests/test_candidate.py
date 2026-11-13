"""The comparator has to be the stock equation, and the candidate one-sided."""
from __future__ import annotations

import numpy as np
import pytest

from inv008.dynisf import isf_v1, quartic_isf
from inv013.candidate import (ANCHOR_GLUCOSE, _QUADRATIC, candidate, curve_isf, curve_k, curve_profile, curve_shape,
                              own_levels, published, stock_v1, trio_log, trio_sigmoid, units_to_target)

BG = np.array([60.0, 80.0, 100.0, 150.0, 200.0, 300.0])


def test_stock_equation_matches_the_replay_with_cap_and_velocity_off():
    assert float(stock_v1(100.0, 40.0, 55.0)) == pytest.approx(1800 / (40 * np.log(100 / 55 + 1)))
    np.testing.assert_allclose(stock_v1(BG, 40.0, 75.0),
                               isf_v1(BG, 40.0, divisor=75.0, velocity=1.0, cap=1e9))


def test_trio_forms_return_the_profile_value_where_their_ratio_is_one():
    assert float(trio_sigmoid(100.0, 50.0)) == pytest.approx(50.0)
    # Unclamped, the profile value cancels and the stock equation is left.
    assert float(trio_log(100.0, 40.0, 50.0, factor=1.0, lo=0.01, hi=100.0)) == pytest.approx(
        float(stock_v1(100.0, 40.0, 55.0)))


def test_candidate_eases_below_the_pivot_and_nowhere_else():
    out = candidate(BG, 50.0, 40.0, 40.0)
    assert np.all(out[BG >= 100.0] == 50.0)
    assert np.all(np.diff(out[BG <= 100.0]) < 0) and out[0] > 50.0


def test_dose_ratio_is_clamped_and_the_high_option_takes_the_stock_shape():
    assert float(candidate(150.0, 50.0, 80.0, 40.0)) == pytest.approx(50.0 / 1.2)
    assert float(candidate(150.0, 50.0, 10.0, 40.0)) == pytest.approx(50.0 / 0.7)
    high = candidate(BG, 50.0, 40.0, 40.0, high_from=200.0)
    assert high[4] == pytest.approx(50.0)
    assert high[5] / high[4] == pytest.approx(float(stock_v1(300.0, 1.0, 55.0) / stock_v1(200.0, 1.0, 55.0)))


def test_curve_is_the_published_quartic_carrying_the_entered_product_at_the_anchor():
    assert ANCHOR_GLUCOSE == 140.0 and float(curve_k(ANCHOR_GLUCOSE)) == pytest.approx(1915.0)
    assert float(curve_profile(ANCHOR_GLUCOSE, 50.0)) == pytest.approx(50.0)
    # The abstract: 35 at 180 mg/dL and 79 at 100, doses about halved at 120 and about doubled above 340.
    assert float(curve_k(100.0) / curve_k(180.0)) == pytest.approx(79.0 / 35.0, rel=0.03)
    assert float(curve_shape(180.0) / curve_shape(120.0)) == pytest.approx(0.5, abs=0.07)
    assert float(curve_shape(180.0) / curve_shape(340.0)) == pytest.approx(2.0, abs=0.1)
    assert float(curve_k(60.0)) == float(curve_k(100.0)) and float(curve_k(450.0)) == float(curve_k(380.0))
    assert np.all(np.diff(curve_k(np.linspace(100.0, 380.0, 100))) < 0)
    # The poster's quadratic gives the abstract's 79 at 100 mg/dL and 65 near 60.
    assert np.polyval(_QUADRATIC[::-1], [60.0, 100.0]) == pytest.approx([65.0, 79.0], abs=0.5)
    assert float(curve_isf(200.0, 40.0, factor=0.5)) == pytest.approx(2 * float(curve_isf(200.0, 40.0)))
    # The two links agree for a person whose profile x daily dose is the entered product.
    assert float(curve_profile(250.0, 1915.0 / 40.0)) == pytest.approx(float(curve_isf(250.0, 40.0)))


def test_units_to_target_reads_the_sensitivity_at_the_current_glucose():
    isf = lambda g: curve_isf(g, 40.0)
    assert units_to_target(isf, 200.0) == pytest.approx(100.0 / float(isf(200.0)))


def test_own_levels_reads_the_positive_insulin_table_and_not_the_one_above_it():
    own = own_levels()
    # Loop at 100 mg/dL is 2132 over every event and 6974 over events with positive insulin used.
    assert own["Loop"] == (842, own["Loop"][1]) and own["Loop"][1][100] == (6974.0, 6659.0, 7194.0)
    assert 250 not in own["Trio"][1] and own["no algorithm (REPLACE-BG)"][1][200][0] == 3283.0
    assert published(100.0) == pytest.approx(79.0, abs=0.1) and published(180.0) == pytest.approx(35.3, abs=0.05)


def test_refit_returns_a_known_quartic_and_the_scale_of_a_scaled_published_curve():
    import pandas as pd
    from inv013.refit import refit
    g = np.arange(100, 310, 10)
    truth = lambda x: 9000 - 60 * x + 0.2 * x ** 2 - 3e-4 * x ** 3 + 2e-7 * x ** 4
    wide = pd.DataFrame(np.tile(truth(g.astype(float)), (12, 1)), columns=g)
    r = refit(wide, n_boot=5)
    assert r["level"][200][0] == pytest.approx(truth(200.0)) and r["miss_refit"] < 1e-6
    assert r["shape"][100][0] == pytest.approx(truth(100.0) / truth(180.0))
    np.testing.assert_allclose(r["coef"][::-1], [9000, -60, 0.2, -3e-4, 2e-7], rtol=1e-4)
    # People who follow the published curve at 60 U a day give that scale back and no miss.
    same = pd.DataFrame(np.tile(60.0 * quartic_isf(g.astype(float)), (12, 1)), columns=g)
    r = refit(same, n_boot=5)
    assert r["scale"][0] == pytest.approx(60.0) and r["miss_published"] < 1e-9


def test_profile_link_gives_the_table_the_device_shadow_code_is_tested_against():
    # The same glucose and sensitivity pairs sit in the AndroidAPS and Trio tests (profile ISF 50).
    table = {40: 79.2753, 60: 79.2753, 100: 79.2753, 120: 62.4565, 140: 50.0,
             180: 34.678, 250: 24.3772, 300: 20.7823, 380: 15.0273, 450: 15.0273}
    for glucose, isf in table.items():
        assert float(curve_profile(glucose, 50.0)) == pytest.approx(isf, abs=1e-4)
