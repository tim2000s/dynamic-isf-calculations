"""The comparator has to be the stock equation, and the candidate one-sided."""
from __future__ import annotations

import numpy as np
import pytest

from inv008.dynisf import isf_v1
from inv013.candidate import (_QUADRATIC, candidate, curve_isf, curve_k, curve_profile, curve_shape,
                              stock_v1, trio_log, trio_sigmoid, units_to_target)

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


def test_curve_is_the_published_quartic_carrying_the_entered_product_at_180():
    assert float(curve_k(180.0)) == pytest.approx(1915.0)
    assert float(curve_isf(180.0, 1915.0 / 35.0)) == pytest.approx(35.0)  # the abstract's reference
    assert float(curve_isf(100.0, 1915.0 / 35.0)) == pytest.approx(79.0, rel=0.03)
    # The abstract: doses about halved at 120 mg/dL and about doubled above 340.
    assert 1.0 / float(curve_shape(120.0)) == pytest.approx(0.5, abs=0.07)
    assert 1.0 / float(curve_shape(340.0)) == pytest.approx(2.0, abs=0.1)
    assert float(curve_k(60.0)) == float(curve_k(100.0)) and float(curve_k(450.0)) == float(curve_k(380.0))
    assert np.all(np.diff(curve_k(np.linspace(100.0, 380.0, 100))) < 0)
    # The poster's quadratic gives the abstract's 79 at 100 mg/dL and 65 near 60.
    assert np.polyval(_QUADRATIC[::-1], [60.0, 100.0]) == pytest.approx([65.0, 79.0], abs=0.5)
    assert float(curve_isf(200.0, 40.0, factor=0.5)) == pytest.approx(2 * float(curve_isf(200.0, 40.0)))
    # The two links agree for a person whose profile x daily dose is the entered product.
    assert float(curve_profile(250.0, 1915.0 / 40.0)) == pytest.approx(float(curve_isf(250.0, 40.0)))


def test_units_to_target_orders_the_three_readings_of_a_falling_curve():
    isf = lambda g: curve_isf(g, 40.0)
    start, mid, path = (units_to_target(isf, 200.0, r) for r in ("start", "midpoint", "path"))
    assert start == pytest.approx(100.0 / float(isf(200.0)))
    assert mid == pytest.approx(100.0 / float(isf(150.0)))
    assert max(mid, path) < start and path == pytest.approx(mid, rel=0.05)
    # A flat sensitivity gives the same answer however it is read.
    flat = lambda g: np.full_like(np.asarray(g, dtype=float), 50.0)
    assert units_to_target(flat, 200.0, "path") == pytest.approx(units_to_target(flat, 200.0))
