"""The comparator has to be the stock equation, and the candidate one-sided."""
from __future__ import annotations

import numpy as np
import pytest

from inv008.dynisf import isf_v1
from inv013.candidate import candidate, stock_v1, trio_log, trio_sigmoid

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
