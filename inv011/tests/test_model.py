import numpy as np
import pytest

from inv011.model import (
    activity_exponential,
    dynamic_isf,
    first_crossing_minutes,
    remaining_exponential,
)
from inv011.run import temp_basal_example


def test_remaining_is_monotone_and_mass_is_conserved():
    t = np.arange(0.0, 361.0)
    for peak in (45.0, 70.0):
        remaining = remaining_exponential(t, 360.0, peak)
        assert remaining[0] == pytest.approx(1.0)
        assert remaining[-1] == pytest.approx(0.0)
        assert np.all(np.diff(remaining) <= 1e-12)
        activity = activity_exponential(t, 360.0, peak)
        assert np.trapezoid(activity, t) == pytest.approx(1.0, abs=2e-4)


def test_activity_peak_matches_parameter():
    t = np.arange(0.0, 361.0, 0.1)
    for peak in (45.0, 70.0):
        activity = activity_exponential(t, 360.0, peak)
        assert t[np.argmax(activity)] == pytest.approx(peak, abs=0.1)


def test_divisor_55_requires_more_insulin_than_75():
    for equation in ("v1", "v2"):
        correct = float(dynamic_isf(equation, 150.0, 40.0, 75.0))
        wrong = float(dynamic_isf(equation, 150.0, 40.0, 55.0))
        assert correct / wrong > 1.0


def test_v2_divisor_error_is_larger_than_v1_near_target():
    v1 = float(dynamic_isf("v1", 99.0, 40.0, 75.0) /
               dynamic_isf("v1", 99.0, 40.0, 55.0))
    v2 = float(dynamic_isf("v2", 99.0, 40.0, 75.0) /
               dynamic_isf("v2", 99.0, 40.0, 55.0))
    assert v1 == pytest.approx(1.2234548052)
    assert v2 == pytest.approx(2.1171450790)
    assert v2 > v1


def test_crossing_interpolates():
    assert first_crossing_minutes([0, 5, 10], [0.8, 0.9, 1.1]) == pytest.approx(7.5)
    assert first_crossing_minutes([0, 5], [1.1, 1.2]) is None


def test_temp_basal_example_conserves_delivery_and_separates_curves():
    result = temp_basal_example()
    assert sum(result["delivery_u"]) == pytest.approx(2.0)
    v1 = result["by_equation"]["v1"]["selected_times"]
    assert v1["120"]["assumed_iob_u"] > v1["120"]["reference_iob_u"]
    assert v1["300"]["remaining_effect_ratio"] > 1.0
