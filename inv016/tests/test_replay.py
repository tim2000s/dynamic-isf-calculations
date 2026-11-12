import numpy as np

from inv016.replay import ARMS, arms, hours_since_meal


def test_arms_at_known_points():
    bg = np.array([180.0, 140.0, 100.0])
    r = arms(bg, np.full(3, 38.3), 50.0, 38.3, trio=False)      # 50 x 38.3 = 1915
    assert set(r) == set(ARMS)
    assert abs(r["Profile link"][1] - 1) < 1e-9                 # the profile holds at the anchor, 140
    assert r["Profile link"][0] > 1 > r["Profile link"][2]      # read at the reading: more above, less below
    assert np.allclose(r["TDD link"], r["Profile link"])         # the links agree at 1915
    assert np.allclose(r["Average-dose link"], r["TDD link"])
    assert r["Profile link"][2] < 0.7 < r["Logarithmic as shipped"][2]
    double = arms(bg, np.full(3, 76.6), 50.0, 38.3, trio=False)
    assert np.allclose(double["TDD link"], 2 * r["TDD link"])    # the dose reaches the TDD arm
    assert np.allclose(double["Profile link"], r["Profile link"])


def test_hours_since_meal():
    ts = np.arange(6) * 1800.0
    hrs, meals = hours_since_meal(ts, [0, 0, 30, 25, np.nan, 28])    # a 3 g rise is no entry
    assert meals == 1 and np.isnan(hrs[:2]).all()
    assert list(hrs[2:]) == [0.0, 0.5, 1.0, 1.5]
