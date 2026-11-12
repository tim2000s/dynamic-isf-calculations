import numpy as np

from inv016.replay import ARMS, arms, hours_since_meal


def test_arms_at_known_points():
    bg = np.array([180.0, 140.0, 100.0])
    r = arms(bg, np.full(3, 38.3), 50.0, 38.3, trio=False)      # 50 x 38.3 = 1915
    assert set(r) == set(ARMS)
    assert abs(r["Profile link, midpoint"][0] - 1) < 1e-9       # midpoint of 180 and 100 is the anchor, 140
    assert abs(r["Profile link, start"][1] - 1) < 1e-9
    assert np.allclose(r["TDD link"], r["Profile link, midpoint"])   # the links agree at 1915
    assert np.allclose(r["Average-dose link"], r["TDD link"])
    assert r["Profile link, midpoint"][2] < 0.7 < r["Logarithmic as shipped"][2]
    double = arms(bg, np.full(3, 76.6), 50.0, 38.3, trio=False)
    assert np.allclose(double["TDD link"], 2 * r["TDD link"])    # the dose reaches the TDD arm
    assert np.allclose(double["Profile link, midpoint"], r["Profile link, midpoint"])


def test_hours_since_meal():
    ts = np.arange(6) * 1800.0
    hrs, meals = hours_since_meal(ts, [0, 0, 30, 25, np.nan, 28])    # a 3 g rise is no entry
    assert meals == 1 and np.isnan(hrs[:2]).all()
    assert list(hrs[2:]) == [0.0, 0.5, 1.0, 1.5]
