import numpy as np

from inv016.replay import ARMS, RBG_COEF, RBG_MIN, arms, hours_since_meal, rbg_shape


def test_arms_at_known_points():
    bg = np.array([260.0, 180.0, 100.0])
    r = arms(bg, np.full(3, 38.3), 50.0, 38.3, trio=False)      # 50 x 38.3 = 1915
    assert set(r) == set(ARMS)
    assert abs(r["Profile link, midpoint"][0] - 1) < 1e-9       # midpoint of 260 and 100 is 180
    assert abs(r["Profile link, start"][1] - 1) < 1e-9
    flat = r["Profile link, start, flat above 200"]
    assert flat[1] == r["Profile link, start"][1] and 1.14 < flat[0] < 1.15 < r["Profile link, start"][0]
    assert np.allclose(r["TDD link"], r["Profile link, midpoint"])   # the links agree at 1915
    assert np.allclose(r["Average-dose link"], r["TDD link"])
    assert r["Profile link, midpoint"][2] < 0.5 < r["Logarithmic as shipped"][2]
    double = arms(bg, np.full(3, 76.6), 50.0, 38.3, trio=False)
    assert np.allclose(double["TDD link"], 2 * r["TDD link"])    # the dose reaches the TDD arm
    assert np.allclose(double["Profile link, midpoint"], r["Profile link, midpoint"])


def test_rbg_shape():
    g = np.arange(100.0, 371.0)
    assert g[np.argmin(np.polyval(RBG_COEF, g))] == RBG_MIN
    # the refit's shape in results/inv013_refit.md: 2.17, 1.28, 0.89, 0.81 and 0.90
    assert np.allclose(rbg_shape(np.array([100.0, 150, 200, 250, 300]), 370.0), [2.17, 1.28, 0.89, 0.81, 0.90], atol=0.01)
    assert rbg_shape(90.0) == rbg_shape(100.0) and rbg_shape(350.0) == rbg_shape(RBG_MIN) < rbg_shape(350.0, 370.0)


def test_hours_since_meal():
    ts = np.arange(6) * 1800.0
    hrs, meals = hours_since_meal(ts, [0, 0, 30, 25, np.nan, 28])    # a 3 g rise is no entry
    assert meals == 1 and np.isnan(hrs[:2]).all()
    assert list(hrs[2:]) == [0.0, 0.5, 1.0, 1.5]
