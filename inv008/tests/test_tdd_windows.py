"""Unit tests for the delivery grid and windowed-TDD reconstruction."""
import numpy as np
import pytest

from inv008.tdd_windows import (TempBasalEvent, build_delivery_grid,
                                profile_rate_at, windowed_tdd)

FLAT_1UH = np.ones(24)  # 1 U/h all day → 24 U/day basal
DAY = 86400


def test_profile_rate_lookup():
    hourly = np.arange(24, dtype=float)
    # 02:30 UTC → hour 2
    assert profile_rate_at(np.array([2 * 3600 + 1800]), hourly)[0] == 2.0


def test_grid_flat_basal_no_events():
    g = build_delivery_grid(0, DAY, np.array([]), np.array([]), [], FLAT_1UH)
    # 5-min bins at 1 U/h → 1/12 U per bin; total ≈ 24 U over the day
    assert g["basal_u"].iloc[0] == pytest.approx(1 / 12)
    assert g["total_u"].iloc[:288].sum() == pytest.approx(24.0)


def test_grid_bolus_lands_in_bin():
    g = build_delivery_grid(0, 3600, np.array([610.0]), np.array([2.5]), [], FLAT_1UH)
    bin_idx = 610 // 300  # bin starting at 600 s
    assert g["bolus_u"].iloc[bin_idx] == pytest.approx(2.5)
    assert g["bolus_u"].sum() == pytest.approx(2.5)


def test_temp_basal_absolute_overrides_profile():
    tb = [TempBasalEvent(ts=600.0, duration_min=30.0, rate_u_h=3.0)]
    g = build_delivery_grid(0, 3600, np.array([]), np.array([]), tb, FLAT_1UH)
    # bins 600–2400 run at 3 U/h (0.25 U/bin); others at 1 U/h
    assert g.loc[g["ts"] == 600, "basal_u"].iloc[0] == pytest.approx(0.25)
    assert g.loc[g["ts"] == 2100, "basal_u"].iloc[0] == pytest.approx(0.25)
    assert g.loc[g["ts"] == 2400, "basal_u"].iloc[0] == pytest.approx(1 / 12)


def test_temp_basal_percent_and_supersede():
    tb = [TempBasalEvent(ts=0.0, duration_min=60.0, percent=50.0),     # 0.5 U/h
          TempBasalEvent(ts=1800.0, duration_min=60.0, rate_u_h=2.0)]  # supersedes at 30 min
    g = build_delivery_grid(0, 7200, np.array([]), np.array([]), tb, FLAT_1UH)
    assert g.loc[g["ts"] == 0, "basal_u"].iloc[0] == pytest.approx(0.5 / 12)
    assert g.loc[g["ts"] == 1800, "basal_u"].iloc[0] == pytest.approx(2.0 / 12)
    # after the second temp expires (1800+3600=5400) → profile
    assert g.loc[g["ts"] == 5400, "basal_u"].iloc[0] == pytest.approx(1 / 12)


def test_windowed_tdd_constant_delivery():
    # 10 full days of flat 1 U/h + one 6 U bolus per day at noon → 30 U/day
    n_days = 10
    bolus_ts = np.array([d * DAY + 12 * 3600 for d in range(n_days)], dtype=float)
    bolus_u = np.full(n_days, 6.0)
    g = build_delivery_grid(0, n_days * DAY - 300, bolus_ts, bolus_u, [], FLAT_1UH)
    w = windowed_tdd(g)
    # late on day 9: all windows defined
    row = w.iloc[-1]
    assert row["tdd_4h"] == pytest.approx(4.0, abs=0.1)        # 4 h of basal only (post-noon... no:
    # last bin is 23:55 day 9 → trailing 4h has basal only = 4 U
    assert row["tdd_8to4h"] == pytest.approx(4.0, abs=0.1)
    assert row["tdd_24h"] == pytest.approx(30.0, abs=0.2)      # one bolus + 24 U basal
    assert row["tdd_1d"] == pytest.approx(30.0, abs=0.01)      # yesterday complete
    assert row["tdd_7d"] == pytest.approx(30.0, abs=0.01)


def test_windowed_tdd_early_rows_nan():
    g = build_delivery_grid(0, 2 * DAY, np.array([3600.0]), np.array([5.0]), [], FLAT_1UH)
    w = windowed_tdd(g)
    assert np.isnan(w["tdd_24h"].iloc[10])   # < 24 h of history
    assert np.isnan(w["tdd_7d"].iloc[10])    # no 7-day history yet


def test_windowed_tdd_zero_bolus_days_invalid():
    # days with no boluses must not poison the 7d average
    n_days = 10
    bolus_ts = np.array([d * DAY + 12 * 3600 for d in range(n_days) if d != 5], dtype=float)
    bolus_u = np.full(len(bolus_ts), 6.0)
    g = build_delivery_grid(0, n_days * DAY - 300, bolus_ts, bolus_u, [], FLAT_1UH)
    w = windowed_tdd(g)
    # day 6: yesterday (day 5) had no boluses → tdd_1d is NaN there
    day6 = w[(w["ts"] >= 6 * DAY) & (w["ts"] < 7 * DAY)]
    assert np.isnan(day6["tdd_1d"]).all()
    # but 7d average still defined (other days valid) and equals 30
    assert day6["tdd_7d"].iloc[-1] == pytest.approx(30.0, abs=0.01)
