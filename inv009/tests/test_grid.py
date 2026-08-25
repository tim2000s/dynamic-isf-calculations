"""The grid has to conserve insulin. If it does not, every sensitivity is wrong."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from inv009 import grid as G


def _streams(cgm, basal=None, bolus=None, carbs=None, sched=None):
    empty = pd.DataFrame(columns=["ts_local"])
    return {"cgm": cgm,
            "basal": basal if basal is not None else empty.assign(rate_u_hr=[]),
            "bolus": bolus if bolus is not None else empty.assign(bolus_u=[], delivery_duration_s=[]),
            "carbs": carbs if carbs is not None else empty.assign(carbs_g=[]),
            "sched": sched if sched is not None else empty.assign(sched_rate_u_hr=[]),
            "wizard": empty}


def _cgm(days=20):
    ts = pd.date_range("2019-01-01", periods=days * 288, freq="5min")
    return pd.DataFrame({"ts_local": ts, "cgm_mgdl": 120.0})


def _terminated(basal, cgm, rate=0.0):
    """Close the pump record at the end of the glucose record.

    build_grid stops at the last insulin row on purpose, because a rate holding
    for ever past the end of an upload is a fiction. A test that wants a full
    length grid has to say when the record ends, exactly as real data does.
    """
    tail = pd.DataFrame({"ts_local": [cgm.ts_local.iloc[-1]], "rate_u_hr": [rate]})
    return pd.concat([basal, tail], ignore_index=True)


def test_epoch_seconds_survive_microsecond_resolution():
    """psycopg2 gives pandas microsecond datetimes, not nanosecond ones."""
    ts = pd.to_datetime(["2019-01-01 00:00:00", "2019-01-01 01:00:00"])
    for unit in ("ns", "us", "ms", "s"):
        s = G.epoch_s(ts.astype(f"datetime64[{unit}]"))
        assert s[1] - s[0] == pytest.approx(3600.0)


def test_constant_basal_integrates_to_rate_times_time():
    cgm = _cgm()
    basal = _terminated(pd.DataFrame({"ts_local": [cgm.ts_local.iloc[0] - pd.Timedelta(hours=1)],
                                      "rate_u_hr": [1.0]}), cgm, rate=1.0)
    g = G.build_grid(_streams(cgm, basal=basal))
    hours = (g.ts.iloc[-1] - g.ts.iloc[0]).total_seconds() / 3600 + 5 / 60
    assert g.basal_u.sum() == pytest.approx(hours * 1.0, rel=1e-9)
    assert g.basal_u.iloc[5] == pytest.approx(5 / 60, rel=1e-9)


def test_basal_step_changes_are_not_samples():
    """A rate holds until the next row. Averaging the rows instead would be wrong."""
    cgm = _cgm()
    t0 = cgm.ts_local.iloc[0]
    basal = _terminated(pd.DataFrame({
        "ts_local": [t0, t0 + pd.Timedelta(hours=1), t0 + pd.Timedelta(hours=2)],
        "rate_u_hr": [2.0, 0.0, 1.0]}), cgm, rate=1.0)
    g = G.build_grid(_streams(cgm, basal=basal))
    first_three_h = g[g.ts < t0 + pd.Timedelta(hours=3)].basal_u.sum()
    assert first_three_h == pytest.approx(2.0 * 1 + 0.0 * 1 + 1.0 * 1, rel=1e-9)


def test_basal_change_off_the_grid_boundary_is_exact():
    cgm = _cgm()
    t0 = cgm.ts_local.iloc[0]
    basal = _terminated(pd.DataFrame({"ts_local": [t0, t0 + pd.Timedelta(minutes=2)],
                                      "rate_u_hr": [6.0, 0.0]}), cgm)
    g = G.build_grid(_streams(cgm, basal=basal))
    assert g.basal_u.iloc[0] == pytest.approx(6.0 * 2 / 60, rel=1e-9)
    assert g.basal_u.iloc[1] == pytest.approx(0.0, abs=1e-12)


def test_standard_bolus_lands_in_one_bin():
    cgm = _cgm()
    t0 = cgm.ts_local.iloc[0]
    bolus = pd.DataFrame({"ts_local": [t0 + pd.Timedelta(minutes=7)],
                          "bolus_u": [4.0], "delivery_duration_s": [0.0]})
    basal = _terminated(pd.DataFrame({"ts_local": [t0], "rate_u_hr": [0.0]}), cgm)
    g = G.build_grid(_streams(cgm, basal=basal, bolus=bolus))
    assert g.bolus_u.sum() == pytest.approx(4.0)
    assert g.bolus_u.iloc[1] == pytest.approx(4.0)


def test_extended_bolus_conserves_mass_and_spreads():
    cgm = _cgm()
    t0 = cgm.ts_local.iloc[0]
    bolus = pd.DataFrame({"ts_local": [t0], "bolus_u": [6.0],
                          "delivery_duration_s": [3600.0]})
    basal = _terminated(pd.DataFrame({"ts_local": [t0], "rate_u_hr": [0.0]}), cgm)
    g = G.build_grid(_streams(cgm, basal=basal, bolus=bolus))
    assert g.bolus_u.sum() == pytest.approx(6.0, rel=1e-9)
    assert (g.bolus_u > 0).sum() == 12                  # one hour of five minute bins
    assert g.bolus_u.iloc[0] == pytest.approx(0.5, rel=1e-9)


def test_grid_is_bounded_by_the_insulin_record():
    """Glucose outside the pump record cannot be attributed to insulin."""
    cgm = _cgm(days=60)
    t0 = cgm.ts_local.iloc[0]
    basal = pd.DataFrame({"ts_local": [t0 + pd.Timedelta(days=10),
                                       t0 + pd.Timedelta(days=30)],
                          "rate_u_hr": [1.0, 1.0]})   # record ends at day 30
    g = G.build_grid(_streams(cgm, basal=basal))
    assert g.ts.iloc[0] >= t0 + pd.Timedelta(days=10)
    assert g.ts.iloc[-1] <= t0 + pd.Timedelta(days=30)


def test_short_records_are_rejected():
    assert G.build_grid(_streams(_cgm(days=3),
                                 basal=pd.DataFrame({"ts_local": [pd.Timestamp("2019-01-01")],
                                                     "rate_u_hr": [1.0]}))) is None


def test_millisecond_durations_are_repaired():
    """Loop, PEDAP and REPLACE-BG store extended bolus durations in milliseconds."""
    fixed = G.fix_duration_units(np.array([0.0, 3599.0, 3_600_000.0, 52_680_000.0, np.nan]))
    assert fixed[0] == 0.0
    assert fixed[1] == 3599.0             # DCLP-style seconds are left alone
    assert fixed[2] == 3600.0             # a one hour square wave
    assert fixed[3] == 0.0                # 14.6 h is beyond any pump, so treat as instant
    assert fixed[4] == 0.0


def test_a_millisecond_duration_does_not_smear_the_record():
    cgm = _cgm(days=20)
    t0 = cgm.ts_local.iloc[0]
    bolus = pd.DataFrame({"ts_local": [t0 + pd.Timedelta(days=1)], "bolus_u": [5.0],
                          "delivery_duration_s": [3_600_000.0]})
    basal = _terminated(pd.DataFrame({"ts_local": [t0], "rate_u_hr": [0.0]}), cgm)
    g = G.build_grid(_streams(cgm, basal=basal, bolus=bolus))
    assert g.bolus_u.sum() == pytest.approx(5.0, rel=1e-9)
    assert (g.bolus_u > 0).sum() == 12          # one hour, not a thousand
