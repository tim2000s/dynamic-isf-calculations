"""Windowed TDD reconstruction: treatments + basal profile → 5-min delivery grid →
the five TDD windows the v1/v2 ISF blend consumes
(TDD_4h, TDD_8to4h, TDD_24h, TDD_1d, TDD_7d-avg).

Pure functions over numpy/pandas — no I/O — so they're unit-testable and safe to call
from multiprocessing workers.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from inv008 import config


@dataclass
class TempBasalEvent:
    ts: float            # epoch seconds
    duration_min: float  # 0 → cancellation
    rate_u_h: float | None = None   # absolute U/h if known
    percent: float | None = None    # % of profile rate if absolute missing


def profile_rate_at(epoch_sec: np.ndarray, hourly_basal: np.ndarray) -> np.ndarray:
    """Profile basal rate (U/h) at each epoch second (UTC hour-of-day lookup)."""
    hours = (np.asarray(epoch_sec, dtype=np.int64) % 86400) // 3600
    return np.asarray(hourly_basal, dtype=float)[hours]


def build_delivery_grid(start_sec: int, end_sec: int,
                        bolus_ts: np.ndarray, bolus_units: np.ndarray,
                        temp_basals: list[TempBasalEvent],
                        hourly_basal: np.ndarray,
                        bin_sec: int = config.GRID_SEC) -> pd.DataFrame:
    """5-min grid of delivered insulin (units/bin) = basal (temp overlay on profile) + boluses.

    Temp basal semantics (Nightscout/oref): an event sets `rate_u_h` (absolute) or
    `percent` of profile from its timestamp for `duration_min`, superseded by any newer
    event; gaps run at profile rate. duration 0 / missing rate+percent = cancel.
    """
    start_sec = int(start_sec) // bin_sec * bin_sec
    end_sec = int(end_sec) // bin_sec * bin_sec
    bins = np.arange(start_sec, end_sec + bin_sec, bin_sec, dtype=np.int64)
    n = len(bins)

    # --- basal rate per bin: start from profile, overlay temp segments ---
    rate = profile_rate_at(bins, hourly_basal)
    if temp_basals:
        evs = sorted(temp_basals, key=lambda e: e.ts)
        for i, ev in enumerate(evs):
            if ev.duration_min <= 0 or (ev.rate_u_h is None and ev.percent is None):
                continue  # cancellation — profile rate already in place
            seg_start = ev.ts
            seg_end = ev.ts + ev.duration_min * 60.0
            if i + 1 < len(evs):
                seg_end = min(seg_end, evs[i + 1].ts)  # superseded by next event
            i0 = np.searchsorted(bins, seg_start, side="left")
            i1 = np.searchsorted(bins, seg_end, side="left")
            if i1 <= i0:
                continue
            if ev.rate_u_h is not None:
                rate[i0:i1] = ev.rate_u_h
            else:
                rate[i0:i1] = profile_rate_at(bins[i0:i1], hourly_basal) * ev.percent / 100.0

    basal_units = rate * bin_sec / 3600.0

    # --- boluses summed into bins ---
    bolus_per_bin = np.zeros(n)
    if len(bolus_ts):
        idx = np.searchsorted(bins, np.asarray(bolus_ts, dtype=np.int64), side="right") - 1
        ok = (idx >= 0) & (idx < n)
        np.add.at(bolus_per_bin, idx[ok], np.asarray(bolus_units, dtype=float)[ok])

    return pd.DataFrame({
        "ts": bins,
        "basal_u": basal_units,
        "bolus_u": bolus_per_bin,
        "total_u": basal_units + bolus_per_bin,
    })


def windowed_tdd(grid: pd.DataFrame, bin_sec: int = config.GRID_SEC,
                 min_days_for_7d: int = config.MIN_DAYS_FOR_7D) -> pd.DataFrame:
    """Per-bin TDD components mirroring the device TDD calculator:

    tdd_4h    = calculateDaily(-4, 0)   → trailing 4 h sum
    tdd_8to4h = calculateDaily(-8, -4)  → sum of the window 8h..4h ago
    tdd_24h   = calculateDaily(-24, 0)  → trailing 24 h sum
    tdd_1d    = averageTDD(calculate(1)) → yesterday's calendar-day total
    tdd_7d    = averageTDD(calculate(7, allowMissingDays=true))
               → mean of valid calendar days among the previous 7 (≥ min_days_for_7d)

    A calendar day is valid if it is complete on the grid and saw ≥ 1 bolus event
    (zero-bolus days indicate upload gaps and would silently deflate the average).
    """
    ts = grid["ts"].to_numpy()
    total = grid["total_u"].to_numpy()
    csum = np.concatenate([[0.0], np.cumsum(total)])

    def trailing(nbins: int, offset: int = 0) -> np.ndarray:
        """Sum of bins (i-offset-nbins, i-offset]."""
        i = np.arange(len(total))
        hi = np.clip(i + 1 - offset, 0, len(total))
        lo = np.clip(i + 1 - offset - nbins, 0, len(total))
        out = csum[hi] - csum[lo]
        # mark windows that extend before the grid start as missing
        out = np.where(i + 1 - offset - nbins < 0, np.nan, out)
        return out

    per_4h = int(4 * 3600 / bin_sec)
    per_24h = int(24 * 3600 / bin_sec)
    tdd_4h = trailing(per_4h)
    tdd_8to4h = trailing(per_4h, offset=per_4h)
    tdd_24h = trailing(per_24h)

    # calendar-day totals
    day = ts // 86400
    df = pd.DataFrame({"day": day, "total_u": total, "bolus_u": grid["bolus_u"].to_numpy()})
    daily = df.groupby("day").agg(total=("total_u", "sum"),
                                  n_bins=("total_u", "size"),
                                  boluses=("bolus_u", lambda s: (s > 0).sum()))
    full_day_bins = 86400 // bin_sec
    daily["valid"] = (daily["n_bins"] == full_day_bins) & (daily["boluses"] >= 1)

    day_index = daily.index.to_numpy()
    day_total = daily["total"].to_numpy()
    day_valid = daily["valid"].to_numpy()

    # yesterday's total / previous-7-day average, evaluated once per day
    tdd_1d_by_day = {}
    tdd_7d_by_day = {}
    for j, d in enumerate(day_index):
        prev1 = np.where(day_index == d - 1)[0]
        tdd_1d_by_day[d] = day_total[prev1[0]] if len(prev1) and day_valid[prev1[0]] else np.nan
        mask = (day_index >= d - 7) & (day_index <= d - 1) & day_valid
        vals = day_total[mask]
        tdd_7d_by_day[d] = float(np.mean(vals)) if len(vals) >= min_days_for_7d else np.nan

    tdd_1d = np.array([tdd_1d_by_day[d] for d in day])
    tdd_7d = np.array([tdd_7d_by_day[d] for d in day])

    return pd.DataFrame({
        "ts": ts,
        "tdd_4h": tdd_4h, "tdd_8to4h": tdd_8to4h, "tdd_24h": tdd_24h,
        "tdd_1d": tdd_1d, "tdd_7d": tdd_7d,
    })
