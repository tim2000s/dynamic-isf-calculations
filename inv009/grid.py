"""One subject's record on a five minute grid.

Five minutes because that is what the CGM reports on. The awkward part is basal:
it arrives as step changes, a rate that holds until the next row rather than a
sample, and treating it as a sample is the standard way to get insulin totals
wrong by a large factor. The same is true of the recovered scheduled basal.

Both are integrated exactly rather than approximately. Delivered units are the
area under the rate, so the cumulative area is piecewise linear with breakpoints
where the rate changed, and reading it at the bin edges gives the units in each
bin with no rounding to grid boundaries at all.

Basal and boluses are kept in separate columns. Downstream needs to know which
insulin a controller chose in response to glucose, and a merged column cannot say.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config

MIN_SPAN_DAYS = 10


def epoch_s(ts) -> np.ndarray:
    """Epoch seconds, whatever datetime resolution the source happens to carry.

    psycopg2 hands back Python datetimes and pandas gives them microsecond
    resolution, not nanosecond, so casting the raw integers and dividing by 1e9
    is off by a thousand. Everything that measures a duration goes through here.
    """
    return np.asarray(ts, dtype="datetime64[ns]").astype("int64") / 1e9


def _step_units(ts_edges: np.ndarray, ev_ts: np.ndarray, rate: np.ndarray,
                end: float) -> np.ndarray:
    """Units delivered in each bin by a rate that holds until the next change.

    Exact: the cumulative delivered volume is piecewise linear in time with a
    breakpoint at every rate change, so interpolating it at the bin edges and
    differencing gives each bin's units however the changes fall against the grid.
    """
    if len(ev_ts) == 0:
        return np.zeros(len(ts_edges) - 1)
    bounds = np.append(ev_ts, end)
    dur_h = np.diff(bounds) / 3600.0
    cum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(rate) * dur_h)])
    u_at_edge = np.interp(ts_edges, bounds, cum)
    return np.diff(u_at_edge)


# No pump delivers a square wave longer than this, so anything beyond it is bad
# data rather than a very long bolus.
MAX_EXTENDED_S = 8 * 3600.0


def fix_duration_units(dur_s: np.ndarray) -> np.ndarray:
    """Repair extended-bolus durations that reached the database in milliseconds.

    studies.bolus.delivery_duration_s is in seconds for DCLP3 and DCLP5 and in
    milliseconds for Loop, PEDAP and REPLACE-BG, where the raw archives store
    Tidepool-style millisecond durations and the loader read them as seconds.
    A one hour square wave therefore arrives as 3,600,000 seconds, which is
    1,000 hours, and spreading a bolus over that smears insulin across most of
    the record: for one REPLACE-BG subject it put a bolus in 41,591 of 46,176
    bins.

    Nothing is guessed. A genuine duration cannot exceed a day, so a value that
    does is milliseconds and is divided by a thousand. Whatever still exceeds
    the longest square wave a pump can deliver is treated as no duration at all,
    which places the dose at its timestamp.
    """
    d = np.asarray(dur_s, dtype=float).copy()
    d = np.where(np.isfinite(d), d, 0.0)
    d = np.where(d > 86400.0, d / 1000.0, d)
    d = np.where(d > MAX_EXTENDED_S, 0.0, d)
    return np.clip(d, 0.0, None)


def _bolus_units(ts_edges: np.ndarray, ev_ts: np.ndarray, units: np.ndarray,
                 dur_s: np.ndarray) -> np.ndarray:
    """Units in each bin, with extended boluses spread over their delivery."""
    n = len(ts_edges) - 1
    out = np.zeros(n)
    dur_s = fix_duration_units(dur_s)
    plain = dur_s <= 0
    if plain.any():
        idx = np.searchsorted(ts_edges, ev_ts[plain], side="right") - 1
        ok = (idx >= 0) & (idx < n)
        np.add.at(out, idx[ok], units[plain][ok])
    # Extended boluses are a few per cent of records, so a loop over them costs
    # nothing and keeps the mass conservation obvious.
    for t, u, d in zip(ev_ts[~plain], units[~plain], dur_s[~plain]):
        lo = np.clip(np.interp(t, ts_edges, np.arange(n + 1)), 0, n)
        hi = np.clip(np.interp(t + d, ts_edges, np.arange(n + 1)), 0, n)
        if hi <= lo:
            i = int(np.clip(lo, 0, n - 1))
            out[i] += u
            continue
        first, last = int(np.floor(lo)), int(np.ceil(hi))
        edges = np.clip(np.arange(first, last + 1, dtype=float), lo, hi)
        frac = np.diff(edges) / (hi - lo)
        out[first:last] += u * frac
    return out


def build_grid(streams: dict[str, pd.DataFrame],
               min_span_days: int = MIN_SPAN_DAYS) -> pd.DataFrame | None:
    """Glucose, insulin and carbohydrate on a five minute grid, or None if unusable.

    Bounded by where insulin records exist: glucose often runs past the end of
    the pump record, and a window with no idea what insulin was given measures
    nothing.
    """
    cgm, basal, bolus = streams["cgm"], streams["basal"], streams["bolus"]
    if cgm.empty or (basal.empty and bolus.empty):
        return None

    ins = [d.ts_local for d in (basal, bolus) if not d.empty]
    lo = max(cgm.ts_local.min(), min(d.min() for d in ins))
    hi = min(cgm.ts_local.max(), max(d.max() for d in ins))
    if (hi - lo) < pd.Timedelta(days=min_span_days):
        return None

    grid = pd.date_range(lo.ceil("5min"), hi.floor("5min"), freq="5min")
    n = len(grid)
    if n < 12 * 24 * min_span_days:
        return None
    edges_dt = grid.union([grid[-1] + pd.Timedelta(minutes=5)])
    edges = epoch_s(edges_dt)
    end = float(edges[-1])

    out = pd.DataFrame({"ts": grid})
    g = cgm.set_index("ts_local").cgm_mgdl
    g = g[~g.index.duplicated(keep="first")]
    out["cgm"] = g.reindex(grid, method="nearest",
                           tolerance=pd.Timedelta("5min")).to_numpy(dtype=float)

    def secs(df):
        return epoch_s(df.ts_local)

    out["basal_u"] = _step_units(edges, secs(basal), basal.rate_u_hr.to_numpy(), end) \
        if not basal.empty else 0.0
    out["bolus_u"] = _bolus_units(edges, secs(bolus), bolus.bolus_u.fillna(0).to_numpy(),
                                  bolus.delivery_duration_s.to_numpy()) \
        if not bolus.empty else 0.0
    out["sched_u"] = _step_units(edges, secs(streams["sched"]),
                                 streams["sched"].sched_rate_u_hr.to_numpy(), end) \
        if not streams["sched"].empty else np.nan

    carbs = np.zeros(n)
    if not streams["carbs"].empty:
        c = streams["carbs"]
        idx = np.searchsorted(edges, secs(c), side="right") - 1
        ok = (idx >= 0) & (idx < n)
        np.add.at(carbs, idx[ok], c.carbs_g.fillna(0).to_numpy()[ok])
    out["carbs_g"] = carbs
    out["total_u"] = out.basal_u + out.bolus_u
    return out


def daily_totals(grid: pd.DataFrame) -> pd.DataFrame:
    """Per calendar day insulin, and whether the day looks complete enough to trust.

    A day missing bins is a day the pump was not uploading, and counting its
    partial total as a daily dose would drag a subject's average down.
    """
    day = grid.ts.dt.normalize()
    agg = grid.groupby(day).agg(total_u=("total_u", "sum"), basal_u=("basal_u", "sum"),
                                bolus_u=("bolus_u", "sum"), n=("ts", "size"),
                                n_cgm=("cgm", "count"))
    agg["complete"] = (agg.n == 288) & (agg.total_u > 0)
    return agg
