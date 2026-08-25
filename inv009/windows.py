"""Overnight fasting windows, and how much insulin acted in each.

One row per window. A window is four hours starting on the hour somewhere in the
night, and it carries the glucose it started and ended at, the insulin action
that happened inside it split by whether that insulin was already committed when
it opened, the dose the person was running, and every screening quantity needed
to decide later whether it was really fasting.

Screens are stored as columns rather than applied here. Which of them a cohort
can afford is a judgement that belongs with the analysis, not with the cache, and
rebuilding a million windows to relax one threshold is a waste of an afternoon.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import action, config, insulin_models as M


def _hours_since(flag: np.ndarray) -> np.ndarray:
    """Hours since the last True, per bin, infinite if there has not been one."""
    n = len(flag)
    idx = np.flatnonzero(flag)
    if len(idx) == 0:
        return np.full(n, np.inf)
    pos = np.searchsorted(idx, np.arange(n), side="right") - 1
    out = np.where(pos >= 0, (np.arange(n) - idx[np.clip(pos, 0, None)]) * config.GRID_MIN / 60.0,
                   np.inf)
    return out


def _centred_median(x: np.ndarray, bins: int) -> np.ndarray:
    return pd.Series(x).rolling(bins, center=True, min_periods=1).median().to_numpy()


def _forward_max(x: np.ndarray, span: int) -> np.ndarray:
    """Maximum of x over the next `span` bins, indexed by the opening bin."""
    s = pd.Series(x[::-1]).rolling(span, min_periods=1).max().to_numpy()[::-1]
    return s


def extract_windows(grid: pd.DataFrame, tddw: pd.DataFrame, subject_id: str, study: str,
                    age: float | None, subj: dict, models: tuple[str, ...]) -> pd.DataFrame:
    """Every candidate overnight window for one subject."""
    n = len(grid)
    h = int(config.HORIZON_MIN / config.GRID_MIN)
    cgm = grid.cgm.to_numpy(dtype=float)
    total_u = grid.total_u.to_numpy(dtype=float)
    bolus_u = grid.bolus_u.to_numpy(dtype=float)
    carbs = grid.carbs_g.to_numpy(dtype=float)

    med_bins = max(int(config.ENDPOINT_MEDIAN_MIN / config.GRID_MIN), 1)
    bg_s = _centred_median(cgm, med_bins)

    # Delivery as a deviation from this person's own routine. Overnight, total
    # delivery barely varies from night to night, so a regression on it is asking
    # a question the data cannot answer; the deviation is the part that moves.
    ref = action.reference_profile(grid.ts, total_u)
    u_net = action.net_delivery(total_u, ref)

    cols: dict[str, np.ndarray] = {}
    for m in models:
        kern = M.kernel(m)
        a_pre, a_in = action.window_action(total_u, kern, h)
        an_pre, an_in = action.window_action(u_net, kern, h)
        cols[f"a_pre_{m}"] = a_pre
        cols[f"a_in_{m}"] = a_in
        cols[f"a_net_pre_{m}"] = an_pre
        cols[f"a_net_in_{m}"] = an_in
        cols[f"iob0_{m}"] = action.iob_series(total_u, kern)

    # Candidate starts: on the hour, in the night hours.
    hour = grid.ts.dt.hour.to_numpy()
    minute = grid.ts.dt.minute.to_numpy()
    ok = np.isin(hour, config.START_HOURS) & (minute == 0)
    ok &= (np.arange(n) + h) < n
    starts = np.flatnonzero(ok)
    if len(starts) == 0:
        return pd.DataFrame()

    # Screening quantities.
    has_carbs = np.isfinite(carbs).any() and (carbs > 0).any()
    carb_free = _hours_since(carbs > 0) if has_carbs else np.full(n, np.nan)
    bolus_free = _hours_since(bolus_u > 0)

    rise30 = np.full(n, np.nan)
    step = max(int(30 / config.GRID_MIN), 1)
    rise30[:n - step] = cgm[step:] - cgm[:n - step]
    max_rise = _forward_max(np.nan_to_num(rise30, nan=-1e9), h)
    min_bg = -_forward_max(-np.nan_to_num(cgm, nan=1e9), h)
    n_cgm = pd.Series(np.isfinite(cgm).astype(float)[::-1]).rolling(h, min_periods=1)\
              .sum().to_numpy()[::-1]

    # Where glucose had been, for windows whose insulin was chosen by a
    # controller watching it. Under closed loop the delivery is close to a
    # function of recent glucose, so conditioning on recent glucose is the only
    # way any independent variation in insulin is left to read.
    lag_cols = {}
    for lag_min in (60, 120):
        lag = int(lag_min / config.GRID_MIN)
        v = np.full(n, np.nan)
        v[lag:] = bg_s[:n - lag]
        lag_cols[f"bg_m{lag_min}"] = v

    pre_slope = np.full(n, np.nan)
    pre = max(int(30 / config.GRID_MIN), 1)
    pre_slope[pre:] = (bg_s[pre:] - bg_s[:n - pre]) / 30.0     # mg/dL per minute

    # A meal-sized bolus, scaled to this person's daily dose, stands in for a
    # meal wherever carbohydrate is not recorded.
    tdd_u = subj.get("tdd_u", np.nan)
    meal_thr = (config.MEAL_BOLUS_FRAC_TDD * tdd_u) if np.isfinite(tdd_u) else np.inf
    meal_free = _hours_since(bolus_u >= meal_thr)

    csum_b = np.concatenate([[0.0], np.cumsum(bolus_u)])
    csum_meal = np.concatenate([[0.0], np.cumsum(np.where(bolus_u >= meal_thr, bolus_u, 0.0))])
    csum_c = np.concatenate([[0.0], np.cumsum(carbs)])

    def carbs_before(hours: float) -> np.ndarray:
        """Grams recorded in the `hours` before each bin.

        Needed to separate two explanations that make the same prediction. A
        person whose recent dose is high has usually been eating, so
        carbohydrate still absorbing at the start of a window would look like
        reduced sensitivity, and would track recent dose while doing it.
        Without this column the dose relationship and the carbohydrate tail
        cannot be told apart.
        """
        k = int(hours * 60 / config.GRID_MIN)
        lo = np.maximum(np.arange(n) - k, 0)
        return csum_c[np.arange(n)] - csum_c[lo]

    csum_bas = np.concatenate([[0.0], np.cumsum(grid.basal_u.to_numpy(dtype=float))])
    first30 = max(int(30 / config.GRID_MIN), 1)

    out = pd.DataFrame({
        "subject_id": subject_id,
        "study": study,
        "t0": grid.ts.to_numpy()[starts],
        "hour": hour[starts],
        "bg0": bg_s[starts],
        "bg_end": bg_s[starts + h],
        "pre_slope": pre_slope[starts],
        "bg_m60": lag_cols["bg_m60"][starts],
        "bg_m120": lag_cols["bg_m120"][starts],
        "n_cgm": n_cgm[starts],
        "min_bg": min_bg[starts],
        "max_rise30": max_rise[starts],
        "carb_free_h": carb_free[starts],
        "bolus_free_h": bolus_free[starts],
        # Measured one bin before the window opens, so that a bolus given at the
        # start of the window does not count as breaking its own quiet period.
        "quiet_before_h": bolus_free[np.maximum(starts - 1, 0)],
        "meal_free_h": meal_free[starts],
        "meal_thr_u": meal_thr,
        "bolus_in_u": csum_b[starts + h] - csum_b[starts],
        "bolus_first30_u": csum_b[starts + first30] - csum_b[starts],
        "meal_in_u": csum_meal[starts + h] - csum_meal[starts],
        "carbs_in_g": csum_c[starts + h] - csum_c[starts],
        "carbs_prev_4h": carbs_before(4)[starts],
        "carbs_prev_8h": carbs_before(8)[starts],
        "carbs_prev_12h": carbs_before(12)[starts],
        "carbs_prev_24h": carbs_before(24)[starts],
        "basal_in_u": csum_bas[starts + h] - csum_bas[starts],
    })
    out["drop"] = out.bg0 - out.bg_end
    for m in models:
        for pref in ("a_pre", "a_in", "a_net_pre", "a_net_in", "iob0"):
            out[f"{pref}_{m}"] = cols[f"{pref}_{m}"][starts]

    for c in ("tdd_4h", "tdd_8to4h", "tdd_1d", "tdd_7d", "tdd_blend"):
        out[c] = tddw[c].to_numpy()[starts]

    out["age"] = age
    out["closed_loop"] = config.COHORTS[study]["closed_loop"]
    out["has_carb_stream"] = has_carbs
    for k, v in subj.items():
        out[k] = v
    return out


def screen(w: pd.DataFrame, strict: bool = False) -> pd.Series:
    """The fasting screen: which windows can be read as a measure of insulin action.

    The primary screen conditions only on things known before the window opened,
    or on the recording itself. That is deliberate. Screening on what glucose
    went on to do would drop the windows where insulin worked least well and keep
    the ones where it worked best, which does not clean the estimate, it inflates
    it. `strict` adds those outcome-conditioned screens so the size of that effect
    can be measured rather than assumed.

    Where carbohydrate is logged a meal is excluded directly. Where it is not, a
    bolus large relative to the person's own daily dose stands in for one, which
    is weaker: it catches announced meals and misses unannounced ones. That gap
    is why the two cohorts that log carbohydrate carry the headline and the rest
    replicate.
    """
    h = config.HORIZON_MIN / config.GRID_MIN
    ok = (w.bg0.between(config.BG0_MIN, config.BG0_MAX)
          & (w.n_cgm >= config.MIN_CGM_FRACTION * h)
          & w.bg_end.notna()
          & w.tdd_blend.notna())
    carbed = w.has_carb_stream.fillna(False).astype(bool)
    ok &= pd.Series(np.where(carbed,
                             (w.carb_free_h.fillna(0) >= config.CARB_FREE_H)
                             & (w.carbs_in_g.fillna(1) <= 0),
                             (w.meal_free_h.fillna(0) >= config.MEAL_FREE_H)
                             & (w.meal_in_u.fillna(1) <= 0)),
                    index=w.index)
    if strict:
        ok &= (w.min_bg >= config.BG_FLOOR_IN_WINDOW) & (w.max_rise30 <= config.RISE_MAX_30)
    return ok.fillna(False)
