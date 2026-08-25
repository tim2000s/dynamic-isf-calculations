"""When somebody changes their sensitivity setting, does the deviation signal notice?

This tests the premise underneath any self-titrating design. A detector like
autosens reports a ratio away from 1 for two quite different reasons: the person
genuinely changed, or their sensitivity factor was wrong all along. Only the
second is something a slow loop should fold into the base, and nothing so far
establishes that the signal carries it.

Loop is the place to ask. It never adjusts sensitivity at runtime, so a setting
that is wrong stays wrong for as long as the person leaves it, and the error sits
in the data undisturbed. Its bolus calculator also records the sensitivity in
force at every dose, so changes are visible.

The test is directional and falsifiable. Someone who RAISES their sensitivity
factor is deciding they were being over-dosed. If the signal reports base error,
it should have been reading below 1 beforehand, and should move toward 1 after.
Someone who LOWERS it should show the opposite. If the readings before a change
are unrelated to its direction, the signal does not carry base error and a slow
loop cannot be built on it.

One confound cannot be removed. Loop applies retrospective correction, which
absorbs part of the discrepancy before it reaches glucose, and the archive does
not say who had the integral variant enabled. That damps whatever is here, so it
works against finding an effect rather than manufacturing one.

    python3 -m inv009.isf_change_test
"""
from __future__ import annotations

import json
import multiprocessing as mp
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps

from . import action, config, db, grid as gridmod, insulin_models as M, tdd as tddmod

MODEL = "oref_6h75"
BIN = config.GRID_MIN
# Defaults are the strict pass. Only 27 of the 196 Loop subjects with a recorded
# sensitivity ever move it by 15%, so a relaxed pass is run as well, and both are
# reported rather than whichever is kinder.
MIN_CHANGE = 0.15          # a setting move smaller than this is not a decision
WINDOW_DAYS = 21           # either side of the change
GAP_DAYS = 3               # skipped either side, so the transition is not scored
MIN_NIGHTS = 40            # usable deviation bins, per side
CARB_ABSORPTION_H = 5.0
AUTOSENS_MIN, AUTOSENS_MAX = 0.7, 1.2


def _entered_series(wiz: pd.DataFrame, ts: pd.Series) -> np.ndarray:
    """The sensitivity in force at each bin, carried forward from the last bolus."""
    w = wiz.dropna(subset=["isf_mgdl_per_u"]).sort_values("ts_local")
    w = w[(w.isf_mgdl_per_u > 5) & (w.isf_mgdl_per_u < 400)]
    if len(w) < 60:
        return None
    idx = np.searchsorted(w.ts_local.to_numpy(), ts.to_numpy(), side="right") - 1
    out = np.where(idx >= 0, w.isf_mgdl_per_u.to_numpy()[np.clip(idx, 0, None)], np.nan)
    return out


def _find_change(entered: pd.DataFrame, min_change: float, window_days: int,
                 gap_days: int):
    """One clear step in the setting, taken as the largest sustained move."""
    # Daily median first, so a single odd bolus entry cannot look like a decision.
    daily = entered.set_index("ts_local").isf_mgdl_per_u.resample("1D").median().dropna()
    if len(daily) < 2 * (window_days + gap_days):
        return None
    best = None
    lo_i, hi_i = window_days + gap_days, len(daily) - (window_days + gap_days)
    for i in range(lo_i, hi_i):
        before = daily.iloc[i - window_days - gap_days:i - gap_days].median()
        after = daily.iloc[i + gap_days:i + gap_days + window_days].median()
        if not (np.isfinite(before) and np.isfinite(after)) or before <= 0:
            continue
        rel = after / before - 1.0
        if abs(rel) < min_change:
            continue
        if best is None or abs(rel) > abs(best[3]):
            best = (daily.index[i], float(before), float(after), float(rel))
    return best


def _autosens(dev: np.ndarray, sens: float, max_daily_basal: float) -> float:
    d = dev[np.isfinite(dev)]
    if len(d) < MIN_NIGHTS:
        return np.nan
    basal_off = float(np.median(d)) * (60.0 / BIN) / sens
    return float(np.clip(1.0 + basal_off / max_daily_basal, AUTOSENS_MIN, AUTOSENS_MAX))


def _rolling_ratios(dev: np.ndarray, sens: float, max_daily_basal: float,
                    hours: float = 24.0) -> np.ndarray:
    """The ratio as it would have been read through the period, not once over it.

    People do not act on an average. They act on seeing the value stuck at a
    bound, so the quantity that matters is how much of the time it spent there,
    which needs the reading reconstructed rolling rather than pooled.
    """
    n_back = int(hours * 60 / BIN)
    s = pd.Series(dev)
    med = s.rolling(n_back, min_periods=MIN_NIGHTS).median().to_numpy()
    basal_off = med * (60.0 / BIN) / sens
    return np.clip(1.0 + basal_off / max_daily_basal, AUTOSENS_MIN, AUTOSENS_MAX)


def analyse(job):
    # Every setting a worker needs travels in the job tuple. macOS spawns rather
    # than forks, so a worker re-imports this module and never sees anything
    # main() assigned to a module global. A first version set them in main() and
    # the strict and relaxed passes returned identical numbers to every decimal,
    # which is the documented signature of a configuration that never arrived.
    subject_id, min_change, window_days, gap_days = job
    try:
        streams = db.streams(subject_id)
        wiz = streams["wizard"]
        if wiz.empty:
            return None
        change = _find_change(wiz, min_change, window_days, gap_days)
        if change is None:
            return None
        t_change, isf_before, isf_after, rel = change

        g = gridmod.build_grid(streams)
        if g is None or g.sched_u.isna().all():
            return None
        subj = tddmod.subject_level(g)
        n = len(g)
        cgm = g.cgm.to_numpy(float)
        total = g.total_u.to_numpy(float)
        sched = np.nan_to_num(g.sched_u.to_numpy(float))

        entered = _entered_series(wiz, g.ts)
        if entered is None:
            return None

        kern = M.kernel(MODEL)
        act_kern = -np.diff(np.append(kern, 0.0))
        acting = np.convolve(total - sched, act_kern)[:n]

        # The deviation is measured against the setting that was in force, which
        # is what makes it a statement about that setting.
        dev = np.full(n, np.nan)
        dev[:-1] = (cgm[1:] - cgm[:-1]) + acting[:-1] * entered[:-1]
        carbs = g.carbs_g.to_numpy(float)
        win = int(CARB_ABSORPTION_H * 60 / BIN)
        absorbing = pd.Series((carbs > 0).astype(float)).rolling(win, min_periods=1).max().to_numpy() > 0
        dev = np.where(absorbing, np.nan, dev)
        dev = np.where((cgm < 80) & (dev > 0), 0.0, dev)

        hourly = pd.Series(g.basal_u.to_numpy(float)).rolling(12, min_periods=12).sum()
        mdb = float(np.nanpercentile(hourly.to_numpy(), 99))
        if not np.isfinite(mdb) or mdb <= 0.05:
            return None

        ts = g.ts
        pre = (ts >= t_change - pd.Timedelta(days=window_days + gap_days)) & \
              (ts < t_change - pd.Timedelta(days=gap_days))
        post = (ts > t_change + pd.Timedelta(days=gap_days)) & \
               (ts <= t_change + pd.Timedelta(days=window_days + gap_days))
        pre, post = pre.to_numpy(), post.to_numpy()
        if pre.sum() < 200 or post.sum() < 200:
            return None

        r_before = _autosens(dev[pre], float(np.nanmedian(entered[pre])), mdb)
        r_after = _autosens(dev[post], float(np.nanmedian(entered[post])), mdb)
        if not (np.isfinite(r_before) and np.isfinite(r_after)):
            return None

        # How often the reading sat at a bound, which is what people report acting on.
        roll = _rolling_ratios(dev, float(np.nanmedian(entered[np.isfinite(entered)])), mdb)
        def sat(mask):
            v = roll[mask]
            v = v[np.isfinite(v)]
            if len(v) < MIN_NIGHTS:
                return np.nan, np.nan
            return (float(np.mean(v >= AUTOSENS_MAX - 1e-6)),
                    float(np.mean(v <= AUTOSENS_MIN + 1e-6)))
        ceil_b, floor_b = sat(pre)
        ceil_a, floor_a = sat(post)
        return dict(subject_id=subject_id, t_change=t_change,
                    isf_before=isf_before, isf_after=isf_after, rel_change=rel,
                    raised=bool(rel > 0), ratio_before=r_before, ratio_after=r_after,
                    moved_toward_one=abs(r_after - 1.0) - abs(r_before - 1.0),
                    ceiling_before=ceil_b, floor_before=floor_b,
                    ceiling_after=ceil_a, floor_after=floor_a,
                    tdd_u=subj["tdd_u"])
    except Exception:
        return None


def main() -> int:
    global MIN_CHANGE, WINDOW_DAYS, GAP_DAYS
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-change", type=float, default=MIN_CHANGE)
    ap.add_argument("--window-days", type=int, default=WINDOW_DAYS)
    ap.add_argument("--gap-days", type=int, default=GAP_DAYS)
    a = ap.parse_args()
    MIN_CHANGE, WINDOW_DAYS, GAP_DAYS = a.min_change, a.window_days, a.gap_days
    print(f"change >= {MIN_CHANGE:.0%}, {WINDOW_DAYS}-day windows, {GAP_DAYS}-day gap")
    config.ensure_dirs()
    subs = db.subjects("Loop").subject_id.tolist()
    jobs = [(s, MIN_CHANGE, WINDOW_DAYS, GAP_DAYS) for s in subs]
    print(f"{len(subs)} Loop subjects, {config.WORKERS} workers")
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        rows = [r for r in pool.imap_unordered(analyse, jobs, chunksize=4) if r]
    R = pd.DataFrame(rows)
    if R.empty:
        print("no usable setting changes found")
        return 0
    tag = f"{int(MIN_CHANGE*100)}pct_{WINDOW_DAYS}d"
    R.to_parquet(config.RESULTS / f"inv009_isf_changes_{tag}.parquet", index=False)

    up, down = R[R.raised], R[~R.raised]
    print(f"\n{len(R)} people made a clear change to their sensitivity factor")
    print(f"  raised it (deciding they were over-dosed): {len(up)}")
    print(f"  lowered it (deciding they were under-dosed): {len(down)}")

    print("\nTHE TEST  what the signal read BEFORE the change")
    print("  if it reports base error, raising should follow a reading below 1")
    print(f"  {'group':>34s} {'n':>5s} {'median reading':>15s}")
    print(f"  {'before raising the setting':>34s} {len(up):5d} {up.ratio_before.median():15.3f}")
    print(f"  {'before lowering the setting':>34s} {len(down):5d} {down.ratio_before.median():15.3f}")
    u = sps.mannwhitneyu(up.ratio_before, down.ratio_before, alternative="less")
    print(f"  one-sided Mann-Whitney, raise below lower: p = {u.pvalue:.3g}")
    rho = sps.spearmanr(R.rel_change, R.ratio_before)
    print(f"  Spearman(size of change, reading before) = {rho.statistic:+.3f}  p = {rho.pvalue:.3g}")

    print("\n  and did the reading move toward 1 afterwards?")
    for name, d in (("raised", up), ("lowered", down), ("all", R)):
        closer = float((d.moved_toward_one < 0).mean())
        print(f"    {name:>8s}  before {d.ratio_before.median():.3f}  after "
              f"{d.ratio_after.median():.3f}  closer to 1 for {100 * closer:.0f}%")
    print("\n  THE BETTER TEST  people act on the reading being STUCK at a bound.")
    print("  Pinned at the ceiling means the loop wants more insulin than the clamp")
    print("  allows, so the setting is too weak and should be lowered.")
    S = R.dropna(subset=["ceiling_before", "floor_before"])
    if len(S) > 10:
        su, sd = S[S.raised], S[~S.raised]
        print(f"  {'group':>34s} {'n':>5s} {'% at ceiling':>13s} {'% at floor':>11s}")
        print(f"  {'before raising the setting':>34s} {len(su):5d} "
              f"{100 * su.ceiling_before.median():13.1f} {100 * su.floor_before.median():11.1f}")
        print(f"  {'before lowering the setting':>34s} {len(sd):5d} "
              f"{100 * sd.ceiling_before.median():13.1f} {100 * sd.floor_before.median():11.1f}")
        cu = sps.mannwhitneyu(sd.ceiling_before, su.ceiling_before, alternative="greater")
        print(f"  one-sided, lowerers more ceiling-bound than raisers: p = {cu.pvalue:.3g}")
        rr = sps.spearmanr(S.rel_change, S.ceiling_before)
        print(f"  Spearman(size of change, time at ceiling) = {rr.statistic:+.3f} "
              f"p = {rr.pvalue:.3g}")
        res_sat = dict(n=int(len(S)),
                       ceiling_before_raised=float(su.ceiling_before.median()),
                       ceiling_before_lowered=float(sd.ceiling_before.median()),
                       floor_before_raised=float(su.floor_before.median()),
                       floor_before_lowered=float(sd.floor_before.median()),
                       mannwhitney_p=float(cu.pvalue),
                       spearman=float(rr.statistic), spearman_p=float(rr.pvalue))
    else:
        res_sat = None

    w = sps.wilcoxon(R.moved_toward_one) if len(R) > 10 else None
    if w is not None:
        print(f"    Wilcoxon on the move toward 1: p = {w.pvalue:.3g}")

    res = dict(n=int(len(R)), n_raised=int(len(up)), n_lowered=int(len(down)),
               ratio_before_raised=float(up.ratio_before.median()),
               ratio_before_lowered=float(down.ratio_before.median()),
               mannwhitney_p=float(u.pvalue),
               spearman_change_vs_before=float(rho.statistic),
               spearman_p=float(rho.pvalue),
               frac_closer_after=float((R.moved_toward_one < 0).mean()),
               wilcoxon_p=float(w.pvalue) if w is not None else None)
    res["saturation_test"] = res_sat
    res["min_change"], res["window_days"] = MIN_CHANGE, WINDOW_DAYS
    (config.RESULTS / f"inv009_isf_changes_{tag}.json").write_text(
        json.dumps(res, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
