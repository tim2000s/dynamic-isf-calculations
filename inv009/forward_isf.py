"""Effective sensitivity from insulin on board at T against glucose change to T+6h.

Stepped every five minutes through every person's record, which is the direct form
of the question a sensitivity factor answers: this much insulin is on board now,
how far does glucose fall before it has finished acting.

    ISF_eff(T) = [ BG(T) - BG(T+6h) ] / IOB(T)

Three things have to hold for that ratio to mean anything, and each is a screen
rather than a correction.

Carbohydrate must be absent, because glucose rising from a meal is subtracted from
the fall and biases the ratio down. No carbohydrate is allowed in the six hours
before T, so nothing is still absorbing, nor in the six hours after.

No further insulin may arrive during the window beyond the programmed basal, or
the denominator understates what acted. Boluses inside the window are excluded
outright, and the net basal deviation across the window is required to be small.

The insulin on board must be net of the programmed basal. Basal exists to offset
hepatic glucose output, so counting it puts insulin in the denominator that is not
there to lower glucose, and the ratio collapses toward zero. Bolus-only insulin on
board is computed alongside, because that is the convention a pump's own display
uses and it makes the two comparable.

The remaining term is the one that cannot be screened away: over six hours glucose
also moves because basal is not exactly right. That is why the estimate is reported
both raw and as a departure from what this person usually does at this time of day,
the second removing any part of the movement that is routine for them at that hour.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp

import numpy as np
import pandas as pd

from . import config, db, grid as gridmod, insulin_models as M

HORIZON_MIN = 360.0            # one insulin duration
STEP_MIN = 5.0                 # every reading
MODEL = "oref_6h75"
CARB_CLEAR_H = 6.0             # before T and inside the window
MIN_IOB_U = 0.30               # below this the ratio is noise over noise
MAX_NET_IN_U = 0.20            # net insulin arriving during the window, above schedule
BG_LO, BG_HI = 70.0, 400.0
MIN_POINTS = 50


def iob_series(u: np.ndarray, kern: np.ndarray) -> np.ndarray:
    """Insulin still to act at each bin, being each dose times what remains of it."""
    return np.convolve(u, kern)[:len(u)]


def _typical_by_halfhour(v: np.ndarray, hh: np.ndarray) -> np.ndarray:
    """What this person usually does at this half hour, as a median."""
    out = np.full(len(v), np.nan)
    ok = np.isfinite(v)
    if not ok.any():
        return out
    s = pd.Series(v[ok]).groupby(hh[ok]).median()
    out[ok] = s.reindex(hh[ok]).to_numpy()
    return out


def analyse(subject_id: str) -> dict | None:
    streams = db.streams(subject_id)
    g = gridmod.build_grid(streams)
    if g is None or g.empty:
        return None
    n = len(g)
    h = int(HORIZON_MIN / config.GRID_MIN)
    if n < 2 * h:
        return None

    kern = M.kernel(MODEL)
    sched = g.sched_u.fillna(0.0).to_numpy(float)
    total = g.total_u.to_numpy(float)
    bolus = g.bolus_u.to_numpy(float)
    net = total - sched

    iob_net = iob_series(net, kern)
    iob_bolus = iob_series(bolus, kern)

    bg = g.cgm.to_numpy(float)
    fwd = np.full(n, np.nan)
    fwd[:n - h] = bg[:n - h] - bg[h:]          # fall from T to T+6h, positive is a fall

    # Insulin arriving during the window, above the programmed basal.
    cnet = np.concatenate([[0.0], np.cumsum(net)])
    net_in = np.full(n, np.nan)
    net_in[:n - h] = cnet[h:n] - cnet[:n - h]
    cbol = np.concatenate([[0.0], np.cumsum(bolus)])
    bol_in = np.full(n, np.nan)
    bol_in[:n - h] = cbol[h:n] - cbol[:n - h]

    # Carbohydrate clearance either side. Studies with no carbohydrate stream get
    # a meal proxy from bolus size, as everywhere else in this package.
    carbs = g.carbs_g.to_numpy(float) if "carbs_g" in g.columns else np.zeros(n)
    if not np.isfinite(carbs).any() or carbs.sum() <= 0:
        tdd = float(np.nansum(total) / max((n * config.GRID_MIN) / 1440.0, 1e-9))
        carbs = np.where(bolus >= config.MEAL_BOLUS_FRAC_TDD * tdd, 1.0, 0.0)
    ccum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(carbs))])
    carb_before = np.full(n, np.nan)
    cb = int(CARB_CLEAR_H * 60 / config.GRID_MIN)
    carb_before[cb:] = ccum[cb:n] - ccum[:n - cb]
    carb_in = np.full(n, np.nan)
    carb_in[:n - h] = ccum[h:n] - ccum[:n - h]

    keep = (np.isfinite(fwd) & np.isfinite(bg) & (bg >= BG_LO) & (bg <= BG_HI)
            & np.isfinite(carb_before) & (carb_before <= 0)
            & np.isfinite(carb_in) & (carb_in <= 0)
            & (bol_in <= 0) & (np.abs(net_in) <= MAX_NET_IN_U)
            & (iob_net >= MIN_IOB_U))
    if keep.sum() < MIN_POINTS:
        return None

    hh = (g.ts.dt.hour.to_numpy() * 2 + (g.ts.dt.minute.to_numpy() >= 30)).astype(int)
    typ_bg = _typical_by_halfhour(np.where(np.isfinite(fwd), fwd, np.nan), hh)
    typ_iob = _typical_by_halfhour(np.where(np.isfinite(iob_net), iob_net, np.nan), hh)

    k = keep
    raw = fwd[k] / iob_net[k]
    d_bg, d_iob = fwd[k] - typ_bg[k], iob_net[k] - typ_iob[k]
    ok = np.isfinite(d_bg) & np.isfinite(d_iob) & (np.abs(d_iob) >= 0.2)
    dep = (d_bg[ok] / d_iob[ok]) if ok.sum() >= 20 else np.array([np.nan])

    tdd = float(np.nansum(total) / max((n * config.GRID_MIN) / 1440.0, 1e-9))
    return dict(
        subject_id=subject_id, study=subject_id.split(":")[0],
        n_points=int(k.sum()), n_dep=int(ok.sum()),
        tdd_u=tdd,
        iob_net_median=float(np.median(iob_net[k])),
        iob_bolus_median=float(np.median(iob_bolus[k])),
        bg_median=float(np.median(bg[k])),
        fall_median=float(np.median(fwd[k])),
        isf_raw=float(np.median(raw)),
        isf_dep=float(np.nanmedian(dep)),
        k_raw=float(np.median(raw)) * tdd,
        k_dep=float(np.nanmedian(dep)) * tdd,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=config.WORKERS)
    args = ap.parse_args()
    config.ensure_dirs()

    studies = [args.study] if args.study else list(config.COHORTS)
    subs = []
    for s in studies:
        t = db.subjects(s)
        subs += t.subject_id.tolist()
    if args.limit:
        subs = subs[:args.limit]
    print(f"{len(subs)} subjects, {args.workers} workers, "
          f"stepping every {STEP_MIN:.0f} min to T+{HORIZON_MIN / 60:.0f}h")

    with mp.Pool(args.workers, maxtasksperchild=8) as pool:
        rows = [r for r in pool.map(analyse, subs) if r]
    r = pd.DataFrame(rows)
    if r.empty:
        print("nothing survived the screen")
        return 1
    r.to_parquet(config.RESULTS / "inv009_forward_isf.parquet", index=False)

    print(f"\n{len(r)} people, {int(r.n_points.sum()):,} five-minute points\n")
    print(f"{'cohort':<12s}{'people':>7s}{'points':>12s}{'IOB U':>7s}{'fall':>7s}"
          f"{'ISF raw':>9s}{'K raw':>7s}{'ISF dep':>9s}{'K dep':>7s}")
    out = []
    for s, d in r.groupby("study"):
        row = dict(cohort=s, n_people=int(len(d)), n_points=int(d.n_points.sum()),
                   iob=float(d.iob_net_median.median()),
                   fall=float(d.fall_median.median()),
                   isf_raw=float(d.isf_raw.median()), k_raw=float(d.k_raw.median()),
                   isf_dep=float(d.isf_dep.median()), k_dep=float(d.k_dep.median()))
        out.append(row)
        print(f"{s:<12s}{row['n_people']:7d}{row['n_points']:12,d}{row['iob']:7.2f}"
              f"{row['fall']:7.1f}{row['isf_raw']:9.1f}{row['k_raw']:7.0f}"
              f"{row['isf_dep']:9.1f}{row['k_dep']:7.0f}")
    print(f"\n{'ALL':<12s}{len(r):7d}{int(r.n_points.sum()):12,d}"
          f"{r.iob_net_median.median():7.2f}{r.fall_median.median():7.1f}"
          f"{r.isf_raw.median():9.1f}{r.k_raw.median():7.0f}"
          f"{r.isf_dep.median():9.1f}{r.k_dep.median():7.0f}")
    (config.RESULTS / "inv009_forward_isf.json").write_text(json.dumps(
        dict(horizon_min=HORIZON_MIN, step_min=STEP_MIN, model=MODEL,
             min_iob_u=MIN_IOB_U, max_net_in_u=MAX_NET_IN_U,
             n_people=int(len(r)), n_points=int(r.n_points.sum()),
             by_cohort=out), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
