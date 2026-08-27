"""Effective sensitivity from insulin action over T to T+6h, stepped every five minutes.

The question a sensitivity factor answers is how far glucose moves per unit of
insulin that acts. Both halves are computable at every reading in the record.

The glucose half is the change from T to T+6h. The insulin half is the action over
that same interval, which follows from the delivery record and the insulin action
curve without needing any part of the record to be quiet:

    action(T, T+6h) = IOB(T) - IOB(T+6h) + delivered(T, T+6h)

Everything on board at T that finishes acting inside the window, plus everything
delivered inside it that acts before the window closes, less whatever is carried
past the end. Boluses and temporary basals both move that quantity, and it is that
variation which identifies the effect. An earlier version of this module screened
out any interval where insulin arrived, which discarded the variation being
measured and left only the tail of the basal rate.

Delivery is taken net of the programmed basal throughout, because basal exists to
offset hepatic glucose output rather than to lower glucose, so counting it puts
insulin in the denominator that is not there to do the job being measured.

Insulin action uses Loop's own model for the Loop cohort, chosen per person against
the insulin on board their app recorded, and the oref exponential elsewhere.

Carbohydrate is the one thing that must be screened rather than modelled, since
glucose rising from a meal is subtracted from the fall. Nothing is allowed in the
six hours before T, so nothing is still absorbing, nor inside the window.

Two estimates are reported per person. The ratio is the median of the per-point
ratio, which is the direct reading of the question. The slope is a regression of
the glucose change on the action with an intercept, which absorbs whatever that
person's glucose does at that time of day irrespective of insulin, and is the
quantity a correction dose is decided with.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp

import numpy as np
import pandas as pd

from . import config, db, grid as gridmod, insulin_models as M

HORIZON_MIN = 360.0            # one insulin duration
MODEL_DEFAULT = "oref_6h75"
LOOP_FALLBACK = "loop_adult"
CARB_CLEAR_H = 6.0
MIN_ACTION_U = 0.50            # below this the ratio is noise over noise
BG_LO, BG_HI = 70.0, 400.0
MIN_POINTS = 100
MIN_ACTION_SD = 0.30           # a person whose action never varies identifies nothing
# Glucose at T decides what a fall can possibly be. From 85 mg/dL there is nowhere
# to go and counter-regulation pushes the other way, so points there return a
# negative sensitivity that is a floor effect rather than a measurement. Pooling
# across all levels averages those against the elevated points where a correction
# is actually given, which halves the estimate. Bands are reported instead, and
# the headline is taken from 200 mg/dL up, where corrections happen.
SGV_BANDS = [(70, 100), (100, 130), (130, 160), (160, 200), (200, 250), (250, 400)]
HEADLINE_BAND = (200, 400)

_LOOP_MODELS: dict[str, str] = {}


def _loop_model(subject_id: str) -> str:
    global _LOOP_MODELS
    if not _LOOP_MODELS:
        p = config.RESULTS / "inv009_loop_model_choice.parquet"
        if p.exists():
            d = pd.read_parquet(p)
            d = d[d.accepted] if "accepted" in d.columns else d
            _LOOP_MODELS = dict(zip(d.subject_id, d.model))
        else:
            _LOOP_MODELS = {"_": LOOP_FALLBACK}
    return _LOOP_MODELS.get(subject_id, LOOP_FALLBACK)


def schedule_units(g: pd.DataFrame, subject_id: str, streams: dict) -> tuple[np.ndarray, str]:
    """Insulin per bin the programme called for, from whichever source records it.

    Three sources in descending order of directness. Loop and REPLACE-BG record a
    scheduled rate against a timestamp, which arrives on the grid already. DCLP3,
    DCLP5 and PEDAP ship the 48 half-hourly programmed rates on the case report
    form, which is a profile by time of day rather than a series. The bionic
    pancreas has no programme at all by design, so the reference is that person's
    own median delivery at each half hour, which is what they typically receive
    rather than what anybody set.

    Getting this wrong is not subtle. With no schedule subtracted the denominator
    becomes total insulin rather than insulin above basal, and six hours of basal
    swamps whatever a correction contributed.
    """
    n = len(g)
    if "sched_u" in g.columns and g.sched_u.notna().any():
        return g.sched_u.fillna(0.0).to_numpy(float), "recorded"
    hh = (g.ts.dt.hour.to_numpy() * 2 + (g.ts.dt.minute.to_numpy() >= 30)).astype(int)
    prof = db.basal_schedule(subject_id)
    if prof is not None and np.isfinite(prof).any():
        rate = np.where(np.isfinite(prof[hh]), prof[hh], np.nanmedian(prof))
        return rate * (config.GRID_MIN / 60.0), "programmed profile"
    tot = g.total_u.to_numpy(float)
    med = pd.Series(tot).groupby(hh).median()
    return med.reindex(hh).to_numpy(dtype=float), "typical delivery"


def analyse(subject_id: str) -> dict | None:
    study = subject_id.split(":")[0]
    model = _loop_model(subject_id) if study == "Loop" else MODEL_DEFAULT

    streams = db.streams(subject_id)
    g = gridmod.build_grid(streams)
    if g is None or g.empty:
        return None
    n = len(g)
    h = int(HORIZON_MIN / config.GRID_MIN)
    if n < 2 * h:
        return None

    kern = M.kernel(model)
    total = g.total_u.to_numpy(float)
    bolus = g.bolus_u.to_numpy(float)
    sched, sched_src = schedule_units(g, subject_id, streams)
    net = total - np.nan_to_num(sched)

    # Insulin still to act at each bin, then action over the window from both the
    # standing load and anything delivered inside it.
    iob = np.convolve(net, kern)[:n]
    cnet = np.concatenate([[0.0], np.cumsum(net)])
    action = np.full(n, np.nan)
    action[:n - h] = iob[:n - h] - iob[h:] + (cnet[h:n] - cnet[:n - h])

    bg = g.cgm.to_numpy(float)
    fall = np.full(n, np.nan)
    fall[:n - h] = bg[:n - h] - bg[h:]

    # Carbohydrate clearance either side; a meal proxy where none is logged.
    carbs = g.carbs_g.to_numpy(float) if "carbs_g" in g.columns else np.zeros(n)
    if not np.isfinite(carbs).any() or np.nansum(carbs) <= 0:
        tdd_est = float(np.nansum(total) / max((n * config.GRID_MIN) / 1440.0, 1e-9))
        carbs = np.where(bolus >= config.MEAL_BOLUS_FRAC_TDD * tdd_est, 1.0, 0.0)
    ccum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(carbs))])
    cb = int(CARB_CLEAR_H * 60 / config.GRID_MIN)
    carb_before = np.full(n, np.nan)
    carb_before[cb:] = ccum[cb:n] - ccum[:n - cb]
    carb_in = np.full(n, np.nan)
    carb_in[:n - h] = ccum[h:n] - ccum[:n - h]

    ok = (np.isfinite(fall) & np.isfinite(action) & np.isfinite(bg)
          & (bg >= BG_LO) & (bg <= BG_HI)
          & np.isfinite(carb_before) & (carb_before <= 0)
          & np.isfinite(carb_in) & (carb_in <= 0))
    if ok.sum() < MIN_POINTS:
        return None

    a, y = action[ok], fall[ok]
    if np.std(a) < MIN_ACTION_SD:
        return None

    # Ratio: the direct reading, on points where enough insulin acted to divide by.
    big = a >= MIN_ACTION_U
    ratio = float(np.median(y[big] / a[big])) if big.sum() >= 30 else np.nan
    bgk = bg[ok]
    bands = {}
    for lo, hi in SGV_BANDS:
        m = big & (bgk >= lo) & (bgk < hi)
        bands[f"isf_{lo}"] = float(np.median(y[m] / a[m])) if m.sum() >= 30 else np.nan
        bands[f"n_{lo}"] = int(m.sum())
    hm = big & (bgk >= HEADLINE_BAND[0]) & (bgk < HEADLINE_BAND[1])
    head = float(np.median(y[hm] / a[hm])) if hm.sum() >= 30 else np.nan

    # Slope: a regression with an intercept, which absorbs whatever this person's
    # glucose does over six hours irrespective of insulin.
    ac = a - a.mean()
    slope = float((ac * (y - y.mean())).sum() / (ac ** 2).sum())

    tdd = float(np.nansum(total) / max((n * config.GRID_MIN) / 1440.0, 1e-9))
    return dict(subject_id=subject_id, study=study, model=model, sched_src=sched_src,
                n_points=int(ok.sum()), n_big=int(big.sum()), tdd_u=tdd,
                action_median=float(np.median(a)), action_sd=float(np.std(a)),
                fall_median=float(np.median(y)), bg_median=float(np.median(bg[ok])),
                isf_ratio=ratio, isf_slope=slope, isf_head=head, n_head=int(hm.sum()),
                k_ratio=ratio * tdd, k_slope=slope * tdd, k_head=head * tdd, **bands)


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
        subs += db.subjects(s).subject_id.tolist()
    if args.limit:
        subs = subs[:args.limit]
    print(f"{len(subs)} subjects, {args.workers} workers, action over "
          f"T to T+{HORIZON_MIN / 60:.0f}h, stepped every {config.GRID_MIN:.0f} min")

    with mp.Pool(args.workers, maxtasksperchild=8) as pool:
        rows = [r for r in pool.map(analyse, subs) if r]
    r = pd.DataFrame(rows)
    if r.empty:
        print("nothing survived the screen")
        return 1
    r.to_parquet(config.RESULTS / "inv009_forward_isf.parquet", index=False)

    print(f"\n{len(r)} people, {int(r.n_points.sum()):,} five-minute points\n")
    print(f"{'cohort':<12s}{'people':>7s}{'points':>12s}{'action U':>9s}"
          f"{'ISF all':>9s}{'K':>7s}{'ISF 200+':>10s}{'K 200+':>8s}{'pts 200+':>10s}")
    out = []
    for s, d in r.groupby("study"):
        row = dict(cohort=s, n_people=int(len(d)), n_points=int(d.n_points.sum()),
                   action=float(d.action_median.median()), fall=float(d.fall_median.median()),
                   isf_ratio=float(d.isf_ratio.median()), k_ratio=float(d.k_ratio.median()),
                   isf_slope=float(d.isf_slope.median()), k_slope=float(d.k_slope.median()))
        row["isf_head"] = float(d.isf_head.median())
        row["k_head"] = float(d.k_head.median())
        row["n_head"] = int(d.n_head.sum())
        row["bands"] = {f"{lo}": float(d[f"isf_{lo}"].median()) for lo, _ in SGV_BANDS}
        out.append(row)
        print(f"{s:<12s}{row['n_people']:7d}{row['n_points']:12,d}{row['action']:9.2f}"
              f"{row['isf_ratio']:9.1f}{row['k_ratio']:7.0f}"
              f"{row['isf_head']:10.1f}{row['k_head']:8.0f}{row['n_head']:10,d}")
    print(f"\n{'ALL':<12s}{len(r):7d}{int(r.n_points.sum()):12,d}"
          f"{r.action_median.median():9.2f}"
          f"{r.isf_ratio.median():9.1f}{r.k_ratio.median():7.0f}"
          f"{r.isf_head.median():10.1f}{r.k_head.median():8.0f}{int(r.n_head.sum()):10,d}")
    print("\nBy glucose at T (mg/dL per acting unit). Negative values below 130 are a")
    print("floor effect: from there glucose cannot fall far and counter-regulation lifts it.\n")
    print(f"{'cohort':<12s}" + "".join(f"{lo}-{hi}".rjust(10) for lo, hi in SGV_BANDS))
    for s_, d in r.groupby("study"):
        print(f"{s_:<12s}" + "".join(
            f"{d[f'isf_{lo}'].median():10.1f}" for lo, _ in SGV_BANDS))
    (config.RESULTS / "inv009_forward_isf.json").write_text(json.dumps(
        dict(horizon_min=HORIZON_MIN, step_min=config.GRID_MIN,
             model_default=MODEL_DEFAULT, min_action_u=MIN_ACTION_U,
             n_people=int(len(r)), n_points=int(r.n_points.sum()),
             by_cohort=out), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
