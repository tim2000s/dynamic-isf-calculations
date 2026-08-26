"""Basal share and carb ratio, measured from behaviour rather than from settings.

The sensitivity constant is not recoverable from this data. INV-009 established
that: the shape of a person's sensitivity survives, the level does not, and the
Walsh sensitivity rule is a statement about the level. Repeating that measurement
adds nothing.

The other two Walsh rules are different, and neither has been measured here.

Basal share needs no model at all. It is delivered basal over delivered total,
straight from the pump record, so it carries no attenuation and covers the two
cohorts with no settings file. That is 446 people every earlier table left out.

Carb ratio can be measured without a sensitivity term. For a meal that was
announced, dosed for, and left glucose where it started, the ratio the person
needed was simply the grams divided by the units. No insulin model, no basal
assumption, no entered setting. Meals that did not end neutral are excluded
rather than corrected, which costs data and keeps the estimator honest.

    python3 -m inv009.measured_basal_cr
"""
from __future__ import annotations

import json
import multiprocessing as mp
import sys

import numpy as np
import pandas as pd

from . import config, data, db, grid as gridmod, stats

WALSH_BASAL, WALSH_CR = 0.50, 500.0
MEAL_MIN_G = 20.0            # smaller entries are corrections or top-ups
BOLUS_WINDOW_MIN = 20.0      # how close the dose must sit to the meal
QUIET_BEFORE_H = 4.0         # no other carbohydrate before
HORIZON_H = 5.0              # how long to let the meal play out
NEUTRAL_MGDL = 20.0          # what counts as returning to where it started
MIN_MEALS = 8


def basal_share() -> pd.DataFrame:
    """Delivered basal as a share of delivered total, per person, every cohort."""
    rows = []
    for f in sorted(config.WINDOW_CACHE.glob("*.parquet")):
        d = pd.read_parquet(f, columns=["subject_id", "study", "tdd_u", "tdd_basal_u",
                                        "basal_frac", "n_days", "age"])
        if d.empty:
            continue
        r = d.iloc[0]
        if not np.isfinite(r.tdd_u) or r.n_days < config.MIN_DAYS_TDD:
            continue
        rows.append(dict(subject_id=r.subject_id, study=r.study, tdd=float(r.tdd_u),
                         basal=float(r.tdd_basal_u), basal_frac=float(r.basal_frac),
                         age=float(r.age) if pd.notna(r.age) else np.nan,
                         n_days=int(r.n_days)))
    R = pd.DataFrame(rows)
    R["age_band"] = stats.band_of(R.age, config.AGE_BANDS)
    return R[R.basal_frac.between(0.05, 0.95)]


def meals_for(job):
    """Announced meals that ended where they started, and the ratio they implied."""
    subject_id, study = job
    try:
        streams = db.streams(subject_id)
        if streams["carbs"].empty or streams["bolus"].empty:
            return None
        g = gridmod.build_grid(streams)
        if g is None:
            return None
        n = len(g)
        cgm = g.cgm.to_numpy(float)
        carbs = g.carbs_g.to_numpy(float)
        bolus = g.bolus_u.to_numpy(float)
        bw = int(BOLUS_WINDOW_MIN / config.GRID_MIN)
        qb = int(QUIET_BEFORE_H * 60 / config.GRID_MIN)
        hz = int(HORIZON_H * 60 / config.GRID_MIN)
        csum_c = np.concatenate([[0.0], np.cumsum(carbs)])
        csum_b = np.concatenate([[0.0], np.cumsum(bolus)])

        out = []
        for i in np.flatnonzero(carbs >= MEAL_MIN_G):
            if i - qb < 0 or i + hz >= n:
                continue
            # Nothing else eaten before it, and nothing else during it.
            if csum_c[i] - csum_c[i - qb] > 0:
                continue
            if csum_c[i + hz + 1] - csum_c[i + 1] > 0:
                continue
            grams = float(carbs[i])
            units = float(csum_b[min(i + bw, n)] - csum_b[max(i - bw, 0)])
            if units <= 0.2:
                continue
            # No further insulin during the meal beyond the dose itself.
            extra = float(csum_b[i + hz + 1] - csum_b[min(i + bw, n)])
            if extra > 0.1 * units:
                continue
            start, end = cgm[i], cgm[i + hz]
            if not (np.isfinite(start) and np.isfinite(end)):
                continue
            excursion = end - start
            out.append(dict(subject_id=subject_id, study=study, grams=grams,
                            units=units, start_bg=start, excursion=excursion,
                            neutral=bool(abs(excursion) <= NEUTRAL_MGDL),
                            cr_implied=grams / units))
        return pd.DataFrame(out) if out else None
    except Exception:
        return None


def boot_ci(x, n_boot=3000, seed=0):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    v = np.array([np.median(x[rng.integers(0, len(x), len(x))]) for _ in range(n_boot)])
    return tuple(float(a) for a in np.percentile(v, [2.5, 97.5]))


def main() -> int:
    config.ensure_dirs()
    res = {}

    print("BASAL SHARE, delivered, needing no model and no settings file\n")
    B = basal_share()
    B.to_parquet(config.RESULTS / "inv009_measured_basal.parquet", index=False)
    print(f"  {'cohort':>28s} {'n':>5s} {'share':>7s} {'95% interval':>16s} {'covers 0.50':>12s}")
    res["basal"] = []
    for label, sub in list(B.groupby("study")) + [("EVERYONE", B)]:
        v = sub.basal_frac.dropna()
        if len(v) < 12:
            continue
        lo, hi = boot_ci(v)
        name = config.COHORTS[label]["label"] if label in config.COHORTS else label
        row = dict(cohort=name, n=int(len(v)), median=float(v.median()),
                   ci=[lo, hi], covers_walsh=bool(lo <= WALSH_BASAL <= hi))
        res["basal"].append(row)
        print(f"  {name:>28s} {row['n']:5d} {row['median']:7.2f} [{lo:.2f}, {hi:.2f}]"
              f" {'yes' if row['covers_walsh'] else 'no':>12s}")

    print("\nCARB RATIO, from meals that ended where they started\n")
    jobs = [(r.subject_id, s) for s in ("Loop", "ReplaceBG")
            for r in db.subjects(s).itertuples()]
    frames = []
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        for out in pool.imap_unordered(meals_for, jobs, chunksize=4):
            if out is not None:
                frames.append(out)
    M = pd.concat(frames, ignore_index=True)
    M.to_parquet(config.RESULTS / "inv009_measured_meals.parquet", index=False)
    neutral = M[M.neutral]
    print(f"  {len(M):,} announced meals from {M.subject_id.nunique()} people; "
          f"{len(neutral):,} ended within {NEUTRAL_MGDL:.0f} mg/dL of where they started")

    per = (neutral.groupby(["study", "subject_id"])
           .agg(cr=("cr_implied", "median"), n=("cr_implied", "size"))
           .reset_index())
    per = per[per.n >= MIN_MEALS]
    tdd = B.set_index("subject_id").tdd
    per["tdd"] = per.subject_id.map(tdd)
    per = per.dropna(subset=["tdd"])
    per["cr_x_tdd"] = per.cr * per.tdd
    per.to_parquet(config.RESULTS / "inv009_measured_cr.parquet", index=False)

    ent = db.entered_isf().set_index("subject_id").cr
    per["cr_entered"] = per.subject_id.map(ent)
    print(f"\n  {'cohort':>28s} {'n':>5s} {'CR x TDD':>10s} {'95% interval':>18s} "
          f"{'covers 500':>11s} {'entered CR':>11s} {'measured':>9s}")
    res["cr"] = []
    for label, sub in list(per.groupby("study")) + [("EVERYONE", per)]:
        v = sub.cr_x_tdd.dropna()
        if len(v) < 12:
            continue
        lo, hi = boot_ci(v)
        name = config.COHORTS[label]["label"] if label in config.COHORTS else label
        e = sub.cr_entered.dropna()
        row = dict(cohort=name, n=int(len(v)), median=float(v.median()), ci=[lo, hi],
                   covers_walsh=bool(lo <= WALSH_CR <= hi),
                   median_cr=float(sub.cr.median()),
                   median_entered=float(e.median()) if len(e) else np.nan)
        res["cr"].append(row)
        print(f"  {name:>28s} {row['n']:5d} {row['median']:10.0f} [{lo:6.0f}, {hi:6.0f}]"
              f" {'yes' if row['covers_walsh'] else 'no':>11s} {row['median_entered']:11.1f}"
              f" {row['median_cr']:9.1f}")

    both = per.dropna(subset=["cr_entered"])
    if len(both) > 20:
        r = (both.cr_entered / both.cr).replace([np.inf, -np.inf], np.nan).dropna()
        res["cr_entered_over_measured"] = float(r.median())
        print(f"\n  entered carb ratio divided by measured: {r.median():.2f} "
              f"across {len(r)} people")
        print(f"  a value near 1 would mean no attenuation, which is the test this")
        print(f"  estimator passes or fails on")

    (config.RESULTS / "inv009_measured_basal_cr.json").write_text(
        json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
