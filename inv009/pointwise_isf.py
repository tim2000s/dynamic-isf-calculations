"""Effective sensitivity at every reading, against what the two equations say there.

Everything else in this package fits a slope across a person's windows. This does
the direct thing instead: at each reading, look back over one insulin duration,
work out how much insulin acted in that time from every dose that contributed,
and ask what the glucose change per unit was. Then evaluate both equations at the
same instant from the same data and put the three numbers side by side.

One correction makes the ratio mean something. Over any six hour lookback the
glucose change is insulin action minus hepatic glucose output, and basal exists
to cancel that output. When basal is right, glucose is flat and the raw ratio is
zero, which is a correct basal rate rather than a sensitivity of zero. Taking
the ratio unadjusted returns 3 to 9 mg/dL per unit against entered settings of 25
to 60, which is that term and not a finding.

So both halves are taken as departures from what this person usually does at this
time of day. How much further glucose moved than usual, divided by how much more
insulin acted than usual. What survives is the marginal effect of the insulin that
was not routine, which is the quantity a correction dose is decided with.

    python3 -m inv009.pointwise_isf
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys

import numpy as np
import pandas as pd

from inv008 import dynisf

from . import action, config, db, grid as gridmod, insulin_models as M, tdd as tddmod

MODEL = "oref_6h75"
LOOKBACK_MIN = 360           # one insulin duration
STRIDE_MIN = 30              # readings are 5 min apart and lookbacks overlap heavily
MIN_EXCESS_U = 0.30          # below this the denominator is noise
CARB_CLEAR_H = 6.0           # clear before the lookback opens, and none inside it
MIN_POINTS = 60
BG_BANDS = [(70, 100), (100, 120), (120, 150), (150, 190), (190, 250), (250, 400)]


def _typical_by_halfhour(values: np.ndarray, half_hour: np.ndarray) -> np.ndarray:
    """This person's usual value at this time of day, as a per-point series."""
    return (pd.DataFrame({"hh": half_hour, "v": values})
            .groupby("hh")["v"].transform("median").to_numpy())


def analyse(job):
    subject_id, study = job
    try:
        streams = db.streams(subject_id)
        g = gridmod.build_grid(streams)
        if g is None:
            return None
        subj = tddmod.subject_level(g)
        if not np.isfinite(subj.get("tdd_u", np.nan)):
            return None
        tw = tddmod.windowed(g)

        n = len(g)
        back = int(LOOKBACK_MIN / config.GRID_MIN)
        kern = M.kernel(MODEL)
        total = g.total_u.to_numpy(float)
        cgm = g.cgm.to_numpy(float)
        carbs = g.carbs_g.to_numpy(float)
        has_carbs = bool((carbs > 0).any())

        # Units that acted during each lookback, from every dose contributing to
        # it, whether given before it opened or inside it.
        a_pre, a_in = action.window_action(total, kern, back)
        acted = a_pre + a_in

        # The same quantity for this person's routine delivery, which is what
        # cancels their hepatic output.
        ref = action.reference_profile(g.ts, total)
        r_pre, r_in = action.window_action(ref, kern, back)
        acted_ref = r_pre + r_in

        half_hour = (g.ts.dt.hour * 2 + g.ts.dt.minute // 30).to_numpy()
        dg = np.full(n, np.nan)
        dg[back:] = cgm[back:] - cgm[:-back]          # change across the lookback
        typical_dg = _typical_by_halfhour(dg, half_hour)

        # Indexed by the END of the lookback, so shift the start-indexed action.
        acted_end = np.full(n, np.nan)
        ref_end = np.full(n, np.nan)
        acted_end[back:] = acted[:-back]
        ref_end[back:] = acted_ref[:-back]
        excess = acted_end - ref_end

        # Carbohydrate anywhere in the lookback, or in the hours before it.
        csum_c = np.concatenate([[0.0], np.cumsum(carbs)])
        idx = np.arange(n)
        clear = int(CARB_CLEAR_H * 60 / config.GRID_MIN)
        carbs_in = np.full(n, np.nan)
        carbs_before = np.full(n, np.nan)
        carbs_in[back:] = csum_c[idx[back:]] - csum_c[idx[back:] - back]
        lo = np.maximum(idx - back - clear, 0)
        carbs_before[back:] = csum_c[np.maximum(idx[back:] - back, 0)] - csum_c[lo[back:]]

        keep = np.zeros(n, bool)
        keep[back::int(STRIDE_MIN / config.GRID_MIN)] = True
        keep &= np.isfinite(dg) & np.isfinite(excess) & (np.abs(excess) >= MIN_EXCESS_U)
        keep &= np.isfinite(cgm)
        bg_start = np.full(n, np.nan)
        bg_start[back:] = cgm[:-back]
        keep &= np.isfinite(bg_start) & (bg_start >= 70) & (bg_start <= 400)
        if has_carbs:
            keep &= (np.nan_to_num(carbs_in, nan=1.0) <= 0) & \
                    (np.nan_to_num(carbs_before, nan=1.0) <= 0)
        if keep.sum() < MIN_POINTS:
            return None

        # The equations are evaluated where they would have been evaluated: at the
        # reading the lookback opens on, from the dose blend available then.
        tdd_blend = tw["tdd_blend"].to_numpy(float)
        tdd_at_start = np.full(n, np.nan)
        tdd_at_start[back:] = tdd_blend[:-back]

        sel = np.flatnonzero(keep)
        bg = bg_start[sel]
        td = tdd_at_start[sel]
        ok = np.isfinite(td) & (td > 0)
        sel, bg, td = sel[ok], bg[ok], td[ok]
        if len(sel) < MIN_POINTS:
            return None

        isf_eff = -(dg[sel] - typical_dg[sel]) / excess[sel]
        return pd.DataFrame({
            "subject_id": subject_id, "study": study,
            "bg": bg, "tdd_blend": td, "tdd_u": subj["tdd_u"],
            "excess_u": excess[sel], "acted_u": acted_end[sel],
            "isf_eff": isf_eff,
            "isf_v1": dynisf.isf_v1(bg, td), "isf_v2": dynisf.isf_v2(bg, td),
            "hour": g.ts.dt.hour.to_numpy()[sel],
        })
    except Exception:
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--studies", nargs="*", default=["Loop", "ReplaceBG"],
                    help="default is the two cohorts that recorded meals")
    args = ap.parse_args()
    config.ensure_dirs()

    jobs = []
    for study in args.studies:
        for r in db.subjects(study).itertuples():
            jobs.append((r.subject_id, study))
    print(f"{len(jobs)} people, {config.WORKERS} workers")
    frames = []
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        for i, out in enumerate(pool.imap_unordered(analyse, jobs, chunksize=4), 1):
            if out is not None:
                frames.append(out)
            if i % 200 == 0:
                print(f"  {i}/{len(jobs)}", flush=True)
    D = pd.concat(frames, ignore_index=True)
    D.to_parquet(config.RESULTS / "inv009_pointwise.parquet", index=False)
    print(f"{len(D):,} readings from {D.subject_id.nunique()} people")

    # Per person first, so that somebody with a long record does not outvote the rest.
    res: dict = {"n_points": int(len(D)), "n_subjects": int(D.subject_id.nunique())}
    per = D.groupby("subject_id").agg(isf_eff=("isf_eff", "median"),
                                      isf_v1=("isf_v1", "median"),
                                      isf_v2=("isf_v2", "median"),
                                      tdd_u=("tdd_u", "first"))
    res["overall"] = {k: float(per[k].median()) for k in ("isf_eff", "isf_v1", "isf_v2")}
    print("\nPer-person median sensitivity, mg/dL per unit")
    print(f"  measured {res['overall']['isf_eff']:7.1f}    "
          f"v1 {res['overall']['isf_v1']:7.1f}    v2 {res['overall']['isf_v2']:7.1f}")

    print("\nBy glucose at the start of the lookback")
    print(f"  {'glucose':>12s} {'people':>7s} {'measured':>9s} {'v1':>8s} {'v2':>8s}")
    res["by_bg"] = []
    for lo, hi in BG_BANDS:
        d = D[(D.bg >= lo) & (D.bg < hi)]
        if d.subject_id.nunique() < 20:
            continue
        pp = d.groupby("subject_id")[["isf_eff", "isf_v1", "isf_v2"]].median()
        row = dict(band=f"{lo}-{hi}", n=int(pp.shape[0]),
                   **{k: float(pp[k].median()) for k in pp.columns})
        res["by_bg"].append(row)
        print(f"  {row['band']:>12s} {row['n']:7d} {row['isf_eff']:9.1f} "
              f"{row['isf_v1']:8.1f} {row['isf_v2']:8.1f}")

    print("\nBy total daily dose")
    print(f"  {'dose':>12s} {'people':>7s} {'measured':>9s} {'v1':>8s} {'v2':>8s}")
    res["by_tdd"] = []
    for lo, hi in config.TDD_BANDS:
        d = D[(D.tdd_u >= lo) & (D.tdd_u < hi)]
        if d.subject_id.nunique() < 20:
            continue
        pp = d.groupby("subject_id")[["isf_eff", "isf_v1", "isf_v2"]].median()
        row = dict(band=f"{lo}-{hi if hi < 999 else '+'}", n=int(pp.shape[0]),
                   **{k: float(pp[k].median()) for k in pp.columns})
        res["by_tdd"].append(row)
        print(f"  {row['band']:>12s} {row['n']:7d} {row['isf_eff']:9.1f} "
              f"{row['isf_v1']:8.1f} {row['isf_v2']:8.1f}")

    # How often each equation lands within a useful distance of the measurement.
    for name in ("isf_v1", "isf_v2"):
        ratio = (D[name] / D.isf_eff).replace([np.inf, -np.inf], np.nan).dropna()
        ratio = ratio[(ratio > 0) & (ratio < 100)]
        res[f"{name}_ratio_median"] = float(ratio.median())
        res[f"{name}_within_30pct"] = float(((ratio > 0.7) & (ratio < 1.3)).mean())
    print(f"\n  v1 sits at {res['isf_v1_ratio_median']:.2f} times the measured value, "
          f"within 30% of it {100 * res['isf_v1_within_30pct']:.0f}% of the time")
    print(f"  v2 sits at {res['isf_v2_ratio_median']:.2f} times the measured value, "
          f"within 30% of it {100 * res['isf_v2_within_30pct']:.0f}% of the time")
    (config.RESULTS / "inv009_pointwise.json").write_text(json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
