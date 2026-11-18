"""Sensitivity variants of the six-hour action-balance outcome proxy of section 4.5.

`inv009/pointwise_isf.py` takes six-hour intervals in Loop and REPLACE-BG,
divides the glucose change beyond the person's usual change at that time of day
by the insulin action beyond a reference, and reports the across-people median by
starting-glucose band. The proxy rose from -1.8 mg/dL per U at 70 to 100 mg/dL to
24.3 at 250 to 400. Review 1, point M2, asked whether that rise survives a
different basal reference and an exposure that leaves out insulin delivered
inside the interval. This script reruns the same estimator, with the same
records, insulin models, carbohydrate screen, stride and 0.30 U minimum, in four
arms:

  A  programmed basal as the reference, all insulin acting in the interval
     (the committed arm, rerun to check that it reproduces the committed counts);
  B  the person's typical half-hourly delivery as the reference, which is the
     fallback `forward_isf.schedule_units` already uses where no programme is
     recoverable: the median delivered units at each half hour of day over the
     whole record;
  C  programmed basal as the reference, with exposure restricted to insulin
     committed before the interval opens: the action, during the interval, of
     doses given before it (`a_pre` of `action.window_action`), less the same
     quantity for the reference;
  D  the typical-delivery reference with the pre-interval exposure of C.

The 0.30 U minimum is applied to each arm's own excess, so the intervals kept
differ between arms. Every band median is given with a 95% percentile interval
from 2,000 bootstrap draws that resample people, and the test of the rise is the
across-people median in the 250 to 400 mg/dL band minus that in the 70 to
100 mg/dL band, with its interval and a verdict.

Expected sign, written before the run. Arm B removes the controller's typical
response from the excess along with basal; in Loop much of that typical delivery is
temporary basal raised in response to glucose, so the excess shrinks most where
glucose is high, which raises the proxy there. The rise is expected to survive and
steepen. Arms C and D remove insulin delivered in response to glucose inside the
interval. Near target that insulin is commonly given after a rise, which is what
drives the proxy negative there, so removing it is expected to raise the low
bands. Pre-interval insulin is larger when the interval opens high, so the high
bands are expected to stay positive. The rise is expected to survive in C and D
with a smaller difference between the top and bottom bands than in A. A
difference whose interval includes zero would mean the rise does not survive
that change.

The proxy is an outcome ratio and its level is not read as a physiological ISF.

    python3 -m inv009.pointwise_variants            # DYNISF_ROOT left unset

Outputs results/inv009_pointwise_variants.{json,md} (aggregates only).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from inv008 import dynisf

from . import action, config, db, grid as gridmod, insulin_models as M, tdd as tddmod
from .forward_isf import _loop_model, schedule_units
from .pointwise_isf import (BG_BANDS, CARB_CLEAR_H, LOOKBACK_MIN, MIN_EXCESS_U, MIN_POINTS,
                            MODEL_DEFAULT, STRIDE_MIN, _typical_by_halfhour)

ARMS = {"A": ("programmed", "all"), "B": ("typical", "all"),
        "C": ("programmed", "pre"), "D": ("typical", "pre")}
ARM_LABEL = {"A": "programmed basal, all insulin (committed)",
             "B": "typical half-hourly delivery, all insulin",
             "C": "programmed basal, insulin committed before the interval",
             "D": "typical half-hourly delivery, insulin committed before the interval"}
N_BOOT, SEED = 2000, 0
OUT = Path(__file__).resolve().parent.parent / "results"


def typical_units(g: pd.DataFrame) -> np.ndarray:
    """Median delivered units per bin at each half hour of day, as schedule_units' fallback."""
    hh = (g.ts.dt.hour.to_numpy() * 2 + (g.ts.dt.minute.to_numpy() >= 30)).astype(int)
    tot = g.total_u.to_numpy(float)
    return pd.Series(tot).groupby(hh).median().reindex(hh).to_numpy(dtype=float)


def end_indexed(x: np.ndarray, back: int) -> np.ndarray:
    out = np.full(len(x), np.nan)
    out[back:] = x[:-back]
    return out


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
        model = _loop_model(subject_id) if study == "Loop" else MODEL_DEFAULT
        kern = M.kernel(model)
        total = g.total_u.to_numpy(float)
        cgm = g.cgm.to_numpy(float)
        carbs = g.carbs_g.to_numpy(float)
        has_carbs = bool((carbs > 0).any())

        a_pre, a_in = action.window_action(total, kern, back)
        prog, src = schedule_units(g, subject_id, streams)
        refs = {"programmed": prog, "typical": typical_units(g)}
        excess = {}
        for arm, (ref_name, expo) in ARMS.items():
            r_pre, r_in = action.window_action(refs[ref_name], kern, back)
            acted = a_pre + a_in if expo == "all" else a_pre
            acted_ref = r_pre + r_in if expo == "all" else r_pre
            excess[arm] = end_indexed(acted, back) - end_indexed(acted_ref, back)

        half_hour = (g.ts.dt.hour * 2 + g.ts.dt.minute // 30).to_numpy()
        dg = np.full(n, np.nan)
        dg[back:] = cgm[back:] - cgm[:-back]
        typical_dg = _typical_by_halfhour(dg, half_hour)

        csum_c = np.concatenate([[0.0], np.cumsum(carbs)])
        idx = np.arange(n)
        clear = int(CARB_CLEAR_H * 60 / config.GRID_MIN)
        carbs_in = np.full(n, np.nan)
        carbs_before = np.full(n, np.nan)
        carbs_in[back:] = csum_c[idx[back:]] - csum_c[idx[back:] - back]
        lo = np.maximum(idx - back - clear, 0)
        carbs_before[back:] = csum_c[np.maximum(idx[back:] - back, 0)] - csum_c[lo[back:]]

        base = np.zeros(n, bool)
        base[back::int(STRIDE_MIN / config.GRID_MIN)] = True
        base &= np.isfinite(dg) & np.isfinite(cgm)
        bg_start = end_indexed(cgm, back)
        base &= np.isfinite(bg_start) & (bg_start >= 70) & (bg_start <= 400)
        if has_carbs:
            base &= (np.nan_to_num(carbs_in, nan=1.0) <= 0) & \
                    (np.nan_to_num(carbs_before, nan=1.0) <= 0)
        tdd_at_start = end_indexed(tw["tdd_blend"].to_numpy(float), back)
        base &= np.isfinite(tdd_at_start) & (tdd_at_start > 0)

        keep = {a: base & np.isfinite(e) & (np.abs(e) >= MIN_EXCESS_U) for a, e in excess.items()}
        # the committed per-person minimum is applied to each arm separately
        keep = {a: (k if k.sum() >= MIN_POINTS else np.zeros(n, bool)) for a, k in keep.items()}
        anyk = np.logical_or.reduce(list(keep.values()))
        if not anyk.any():
            return None
        sel = np.flatnonzero(anyk)
        bg, td = bg_start[sel], tdd_at_start[sel]
        out = pd.DataFrame({"subject_id": subject_id, "study": study, "bg": bg,
                            "isf_v1": dynisf.isf_v1(bg, td), "isf_v2": dynisf.isf_v2(bg, td)})
        for a in ARMS:
            out[f"keep_{a}"] = keep[a][sel]
            with np.errstate(divide="ignore", invalid="ignore"):
                out[f"proxy_{a}"] = -(dg[sel] - typical_dg[sel]) / excess[a][sel]
            out[f"excess_{a}"] = excess[a][sel]
        return out
    except Exception:
        return None


def band_medians(D: pd.DataFrame, arm: str) -> pd.DataFrame:
    """Per person, per band: median proxy, V1 and V2 over that arm's intervals."""
    d = D[D[f"keep_{arm}"]]
    band = np.full(len(d), -1)
    for k, (lo, hi) in enumerate(BG_BANDS):
        band[(d.bg >= lo).to_numpy() & (d.bg < hi).to_numpy()] = k
    d = d.assign(band=band)[band >= 0]
    return d.groupby(["subject_id", "band"]).agg(proxy=(f"proxy_{arm}", "median"),
                                                 v1=("isf_v1", "median"),
                                                 v2=("isf_v2", "median")).reset_index()


def summarise(D: pd.DataFrame, arm: str, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    pb = band_medians(D, arm)
    wide = pb.pivot(index="subject_id", columns="band", values="proxy")
    people = wide.index.to_numpy()
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(people), (n_boot, len(people)))
    W = wide.to_numpy(float)
    res = {"n_people": int(D.loc[D[f"keep_{arm}"], "subject_id"].nunique()),
           "n_intervals": int(D[f"keep_{arm}"].sum()), "bands": []}
    for k, (lo, hi) in enumerate(BG_BANDS):
        if k not in wide.columns:
            continue
        col = list(wide.columns).index(k)
        v = W[:, col]
        meds = np.nanmedian(np.where(np.isfinite(v[idx]), v[idx], np.nan), axis=1)
        sub = pb[pb.band == k]
        res["bands"].append(dict(band=f"{lo}-{hi}", n=int(np.isfinite(v).sum()),
                                 proxy=float(np.nanmedian(v)),
                                 lo=float(np.percentile(meds, 2.5)), hi=float(np.percentile(meds, 97.5)),
                                 v1=float(sub.v1.median()), v2=float(sub.v2.median())))
    c0, ct = list(wide.columns).index(0), list(wide.columns).index(len(BG_BANDS) - 1)
    d = np.nanmedian(W[idx, ct], axis=1) - np.nanmedian(W[idx, c0], axis=1)
    point = float(np.nanmedian(W[:, ct]) - np.nanmedian(W[:, c0]))
    dlo, dhi = np.percentile(d, [2.5, 97.5])
    res["top_minus_lowest"] = dict(diff=point, lo=float(dlo), hi=float(dhi),
                                   verdict=("interval above zero" if dlo > 0 else
                                            "interval below zero" if dhi < 0 else
                                            "interval includes zero"))
    k = D[D[f"keep_{arm}"]]
    for name in ("isf_v1", "isf_v2"):
        r = (k[name] / k[f"proxy_{arm}"]).replace([np.inf, -np.inf], np.nan).dropna()
        r = r[(r > 0) & (r < 100)]
        res[f"{name}_ratio_median"] = float(r.median())
    return res


def self_check() -> None:
    """Pre-only exposure must exclude a dose given inside the interval."""
    kern = M.kernel(MODEL_DEFAULT)
    u = np.zeros(400)
    u[150] = 2.0                                  # one dose, inside an interval opening at 100
    a_pre, a_in = action.window_action(u, kern, 72)
    assert abs(a_pre[100]) < 1e-12 and a_in[100] > 0.5
    assert a_pre[160] > 0                         # but before an interval opening at 160
    g = pd.DataFrame({"ts": pd.date_range("2020-01-01", periods=576, freq="5min"),
                      "total_u": np.tile(np.r_[np.full(144, 0.1), np.full(144, 0.05)], 2)})
    t = typical_units(g)
    assert np.allclose(t[:144], 0.1) and np.allclose(t[144:288], 0.05)


def markdown(res: dict) -> str:
    L = ["# Six-hour action-balance proxy under other references and exposures", "",
         "Generated by `inv009/pointwise_variants.py`. Across-people median of each person's median "
         "proxy in each starting-glucose band, mg/dL per U, with a 95% percentile interval from "
         f"{N_BOOT} bootstrap draws resampling people. The proxy is an outcome ratio, not a "
         "measured ISF. V1 and V2 are across-people medians at the same intervals.", ""]
    for arm, r in res["arms"].items():
        L += [f"## Arm {arm}: {ARM_LABEL[arm]}", "",
              f"{r['n_people']} people, {r['n_intervals']:,} intervals. Pooled point-level ratio "
              f"V1/proxy {r['isf_v1_ratio_median']:.2f}, V2/proxy {r['isf_v2_ratio_median']:.2f}.", "",
              "| Starting glucose | People | Proxy (95% interval) | V1 | V2 |", "|---|---:|---|---:|---:|"]
        for b in r["bands"]:
            L.append(f"| {b['band']} mg/dL | {b['n']} | {b['proxy']:.1f} ({b['lo']:.1f} to {b['hi']:.1f}) | "
                     f"{b['v1']:.1f} | {b['v2']:.1f} |")
        t = r["top_minus_lowest"]
        L += ["", f"250 to 400 minus 70 to 100: {t['diff']:+.1f} mg/dL per U "
                  f"({t['lo']:+.1f} to {t['hi']:+.1f}), {t['verdict']}.", ""]
        for s, rs in r.get("by_study", {}).items():
            t = rs["top_minus_lowest"]
            L.append(f"- {s}: {rs['n_people']} people; 70 to 100 {rs['bands'][0]['proxy']:.1f}, "
                     f"250 to 400 {rs['bands'][-1]['proxy']:.1f}; difference {t['diff']:+.1f} "
                     f"({t['lo']:+.1f} to {t['hi']:+.1f}), {t['verdict']}.")
        L.append("")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--studies", nargs="*", default=["Loop", "ReplaceBG"])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--raw", type=Path, default=None,
                    help="where to write the per-interval table (kept out of git)")
    a = ap.parse_args()
    self_check()
    jobs = []
    for study in a.studies:
        ids = db.subjects(study).subject_id.tolist()
        jobs += [(s, study) for s in (ids[:a.limit] if a.limit else ids)]
    t0 = time.time()
    frames = []
    with mp.Pool(min(a.workers, 7), maxtasksperchild=8) as pool:
        for i, out in enumerate(pool.imap_unordered(analyse, jobs, chunksize=2), 1):
            if out is not None:
                frames.append(out)
            if i % 100 == 0:
                print(f"  {i}/{len(jobs)} {time.time() - t0:.0f} s", flush=True)
    D = pd.concat(frames, ignore_index=True)
    if a.raw:
        D.to_parquet(a.raw, index=False)          # raw output first (REVIEW.md item 16)
    res = {"n_jobs": len(jobs), "arms": {}}
    for arm in ARMS:
        r = summarise(D, arm)
        r["by_study"] = {s: summarise(D[D.study == s], arm) for s in a.studies}
        res["arms"][arm] = r
        print(f"  arm {arm}: {r['n_people']} people, {r['n_intervals']:,} intervals", flush=True)
    # arms that agree to every decimal have usually not received their parameter
    for x, y in (("A", "B"), ("A", "C"), ("C", "D")):
        same = np.allclose(D[f"excess_{x}"], D[f"excess_{y}"], equal_nan=True)
        res[f"excess_identical_{x}_{y}"] = bool(same)
        assert not same, f"arms {x} and {y} have identical excess"
    a.out.mkdir(exist_ok=True)
    (a.out / "inv009_pointwise_variants.json").write_text(json.dumps(res, indent=1, default=float))
    md = markdown(res)
    (a.out / "inv009_pointwise_variants.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
