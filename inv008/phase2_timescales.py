#!/usr/bin/env python3
"""Phase 2 of the data-derived-ISF study: decompose within-user sensitivity variance
by timescale, to see what is a stable baseline, what is deterministic-trackable
(circadian/weekly), and what is irreducible stochastic — and whether the fast
variation is persistent (trackable) or white (not), real or artefact.

Per clean fasting window (same machinery as Phase 1) with ΔIOB ≥ 0.3 U we form a
per-window "local ISF" = -(ΔBG - c·trend)/ΔIOB using the user's own fitted trend
coefficient. We centre log(local ISF) per user (removes the stable baseline) and
decompose the remaining variance:

  - circadian   : hour-of-day (DB `hour` column — no anchor needed)
  - weekly      : relative day-of-week (record-anchored)
  - slow        : ~monthly index (record-relative)
  - residual    : everything else (true day-to-day stochastic + measurement noise)

Then per-user autocorrelation of the *daily* sensitivity series at lag 1/7/28
separates persistent (trackable by a fast EMA, i.e. autosens) from white noise.
A stricter-ΔIOB pass separates real physiology (circadian survives) from artefact
(residual shrinks).

Output: results/phase2_timescales.{json,md}
Run: python -m inv008.phase2_timescales
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
sys.path.insert(0, str(ROOT))
OUT = ROOT / "results"
DSN = "dbname=oref"
N_WORKERS = min(12, mp.cpu_count())

from inv008.phase1_convergence import _compute_rows, COL_MAP  # reuse window mechanics

DIOB_FLOOR = 0.3
ISF_LO, ISF_HI = 3.0, 600.0
MIN_WINDOWS = 150
HOUR_BUCKETS = {**{h: "night" for h in range(0, 6)}, **{h: "dawn" for h in range(6, 10)},
                **{h: "midday" for h in range(10, 14)}, **{h: "afternoon" for h in range(14, 18)},
                **{h: "evening" for h in range(18, 22)}, **{h: "late" for h in range(22, 24)}}


def _local_isf(df, table, diob_floor):
    ts, keep, diob, dbg, trend = _compute_rows(df, table)
    hour = df["hour"].values.astype(float)
    m = keep & (diob >= diob_floor)
    if m.sum() < MIN_WINDOWS:
        return None
    # user baseline regression to get the trend coefficient c
    X = np.column_stack([np.ones(m.sum()), diob[m], trend[m]])
    beta, *_ = np.linalg.lstsq(X, dbg[m], rcond=None)
    c = beta[2]
    local = -(dbg[m] - c * trend[m]) / diob[m]
    ok = np.isfinite(local) & (local >= ISF_LO) & (local <= ISF_HI)
    out = pd.DataFrame({
        "ts": ts[m][ok], "hour": hour[m][ok], "isf": local[ok],
    })
    if len(out) < MIN_WINDOWS:
        return None
    out["logisf"] = np.log(out.isf)
    out["dev"] = out.logisf - out.logisf.mean()           # per-user centred (baseline removed)
    out["bucket"] = out.hour.astype(int).map(HOUR_BUCKETS)
    day = ((out.ts - out.ts.min()) // 86400).astype(int)
    out["day"] = day
    out["dow"] = (day % 7)
    out["month"] = (day // 30)
    return out


def _eta2(df, factor):
    """Fraction of dev variance explained by group means of `factor` (one-way η²)."""
    ss_tot = float((df.dev ** 2).sum())
    if ss_tot <= 0:
        return 0.0
    ss_between = float(df.groupby(factor).dev.apply(lambda s: len(s) * s.mean() ** 2).sum())
    return ss_between / ss_tot


def _acf(x, lag):
    x = np.asarray(x, dtype=float)
    if len(x) <= lag + 3 or np.std(x) == 0:
        return None
    z = x - x.mean()
    return float((z[:-lag] * z[lag:]).sum() / (z * z).sum())


def run_user(args):
    user_id, table = args
    cm = COL_MAP[table]
    sql = (f"SELECT ts_relative_sec, cgm_mgdl, iob_iob, {cm['cob']} AS cob, "
           f"sug_smb_units, hour FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL "
           f"AND iob_iob IS NOT NULL ORDER BY ts_relative_sec")
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(sql, conn, params=(user_id,))
    finally:
        conn.close()
    if len(df) < 1000:
        return {"user_id": user_id, "skipped": "few_rows"}

    rec = {"user_id": user_id, "table": table}
    for tag, floor in (("base", DIOB_FLOOR), ("strict", 0.6)):
        w = _local_isf(df, table, floor)
        if w is None:
            rec[tag] = None
            continue
        # daily median sensitivity series for autocorrelation
        daily = w.groupby("day").agg(dev=("dev", "median"), n=("dev", "size"))
        daily = daily[daily.n >= 5]["dev"]
        rec[tag] = {
            "n_windows": int(len(w)),
            "total_var": float((w.dev ** 2).mean()),
            "eta2_circadian": round(_eta2(w, "bucket"), 4),
            "eta2_weekly": round(_eta2(w, "dow"), 4),
            "eta2_slow": round(_eta2(w, "month"), 4),
            "n_days": int(len(daily)),
            "daily_acf1": _acf(daily.values, 1),
            "daily_acf7": _acf(daily.values, 7),
            # circadian shape: mean dev by bucket (the deterministic curve)
            "circadian_profile": {b: round(float(v), 3)
                                  for b, v in w.groupby("bucket").dev.mean().items()},
        }
    return rec


def main():
    from canonical_cohort import load_canonical_cohort
    coh = load_canonical_cohort()
    coh = coh[coh["in_cohort"]]
    s2t = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}
    work = [(r["user_id"], s2t[r["cohort"]]) for _, r in coh.iterrows()]
    print(f"Phase 2 timescales: {len(work)} users on {N_WORKERS} workers")
    with mp.Pool(N_WORKERS) as pool:
        res = pool.map(run_user, work, chunksize=2)
    OUT.mkdir(exist_ok=True)
    (OUT / "phase2_timescales.json").write_text(json.dumps(res, indent=1))

    ok = [r for r in res if r.get("base")]
    B = [r["base"] for r in ok]

    def med(key):
        v = np.array([b[key] for b in B if b.get(key) is not None], dtype=float)
        return np.median(v), v

    md = ["# Phase 2 — sensitivity variance decomposition by timescale\n",
          f"{len(ok)} users with ≥{MIN_WINDOWS} clean windows (ΔIOB ≥ {DIOB_FLOOR} U).",
          "Per-window log local-ISF, centred per user (baseline removed); variance attributed "
          "to each timescale (one-way η²); daily-series autocorrelation separates persistent "
          "(trackable) from white.\n",
          "## Variance budget (median across users, fraction of within-user log-sensitivity variance)\n"]
    for key, label in (("eta2_circadian", "circadian (hour-of-day)"),
                       ("eta2_weekly", "weekly (day-of-week)"),
                       ("eta2_slow", "slow (~monthly)")):
        m, v = med(key)
        md.append(f"- **{label}**: {100*m:.0f}%  [Q1 {100*np.quantile(v,.25):.0f}%, "
                  f"Q3 {100*np.quantile(v,.75):.0f}%]")
    circ = np.array([b["eta2_circadian"] for b in B])
    wk = np.array([b["eta2_weekly"] for b in B])
    sl = np.array([b["eta2_slow"] for b in B])
    resid = 1 - (circ + wk + sl)
    md.append(f"- **residual (day-to-day stochastic + measurement noise)**: "
              f"~{100*np.median(resid):.0f}%")

    md.append("\n## Is the day-to-day residual trackable or white?\n")
    a1, a1v = med("daily_acf1")
    a7, a7v = med("daily_acf7")
    md.append(f"- median daily lag-1 autocorrelation: **{a1:+.2f}** "
              f"({'persistent → trackable by a fast EMA (autosens)' if a1 > 0.2 else 'weak/none → largely white at daily scale'})")
    md.append(f"- median daily lag-7 autocorrelation: {a7:+.2f}")

    md.append("\n## Circadian shape (population mean dev by bucket, log units)\n")
    order = ["night", "dawn", "midday", "afternoon", "evening", "late"]
    pooled = {b: [] for b in order}
    for b in B:
        for k, v in b["circadian_profile"].items():
            pooled[k].append(v)
    md.append("| bucket | mean log-dev | ≈ relative ISF |")
    md.append("|---|---|---|")
    for b in order:
        if pooled[b]:
            mdev = np.median(pooled[b])
            md.append(f"| {b} | {mdev:+.3f} | ×{np.exp(mdev):.2f} |")

    md.append("\n## Real vs artefact (stricter ΔIOB ≥ 0.6 U)\n")
    S = [r["strict"] for r in res if r.get("strict")]
    if S:
        cs = np.median([s["eta2_circadian"] for s in S])
        rs = np.median([1 - (s["eta2_circadian"] + s["eta2_weekly"] + s["eta2_slow"]) for s in S])
        md.append(f"- circadian η²: {100*np.median(circ):.0f}% (base) → {100*cs:.0f}% (strict)")
        md.append(f"- residual: {100*np.median(resid):.0f}% (base) → {100*rs:.0f}% (strict)")
        md.append(f"- → {'residual shrinks under stricter filter ⇒ part of it is estimator artefact' if rs < np.median(resid) - 0.03 else 'residual stable ⇒ largely real day-to-day variability'}; "
                  f"circadian {'survives' if cs > 0.5*np.median(circ) else 'weakens'} ⇒ "
                  f"{'real physiology' if cs > 0.5*np.median(circ) else 'check'}.")

    (OUT / "phase2_timescales.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
