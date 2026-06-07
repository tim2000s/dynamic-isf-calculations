#!/usr/bin/env python3
"""Phase 6 pivot: can OVERNIGHT data give a cleaner, per-user sensitivity than the
all-day fasting estimator that Phase 5 showed is hypo-biased?

Overnight (00:00–06:00) is the classic clean window: genuinely post-absorptive, no
meals/activity, stable. We restrict the same ΔIOB regression to nocturnal windows and ask:

  1. Reproducibility — is the overnight estimate steadier block-to-block than the all-day
     ±34% (Phase 1)?
  2. The decisive bias test — re-run Phase 5's directional outcome test on the overnight ISF:
     R = entered/overnight; if overnight sensitivity is unbiased, the +0.38 logR-vs-TBR
     hypo-signal should vanish (and ideally TAR should rise with R). If it persists, overnight
     does not fix the bias.

The estimator's pre-window BG trend covariate absorbs part of the dawn drift.

Output: results/phase6_overnight.{json,md}
Run: python -m inv008.phase6_overnight
"""
from __future__ import annotations

import json
import math
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

from inv008.phase1_convergence import _compute_rows, _fit, COL_MAP

NIGHT = set(range(0, 6))           # 00:00–06:00 (window start hour)
DIOB_FLOOR = 0.25
MIN_NIGHT = 60
BLOCK_D = 14
TABLES = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}


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
        return {"user_id": user_id, "skipped": True}
    ts, keep, diob, dbg, trend = _compute_rows(df, table)
    hour = df["hour"].values.astype(int)
    night = keep & np.isin(hour, list(NIGHT)) & (diob >= DIOB_FLOOR)
    if night.sum() < MIN_NIGHT:
        return {"user_id": user_id, "skipped": True, "n_night": int(night.sum())}

    base = _fit(night, diob, dbg, trend)
    if not base or not (5 < base["isf"] < 500):
        return {"user_id": user_id, "skipped": True}

    # overnight test-retest across 14-day blocks
    t0, t1 = ts.min(), ts.max()
    blocks = []
    b = t0
    while b < t1:
        m = night & (ts >= b) & (ts < b + BLOCK_D * 86400)
        r = _fit(m, diob, dbg, trend)
        if r and 5 < r["isf"] < 500:
            blocks.append(r["isf"])
        b += BLOCK_D * 86400
    blocks = np.array(blocks)
    rec = {"user_id": user_id, "table": table, "n_night": int(night.sum()),
           "overnight_isf": round(base["isf"], 1), "r2": round(base["r2"], 3)}
    if len(blocks) >= 3:
        rec["block_cv"] = float(np.std(blocks) / np.mean(blocks))
    return rec


def main():
    from canonical_cohort import load_canonical_cohort
    from inv008.phase5_outcomes import outcomes
    coh = load_canonical_cohort()
    coh = coh[coh["in_cohort"]].copy()
    work = [(r["user_id"], TABLES[r["cohort"]]) for _, r in coh.iterrows()]
    print(f"Phase 6 overnight: {len(work)} users on {N_WORKERS} workers")
    with mp.Pool(N_WORKERS) as pool:
        res = [r for r in pool.map(run_user, work, chunksize=2) if not r.get("skipped")]
        oc = pd.DataFrame([r for r in pool.map(outcomes, work, chunksize=2) if not r.get("skipped")])
    nightdf = pd.DataFrame(res)

    emp = pd.DataFrame(json.loads((ROOT / "empirical_isf_v5.json").read_text()))[
        ["user_id", "empirical_isf", "r2"]].rename(columns={"r2": "r2_allday"})
    df = (nightdf.merge(coh[["user_id", "isf", "tdd"]], on="user_id")
                 .merge(emp, on="user_id").merge(oc, on="user_id"))
    df = df[df.isf.notna() & (df.isf > 0) & df.overnight_isf.between(5, 500)]
    df["R_night"] = df.isf / df.overnight_isf
    df["logR_night"] = np.log(df.R_night)
    n = len(df)

    def sp(a, b, d=df):
        return float(d[[a, b]].corr(method="spearman").iloc[0, 1])

    # reproducibility
    cvs = nightdf["block_cv"].dropna()
    # bias test (overnight) vs the all-day numbers from Phase 5
    tests = {
        "overnight logR_vs_TAR (expect +)": sp("logR_night", "tar"),
        "overnight logR_vs_TBR (expect -)": sp("logR_night", "tbr"),
        "overnight logR_vs_TIR": sp("logR_night", "tir"),
    }
    # overnight vs all-day ISF level
    df["night_over_allday"] = df.overnight_isf / df.empirical_isf
    df["night_over_entered"] = df.overnight_isf / df.isf

    summary = {
        "n": n,
        "median_overnight_isf": round(float(df.overnight_isf.median()), 1),
        "median_entered_isf": round(float(df.isf.median()), 1),
        "median_allday_empirical": round(float(df.empirical_isf.median()), 1),
        "median_R_entered_over_overnight": round(float(df.R_night.median()), 2),
        "median_overnight_over_entered": round(float(df.night_over_entered.median()), 2),
        "overnight_block_cv_median": round(float(cvs.median()), 3) if len(cvs) else None,
        "allday_block_cv_ref": 0.34,
        "spearman": {k: round(v, 3) for k, v in tests.items()},
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "phase6_overnight.json").write_text(json.dumps(summary, indent=1))

    md = ["# Phase 6 — overnight-derived per-user sensitivity\n",
          f"{n} users with ≥{MIN_NIGHT} nocturnal (00:00–06:00) clean windows.\n",
          "## Level\n",
          f"- median overnight ISF **{summary['median_overnight_isf']}** vs entered "
          f"{summary['median_entered_isf']} vs all-day empirical "
          f"{summary['median_allday_empirical']} mg/dL/U",
          f"- overnight/entered ratio median **{summary['median_overnight_over_entered']}** "
          f"(all-day empirical/entered was ~0.41 — i.e. overnight sits "
          f"{'closer to' if summary['median_overnight_over_entered']>0.6 else 'still below'} the entered value)",
          "\n## Reproducibility (14-day block test-retest)\n",
          f"- overnight block-to-block CV: **±{100*summary['overnight_block_cv_median']:.0f}%** "
          f"(all-day was ±34%) — {'tighter' if summary['overnight_block_cv_median']<0.30 else 'similar'}",
          "\n## Decisive bias test (Phase-5 directional test on the OVERNIGHT ISF)\n",
          "| association | ρ (overnight) | all-day (Phase 5) | expected |",
          "|---|---|---|---|",
          f"| logR vs TAR | {tests['overnight logR_vs_TAR (expect +)']:+.2f} | +0.00 | + |",
          f"| logR vs TBR | {tests['overnight logR_vs_TBR (expect -)']:+.2f} | +0.38 | − |"]
    tbr_sp = tests["overnight logR_vs_TBR (expect -)"]
    md.append("\n## Reading\n")
    if tbr_sp < 0.1:
        md.append("- The all-day **hypo-bias is gone** overnight (logR-vs-TBR no longer positive): "
                  "overnight sensitivity is not entangled with hypo tendency the way the all-day "
                  "estimate was. Overnight is a cleaner per-user sensitivity probe and a safer "
                  "basis for a derived level.")
    else:
        md.append(f"- The hypo-bias **persists overnight** (logR-vs-TBR {tbr_sp:+.2f}): restricting "
                  "to nocturnal windows does not de-bias the sensitivity estimate; overnight alone "
                  "does not solve the Phase-5 problem.")
    md.append("\n*Caveat:* dawn drift partly absorbed by the trend covariate but not eliminated; "
              "overnight ΔIOB is smaller (post-dinner tail), so per-user n is lower; observational.")
    (OUT / "phase6_overnight.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
