#!/usr/bin/env python3
"""Phase 3 (the crux): the baseline empirical-Bayes shrinkage estimator.

The derivable target is each user's stable BASELINE sensitivity. From limited data we
estimate it by shrinking the user's own W-day empirical ISF toward the population prior
K_pop/√TDD. This quantifies, as a function of how long a user has run:

  - prior-only error  : |log(K_pop/√TDD) − log(baseline)|         (no own data)
  - own-only error    : |log(W-day empirical ISF) − log(baseline)| (no prior)
  - shrinkage error   : own/prior blended at the empirical-Bayes weight w = τ²/(τ²+σ²_W)

Honest noise: σ²_W is MEASURED as the spread of a user's W-day window estimates around their
own long-run baseline (Phase 1 showed the regression SE understates this ~3.3×), so the
shrinkage weight and errors reflect real reproducibility, not the optimistic regression CI.

baseline_i  = empirical ISF over the user's full clean-window history
prior_i     = K_pop / √TDD_i,  K_pop = leave-one-out median of (baseline·√TDD)
τ²          = between-user variance of log(baseline·√TDD / K_pop)  (prior's irreducible spread)

Output: results/phase3_baseline.{json,md}
Run: python -m inv008.phase3_baseline
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

from inv008.phase1_convergence import _compute_rows, _fit, COL_MAP, MIN_W_FIT

WINDOWS_D = [7, 14, 30, 60, 90]
ISF_LO, ISF_HI = 5.0, 500.0


def run_user(args):
    user_id, table, tdd = args
    cm = COL_MAP[table]
    sql = (f"SELECT ts_relative_sec, cgm_mgdl, iob_iob, {cm['cob']} AS cob, "
           f"sug_smb_units FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL "
           f"AND iob_iob IS NOT NULL ORDER BY ts_relative_sec")
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(sql, conn, params=(user_id,))
    finally:
        conn.close()
    if len(df) < 1000 or tdd is None or tdd <= 0:
        return {"user_id": user_id, "skipped": "insufficient"}
    ts, keep, diob, dbg, trend = _compute_rows(df, table)
    if keep.sum() < 200:
        return {"user_id": user_id, "skipped": "few_windows"}

    base = _fit(keep, diob, dbg, trend)        # full-history baseline
    if not base or not (ISF_LO < base["isf"] < ISF_HI):
        return {"user_id": user_id, "skipped": "bad_baseline"}
    log_base = math.log(base["isf"])
    t0, t1 = ts.min(), ts.max()

    # per-W: non-overlapping window estimates → log-deviations from baseline
    win = {}
    for W in WINDOWS_D:
        devs, ests = [], []
        b = t0
        while b < t1:
            m = keep & (ts >= b) & (ts < b + W * 86400)
            r = _fit(m, diob, dbg, trend)
            if r and ISF_LO < r["isf"] < ISF_HI:
                ests.append(math.log(r["isf"])); devs.append(math.log(r["isf"]) - log_base)
            b += W * 86400
        win[W] = ests if len(ests) >= 1 else None
    return {"user_id": user_id, "table": table, "tdd": float(tdd),
            "baseline_isf": round(base["isf"], 1), "log_base": log_base,
            "win_log_est": {str(W): win[W] for W in WINDOWS_D}}


def main():
    from canonical_cohort import load_canonical_cohort
    coh = load_canonical_cohort()
    coh = coh[coh["in_cohort"]]
    s2t = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}
    work = [(r["user_id"], s2t[r["cohort"]], r.get("tdd")) for _, r in coh.iterrows()]
    print(f"Phase 3 baseline: {len(work)} users on {N_WORKERS} workers")
    with mp.Pool(N_WORKERS) as pool:
        res = pool.map(run_user, work, chunksize=2)
    OUT.mkdir(exist_ok=True)

    ok = [r for r in res if "log_base" in r]
    logK = np.array([r["log_base"] + 0.5 * math.log(r["tdd"]) for r in ok])  # log(baseline·√TDD)
    Kpop_log = float(np.median(logK))
    tau2 = float(np.var(logK))                                # prior's between-user spread
    # prior error per user = log(prior) − log(baseline) = Kpop_log − 0.5 ln(tdd) − log_base = -(logK_i - Kpop_log)
    prior_err = np.abs(logK - Kpop_log)

    def pct(x):  # median |log| → ± percent
        return 100 * (math.exp(np.median(x)) - 1)

    rows = []
    for W in WINDOWS_D:
        own_devs, shrink_errs = [], []
        ws = []
        for r in ok:
            ests = r["win_log_est"].get(str(W))
            if not ests or len(ests) < 2:
                continue
            ests = np.array(ests)
            sigma2 = float(np.var(ests - r["log_base"]))      # honest within-user W-day error²
            w = tau2 / (tau2 + sigma2) if (tau2 + sigma2) > 0 else 0.0
            ws.append(w)
            log_prior = Kpop_log - 0.5 * math.log(r["tdd"])
            for e in ests:
                own_devs.append(abs(e - r["log_base"]))
                shrink = w * e + (1 - w) * log_prior
                shrink_errs.append(abs(shrink - r["log_base"]))
        if not own_devs:
            continue
        rows.append({
            "W": W, "n_users": len(ws),
            "w_median": round(float(np.median(ws)), 2),
            "prior_only_pct": round(pct(prior_err), 0),
            "own_only_pct": round(pct(np.array(own_devs)), 0),
            "shrinkage_pct": round(pct(np.array(shrink_errs)), 0),
        })

    summary = {"n_users": len(ok), "K_pop": round(math.exp(Kpop_log), 1),
               "tau_log": round(math.sqrt(tau2), 3),
               "prior_only_pct": round(pct(prior_err), 0), "by_window": rows}
    (OUT / "phase3_baseline.json").write_text(json.dumps(summary, indent=1))

    md = ["# Phase 3 — baseline empirical-Bayes shrinkage estimator\n",
          f"{len(ok)} users. Population prior K_pop = **{summary['K_pop']}** (ISF·√TDD), "
          f"between-user spread τ = {summary['tau_log']} in log "
          f"(prior-only baseline error **±{summary['prior_only_pct']:.0f}%**).\n",
          "Error = median |log(estimate) − log(full-history baseline)|, as ± percent. "
          "Own-data error uses MEASURED window-to-baseline spread (not the optimistic "
          "regression SE). Shrinkage weight w = τ²/(τ²+σ²_W).\n",
          "| trailing window | n users | shrink weight w | prior-only | own-only | **shrinkage** |",
          "|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['W']}d | {r['n_users']} | {r['w_median']:.2f} | "
                  f"±{r['prior_only_pct']:.0f}% | ±{r['own_only_pct']:.0f}% | "
                  f"**±{r['shrinkage_pct']:.0f}%** |")
    md.append("")
    # crossover + floor commentary
    cross = next((r["W"] for r in rows if r["own_only_pct"] < r["prior_only_pct"]), None)
    best = min(rows, key=lambda r: r["shrinkage_pct"]) if rows else None
    md.append(f"- Own-data overtakes the prior at **~{cross} days**." if cross
              else "- Own-data never clearly beats the prior in this cohort.")
    if best:
        md.append(f"- Best achievable baseline error here: **±{best['shrinkage_pct']:.0f}%** "
                  f"at {best['W']}-day windows (shrinkage).")
    md.append(f"- Prior-only (zero own-data, cold start) baseline error: "
              f"**±{summary['prior_only_pct']:.0f}%** — the floor a brand-new user starts at.")
    (OUT / "phase3_baseline.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
