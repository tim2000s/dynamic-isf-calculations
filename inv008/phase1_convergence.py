#!/usr/bin/env python3
"""Phase 1 of the data-derived-ISF feasibility study: how long until a per-user
empirical-ISF estimate is dosing-grade, and is a short-window estimate reproducible?

Reuses the empirical_isf_v5 window construction and ΔBG = a + b·ΔIOB + c·trend
regression verbatim, refactored to fit on an arbitrary time-subset of a user's history.

Two measures, per user:
  (1) Growing trailing window — empirical ISF on the most recent W days for
      W in {7,14,30,60,90}; record estimate, regression SE, n windows, relative
      95% CI half-width (1.96·se/estimate). Within-window precision.
  (2) Test-retest — empirical ISF on consecutive non-overlapping 14-day blocks;
      between-block coefficient of variation. This is the operational truth for a
      weekly recalibration: does the number actually hold steady, or bounce?

"Dosing-grade" is reported at several relative-precision thresholds rather than one.

Output: results/phase1_convergence.{json,md}
Run: python -m inv008.phase1_convergence
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

# --- identical constants to empirical_isf_v5 ---
COL_MAP = {
    "oref_v5": {"cob": '"sug_COB"'}, "oref_v6": {"cob": "sug_cob"},
    "oref_v7": {"cob": "sug_cob"},
}
WINDOW_S, PRE_S, COB_LEAD_S, SMB_LEAD_S = 1800, 1800, 5400, 1800
COB_TOL, SMB_TOL = 1.0, 0.05
CGM_LO, CGM_HI, DBG_MAX, DIOB_MIN, DIOB_MAX = 70, 250, 80.0, 0.0, 2.0

WINDOWS_D = [7, 14, 30, 60, 90]
BLOCK_D = 14
MIN_W_FIT = 15          # minimum kept windows to attempt a fit on a subset
REL_THRESHOLDS = [0.10, 0.15, 0.20, 0.25]


def _compute_rows(df: pd.DataFrame, table: str):
    """Full-series window mechanics → per-row (ts, keep, diob, dbg, trend)."""
    cm = COL_MAP[table]
    ts = df["ts_relative_sec"].values.astype(float)
    bg = df["cgm_mgdl"].values.astype(float)
    iob = df["iob_iob"].values.astype(float)
    cob = df["cob"].fillna(0).values.astype(float)
    smb = df["sug_smb_units"].fillna(0).values.astype(float)
    n = len(df)

    ends = np.searchsorted(ts, ts + WINDOW_S)
    valid_end = ends < n
    end_idx = np.where(valid_end, ends, 0)
    valid_end &= np.abs(ts[end_idx] - (ts + WINDOW_S)) <= 300

    pre_start = np.searchsorted(ts, ts - PRE_S)
    valid_pre = (pre_start < np.arange(n)) & (
        np.abs(ts[np.clip(pre_start, 0, n - 1)] - (ts - PRE_S)) <= 300)
    pre_idx = np.clip(pre_start, 0, n - 1)

    diob = iob - iob[end_idx]
    dbg = bg[end_idx] - bg

    rolling_cob_max = np.zeros(n)
    j = 0
    for i in range(n):
        while j < i and ts[j] < ts[i] - COB_LEAD_S:
            j += 1
        rolling_cob_max[i] = cob[j:i + 1].max() if i >= j else 0.0

    psum_smb = np.concatenate([[0.0], np.cumsum(smb)])
    rolling_smb_sum = psum_smb[np.arange(n)] - psum_smb[np.searchsorted(ts, ts - SMB_LEAD_S)]
    cob_pos = (cob > COB_TOL).astype(np.int32)
    psum_cp = np.concatenate([[0], np.cumsum(cob_pos)])
    any_cob = psum_cp[end_idx] - psum_cp[np.arange(n)]

    pre_dt = (ts - ts[pre_idx]) / 60.0
    trend = (bg - bg[pre_idx]) / np.where(pre_dt > 0, pre_dt, np.nan)

    keep = (valid_end & valid_pre & (any_cob == 0) & (rolling_cob_max <= COB_TOL)
            & (rolling_smb_sum <= SMB_TOL) & (bg >= CGM_LO) & (bg <= CGM_HI)
            & (bg[end_idx] >= 50) & (bg[end_idx] <= 300) & (np.abs(dbg) <= DBG_MAX)
            & (diob > DIOB_MIN) & (diob < DIOB_MAX) & np.isfinite(trend))
    return ts, keep, diob, dbg, trend


def _fit(mask, diob, dbg, trend):
    """OLS ΔBG = a + b·ΔIOB + c·trend on the masked rows → (isf, se, r2, n)."""
    y = dbg[mask]
    if len(y) < MIN_W_FIT:
        return None
    X = np.column_stack([np.ones(len(y)), diob[mask], trend[mask]])
    try:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    except np.linalg.LinAlgError:
        return None
    pred = X @ beta
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    sigma2 = ss_res / max(len(y) - 3, 1)
    se = float(np.sqrt(sigma2 * np.linalg.pinv(X.T @ X)[1, 1]))
    return dict(isf=float(-beta[1]), se=se, r2=float(r2), n=int(len(y)))


def run_user(args):
    user_id, table = args
    cm = COL_MAP[table]
    sql = (f"SELECT ts_relative_sec, cgm_mgdl, iob_iob, {cm['cob']} AS cob, "
           f"sug_smb_units FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL "
           f"AND iob_iob IS NOT NULL ORDER BY ts_relative_sec")
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(sql, conn, params=(user_id,))
    finally:
        conn.close()
    if len(df) < 500:
        return {"user_id": user_id, "skipped": "few_rows"}
    ts, keep, diob, dbg, trend = _compute_rows(df, table)
    if keep.sum() < MIN_W_FIT:
        return {"user_id": user_id, "skipped": "few_windows", "n_kept": int(keep.sum())}

    t0, t1 = ts.min(), ts.max()
    span_d = (t1 - t0) / 86400.0
    start_ts = ts  # window assignment by window START time

    # (1) growing trailing windows
    trailing = {}
    for W in WINDOWS_D:
        m = keep & (start_ts >= t1 - W * 86400)
        r = _fit(m, diob, dbg, trend)
        if r and r["isf"] > 0:
            r["rel_ci"] = 1.96 * r["se"] / r["isf"]
        trailing[W] = r

    # (2) test-retest over consecutive 14-day blocks
    block_isf = []
    block_se = []
    b = t0
    while b < t1:
        m = keep & (start_ts >= b) & (start_ts < b + BLOCK_D * 86400)
        r = _fit(m, diob, dbg, trend)
        if r and 5 < r["isf"] < 500:
            block_isf.append(r["isf"]); block_se.append(r["se"])
        b += BLOCK_D * 86400
    block_isf = np.array(block_isf); block_se = np.array(block_se)

    rec = {"user_id": user_id, "table": table, "span_days": round(span_d, 1),
           "n_kept_total": int(keep.sum()),
           "trailing": {str(W): trailing[W] for W in WINDOWS_D},
           "n_blocks": int(len(block_isf))}
    if len(block_isf) >= 3:
        rec["block_median_isf"] = float(np.median(block_isf))
        rec["block_cv"] = float(np.std(block_isf) / np.mean(block_isf))      # whole-history spread
        rec["block_within_se_rel"] = float(np.median(block_se) / np.median(block_isf))
        # adjacent-block relative change = the real fortnight-to-fortnight jitter
        adj = np.abs(np.diff(block_isf)) / ((block_isf[:-1] + block_isf[1:]) / 2)
        rec["adj_change"] = float(np.median(adj))
        # lag-1 autocorrelation: high → variation is slow drift (trackable);
        # near 0 → uncorrelated noise (untrackable)
        if len(block_isf) >= 5 and np.std(block_isf) > 0:
            z = block_isf - block_isf.mean()
            rec["block_acf1"] = float((z[:-1] * z[1:]).sum() / (z * z).sum())
    return rec


def main():
    from canonical_cohort import load_canonical_cohort
    coh = load_canonical_cohort()
    coh = coh[coh["in_cohort"]]
    s2t = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}
    work = [(r["user_id"], s2t[r["cohort"]]) for _, r in coh.iterrows()]
    print(f"Phase 1 convergence: {len(work)} users on {N_WORKERS} workers")
    with mp.Pool(N_WORKERS) as pool:
        res = pool.map(run_user, work, chunksize=2)
    OUT.mkdir(exist_ok=True)
    (OUT / "phase1_convergence.json").write_text(json.dumps(res, indent=1))

    ok = [r for r in res if "trailing" in r]
    md = ["# Phase 1 — empirical-ISF convergence & reproducibility\n",
          f"{len(ok)} users with a fit (of {len(work)} canonical).\n",
          "## (1) Growing trailing window — within-window precision\n",
          "| window | n users fit | median rel. 95% CI half-width | "
          + " | ".join(f"≤±{int(t*100)}%" for t in REL_THRESHOLDS) + " |",
          "|---|---|---|" + "---|" * len(REL_THRESHOLDS)]
    for W in WINDOWS_D:
        vals = [r["trailing"][str(W)]["rel_ci"] for r in ok
                if r["trailing"].get(str(W)) and "rel_ci" in r["trailing"][str(W)]]
        vals = np.array([v for v in vals if np.isfinite(v)])
        if not len(vals):
            md.append(f"| {W}d | 0 | – |" + " |" * len(REL_THRESHOLDS)); continue
        fr = [f"{100*np.mean(vals <= t):.0f}%" for t in REL_THRESHOLDS]
        md.append(f"| {W}d | {len(vals)} | ±{100*np.median(vals):.0f}% | " + " | ".join(fr) + " |")

    md.append("\n## (2) Test-retest — reproducibility of a 14-day estimate\n")
    cvs = np.array([r["block_cv"] for r in ok if "block_cv" in r])
    wse = np.array([r["block_within_se_rel"] for r in ok if "block_within_se_rel" in r])
    if len(cvs):
        md.append(f"- users with ≥3 fourteen-day blocks: **{len(cvs)}**")
        md.append(f"- median block-to-block CV (true reproducibility): "
                  f"**±{100*np.median(cvs):.0f}%** [Q1 {100*np.quantile(cvs,.25):.0f}%, "
                  f"Q3 {100*np.quantile(cvs,.75):.0f}%]")
        md.append(f"- median within-block regression SE (the model's own claim): "
                  f"±{100*np.median(wse):.0f}%")
        md.append(f"- → the regression SE {'understates' if np.median(cvs)>np.median(wse) else 'matches'} "
                  f"true variability by ~{np.median(cvs)/max(np.median(wse),1e-9):.1f}×")
        for t in REL_THRESHOLDS:
            md.append(f"- fraction of users reproducible within ±{int(t*100)}% "
                      f"(whole-history CV): {100*np.mean(cvs <= t):.0f}%")
        adj = np.array([r["adj_change"] for r in ok if "adj_change" in r])
        acf = np.array([r["block_acf1"] for r in ok if "block_acf1" in r])
        md.append("\n### Drift vs jitter (adjacent fortnight-to-fortnight)\n")
        md.append(f"- median adjacent-block change: **±{100*np.median(adj):.0f}%** "
                  f"(vs whole-history ±{100*np.median(cvs):.0f}%)")
        md.append(f"- median lag-1 autocorrelation of block estimates: {np.median(acf):+.2f} "
                  f"({'mostly slow drift — trackable' if np.median(acf) > 0.2 else 'mostly uncorrelated noise — not trackable at this resolution'})")
        for t in REL_THRESHOLDS:
            md.append(f"- fraction with adjacent change ≤ ±{int(t*100)}%: "
                      f"{100*np.mean(adj <= t):.0f}%")
    (OUT / "phase1_convergence.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
