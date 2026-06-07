#!/usr/bin/env python3
"""Empirical-ISF v5 — uses ΔIOB (unambiguous: U absorbed in window) instead of
the unit-uncertain `iob_activity`.

Per 30-min fasting window:
  • ΔIOB = iob_iob[start] − iob_iob[end]
  • ΔBG  = cgm[end] − cgm[start]
  • Pre-window BG trend (mg/dL/min over the preceding 30 min)
  • Regression: ΔBG = a + b · ΔIOB + c · pre_trend
  • Empirical ISF = b   (mg/dL per U absorbed; positive because we expect
    BG drop when ΔIOB > 0, i.e. insulin absorbed → BG falls)

  Fasting window definition (same as v3/v4):
    sug_COB == 0 throughout, rolling 90-min preceding COB max ≤ 1g, pre-window
    SMB sum ≤ 0.05 U.

Reads canonical cohort from `canonical_cohort.load_canonical_cohort()` so
DynISF group / entered ISF / entered TDD come from the same authoritative source.

Output:
  empirical_isf_v5.{md,json}

Parallelism: multiprocessing.Pool (capped 12).
"""
from __future__ import annotations

import json
import multiprocessing as mp
import sys
import warnings
import os
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
sys.path.insert(0, str(ROOT))

OUT_JSON = ROOT / "empirical_isf_v5.json"
OUT_MD   = ROOT / "empirical_isf_v5.md"
DSN = "dbname=oref"
N_WORKERS = min(12, mp.cpu_count())

COL_MAP = {
    "oref_v5": {"cob": '"sug_COB"', "tgt": "sug_current_target",
                "sens": '"sug_sensitivityRatio"'},
    "oref_v6": {"cob": "sug_cob", "tgt": "sug_current_target",
                "sens": "sug_sensitivityratio"},
    "oref_v7": {"cob": "sug_cob", "tgt": "sug_current_target",
                "sens": "sug_sensitivityratio"},
}

WINDOW_S = 1800
PRE_S    = 1800
COB_LEAD_S = 5400
SMB_LEAD_S = 1800
COB_TOL = 1.0
SMB_TOL = 0.05
CGM_LO, CGM_HI = 70, 250
DBG_MAX = 80.0
DIOB_MIN = 0.0
DIOB_MAX = 2.0
MIN_WINDOWS = 80


def fit_user(args):
    user_id, table = args
    cm = COL_MAP[table]
    sql = f"""
        SELECT ts_relative_sec, cgm_mgdl, iob_iob,
               {cm['cob']} AS cob, {cm['tgt']} AS target,
               has_dynisf, {cm['sens']} AS sens, sug_smb_units
        FROM {table}
        WHERE user_id = %s
          AND cgm_mgdl IS NOT NULL
          AND iob_iob IS NOT NULL
        ORDER BY ts_relative_sec
    """
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(sql, conn, params=(user_id,))
    finally:
        conn.close()
    n = len(df)
    if n < 200:
        return {"user_id": user_id, "table": table, "skipped": "few_rows"}

    ts = df["ts_relative_sec"].values.astype(float)
    bg = df["cgm_mgdl"].values.astype(float)
    iob = df["iob_iob"].values.astype(float)
    cob = df["cob"].fillna(0).values.astype(float)
    smb = df["sug_smb_units"].fillna(0).values.astype(float)
    has_dyn = df["has_dynisf"].fillna(0).values.astype(float)
    sens = df["sens"].fillna(np.nan).values.astype(float)

    ends = np.searchsorted(ts, ts + WINDOW_S)
    valid_end = (ends < n)
    end_idx = np.where(valid_end, ends, 0)
    valid_end &= (np.abs(ts[end_idx] - (ts + WINDOW_S)) <= 300)

    pre_start = np.searchsorted(ts, ts - PRE_S)
    valid_pre = (pre_start < np.arange(n)) & (np.abs(ts[np.clip(pre_start, 0, n-1)] - (ts - PRE_S)) <= 300)
    pre_idx = np.clip(pre_start, 0, n-1)

    diob = iob - iob[end_idx]
    dbg  = bg[end_idx] - bg

    # Rolling COB max (preceding 90 min)
    rolling_cob_max = np.zeros(n)
    j = 0
    for i in range(n):
        while j < i and ts[j] < ts[i] - COB_LEAD_S:
            j += 1
        rolling_cob_max[i] = cob[j:i + 1].max() if i >= j else 0.0

    # Pre-window SMB sum (preceding 30 min)
    psum_smb = np.concatenate([[0.0], np.cumsum(smb)])
    smb_starts = np.searchsorted(ts, ts - SMB_LEAD_S)
    rolling_smb_sum = psum_smb[np.arange(n)] - psum_smb[smb_starts]

    # In-window COB > 1
    cob_pos = (cob > COB_TOL).astype(np.int32)
    psum_cp = np.concatenate([[0], np.cumsum(cob_pos)])
    any_cob_in_window = psum_cp[end_idx] - psum_cp[np.arange(n)]

    pre_dt_min = (ts - ts[pre_idx]) / 60.0
    pre_dt_min_safe = np.where(pre_dt_min > 0, pre_dt_min, np.nan)
    bg_pre_trend = (bg - bg[pre_idx]) / pre_dt_min_safe

    keep = (
        valid_end & valid_pre
        & (any_cob_in_window == 0)
        & (rolling_cob_max <= COB_TOL)
        & (rolling_smb_sum <= SMB_TOL)
        & (bg >= CGM_LO) & (bg <= CGM_HI)
        & (bg[end_idx] >= 50) & (bg[end_idx] <= 300)
        & (np.abs(dbg) <= DBG_MAX)
        & (diob > DIOB_MIN) & (diob < DIOB_MAX)
        & np.isfinite(bg_pre_trend)
    )
    n_w = int(keep.sum())
    if n_w < MIN_WINDOWS:
        return {"user_id": user_id, "table": table,
                "skipped": "few_windows", "n_windows": n_w}

    x_diob  = diob[keep]
    x_trend = bg_pre_trend[keep]
    y       = dbg[keep]
    n_obs = len(y)
    X = np.column_stack([np.ones(n_obs), x_diob, x_trend])
    try:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    except np.linalg.LinAlgError:
        return {"user_id": user_id, "table": table, "skipped": "lstsq_failed"}

    pred = X @ beta
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # SE on the diob coefficient (per-user uncertainty)
    sigma2 = ss_res / max(n_obs - 3, 1)
    XtX_inv = np.linalg.pinv(X.T @ X)
    se_diob = float(np.sqrt(sigma2 * XtX_inv[1, 1]))

    # Empirical ISF = -slope_on_diob.  ΔBG is positive when BG rose, negative when it fell.
    # If ΔIOB > 0 (insulin absorbed) ⇒ BG should fall (ΔBG < 0) ⇒ slope is negative ⇒ ISF > 0.
    empirical_isf = float(-beta[1])
    se_isf = se_diob

    dyn_frac = float(has_dyn[keep].mean())
    dyn_on_mask = (has_dyn[keep] > 0)
    sens_when_on = float(np.nanmean(sens[keep][dyn_on_mask])) if dyn_on_mask.sum() > 30 else None

    return {
        "user_id": user_id, "table": table,
        "n_windows": n_w,
        "empirical_isf": empirical_isf,
        "se_isf": se_isf,
        "ci_low_isf":  float(empirical_isf - 1.96 * se_isf),
        "ci_high_isf": float(empirical_isf + 1.96 * se_isf),
        "intercept": float(beta[0]),
        "bg_trend_coef": float(beta[2]),
        "r2": float(r2),
        "dynisf_frac": dyn_frac,
        "mean_sens_when_dyn_on": sens_when_on,
    }


def main():
    from canonical_cohort import load_canonical_cohort
    cohort = load_canonical_cohort()
    incohort = cohort[cohort["in_cohort"]]
    work = []
    src_to_table = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}
    for _, r in incohort.iterrows():
        work.append((r["user_id"], src_to_table[r["cohort"]]))
    print(f"v5 (ΔIOB-based) inference: {len(work)} canonical users on {N_WORKERS} workers")

    with mp.Pool(N_WORKERS) as pool:
        results = pool.map(fit_user, work, chunksize=2)
    df = pd.DataFrame(results)
    df.to_json(OUT_JSON, orient="records", indent=2)

    df = df.merge(incohort[["user_id", "isf", "tdd", "cohort", "group"]],
                  on="user_id", how="left")
    df["entered_isf"] = df["isf"]
    df["ratio"] = df["empirical_isf"] / df["entered_isf"]
    valid = df[df["empirical_isf"].between(5, 500) & (df["r2"] >= 0.10)
               & df["entered_isf"].notna()]

    md = []
    md.append("# Empirical-ISF v5 — ΔIOB-based regression on canonical cohort\n")
    md.append("Replaces `iob_activity` (unit ambiguous) with `ΔIOB` (unambiguous: U "
              "absorbed in 30-min window).  Same fasting filter as v4. Per-user "
              "OLS: `ΔBG = a + b·ΔIOB + c·BG_pre_trend`. Empirical ISF = -b "
              "(mg/dL per U absorbed). Per-user 95 % CI from regression SE.\n")
    md.append(f"## Coverage\n")
    md.append(f"- Canonical cohort users in the analysis: {len(work)}")
    md.append(f"- Users with valid empirical_isf in [5,500] AND R² ≥ 0.10: **{len(valid)}**")
    if "skipped" in df.columns:
        for r, n_s in df["skipped"].dropna().value_counts().items():
            md.append(f"- skipped ({r}): {n_s}")
    md.append("")

    if len(valid):
        md.append("## Cohort-level results\n")
        md.append(f"- Median empirical ISF: **{valid['empirical_isf'].median():.0f}** "
                  f"[Q1 {valid['empirical_isf'].quantile(0.25):.0f}, "
                  f"Q3 {valid['empirical_isf'].quantile(0.75):.0f}] mg/dL/U")
        md.append(f"- Median entered ISF: **{valid['entered_isf'].median():.0f}** mg/dL/U")
        md.append(f"- Median empirical/entered ratio: **{valid['ratio'].median():.2f}** "
                  f"[Q1 {valid['ratio'].quantile(0.25):.2f}, "
                  f"Q3 {valid['ratio'].quantile(0.75):.2f}]")
        md.append(f"- Median R²: {valid['r2'].median():.2f}")
        md.append(f"- Median windows per user: {valid['n_windows'].median():.0f}")
        md.append("")

        md.append("## Stratified by DynISF group\n")
        md.append("| Group | n | Median entered | Median empirical | Median ratio | Median R² | Median per-user 95% CI half-width |")
        md.append("|---|---|---|---|---|---|---|")
        for g in ("no_dynisf", "dynisf_sigmoid", "dynisf_log"):
            s = valid[valid["group"] == g]
            if len(s) < 5:
                md.append(f"| **{g}** | {len(s)} | (n<5) | | | | |"); continue
            ci_hw = (s["ci_high_isf"] - s["ci_low_isf"]).median() / 2
            md.append(f"| **{g}** | {len(s)} | "
                      f"{s['entered_isf'].median():.0f} | "
                      f"{s['empirical_isf'].median():.0f} | "
                      f"{s['ratio'].median():.2f} | "
                      f"{s['r2'].median():.2f} | "
                      f"±{ci_hw:.0f} mg/dL/U |")
        md.append("")

        md.append("**Per-user uncertainty caveat**: median per-user 95% CI half-width "
                  "shows the regression noise per individual estimate. Cohort-level "
                  "medians are robust because errors cancel; individual user numbers "
                  "have wide CIs and should not be reported as 'this user's empirical ISF "
                  "is X' without the band.\n")

    OUT_MD.write_text("\n".join(md))
    df.to_json(OUT_JSON, orient="records", indent=2)
    print(f"Wrote {OUT_MD}")
    print("\n".join(md[-25:]))


if __name__ == "__main__":
    main()
