#!/usr/bin/env python3
"""Quantify the glucose double-count in the Tier-2 anchor.

The Tier-2 anchor is the measured ISF — the slope of ΔBG on ΔIOB across a person's 30-min
fasting windows (empirical_isf_v5). That regression has no glucose-level term, so its slope is
a ΔIOB²-weighted average of the true ISF over whatever BG the windows sat at:

    measured_ISF  ≈  anchor_at_target · E_w[g(BG_window)]

Tier-2 then uses measured_ISF as the at-target anchor and multiplies by g(BG) again, so its
ISF is biased by the factor E_w[g(BG_window)] at every glucose. This script re-extracts the
exact fasting windows, records each window's BG, and computes that mean-g factor per person —
giving a corrected at-target anchor = measured_ISF / E_w[g], and the size of the bias.

Windows use the identical definition as empirical_isf_v5.py. Representative window glucose is
the mean of the start and end CGM (the BG over which the insulin acted).

Output: results/tier2_anchor_debias.{json,md}, charts/inv008/fig_tier2_anchor_debias.png
Run: python -m inv008.tier2_anchor_debias
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")

from inv008 import config
from inv008.dynisf import g_quartic

ROOT = config.ROOT
EMPIRICAL = ROOT / "empirical_isf_v5.json"
COHORT = ROOT / "canonical_cohort.json"
OUT = ROOT / "results"
CHART = ROOT / "charts" / "inv008"
DSN = "dbname=oref"
TARGET = config.NORMAL_TARGET
R2_MIN = 0.10

# window constants (identical to empirical_isf_v5.py)
WINDOW_S, PRE_S, COB_LEAD_S, SMB_LEAD_S = 1800, 1800, 5400, 1800
COB_TOL, SMB_TOL = 1.0, 0.05
CGM_LO, CGM_HI, DBG_MAX, DIOB_MIN, DIOB_MAX = 70, 250, 80.0, 0.0, 2.0
COL_MAP = {
    "oref_v5": {"cob": '"sug_COB"'}, "oref_v6": {"cob": "sug_cob"}, "oref_v7": {"cob": "sug_cob"},
}


def window_bg_diob(user_id, table):
    """Return (bg_window, diob) for the kept fasting windows — same filter as empirical_isf_v5."""
    cob_col = COL_MAP[table]["cob"]
    sql = f"""
        SELECT ts_relative_sec, cgm_mgdl, iob_iob, {cob_col} AS cob, sug_smb_units
        FROM {table} WHERE user_id = %s AND cgm_mgdl IS NOT NULL AND iob_iob IS NOT NULL
        ORDER BY ts_relative_sec
    """
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(sql, conn, params=(user_id,))
    finally:
        conn.close()
    n = len(df)
    if n < 200:
        return None
    ts = df["ts_relative_sec"].values.astype(float)
    bg = df["cgm_mgdl"].values.astype(float)
    iob = df["iob_iob"].values.astype(float)
    cob = df["cob"].fillna(0).values.astype(float)
    smb = df["sug_smb_units"].fillna(0).values.astype(float)

    ends = np.searchsorted(ts, ts + WINDOW_S)
    valid_end = ends < n
    end_idx = np.where(valid_end, ends, 0)
    valid_end &= np.abs(ts[end_idx] - (ts + WINDOW_S)) <= 300
    pre_start = np.searchsorted(ts, ts - PRE_S)
    valid_pre = (pre_start < np.arange(n)) & (np.abs(ts[np.clip(pre_start, 0, n-1)] - (ts - PRE_S)) <= 300)

    diob = iob - iob[end_idx]
    dbg = bg[end_idx] - bg

    rolling_cob_max = np.zeros(n); j = 0
    for i in range(n):
        while j < i and ts[j] < ts[i] - COB_LEAD_S:
            j += 1
        rolling_cob_max[i] = cob[j:i + 1].max() if i >= j else 0.0
    psum_smb = np.concatenate([[0.0], np.cumsum(smb)])
    rolling_smb_sum = psum_smb[np.arange(n)] - psum_smb[np.searchsorted(ts, ts - SMB_LEAD_S)]
    cob_pos = (cob > COB_TOL).astype(np.int32)
    psum_cp = np.concatenate([[0], np.cumsum(cob_pos)])
    any_cob_in_window = psum_cp[end_idx] - psum_cp[np.arange(n)]
    pre_idx = np.clip(pre_start, 0, n - 1)
    pre_dt_min = (ts - ts[pre_idx]) / 60.0
    bg_pre_trend = (bg - bg[pre_idx]) / np.where(pre_dt_min > 0, pre_dt_min, np.nan)

    keep = (valid_end & valid_pre & (any_cob_in_window == 0) & (rolling_cob_max <= COB_TOL)
            & (rolling_smb_sum <= SMB_TOL) & (bg >= CGM_LO) & (bg <= CGM_HI)
            & (bg[end_idx] >= 50) & (bg[end_idx] <= 300) & (np.abs(dbg) <= DBG_MAX)
            & (diob > DIOB_MIN) & (diob < DIOB_MAX) & np.isfinite(bg_pre_trend))
    if keep.sum() < 80:
        return None
    bg_w = 0.5 * (bg[keep] + bg[end_idx][keep])     # mean BG over the window
    return bg_w, diob[keep]


def analyse(args):
    user_id, table, measured_isf, profile_isf = args
    out = window_bg_diob(user_id, table)
    if out is None:
        return None
    bg_w, diob = out
    g = np.asarray(g_quartic(bg_w))
    w = diob ** 2                                    # OLS implicit weight on ISF(BG_w)
    mean_g_w = float(np.sum(w * g) / np.sum(w))      # ΔIOB²-weighted (matches the regression)
    mean_g_simple = float(np.mean(g))
    corrected_anchor = measured_isf / mean_g_w
    return {
        "user": user_id, "n_windows": int(len(bg_w)),
        "median_window_bg": round(float(np.median(bg_w)), 1),
        "frac_windows_above_target": round(float(np.mean(bg_w > TARGET)), 3),
        "measured_isf": round(measured_isf, 1), "profile_isf": round(profile_isf, 1),
        "mean_g_weighted": round(mean_g_w, 3), "mean_g_simple": round(mean_g_simple, 3),
        "corrected_anchor": round(corrected_anchor, 1),
    }


def to_mgdl(v):
    v = float(v); return v * 18.018 if v < 20 else v


def main():
    emp = {e["user_id"]: e for e in json.load(open(EMPIRICAL))}
    coh = {r["user_id"]: r for r in json.load(open(COHORT))}
    jobs = []
    for u, e in emp.items():
        if (u in coh and e.get("r2", 0) >= R2_MIN and 5 <= e.get("empirical_isf", 0) <= 500
                and (config.REPLAY_DIR / f"{u}.parquet").exists()):
            jobs.append((u, e["table"], float(e["empirical_isf"]), to_mgdl(coh[u]["isf"])))
    n_workers = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"Tier-2 anchor de-bias: {len(jobs)} users on {n_workers} workers")
    with mp.Pool(n_workers) as pool:
        res = [r for r in pool.map(analyse, jobs) if r]
    n = len(res)

    mg = np.array([r["mean_g_weighted"] for r in res])
    wbg = np.array([r["median_window_bg"] for r in res])
    shift = np.array([r["corrected_anchor"] / r["measured_isf"] for r in res])  # = 1/mean_g

    summary = {
        "n_patients": n,
        "median_window_bg": round(float(np.median(wbg)), 1),
        "median_frac_windows_above_target": round(float(np.median(
            [r["frac_windows_above_target"] for r in res])), 3),
        "mean_g_factor": {"median": round(float(np.median(mg)), 3),
                          "p25": round(float(np.percentile(mg, 25)), 3),
                          "p75": round(float(np.percentile(mg, 75)), 3),
                          "min": round(float(mg.min()), 3), "max": round(float(mg.max()), 3)},
        "interpretation": ("Tier-2 ISF is biased by the mean_g factor at every glucose: "
                           ">1 means the current anchor makes Tier-2 too weak (corrected anchor "
                           "lower / more aggressive); <1 means too aggressive."),
        "anchor_shift_corrected_over_current": {
            "median": round(float(np.median(shift)), 3),
            "p25": round(float(np.percentile(shift, 25)), 3),
            "p75": round(float(np.percentile(shift, 75)), 3)},
        "per_person": res,
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    (OUT / "tier2_anchor_debias.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].hist(wbg, bins=24, color="#6a6ad6", alpha=0.85)
    ax[0].axvline(TARGET, color="r", ls="--", lw=1.5, label=f"target {TARGET:.0f}")
    ax[0].axvline(np.median(wbg), color="k", ls=":", lw=1.5, label=f"median {np.median(wbg):.0f}")
    ax[0].set_xlabel("per-person median fasting-window BG (mg/dL)")
    ax[0].set_ylabel("people"); ax[0].set_title("Where the measured-ISF windows sit")
    ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    ax[1].hist(mg, bins=24, color="#d62728", alpha=0.8)
    ax[1].axvline(1.0, color="k", ls="--", lw=1.5, label="1.0 (no bias)")
    ax[1].axvline(np.median(mg), color="b", ls=":", lw=1.5, label=f"median {np.median(mg):.2f}")
    ax[1].set_xlabel("mean g(BG) over the windows  =  Tier-2 glucose double-count factor")
    ax[1].set_ylabel("people")
    ax[1].set_title("How far the Tier-2 anchor is off\n(ISF biased by this factor at every BG)")
    ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_tier2_anchor_debias.png", dpi=150); plt.close(fig)

    mgm = summary["mean_g_factor"]["median"]
    md = ["# Tier-2 anchor — quantifying the glucose double-count\n",
          f"{n} people. The measured ISF is a ΔIOB²-weighted average of the true ISF over each "
          "person's fasting windows; with no glucose-level term in the regression, it carries the "
          "average g(BG) of those windows. Tier-2 then re-applies g(BG), so its ISF is off by the "
          "mean-g factor at every glucose.\n",
          "## Result\n",
          f"- Fasting windows sit at a median BG of **{summary['median_window_bg']:.0f} mg/dL** "
          f"(target {TARGET:.0f}); a median **{100*summary['median_frac_windows_above_target']:.0f}%** "
          "of each person's windows are above target.",
          f"- Mean g(BG) over the windows (the bias factor): median **{mgm}** "
          f"[IQR {mgm and summary['mean_g_factor']['p25']}–{summary['mean_g_factor']['p75']}, "
          f"range {summary['mean_g_factor']['min']}–{summary['mean_g_factor']['max']}].",
          f"- So the corrected at-target anchor = measured_ISF / {mgm} → a median "
          f"**{summary['anchor_shift_corrected_over_current']['median']}×** the current anchor "
          f"(IQR {summary['anchor_shift_corrected_over_current']['p25']}–"
          f"{summary['anchor_shift_corrected_over_current']['p75']}×).",
          "\n![Tier-2 anchor de-bias](charts/inv008/fig_tier2_anchor_debias.png)\n",
          "## Reading\n",
          (f"The double-count is real but modest: a median factor of {mgm}. Because fasting "
           "windows sit close to target for most people, g(BG) over them is near 1, so the bias "
           "is small on average — though it runs both ways across the cohort "
           f"({summary['mean_g_factor']['min']}–{summary['mean_g_factor']['max']}) depending on "
           "whether a person's fasting glucose ran above or below target. Where windows sat above "
           "target (g<1) the current anchor is too low and Tier-2 over-doses; where below (g>1) "
           "it under-doses. The clean fix is to divide the measured slope by this mean-g factor "
           "before using it as the at-target anchor, or to fit the slope on near-target windows "
           "only."),
          "\n## Per-person\n",
          "| person | n win | median win BG | % win >target | measured ISF | mean g | corrected anchor |",
          "|---|---|---|---|---|---|---|"]
    for r in sorted(res, key=lambda x: x["mean_g_weighted"]):
        md.append(f"| {r['user']} | {r['n_windows']} | {r['median_window_bg']:.0f} | "
                  f"{100*r['frac_windows_above_target']:.0f}% | {r['measured_isf']:.0f} | "
                  f"{r['mean_g_weighted']:.2f} | {r['corrected_anchor']:.0f} |")
    (OUT / "tier2_anchor_debias.md").write_text("\n".join(md))
    print("\n".join(md[:12]))


if __name__ == "__main__":
    main()
