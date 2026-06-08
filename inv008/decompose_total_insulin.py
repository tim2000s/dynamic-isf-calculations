#!/usr/bin/env python3
"""Re-do the overnight sensitivity decomposition using TOTAL delivered insulin.

iob_iob is net of scheduled basal, so it misses the scheduled-basal insulin that (we showed)
drives the overnight glucose fall. Here the dose is the actual insulin delivered over the
window — scheduled basal + temp overlay + boluses — taken from the stage-1 delivery grid
(total_u per 5-min bin, which already reflects temp reductions and suspensions). CGM is aligned
to the grid via the same recovered anchor stage 2 used.

Model per overnight carb-screened window (v6/v7 only — these have the reconstructed grid):
    drop(T→T+4h) = a + b·INS + c·(BG−100) + d·INS·(BG−100)
where INS = total insulin delivered over the window.
    b   = sensitivity at BG 100 (mg/dL per U)         → Q2 anchor
    −a  = EGP rise over 4h with zero insulin           (should be positive: glucose rises)
    c   = residual glucose term; if total insulin truly explains the fall, c shrinks toward 0
    d   = how sensitivity changes per mg/dL above target → Q1 shape

Output: results/decompose_total_insulin.{json,md}, charts/inv008/fig_decompose_total_insulin.png
Run: python -m inv008.decompose_total_insulin
"""
from __future__ import annotations

import json
import multiprocessing as mp
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")
from inv008 import config, sources
from inv008.tdd_windows import build_delivery_grid
from inv008.stage2_replay import choose_anchor

ROOT = config.ROOT
OUT = ROOT / "results"; CHART = ROOT / "charts" / "inv008"
TARGET = 100.0
START_HOURS = {23, 0, 1, 2}
HZ, TOL = 4 * 3600, 300
RISE_MAX = 2.0
TBL = {"v6": "oref_v6", "v7": "oref_v7"}


def analyse(args):
    user_id, platform, raw_dir = args
    if platform not in TBL:
        return None
    hourly = sources.load_hourly_basal(user_id)
    if hourly is None:
        return None
    try:
        if platform == "v7":
            bolus_ts, bolus_u, temps = sources.load_v7_events(raw_dir)
        else:
            bolus_ts, bolus_u, temps = sources.load_v6_events(raw_dir)
    except Exception:
        return None
    if len(bolus_ts) < 50:
        return None
    t0, t1 = float(bolus_ts.min()), float(bolus_ts.max())
    if temps:
        t0 = min(t0, temps[0].ts); t1 = max(t1, temps[-1].ts)
    anchor_start = sources.recover_anchor(platform, raw_dir)
    if anchor_start:
        t0 = min(t0, anchor_start)
    grid = build_delivery_grid(int(t0), int(t1), bolus_ts, bolus_u, temps, hourly)
    gts = grid["ts"].to_numpy(float)
    gcum = np.concatenate([[0.0], np.cumsum(grid["total_u"].to_numpy(float))])

    tcache = config.TDD_DIR / f"{user_id}.parquet"
    if not tcache.exists():
        return None
    win = pd.read_parquet(tcache)
    conn = psycopg2.connect("dbname=oref")
    try:
        df = pd.read_sql(f"""SELECT ts_relative_sec, cgm_mgdl, hour, iob_iob
                             FROM {TBL[platform]} WHERE user_id=%s AND cgm_mgdl IS NOT NULL
                             ORDER BY ts_relative_sec""", conn, params=(user_id,))
    finally:
        conn.close()
    if len(df) < 500:
        return None
    try:
        _, ameta = choose_anchor(df, win, anchor_start)
    except Exception:
        return None
    anchor = ameta["anchor_epoch_sec"]

    ts = df.ts_relative_sec.values.astype(float)
    bg = df.cgm_mgdl.values.astype(float)
    hr = df.hour.values.astype(float)
    n = len(df)
    abs_ts = anchor + ts
    end4 = np.searchsorted(ts, ts + HZ); end4 = np.where(end4 < n, end4, n - 1)
    ok_end = np.abs(ts[end4] - (ts + HZ)) <= TOL
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    sl = np.where(np.abs(ts[p15] - (ts + 900)) <= TOL, (bg[p15] - bg) / 3.0, 0.0)

    rows = []
    for i in range(n):
        if hr[i] not in START_HOURS or not ok_end[i] or bg[i] < 80 or bg[i] > 260:
            continue
        j = end4[i]
        if np.nanmax(sl[i:j + 1]) > RISE_MAX:
            continue
        a0 = abs_ts[i]; a1 = a0 + HZ
        if a1 > gts[-1] or a0 < gts[0]:
            continue
        ins = float(gcum[np.searchsorted(gts, a1)] - gcum[np.searchsorted(gts, a0)])
        if ins <= 0:
            continue
        rows.append((bg[i], ins, bg[i] - bg[j]))
    if len(rows) < 80:
        return None
    arr = np.array(rows)
    return user_id, arr[:, 0], arr[:, 1], arr[:, 2]


def reg(bg, ins, drop):
    g = bg - TARGET
    X = np.column_stack([np.ones_like(ins), ins, g, ins * g])
    beta, *_ = np.linalg.lstsq(X, drop, rcond=None)
    return beta  # a, b(ISF@100), c, d


def main():
    maps = sources.load_mappings()
    jobs = [(u, p, r) for u, (p, r) in maps.items() if p in TBL]
    nw = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"total-insulin decomposition: {len(jobs)} v6/v7 users on {nw} workers")
    with mp.Pool(nw) as pool:
        res = [r for r in pool.map(analyse, jobs) if r]

    BG = np.concatenate([r[1] for r in res])
    INS = np.concatenate([r[2] for r in res])
    DROP = np.concatenate([r[3] for r in res])
    a, b, c, d = reg(BG, INS, DROP)
    egp_h = -a / 4.0

    # diagnostic: at the same glucose, do high-insulin windows drop more now?
    bands = [(90, 110), (110, 130), (130, 155), (155, 200)]
    lo_d, hi_d, ctr = [], [], []
    for a0, a1 in bands:
        m = (BG >= a0) & (BG < a1)
        if m.sum() < 200:
            lo_d.append(np.nan); hi_d.append(np.nan); ctr.append((a0+a1)/2); continue
        mi = np.median(INS[m]); lo = m & (INS <= mi); hi = m & (INS > mi)
        lo_d.append(float(np.median(DROP[lo]))); hi_d.append(float(np.median(DROP[hi]))); ctr.append((a0+a1)/2)
    ins_gap = float(np.nanmedian(np.array(hi_d) - np.array(lo_d)))

    # per-user ISF at 100
    per_user = []
    for u, bg, ins, drop in res:
        if len(bg) < 120:
            continue
        adj = drop - d * ins * (bg - TARGET)
        Xu = np.column_stack([np.ones_like(ins), ins, bg - TARGET])
        bu, *_ = np.linalg.lstsq(Xu, adj, rcond=None)
        per_user.append({"user": u, "n": int(len(bg)), "isf100": float(bu[1]),
                         "egp_h": float(-bu[0] / 4.0)})
    bu_arr = np.array([p["isf100"] for p in per_user])

    summary = {
        "n_users": len(res), "n_windows": int(len(BG)),
        "pooled": {"isf_at_100": round(float(b), 1),
                   "egp_rise_mgdl_per_h_zero_insulin": round(float(egp_h), 1),
                   "residual_c_per_mgdl": round(float(c), 3),
                   "interaction_d": round(float(d), 3),
                   "frac_change_per_mgdl": round(float(d / b), 4) if b else None},
        "high_minus_low_insulin_drop_gap_mgdl": round(ins_gap, 1),
        "per_user_isf_at_100": {"n": len(per_user), "median": round(float(np.median(bu_arr)), 1),
                                "p25": round(float(np.percentile(bu_arr, 25)), 1),
                                "p75": round(float(np.percentile(bu_arr, 75)), 1),
                                "frac_positive": round(float(np.mean(bu_arr > 0)), 3)},
        "diagnostic_bands": {"centres": ctr, "low_ins_drop": lo_d, "high_ins_drop": hi_d},
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    summary["per_user"] = per_user
    (OUT / "decompose_total_insulin.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].plot(ctr, lo_d, "o-", color="#1f77b4", lw=2, label="low total insulin")
    ax[0].plot(ctr, hi_d, "s-", color="#d62728", lw=2, label="high total insulin")
    ax[0].set_xlabel("starting glucose (mg/dL)"); ax[0].set_ylabel("median 4h drop (mg/dL)")
    ax[0].set_title(f"Now does insulin separate?\nhigh− low-insulin drop gap ≈ {ins_gap:.0f} mg/dL")
    ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    ax[1].hist(bu_arr[(bu_arr > -20) & (bu_arr < 120)], bins=26, color="#1f77b4", alpha=0.85)
    ax[1].axvline(np.median(bu_arr), color="k", ls="--", lw=1.5, label=f"median {np.median(bu_arr):.0f}")
    ax[1].set_xlabel("per-user ISF at 100 (mg/dL per U, total-insulin model)"); ax[1].set_ylabel("people")
    ax[1].set_title(f"Q2: per-user ISF at 100 ({len(per_user)} people)\nEGP rise {egp_h:.0f} mg/dL/h at zero insulin")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_decompose_total_insulin.png", dpi=150); plt.close(fig)

    p = summary["pooled"]; q2 = summary["per_user_isf_at_100"]
    md = ["# Overnight sensitivity on TOTAL delivered insulin (basal + temp + bolus)\n",
          f"{len(res)} v6/v7 users, {len(BG):,} overnight carb-screened windows. Dose = insulin "
          "delivered over the window from the reconstructed grid (includes scheduled basal). "
          "Model: drop = a + b·INS + c·(BG−100) + d·INS·(BG−100).\n",
          "## Result\n",
          f"- **ISF at 100 (pooled): {p['isf_at_100']} mg/dL per U** — "
          f"{'physiological' if 10 < p['isf_at_100'] < 120 else 'still off'}.",
          f"- EGP rise with zero insulin: **{p['egp_rise_mgdl_per_h_zero_insulin']} mg/dL/h** "
          "(glucose rises without insulin — as it should).",
          f"- Residual glucose term c = {p['residual_c_per_mgdl']} (was the 'reversion'; smaller "
          "now that scheduled basal is counted).",
          f"- Glucose interaction d = {p['interaction_d']} → {p['frac_change_per_mgdl']*100 if p['frac_change_per_mgdl'] else 0:.2f}%/mg/dL "
          f"({'rises' if p['interaction_d']>0 else 'falls'} with glucose).",
          f"- At the same glucose, high- vs low-insulin windows now differ by "
          f"**{summary['high_minus_low_insulin_drop_gap_mgdl']:.0f} mg/dL** (was ~7 with iob_iob) "
          "— insulin is starting to separate from the glucose level.",
          "\n## Q2 — per-user ISF at 100\n",
          f"- Median **{q2['median']}** mg/dL per U [IQR {q2['p25']}–{q2['p75']}], "
          f"{100*q2['frac_positive']:.0f}% positive ({q2['n']} people).",
          "\n![Total-insulin decomposition](charts/inv008/fig_decompose_total_insulin.png)\n",
          "*v6/v7 only (they have the reconstructed delivery grid; v5/Trio excluded). The grid's "
          "total_u reflects actual delivery including temp reductions and suspensions. CGM aligned "
          "to the grid by the stage-2 recovered anchor.*"]
    (OUT / "decompose_total_insulin.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
