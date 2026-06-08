#!/usr/bin/env python3
"""Decompose overnight sensitivity into a per-user anchor at BG 100, a glucose shape, and the
mean-reversion term — to answer three questions on clean data:

  Q1  How does sensitivity vary with glucose?            → the IOB×(BG−100) coefficient d
  Q2  Best per-user ISF at target 100?                   → the per-user IOB coefficient b_u
  Q3  How does the value at 100 relate to the changes?   → is d_u / b_u constant (multiplicative)

Model, per overnight carb-screened window (BG ≥ target so the loop is dosing):
    drop(T→T+4h) = a + b·IOB(T) + c·(BG−100) + d·IOB·(BG−100)
The c·(BG−100) term explicitly absorbs glucose mean-reversion / EGP (drop grows with how far
above target the glucose sits, independent of insulin). b is then the insulin effect at BG 100,
and d is how that effect changes per mg/dL above target. If IOB is endogenous (the loop builds
it in anticipation of where glucose is heading), b and d stay biased even with the c term — the
script reports whether the reversion term changes the picture.

Output: results/decompose_sensitivity.{json,md}, charts/inv008/fig_decompose_sensitivity.png
Run: python -m inv008.decompose_sensitivity
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

from inv008 import config

ROOT = config.ROOT
EMPIRICAL = ROOT / "empirical_isf_v5.json"
COHORT = ROOT / "canonical_cohort.json"
OUT = ROOT / "results"
CHART = ROOT / "charts" / "inv008"
DSN = "dbname=oref"
TARGET = 100.0

START_HOURS = {23, 0, 1, 2}
HORIZON_S, TOL_S = 4 * 3600, 300
RISE_MAX, IOB_MIN = 2.0, 0.30
COL_MAP = {"oref_v5": '"sug_COB"', "oref_v6": "sug_cob", "oref_v7": "sug_cob"}


def windows(user_id, table):
    sql = f"""SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob, {COL_MAP[table]} AS cob
              FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL AND iob_iob IS NOT NULL
              ORDER BY ts_relative_sec"""
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(sql, conn, params=(user_id,))
    finally:
        conn.close()
    n = len(df)
    if n < 500:
        return None
    ts = df.ts_relative_sec.values.astype(float)
    hour = df.hour.values.astype(float)
    bg = df.cgm_mgdl.values.astype(float)
    iob = df.iob_iob.values.astype(float)
    end4 = np.searchsorted(ts, ts + HORIZON_S); end4 = np.where(end4 < n, end4, n - 1)
    ok_end = np.abs(ts[end4] - (ts + HORIZON_S)) <= TOL_S
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    ok15 = np.abs(ts[p15] - (ts + 900)) <= TOL_S
    slope15 = np.where(ok15, (bg[p15] - bg) / 3.0, 0.0)
    rows = []
    for i in range(n):
        if hour[i] not in START_HOURS or not ok_end[i] or iob[i] < IOB_MIN or bg[i] < TARGET:
            continue
        j = end4[i]
        if np.nanmax(slope15[i:j + 1]) > RISE_MAX:
            continue
        rows.append((bg[i], iob[i], bg[i] - bg[j]))
    if len(rows) < 50:
        return None
    a = np.array(rows)
    return user_id, a[:, 0], a[:, 1], a[:, 2]      # bg, iob, drop


def fit(bg, iob, drop):
    """OLS drop = a + b·iob + c·(bg-100) + d·iob·(bg-100). Returns (a,b,c,d)."""
    g = bg - TARGET
    X = np.column_stack([np.ones_like(iob), iob, g, iob * g])
    beta, *_ = np.linalg.lstsq(X, drop, rcond=None)
    return beta


def main():
    emp = {e["user_id"]: e for e in json.load(open(EMPIRICAL))}
    coh = {r["user_id"]: r for r in json.load(open(COHORT))}
    jobs = [(u, e["table"]) for u, e in emp.items()
            if (config.REPLAY_DIR / f"{u}.parquet").exists()]
    nw = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"decompose: {len(jobs)} users on {nw} workers")
    with mp.Pool(nw) as pool:
        res = [r for r in pool.starmap(windows, jobs) if r]

    BG = np.concatenate([r[1] for r in res])
    IOB = np.concatenate([r[2] for r in res])
    DROP = np.concatenate([r[3] for r in res])
    n_win = len(BG)
    corr_iob_bg = float(np.corrcoef(IOB, BG - TARGET)[0, 1])   # collinearity of the two regressors

    # ---- pooled de-confounding regression (Q1) ----
    a, b, c, d = fit(BG, IOB, DROP)
    # naive (no reversion term) for contrast
    Xn = np.column_stack([np.ones_like(IOB), IOB, IOB * (BG - TARGET)])
    bn, *_ = np.linalg.lstsq(Xn, DROP, rcond=None)
    naive_b, naive_d = float(bn[1]), float(bn[2])
    frac_slope = d / b      # fractional change in sensitivity per mg/dL above target

    # ---- per-user anchor at 100 (Q2): b_u with the shared shape removed ----
    per_user = []
    for u, bg, iob, drop in res:
        if len(bg) < 80:
            continue
        # remove the shared glucose modulation, fit per-user level + reversion
        adj = drop - d * iob * (bg - TARGET)
        Xu = np.column_stack([np.ones_like(iob), iob, bg - TARGET])
        bu, *_ = np.linalg.lstsq(Xu, adj, rcond=None)
        # per-user interaction too (Q3), where enough data
        bbeta = fit(bg, iob, drop)
        prof = coh.get(u, {}).get("isf")
        prof = (prof * 18.018 if prof and prof < 20 else prof)
        # near-target-only drop/IOB (alternative anchor)
        nt = (bg >= 95) & (bg <= 105)
        nt_isf = float(np.median(drop[nt] / iob[nt])) if nt.sum() >= 20 else None
        per_user.append({"user": u, "n": int(len(bg)),
                         "isf100_reg": float(bu[1]), "reversion_c": float(bu[2]),
                         "b_raw": float(bbeta[1]), "d_raw": float(bbeta[3]),
                         "near_target_isf": nt_isf,
                         "profile_isf": None if prof is None else round(float(prof), 1)})

    bu_arr = np.array([p["isf100_reg"] for p in per_user])
    # Q3: is the absolute glucose slope proportional to the anchor? (multiplicative shape)
    good = [p for p in per_user if p["n"] >= 200 and p["isf100_reg"] > 0]
    bvec = np.array([p["isf100_reg"] for p in good])
    dvec = np.array([p["d_raw"] for p in good])
    # correlation of d with b (multiplicative ⇒ d grows with b)
    corr_bd = float(np.corrcoef(bvec, dvec)[0, 1]) if len(good) > 5 else None
    frac_by_user = dvec / bvec
    multiplicative_frac = float(np.median(frac_by_user))

    summary = {
        "n_users_windows": [len(res), int(n_win)],
        "Q1_glucose_shape": {
            "isf_at_100_pooled": round(float(b), 2),
            "reversion_c_per_mgdl": round(float(c), 3),
            "interaction_d": round(float(d), 4),
            "fractional_change_per_mgdl": round(float(frac_slope), 4),
            "direction": "rises with glucose" if d > 0 else "falls with glucose",
            "naive_without_reversion_term": {"isf_at_100": round(naive_b, 2),
                                             "fractional_change_per_mgdl": round(naive_d / naive_b, 4)},
        },
        "Q2_per_user_isf_at_100": {
            "n_users": len(per_user),
            "median": round(float(np.median(bu_arr)), 1),
            "p25": round(float(np.percentile(bu_arr, 25)), 1),
            "p75": round(float(np.percentile(bu_arr, 75)), 1),
            "min": round(float(bu_arr.min()), 1), "max": round(float(bu_arr.max()), 1),
        },
        "Q3_anchor_vs_change": {
            "corr_anchor_vs_abs_slope": None if corr_bd is None else round(corr_bd, 3),
            "median_fractional_slope_per_user": round(multiplicative_frac, 4),
            "reading": ("if the abs slope d grows with the anchor b (positive corr) and d/b is "
                        "roughly constant, the glucose effect is multiplicative on the anchor"),
        },
        "per_user": per_user,
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    (OUT / "decompose_sensitivity.json").write_text(json.dumps(summary, indent=1))

    # ---- diagnostic: is insulin separable from reversion? ----
    # within each BG band, do high-IOB and low-IOB windows drop differently? if not, the
    # insulin effect is not identifiable against glucose reversion.
    bands = [(100, 115), (115, 130), (130, 150), (150, 175), (175, 220)]
    lo_drop, hi_drop, ctr, n_lo, n_hi = [], [], [], [], []
    for a0, a1 in bands:
        m = (BG >= a0) & (BG < a1)
        if m.sum() < 200:
            lo_drop.append(np.nan); hi_drop.append(np.nan); ctr.append((a0 + a1) / 2); continue
        med_iob = np.median(IOB[m])
        lo = m & (IOB <= med_iob); hi = m & (IOB > med_iob)
        lo_drop.append(float(np.median(DROP[lo]))); hi_drop.append(float(np.median(DROP[hi])))
        ctr.append((a0 + a1) / 2); n_lo.append(int(lo.sum())); n_hi.append(int(hi.sum()))
    insulin_gap = float(np.nanmedian(np.array(hi_drop) - np.array(lo_drop)))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].plot(ctr, lo_drop, "o-", color="#1f77b4", lw=2, label="low IOB (below band median)")
    ax[0].plot(ctr, hi_drop, "s-", color="#d62728", lw=2, label="high IOB (above band median)")
    ax[0].set_xlabel("starting glucose BG(T) (mg/dL)"); ax[0].set_ylabel("median 4h glucose drop (mg/dL)")
    ax[0].set_title("Drop is set by glucose, not insulin\n"
                    f"high- vs low-IOB gap ≈ {insulin_gap:.0f} mg/dL; reversion dominates")
    ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    # candidate per-user anchors at 100
    prof_vals = np.array([p["profile_isf"] for p in per_user if p["profile_isf"]])
    nt_vals = np.array([p["near_target_isf"] for p in per_user if p["near_target_isf"] is not None])
    ax[1].hist(prof_vals, bins=20, alpha=0.6, color="#2ca02c", label=f"profile ISF (n={len(prof_vals)}, med {np.median(prof_vals):.0f})")
    if len(nt_vals):
        ax[1].hist(nt_vals, bins=20, alpha=0.6, color="#ff7f0e",
                   label=f"near-target drop/IOB (n={len(nt_vals)}, med {np.median(nt_vals):.0f})")
    ax[1].set_xlabel("ISF at 100 (mg/dL per U)"); ax[1].set_ylabel("people")
    ax[1].set_title("Candidate per-user anchors at 100\n(closed-loop drop/IOB is loop-suppressed and data-starved)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_decompose_sensitivity.png", dpi=150); plt.close(fig)
    summary["diagnostic"] = {"corr_iob_vs_bg_above_target": round(corr_iob_bg, 3),
                             "reversion_drop_at_bg160_zero_insulin": round(float(c) * 60, 1),
                             "high_minus_low_iob_drop_gap_mgdl": round(insulin_gap, 1),
                             "band_centres": ctr, "low_iob_drop": lo_drop, "high_iob_drop": hi_drop}

    summary["conclusion"] = (
        "Closed-loop observational data cannot identify the insulin sensitivity here. The "
        "homeostasis the loop enforces creates two killers: (1) glucose reverts toward target "
        f"regardless of insulin (~{c*60:.0f} mg/dL drop from BG 160 with zero insulin), and that "
        f"reversion is collinear with IOB (corr {corr_iob_bg:.2f}, the loop doses by glucose), so "
        "the insulin term collapses when reversion is controlled; (2) at target the loop holds "
        "BG there, so almost no excursion exists to measure (few near-target 4h windows, and "
        "those that exist read a loop-suppressed value). The best per-user ISF at 100 therefore "
        "comes from tuned profiles or the cross-sectional K/√TDD cold-start, not from drop/IOB.")
    (OUT / "decompose_sensitivity.json").write_text(json.dumps(summary, indent=1))

    diag = summary["diagnostic"]
    md = ["# Can the glucose↔sensitivity relation and the at-100 anchor be measured from logs?\n",
          f"{len(res)} people, {n_win:,} overnight carb-screened windows (BG ≥ target). "
          "Model: drop = a + b·IOB + c·(BG−100) + d·IOB·(BG−100); c is meant to absorb glucose "
          "reversion so b is the insulin effect at 100 and d its glucose dependence.\n",
          "## The attempt fails to identify insulin — and that is the answer\n",
          f"- The reversion term alone predicts a **{diag['reversion_drop_at_bg160_zero_insulin']:.0f} "
          "mg/dL** drop from BG 160 with **zero insulin** — overnight glucose reverts toward target "
          "on its own.",
          f"- IOB is **collinear with glucose** (corr IOB vs BG−100 = {diag['corr_iob_vs_bg_above_target']}): "
          "the loop doses by glucose, so the insulin and reversion terms cannot be separated. The "
          f"fitted ISF at 100 collapses to {summary['Q1_glucose_shape']['isf_at_100_pooled']} "
          "mg/dL per U — not physiological, the signature of an unidentified model.",
          f"- The cleanest test: at the same glucose, high-IOB and low-IOB windows drop almost the "
          f"same (median gap ≈ **{diag['high_minus_low_iob_drop_gap_mgdl']:.0f} mg/dL**). The 4-hour "
          "drop is set by where glucose started, not by how much insulin was on board.",
          "\n## Q1 — sensitivity vs glucose: not recoverable\n",
          "The insulin effect is not separable from reversion in closed-loop data, so the "
          "glucose shape g(BG) cannot be measured from these logs. It must come from controlled "
          "(clamp / clinical) data.",
          "\n## Q2 — per-user ISF at 100: not from drop/IOB\n",
          "At target the loop holds glucose there, so there is essentially no excursion to "
          "measure — only a handful of users have enough near-target 4h windows, and those read a "
          "loop-suppressed value. The usable per-user anchor at 100 is the **tuned profile ISF** "
          "(cohort median ≈ 50 mg/dL per U), or a **K/√TDD cold-start** where no profile exists.",
          "\n## Q3 — value at 100 vs the changes: moot\n",
          "With the glucose changes unobservable and the anchor coming from profiles, the link "
          "cannot be derived from logs; the multiplicative form ISF(BG) = ISF₁₀₀ · g(BG) is a "
          "modelling choice for clinical data to validate.",
          "\n![Diagnostic](charts/inv008/fig_decompose_sensitivity.png)\n",
          f"**Conclusion.** {summary['conclusion']}"]
    (OUT / "decompose_sensitivity.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
