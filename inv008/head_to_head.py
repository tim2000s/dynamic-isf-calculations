#!/usr/bin/env python3
"""Same-window head-to-head: which ISF predicts the realised drop best — static or dynamic?

The loop's IOB-based BG prediction is linear in ISF (predicted drop = ISF × activity-integral),
so on any window we can take the loop's own prediction (made with the ISF it ran, sug_isf) and
rescale it to any candidate ISF, then compare to the actual outcome — on the identical window.
This removes the between-user confound: every ISF form is tested on every window.

Per overnight carb-screened window (full-action 4h horizon), candidate ISFs (all in mg/dL/U):
    static    = the user's tuned profile ISF (units-cleaned, constant)
    v1        = dynamic 1800/(TDD·ln(BG/75+1))           (from the replay cache)
    v2        = dynamic 115000/(TDD²·ln(BG_floored/75))
    loop      = what the loop actually used (sug_isf) — the calibration baseline
predicted_end(ISF) = BG − (BG − reason_iobpredbg)·(ISF / sug_isf)
error(ISF)         = actual_end − predicted_end(ISF)

Also derive the per-window realised ISF (the ISF that would have made the prediction exact):
    realised_isf = sug_isf · (BG − actual_end) / (BG − reason_iobpredbg)
saved with features (bg, tdd, iob, hour) as the dataset for the ML pattern step.

Units: per user, if median sug_isf < 20 it is mmol/L per U → ×18.018 to mg/dL (same for profile).
v6 is excluded (its iob/isf accounting did not reconcile). Output:
results/head_to_head.{json,md}, results/head_to_head_windows.parquet,
charts/inv008/fig_head_to_head.png
Run: python -m inv008.head_to_head
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
OUT = ROOT / "results"; CHART = ROOT / "charts" / "inv008"
TBL = {"v5_trio": "oref_v5", "v7_oref0": "oref_v7"}    # v6 excluded (iob/isf reconcile issue)
PREDC = {"oref_v5": '"reason_IOBpredBG"', "oref_v7": "reason_iobpredbg"}
ISFC = {"oref_v5": '"sug_ISF"', "oref_v7": "sug_isf"}
START_HOURS = {23, 0, 1, 2}
HZ, TOL, RISE_MAX = 4 * 3600, 300, 2.0
MIN_PRED_DROP = 8.0     # only windows where the loop predicted a real drop (ratio stable)
BG_BANDS = [(80, 100), (100, 120), (120, 145), (145, 175), (175, 230)]


def to_mgdl_scalar(v):
    return v * 18.018 if (v is not None and v < 20) else v


def analyse(args):
    user_id, table, profile_isf = args
    rp = config.REPLAY_DIR / f"{user_id}.parquet"
    if not rp.exists():
        return None
    rep = pd.read_parquet(rp, columns=["ts_relative_sec", "bg", "tdd", "isf_v1", "isf_v2"])
    conn = psycopg2.connect("dbname=oref")
    try:
        db = pd.read_sql(f"""SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob,
                                    {PREDC[table]} AS pred, {ISFC[table]} AS sugisf
                             FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL
                               AND {PREDC[table]} IS NOT NULL ORDER BY ts_relative_sec""",
                         conn, params=(user_id,))
    finally:
        conn.close()
    if len(db) < 500:
        return None
    d = rep.merge(db, on="ts_relative_sec", how="inner")
    if len(d) < 500:
        return None
    d = d.sort_values("ts_relative_sec").reset_index(drop=True)
    ts = d.ts_relative_sec.values.astype(float)
    bg = d.bg.values.astype(float)
    hr = d.hour.values.astype(float)
    pred = pd.to_numeric(d["pred"], errors="coerce").values.astype(float)
    sugisf = pd.to_numeric(d["sugisf"], errors="coerce").values.astype(float)
    iob = pd.to_numeric(d["iob_iob"], errors="coerce").values.astype(float)
    v1 = d.isf_v1.values.astype(float); v2 = d.isf_v2.values.astype(float)
    tdd = d.tdd.values.astype(float)

    # unit cleaning: scale sug_isf (and profile) to mg/dL if mmol-scale
    med_sug = float(np.nanmedian(sugisf[sugisf > 0])) if np.any(sugisf > 0) else np.nan
    scale = 18.018 if (np.isfinite(med_sug) and med_sug < 20) else 1.0
    sugisf = sugisf * scale
    pisf = to_mgdl_scalar(profile_isf)

    n = len(d)
    end4 = np.searchsorted(ts, ts + HZ); end4 = np.where(end4 < n, end4, n - 1)
    oke = np.abs(ts[end4] - (ts + HZ)) <= TOL
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    sl = np.where(np.abs(ts[p15] - (ts + 900)) <= TOL, (bg[p15] - bg) / 3.0, 0.0)

    rows = []
    for i in range(n):
        if hr[i] not in START_HOURS or not oke[i] or bg[i] < 80 or bg[i] > 260:
            continue
        if not (np.isfinite(pred[i]) and np.isfinite(sugisf[i]) and sugisf[i] > 0):
            continue
        pdl = bg[i] - pred[i]                  # loop predicted drop (with sug_isf)
        if pdl < MIN_PRED_DROP:
            continue
        j = end4[i]
        if np.nanmax(sl[i:j + 1]) > RISE_MAX:
            continue
        act_end = bg[j]
        act_drop = bg[i] - act_end
        ai = pdl / sugisf[i]                   # activity-integral (ISF-independent)
        # predicted end under each candidate ISF
        def end_for(isf):
            return bg[i] - ai * isf
        e_static = act_end - end_for(pisf) if pisf else np.nan
        e_v1 = act_end - end_for(v1[i]) if np.isfinite(v1[i]) else np.nan
        e_v2 = act_end - end_for(v2[i]) if np.isfinite(v2[i]) else np.nan
        e_loop = act_end - pred[i]
        realised = sugisf[i] * act_drop / pdl
        rows.append((bg[i], tdd[i], iob[i], hr[i], sugisf[i], pisf or np.nan,
                     v1[i], v2[i], realised, e_static, e_v1, e_v2, e_loop))
    if len(rows) < 40:
        return None
    cols = ["bg", "tdd", "iob", "hour", "sug_isf", "profile_isf", "isf_v1", "isf_v2",
            "realised_isf", "err_static", "err_v1", "err_v2", "err_loop"]
    df = pd.DataFrame(rows, columns=cols); df["user"] = user_id; df["table"] = table
    return df


def main():
    coh = {r["user_id"]: r for r in json.load(open(ROOT / "canonical_cohort.json"))}
    jobs = [(u, TBL[r["cohort"]], r.get("isf")) for u, r in coh.items()
            if r.get("cohort") in TBL and r.get("isf")]
    nw = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"head-to-head: {len(jobs)} v5/v7 users on {nw} workers")
    with mp.Pool(nw) as pool:
        parts = [p for p in pool.map(analyse, jobs) if p is not None]
    D = pd.concat(parts, ignore_index=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    D.to_parquet(OUT / "head_to_head_windows.parquet", index=False)

    forms = {"static": "err_static", "v1": "err_v1", "v2": "err_v2", "loop": "err_loop"}

    def mae_bias(col, mask=None):
        e = D[col] if mask is None else D[col][mask]
        e = e.dropna()
        return {"mae": round(float(e.abs().median()), 1), "bias": round(float(e.median()), 1),
                "n": int(len(e))}

    overall = {f: mae_bias(c) for f, c in forms.items()}
    by_bg = {}
    for a, b in BG_BANDS:
        m = (D.bg >= a) & (D.bg < b)
        by_bg[f"{a}-{b}"] = {f: mae_bias(c, m) for f, c in forms.items()}
    # per-user: which form wins (lowest MAE)?
    win = {f: 0 for f in forms}
    nuser = 0
    for u, g in D.groupby("user"):
        if len(g) < 60:
            continue
        nuser += 1
        maes = {f: g[c].abs().median() for f, c in forms.items()}
        win[min(maes, key=maes.get)] += 1

    summary = {
        "n_users": int(D.user.nunique()), "n_windows": int(len(D)),
        "overall_mae_bias": overall,
        "by_bg_band": by_bg,
        "per_user_best_form_counts": win, "n_users_scored": nuser,
        "note": "error = actual_end − predicted_end(ISF); MAE & bias in mg/dL. Lower MAE = better "
                "predictor of the realised outcome on the same windows.",
    }
    (OUT / "head_to_head.json").write_text(json.dumps(summary, indent=1))

    # figure: MAE by form overall + bias by BG band
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    fl = ["static", "v1", "v2", "loop"]
    ax[0].bar(fl, [overall[f]["mae"] for f in fl],
              color=["#2ca02c", "#888", "#d62728", "#1f77b4"])
    ax[0].set_ylabel("median |error| (mg/dL)"); ax[0].set_title("Prediction error by ISF form (same windows)")
    ax[0].grid(alpha=0.3, axis="y")
    ctr = [(a + b) / 2 for a, b in BG_BANDS]
    for f, col in [("static", "#2ca02c"), ("v1", "#888"), ("v2", "#d62728")]:
        ax[1].plot(ctr, [by_bg[f"{a}-{b}"][f]["bias"] for a, b in BG_BANDS], "o-", color=col, lw=2, label=f)
    ax[1].axhline(0, color="k", lw=1, ls="--")
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("median error / bias (mg/dL)")
    ax[1].set_title("Bias vs glucose by ISF form\n(>0: ISF over-predicts the drop)")
    ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_head_to_head.png", dpi=150); plt.close(fig)

    md = ["# Same-window head-to-head: static vs dynamic ISF as a drop predictor\n",
          f"{summary['n_users']} v5/v7 users, {len(D):,} overnight carb-screened windows. Each "
          "ISF form is tested on the *same* windows by rescaling the loop's IOB prediction. "
          "error = actual end BG − predicted end BG; lower MAE = better.\n",
          "## Overall (median |error|, mg/dL)\n",
          "| form | MAE | bias | n |", "|---|---|---|---|"]
    for f in fl:
        o = overall[f]
        md.append(f"| {f} | {o['mae']} | {o['bias']} | {o['n']:,} |")
    md += ["\n## Bias by glucose (mg/dL; >0 = over-predicts the drop)\n",
           "| BG band | static | v1 | v2 | loop |", "|---|---|---|---|---|"]
    for a, b in BG_BANDS:
        r = by_bg[f"{a}-{b}"]
        md.append(f"| {a}-{b} | {r['static']['bias']} | {r['v1']['bias']} | {r['v2']['bias']} | {r['loop']['bias']} |")
    md += [f"\n## Per-user best predictor (lowest MAE), {nuser} users\n",
           "| form | users where best |", "|---|---|"]
    for f in fl:
        md.append(f"| {f} | {win[f]} |")
    md += ["\n![Head-to-head](charts/inv008/fig_head_to_head.png)\n",
           "*Per-window dataset (features + realised ISF) saved to "
           "`results/head_to_head_windows.parquet` for the pattern/ML step. v6 excluded "
           "(iob/isf accounting did not reconcile). Units cleaned (mmol-scale ISF ×18.018).*"]
    (OUT / "head_to_head.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
