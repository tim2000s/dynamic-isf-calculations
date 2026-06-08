#!/usr/bin/env python3
"""Does the ISF the loop uses match the actual effect of the insulin it delivered?

This does not try to measure ISF from scratch (we showed that is not identifiable in a closed
loop). Instead it audits the loop's own forecast: each cycle the loop predicts where glucose
will settle from IOB, using the ISF it chose (reason_iobpredbg). We compare that prediction to
the realised glucose at full insulin action, in fasting overnight windows, and ask where the
prediction is biased — i.e. where the loop's ISF over- or under-states what insulin actually did.

Per overnight carb-screened window (start hour 23–02):
    pred_end   = reason_iobpredbg(T)      # loop's IOB-only predicted eventual BG (uses its ISF)
    actual_end = CGM(T + 4h)
    error      = actual_end − pred_end
        error > 0  → ended higher than predicted → insulin did LESS than the loop's ISF expected
                     → ISF too aggressive (too low)
        error < 0  → ended lower  than predicted → insulin did MORE → ISF too weak (too high)

We break the error down by glucose at decision, ISF formula (no_dynisf / sigmoid / log), TDD
band, and the loop's own autosens ratio at the time — to see where and why ISF diverges.

Output: results/calibration_audit.{json,md}, charts/inv008/fig_calibration_audit.png
Run: python -m inv008.calibration_audit
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
TBL = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}
# column names differ in case between tables (v5 quoted mixed-case)
COLS = {
    "oref_v5": ('"reason_IOBpredBG"', '"sug_ISF"', '"sug_sensitivityRatio"'),
    "oref_v6": ("reason_iobpredbg", "sug_isf", "sug_sensitivityratio"),
    "oref_v7": ("reason_iobpredbg", "sug_isf", "sug_sensitivityratio"),
}
START_HOURS = {23, 0, 1, 2}
HZ, TOL, RISE_MAX = 4 * 3600, 300, 2.0
BG_BANDS = [(80, 100), (100, 120), (120, 145), (145, 175), (175, 230)]


def analyse_user(args):
    user_id, table, group, tdd = args
    predc, isfc, src = COLS[table]
    conn = psycopg2.connect("dbname=oref")
    try:
        df = pd.read_sql(f"""SELECT ts_relative_sec, hour, cgm_mgdl, {predc} AS pred,
                                    {isfc} AS isf, {src} AS sr
                             FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL
                               AND {predc} IS NOT NULL ORDER BY ts_relative_sec""",
                         conn, params=(user_id,))
    finally:
        conn.close()
    n = len(df)
    if n < 500:
        return None
    ts = df.ts_relative_sec.values.astype(float); hr = df.hour.values.astype(float)
    bg = df.cgm_mgdl.values.astype(float)
    pred = pd.to_numeric(df["pred"], errors="coerce").values.astype(float)
    isf = pd.to_numeric(df["isf"], errors="coerce").values.astype(float)
    sr = pd.to_numeric(df["sr"], errors="coerce").values.astype(float)
    end4 = np.searchsorted(ts, ts + HZ); end4 = np.where(end4 < n, end4, n - 1)
    oke = np.abs(ts[end4] - (ts + HZ)) <= TOL
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    sl = np.where(np.abs(ts[p15] - (ts + 900)) <= TOL, (bg[p15] - bg) / 3.0, 0.0)
    rows = []
    for i in range(n):
        if hr[i] not in START_HOURS or not oke[i] or not np.isfinite(pred[i]) or bg[i] < 80 or bg[i] > 260:
            continue
        j = end4[i]
        if np.nanmax(sl[i:j + 1]) > RISE_MAX:
            continue
        rows.append((bg[i], pred[i], bg[j], isf[i], sr[i]))
    if len(rows) < 40:
        return None
    a = np.array(rows)
    bg0, pred_, end_, isf_, sr_ = a.T
    err = end_ - pred_
    return {"user": user_id, "group": group, "tdd": tdd, "n": int(len(err)),
            "median_error": float(np.median(err)),
            "median_isf": float(np.nanmedian(isf_)),
            "median_sr": float(np.nanmedian(sr_)),
            "bg": bg0.tolist(), "err": err.tolist(), "sr": sr_.tolist()}


def main():
    coh = json.load(open(ROOT / "canonical_cohort.json"))
    jobs = [(r["user_id"], TBL[r["cohort"]], r.get("group", "?"), float(r.get("tdd", 0)))
            for r in coh if r.get("cohort") in TBL and r.get("in_cohort", True)]
    nw = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"calibration audit: {len(jobs)} users on {nw} workers")
    with mp.Pool(nw) as pool:
        res = [r for r in pool.map(analyse_user, jobs) if r]

    BG = np.concatenate([np.array(r["bg"]) for r in res])
    ERR = np.concatenate([np.array(r["err"]) for r in res])
    SR = np.concatenate([np.array(r["sr"]) for r in res])
    grp = np.concatenate([[r["group"]] * r["n"] for r in res])

    def band_stats(mask):
        e = ERR[mask]
        return {"n": int(mask.sum()), "median_error": round(float(np.median(e)), 1),
                "p25": round(float(np.percentile(e, 25)), 1),
                "p75": round(float(np.percentile(e, 75)), 1)} if mask.sum() >= 100 else None

    by_bg = {f"{a}-{b}": band_stats((BG >= a) & (BG < b)) for a, b in BG_BANDS}
    by_group = {g: band_stats(grp == g) for g in ("no_dynisf", "dynisf_sigmoid", "dynisf_log")}
    # per-user median error vs TDD band
    tdds = np.array([r["tdd"] for r in res]); permed = np.array([r["median_error"] for r in res])
    by_tdd = {}
    for lab, lo, hi in [("<25", 0, 25), ("25-45", 25, 45), ("45-70", 45, 70), ("70+", 70, 1e9)]:
        m = (tdds >= lo) & (tdds < hi)
        by_tdd[lab] = round(float(np.median(permed[m])), 1) if m.sum() >= 3 else None
    # autosens: when the loop set SR != 1, did it help? correlation of error with (sr-1)
    oksr = np.isfinite(SR)
    corr_err_sr = float(np.corrcoef(ERR[oksr], SR[oksr] - 1)[0, 1])

    summary = {
        "n_users": len(res), "n_windows": int(len(ERR)),
        "overall_median_error_mgdl": round(float(np.median(ERR)), 1),
        "interpretation": ("error = actual_end − loop_predicted_end. >0: insulin did less than the "
                           "loop's ISF expected (ISF too aggressive); <0: did more (ISF too weak)."),
        "by_glucose_at_decision": by_bg,
        "by_isf_formula": by_group,
        "by_tdd_band_per_user": by_tdd,
        "corr_error_with_autosens_minus_1": round(corr_err_sr, 3),
        "median_autosens_ratio": round(float(np.nanmedian(SR)), 3),
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    summary["per_user"] = [{k: r[k] for k in ("user", "group", "tdd", "n", "median_error",
                                              "median_isf", "median_sr")} for r in res]
    (OUT / "calibration_audit.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ctr = [(a + b) / 2 for a, b in BG_BANDS]
    med = [by_bg[f"{a}-{b}"]["median_error"] if by_bg[f"{a}-{b}"] else np.nan for a, b in BG_BANDS]
    p25 = [by_bg[f"{a}-{b}"]["p25"] if by_bg[f"{a}-{b}"] else np.nan for a, b in BG_BANDS]
    p75 = [by_bg[f"{a}-{b}"]["p75"] if by_bg[f"{a}-{b}"] else np.nan for a, b in BG_BANDS]
    ax[0].plot(ctr, med, "o-", color="#1f77b4", lw=2.2)
    ax[0].fill_between(ctr, p25, p75, alpha=0.15, color="#1f77b4")
    ax[0].axhline(0, color="k", lw=1, ls="--", label="perfectly calibrated")
    ax[0].set_xlabel("glucose at decision (mg/dL)")
    ax[0].set_ylabel("actual − predicted end BG (mg/dL)")
    ax[0].set_title("Loop ISF calibration vs glucose\n>0: insulin did less than ISF predicted (too aggressive)")
    ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    gl = [g for g in ("no_dynisf", "dynisf_sigmoid", "dynisf_log") if by_group[g]]
    vals = [by_group[g]["median_error"] for g in gl]
    ax[1].bar(range(len(gl)), vals, color=["#7f7f7f", "#1f77b4", "#2ca02c"][:len(gl)], alpha=0.85)
    ax[1].axhline(0, color="k", lw=1)
    ax[1].set_xticks(range(len(gl))); ax[1].set_xticklabels([g.replace("dynisf_", "") for g in gl])
    ax[1].set_ylabel("median actual − predicted end BG (mg/dL)")
    ax[1].set_title("Calibration by ISF formula")
    ax[1].grid(alpha=0.3, axis="y")
    fig.tight_layout(); fig.savefig(CHART / "fig_calibration_audit.png", dpi=150); plt.close(fig)

    s = summary
    md = ["# Does the loop's ISF match the actual effect of insulin? (calibration audit)\n",
          f"{len(res)} users, {len(ERR):,} fasting overnight windows. error = actual end BG − the "
          "loop's IOB-predicted end BG (which uses the ISF it chose). >0 means insulin did less "
          "than the loop expected (ISF too aggressive); <0 means it did more (ISF too weak).\n",
          f"## Overall: median error **{s['overall_median_error_mgdl']} mg/dL**\n",
          "## By glucose at decision\n", "| BG band | n | median error | IQR |", "|---|---|---|---|"]
    for a, b in BG_BANDS:
        st = by_bg[f"{a}-{b}"]
        md.append(f"| {a}-{b} | {st['n']:,} | {st['median_error']} | {st['p25']}–{st['p75']} |" if st else f"| {a}-{b} | – | – | – |")
    md += ["\n## By ISF formula\n", "| formula | n | median error | IQR |", "|---|---|---|---|"]
    for g in ("no_dynisf", "dynisf_sigmoid", "dynisf_log"):
        st = by_group[g]
        md.append(f"| {g} | {st['n']:,} | {st['median_error']} | {st['p25']}–{st['p75']} |" if st else f"| {g} | – | – | – |")
    md += ["\n## By TDD band (per-user median error)\n",
           "| TDD | median error |", "|---|---|"]
    for k, v in by_tdd.items():
        md.append(f"| {k} | {v} |")
    md += [f"\n- Loop's own autosens ratio: median {s['median_autosens_ratio']}; correlation of "
           f"prediction error with (autosens−1) = {s['corr_error_with_autosens_minus_1']} "
           "(if autosens were fully correcting, residual error would be ~0 and uncorrelated).",
           "\n![Calibration audit](charts/inv008/fig_calibration_audit.png)\n",
           "*Fasting overnight, carb-screened, full-action (4h) horizon. The loop keeps dosing "
           "over the window, so the error is the net calibration of the loop's ISF-based forecast "
           "against the realised outcome, not a pure open-loop insulin response.*"]
    (OUT / "calibration_audit.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
