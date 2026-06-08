#!/usr/bin/env python3
"""Overnight insulin sensitivity at a full 4-hour horizon (per spec).

For each 5-minute CGM reading whose hour-of-day is in 23:00–02:55 (so the 4-hour measurement
ends by 07:00, before the dawn rise), in a period we are confident is fasting:

  1. Carb screen: across the T → T+4h window, require the 15-minute forward slope to stay
     ≤ 2 mg/dL per 5 min everywhere. If glucose rises faster than that at any point, assume
     carbs on board and discard the window.
  2. IOB(T): insulin on board at T (require ≥ IOB_MIN so the ratio is meaningful).
  3. Glucose drop over the next 4 h: ΔBG = BG(T) − BG(T+4h).
  4. Observed sensitivity at T = ΔBG / IOB(T)   (mg/dL per unit).

This lets the insulin on board fully act (4 h ≈ insulin duration) instead of the 30-min
ΔIOB window, which is why it is far less exposed to short-horizon glucose mean-reversion.

We report the distribution of observed sensitivity per person and across the population, and
sensitivity binned by the starting glucose BG(T) — overlaid with what each equation calculates
— to see whether, measured this way, sensitivity falls with glucose.

Output: results/overnight_sensitivity.{json,md},
        charts/inv008/fig_overnight_sensitivity.png
Run: python -m inv008.overnight_sensitivity
"""
from __future__ import annotations

import json
import math
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
from inv008.dynisf import quartic_isf

ROOT = config.ROOT
EMPIRICAL = ROOT / "empirical_isf_v5.json"
OUT = ROOT / "results"
CHART = ROOT / "charts" / "inv008"
DSN = "dbname=oref"
TARGET, DIV = config.NORMAL_TARGET, config.INSULIN_DIVISOR

START_HOURS = {23, 0, 1, 2}        # 11pm–3am start
HORIZON_S = 4 * 3600               # measure drop at T + 4h
TOL_S = 300                        # ±5 min tolerance on the horizon point
RISE_MAX = 2.0                     # mg/dL per 5 min; above this over 15 min ⇒ assume COB
IOB_MIN = 0.30                     # U; below this the drop/IOB ratio is too noisy
SENS_RANGE = (-50, 400)            # plausible mg/dL per U for the binned-shape view

BG_EDGES = np.array([99, 112, 125, 145, 175, 220], float)   # at/above target only
BG_CENTRES = 0.5 * (BG_EDGES[:-1] + BG_EDGES[1:])
TARGET_BIN = 0                                              # first band sits at target


def analyse(user_id, table):
    sql = f"""
        SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob
        FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL AND iob_iob IS NOT NULL
        ORDER BY ts_relative_sec
    """
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

    # index at T+4h (within tolerance) and at T+15min (for the rolling slope)
    end4 = np.searchsorted(ts, ts + HORIZON_S)
    end4 = np.where(end4 < n, end4, n - 1)
    ok_end = np.abs(ts[end4] - (ts + HORIZON_S)) <= TOL_S
    p15 = np.searchsorted(ts, ts + 900)
    p15 = np.where(p15 < n, p15, n - 1)
    ok15 = np.abs(ts[p15] - (ts + 900)) <= TOL_S
    slope15 = np.where(ok15, (bg[p15] - bg) / 3.0, 0.0)     # mg/dL per 5 min over next 15 min

    rows = []
    n_subtarget = 0
    for i in range(n):
        if hour[i] not in START_HOURS or not ok_end[i] or iob[i] < IOB_MIN:
            continue
        j = end4[i]
        # exclude sub-target starting points: below target the loop withholds insulin (and IOB
        # can go negative), so drop/IOB there does not reflect insulin sensitivity.
        if bg[i] < TARGET:
            n_subtarget += 1
            continue
        # carb screen: no 15-min forward slope above RISE_MAX anywhere in [T, T+4h]
        seg = slope15[i:j + 1]
        if seg.size == 0 or np.nanmax(seg) > RISE_MAX:
            continue
        drop = bg[i] - bg[j]                # positive = glucose fell
        sens = drop / iob[i]                # mg/dL per U
        rows.append((bg[i], iob[i], drop, sens))
    if len(rows) < 30:
        return None
    arr = np.array(rows)
    bg0, iob0, drop, sens = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
    binned = {}
    b = np.digitize(bg0, BG_EDGES) - 1
    for k in range(len(BG_CENTRES)):
        m = (b == k) & (sens > SENS_RANGE[0]) & (sens < SENS_RANGE[1])
        binned[k] = float(np.median(sens[m])) if m.sum() >= 10 else None
    return {
        "user": user_id, "n_windows": int(len(rows)),
        "median_bg_start": round(float(np.median(bg0)), 1),
        "median_iob": round(float(np.median(iob0)), 2),
        "median_drop": round(float(np.median(drop)), 1),
        "median_sens": round(float(np.median(sens)), 1),
        "sens_p25": round(float(np.percentile(sens, 25)), 1),
        "sens_p75": round(float(np.percentile(sens, 75)), 1),
        "binned_sens": binned,
        # keep a capped sample for pooled views
        "sens_sample": sens[(sens > SENS_RANGE[0]) & (sens < SENS_RANGE[1])].tolist(),
        "bg_sample": bg0[(sens > SENS_RANGE[0]) & (sens < SENS_RANGE[1])].tolist(),
    }


def main():
    emp = {e["user_id"]: e for e in json.load(open(EMPIRICAL))}
    jobs = [(u, e["table"]) for u, e in emp.items()]
    nw = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"overnight sensitivity: {len(jobs)} users on {nw} workers")
    with mp.Pool(nw) as pool:
        res = [r for r in pool.starmap(analyse, jobs) if r]
    n = len(res)

    per_med = np.array([r["median_sens"] for r in res])
    all_sens = np.concatenate([np.array(r["sens_sample"]) for r in res])
    all_bg = np.concatenate([np.array(r["bg_sample"]) for r in res])

    # population sensitivity vs BG(T): pooled median per band + the equation curves
    bsd = np.digitize(all_bg, BG_EDGES) - 1
    obs_band = []
    for k in range(len(BG_CENTRES)):
        m = bsd == k
        obs_band.append(float(np.median(all_sens[m])) if m.sum() >= 200 else None)
    ref = obs_band[TARGET_BIN]
    obs_norm = [None if v is None else v / ref for v in obs_band]

    def g_quartic_n(x): return float(quartic_isf(x) / quartic_isf(TARGET))
    def v1_s(x):
        xc = x if x <= 210 else 210 + (x - 210) / 3
        return math.log(TARGET / DIV + 1) / math.log(xc / DIV + 1)
    def v2_s(x): return math.log(TARGET / DIV) / math.log(max(x, DIV + 1) / DIV)
    gq = [g_quartic_n(c) for c in BG_CENTRES]
    g1 = [v1_s(c) for c in BG_CENTRES]
    g2 = [v2_s(c) for c in BG_CENTRES]

    okb = [i for i, v in enumerate(obs_norm) if v is not None and v > 0]
    k_obs = float(np.polyfit(np.log(TARGET / BG_CENTRES[okb]),
                             np.log([obs_norm[i] for i in okb]), 1)[0]) if len(okb) >= 3 else None

    summary = {
        "n_patients": n, "total_windows": int(sum(r["n_windows"] for r in res)),
        "method": "overnight 11pm-3am start, drop over T+4h / IOB(T), carb-screened",
        "population_median_sensitivity": round(float(np.median(all_sens)), 1),
        "population_p25_p75": [round(float(np.percentile(all_sens, 25)), 1),
                               round(float(np.percentile(all_sens, 75)), 1)],
        "per_person_median_sensitivity": {
            "median": round(float(np.median(per_med)), 1),
            "p25": round(float(np.percentile(per_med, 25)), 1),
            "p75": round(float(np.percentile(per_med, 75)), 1),
            "min": round(float(per_med.min()), 1), "max": round(float(per_med.max()), 1)},
        "sensitivity_vs_bg_pooled_median": {
            f"{BG_CENTRES[k]:.0f}": (None if obs_band[k] is None else round(obs_band[k], 1))
            for k in range(len(BG_CENTRES))},
        "sensitivity_vs_bg_normalised": {
            f"{BG_CENTRES[k]:.0f}": (None if obs_norm[k] is None else round(obs_norm[k], 3))
            for k in range(len(BG_CENTRES))},
        "calculated_normalised": {"vnext_quartic": [round(v, 3) for v in gq],
                                  "v1_log": [round(v, 3) for v in g1],
                                  "v2_log": [round(v, 3) for v in g2]},
        "observed_powerlaw_exponent": None if k_obs is None else round(k_obs, 2),
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    summary["per_person"] = [{k: r[k] for k in
                              ("user", "n_windows", "median_bg_start", "median_iob",
                               "median_drop", "median_sens", "sens_p25", "sens_p75")} for r in res]
    (OUT / "overnight_sensitivity.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].hist(per_med, bins=24, color="#2ca02c", alpha=0.85)
    ax[0].axvline(np.median(per_med), color="k", ls="--", lw=1.5,
                  label=f"median {np.median(per_med):.0f}")
    ax[0].set_xlabel("per-person median overnight sensitivity (mg/dL per U)")
    ax[0].set_ylabel("people"); ax[0].set_title(f"Observed overnight sensitivity ({n} people)")
    ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    xc = BG_CENTRES
    yv = np.array([v if v is not None else np.nan for v in obs_norm])
    ax[1].plot(xc, yv, "o-", color="#111", ms=7, lw=1.6, label="observed (overnight, normalised)", zorder=10)
    ax[1].plot(xc, gq, "-", color="#1f77b4", lw=2.2, label="v-next g(BG) quartic")
    ax[1].plot(xc, g1, "--", color="#888", lw=1.6, label="v1 log")
    ax[1].plot(xc, g2, ":", color="#d62728", lw=1.6, label="v2 log")
    ax[1].axhline(1, color="k", lw=0.6, alpha=0.4); ax[1].axvline(TARGET, color="k", lw=0.6, alpha=0.4)
    ax[1].set_xlabel("starting glucose BG(T) (mg/dL)"); ax[1].set_ylabel("sensitivity relative to target")
    ttl = f"Sensitivity vs glucose (overnight, 4 h horizon)"
    if k_obs is not None:
        ttl += f"\nobserved exponent k≈{k_obs:.1f}  (quartic +1.3)"
    ax[1].set_title(ttl); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    ax[1].set_ylim(0, max(2.0, np.nanmax(yv) * 1.15))
    fig.tight_layout(); fig.savefig(CHART / "fig_overnight_sensitivity.png", dpi=150); plt.close(fig)

    md = ["# Overnight insulin sensitivity at a 4-hour horizon\n",
          f"{n} people, {summary['total_windows']:,} carb-screened overnight windows "
          "(11pm–3am start). Sensitivity = (BG(T) − BG(T+4h)) / IOB(T), mg/dL per U.\n",
          "## Population\n",
          f"- Per-person median sensitivity: median **{summary['per_person_median_sensitivity']['median']}** "
          f"mg/dL/U [IQR {summary['per_person_median_sensitivity']['p25']}–"
          f"{summary['per_person_median_sensitivity']['p75']}, range "
          f"{summary['per_person_median_sensitivity']['min']}–{summary['per_person_median_sensitivity']['max']}].",
          f"- Pooled per-window: median {summary['population_median_sensitivity']} "
          f"[IQR {summary['population_p25_p75'][0]}–{summary['population_p25_p75'][1]}].",
          f"- Observed sensitivity-vs-glucose exponent: **k ≈ {summary['observed_powerlaw_exponent']}** "
          "(v-next quartic is +1.3; positive = falls with glucose).",
          "\n## Sensitivity vs starting glucose (normalised to target)\n",
          "| BG(T) | observed | v-next g(BG) | v1 | v2 |", "|---|---|---|---|---|"]
    for k in range(len(BG_CENTRES)):
        o = obs_norm[k]
        md.append(f"| {BG_CENTRES[k]:.0f} | {'–' if o is None else f'{o:.2f}'} | "
                  f"{gq[k]:.2f} | {g1[k]:.2f} | {g2[k]:.2f} |")
    md += ["\n![Overnight sensitivity](charts/inv008/fig_overnight_sensitivity.png)\n",
           "## Per-person (first 30 by median sensitivity)\n",
           "| user | n win | median BG(T) | median IOB | median drop | median sens | IQR |",
           "|---|---|---|---|---|---|---|"]
    for r in sorted(res, key=lambda x: x["median_sens"])[:30]:
        md.append(f"| {r['user']} | {r['n_windows']} | {r['median_bg_start']:.0f} | "
                  f"{r['median_iob']:.2f} | {r['median_drop']:.0f} | {r['median_sens']:.0f} | "
                  f"{r['sens_p25']:.0f}–{r['sens_p75']:.0f} |")
    md.append("\n*Caveat: basal continues over the 4 h and roughly offsets endogenous glucose in "
              "a fasting state, so drop/IOB(T) is an approximation of sensitivity; residual dawn "
              "effect, basal mis-set, and counter-regulation still bias the tails.*")
    (OUT / "overnight_sensitivity.md").write_text("\n".join(md))
    print("\n".join(md[:10]))
    print("observed k =", summary["observed_powerlaw_exponent"])


if __name__ == "__main__":
    main()
