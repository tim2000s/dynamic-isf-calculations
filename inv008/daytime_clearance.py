#!/usr/bin/env python3
"""Phase C: does the clearance↔resistance cancellation hold by DAY, or does net resistance emerge?

Overnight (23:00-02:00), insulin resistance at high BG is real but insulin-independent clearance rises
with BG and cancels it, leaving the NET realised-ISF flat. The power law steepened in DAYTIME data
(k→4) and daytime errors are ~2× overnight (activity/stress/hormones). So the open question for a
glucose term is the daytime regime: if clearance still offsets resistance, the net stays flat and no
glucose term is warranted; if it does NOT, the net realised ISF falls with BG by day and a (gentle)
daytime glucose term is justified.

Same method as `clearance_corrected_isf` but daytime start hours (09:00-16:00, carb-screened so meal
windows are excluded — the rise screen drops windows where glucose climbs). Compares daytime vs the
saved overnight curves. The decisive metric is the RAW (net) ratio by BG: flat ⇒ clearance still
offsets ⇒ no glucose term; falling ⇒ net daytime resistance ⇒ a glucose term is warranted by day.

Caveat: daytime between-meal fasting windows are fewer and more contaminated (post-meal decay can
leak past a rise-only screen). Single-process. Output: results/daytime_clearance.{json,md},
charts/inv008/fig_daytime_clearance.png. Run: python -m inv008.daytime_clearance
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config
from inv008.clearance_corrected_isf import (BANDS, LBL, CTR, PREDC, TBL, band_of,
                                            HZ, TOL, RISE_MAX, POWERLAW_RATIO_HI)
import psycopg2

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
DAY_HOURS = {9, 10, 11, 12, 13, 14, 15, 16}


def user_windows(user_id, table, hours):
    conn = psycopg2.connect("dbname=oref")
    try:
        d = pd.read_sql(f"""SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob,
                                   {PREDC[table]} AS pred
                            FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL
                              AND {PREDC[table]} IS NOT NULL ORDER BY ts_relative_sec""",
                        conn, params=(user_id,))
    finally:
        conn.close()
    if len(d) < 500:
        return None
    ts = d.ts_relative_sec.values.astype(float); bg = d.cgm_mgdl.values.astype(float)
    hr = d.hour.values.astype(float)
    pred = pd.to_numeric(d.pred, errors="coerce").values.astype(float)
    n = len(d)
    end4 = np.searchsorted(ts, ts + HZ); end4 = np.where(end4 < n, end4, n - 1)
    oke = np.abs(ts[end4] - (ts + HZ)) <= TOL
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    sl = np.where(np.abs(ts[p15] - (ts + 900)) <= TOL, (bg[p15] - bg) / 3.0, 0.0)
    rows = []
    for i in range(n):
        if hr[i] not in hours or not oke[i] or bg[i] < 80 or bg[i] > 260 or not np.isfinite(pred[i]):
            continue
        j = end4[i]
        if np.nanmax(sl[i:j + 1]) > RISE_MAX:
            continue
        rows.append((bg[i], bg[i] - bg[j], bg[i] - pred[i]))
    if len(rows) < 30:
        return None
    df = pd.DataFrame(rows, columns=["bg", "drop", "pred_drop"]); df["user"] = user_id
    return df


def analyse(D):
    bidx = band_of(D.bg.values)
    low = np.abs(D.pred_drop.values) < 5.0
    nonins = {}
    for kk in range(len(BANDS)):
        m = low & (bidx == kk)
        nonins[kk] = float(np.median(D["drop"].values[m])) if m.sum() >= 40 else np.nan
    act = D[D.pred_drop >= 8].copy()
    ab = band_of(act.bg.values)
    act["raw"] = (act["drop"] / act.pred_drop).clip(-1, 3)
    act["corr"] = ((act["drop"] - np.array([nonins.get(b, np.nan) for b in ab])) / act.pred_drop).clip(-1, 3)

    def curve(col):
        a = act; ab = band_of(a.bg.values)
        return [round(float(a[ab == kk].groupby("user")[col].median().median()), 2)
                if (ab == kk).sum() >= 100 and a[ab == kk].user.nunique() >= 6 else None
                for kk in range(len(BANDS))]
    return {"nonins": [round(nonins[k], 1) if np.isfinite(nonins[k]) else None for k in range(len(BANDS))],
            "raw": curve("raw"), "corr": curve("corr"), "n": int(len(act))}


def slope_k(c):
    pts = [(CTR[i], c[i]) for i in range(len(c)) if c[i] and CTR[i] >= 120 and c[i] > 0.02]
    if len(pts) < 3:
        return None
    x = np.log(np.array([p[0] for p in pts]) / 100.0); y = np.log([p[1] for p in pts])
    return round(-float(np.polyfit(x, y, 1)[0]), 2)


def main():
    coh = {r["user_id"]: r for r in json.load(open(config.ROOT / "canonical_cohort.json"))}
    jobs = [(u, TBL[r["cohort"]]) for u, r in coh.items()
            if r.get("cohort") in TBL and r.get("isf")]
    print(f"daytime: {len(jobs)} users, single-process")
    parts = []
    for k, (u, table) in enumerate(jobs):
        df = user_windows(u, table, DAY_HOURS)
        if df is not None:
            parts.append(df)
        if k % 25 == 0:
            print(f"  {k}/{len(jobs)}")
    D = pd.concat(parts, ignore_index=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    day = analyse(D)

    # overnight comparison (from the saved test #1 run)
    ov = json.load(open(OUT / "clearance_corrected_isf.json"))
    ov_raw = [ov["raw_ratio_by_bg"].get(l) for l in LBL]
    ov_corr = [ov["clearance_corrected_ratio_by_bg"].get(l) for l in LBL]
    ov_non = [ov["data_derived_nonInsulin_flux_mgdl_by_bg"].get(l) for l in LBL]

    k_day_raw, k_day_corr = slope_k(day["raw"]), slope_k(day["corr"])
    k_ov_raw, k_ov_corr = slope_k(ov_raw), slope_k(ov_corr)
    net_resistance_day = bool(k_day_raw and k_day_raw > 0.2)

    summary = {
        "day_hours": sorted(DAY_HOURS), "n_users": int(D.user.nunique()), "n_day_active": day["n"],
        "daytime": {"nonInsulin_flux": dict(zip(LBL, day["nonins"])),
                    "raw_net_ratio": dict(zip(LBL, day["raw"])),
                    "clearance_corrected_ratio": dict(zip(LBL, day["corr"]))},
        "overnight_ref": {"raw_net_ratio": dict(zip(LBL, ov_raw)),
                          "clearance_corrected_ratio": dict(zip(LBL, ov_corr)),
                          "nonInsulin_flux": dict(zip(LBL, ov_non))},
        "implied_k": {"day_raw_net": k_day_raw, "day_corr_insulin": k_day_corr,
                      "overnight_raw_net": k_ov_raw, "overnight_corr_insulin": k_ov_corr},
        "verdict": ("DAYTIME NET RESISTANCE — raw net ratio falls with BG by day (clearance does NOT "
                    "fully offset); a gentle daytime glucose term is warranted" if net_resistance_day
                    else "clearance still offsets by day — net ratio flat; NO daytime glucose term "
                    "warranted (same as overnight)"),
        "caveat": "daytime between-meal fasting windows are fewer and more carb-contaminated (post-meal "
                  "decay can leak past a rise-only screen); treat as indicative.",
    }
    (OUT / "daytime_clearance.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].plot(CTR, day["raw"], "o-", color="#d62728", lw=2.5, label=f"DAY net (k={k_day_raw})")
    ax[0].plot(CTR, ov_raw, "s--", color="#1f77b4", lw=2, label=f"overnight net (k={k_ov_raw})")
    ax[0].axhline(1, color="k", ls="--", lw=1)
    ax[0].set_xlabel("glucose (mg/dL)"); ax[0].set_ylabel("raw (net) realised ÷ profile ISF")
    ax[0].set_title("NET realised ISF: day vs overnight\n(day falling ⇒ glucose term warranted)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(CTR, day["nonins"], "o-", color="#9467bd", lw=2, label="day non-insulin flux")
    ax[1].plot(CTR, ov_non, "s--", color="#888", lw=2, label="overnight non-insulin flux")
    ax[1].axvline(180, color="#888", ls=":", lw=1)
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("non-insulin flux (mg/dL/4h)")
    ax[1].set_title("Clearance: day vs overnight"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_daytime_clearance.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Phase C: does the clearance↔resistance cancellation hold by day?\n",
          f"{summary['n_users']} users, {day['n']:,} daytime insulin-active windows "
          f"(start hours {sorted(DAY_HOURS)}, carb-screened). Compared to the saved overnight run.\n",
          "## NET realised ÷ profile ISF by glucose (the decisive metric)\n",
          "| BG band | daytime net | overnight net |", "|---|---|---|"]
    for i in range(len(BANDS)):
        md.append(f"| {LBL[i]} | {day['raw'][i]} | {ov_raw[i]} |")
    md += [f"\nImplied k (net): **day {k_day_raw}**, overnight {k_ov_raw} (k>0 = net resistance).\n",
           "## Clearance-corrected (insulin-only) ratio\n", "| BG band | daytime | overnight |",
           "|---|---|---|"]
    for i in range(len(BANDS)):
        md.append(f"| {LBL[i]} | {day['corr'][i]} | {ov_corr[i]} |")
    md += [f"\nImplied k (insulin-only): day {k_day_corr}, overnight {k_ov_corr}.\n",
           "## Non-insulin clearance flux (mg/dL/4h)\n", "| BG band | daytime | overnight |", "|---|---|---|"]
    for i in range(len(BANDS)):
        md.append(f"| {LBL[i]} | {day['nonins'][i]} | {ov_non[i]} |")
    md += [f"\n**Verdict: {summary['verdict']}.**\n",
           "![daytime](charts/inv008/fig_daytime_clearance.png)\n", "*" + summary["caveat"] + "*"]
    (OUT / "daytime_clearance.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
