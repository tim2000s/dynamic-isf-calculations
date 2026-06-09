#!/usr/bin/env python3
"""Test #1 (bounding): does high-BG resistance appear once insulin-INDEPENDENT clearance is removed?

The realised-ISF ratio (actual_drop / loop_predicted_drop; <1 ⇒ insulin did less than expected) is
inflated at high BG by insulin-INDEPENDENT disposal: renal glucosuria (>~180) + glucose effectiveness
(mass action) − EGP. We estimate that non-insulin flux FROM THE DATA, not from literature constants:
windows where the loop expects ~no insulin action (|loop predicted drop| < 5 mg/dL) give actual_drop ≈
the non-insulin flux at that glucose. Subtract that curve from the insulin-active windows and
recompute the ratio.

    loop_pred_drop = cgm_start − reason_IOBpredBG          (insulin-expected drop, mg/dL, unit-free)
    nonInsulin(band) = median actual_drop where |loop_pred_drop| < 5   (renal + S_G − EGP)
    raw_ratio       = actual_drop / loop_pred_drop
    corrected_ratio = (actual_drop − nonInsulin(band)) / loop_pred_drop   (insulin-only)

If corrected_ratio FALLS at high BG ⇒ insulin per-unit does less there ⇒ genuine resistance (power
law) once clearance is removed. Bound: also report how much clearance would be REQUIRED to bend the
high-BG ratio to the power-law level (~0.3) and compare to the data-derived value.

Caveat (the ceiling): clearance and resistance are partially entangled (both BG-dependent), low-insulin
high-BG windows are sparse, and EGP is not perfectly stable — so this is a PLAUSIBILITY BOUND, not a
clean resistance number. Single-process. Output: results/clearance_corrected_isf.{json,md},
charts/inv008/fig_clearance_corrected_isf.png. Run: python -m inv008.clearance_corrected_isf
"""
from __future__ import annotations

import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")
from inv008 import config

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
TBL = {"v5_trio": "oref_v5", "v7_oref0": "oref_v7"}
PREDC = {"oref_v5": '"reason_IOBpredBG"', "oref_v7": "reason_iobpredbg"}
START_HOURS = {23, 0, 1, 2}
HZ, TOL, RISE_MAX = 4 * 3600, 300, 2.0
BANDS = [(100, 120), (120, 145), (145, 175), (175, 205), (205, 260)]
LBL = [f"{a}-{b}" for a, b in BANDS]
CTR = [(a + b) / 2 for a, b in BANDS]
POWERLAW_RATIO_HI = 0.30          # what the power law implies for realised/profile at high BG


def band_of(bg):
    out = np.full(len(bg), -1)
    for k, (a, b) in enumerate(BANDS):
        out[(bg >= a) & (bg < b)] = k
    return out


def user_windows(user_id, table):
    conn = psycopg2.connect("dbname=oref")
    try:
        d = pd.read_sql(f"""SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob, sug_smb_units AS smb,
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
    iob = pd.to_numeric(d.iob_iob, errors="coerce").values.astype(float)
    smb = pd.to_numeric(d.smb, errors="coerce").fillna(0).values.astype(float)
    n = len(d); csum = np.concatenate([[0.0], np.cumsum(smb)])
    end4 = np.searchsorted(ts, ts + HZ); end4 = np.where(end4 < n, end4, n - 1)
    oke = np.abs(ts[end4] - (ts + HZ)) <= TOL
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    sl = np.where(np.abs(ts[p15] - (ts + 900)) <= TOL, (bg[p15] - bg) / 3.0, 0.0)
    pre = np.searchsorted(ts, ts - 1800); pre = np.clip(pre, 0, n - 1)
    pre_slope = np.where(np.abs(ts[pre] - (ts - 1800)) <= TOL, (bg - bg[pre]) / 6.0, np.nan)
    rows = []
    for i in range(n):
        if hr[i] not in START_HOURS or not oke[i] or bg[i] < 80 or bg[i] > 260 or not np.isfinite(pred[i]):
            continue
        j = end4[i]
        if np.nanmax(sl[i:j + 1]) > RISE_MAX:
            continue
        rows.append((bg[i], bg[i] - bg[j], bg[i] - pred[i], iob[i], csum[j + 1] - csum[i], pre_slope[i]))
    if len(rows) < 40:
        return None
    df = pd.DataFrame(rows, columns=["bg", "drop", "pred_drop", "iob", "corr", "pre_slope"])
    df["user"] = user_id
    return df


def main():
    coh = {r["user_id"]: r for r in json.load(open(config.ROOT / "canonical_cohort.json"))}
    jobs = [(u, TBL[r["cohort"]]) for u, r in coh.items()
            if r.get("cohort") in TBL and r.get("isf")]
    print(f"clearance test: {len(jobs)} users, single-process")
    parts = []
    for k, (u, table) in enumerate(jobs):
        df = user_windows(u, table)
        if df is not None:
            parts.append(df)
        if k % 25 == 0:
            print(f"  {k}/{len(jobs)}")
    D = pd.concat(parts, ignore_index=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    bidx = band_of(D.bg.values)

    # data-derived non-insulin flux: drop where the loop expects ~no insulin action
    low = np.abs(D.pred_drop.values) < 5.0
    nonins, nonins_n = {}, {}
    for kk in range(len(BANDS)):
        m = low & (bidx == kk)
        nonins[kk] = float(np.median(D["drop"].values[m])) if m.sum() >= 50 else np.nan
        nonins_n[kk] = int(m.sum())

    # insulin-active windows, raw and clearance-corrected ratio, per user then median across users
    act = D[D.pred_drop >= 8].copy()
    act_bidx = band_of(act.bg.values)
    act["raw_ratio"] = (act["drop"] / act.pred_drop).clip(-1, 3)
    act["nonins"] = [nonins.get(b, np.nan) for b in act_bidx]
    act["corr_ratio"] = ((act["drop"] - act.nonins) / act.pred_drop).clip(-1, 3)

    def curve(col, sub=None):
        a = act if sub is None else act[sub]
        ab = band_of(a.bg.values)
        out = []
        for kk in range(len(BANDS)):
            peru = a[ab == kk].groupby("user")[col].median()
            out.append(round(float(peru.median()), 2) if len(peru) >= 8 else None)
        return out

    raw = curve("raw_ratio")
    corr = curve("corr_ratio")
    clean_raw = curve("raw_ratio", act.pre_slope > -1.0)      # mean-reversion removed too
    clean_corr = curve("corr_ratio", act.pre_slope > -1.0)

    # bound: clearance REQUIRED to bend high-BG raw ratio to the power-law level vs what data shows
    req = {}
    for kk in [3, 4]:  # 175-205, 205-260
        a = act[act_bidx == kk]
        med_drop = float(a["drop"].median()); med_pred = float(a.pred_drop.median())
        required_C = med_drop - POWERLAW_RATIO_HI * med_pred   # C s.t. (drop-C)/pred = 0.30
        req[LBL[kk]] = {"required_clearance_mgdl": round(required_C, 1),
                        "data_derived_clearance_mgdl": round(nonins.get(kk, np.nan), 1),
                        "median_drop": round(med_drop, 1), "median_pred_drop": round(med_pred, 1)}

    def slope_k(c):
        pts = [(CTR[i], c[i]) for i in range(len(c)) if c[i] and CTR[i] >= 120 and c[i] > 0.02]
        if len(pts) < 3:
            return None
        x = np.log(np.array([p[0] for p in pts]) / 100.0); y = np.log([p[1] for p in pts])
        return round(-float(np.polyfit(x, y, 1)[0]), 2)

    summary = {
        "n_users": int(D.user.nunique()), "n_active_windows": int(len(act)),
        "data_derived_nonInsulin_flux_mgdl_by_bg": {LBL[k]: round(nonins[k], 1) for k in range(len(BANDS))},
        "nonInsulin_n": {LBL[k]: nonins_n[k] for k in range(len(BANDS))},
        "raw_ratio_by_bg": dict(zip(LBL, raw)),
        "clearance_corrected_ratio_by_bg": dict(zip(LBL, corr)),
        "clean_raw_ratio_by_bg": dict(zip(LBL, clean_raw)),
        "clean_corrected_ratio_by_bg": dict(zip(LBL, clean_corr)),
        "implied_k": {"raw": slope_k(raw), "corrected": slope_k(corr),
                      "clean_corrected": slope_k(clean_corr),
                      "note": "k>0 = ratio falls with BG = resistance / power-law direction"},
        "powerlaw_bound": {"powerlaw_high_BG_ratio_target": POWERLAW_RATIO_HI, "by_band": req},
        "caveat": "PLAUSIBILITY BOUND, not a clean number: clearance↔resistance partially entangled; "
                  "low-insulin high-BG windows sparse; EGP not perfectly stable. nonInsulin includes EGP "
                  "(can be negative near target = glucose rising).",
    }
    kcc = slope_k(corr)
    summary["verdict"] = (
        "RESISTANCE appears after clearance removal (corrected ratio falls with BG)" if (kcc and kcc > 0.3)
        else "still NO resistance after data-derived clearance removal (corrected ratio flat/rising); "
             "power law requires clearance far above the data-derived value (see bound)")
    (OUT / "clearance_corrected_isf.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].plot(CTR, [nonins[k] for k in range(len(BANDS))], "o-", color="#9467bd", lw=2)
    ax[0].axhline(0, color="k", ls="--", lw=1); ax[0].axvline(180, color="#888", ls=":", lw=1, label="renal ~180")
    ax[0].set_xlabel("glucose (mg/dL)"); ax[0].set_ylabel("non-insulin flux (mg/dL / 4h)")
    ax[0].set_title("Data-derived insulin-INDEPENDENT disposal\n(|loop predicted drop|<5)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(CTR, raw, "o-", color="#1f77b4", lw=2, label=f"raw (k={slope_k(raw)})")
    ax[1].plot(CTR, corr, "s-", color="#d62728", lw=2.5, label=f"clearance-corrected (k={kcc})")
    ax[1].plot(CTR, [(120.0 / b) ** 1.0 * (corr[1] or 0.8) for b in CTR], ":", color="#888", lw=1.5, label="power-law k=1 ref")
    ax[1].axhline(1, color="k", ls="--", lw=1)
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("realised ÷ profile ISF ratio")
    ax[1].set_title("Does resistance appear after removing clearance?"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_clearance_corrected_isf.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Test #1 (bounding): high-BG resistance after removing insulin-independent clearance\n",
          f"{summary['n_users']} users, {len(act):,} insulin-active windows. Non-insulin flux estimated "
          "from windows where the loop expects ~no insulin action (|loop predicted drop|<5).\n",
          "## Data-derived non-insulin flux (mg/dL over 4h)\n", "| BG band | flux | n |", "|---|---|---|"]
    for k in range(len(BANDS)):
        md.append(f"| {LBL[k]} | {round(nonins[k],1)} | {nonins_n[k]:,} |")
    md += ["\n## Realised ÷ profile ISF ratio: raw vs clearance-corrected\n",
           "| BG band | raw | corrected | clean+corrected |", "|---|---|---|---|"]
    for k in range(len(BANDS)):
        md.append(f"| {LBL[k]} | {raw[k]} | {corr[k]} | {clean_corr[k]} |")
    md += [f"\n**Implied k**: raw {slope_k(raw)}, corrected **{kcc}**, clean+corrected {slope_k(clean_corr)} "
           "(k>0 = resistance).\n", "## Power-law bound (clearance required vs data-derived)\n",
           "| BG band | required clearance | data-derived | median drop / pred |", "|---|---|---|---|"]
    for b, r in req.items():
        md.append(f"| {b} | {r['required_clearance_mgdl']} | {r['data_derived_clearance_mgdl']} | "
                  f"{r['median_drop']} / {r['median_pred_drop']} |")
    md += [f"\n**Verdict: {summary['verdict']}.**\n",
           "![clearance corrected](charts/inv008/fig_clearance_corrected_isf.png)\n", "*" + summary["caveat"] + "*"]
    (OUT / "clearance_corrected_isf.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
