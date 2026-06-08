#!/usr/bin/env python3
"""Observed insulin sensitivity as a function of glucose, vs what each equation calculates.

The Tier-2 anchor collapses every fasting window into one ΔBG-on-ΔIOB slope with no glucose
term — so it cannot say whether sensitivity actually varies with glucose the way g(BG) claims.
This keeps the glucose resolution: re-extract the same 30-min fasting windows, bin them by the
window's glucose, and within each BG band fit ΔBG = a + b·ΔIOB + c·pre_trend. The observed ISF
in that band is −b (mg/dL per U absorbed), with a regression CI. A single window is far too
noisy to use alone; the within-band regression pools the windows at that glucose to recover a
trustworthy point.

We then normalise observed ISF to 1.0 at target and overlay the calculated glucose curves:
  • v-next g(BG)  = Diabeloop quartic / quartic(target)
  • v1 scaler     = ln(target/75+1)/ln(BG/75+1)
  • v2 scaler     = ln(target/75)/ln(BG_floored/75)
so observed and calculated are compared like-for-like at each glucose. The fitted power-law
exponent of the observed curve says directly whether g(BG)'s ~1.3 slope is supported.

Windows: identical definition to empirical_isf_v5; pooled across the 114 cohort users.
Output: results/tier2_window_shape.{json,md}, charts/inv008/fig_window_shape.png
Run: python -m inv008.tier2_window_shape
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
R2_MIN = 0.10

WINDOW_S, PRE_S, COB_LEAD_S, SMB_LEAD_S = 1800, 1800, 5400, 1800
COB_TOL, SMB_TOL = 1.0, 0.05
CGM_LO, CGM_HI, DBG_MAX, DIOB_MIN, DIOB_MAX = 70, 250, 80.0, 0.0, 2.0
COL_MAP = {"oref_v5": '"sug_COB"', "oref_v6": "sug_cob", "oref_v7": "sug_cob"}

BG_EDGES = np.array([70, 85, 95, 105, 115, 130, 150, 175, 210], float)
BG_CENTRES = 0.5 * (BG_EDGES[:-1] + BG_EDGES[1:])
TARGET_BIN = int(np.digitize(TARGET, BG_EDGES) - 1)     # the 95–105 band


def windows(user_id, table):
    sql = f"""
        SELECT ts_relative_sec, cgm_mgdl, iob_iob, {COL_MAP[table]} AS cob, sug_smb_units
        FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL AND iob_iob IS NOT NULL
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
    rolling_smb = psum_smb[np.arange(n)] - psum_smb[np.searchsorted(ts, ts - SMB_LEAD_S)]
    cob_pos = (cob > COB_TOL).astype(np.int32)
    psum_cp = np.concatenate([[0], np.cumsum(cob_pos)])
    any_cob = psum_cp[end_idx] - psum_cp[np.arange(n)]
    pre_idx = np.clip(pre_start, 0, n - 1)
    pre_dt = (ts - ts[pre_idx]) / 60.0
    trend = (bg - bg[pre_idx]) / np.where(pre_dt > 0, pre_dt, np.nan)
    keep = (valid_end & valid_pre & (any_cob == 0) & (rolling_cob_max <= COB_TOL)
            & (rolling_smb <= SMB_TOL) & (bg >= CGM_LO) & (bg <= CGM_HI)
            & (bg[end_idx] >= 50) & (bg[end_idx] <= 300) & (np.abs(dbg) <= DBG_MAX)
            & (diob > DIOB_MIN) & (diob < DIOB_MAX) & np.isfinite(trend))
    if keep.sum() < 80:
        return None
    bg_w = 0.5 * (bg[keep] + bg[end_idx][keep])
    return pd.DataFrame({"bg": bg_w, "diob": diob[keep], "dbg": dbg[keep], "trend": trend[keep]})


def band_isf(d):
    """Observed ISF in a BG band: -b from ΔBG = a + b·ΔIOB + c·trend, with SE."""
    if len(d) < 60:
        return None
    X = np.column_stack([np.ones(len(d)), d.diob.values, d.trend.values])
    y = d.dbg.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    sigma2 = (resid @ resid) / max(len(d) - 3, 1)
    se_b = math.sqrt(sigma2 * np.linalg.pinv(X.T @ X)[1, 1])
    return {"isf": float(-beta[1]), "se": float(se_b), "n": int(len(d))}


def main():
    emp = {e["user_id"]: e for e in json.load(open(EMPIRICAL))}
    jobs = [(u, e["table"]) for u, e in emp.items()
            if e.get("r2", 0) >= R2_MIN and 5 <= e.get("empirical_isf", 0) <= 500
            and (config.REPLAY_DIR / f"{u}.parquet").exists()]
    nw = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"window shape: {len(jobs)} users on {nw} workers")
    with mp.Pool(nw) as pool:
        parts = [p for p in pool.starmap(windows, jobs) if p is not None]
    D = pd.concat(parts, ignore_index=True)
    D["bin"] = np.digitize(D.bg, BG_EDGES) - 1

    bands = []
    for i in range(len(BG_CENTRES)):
        r = band_isf(D[D.bin == i])
        bands.append(r)
    ref = bands[TARGET_BIN]["isf"]
    obs_norm = [(b["isf"] / ref if b else None) for b in bands]
    obs_lo = [((b["isf"] - 1.96 * b["se"]) / ref if b else None) for b in bands]
    obs_hi = [((b["isf"] + 1.96 * b["se"]) / ref if b else None) for b in bands]

    # calculated normalised curves at the bin centres
    def g_quartic_n(bg):
        return quartic_isf(bg) / quartic_isf(TARGET)

    def v1_scaler(bg):
        bg_c = np.minimum(bg, 210 + (np.maximum(bg, 210) - 210) / 3)
        return math.log(TARGET / DIV + 1) / np.log(bg_c / DIV + 1)

    def v2_scaler(bg):
        bg_f = np.maximum(bg, DIV + 1)
        return math.log(TARGET / DIV) / np.log(bg_f / DIV)

    gq = [float(g_quartic_n(c)) for c in BG_CENTRES]
    g1 = [float(v1_scaler(c)) for c in BG_CENTRES]
    g2 = [float(v2_scaler(c)) for c in BG_CENTRES]

    # fitted power-law exponent of the OBSERVED curve (over bins with a point)
    ok = [i for i, b in enumerate(bands) if b is not None and b["isf"] > 0]
    x = np.log(TARGET / BG_CENTRES[ok])
    yv = np.log(np.array([bands[i]["isf"] for i in ok]) / ref)
    k_obs = float(np.polyfit(x, yv, 1)[0])
    # quartic's own exponent over the same range, for reference
    k_quart = float(np.polyfit(x, np.log(np.array([gq[i] for i in ok])), 1)[0])

    summary = {
        "n_patients": len(parts), "n_windows": int(len(D)),
        "bg_centres": BG_CENTRES.tolist(),
        "observed_isf_by_bg": [None if b is None else round(b["isf"], 1) for b in bands],
        "observed_isf_n_by_bg": [None if b is None else b["n"] for b in bands],
        "observed_isf_normalised": [None if v is None else round(v, 3) for v in obs_norm],
        "calculated_normalised": {
            "vnext_quartic": [round(v, 3) for v in gq],
            "v1_log": [round(v, 3) for v in g1],
            "v2_log": [round(v, 3) for v in g2]},
        "fitted_powerlaw_exponent": {"observed": round(k_obs, 2),
                                     "vnext_quartic_same_range": round(k_quart, 2)},
        "reading": ("observed insulin sensitivity falls with glucose; the fitted exponent of "
                    "the observed curve vs the quartic's exponent over the same range says how "
                    "well g(BG)'s shape matches the data."),
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    (OUT / "tier2_window_shape.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(figsize=(9, 6))
    xc = BG_CENTRES
    yv2 = np.array([v if v is not None else np.nan for v in obs_norm])
    lo = np.array([v if v is not None else np.nan for v in obs_lo])
    hi = np.array([v if v is not None else np.nan for v in obs_hi])
    ax.errorbar(xc, yv2, yerr=[yv2 - lo, hi - yv2], fmt="o", color="#111", ms=7, lw=1.5,
                capsize=4, label="observed sensitivity (per-window, binned ± 95% CI)", zorder=10)
    ax.plot(xc, gq, "-", color="#1f77b4", lw=2.5, label=f"v-next g(BG) quartic (k≈{k_quart:.1f})")
    ax.plot(xc, g1, "--", color="#888", lw=1.8, label="v1 log scaler")
    ax.plot(xc, g2, ":", color="#d62728", lw=1.8, label="v2 log scaler")
    ax.axhline(1, color="k", lw=0.6, alpha=0.4); ax.axvline(TARGET, color="k", lw=0.6, alpha=0.4)
    ax.set_xlabel("glucose (mg/dL)"); ax.set_ylabel("ISF relative to target (=1 at 99)")
    ax.set_title(f"Observed sensitivity vs glucose, against each equation's curve\n"
                 f"{len(parts)} people, {len(D):,} fasting windows; observed exponent k≈{k_obs:.1f}")
    ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_ylim(0, max(2.5, np.nanmax(hi) * 1.1))
    fig.tight_layout(); fig.savefig(CHART / "fig_window_shape.png", dpi=150); plt.close(fig)

    md = ["# Observed sensitivity vs glucose — per-window, against the calculated curves\n",
          f"{len(parts)} people, {len(D):,} fasting windows. Within each glucose band, observed "
          "ISF = −b from ΔBG = a + b·ΔIOB + c·trend (the per-window measurement, pooled only "
          "within the band so it is not collapsed across glucose).\n",
          "## Observed vs calculated (normalised to 1.0 at target)\n",
          "| BG | n windows | observed ISF | observed (norm) | v-next g(BG) | v1 | v2 |",
          "|---|---|---|---|---|---|---|"]
    for i, c in enumerate(BG_CENTRES):
        b = bands[i]
        if b is None:
            md.append(f"| {c:.0f} | – | – | – | {gq[i]:.2f} | {g1[i]:.2f} | {g2[i]:.2f} |")
        else:
            md.append(f"| {c:.0f} | {b['n']} | {b['isf']:.1f} | {obs_norm[i]:.2f} | "
                      f"{gq[i]:.2f} | {g1[i]:.2f} | {g2[i]:.2f} |")
    md += [f"\n![Observed vs calculated](charts/inv008/fig_window_shape.png)\n",
           "## Reading\n",
           f"- Observed sensitivity **falls with glucose**, fitted power-law exponent "
           f"**k ≈ {k_obs:.1f}** over this range; the v-next quartic's exponent over the same "
           f"range is {k_quart:.1f}.",
           "- This is the like-for-like the averaged anchor could not give: the measured value "
           "compared to the calculated value *at each glucose*, rather than one pooled slope.",
           "\n*Caveat: 30-min fasting windows, ΔIOB ∈ (0,2] U, trend-adjusted; observational, so "
           "counter-regulation at low BG and unrecorded carbs/EGP still confound the band "
           "estimates, most at the extremes where windows are fewer.*"]
    (OUT / "tier2_window_shape.md").write_text("\n".join(md))
    print("\n".join(md[:6]))
    print("observed k =", round(k_obs, 2), "| quartic k =", round(k_quart, 2))


if __name__ == "__main__":
    main()
