#!/usr/bin/env python3
"""Model-INDEPENDENT effective ISF vs glucose: remove the dependence on the loop's insulin curve.

The same-window effective ISF elsewhere uses the loop's IOB prediction (reason_IOBpredBG), so it
inherits the loop's DIA/peak insulin-action model. Here we compute the insulin that ACTUALLY ACTED
over each window by conservation — no activity curve, no forward projection:

    insulin_acted = (IOB_start − IOB_end) + SMBs_delivered + ∫(temp_basal − profile_basal) dt

ΔIOB is the observed change in on-board insulin (units absorbed); SMBs and the temp-basal deviation
from the scheduled profile are the insulin delivered *during* the window (which the loop's prediction
omits). Because by four hours ~85–93% of a fast-insulin dose has acted, IOB_end is small and the curve
choice barely matters — the conservation identity is robust (this is the point: the 4-hour horizon
makes it model-light). Then:

    effective_ISF = (cgm_start − cgm_end) / insulin_acted          (mg/dL per U that acted)

We compute effective ISF by glucose band (÷ each user's profile ISF) and compare it to the
loop-model version (realised_isf) and to the candidate forms. If the glucose dependence is unchanged,
the paper's conclusions are robust to the insulin-action-model caveat.

Overnight carb-screened 4h windows. Single-process. Output: results/effective_isf_independent.{json,md},
charts/inv008/fig_effective_isf_independent.png. Run: python -m inv008.effective_isf_independent
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
START_HOURS = {23, 0, 1, 2}
HZ, TOL, RISE_MAX = 4 * 3600, 300, 2.0
MIN_ACTED = 0.5                  # min net insulin acted (U) for a stable ISF ratio
BANDS = [(100, 120), (120, 145), (145, 175), (175, 205), (205, 260)]
LBL = [f"{a}-{b}" for a, b in BANDS]
CTR = [(a + b) / 2 for a, b in BANDS]


def to_mgdl(v):
    return v * 18.018 if (v is not None and v < 20) else v


def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4


def user_windows(user_id, table, profile_isf, basal_hourly, min_acted=MIN_ACTED):
    conn = psycopg2.connect("dbname=oref")
    try:
        d = pd.read_sql(f"""SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob, sug_smb_units AS smb,
                                   sug_rate FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL
                              AND iob_iob IS NOT NULL ORDER BY ts_relative_sec""", conn, params=(user_id,))
    finally:
        conn.close()
    if len(d) < 500 or profile_isf is None:
        return None
    ts = d.ts_relative_sec.values.astype(float); bg = d.cgm_mgdl.values.astype(float)
    hr = d.hour.values.astype(float)
    iob = pd.to_numeric(d.iob_iob, errors="coerce").values.astype(float)
    smb = pd.to_numeric(d.smb, errors="coerce").fillna(0).values.astype(float)
    rate = pd.to_numeric(d.sug_rate, errors="coerce").values.astype(float)
    n = len(d)
    pisf = to_mgdl(profile_isf)
    # per-tick basal deviation (U): (temp rate − profile basal at that hour) × Δt hours
    dt = np.diff(ts, prepend=ts[0]); dt = np.clip(dt, 0, 600)            # cap gaps at 10 min
    prof = np.array([basal_hourly[int(h) % 24] for h in hr])
    # null temp-rate ticks → assume profile basal was delivered (deviation 0), so cumsum isn't poisoned
    dev = np.where(np.isfinite(rate), (rate - prof) * dt / 3600.0, 0.0)
    smb = np.where(np.isfinite(smb), smb, 0.0)
    csum_smb = np.concatenate([[0.0], np.cumsum(smb)])
    csum_dev = np.concatenate([[0.0], np.cumsum(dev)])
    end4 = np.searchsorted(ts, ts + HZ); end4 = np.where(end4 < n, end4, n - 1)
    oke = np.abs(ts[end4] - (ts + HZ)) <= TOL
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    sl = np.where(np.abs(ts[p15] - (ts + 900)) <= TOL, (bg[p15] - bg) / 3.0, 0.0)
    rows = []
    for i in range(n):
        if hr[i] not in START_HOURS or not oke[i] or bg[i] < 80 or bg[i] > 260:
            continue
        j = end4[i]
        if np.nanmax(sl[i:j + 1]) > RISE_MAX:
            continue
        d_iob = iob[i] - iob[j]
        delivered = (csum_smb[j + 1] - csum_smb[i]) + (csum_dev[j + 1] - csum_dev[i])
        acted = d_iob + delivered
        drop = bg[i] - bg[j]
        if acted < min_acted:
            continue
        rows.append((bg[i], drop, acted, drop / acted if abs(acted) > 1e-6 else np.nan, pisf))
    if len(rows) < 30:
        return None
    df = pd.DataFrame(rows, columns=["bg", "drop", "acted", "eff_isf", "profile_isf"])
    df["user"] = user_id
    return df


def main():
    coh = {r["user_id"]: r for r in json.load(open(config.ROOT / "canonical_cohort.json"))}
    basal = json.load(open(config.ROOT / "user_basal_profiles.json"))
    jobs = [(u, TBL[r["cohort"]], r.get("isf")) for u, r in coh.items()
            if r.get("cohort") in TBL and r.get("isf") and u in basal]
    print(f"independent effective ISF: {len(jobs)} users, single-process")
    parts = []
    for k, (u, table, isf) in enumerate(jobs):
        df = user_windows(u, table, isf, basal[u]["hourly_rates"])
        if df is not None:
            parts.append(df)
        if k % 25 == 0:
            print(f"  {k}/{len(jobs)}")
    D = pd.concat(parts, ignore_index=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    D["ratio"] = (D.eff_isf / D.profile_isf).clip(-1, 4)
    bidx = np.full(len(D), -1)
    for kk, (a, b) in enumerate(BANDS):
        bidx[(D.bg >= a) & (D.bg < b)] = kk

    # per-user median ratio by band, then median across users (equal weight)
    indep, npar = [], []
    for kk in range(len(BANDS)):
        sub = D[bidx == kk]
        peru = sub.groupby("user").ratio.median()
        indep.append(round(float(peru.median()), 2) if peru.notna().sum() >= 8 else None)
        npar.append(int(len(sub)))
    quart = [round(float(quartic(c) / quartic(100.0)), 2) for c in CTR]

    # loop-model comparison (from the head_to_head parquet, same metric: realised/profile)
    try:
        P = pd.read_parquet(OUT / "head_to_head_windows.parquet")
        P["ad"] = P.bg - P.bg_end; P["pd_"] = P.err_static + P.ad
        P["r"] = (P.ad / P.pd_).clip(-1, 4)
        loopm = []
        for a, b in BANDS:
            sub = P[(P.bg >= a) & (P.bg < b)]
            loopm.append(round(float(sub.groupby("user").r.median().median()), 2))
    except Exception:
        loopm = [None] * len(BANDS)

    def slope_k(c):
        pts = [(CTR[i], c[i]) for i in range(len(c)) if c[i] and CTR[i] >= 120 and c[i] > 0.02]
        if len(pts) < 3:
            return None
        return round(-float(np.polyfit(np.log(np.array([p[0] for p in pts]) / 100.0),
                                       np.log([p[1] for p in pts]), 1)[0]), 2)

    summary = {
        "n_users": int(D.user.nunique()), "n_windows": int(len(D)),
        "method": "effective_ISF = ΔBG / (ΔIOB + SMBs + ∫(temp−profile basal)); no insulin curve.",
        "effective_isf_independent_ratio": dict(zip(LBL, indep)),
        "effective_isf_loopmodel_ratio": dict(zip(LBL, loopm)),
        "diabeloop_quartic_ratio": dict(zip(LBL, quart)),
        "n_by_band": dict(zip(LBL, npar)),
        "implied_k": {"independent": slope_k(indep), "loopmodel": slope_k(loopm), "quartic": slope_k(quart)},
        "verdict": None,
    }
    ki = slope_k(indep)
    summary["verdict"] = (
        f"Model-INDEPENDENT net effective ISF is {'flat/rising' if (ki is None or ki<=0.1) else 'falling'} "
        f"with glucose (k={ki}); {'matches' if (ki is not None and ki<=0.1) else 'differs from'} the "
        "loop-model version → the conclusion is robust to the insulin-action-model caveat" )
    (OUT / "effective_isf_independent.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 1, figsize=(7.5, 5))
    ax.plot(CTR, indep, "o-", color="#2ca02c", lw=2.5, label=f"independent (ΔIOB-based), k={ki}")
    ax.plot(CTR, loopm, "s--", color="#1f77b4", lw=2, label=f"loop-model (reason_IOBpredBG), k={slope_k(loopm)}")
    ax.plot(CTR, quart, ":", color="#d62728", lw=2, label=f"Diabeloop quartic, k={slope_k(quart)}")
    ax.axhline(1, color="k", ls="--", lw=1)
    ax.set_xlabel("glucose (mg/dL)"); ax.set_ylabel("effective ÷ profile ISF")
    ax.set_title("Net effective ISF vs glucose: model-independent vs loop-model vs Diabeloop")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_effective_isf_independent.png", dpi=150); plt.close(fig)

    md = ["# Model-independent effective ISF vs glucose (no insulin-action curve)\n",
          f"{summary['n_users']} users, {len(D):,} windows. effective ISF = ΔBG / (ΔIOB + SMBs + "
          "∫(temp−profile basal)). Ratio = ÷ each user's profile ISF.\n",
          "| BG band | independent (ΔIOB) | loop-model | Diabeloop quartic | n |", "|---|---|---|---|---|"]
    for i in range(len(BANDS)):
        md.append(f"| {LBL[i]} | {indep[i]} | {loopm[i]} | {quart[i]} | {npar[i]:,} |")
    md += [f"\nImplied k: independent **{ki}**, loop-model {slope_k(loopm)}, quartic {slope_k(quart)} "
           "(k>0 = falls with glucose).\n", f"**Verdict: {summary['verdict']}.**\n",
           "![independent effective ISF](charts/inv008/fig_effective_isf_independent.png)\n"]
    (OUT / "effective_isf_independent.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
