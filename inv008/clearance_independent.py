#!/usr/bin/env python3
"""Model-INDEPENDENT clearance decomposition: does the resistance↔clearance reconciliation survive?

Re-runs the clearance decomposition (test #1) without the loop's insulin-action curve. Insulin that
acted is computed by conservation (inv008.effective_isf_independent): acted = ΔIOB + SMBs +
∫(temp−profile basal). Non-insulin flux is estimated from windows where ~no net insulin acted
(|acted| < 0.3 U) — there the 4h glucose change is renal/mass-action clearance − EGP. We then split
the insulin-active windows (acted ≥ 0.5 U) into net and clearance-corrected (insulin-only) effective
ISF, by glucose.

    net ratio        = (drop / acted) / profile_isf
    corrected ratio  = ((drop − nonInsulin(BG)) / acted) / profile_isf

If — as in the loop-model version — the NET ratio is flat/rising but the CORRECTED (insulin-only)
ratio FALLS with glucose, the reconciliation holds model-independently: resistance is real but offset
by clearance in the net. Single-process. Output: results/clearance_independent.{json,md},
charts/inv008/fig_clearance_independent.png. Run: python -m inv008.clearance_independent
"""
from __future__ import annotations

import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
from inv008 import config
from inv008.effective_isf_independent import user_windows, TBL, BANDS, LBL, CTR

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
ZERO_ACTED = 0.3        # |acted| below this ⇒ ~no insulin → drop is non-insulin flux
ACTIVE = 0.5            # acted at/above this ⇒ insulin-active window


def band_of(bg):
    out = np.full(len(bg), -1)
    for k, (a, b) in enumerate(BANDS):
        out[(bg >= a) & (bg < b)] = k
    return out


def slope_k(c):
    pts = [(CTR[i], c[i]) for i in range(len(c)) if c[i] and CTR[i] >= 120 and c[i] > 0.02]
    if len(pts) < 3:
        return None
    return round(-float(np.polyfit(np.log(np.array([p[0] for p in pts]) / 100.0),
                                   np.log([p[1] for p in pts]), 1)[0]), 2)


def main():
    coh = {r["user_id"]: r for r in json.load(open(config.ROOT / "canonical_cohort.json"))}
    basal = json.load(open(config.ROOT / "user_basal_profiles.json"))
    jobs = [(u, TBL[r["cohort"]], r.get("isf")) for u, r in coh.items()
            if r.get("cohort") in TBL and r.get("isf") and u in basal]
    print(f"independent clearance: {len(jobs)} users, single-process")
    parts = []
    for k, (u, table, isf) in enumerate(jobs):
        df = user_windows(u, table, isf, basal[u]["hourly_rates"], min_acted=-1e9)  # keep ALL windows
        if df is not None:
            parts.append(df)
        if k % 25 == 0:
            print(f"  {k}/{len(jobs)}")
    D = pd.concat(parts, ignore_index=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    bidx = band_of(D.bg.values)

    # non-insulin flux from near-zero-insulin windows
    low = D.acted.abs().values < ZERO_ACTED
    nonins, nonins_n = {}, {}
    for kk in range(len(BANDS)):
        m = low & (bidx == kk)
        nonins[kk] = float(np.median(D["drop"].values[m])) if m.sum() >= 40 else np.nan
        nonins_n[kk] = int(m.sum())

    act = D[D.acted >= ACTIVE].copy()
    ab = band_of(act.bg.values)
    act["net_ratio"] = (act["drop"] / act.acted / act.profile_isf).clip(-1, 4)
    act["corr_ratio"] = ((act["drop"] - np.array([nonins.get(b, np.nan) for b in ab])) / act.acted
                         / act.profile_isf).clip(-1, 4)

    def curve(col):
        a = act; ab = band_of(a.bg.values)
        return [round(float(a[ab == kk].groupby("user")[col].median().median()), 2)
                if (ab == kk).sum() >= 100 and a[ab == kk].user.nunique() >= 6 else None
                for kk in range(len(BANDS))]

    net, corr = curve("net_ratio"), curve("corr_ratio")
    k_net, k_corr = slope_k(net), slope_k(corr)
    nonins_curve = [round(nonins[k], 1) if np.isfinite(nonins[k]) else None for k in range(len(BANDS))]

    summary = {
        "n_users": int(D.user.nunique()), "n_active": int(len(act)),
        "method": "model-independent: acted = ΔIOB + SMBs + ∫(temp−profile basal); non-insulin flux "
                  "from |acted|<0.3 U windows.",
        "nonInsulin_flux_mgdl": dict(zip(LBL, nonins_curve)),
        "nonInsulin_n": dict(zip(LBL, [nonins_n[k] for k in range(len(BANDS))])),
        "net_ratio_by_bg": dict(zip(LBL, net)),
        "clearance_corrected_ratio_by_bg": dict(zip(LBL, corr)),
        "implied_k": {"net": k_net, "corrected": k_corr,
                      "note": "k>0 = ratio falls with BG = resistance"},
        "verdict": None,
    }
    net_flat = (k_net is None) or (k_net <= 0.15)
    resist = (k_corr is not None) and (k_corr > 0.2)
    summary["verdict"] = (
        "RECONCILIATION HOLDS model-independently — net effective ISF flat/rising "
        f"(k={k_net}) while clearance-corrected insulin-only ISF falls with glucose (k={k_corr}); "
        "resistance is real but offset by clearance in the net"
        if (net_flat and resist) else
        f"net k={k_net}, corrected k={k_corr} — does not cleanly reproduce the loop-model pattern; inspect")
    (OUT / "clearance_independent.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].plot(CTR, nonins_curve, "o-", color="#9467bd", lw=2)
    ax[0].axhline(0, color="k", ls="--", lw=1); ax[0].axvline(180, color="#888", ls=":", lw=1, label="renal ~180")
    ax[0].set_xlabel("glucose (mg/dL)"); ax[0].set_ylabel("non-insulin flux (mg/dL/4h)")
    ax[0].set_title("Model-independent clearance\n(|acted|<0.3 U windows)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(CTR, net, "o-", color="#1f77b4", lw=2, label=f"net (k={k_net})")
    ax[1].plot(CTR, corr, "s-", color="#d62728", lw=2.5, label=f"clearance-corrected insulin-only (k={k_corr})")
    ax[1].axhline(1, color="k", ls="--", lw=1)
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("effective ÷ profile ISF")
    ax[1].set_title("Resistance↔clearance, model-independent"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_clearance_independent.png", dpi=150); plt.close(fig)

    md = ["# Model-independent clearance decomposition (no insulin-action curve)\n",
          f"{summary['n_users']} users, {len(act):,} insulin-active windows. acted = ΔIOB + SMBs + "
          "∫(temp−profile basal); non-insulin flux from |acted|<0.3 U windows.\n",
          "## Non-insulin flux (mg/dL/4h)\n", "| BG band | flux | n |", "|---|---|---|"]
    for k in range(len(BANDS)):
        md.append(f"| {LBL[k]} | {nonins_curve[k]} | {nonins_n[k]:,} |")
    md += ["\n## Effective ISF ratio: net vs clearance-corrected (insulin-only)\n",
           "| BG band | net | clearance-corrected |", "|---|---|---|"]
    for k in range(len(BANDS)):
        md.append(f"| {LBL[k]} | {net[k]} | {corr[k]} |")
    md += [f"\nImplied k: net **{k_net}**, corrected **{k_corr}** (k>0 = falls with glucose = resistance).\n",
           f"**Verdict: {summary['verdict']}.**\n",
           "![clearance independent](charts/inv008/fig_clearance_independent.png)\n"]
    (OUT / "clearance_independent.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
