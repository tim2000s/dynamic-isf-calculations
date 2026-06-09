#!/usr/bin/env python3
"""(c) The IOB-model magnitude bias: is it an ISF-expressible, deployable correction or not?

The loop over-predicts the drop, worse the larger the predicted drop. The decisive question for an
oref-deployable ISF algorithm is the SHAPE of actual_drop vs predicted_drop:

  actual_drop = a + b·pred_drop                 (affine)
    intercept a  → baseline glucose drift (basal/endogenous) — NOT an ISF effect
    slope    b<1 → constant insulin OVER-scale  — a single effective-ISF / activity multiplier
    A positive intercept ALONE makes the ratio actual/pred shrink with magnitude WITHOUT any
    nonlinearity — so the "shrink with pred_drop" can be fully affine.
  + curvature (pred_drop²)                       → genuine SATURATION → needs an IOB-conditioned ISF

Decision: if the relationship is affine (no curvature), the deployable fix is a constant
activity/ISF rescale (×b) plus a baseline term — no IOB-aware ISF needed, and it can live in the
insulin-action model rather than overloading ISF. If curvature is real and IOB-driven, an
IOB-conditioned ISF is warranted. Either way it MUST transfer across platforms (Trio v5 ↔ oref0 v7)
to be addable to oref instantiations — that is tested explicitly.

Parquet-only (no DB, no Pool). Output: results/magnitude_bias.{json,md},
charts/inv008/fig_magnitude_bias.png.  Run: python -m inv008.magnitude_bias
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config, err_common as ec
from inv008.bridge_diabeloop import louo_multiplier

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"


def per_user_fits(d):
    """Per-user affine (a, b) and quadratic curvature c of actual_drop vs pred_drop (scaled /100)."""
    rec = {}
    for u, g in d.groupby("user"):
        if len(g) < 30:
            continue
        p = g.pred_drop.values / 100.0          # scale x for conditioning; convert coeffs back below
        y = g.actual_drop.values
        lin = ec.ols(y, np.column_stack([np.ones(len(g)), p]))
        quad = ec.ols(y, np.column_stack([np.ones(len(g)), p, p**2]))
        if lin is None or quad is None:
            continue
        # convert to real units: slope b is per (pred/100) → /100; curvature c per (pred/100)² → /1e4
        rec[u] = {"a": lin[0][0], "a_se": lin[1][0], "b": lin[0][1] / 100, "b_se": lin[1][1] / 100,
                  "c": quad[0][2] / 1e4, "c_se": quad[1][2] / 1e4, "platform": g.table.iloc[0]}
    return pd.DataFrame(rec).T


def pool(F, coef):
    r = ec.dersimonian_laird(F[coef].astype(float).values, F[coef + "_se"].astype(float).values)
    return {k: round(v, 4) for k, v in r.items()} if r else None


def main():
    d = ec.load_windows().copy()
    d["actual_drop"] = d.bg - d.bg_end
    d["pred_drop"] = d.err_static + d.actual_drop
    d = d[d.pred_drop > 0].reset_index(drop=True)
    d["activity"] = d.pred_drop / d.profile_isf          # activity-integral proxy (IOB-action)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)

    F = per_user_fits(d)
    pooled = {"intercept_a_mgdl": pool(F, "a"), "slope_b": pool(F, "b"),
              "curvature_c_per100": pool(F, "c")}

    # shrink curve: median actual/pred ratio by pred_drop decile and by IOB decile
    d["ratio"] = (d.actual_drop / d.pred_drop).clip(-1, 2)
    def decile_curve(col):
        q = pd.qcut(d[col], 8, duplicates="drop")
        g = d.groupby(q, observed=True)
        return [{"x": round(float(d[col][g.groups[k]].median()), 1),
                 "ratio": round(float(d.ratio[g.groups[k]].median()), 2),
                 "n": int(len(g.groups[k]))} for k in g.groups]
    shrink = {"by_pred_drop": decile_curve("pred_drop"), "by_iob": decile_curve("iob")}

    # which conditioning variable captures the bias, out-of-user (LOUO MAE of m·pred_drop)
    def louo_mae(cols):
        e = louo_multiplier(d, cols)
        s = pd.DataFrame({"u": d.user.values, "e": np.abs(e)}).groupby("u").e.median()
        return round(float(s.median()), 2)
    static_mae = round(float(pd.DataFrame({"u": d.user.values, "e": (d.actual_drop - d.pred_drop).abs()})
                              .groupby("u").e.median().median()), 2)
    conditioning = {"static": static_mae,
                    "m(pred_drop)": louo_mae(["pred_drop"]),
                    "m(iob)": louo_mae(["iob"]),
                    "m(activity)": louo_mae(["activity"]),
                    "m(iob,pred_drop)": louo_mae(["iob", "pred_drop"])}

    # cross-platform transfer: fit pooled affine per platform, apply each to the other
    plat = {}
    for p in ["oref_v5", "oref_v7"]:
        Fp = F[F.platform == p]
        plat[p] = {"a": float(Fp.a.astype(float).median()), "b": float(Fp.b.astype(float).median()),
                   "n_users": int(len(Fp))}
    def apply_model(src, tgt):
        a, b = plat[src]["a"], plat[src]["b"]
        m = d.table == tgt
        e = (a + b * d.pred_drop[m]) - d.actual_drop[m]
        s = pd.DataFrame({"u": d.user[m].values, "e": e.abs()}).groupby("u").e.median()
        return round(float(s.median()), 2)
    cross = {"v7_model_on_v5": apply_model("oref_v7", "oref_v5"),
             "v5_model_on_v7": apply_model("oref_v5", "oref_v7"),
             "v7_on_v7": apply_model("oref_v7", "oref_v7"),
             "v5_on_v5": apply_model("oref_v5", "oref_v5")}

    summary = {
        "n_users": int(d.user.nunique()), "n_windows": int(len(d)),
        "pooled_affine": pooled,
        "interpretation": {
            "slope_b": "actual_drop ≈ b·pred_drop; b<1 = loop over-scales insulin action by 1/b "
                       "(a single constant effective-ISF/activity multiplier could fix it).",
            "intercept_a": "baseline drift independent of the prediction (basal/endogenous, NOT ISF). "
                           "A positive a alone makes actual/pred shrink with magnitude affinely.",
            "curvature_c": "coef of pred_drop²; ~0 (CI spans 0) ⇒ no saturation ⇒ the bias is AFFINE "
                           "⇒ a constant rescale + baseline, NO IOB-conditioned ISF needed. "
                           "Clearly negative ⇒ real saturation ⇒ IOB-dependent correction warranted.",
        },
        "shrink_curve": shrink,
        "conditioning_louo_mae": conditioning,
        "iob_adds_over_pred_drop": round(conditioning["m(pred_drop)"] - conditioning["m(iob,pred_drop)"], 2),
        "cross_platform": cross, "per_platform_affine": plat,
    }
    # curvature consistency across users (sign agreement); pooled-significant ≠ universal
    cc = pooled["curvature_c_per100"]
    curv_consistent = bool(cc and cc["frac_same_sign"] >= 0.8 and cc["I2_pct"] < 90)
    summary["curvature_consistent_across_users"] = curv_consistent
    summary["verdict"] = ("AFFINE (no consistent saturation) — constant rescale captures it"
                          if not curv_consistent else "CURVED (consistent saturation) — IOB-dependence warranted")
    (OUT / "magnitude_bias.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    # 1: actual vs predicted drop with the pooled affine line (binned medians)
    qb = pd.qcut(d.pred_drop, 12, duplicates="drop")
    bx = d.groupby(qb, observed=True).pred_drop.median()
    byy = d.groupby(qb, observed=True).actual_drop.median()
    ax[0].plot(bx, byy, "o", color="#1f77b4", label="binned median")
    xs = np.linspace(d.pred_drop.min(), d.pred_drop.quantile(0.99), 50)
    a, b = pooled["intercept_a_mgdl"]["b_re"], pooled["slope_b"]["b_re"]
    ax[0].plot(xs, a + b * xs, "-", color="#d62728", lw=2, label=f"affine a={a:.0f}, b={b:.2f}")
    ax[0].plot(xs, xs, "--", color="#999", label="y=x (perfect)")
    ax[0].set_xlabel("predicted drop (mg/dL)"); ax[0].set_ylabel("actual drop (mg/dL)")
    ax[0].set_title("Actual vs predicted drop\n(slope b<1 = constant over-scale)"); ax[0].legend(fontsize=8)
    # 2: shrink ratio by pred_drop and by IOB
    sp = shrink["by_pred_drop"]; si = shrink["by_iob"]
    ax[1].plot([p["x"] for p in sp], [p["ratio"] for p in sp], "o-", color="#1f77b4", label="vs pred_drop")
    ax2 = ax[1].twiny()
    ax2.plot([p["x"] for p in si], [p["ratio"] for p in si], "s-", color="#2ca02c", label="vs IOB")
    ax[1].axhline(1, color="k", ls="--", lw=1); ax[1].set_ylabel("actual/predicted ratio")
    ax[1].set_xlabel("predicted drop (mg/dL)", color="#1f77b4"); ax2.set_xlabel("IOB (U)", color="#2ca02c")
    ax[1].set_title("Realised/predicted shrink"); ax[1].legend(fontsize=8, loc="upper right")
    # 3: cross-platform affine + conditioning MAE
    names = list(conditioning); ax[2].bar(range(len(names)), [conditioning[n] for n in names],
                                          color=["#999", "#1f77b4", "#2ca02c", "#9467bd", "#d62728"])
    ax[2].set_xticks(range(len(names))); ax[2].set_xticklabels(names, rotation=40, ha="right", fontsize=8)
    ax[2].set_ylabel("out-of-user MAE (mg/dL)")
    ax[2].set_title(f"Conditioning variable\n(cross-plat: v7→v5 {cross['v7_model_on_v5']}, "
                    f"v5→v7 {cross['v5_model_on_v7']})")
    fig.tight_layout(); fig.savefig(CHART / "fig_magnitude_bias.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    cc = pooled["curvature_c_per100"]
    affine = not curv_consistent
    md = ["# (c) The magnitude bias: affine rescale or nonlinear saturation? Deployable in oref?\n",
          f"{summary['n_users']} users, {len(d):,} windows. Fit actual_drop = a + b·pred_drop (+ "
          "curvature). Out-of-user MAEs are per-user median |err|, median across users.\n",
          "## Pooled shape (random-effects)\n", "| term | estimate | p | between-user τ |",
          "|---|---|---|---|",
          f"| intercept a (mg/dL drift) | {pooled['intercept_a_mgdl']['b_re']} | "
          f"{pooled['intercept_a_mgdl']['p']:.1e} | {pooled['intercept_a_mgdl']['tau']} |",
          f"| slope b (insulin scale) | {pooled['slope_b']['b_re']} | {pooled['slope_b']['p']:.1e} | "
          f"{pooled['slope_b']['tau']} |",
          f"| curvature c (per mg/dL²) | {cc['b_re']:.2e} | {cc['p']:.1e} | {cc['tau']:.2e} |",
          f"\n**Shape verdict: {'AFFINE — saturation not consistent across users' if affine else 'CURVED — consistent saturation'}** "
          f"(curvature same-sign across users {cc['frac_same_sign']:.0%}, I²={cc['I2_pct']:.0f}%). "
          f"slope b={pooled['slope_b']['b_re']:.2f} ⇒ loop over-scales insulin action by ~{1/pooled['slope_b']['b_re']:.1f}×; "
          f"intercept a={pooled['intercept_a_mgdl']['b_re']:.0f} mg/dL is baseline drift (NOT ISF).\n",
          "## What captures it, out-of-user (MAE mg/dL)\n", "| conditioning | MAE |", "|---|---|"]
    for n, v in conditioning.items():
        md.append(f"| {n} | {v} |")
    md += [f"\nIOB adds over pred_drop: {summary['iob_adds_over_pred_drop']} mg/dL.\n",
           "## Cross-platform transfer (Trio v5 ↔ oref0 v7)\n",
           f"per-platform affine: v5 a={plat['oref_v5']['a']:.0f} b={plat['oref_v5']['b']:.2f} "
           f"(n={plat['oref_v5']['n_users']}); v7 a={plat['oref_v7']['a']:.0f} b={plat['oref_v7']['b']:.2f} "
           f"(n={plat['oref_v7']['n_users']}).\n",
           f"v7-model on v5: MAE {cross['v7_model_on_v5']} (vs v5-on-v5 {cross['v5_on_v5']}); "
           f"v5-model on v7: {cross['v5_model_on_v7']} (vs v7-on-v7 {cross['v7_on_v7']}).\n",
           "![magnitude bias](charts/inv008/fig_magnitude_bias.png)\n",
           "## Deployability read\n",
           ("- The bias is **affine** → a **constant effective-insulin rescale (×b) plus a baseline "
            "term** captures it; **no IOB-conditioned ISF is required**. Cleanest home is the "
            "insulin-action model (scale activity by b — affects prediction only), not ISF (which is "
            "overloaded: it also sets correction dose). " if affine else
            "- The bias is **nonlinear (saturation)** → an **IOB-conditioned correction** is "
            "warranted; a constant ISF cannot capture it. ") +
           ("If the per-platform b/a and the cross-applied MAEs match, the correction is "
            "**platform-invariant → addable to all oref instantiations**; if they diverge it is "
            "algorithm-specific.")]
    (OUT / "magnitude_bias.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
