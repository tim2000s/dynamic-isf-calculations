#!/usr/bin/env python3
"""Best-fit individualised glucose-ISF: let the data choose the glucose-dependence per user.

Earlier bridging IMPOSED fixed steep shapes anchored at profile@100 → they failed. That tested the
wrong thing. Here we FIT, the way the Diabeloop curve is meant to be used: a shared glucose SHAPE
g(BG;k) with a per-user SCALE (and baseline), chosen by optimisation. k=0 is flat (static) as a
special case, so the data decides how much glucose-dependence each person wants.

Model (linear in the per-user nuisance params given k, so the scale/baseline profile out and we
search the shared shape steepness k — a separable nonlinear least squares / "gradient" fit):

    actual_drop ≈ a_u + s_u · ( pred_drop · g(BG; k) )        g(BG;k) = (100/BG)^k,  g(·;0)=1

  a_u = per-user baseline (drift/basal), s_u = per-user effective-sensitivity scale (individualises
  the curve, the Diabeloop multiplicative factor), k = SHARED glucose steepness (the scientific
  question). k=0 reduces to the affine magnitude model; k>0 adds glucose shape.

Honesty: (a_u, s_u) are fit by WITHIN-USER 5-fold CV at every k, so a per-user level fit cannot
overfit; the shared k is the only thing read across users (1 dof / cohort). We compare:
    static        pred_drop − actual           (no fit)
    magnitude     k=0 (a_u,s_u fit)             (the strong baseline to beat)
    best-fit k*   argmin over k                 (does a glucose shape help BEYOND magnitude?)
    diabeloop     a_u + s_u·pred_drop·H(BG)/H(100)   (the literal Diabeloop shape, per-user scaled)
    cold-start    s_u=1,a_u=0 swept over k      (population shape, NO adaptation — deployability floor)
Also: per-user optimal k distribution, and whether s_u is predictable from TDD/profile (cold-start).
Parquet-only. Output: results/gradient_isf_fit.{json,md}, charts/inv008/fig_gradient_isf_fit.png
Run: python -m inv008.gradient_isf_fit
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import KFold

from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
KGRID = np.round(np.arange(0.0, 4.01, 0.25), 2)
MIN_W = 50


def hybrid_shape(bg):
    H = np.where(bg >= 105.0, 272 - 3.121 * bg + 0.01511 * bg**2 - 3.305e-5 * bg**3 + 2.69e-8 * bg**4,
                 75.8 * (105.0 / bg) ** 3.5)
    H100 = 75.8 * (105.0 / 100.0) ** 3.5
    return H / H100


def user_cv_err(y, z, seed=0):
    """Within-user 5-fold OOF error for actual ≈ a + s·z (z = pred_drop·shape)."""
    err = np.full(len(y), np.nan)
    kf = KFold(5, shuffle=True, random_state=seed)
    for tr, te in kf.split(z):
        X = np.column_stack([np.ones(len(tr)), z[tr]])
        beta, *_ = np.linalg.lstsq(X, y[tr], rcond=None)
        err[te] = (beta[0] + beta[1] * z[te]) - y[te]
    return err


def agg_mae(d, zfn):
    """median across users of within-user-CV median|err|, for feature z = pred_drop·zfn(bg)."""
    per = []
    for _, g in d.groupby("user"):
        z = g.pred_drop.values * zfn(g.bg.values)
        e = user_cv_err(g.actual_drop.values, z)
        per.append(np.median(np.abs(e)))
    return float(np.median(per)), np.array(per)


def main():
    d = ec.load_windows().copy()
    d["actual_drop"] = d.bg - d.bg_end
    d["pred_drop"] = d.err_static + d.actual_drop
    d = d[d.pred_drop > 0]
    cnt = d.groupby("user").bg.transform("size")
    d = d[cnt >= MIN_W].reset_index(drop=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)

    # static (no fit) and magnitude (k=0) and the k-sweep, all adaptive within-user CV
    static_per = d.groupby("user").apply(
        lambda g: np.median(np.abs((g.pred_drop - g.actual_drop).values)), include_groups=False)
    static_mae = float(static_per.median())

    curve = {}
    per_user_by_k = {}
    for k in KGRID:
        m, per = agg_mae(d, lambda bg, k=k: (100.0 / bg) ** k)
        curve[float(k)] = round(m, 2)
        per_user_by_k[float(k)] = per
    k_star = min(curve, key=curve.get)
    mag_mae = curve[0.0]

    # diabeloop shape, per-user scaled (adaptive)
    dia_mae, _ = agg_mae(d, hybrid_shape)

    # cold-start: s=1, a=0, population shape swept over k (no per-user adaptation)
    cold = {}
    for k in KGRID:
        e = d.pred_drop.values * (100.0 / d.bg.values) ** k - d.actual_drop.values
        per = pd.DataFrame({"u": d.user.values, "e": np.abs(e)}).groupby("u").e.median()
        cold[float(k)] = round(float(per.median()), 2)
    cold_kstar = min(cold, key=cold.get)

    # per-user optimal k (which steepness each person wants)
    opt_k = []
    for ui, u in enumerate(sorted(d.user.unique())):
        errs = {k: per_user_by_k[k][ui] for k in [float(x) for x in KGRID]}
        opt_k.append(min(errs, key=errs.get))
    opt_k = np.array(opt_k)

    # is the per-user scale s_u predictable from TDD/profile (cold-start deployability)?
    s_rows = []
    for u, g in d.groupby("user"):
        z = g.pred_drop.values * (100.0 / g.bg.values) ** k_star
        X = np.column_stack([np.ones(len(g)), z])
        beta, *_ = np.linalg.lstsq(X, g.actual_drop.values, rcond=None)
        s_rows.append({"user": u, "s": beta[1], "a": beta[0],
                       "tdd": float(np.nanmedian(g.tdd)), "profile_isf": float(np.nanmedian(g.profile_isf))})
    S = pd.DataFrame(s_rows)
    rho_tdd = stats.spearmanr(S.tdd, S.s, nan_policy="omit")
    rho_isf = stats.spearmanr(S.profile_isf, S.s, nan_policy="omit")

    summary = {
        "n_users": int(d.user.nunique()), "n_windows": int(len(d)),
        "static_mae": round(static_mae, 2),
        "magnitude_mae_k0": mag_mae,
        "best_fit_k": k_star, "best_fit_mae": curve[k_star],
        "glucose_adds_over_magnitude": round(mag_mae - curve[k_star], 2),
        "diabeloop_shape_scaled_mae": round(dia_mae, 2),
        "mae_vs_k_adaptive": curve,
        "cold_start_best_k": cold_kstar, "cold_start_best_mae": cold[cold_kstar],
        "cold_start_mae_vs_k": cold,
        "per_user_optimal_k": {"median": float(np.median(opt_k)), "iqr": [float(np.percentile(opt_k, 25)),
                               float(np.percentile(opt_k, 75))], "share_k0": float(np.mean(opt_k == 0)),
                               "share_k_ge1": float(np.mean(opt_k >= 1))},
        "scale_s_predictable": {"spearman_vs_tdd": round(float(rho_tdd.statistic), 2),
                                "p_tdd": round(float(rho_tdd.pvalue), 4),
                                "spearman_vs_profile_isf": round(float(rho_isf.statistic), 2),
                                "p_isf": round(float(rho_isf.pvalue), 4)},
        "read": "k=0 is flat (=magnitude model). best_fit_k>0 with lower MAE ⇒ an individualised "
                "glucose shape helps BEYOND magnitude. cold-start = population shape with no per-user "
                "adaptation. per-user optimal k shows whether people want different steepness.",
    }
    (OUT / "gradient_isf_fit.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    ks = [float(x) for x in KGRID]
    ax[0].plot(ks, [curve[k] for k in ks], "o-", color="#1f77b4", label="adaptive (per-user scale)")
    ax[0].plot(ks, [cold[k] for k in ks], "s--", color="#d62728", label="cold-start (no adaptation)")
    ax[0].axhline(static_mae, color="#999", ls=":", label=f"static {static_mae:.1f}")
    ax[0].axvline(k_star, color="#1f77b4", ls=":", lw=1)
    ax[0].set_xlabel("glucose steepness k  (0 = flat)"); ax[0].set_ylabel("out-of-user MAE (mg/dL)")
    ax[0].set_title(f"Best-fit glucose-dependence\nk*={k_star} (adaptive)"); ax[0].legend(fontsize=8)
    ax[1].hist(opt_k, bins=np.arange(-0.125, 4.2, 0.25), color="#9ecae1", edgecolor="#3182bd")
    ax[1].axvline(float(np.median(opt_k)), color="#d62728", lw=2, label=f"median {np.median(opt_k):.2f}")
    ax[1].set_xlabel("per-user optimal k"); ax[1].set_ylabel("users")
    ax[1].set_title("Steepness each person 'wants'"); ax[1].legend(fontsize=8)
    ax[2].scatter(S.tdd, S.s, s=22, color="#2ca02c", alpha=0.7)
    ax[2].set_xlabel("user median TDD (U/day)"); ax[2].set_ylabel(f"fitted per-user scale s (at k={k_star})")
    ax[2].set_title(f"Is the scale cold-start-predictable?\nρ(TDD)={rho_tdd.statistic:.2f}")
    ax[2].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_gradient_isf_fit.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Best-fit individualised glucose-ISF (gradient/separable-NLLS over a shared shape)\n",
          f"{summary['n_users']} users, {len(d):,} windows. Model: actual_drop ≈ a_u + s_u·(pred_drop·"
          "(100/BG)^k). Per-user (a_u,s_u) fit by within-user 5-fold CV; shared k searched. **k=0 = "
          "flat magnitude model; does k>0 help?**\n",
          "## Headline\n", "| model | out-of-user MAE |", "|---|---|",
          f"| static (no fit) | {static_mae:.2f} |",
          f"| magnitude (k=0, per-user scale) | {mag_mae} |",
          f"| **best-fit glucose k*={k_star} (per-user scale)** | **{curve[k_star]}** |",
          f"| diabeloop shape, per-user scaled | {round(dia_mae,2)} |",
          f"| cold-start best (no adaptation, k={cold_kstar}) | {cold[cold_kstar]} |",
          f"\n**Glucose shape adds {summary['glucose_adds_over_magnitude']} mg/dL beyond magnitude** "
          f"(k*={k_star}). Per-user optimal k: median {np.median(opt_k):.2f}, "
          f"{summary['per_user_optimal_k']['share_k0']:.0%} want k=0, "
          f"{summary['per_user_optimal_k']['share_k_ge1']:.0%} want k≥1.\n",
          f"Per-user scale s predictable from TDD: ρ={rho_tdd.statistic:.2f} (p={rho_tdd.pvalue:.3f}); "
          f"from profile ISF: ρ={rho_isf.statistic:.2f} (p={rho_isf.pvalue:.3f}).\n",
          "![gradient fit](charts/inv008/fig_gradient_isf_fit.png)\n", "*" + summary["read"] + "*"]
    (OUT / "gradient_isf_fit.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
