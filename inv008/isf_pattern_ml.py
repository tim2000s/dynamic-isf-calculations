#!/usr/bin/env python3
"""Joint-control test: does a glucose→ISF effect survive once the confounds are in the model?

The univariate glucose gradients in this data are confound mirages — each observational metric
manufactures a different sign (min-drop selection near target; mean-reversion / endogenous
momentum at high BG). This script does the adjudication: a gradient-boosted model predicts the
static-ISF error from glucose AND the physiological confounders together, under leave-one-USER-out
CV, and we ask whether glucose keeps predictive weight once the model can already see momentum,
IOB and the predicted-drop magnitude.

Physiological reading of the features (the realised fasting drop = insulin action − endogenous
glucose appearance + insulin-independent/renal clearance, all gated by counterregulation):

    bg          glucose level — the question. Net of the others: mass-action/renal clearance
                (raises sensitivity with BG) vs glucotoxic resistance (lowers it) vs
                counterregulatory blunting near target.
    pre_slope   30-min momentum BEFORE the window — proxy for endogenous flux the loop dosed
                into (dawn / hepatic output / carb tail), ~uncorrelated with BG (the clean one).
    momentum_endo  pre_slope orthogonalised to IOB within user — the part of entry momentum NOT
                explained by insulin on board → closer to purely endogenous drive.
    start_slope post-decision 15-min slope — contains the insulin-driven fall (a mediator; kept
                only to show how much it over-controls vs pre_slope).
    iob         insulin on board — dose magnitude / pharmacodynamics.
    pred_drop   the drop the loop predicted — the selection variable that flipped the low-BG sign.
    hour        23→02; proxy for the dawn phenomenon (hepatic glucose output rises pre-dawn).
    tdd         near-constant within user (we centre within user) → expected to add ~nothing here.

Targets, both centred WITHIN user so the model can only earn credit for shared, actionable
structure (the per-user level — the thing static ISF already captures — is removed):
    y_err   err_static − user median   (bounded mg/dL; outcome-relevant)
    y_frac  frac_overpred − user median (effective-sensitivity gap; physiological, de-scaled)

Decision rule: ship a glucose term only if FULL beats NO_BG on out-of-user error AND the SHAP
curve of bg is physiologically coherent. Output: results/isf_pattern_ml.{json,md},
charts/inv008/fig_isf_pattern_ml.png.  Run: python -m inv008.isf_pattern_ml
"""
from __future__ import annotations

import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.model_selection import GroupKFold

warnings.filterwarnings("ignore")
from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"

FEATURES = ["bg", "pre_slope", "momentum_endo", "start_slope", "iob", "pred_drop", "hour", "tdd"]
N_FOLDS = 5
LGB_PARAMS = dict(objective="huber", alpha=0.9, n_estimators=400, learning_rate=0.03,
                  num_leaves=31, min_child_samples=200, subsample=0.8, subsample_freq=1,
                  colsample_bytree=0.8, reg_lambda=1.0, verbosity=-1)


def prepare(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["pred_drop"] = d.err_static + (d.bg - d.bg_end)
    # endogenous momentum: pre_slope with the IOB-explained part removed, per user
    d["momentum_endo"] = np.nan
    for u, g in d.groupby("user"):
        m = g.pre_slope.notna() & g.iob.notna()
        if m.sum() >= 20:
            X = np.column_stack([np.ones(m.sum()), g.iob[m].values])
            beta, *_ = np.linalg.lstsq(X, g.pre_slope[m].values, rcond=None)
            d.loc[g.index[m], "momentum_endo"] = g.pre_slope[m].values - X @ beta
    # within-user-centred targets (per-user level removed)
    d["y_err"] = d.err_static - d.groupby("user").err_static.transform("median")
    d["y_frac"] = d.frac_overpred - d.groupby("user").frac_overpred.transform("median")
    d["w"] = 1.0 / d.groupby("user").bg.transform("size")          # equal weight per user
    return d


def oof(d: pd.DataFrame, target: str, feats: list[str]):
    """Leave-one-user-out (grouped) out-of-fold predictions for a feature set."""
    X = d[feats].values
    y = d[target].values
    g = d.user.values
    w = d.w.values
    pred = np.full(len(d), np.nan)
    shap = np.zeros((len(d), len(feats)))
    gkf = GroupKFold(n_splits=N_FOLDS)
    for tr, te in gkf.split(X, y, groups=g):
        m = lgb.LGBMRegressor(**LGB_PARAMS)
        m.fit(X[tr], y[tr], sample_weight=w[tr])
        pred[te] = m.predict(X[te])
        contrib = m.predict(X[te], pred_contrib=True)          # [n, nfeat+1] incl. base value
        shap[te] = contrib[:, :len(feats)]
    return pred, shap


def wmae(y, yhat, w):
    e = np.abs(y - yhat)
    ok = np.isfinite(e) & np.isfinite(w)
    return float(np.median(e[ok]))                              # median is already robust; w≈equal within band


def main():
    d = prepare(ec.load_windows())
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    results = {}

    for target, lbl in [("y_err", "err_static (mg/dL, outcome)"),
                        ("y_frac", "frac_overpred (effective-sensitivity)")]:
        dd = d[d[target].notna()].copy()
        y, w = dd[target].values, dd.w.values
        base = wmae(y, np.zeros(len(y)), w)                     # predict the per-user mean (=0)
        full_pred, full_shap = oof(dd, target, FEATURES)
        nobg_pred, _ = oof(dd, target, [f for f in FEATURES if f != "bg"])
        bgonly_pred, _ = oof(dd, target, ["bg"])
        # SHAP curve of bg (its learned effect net of the other features), by BG band
        bidx = ec.band_of(dd.bg.values)
        shap_bg = full_shap[:, FEATURES.index("bg")]
        curve = []
        for k in range(len(ec.BG_BANDS)):
            v = shap_bg[bidx == k]
            curve.append({"band": ec.BAND_LBL[k], "n": int(len(v)),
                          "mean_shap_bg": round(float(np.mean(v)), 2) if len(v) else None})
        imp = {f: round(float(np.mean(np.abs(full_shap[:, i]))), 3) for i, f in enumerate(FEATURES)}
        results[target] = {
            "label": lbl, "n": int(len(dd)),
            "baseline_mae_predict_user_mean": round(base, 2),
            "full_mae": round(wmae(y, full_pred, w), 2),
            "no_bg_mae": round(wmae(y, nobg_pred, w), 2),
            "bg_only_mae": round(wmae(y, bgonly_pred, w), 2),
            "bg_adds_vs_confounds": round(wmae(y, nobg_pred, w) - wmae(y, full_pred, w), 3),
            "full_improves_vs_baseline": round(base - wmae(y, full_pred, w), 3),
            "mean_abs_shap": dict(sorted(imp.items(), key=lambda kv: -kv[1])),
            "bg_shap_curve": curve,
        }

    (OUT / "isf_pattern_ml.json").write_text(json.dumps(results, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    r = results["y_frac"]
    feats = list(r["mean_abs_shap"]); vals = [r["mean_abs_shap"][f] for f in feats]
    ax[0].barh(feats[::-1], vals[::-1], color="#4c78a8")
    ax[0].set_xlabel("mean |SHAP|  (effective-sensitivity target)")
    ax[0].set_title("Feature importance, jointly controlled\n(does bg survive beside momentum/IOB?)")
    x = ec.BAND_CTR
    for tgt, col, lab in [("y_frac", "#d62728", "effective-sensitivity"),
                          ("y_err", "#1f77b4", "err mg/dL (scaled)")]:
        c = results[tgt]["bg_shap_curve"]
        ax[1].plot(x, [p["mean_shap_bg"] for p in c], "o-", color=col, lw=2, label=lab)
    ax[1].axhline(0, color="k", lw=1, ls="--")
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("mean SHAP of bg (learned effect)")
    ax[1].set_title("Glucose effect on ISF error, net of confounds\n(SHAP curve of bg)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_isf_pattern_ml.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Does a glucose→ISF effect survive joint control? (leave-one-user-out GBM)\n",
          f"{results['y_err']['n']:,} windows. Target centred within user, so only *shared* "
          "within-user structure can be learned. Lower MAE = better. Decision: a glucose term is "
          "justified only if **full beats no-bg** and the bg SHAP curve is physiologically coherent.\n"]
    for tgt in ["y_err", "y_frac"]:
        r = results[tgt]
        md += [f"## Target: {r['label']}\n",
               "| model | OOF MAE |", "|---|---|",
               f"| baseline (per-user mean) | {r['baseline_mae_predict_user_mean']} |",
               f"| bg only | {r['bg_only_mae']} |",
               f"| confounds, no bg | {r['no_bg_mae']} |",
               f"| **full (bg + confounds)** | **{r['full_mae']}** |",
               f"\n*bg adds vs confounds: **{r['bg_adds_vs_confounds']}** · full vs baseline: "
               f"{r['full_improves_vs_baseline']} (MAE units).*\n",
               "mean |SHAP|: " + ", ".join(f"{k} {v}" for k, v in r["mean_abs_shap"].items()) + "\n",
               "bg SHAP by glucose: " + ", ".join(
                   f"{p['band']}:{p['mean_shap_bg']}" for p in r["bg_shap_curve"]) + "\n"]
    md += ["![pattern](charts/inv008/fig_isf_pattern_ml.png)\n"]
    (OUT / "isf_pattern_ml.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
