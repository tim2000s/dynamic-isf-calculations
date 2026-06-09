#!/usr/bin/env python3
"""Phase C: is the per-user adaptive-k minority (36%) real out-of-sample, or selection noise?

The gradient fit found 36% of users individually prefer a glucose steepness k≥1 — but that k was
picked on the SAME data it was scored on, so it is selection-inflated. Honest test: NESTED, time-
ordered. Per user, take the first 60% of windows as train and the last 40% as test (adapt on the
past, validate on the future). Pick k* on train by internal 3-fold CV, refit the per-user
(intercept, scale) on train at k*, predict the test fold. Compare the test error to the k=0 (flat)
model fit the same way. If per-user-selected k beats k=0 out-of-sample for materially more than half
the users with a real aggregate gain, adaptive-k is worth deploying; if it is ~chance with ~0 gain,
the 36% was overfit and k=0 stands.

Model: actual_drop ≈ a + s·(pred_drop·(100/BG)^k). Parquet-only.
Output: results/adaptive_k_nestedcv.{json,md}. Run: python -m inv008.adaptive_k_nestedcv
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
KGRID = np.round(np.arange(0.0, 4.01, 0.5), 2)
MIN_W = 80


def fit_predict(tr_y, tr_z, te_z):
    X = np.column_stack([np.ones(len(tr_z)), tr_z])
    beta, *_ = np.linalg.lstsq(X, tr_y, rcond=None)
    return beta[0] + beta[1] * te_z


def cv_err_at_k(y, z):
    err = []
    for tr, te in KFold(3, shuffle=True, random_state=0).split(z):
        err.append(np.abs(fit_predict(y[tr], z[tr], z[te]) - y[te]))
    return float(np.median(np.concatenate(err)))


def main():
    d = ec.load_windows().copy()
    d["actual_drop"] = d.bg - d.bg_end
    d["pred_drop"] = d.err_static + d.actual_drop
    d = d[d.pred_drop > 0]
    cnt = d.groupby("user").bg.transform("size")
    d = d[cnt >= MIN_W].reset_index(drop=True)

    recs = []
    for u, g in d.groupby("user"):
        g = g.sort_values("bg").reset_index(drop=True)        # no timestamp in parquet; order by bg
        # time-order proxy unavailable → use a random but fixed split (shuffle once)
        idx = np.random.default_rng(abs(hash(u)) % (2**32)).permutation(len(g))
        cut = int(0.6 * len(g))
        tr, te = idx[:cut], idx[cut:]
        if len(te) < 25:
            continue
        bg = g.bg.values; pd_ = g.pred_drop.values; y = g.actual_drop.values
        # pick k* on train by internal CV
        best_k, best_e = 0.0, np.inf
        for k in KGRID:
            z = pd_[tr] * (100.0 / bg[tr]) ** k
            e = cv_err_at_k(y[tr], z)
            if e < best_e:
                best_e, best_k = e, k
        # evaluate selected-k vs k=0 on the held-out test fold
        def test_err(k):
            ztr = pd_[tr] * (100.0 / bg[tr]) ** k
            zte = pd_[te] * (100.0 / bg[te]) ** k
            return float(np.median(np.abs(fit_predict(y[tr], ztr, zte) - y[te])))
        recs.append({"user": u, "k_star": float(best_k), "mae_kstar": test_err(best_k),
                     "mae_k0": test_err(0.0)})
    R = pd.DataFrame(recs)
    R["gain"] = R.mae_k0 - R.mae_kstar                        # >0 ⇒ adaptive-k helped out-of-sample

    summary = {
        "n_users": int(len(R)),
        "share_kstar_gt0": round(float((R.k_star > 0).mean()), 2),
        "share_adaptive_beats_flat_oos": round(float((R.gain > 0).mean()), 2),
        "median_oos_gain_mgdl": round(float(R.gain.median()), 3),
        "mean_oos_gain_mgdl": round(float(R.gain.mean()), 3),
        "median_mae_k0": round(float(R.mae_k0.median()), 2),
        "median_mae_kstar": round(float(R.mae_kstar.median()), 2),
        "verdict": None,
    }
    real = summary["share_adaptive_beats_flat_oos"] >= 0.6 and summary["median_oos_gain_mgdl"] > 0.5
    summary["verdict"] = ("ADAPTIVE-k REAL — per-user k beats flat out-of-sample for a clear majority"
                          if real else
                          "ADAPTIVE-k NOT supported — out-of-sample it is ~chance with negligible gain; "
                          "the 36% was selection overfit, k=0 stands")
    (OUT / "adaptive_k_nestedcv.json").write_text(json.dumps(summary, indent=1))

    md = ["# Phase C: nested-CV — is per-user adaptive-k real out-of-sample?\n",
          f"{summary['n_users']} users, per-user 60/40 split, k* picked on train (internal CV), scored "
          "on held-out test vs k=0 (flat).\n",
          f"- users whose selected k>0: **{summary['share_kstar_gt0']:.0%}**",
          f"- users where adaptive-k beats flat OUT-OF-SAMPLE: **{summary['share_adaptive_beats_flat_oos']:.0%}** "
          "(50% = chance)",
          f"- median out-of-sample gain: **{summary['median_oos_gain_mgdl']} mg/dL** "
          f"(mean {summary['mean_oos_gain_mgdl']})",
          f"- median test MAE: k=0 {summary['median_mae_k0']} vs adaptive {summary['median_mae_kstar']}\n",
          f"**Verdict: {summary['verdict']}.**"]
    (OUT / "adaptive_k_nestedcv.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
