#!/usr/bin/env python3
"""Bridging experiment: do the Diabeloop / power-law glucose SHAPES transfer to the oref cohort,
and does any benefit survive controlling for correction magnitude (pred_drop)?

The Diabeloop curve (n≈1000s) and the N=1 power-law are ALGORITHM-FIT, so their glucose steepness
conflates real glucose-sensitivity with the IOB model's over-trust of large corrections (bg and
pred_drop are coupled, corr≈0.59). This script adjudicates, parquet-only (no DB, no Pool):

Every candidate ISF is applied as a SHAPE s(BG) anchored to each user's profile ISF at 100 mg/dL,
so candidates differ only in curvature (the level — the √TDD question — is held at profile). Because
the loop's prediction is linear in ISF, the per-window error is simply:

    pred_drop_static = err_static + (bg − bg_end)          # = activity-integral × profile_ISF
    err_candidate    = pred_drop_static · s(bg) − actual_drop      actual_drop = bg − bg_end

Shapes (normalised s(100)=1):
    static     1
    power-law  (100/BG)^k     k = 2, 3, 3.5
    diabeloop  H(BG)/H(100),  H = quartic (≥105) / 75.8·(105/BG)^3.5 (<105)   [the hybrid]

Controls, fit leave-one-USER-out (GroupKFold), glucose-BLIND vs glucose-aware optimal multipliers
on the predicted drop (m·pred_drop ≈ actual_drop), the decisive test:
    magnitude-only  m*(pred_drop)        — corrects the over-trust of big corrections, no glucose
    bg-only         m*(bg)               — empirical CEILING for any glucose-shaped ISF
    bg+pred_drop    m*(bg, pred_drop)     — does bg add anything beyond magnitude?

Verdict: a glucose shape is transferable physiology only if it beats static AND beats the
magnitude-only control. If it beats static but not magnitude-only, its gain was algorithm-magnitude
compensation → fix the IOB/activity model, don't bend ISF by glucose.
Output: results/bridge_diabeloop.{json,md}, charts/inv008/fig_bridge_diabeloop.png
Run: python -m inv008.bridge_diabeloop
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
ANCHOR = 100.0
N_FOLDS = 5


def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4


def hybrid_isf(g):
    g = np.asarray(g, float)
    return np.where(g >= 105.0, quartic(g), 75.8 * (105.0 / g) ** 3.5)


SHAPES = {
    "static": lambda bg: np.ones_like(bg),
    "powerlaw_k2": lambda bg: (ANCHOR / bg) ** 2,
    "powerlaw_k3": lambda bg: (ANCHOR / bg) ** 3,
    "powerlaw_k3.5": lambda bg: (ANCHOR / bg) ** 3.5,
    "diabeloop_hybrid": lambda bg: hybrid_isf(bg) / hybrid_isf(np.array([ANCHOR]))[0],
}


def louo_multiplier(d, cols, nbins=10):
    """Leave-one-user-out optimal multiplier m on pred_drop, binned on `cols`.

    For each held-out fold we learn, from the training users, the median realised/predicted ratio
    in each bin of `cols`, then apply it to the held-out users. This is the best a glucose-blind
    (cols=['pred_drop']) or glucose-aware (cols incl. 'bg') correction of the predicted drop can do
    out-of-user. Returns per-window err = m·pred_drop − actual_drop.
    """
    r = (d.actual_drop / d.pred_drop).clip(-1.0, 3.0).values
    err = np.full(len(d), np.nan)
    gkf = GroupKFold(n_splits=N_FOLDS)
    Xb = {c: d[c].values for c in cols}
    for tr, te in gkf.split(d, groups=d.user.values):
        # bin edges from training fold (quantiles per column)
        edges = {c: np.unique(np.quantile(Xb[c][tr], np.linspace(0, 1, nbins + 1))) for c in cols}
        def keys(idx):
            return tuple(np.clip(np.digitize(Xb[c][idx], edges[c][1:-1]), 0, nbins) for c in cols)
        tr_keys = np.stack(keys(tr), axis=1)
        te_keys = np.stack(keys(te), axis=1)
        # median ratio per training bin
        tab = {}
        dfk = pd.DataFrame(tr_keys); dfk["r"] = r[tr]
        for key, grp in dfk.groupby(list(dfk.columns[:-1])):
            tab[key if isinstance(key, tuple) else (key,)] = float(np.median(grp.r))
        glob = float(np.median(r[tr]))
        m = np.array([tab.get(tuple(k), glob) for k in te_keys])
        err[te] = m * d.pred_drop.values[te] - d.actual_drop.values[te]
    return err


def score(err, d):
    """Out-of-user outcome: per-user median |err|, then median across users (equal user weight)."""
    s = pd.DataFrame({"u": d.user.values, "e": np.abs(err), "b": err})
    per_user_mae = s.groupby("u").e.median()
    per_user_bias = s.groupby("u").b.median()
    return {"mae": round(float(per_user_mae.median()), 2),
            "bias": round(float(per_user_bias.median()), 2),
            "n_users": int(per_user_mae.notna().sum())}


def main():
    d = ec.load_windows().copy()
    d["actual_drop"] = d.bg - d.bg_end
    d["pred_drop"] = d.err_static + d.actual_drop          # = activity-integral × profile ISF
    d = d[d.pred_drop > 0].reset_index(drop=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)

    res = {}
    errs = {}
    for name, s in SHAPES.items():
        e = d.pred_drop.values * s(d.bg.values) - d.actual_drop.values
        errs[name] = e
        res[name] = score(e, d)
    for name, cols in [("ctrl_magnitude_only", ["pred_drop"]),
                       ("ctrl_bg_only", ["bg"]),
                       ("ctrl_bg+pred_drop", ["bg", "pred_drop"])]:
        e = louo_multiplier(d, cols)
        errs[name] = e
        res[name] = score(e, d)

    static_mae = res["static"]["mae"]
    mag_mae = res["ctrl_magnitude_only"]["mae"]
    verdict = {}
    for shape in ["powerlaw_k2", "powerlaw_k3", "powerlaw_k3.5", "diabeloop_hybrid"]:
        m = res[shape]["mae"]
        verdict[shape] = {
            "beats_static": m < static_mae,
            "beats_magnitude_only": m < mag_mae,
            "vs_static": round(static_mae - m, 2),
            "vs_magnitude_only": round(mag_mae - m, 2),
        }

    # per-BG-band MAE for the key candidates (where do shapes help/hurt?)
    bidx = ec.band_of(d.bg.values)
    by_band = {}
    for shape in ["static", "powerlaw_k3.5", "diabeloop_hybrid", "ctrl_magnitude_only"]:
        by_band[shape] = {}
        for k in range(len(ec.BG_BANDS)):
            m = bidx == k
            by_band[shape][ec.BAND_LBL[k]] = round(float(np.median(np.abs(errs[shape][m]))), 1) if m.sum() else None

    summary = {
        "n_users": int(d.user.nunique()), "n_windows": int(len(d)),
        "design": "candidate ISF applied as a shape anchored to each user's profile ISF at 100 "
                  "mg/dL; err = pred_drop·shape(bg) − actual_drop; out-of-user MAE = per-user median "
                  "|err|, median across users. Controls are LOUO optimal multipliers on pred_drop.",
        "scores": res,
        "verdict_vs_controls": verdict,
        "mae_by_bg_band": by_band,
        "read": "A glucose shape is transferable physiology only if it beats BOTH static and the "
                "magnitude-only control. Beats static but not magnitude-only ⇒ the gain is "
                "correction-magnitude compensation (fix the IOB model), not a glucose effect. "
                "ctrl_bg_only is the empirical ceiling for any glucose-shaped ISF; "
                "ctrl_bg+pred_drop − ctrl_magnitude_only = what bg adds beyond magnitude.",
    }
    (OUT / "bridge_diabeloop.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    order = ["static", "powerlaw_k2", "powerlaw_k3", "powerlaw_k3.5", "diabeloop_hybrid",
             "ctrl_magnitude_only", "ctrl_bg_only", "ctrl_bg+pred_drop"]
    cols = ["#999", "#fdae6b", "#fd8d3c", "#e6550d", "#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
    ax[0].bar(range(len(order)), [res[o]["mae"] for o in order], color=cols)
    ax[0].axhline(static_mae, color="#999", ls="--", lw=1)
    ax[0].axhline(mag_mae, color="#1f77b4", ls="--", lw=1, label="magnitude-only control")
    ax[0].set_xticks(range(len(order))); ax[0].set_xticklabels(order, rotation=40, ha="right", fontsize=8)
    ax[0].set_ylabel("out-of-user MAE (mg/dL)"); ax[0].legend(fontsize=8)
    ax[0].set_title("Candidate ISF shapes vs controls\n(lower = better; must beat the blue line)")
    x = ec.BAND_CTR
    for shape, c in [("static", "#999"), ("powerlaw_k3.5", "#e6550d"),
                     ("diabeloop_hybrid", "#d62728"), ("ctrl_magnitude_only", "#1f77b4")]:
        ax[1].plot(x, [by_band[shape][lbl] for lbl in ec.BAND_LBL], "o-", color=c, lw=2, label=shape)
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("median |err| (mg/dL)")
    ax[1].set_title("Where shapes help or hurt, by glucose"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_bridge_diabeloop.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Bridging experiment: do the Diabeloop / power-law glucose shapes transfer to oref?\n",
          f"{summary['n_users']} users, {len(d):,} windows. Each shape anchored to each user's profile "
          "ISF at 100 mg/dL (curvature only). Out-of-user MAE = per-user median |err|, median across "
          "users. **A glucose shape is real, transferable physiology only if it beats BOTH `static` "
          "and `ctrl_magnitude_only`.**\n",
          "## Out-of-user scores (mg/dL)\n",
          "| candidate | MAE | bias | vs static | vs magnitude-only |", "|---|---|---|---|---|"]
    for o in order:
        v = verdict.get(o, {})
        vs_s = f"{v['vs_static']:+}" if v else "—"
        vs_m = f"{v['vs_magnitude_only']:+}" if v else "—"
        md.append(f"| {o} | {res[o]['mae']} | {res[o]['bias']} | {vs_s} | {vs_m} |")
    md += ["\n## MAE by glucose band (mg/dL)\n",
           "| candidate | " + " | ".join(ec.BAND_LBL) + " |",
           "|---|" + "---|" * len(ec.BAND_LBL)]
    for shape in ["static", "powerlaw_k3.5", "diabeloop_hybrid", "ctrl_magnitude_only"]:
        md.append(f"| {shape} | " + " | ".join(str(by_band[shape][l]) for l in ec.BAND_LBL) + " |")
    md += ["\n![bridge](charts/inv008/fig_bridge_diabeloop.png)\n", "*" + summary["read"] + "*"]
    (OUT / "bridge_diabeloop.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
