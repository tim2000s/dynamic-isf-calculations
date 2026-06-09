#!/usr/bin/env python3
"""Is the static-ISF error-vs-glucose pattern consistent across people, or idiosyncratic?

err_curve.py draws the average curve; this script asks whether it is a shared law or a pooled
average of disagreeing individuals. Per user we fit two shapes and pool them with a random-effects
meta-analysis (DerSimonian-Laird), which separates the genuine population effect from the
between-user spread τ:

  monotone gradient   frac_overpred ~ BG            slope < 0 ⇒ over-prediction eases as BG rises
                                                    (the low-glucose effective-sensitivity loss)
  U-shape / high-BG   err_static  ~ BG + BG²         curvature > 0 ⇒ the high-glucose up-turn

Each is fitted raw and again with confound controls (start_slope, IOB, late-hour) so we can see
how much of each effect is really glucose rather than dawn/carbs. Finally we test what predicts a
person's high-BG over-prediction (their TDD, profile ISF) to explain the spread.
Output: results/err_consistency.{json,md}, charts/inv008/fig_err_consistency.png
Run: python -m inv008.err_consistency
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
MIN_FIT = 25            # min windows for a stable per-user regression


def per_user_fits(d: pd.DataFrame):
    """Per-user coefficients for the monotone (fractional) and curvature (absolute) shapes.

    Returns a DataFrame with, per user: fractional BG-slope and its SE (raw + adjusted), and the
    absolute BG² curvature and its SE (raw + adjusted), plus user-level moderators.
    """
    recs = []
    for u, g in d.groupby("user"):
        if len(g) < MIN_FIT:
            continue
        bg = g.bg.values
        bgc = (bg - bg.mean()) / 100.0                 # per 100 mg/dL, centred
        z = np.column_stack([g.start_slope.values, g.iob.values, g.late_hour.values])
        rec = {"user": u, "n": len(g), "tdd": float(np.nanmedian(g.tdd)),
               "profile_isf": float(np.nanmedian(g.profile_isf)),
               "hi_overpred": float(np.nanmedian(g.err_static[bg >= 175])) if (bg >= 175).sum() >= 10 else np.nan}

        # monotone: fractional ~ BG   (slope, index 1)
        for tag, X in [("frac_slope", np.column_stack([np.ones(len(g)), bgc])),
                       ("frac_slope_adj", np.column_stack([np.ones(len(g)), bgc, z]))]:
            f = ec.ols(g.frac_overpred.values, X)
            rec[tag], rec[tag + "_se"] = (f[0][1], f[1][1]) if f else (np.nan, np.nan)

        # curvature: err_static ~ BG + BG²   (quadratic coef, index 2)
        for tag, X in [("curv", np.column_stack([np.ones(len(g)), bgc, bgc**2])),
                       ("curv_adj", np.column_stack([np.ones(len(g)), bgc, bgc**2, z]))]:
            f = ec.ols(g.err_static.values, X)
            rec[tag], rec[tag + "_se"] = (f[0][2], f[1][2]) if f else (np.nan, np.nan)
        recs.append(rec)
    return pd.DataFrame(recs)


def pool(F: pd.DataFrame, coef: str):
    r = ec.dersimonian_laird(F[coef].values, F[coef + "_se"].values)
    if r is None:
        return None
    return {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}


def main():
    d = ec.load_windows()
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    F = per_user_fits(d)

    pooled = {
        "frac_slope_raw": pool(F, "frac_slope"),
        "frac_slope_adjusted": pool(F, "frac_slope_adj"),
        "curvature_raw": pool(F, "curv"),
        "curvature_adjusted": pool(F, "curv_adj"),
    }

    # what predicts a person's high-BG over-prediction?
    moder = {}
    sub = F.dropna(subset=["hi_overpred"])
    for m in ["tdd", "profile_isf"]:
        s = sub.dropna(subset=[m])
        rho, p = stats.spearmanr(s[m], s.hi_overpred)
        moder[m] = {"spearman_rho": round(float(rho), 3), "p": round(float(p), 4), "n": int(len(s))}

    summary = {
        "n_users_fitted": int(len(F)), "n_windows": int(len(d)),
        "interpretation": {
            "frac_slope": "fractional over-prediction vs BG; NEGATIVE = static ISF over-predicts "
                          "the drop more at LOW glucose (effective sensitivity falls near normal).",
            "curvature": "quadratic coef of err_static vs BG; POSITIVE = U-shaped, i.e. a high-"
                         "glucose up-turn in the absolute error on top of the low-glucose effect.",
            "pooling": "random-effects (DerSimonian-Laird): b_re=population effect, tau=between-"
                       "user SD, I2=% variance between users, frac_same_sign=share of users "
                       "agreeing with the pooled sign, p=is the population effect ≠ 0.",
        },
        "pooled_effects": pooled,
        "spread_moderators_vs_high_bg_overpred": moder,
    }
    (OUT / "err_consistency.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    # left: per-user fractional BG-slope distribution + pooled estimate
    sl = F.frac_slope.dropna().values
    ax[0].hist(sl, bins=22, color="#9ecae1", edgecolor="#3182bd")
    pr = pooled["frac_slope_raw"]
    ax[0].axvline(0, color="k", lw=1, ls="--")
    if pr:
        ax[0].axvline(pr["b_re"], color="#d62728", lw=2.5,
                      label=f"pooled {pr['b_re']:+.2f}  (p={pr['p']:.1e}, I²={pr['I2_pct']:.0f}%)")
    ax[0].set_xlabel("per-user fractional slope vs BG (per 100 mg/dL)")
    ax[0].set_ylabel("users"); ax[0].legend(fontsize=8)
    ax[0].set_title("Effective-sensitivity gradient by user\n(<0: worse over-prediction at low glucose)")
    # right: per-user high-BG over-prediction vs TDD
    s = sub.dropna(subset=["tdd"])
    ax[1].scatter(s.tdd, s.hi_overpred, s=22, color="#2ca02c", alpha=0.7)
    ax[1].axhline(0, color="k", lw=1, ls="--")
    rho = moder.get("tdd", {}).get("spearman_rho")
    ax[1].set_xlabel("user median TDD (U/day)")
    ax[1].set_ylabel("median err_static at BG≥175 (mg/dL)")
    ax[1].set_title(f"High-glucose over-prediction vs TDD\n(Spearman ρ={rho})")
    ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_err_consistency.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    def row(name, r):
        if not r:
            return f"| {name} | – | – | – | – | – |"
        return (f"| {name} | {r['b_re']:+.3f} | {r['p']:.1e} | {r['tau']:.3f} | "
                f"{r['I2_pct']:.0f}% | {r['frac_same_sign']:.0%} |")

    md = ["# Is the static-ISF error pattern consistent across users?\n",
          f"{len(F)} users fitted, {len(d):,} windows. Per-user shapes pooled by random-effects "
          "meta-analysis. `b_re` = population effect, `τ` = between-user SD, `I²` = share of "
          "variance that is between-user, `same-sign` = users agreeing with the pooled direction.\n",
          "## Pooled effects\n",
          "| effect | b_re | p | τ | I² | same-sign |", "|---|---|---|---|---|---|",
          row("fractional slope vs BG (raw)", pooled["frac_slope_raw"]),
          row("fractional slope vs BG (confound-adj)", pooled["frac_slope_adjusted"]),
          row("curvature err~BG² (raw)", pooled["curvature_raw"]),
          row("curvature err~BG² (confound-adj)", pooled["curvature_adjusted"]),
          "\n*Negative fractional slope ⇒ over-prediction is worse at LOW glucose (effective "
          "sensitivity falls near normal). Positive curvature ⇒ a high-glucose up-turn on top. "
          "High I² / non-trivial τ ⇒ the effect, even if real on average, varies a lot between "
          "people; `same-sign` is the cleanest read on consistency.*\n",
          "## What predicts a person's high-glucose over-prediction?\n",
          "| moderator | Spearman ρ | p | n |", "|---|---|---|---|"]
    for m, r in moder.items():
        md.append(f"| {m} | {r['spearman_rho']} | {r['p']} | {r['n']} |")
    md += ["\n![Consistency](charts/inv008/fig_err_consistency.png)\n"]
    (OUT / "err_consistency.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
