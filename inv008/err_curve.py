#!/usr/bin/env python3
"""Characterise the static-ISF prediction error as a function of glucose (the full curve).

For every overnight window we know how far the person's static profile ISF was from the ISF
that would have made the loop's drop prediction exact. This script maps that error across the
glucose range, in two views (see err_common):

    absolute    err_static (mg/dL)         — over-prediction the controller feels
    fractional  1 − realised/profile        — effective-sensitivity gap (mechanically de-biased)

It builds the population curve (median across users, bootstrapped over users so one heavy user
cannot dominate), overlays the per-user curves to show spread, and reports a confound-adjusted
curve (regressing out start-slope / IOB / hour within user) so dawn-and-carb windows are not
charged to glucose. Output: results/err_curve.{json,md}, charts/inv008/fig_err_curve.png
Run: python -m inv008.err_curve
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"


def per_user_band_curve(d: pd.DataFrame, value: str):
    """Per-user median of `value` in each BG band → matrix [user x band] (NaN if too few)."""
    users = sorted(d.user.unique())
    M = np.full((len(users), len(ec.BG_BANDS)), np.nan)
    bidx = ec.band_of(d.bg.values)
    for ui, u in enumerate(users):
        gu = d.user.values == u
        for k in range(len(ec.BG_BANDS)):
            v = d[value].values[gu & (bidx == k)]
            if len(v) >= ec.MIN_BAND_WINDOWS:
                M[ui, k] = np.nanmedian(v)
    return users, M


def pooled_curve(M: np.ndarray, seed: int = 1):
    """Population curve = median across users per band, with a bootstrap-over-users 90% band."""
    out = []
    for k in range(M.shape[1]):
        col = M[:, k]
        col = col[np.isfinite(col)]
        med, lo, hi = ec.boot_median_ci(col, seed=seed + k)
        share_pos = float(np.mean(col > 0)) if len(col) else np.nan
        out.append({"band": ec.BAND_LBL[k], "n_users": int(len(col)),
                    "median": round(med, 2), "ci_lo": round(lo, 2), "ci_hi": round(hi, 2),
                    "share_users_overpred": round(share_pos, 2),
                    "iqr_lo": round(float(np.nanpercentile(col, 25)), 2) if len(col) else None,
                    "iqr_hi": round(float(np.nanpercentile(col, 75)), 2) if len(col) else None})
    return out


def confound_adjusted_curve(d: pd.DataFrame):
    """Within-user residual err_static after removing start_slope, IOB, late-hour; then re-band.

    For each user we regress err_static on the confounds (no BG term), and keep the residual +
    the user's mean. What survives is the part of the error that the confounds do NOT explain;
    re-banding it by BG shows the glucose curve net of dawn/carb/IOB structure.
    """
    d = d.copy()
    d["adj"] = np.nan
    for u, g in d.groupby("user"):
        X = np.column_stack([np.ones(len(g)), g.start_slope.values, g.iob.values, g.late_hour.values])
        fit = ec.ols(g.err_static.values, X)
        if fit is None:
            continue
        beta, _, _ = fit
        # residual about the confound model, recentred on the user's mean error
        pred = X @ beta
        d.loc[g.index, "adj"] = g.err_static.values - pred + float(np.nanmean(g.err_static.values))
    _, Madj = per_user_band_curve(d.dropna(subset=["adj"]), "adj")
    return pooled_curve(Madj, seed=50)


def main():
    d = ec.load_windows()
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)

    users_a, Mabs = per_user_band_curve(d, "err_static")
    users_f, Mfrac = per_user_band_curve(d, "frac_overpred")
    abs_curve = pooled_curve(Mabs, seed=1)
    frac_curve = pooled_curve(Mfrac, seed=20)
    adj_curve = confound_adjusted_curve(d)

    summary = {
        "n_users": int(d.user.nunique()), "n_windows": int(len(d)),
        "note": "err_static>0 ⇒ static profile ISF over-predicted the drop. absolute = mg/dL; "
                "fractional = 1−realised/profile (effective-sensitivity gap). Population value = "
                "median across users; CI = bootstrap over users (90%). adjusted = absolute curve "
                "after regressing out start_slope/IOB/late-hour within user.",
        "absolute_mgdl": abs_curve,
        "fractional_pct": [{**r, "median": round(r["median"] * 100, 1),
                            "ci_lo": round(r["ci_lo"] * 100, 1), "ci_hi": round(r["ci_hi"] * 100, 1)}
                           for r in frac_curve],
        "absolute_confound_adjusted_mgdl": adj_curve,
    }
    (OUT / "err_curve.json").write_text(json.dumps(summary, indent=1))

    # ---- figure: 3 panels ----
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    x = ec.BAND_CTR

    def band_plot(a, M, curve, ylab, title, pct=False):
        for row in M:                                   # faint per-user curves
            a.plot(x, row * (100 if pct else 1), color="#bbb", lw=0.6, alpha=0.5)
        med = [r["median"] for r in curve]
        lo = [r["ci_lo"] for r in curve]; hi = [r["ci_hi"] for r in curve]
        a.plot(x, med, "o-", color="#1f77b4", lw=2.5, label="population median")
        a.fill_between(x, lo, hi, color="#1f77b4", alpha=0.18, label="90% CI (over users)")
        a.axhline(0, color="k", lw=1, ls="--")
        a.set_xlabel("glucose (mg/dL)"); a.set_ylabel(ylab); a.set_title(title)
        a.grid(alpha=0.3); a.legend(fontsize=8)

    band_plot(ax[0], Mabs, abs_curve, "median err_static (mg/dL)",
              "Absolute over-prediction vs glucose\n(>0: static ISF over-predicts the drop)")
    band_plot(ax[1], Mfrac * 100, summary["fractional_pct"], "effective-sensitivity gap (%)",
              "Fractional over-prediction vs glucose\n(mechanically de-biased)", pct=True)
    # panel 3: raw vs confound-adjusted absolute curve
    ax[2].plot(x, [r["median"] for r in abs_curve], "o-", color="#1f77b4", lw=2.5, label="raw")
    ax[2].plot(x, [r["median"] for r in adj_curve], "s--", color="#d62728", lw=2.5,
               label="confound-adjusted")
    ax[2].axhline(0, color="k", lw=1, ls="--")
    ax[2].set_xlabel("glucose (mg/dL)"); ax[2].set_ylabel("median err_static (mg/dL)")
    ax[2].set_title("Raw vs confound-adjusted\n(start-slope / IOB / hour removed)")
    ax[2].grid(alpha=0.3); ax[2].legend(fontsize=9)
    fig.tight_layout(); fig.savefig(CHART / "fig_err_curve.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Static-ISF prediction error vs glucose (full curve)\n",
          f"{summary['n_users']} users, {len(d):,} overnight windows. `err_static > 0` means the "
          "person's static profile ISF *over-predicted the drop* (the real fall was smaller). "
          "Population value is the median across users; CI is bootstrapped over users.\n",
          "## Absolute over-prediction (mg/dL)\n",
          "| BG band | median | 90% CI | users over-pred | users |",
          "|---|---|---|---|---|"]
    for r in abs_curve:
        md.append(f"| {r['band']} | {r['median']} | [{r['ci_lo']}, {r['ci_hi']}] | "
                  f"{r['share_users_overpred']:.0%} | {r['n_users']} |")
    md += ["\n## Fractional over-prediction (effective-sensitivity gap, %)\n",
           "| BG band | median | 90% CI | users |", "|---|---|---|---|"]
    for r in summary["fractional_pct"]:
        md.append(f"| {r['band']} | {r['median']}% | [{r['ci_lo']}%, {r['ci_hi']}%] | {r['n_users']} |")
    md += ["\n## Absolute curve, raw vs confound-adjusted (mg/dL)\n",
           "| BG band | raw | adjusted |", "|---|---|---|"]
    for r0, r1 in zip(abs_curve, adj_curve):
        md.append(f"| {r0['band']} | {r0['median']} | {r1['median']} |")
    md += ["\n![Error curve](charts/inv008/fig_err_curve.png)\n",
           "*Absolute and fractional disagree by construction: corrections are bigger at high "
           "glucose, so a flat proportional mismatch shows up as a larger mg/dL error there. The "
           "fractional panel is the physiological view; the confound-adjusted panel shows how much "
           "of the curve survives removing dawn/IOB/hour structure.*"]
    (OUT / "err_curve.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
