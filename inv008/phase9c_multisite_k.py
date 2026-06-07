#!/usr/bin/env python3
"""Phase 9c: fit the glucose exponent k per site at full insulin action, across the
12-site multisite 4h cache — to firm up the single-patient k≈2.25 and see its spread,
and to read the TDD exponent across sites (TDD 5.6-83.8 U/day).

Prediction-error design per site: ISF_cand = α·(target/BG)^k; scale the loop's predicted
drop by ISF_cand/ISF_loop; compare to actual_bg_end; α and k fit to minimise MAE.

Output: results/phase9c_multisite_k.{json,md}
Run: python -m inv008.phase9c_multisite_k
"""
from __future__ import annotations

import json
import math
import os
import pickle
from pathlib import Path

import numpy as np

CACHE = Path(os.environ.get("MULTISITE_4H_CACHE",
    "/Users/timstreet/Library/CloudStorage/GoogleDrive-tim.street@liveintheirshoes.com/"
    "My Drive/Dynamic ISF data/multisite_4h_sample_cache.pkl"))
OUT = Path(os.environ.get("DYNISF_ROOT", Path.cwd())) / "results"
TARGET, D = 99.0, 75.0
KGRID = np.arange(0.5, 6.01, 0.25)


def fit_site(bg, isf_loop, pred_drop, actual_end):
    m = (np.isfinite(bg) & np.isfinite(isf_loop) & np.isfinite(pred_drop)
         & np.isfinite(actual_end) & (isf_loop > 0) & (bg > 0) & (np.abs(pred_drop) >= 3))
    bg, isf_loop, pred_drop, actual_end = bg[m], isf_loop[m], pred_drop[m], actual_end[m]
    if len(bg) < 30:
        return None

    def score(isf):
        return float(np.abs(actual_end - (bg - pred_drop * (isf / isf_loop))).mean())

    def fit_scaled(base):
        a0 = np.median(isf_loop) / np.median(base)
        best = None
        for a in a0 * np.linspace(0.3, 2.5, 121):
            mae = score(a * base)
            if best is None or mae < best[0]:
                best = (mae, a)
        return best

    # power-law: sweep k
    bestpl = None
    for k in KGRID:
        mae, a = fit_scaled((TARGET / bg) ** k)
        if bestpl is None or mae < bestpl[0]:
            bestpl = (mae, k, a)
    # log scaler baseline
    logt = math.log(TARGET / D + 1.0) / np.log(bg / D + 1.0)
    mae_log, _ = fit_scaled(logt)
    mae_loop = score(isf_loop)
    return {"n": int(len(bg)), "best_k": float(bestpl[1]), "isf_at_target": float(bestpl[2]),
            "mae_pl": round(bestpl[0], 2), "mae_log": round(mae_log, 2), "mae_loop": round(mae_loop, 2)}


def main():
    sites = pickle.load(open(CACHE, "rb"))
    rows = []
    for s in sites:
        r = fit_site(np.asarray(s["bg"], float), np.asarray(s["isf_actual"], float),
                     np.asarray(s["pred_drop"], float), np.asarray(s["actual_bg_end"], float))
        if r:
            r.update(name=s["name"], model=s["model"], tdd=float(s["tdd_median"]),
                     horizon_min=float(s.get("median_horizon_min", float("nan"))))
            rows.append(r)

    ks = np.array([r["best_k"] for r in rows])
    pl_wins = sum(r["mae_pl"] <= r["mae_log"] for r in rows)
    # between-site TDD exponent: log(ISF_at_target) vs log(TDD)
    tdd = np.array([r["tdd"] for r in rows]); anch = np.array([r["isf_at_target"] for r in rows])
    good = (tdd > 0) & (anch > 0)
    p_slope = float(np.polyfit(np.log(tdd[good]), np.log(anch[good]), 1)[0])

    summary = {
        "n_sites": len(rows),
        "best_k_median": round(float(np.median(ks)), 2),
        "best_k_IQR": [round(float(np.quantile(ks, .25)), 2), round(float(np.quantile(ks, .75)), 2)],
        "best_k_range": [round(float(ks.min()), 2), round(float(ks.max()), 2)],
        "n_weighted_k": round(float(np.average(ks, weights=[r["n"] for r in rows])), 2),
        "powerlaw_beats_or_ties_log_sites": f"{pl_wins}/{len(rows)}",
        "tdd_exponent_across_sites": round(-p_slope, 3),
        "tdd_exponent_closer_to": "0.5 (√TDD)" if abs(-p_slope - 0.5) < abs(-p_slope - 1.0) else "1.0 (1/TDD)",
        "sites": rows,
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "phase9c_multisite_k.json").write_text(json.dumps(summary, indent=1))

    md = ["# Phase 9c — glucose exponent k across 12 sites (full insulin action)\n",
          f"Per-site prediction-error fit of ISF = α·(target/BG)^k at the end-of-insulin-action "
          f"horizon (~150–225 min). TDD spans {tdd.min():.0f}–{tdd.max():.0f} U/day.\n",
          f"- **median best k = {summary['best_k_median']}** "
          f"[IQR {summary['best_k_IQR'][0]}–{summary['best_k_IQR'][1]}, "
          f"range {summary['best_k_range'][0]}–{summary['best_k_range'][1]}], "
          f"n-weighted {summary['n_weighted_k']}",
          f"- power-law ≥ log for **{summary['powerlaw_beats_or_ties_log_sites']}** sites",
          f"- between-site TDD exponent (ISF∝1/TDD^p): **p = {summary['tdd_exponent_across_sites']}** "
          f"→ closer to **{summary['tdd_exponent_closer_to']}**",
          "\n| site | model | n | TDD | best k | MAE pl | MAE log | MAE loop |",
          "|---|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda x: x["tdd"]):
        md.append(f"| {r['name'][:18]} | {r['model']} | {r['n']} | {r['tdd']:.0f} | "
                  f"**{r['best_k']:.2f}** | {r['mae_pl']} | {r['mae_log']} | {r['mae_loop']} |")
    md.append(f"\n- Single-patient (boost cache) gave k≈2.25 at ~3.17h; the multisite median "
              f"({summary['best_k_median']}) {'agrees' if abs(summary['best_k_median']-2.25)<1 else 'differs'}.")
    md.append("\n*N small per site (40–2819 windows); k is noisy per site but the central "
              "tendency firms up the exponent. Prediction-error design; mixed sigmoid/log loop "
              "formulas (the scaling is formula-agnostic).*")
    (OUT / "phase9c_multisite_k.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
