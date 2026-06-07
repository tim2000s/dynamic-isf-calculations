#!/usr/bin/env python3
"""Phase 9: backtest the consolidated equation  ISF = (C/TDD^p)·(target/BG)^k  on the
detailed N=1 overnight prediction dataset (ns_backtest_overnight.csv).

This is the prediction-error design — the *only* valid way to test the glucose curve g(BG)
(the cohort ISF-vs-BG regression is artefactually wrong-signed; Phase 8). For each overnight
fasting cycle the loop logged a 2-hour IOB prediction (pred_iob_24) computed with its own ISF
(isf_v1). For a candidate ISF, the predicted 2-h drop scales by isf_cand/isf_loop (drop ∝ ISF
for the insulin actually delivered); we compare the candidate's predicted endpoint to the
actual 2-h BG and score MAE / bias.

Candidates: loop (baseline), log scaler (7d-TDD), power-law (target/BG)^k for swept k, each
with TDD exponent p∈{1, 0.5}. tdd_7day varies ~3.5× within this record, giving limited
within-patient leverage on p.

Caveat: N=1 (the author's own closed-loop). Validates the glucose-curve SHAPE and this
patient's constant; the population TDD exponent comes from the cohort, not here.

Output: results/phase9_powerlaw_backtest.{json,md}
Run: python -m inv008.phase9_powerlaw_backtest
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(os.environ.get("NS_BACKTEST_CSV",
    "/Users/timstreet/Library/CloudStorage/GoogleDrive-tim.street@liveintheirshoes.com/"
    "My Drive/Dynamic ISF data/ns_backtest_overnight.csv"))
OUT = Path(os.environ.get("DYNISF_ROOT", Path.cwd())) / "results"
TARGET, D = 99.0, 82.0          # normal target, insulin divisor (prior work)


def main():
    df = pd.read_csv(DATA)
    df = df[df.pred_iob_24.notna() & df.actual_bg_2h.notna() & df.tdd_7day.notna()
            & df.isf_v1.notna() & (df.isf_v1 > 0) & (df.bg > 0)].copy()
    df["drop_pred"] = df.bg - df.pred_iob_24            # loop's predicted 2h drop (uses isf_v1)
    df = df[df.drop_pred.abs() >= 3]                    # need a real predicted move to scale
    n = len(df)

    bg, tdd, isf_loop, actual = df.bg.values, df.tdd_7day.values, df.isf_v1.values, df.actual_bg_2h.values
    drop_pred = df.drop_pred.values
    logterm = math.log(TARGET / D + 1.0) / np.log(bg / D + 1.0)   # log glucose factor
    plterm = lambda k: (TARGET / bg) ** k                         # power-law glucose factor

    def score(isf_cand):
        pred = bg - drop_pred * (isf_cand / isf_loop)
        err = actual - pred
        return float(np.abs(err).mean()), float(err.mean())

    def fit_C(glucose_factor, p):
        """Best multiplicative scale α for ISF = α · (glucose_factor / tdd^p), so each
        formula is compared at its own optimal level (shape-fair). α0 matches the loop's
        median ISF, then a wide relative search; returns (MAE, bias, α)."""
        base = glucose_factor / (tdd ** p)
        a0 = np.median(isf_loop) / np.median(base)
        best = None
        for a in a0 * np.linspace(0.3, 2.5, 221):
            mae, bias = score(a * base)
            if best is None or mae < best[0]:
                best = (mae, bias, a)
        return best

    results = {}
    # baselines
    results["loop (isf_v1)"] = (*score(isf_loop), None, None)
    mae_log, bias_log, C_log = fit_C(logterm, 1.0)
    results["log scaler /TDD"] = (mae_log, bias_log, C_log, "log, p=1")

    # power-law: sweep k, both TDD exponents
    grid = []
    for p in (1.0, 0.5):
        for k in np.arange(1.0, 6.01, 0.25):
            mae, bias, C = fit_C(plterm(k), p)
            grid.append({"p": p, "k": round(k, 2), "C": round(C, 0),
                         "MAE": round(mae, 2), "Bias": round(bias, 2)})
    g = pd.DataFrame(grid)
    best_p1 = g[g.p == 1.0].loc[g[g.p == 1.0].MAE.idxmin()]
    best_p05 = g[g.p == 0.5].loc[g[g.p == 0.5].MAE.idxmin()]

    loop_mae = results["loop (isf_v1)"][0]
    # per-BG-band MAE for the winner (power-law, p=1, best k) vs log vs loop
    kbest, Cbest = float(best_p1.k), float(best_p1.C)
    isf_pl = Cbest * plterm(kbest) / (tdd ** 1.0)
    isf_log = C_log * logterm / tdd
    df["err_pl"] = actual - (bg - drop_pred * (isf_pl / isf_loop))
    df["err_log"] = actual - (bg - drop_pred * (isf_log / isf_loop))
    df["err_loop"] = actual - (bg - drop_pred)
    df["band"] = pd.cut(df.bg, [70, 90, 105, 120, 150, 200])
    bands = df.groupby("band").apply(lambda x: pd.Series({
        "n": len(x), "loop": np.abs(x.err_loop).mean(),
        "log": np.abs(x.err_log).mean(), "powerlaw": np.abs(x.err_pl).mean()})).round(1)

    summary = {
        "n_cycles": n, "tdd_range": [round(float(tdd.min()), 1), round(float(tdd.max()), 1)],
        "loop_MAE": round(loop_mae, 2),
        "log_MAE": round(mae_log, 2), "log_bias": round(bias_log, 2),
        "best_powerlaw_p1": {"k": kbest, "C": Cbest, "MAE": round(float(best_p1.MAE), 2),
                             "bias": round(float(best_p1.Bias), 2)},
        "best_powerlaw_sqrtTDD": {"k": float(best_p05.k), "C": float(best_p05.C),
                                  "MAE": round(float(best_p05.MAE), 2)},
        "powerlaw_vs_log_pct": round(100 * (mae_log - float(best_p1.MAE)) / mae_log, 1),
        "powerlaw_vs_loop_pct": round(100 * (loop_mae - float(best_p1.MAE)) / loop_mae, 1),
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "phase9_powerlaw_backtest.json").write_text(json.dumps(summary, indent=1))

    s = summary
    md = ["# Phase 9 — consolidated-equation backtest (N=1 overnight prediction data)\n",
          f"{n} overnight fasting cycles; tdd_7day range {s['tdd_range'][0]}–{s['tdd_range'][1]} U/day. "
          "Prediction-error design (scale the loop's 2h drop by isf_cand/isf_loop, compare to actual).\n",
          "## Overall MAE (mg/dL, lower better)\n",
          "| formula | MAE | bias |",
          "|---|---|---|",
          f"| loop (its own ISF) | {s['loop_MAE']} | – |",
          f"| log scaler /TDD | {s['log_MAE']} | {s['log_bias']:+} |",
          f"| **power-law (target/BG)^{s['best_powerlaw_p1']['k']}, /TDD** | "
          f"**{s['best_powerlaw_p1']['MAE']}** | {s['best_powerlaw_p1']['bias']:+} |",
          f"| power-law, /√TDD (best k={s['best_powerlaw_sqrtTDD']['k']}) | "
          f"{s['best_powerlaw_sqrtTDD']['MAE']} | – |",
          "",
          f"- power-law beats log by **{s['powerlaw_vs_log_pct']}%**, and the loop by "
          f"{s['powerlaw_vs_loop_pct']}%.",
          f"- best glucose exponent **k = {s['best_powerlaw_p1']['k']}** (prior work found ≈3.5).",
          f"- √TDD vs 1/TDD: MAE {s['best_powerlaw_sqrtTDD']['MAE']} vs "
          f"{s['best_powerlaw_p1']['MAE']} — within-patient TDD range is small, so this N=1 "
          "data barely distinguishes the TDD exponent (that is the cohort's job; cohort = √TDD).",
          "\n## Per-BG-band MAE\n",
          "| BG band | n | loop | log | power-law |",
          "|---|---|---|---|---|"]
    for b, r in bands.iterrows():
        md.append(f"| {b} | {int(r.n)} | {r.loop} | {r['log']} | {r.powerlaw} |")
    md.append("\n*N=1 (author's own closed-loop). Validates the glucose-curve shape and this "
              "patient's constant; the population TDD exponent (√TDD) comes from the cohort.*")
    (OUT / "phase9_powerlaw_backtest.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
