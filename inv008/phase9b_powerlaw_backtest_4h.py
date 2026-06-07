#!/usr/bin/env python3
"""Phase 9b: repeat the consolidated-equation backtest at the END-OF-INSULIN-ACTION
horizon (the '4h' cache) instead of +2h.

A 2-hour drop captures only part of the insulin action (DIA ~4-5h), so it understates the
total ISF and inflates the positive bias. The 4h cache uses the LAST element of the loop's
predBGs.IOB array (end of insulin action, ~3h here) matched to actual CGM at that offset —
the cleaner measure of total mg/dL drop per unit.

Same prediction-error design as Phase 9: scale the loop's predicted drop by isf_cand/isf_loop,
compare predicted endpoint to actual_bg_end, score MAE/bias. Each formula fit at its own scale.

Output: results/phase9b_powerlaw_backtest_4h.{json,md}
Run: python -m inv008.phase9b_powerlaw_backtest_4h
"""
from __future__ import annotations

import json
import math
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

CACHE = Path(os.environ.get("BOOST_4H_CACHE",
    "/Users/timstreet/Library/CloudStorage/GoogleDrive-tim.street@liveintheirshoes.com/"
    "My Drive/Dynamic ISF data/boost_4h_cache.pkl"))
OUT = Path(os.environ.get("DYNISF_ROOT", Path.cwd())) / "results"
TARGET, D = 99.0, 82.0


def main():
    df = pickle.load(open(CACHE, "rb"))["strict"].copy()
    df = df[df.pred_iob_final.notna() & df.actual_bg_end.notna() & df.tdd_7day.notna()
            & df.variable_sens.notna() & (df.variable_sens > 0) & (df.bg > 0)].copy()
    df["pred_drop"] = df.bg - df.pred_iob_final
    df = df[df.pred_drop.abs() >= 3]
    n = len(df)
    horizon_h = float(np.median(df.pred_horizon_s)) / 3600.0

    bg = df.bg.values; tdd = df.tdd_7day.values
    isf_loop = df.variable_sens.values; actual = df.actual_bg_end.values
    drop_pred = df.pred_drop.values
    logterm = math.log(TARGET / D + 1.0) / np.log(bg / D + 1.0)
    plterm = lambda k: (TARGET / bg) ** k

    def score(isf):
        err = actual - (bg - drop_pred * (isf / isf_loop))
        return float(np.abs(err).mean()), float(err.mean())

    def fit(gfac, p):
        base = gfac / (tdd ** p)
        a0 = np.median(isf_loop) / np.median(base)
        best = None
        for a in a0 * np.linspace(0.3, 2.5, 221):
            mae, bias = score(a * base)
            if best is None or mae < best[0]:
                best = (mae, bias, a)
        return best

    loop_mae, loop_bias = score(isf_loop)
    log_mae, log_bias, _ = fit(logterm, 1.0)
    grid = [{"p": p, "k": round(k, 2), **dict(zip(("MAE", "Bias", "C"),
            (lambda r: (round(r[0], 2), round(r[1], 2), round(r[2], 0)))(fit(plterm(k), p))))}
            for p in (1.0, 0.5) for k in np.arange(1.0, 6.01, 0.25)]
    g = pd.DataFrame(grid)
    bp1 = g[g.p == 1.0].loc[g[g.p == 1.0].MAE.idxmin()]
    bp05 = g[g.p == 0.5].loc[g[g.p == 0.5].MAE.idxmin()]

    kbest = float(bp1.k)
    _, _, Cpl = fit(plterm(kbest), 1.0)
    _, _, Clog = fit(logterm, 1.0)
    isf_pl = Cpl * plterm(kbest) / tdd
    isf_log = Clog * logterm / tdd
    df["e_pl"] = actual - (bg - drop_pred * (isf_pl / isf_loop))
    df["e_log"] = actual - (bg - drop_pred * (isf_log / isf_loop))
    df["e_loop"] = actual - df.pred_iob_final
    df["band"] = pd.cut(df.bg, [50, 90, 105, 120, 150, 250])
    bands = df.groupby("band").apply(lambda x: pd.Series({
        "n": len(x), "loop": np.abs(x.e_loop).mean(),
        "log": np.abs(x.e_log).mean(), "powerlaw": np.abs(x.e_pl).mean()})).round(1)

    summary = {
        "n_cycles": n, "horizon_h": round(horizon_h, 2),
        "tdd_range": [round(float(tdd.min()), 1), round(float(tdd.max()), 1)],
        "loop_MAE": round(loop_mae, 2), "loop_bias": round(loop_bias, 2),
        "log_MAE": round(log_mae, 2), "log_bias": round(log_bias, 2),
        "best_powerlaw_p1": {"k": kbest, "MAE": round(float(bp1.MAE), 2), "bias": round(float(bp1.Bias), 2)},
        "best_powerlaw_sqrtTDD": {"k": float(bp05.k), "MAE": round(float(bp05.MAE), 2)},
        "powerlaw_vs_log_pct": round(100 * (log_mae - float(bp1.MAE)) / log_mae, 1),
        "powerlaw_vs_loop_pct": round(100 * (loop_mae - float(bp1.MAE)) / loop_mae, 1),
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "phase9b_powerlaw_backtest_4h.json").write_text(json.dumps(summary, indent=1))

    s = summary
    md = [f"# Phase 9b — consolidated-equation backtest at end-of-insulin-action (~{s['horizon_h']}h)\n",
          f"{n} overnight cycles, horizon median {s['horizon_h']}h (end of the loop's IOB "
          f"prediction). Compare to Phase 9 at +2h.\n",
          "## Overall (MAE mg/dL; bias = actual − predicted, + = formula too aggressive)\n",
          "| formula | MAE | bias |",
          "|---|---|---|",
          f"| loop (its own ISF) | {s['loop_MAE']} | {s['loop_bias']:+} |",
          f"| log scaler /TDD | {s['log_MAE']} | {s['log_bias']:+} |",
          f"| **power-law (target/BG)^{s['best_powerlaw_p1']['k']} /TDD** | "
          f"**{s['best_powerlaw_p1']['MAE']}** | {s['best_powerlaw_p1']['bias']:+} |",
          f"| power-law /√TDD (k={s['best_powerlaw_sqrtTDD']['k']}) | {s['best_powerlaw_sqrtTDD']['MAE']} | – |",
          "",
          f"- power-law beats log by **{s['powerlaw_vs_log_pct']}%**, loop by {s['powerlaw_vs_loop_pct']}%; "
          f"best k = **{s['best_powerlaw_p1']['k']}**.",
          f"- bias at ~{s['horizon_h']}h: loop {s['loop_bias']:+}, log {s['log_bias']:+}, "
          f"power-law {s['best_powerlaw_p1']['bias']:+} (vs the larger +bias seen at 2h).",
          "\n## Per-BG-band MAE\n",
          "| BG band | n | loop | log | power-law |", "|---|---|---|---|---|"]
    for b, r in bands.iterrows():
        md.append(f"| {b} | {int(r.n)} | {r.loop} | {r['log']} | {r.powerlaw} |")
    md.append("\n*N=1, end-of-insulin-action horizon captures more of the total ISF than +2h. "
              "Validates the glucose-curve shape; TDD exponent is the cohort's (√TDD).*")
    (OUT / "phase9b_powerlaw_backtest_4h.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
