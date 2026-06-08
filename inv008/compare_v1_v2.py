#!/usr/bin/env python3
"""Compare the v1 and v2 dynamic-ISF equations across the cohort.

    v1:  ISF(BG) = 1800   / ( TDD  · ln(BG_capped/75 + 1) )
    v2:  ISF(BG) = 115000 / ( TDD² · ln(BG_floored/75) )      BG floored at 76

v1 keeps the +1 in the glucose log; v2 drops it and floors glucose at divisor+1. The two
glucose terms differ, so the v2/v1 ratio depends on glucose as well as TDD: v2 gives a much
higher ISF (weaker correction) at low glucose, easing to a modest margin when high.

Reads the cached per-tick replay (bg, tdd, isf_v1) and evaluates v2 from bg and tdd.
Output: charts/inv008/fig_v1_v2.png, results/v1_v2_comparison.{json,md}
Run: python -m inv008.compare_v1_v2
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config
from inv008.dynisf import isf_v1, isf_v2_updated

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"


def main():
    rows, per_user = [], []
    for f in sorted(glob.glob(str(config.REPLAY_DIR / "*.parquet"))):
        d = pd.read_parquet(f, columns=["bg", "tdd", "isf_v1"]).dropna()
        d = d[(d.bg > 0) & (d.tdd > 0) & (d.isf_v1 > 0)]
        if len(d) < 500:
            continue
        bg, tdd, v1 = d.bg.to_numpy(), d.tdd.to_numpy(), d.isf_v1.to_numpy()
        v2 = isf_v2_updated(bg, tdd)
        ok = np.isfinite(v2) & np.isfinite(v1) & (v1 > 0)
        bg, tdd, v1, v2 = bg[ok], tdd[ok], v1[ok], v2[ok]
        rows.append(pd.DataFrame({"bg": bg, "tdd": tdd, "ratio": v2 / v1}))
        per_user.append({"user": Path(f).stem, "median_tdd": float(np.median(tdd)),
                         "median_ratio": float(np.median(v2 / v1))})
    D = pd.concat(rows, ignore_index=True)
    pu = pd.DataFrame(per_user)

    D["bgband"] = pd.cut(D.bg, [40, 80, 100, 120, 150, 200, 360])
    by_bg = D.groupby("bgband", observed=True).agg(
        n=("ratio", "size"), median_ratio=("ratio", "median")).round(2)
    frac_weaker = float((D.ratio > 1).mean())

    summary = {
        "n_users": int(len(pu)), "n_ticks": int(len(D)),
        "median_ratio_v2_over_v1": round(float(D.ratio.median()), 2),
        "frac_ticks_v2_weaker_than_v1": round(frac_weaker, 3),
        "ratio_by_BG_band": {str(b): {"n": int(r.n), "median_ratio": float(r.median_ratio)}
                             for b, r in by_bg.iterrows()},
        "note": "ratio = ISF_v2/ISF_v1; >1 means v2 gives a higher ISF (weaker correction). "
                "Depends on glucose because the two equations use different glucose terms.",
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    (OUT / "v1_v2_comparison.json").write_text(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    s = D.sample(min(40000, len(D)), random_state=0)
    ax[0].scatter(s.bg, s.ratio, s=2, alpha=0.05, color="#7f7f7f")
    bm = D.groupby(pd.cut(D.bg, np.arange(60, 320, 10)), observed=True).ratio.median()
    ax[0].plot([iv.mid for iv in bm.index], bm.values, "r-", lw=2, label="median ratio")
    ax[0].axhline(1, color="k", lw=0.8, ls=":")
    ax[0].set_xlabel("glucose (mg/dL)"); ax[0].set_ylabel("ISF_v2 / ISF_v1")
    ax[0].set_title("v2 vs v1 by glucose\n(>1: v2 weaker / more protective)")
    ax[0].set_ylim(0, 8); ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    ax[1].scatter(pu.median_tdd, pu.median_ratio, s=22, alpha=0.7, color="#1f77b4")
    ax[1].axhline(1, color="k", lw=0.8, ls=":")
    ax[1].set_xscale("log"); ax[1].set_xlabel("median TDD (U/day)")
    ax[1].set_ylabel("median ISF_v2 / ISF_v1 per user")
    ax[1].set_title("Per-user median ratio vs TDD"); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(CHART / "fig_v1_v2.png", dpi=150); plt.close(fig)

    md = ["# v1 vs v2 comparison\n",
          "v1 `ISF = 1800/(TDD·ln(BG/75+1))`; v2 `ISF = 115000/(TDD²·ln(BG_floored/75))`. "
          f"{summary['n_users']} users, {summary['n_ticks']:,} ticks (cached replay).\n",
          f"- v2 gives a median **{summary['median_ratio_v2_over_v1']}× the ISF of v1** — "
          f"weaker corrections — and is weaker on **{100*summary['frac_ticks_v2_weaker_than_v1']:.0f}%** "
          "of readings.",
          "- The margin depends on glucose (the equations use different glucose terms):",
          "\n| glucose band | n | median ISF_v2 / ISF_v1 |", "|---|---|---|"]
    for b, r in by_bg.iterrows():
        md.append(f"| {b} | {int(r.n):,} | {r.median_ratio} |")
    md += ["\nThe v2 log approaches zero as glucose nears its floor, so v2's ISF climbs steeply "
           "below ~100 mg/dL (near-zero correction when low) and settles to roughly 1.5× v1 when "
           "high. v2 is the gentler equation everywhere, markedly so at low glucose.\n",
           "*Counterfactual replay on the open-source cohort; high cap 210, low floor 76.*"]
    (OUT / "v1_v2_comparison.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
