#!/usr/bin/env python3
"""Rerun the v1 vs v2 comparison with the UPDATED v2 equation.

Updated v2:  sensNormalTarget = 2300 / (ln(target/divisor) · TDD² · 0.02)   [no +1]
             BG floored at divisor+1 (so ln(BG/divisor) > 0)
             → ISF = 115000 / (TDD² · ln(BG_floored/divisor))

vs v1:       ISF = 1800 / (TDD · ln(BG/divisor + 1))

The glucose terms no longer cancel, so the ratio is BG-DEPENDENT:
   ISF_v2u / ISF_v1 = (63.888/TDD) · ln(BG/divisor+1) / ln(BG_floored/divisor)

Reuses the cached per-tick replay (bg, tdd, isf_v1, isf_v2-old) — no DB rerun.
Output: charts/inv008/fig_v2updated.png, results/v1_v2updated_comparison.{json,md}
Run: python -m inv008.compare_v1_v2updated
"""
from __future__ import annotations

import glob
import json
import math
import os
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
    rows = []
    per_user = []
    for f in sorted(glob.glob(str(config.REPLAY_DIR / "*.parquet"))):
        d = pd.read_parquet(f, columns=["bg", "tdd", "isf_v1", "isf_v2"]).dropna()
        d = d[(d.bg > 0) & (d.tdd > 0) & (d.isf_v1 > 0)]
        if len(d) < 500:
            continue
        bg, tdd = d.bg.to_numpy(), d.tdd.to_numpy()
        v1 = d.isf_v1.to_numpy()
        v2_old = d.isf_v2.to_numpy()
        v2u = isf_v2_updated(bg, tdd)
        ok = np.isfinite(v2u) & np.isfinite(v1) & (v1 > 0)
        bg, tdd, v1, v2_old, v2u = bg[ok], tdd[ok], v1[ok], v2_old[ok], v2u[ok]
        rows.append(pd.DataFrame({"bg": bg, "tdd": tdd,
                                  "r_v2u_v1": v2u / v1, "r_v2u_v2old": v2u / v2_old}))
        per_user.append({"user": Path(f).stem, "median_tdd": float(np.median(tdd)),
                         "median_ratio_v2u_v1": float(np.median(v2u / v1)),
                         "median_ratio_v2u_v2old": float(np.median(v2u / v2_old))})
    D = pd.concat(rows, ignore_index=True)
    pu = pd.DataFrame(per_user)

    # ratio by BG band (the new BG dependence)
    D["bgband"] = pd.cut(D.bg, [40, 80, 100, 120, 150, 200, 360])
    by_bg = D.groupby("bgband").agg(n=("r_v2u_v1", "size"),
                                    median_ratio=("r_v2u_v1", "median")).round(2)

    frac_weaker = float((D.r_v2u_v1 > 1).mean())
    summary = {
        "n_users": int(len(pu)), "n_ticks": int(len(D)),
        "median_ratio_v2updated_over_v1": round(float(D.r_v2u_v1.median()), 2),
        "median_ratio_v2updated_over_v2old": round(float(D.r_v2u_v2old.median()), 2),
        "frac_ticks_v2updated_weaker_than_v1": round(frac_weaker, 3),
        "ratio_by_BG_band": {str(b): {"n": int(r.n), "median_ratio": float(r.median_ratio)}
                             for b, r in by_bg.iterrows()},
        "note": "ratio = ISF_v2updated/ISF_v1; >1 means v2updated gives higher ISF (weaker "
                "corrections). Now BG-dependent (was flat 63.9/TDD for old v2).",
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "v1_v2updated_comparison.json").write_text(json.dumps(summary, indent=1))

    # figure: ratio vs BG (the new dependence) + per-user ratio vs TDD
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    # left: ratio vs BG, hexbin-ish via band medians + scatter sample
    s = D.sample(min(40000, len(D)), random_state=0)
    ax[0].scatter(s.bg, s.r_v2u_v1, s=2, alpha=0.05, color="#7f7f7f")
    bm = D.groupby(pd.cut(D.bg, np.arange(60, 320, 10))).r_v2u_v1.median()
    cx = [iv.mid for iv in bm.index]
    ax[0].plot(cx, bm.values, "r-", lw=2, label="median ratio")
    ax[0].axhline(1, color="k", lw=0.8, ls=":")
    ax[0].set_xlabel("glucose (mg/dL)"); ax[0].set_ylabel("ISF_v2updated / ISF_v1")
    ax[0].set_title("Updated v2 vs v1 — now BG-dependent\n(>1: v2updated weaker / more protective)")
    ax[0].set_ylim(0, 8); ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
    # right: per-user median ratio vs TDD
    ax[1].scatter(pu.median_tdd, pu.median_ratio_v2u_v1, s=22, alpha=0.7, color="#1f77b4")
    tt = np.linspace(pu.median_tdd.min(), pu.median_tdd.max(), 100)
    ax[1].plot(tt, 63.888 / tt, "k--", lw=1, label="old v2: 63.9/TDD (BG-indep.)")
    ax[1].axhline(1, color="k", lw=0.8, ls=":")
    ax[1].set_xscale("log"); ax[1].set_xlabel("median TDD (U/day)")
    ax[1].set_ylabel("median ISF_v2updated / ISF_v1 per user")
    ax[1].set_title("Per-user median ratio vs TDD"); ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(CHART / "fig_v2updated.png", dpi=150); plt.close(fig)

    md = ["# v1 vs v2 (UPDATED) comparison\n",
          "Updated v2: `ISF = 115000/(TDD²·ln(BG_floored/divisor))` (anchor "
          "`2300/(ln(target/divisor)·TDD²·0.02)`, no +1; BG floored at divisor+1). "
          f"{summary['n_users']} users, {summary['n_ticks']:,} ticks (cached replay).\n",
          "## Headline\n",
          f"- updated v2 gives **{summary['median_ratio_v2updated_over_v1']}× the ISF of v1** "
          f"(median) — i.e. much weaker corrections; weaker than v1 on "
          f"**{100*summary['frac_ticks_v2updated_weaker_than_v1']:.0f}%** of ticks.",
          f"- updated v2 is **{summary['median_ratio_v2updated_over_v2old']}× the OLD v2** "
          "(median) — the dropped +1 raises ISF substantially.",
          "- the ratio is now **BG-dependent** (old v2 was a flat 63.9/TDD):",
          "\n| BG band | n | median ISF_v2updated/ISF_v1 |", "|---|---|---|"]
    for b, r in by_bg.iterrows():
        md.append(f"| {b} | {int(r.n):,} | {r.median_ratio} |")
    md.append("\n## Reading\n")
    md.append("- The dropped +1 makes ISF blow up as BG approaches the divisor floor → "
              "**very high ISF (near-zero correction) at low BG** — strong hypo protection — "
              "tapering toward ~1.6× v1 at high BG.")
    md.append("- Net effect vs v1: updated v2 is *less* aggressive everywhere, dramatically so "
              "below ~100 mg/dL. Versus the old v2 it is uniformly higher-ISF (the +1 removal "
              "≈3× at target).")
    md.append("\n*Caveat: counterfactual replay on plain-oref cohort; v1 unchanged "
              "(1800/(TDD·ln(BG/div+1))); high cap 210 retained, low floor divisor+1 added.*")
    (OUT / "v1_v2updated_comparison.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
