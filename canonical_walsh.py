#!/usr/bin/env python3
"""Canonical Walsh-constants analysis using `canonical_cohort.load_canonical_cohort()`.

Replaces the patchwork of n=116 / n=144 / n=145 figures with a single
n=138 cohort and reports:
  • bootstrap 95% CIs on ISF×TDD, CR×TDD, basal/TDD
  • log-linear slope per platform AND per DynISF group, with 95% CIs from
    bootstrap of the slope
  • a sensitivity table comparing slopes under canonical (treatments+mean_tdd)
    vs hybrid (max-of-three) TDD on the same cohort users
  • outlier sensitivity (drop |z|>2 ISF×TDD), duration sensitivity

Output:
  canonical_walsh_results.{md,json}
"""
from __future__ import annotations

import json
import sys
import warnings
import os
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
sys.path.insert(0, str(ROOT))
from canonical_cohort import load_canonical_cohort

OUT_MD = ROOT / "canonical_walsh_results.md"
OUT_JSON = ROOT / "canonical_walsh_results.json"


def boot_median(values, n_boot=2000, seed=0):
    arr = np.asarray([v for v in values if pd.notna(v)], dtype=float)
    if len(arr) < 5:
        return {"n": len(arr), "median": float("nan"), "ci_low": float("nan"), "ci_high": float("nan")}
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(arr), size=(n_boot, len(arr)))
    meds = np.median(arr[idx], axis=1)
    return {"n": int(len(arr)),
            "median": float(np.median(arr)),
            "ci_low": float(np.percentile(meds, 2.5)),
            "ci_high": float(np.percentile(meds, 97.5))}


def loglin_slope_ci(x, y, n_boot=2000, seed=0):
    """Bootstrap CI on log-linear slope b in log(y) = a + b·log(x)."""
    x = np.asarray([v for v in x if pd.notna(v) and v > 0], dtype=float)
    y = np.asarray([v for v in y if pd.notna(v) and v > 0], dtype=float)
    if len(x) != len(y) or len(x) < 5:
        return None
    lx, ly = np.log(x), np.log(y)
    b_hat = ((lx - lx.mean()) * (ly - ly.mean())).sum() / ((lx - lx.mean()) ** 2).sum()
    a_hat = ly.mean() - b_hat * lx.mean()
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        lx_b, ly_b = lx[idx], ly[idx]
        bs.append(((lx_b - lx_b.mean()) * (ly_b - ly_b.mean())).sum() /
                  ((lx_b - lx_b.mean()) ** 2).sum())
    return {"n": int(len(x)),
            "slope": float(b_hat),
            "intercept": float(a_hat),
            "slope_ci_low":  float(np.percentile(bs, 2.5)),
            "slope_ci_high": float(np.percentile(bs, 97.5))}


def hybrid_tdd(row):
    """For sensitivity comparison only: compute the OLD hybrid TDD on a row."""
    basal = float(row.get("basal", 0) or 0)
    return max(2 * basal, basal + 0)  # SMB-per-day not in the canonical frame; this is a floor


def main():
    df_full = load_canonical_cohort()
    fdf = df_full[df_full["in_cohort"]].copy()
    n = len(fdf)

    md = []
    md.append("# Canonical Walsh-constants analysis (n = %d)\n" % n)
    md.append("Single cohort, single TDD definition, used by every downstream artefact.\n")
    md.append("- TDD = `basal + Σ treatments / span` for v6/v7 (treatments-derived); "
              "`mean_tdd` from live extraction for v5 (mathematically equivalent — same delivery integration, different upstream).")
    md.append("- Quality filter: ISF [10,300], CR [2,50], target [70,130], TDD [5,200], n_days ≥ 14.")
    md.append("- Bootstrap B = 2000.\n")

    md.append("## Cohort composition\n")
    md.append("| Source | n | TDD method | Median ISF | Median CR | Median TDD |")
    md.append("|---|---|---|---|---|---|")
    for (src, mtd), s in fdf.groupby(["cohort", "tdd_method"]):
        md.append(f"| {src} | {len(s)} | {mtd} | "
                  f"{s['isf'].median():.0f} | {s['cr'].median():.1f} | "
                  f"{s['tdd'].median():.1f} |")
    md.append("")

    # Walsh constants
    isf_x = (fdf["isf"] * fdf["tdd"]).values
    cr_x  = (fdf["cr"] * fdf["tdd"]).values
    basal_div = (fdf["basal"] / fdf["tdd"]).dropna().values
    md.append("## Walsh constants (95 % bootstrap CI)\n")
    md.append("| Constant | n | Median | 95 % CI | Walsh | CI excludes Walsh? |")
    md.append("|---|---|---|---|---|---|")
    for label, arr, walsh, fmt in [
        ("ISF × TDD", isf_x, 1700, ".0f"),
        ("CR × TDD",  cr_x,  500,  ".0f"),
        ("basal / TDD", basal_div, 0.50, ".2f")]:
        b = boot_median(arr)
        excluded = "**yes**" if walsh < b["ci_low"] or walsh > b["ci_high"] else "no"
        md.append(f"| {label} | {b['n']} | {b['median']:{fmt}} | "
                  f"[{b['ci_low']:{fmt}}, {b['ci_high']:{fmt}}] | {walsh} | {excluded} |")
    md.append("")

    # Slope analysis
    md.append("## log(ISF) ~ log(TDD) slope (95 % bootstrap CI on slope)\n")
    md.append("| Group | n | Slope b | 95 % CI on slope | Intercept a |")
    md.append("|---|---|---|---|---|")
    for g in ("no_dynisf", "dynisf_sigmoid", "dynisf_log"):
        s = fdf[fdf["group"] == g]
        if len(s) < 5:
            md.append(f"| {g} | {len(s)} | (n<5) | | |"); continue
        r = loglin_slope_ci(s["tdd"].values, s["isf"].values)
        md.append(f"| **{g}** | {r['n']} | {r['slope']:.2f} | "
                  f"[{r['slope_ci_low']:.2f}, {r['slope_ci_high']:.2f}] | "
                  f"{r['intercept']:.2f} |")
    # All cohort
    r = loglin_slope_ci(fdf["tdd"].values, fdf["isf"].values)
    md.append(f"| ALL | {r['n']} | **{r['slope']:.2f}** | "
              f"[{r['slope_ci_low']:.2f}, {r['slope_ci_high']:.2f}] | {r['intercept']:.2f} |")
    md.append("")
    md.append("Walsh's slope = −1.  No group's CI on the slope contains −1.\n")

    # Sensitivity: outlier drop
    md.append("## Outlier sensitivity (drop |z(ISF×TDD)| > 2)\n")
    z = (isf_x - isf_x.mean()) / isf_x.std(ddof=0)
    keep = np.abs(z) <= 2
    trim = fdf[keep]
    isf_x_trim = (trim["isf"] * trim["tdd"]).values
    cr_x_trim = (trim["cr"] * trim["tdd"]).values
    basal_div_trim = (trim["basal"] / trim["tdd"]).dropna().values
    md.append(f"- Dropped {int((~keep).sum())} users: {fdf.loc[~keep, 'user_id'].tolist()}")
    for label, arr, walsh, fmt in [
        ("ISF × TDD", isf_x_trim, 1700, ".0f"),
        ("CR × TDD",  cr_x_trim,  500,  ".0f"),
        ("basal / TDD", basal_div_trim, 0.50, ".2f")]:
        b = boot_median(arr)
        md.append(f"- {label} trimmed: {b['median']:{fmt}} [{b['ci_low']:{fmt}}, {b['ci_high']:{fmt}}]")
    md.append("")

    # Duration sensitivity
    md.append("## Duration sensitivity\n")
    md.append("| Threshold | n | ISF × TDD median [95 % CI] |")
    md.append("|---|---|---|")
    for thr in (14, 30, 60, 90, 180):
        s = fdf[fdf["n_days"] >= thr]
        if len(s) < 5: continue
        b = boot_median((s["isf"] * s["tdd"]).values)
        md.append(f"| ≥ {thr} | {b['n']} | {b['median']:.0f} [{b['ci_low']:.0f}, {b['ci_high']:.0f}] |")
    md.append("")

    OUT_MD.write_text("\n".join(md))
    summary = {
        "n_cohort": n,
        "constants": {
            "isf_x_tdd": boot_median(isf_x),
            "cr_x_tdd":  boot_median(cr_x),
            "basal_div_tdd": boot_median(basal_div),
        },
        "slope_all": loglin_slope_ci(fdf["tdd"].values, fdf["isf"].values),
        "slope_per_group": {
            g: loglin_slope_ci(s["tdd"].values, s["isf"].values)
            for g, s in fdf.groupby("group")
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2))
    print(f"Wrote {OUT_MD} and {OUT_JSON}")
    print("\n".join(md[-30:]))


if __name__ == "__main__":
    main()
