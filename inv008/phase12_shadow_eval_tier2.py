#!/usr/bin/env python3
"""Phase 12: cohort shadow evaluation of Tier-2 v-next (sensitivity-anchored).

Tier-2 anchors the per-user constant to the person's *measured* sensitivity rather than
their profile ISF:
    K_user = measured_ISF · √(median TDD)        (Tier-1 uses profile_ISF here)
so the equation doses toward observed insulin effect. measured_ISF is the ΔIOB fasting-window
regression sensitivity (empirical_isf_v5.json), quality-gated to r² ≥ 0.10 and 5–500 mg/dL/U.

This sweep computes the Tier-2 ISF across each person's real per-tick history and reports how
much it would change a correction dose versus their static profile and versus Tier-1, both
with and without the §8.2 level clamp (floor the level at profile_ISF / 1.5, applied to the
level term only). The clamp is central here: across the cohort measured sensitivity is far
stronger than tuned profiles, so the unclamped Tier-2 level doses well above profile.

Safety context (not re-derived here): the data-derived study (Phases 5–6) found dosing to
measured sensitivity is hypo-biased — people who run low read as "very sensitive" and would
get the most insulin — which is why Tier-2 is gated behind the clamp and forward validation.
This shadow quantifies the dosing magnitude that gating has to contain.

Patients: 114 users with measured ISF + replay + profile ISF. Parallel over patients.
Output: results/phase12_shadow_eval_tier2.{json,md}, charts/inv008/fig_shadow_eval_tier2.png
Run: python -m inv008.phase12_shadow_eval_tier2
"""
from __future__ import annotations

import glob
import json
import multiprocessing as mp
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config
from inv008.dynisf import g_quartic, k_user_tier1

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
COHORT = config.ROOT / "canonical_cohort.json"
EMPIRICAL = config.ROOT / "empirical_isf_v5.json"

BG_BANDS = [54, 80, 100, 120, 150, 200, 260]
BG_LABELS = ["54-80", "80-100", "100-120", "120-150", "150-200", "200-260"]
TDD_BANDS = [0, 15, 25, 40, 65, 1e9]
TDD_LABELS = ["<15", "15-25", "25-40", "40-65", "65+"]
CLAMP = 1.5            # §8.2: level never more than 1.5x stronger than profile ISF
R2_MIN = 0.10          # fit-quality gate on the measured-sensitivity regression


def to_mgdl(v):
    v = float(v)
    return v * 18.018 if v < 20 else v


def shadow_user(args):
    fpath, profile_isf, measured_isf = args
    d = pd.read_parquet(fpath, columns=["bg", "tdd", "isf_v1"]).dropna()
    d = d[(d.bg > 0) & (d.tdd > 0) & (d.isf_v1 > 0)]
    if len(d) < 500:
        return None
    bg, tdd = d.bg.to_numpy(), d.tdd.to_numpy()
    med_tdd = float(np.median(tdd))
    g = np.asarray(g_quartic(bg))
    f_tdd = np.sqrt(med_tdd / tdd)

    # Tier-2 level = measured_ISF * sqrt(median/tdd); Tier-1 level = profile_ISF * sqrt(...)
    level_t2 = measured_isf * f_tdd
    level_t2_clamped = np.maximum(level_t2, profile_isf / CLAMP)       # §8.2 level floor
    isf_t2 = level_t2 * g
    isf_t2_clamped = level_t2_clamped * g
    isf_t1 = profile_isf * f_tdd * g

    r_uncl = isf_t2 / profile_isf                # sensitivity ratio vs profile (unclamped)
    r_cl = isf_t2_clamped / profile_isf          # ... clamped
    ok = np.isfinite(r_uncl) & (r_uncl > 0)
    bg, r_uncl, r_cl, f_tdd = bg[ok], r_uncl[ok], r_cl[ok], f_tdd[ok]
    isf_t2, isf_t2_clamped, isf_t1 = isf_t2[ok], isf_t2_clamped[ok], isf_t1[ok]

    def band_median(values, x):
        idx = np.digitize(x, BG_BANDS) - 1
        return {i: (float(np.median(values[idx == i])) if (idx == i).sum() >= 20 else None)
                for i in range(len(BG_LABELS))}

    return {
        "user": Path(fpath).stem, "n": int(len(r_uncl)), "median_tdd": round(med_tdd, 1),
        "profile_isf": round(profile_isf, 1), "measured_isf": round(measured_isf, 1),
        "measured_over_profile": float(measured_isf / profile_isf),
        # dose change vs profile (units ∝ 1/ISF): unclamped vs clamped
        "med_abs_dose_change_uncl": float(np.median(np.abs(1.0 / r_uncl - 1.0))),
        "med_abs_dose_change_cl": float(np.median(np.abs(1.0 / r_cl - 1.0))),
        # fraction of ticks the level clamp binds (level would be >1.5x stronger than profile)
        "frac_level_clamp_binds": float((f_tdd * measured_isf < profile_isf / CLAMP).mean()),
        # fraction of ticks dosing >1.5x stronger than profile, incl. glucose (unclamped/clamped)
        "frac_stronger_15x_uncl": float((r_uncl < 1.0 / CLAMP).mean()),
        "frac_stronger_15x_cl": float((r_cl < 1.0 / CLAMP).mean()),
        # vs Tier-1 (median over ticks)
        "med_t2_over_t1_uncl": float(np.median(isf_t2 / isf_t1)),
        "med_t2_over_t1_cl": float(np.median(isf_t2_clamped / isf_t1)),
        "r_by_bg_uncl": band_median(r_uncl, bg),
        "r_by_bg_cl": band_median(r_cl, bg),
    }


def main():
    cohort = {r["user_id"]: r for r in json.load(open(COHORT))}
    measured = {e["user_id"]: e for e in json.load(open(EMPIRICAL))}
    jobs = []
    for f in sorted(glob.glob(str(config.REPLAY_DIR / "*.parquet"))):
        u = Path(f).stem
        e = measured.get(u)
        if u in cohort and e and e.get("r2", 0) >= R2_MIN and 5 <= e.get("empirical_isf", 0) <= 500:
            jobs.append((f, to_mgdl(cohort[u]["isf"]), float(e["empirical_isf"])))
    n_workers = min(config.DEFAULT_WORKERS, mp.cpu_count())
    print(f"Phase 12 (Tier-2): {len(jobs)} users with measured ISF on {n_workers} workers")
    with mp.Pool(n_workers) as pool:
        res = [r for r in pool.map(shadow_user, jobs) if r]
    n = len(res)

    def xmed(key):
        return round(float(np.median([r[key] for r in res if r.get(key) is not None])), 3)

    def band_xmed(key, i):
        vals = [r[key][i] for r in res if r[key].get(i) is not None]
        return round(float(np.median(vals)), 3) if vals else None

    summary = {
        "n_patients": n, "total_ticks": int(sum(r["n"] for r in res)),
        "tier": 2, "anchor": "measured_ISF * sqrt(median TDD)", "r2_gate": R2_MIN,
        "median_measured_over_profile_isf": xmed("measured_over_profile"),
        "implied_dose_vs_profile_at_level": round(1.0 / xmed("measured_over_profile"), 2),
        "median_abs_dose_change_pct_unclamped": round(xmed("med_abs_dose_change_uncl") * 100, 1),
        "median_abs_dose_change_pct_clamped": round(xmed("med_abs_dose_change_cl") * 100, 1),
        "clamp_threshold_x": CLAMP,
        "median_frac_ticks_level_clamp_binds": xmed("frac_level_clamp_binds"),
        "median_frac_ticks_stronger_1p5x_unclamped": xmed("frac_stronger_15x_uncl"),
        "median_frac_ticks_stronger_1p5x_clamped": xmed("frac_stronger_15x_cl"),
        "median_t2_over_t1_unclamped": xmed("med_t2_over_t1_uncl"),
        "median_t2_over_t1_clamped": xmed("med_t2_over_t1_cl"),
        "sensitivity_ratio_vs_profile_by_bg_unclamped":
            {BG_LABELS[i]: band_xmed("r_by_bg_uncl", i) for i in range(len(BG_LABELS))},
        "sensitivity_ratio_vs_profile_by_bg_clamped":
            {BG_LABELS[i]: band_xmed("r_by_bg_cl", i) for i in range(len(BG_LABELS))},
    }

    # by per-user median TDD
    tdds = np.array([r["median_tdd"] for r in res])
    tband = np.digitize(tdds, TDD_BANDS) - 1
    by_tdd = []
    for i, lab in enumerate(TDD_LABELS):
        grp = [r for r, b in zip(res, tband) if b == i]
        if not grp:
            continue
        by_tdd.append({"tdd_band": lab, "n_users": len(grp),
                       "median_measured_over_profile": round(float(np.median(
                           [r["measured_over_profile"] for r in grp])), 3),
                       "median_frac_clamp_binds": round(float(np.median(
                           [r["frac_level_clamp_binds"] for r in grp])), 3)})
    summary["by_tdd_band"] = by_tdd

    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    (OUT / "phase12_shadow_eval_tier2.json").write_text(json.dumps(summary, indent=1))

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(BG_LABELS))
    ru = [summary["sensitivity_ratio_vs_profile_by_bg_unclamped"][l] for l in BG_LABELS]
    rc = [summary["sensitivity_ratio_vs_profile_by_bg_clamped"][l] for l in BG_LABELS]
    ax[0].plot(x, ru, "o-", color="#d62728", lw=2, label="Tier-2, no clamp")
    ax[0].plot(x, rc, "s--", color="#1f77b4", lw=2, label="Tier-2, level clamp")
    ax[0].axhline(1, color="k", lw=0.8, ls=":")
    ax[0].axhline(1 / CLAMP, color="green", lw=0.8, ls="--", alpha=0.6,
                  label=f"profile/{CLAMP} (clamp floor)")
    ax[0].set_xticks(x); ax[0].set_xticklabels(BG_LABELS, rotation=30)
    ax[0].set_xlabel("glucose band (mg/dL)"); ax[0].set_ylabel("ISF Tier-2 / profile ISF")
    ax[0].set_title("Tier-2 sensitivity ratio vs profile\n(<1 = more aggressive / more insulin)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    mop = np.array([r["measured_over_profile"] for r in res])
    ax[1].scatter(tdds, mop, s=24, alpha=0.7, color="#d62728")
    ax[1].axhline(1, color="k", lw=0.8, ls=":")
    ax[1].axhline(1 / CLAMP, color="green", lw=0.8, ls="--", alpha=0.6, label=f"1/{CLAMP} clamp floor")
    ax[1].set_xscale("log"); ax[1].set_xlabel("per-user median TDD (U/day)")
    ax[1].set_ylabel("measured ISF / profile ISF")
    ax[1].set_title("Measured sensitivity vs profile (the Tier-2 level)\nbelow the green line, the clamp binds")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(CHART / "fig_shadow_eval_tier2.png", dpi=150); plt.close(fig)

    md = ["# Phase 12 — cohort shadow evaluation of Tier-2 v-next (sensitivity-anchored)\n",
          f"{n} people with a usable measured-sensitivity fit (r² ≥ {R2_MIN}), "
          f"{summary['total_ticks']:,} per-tick readings. Tier-2 anchors "
          "K_user = measured_ISF · √(median TDD); the glucose curve g(BG) is unchanged.\n",
          "## Headline\n",
          f"- Measured sensitivity is a median **{summary['median_measured_over_profile_isf']}× "
          f"the profile ISF**, so at the level Tier-2 doses about "
          f"**{summary['implied_dose_vs_profile_at_level']}× the insulin** of the person's "
          "current setting — much more aggressive than Tier-1, which preserves the level.",
          f"- **Without the clamp**, a correction dose changes by a median "
          f"**{summary['median_abs_dose_change_pct_unclamped']}%** vs profile, and "
          f"**{100*summary['median_frac_ticks_stronger_1p5x_unclamped']:.0f}%** of readings "
          "dose more than 1.5× stronger than profile.",
          f"- **With the §8.2 level clamp** (floor the level at profile/{CLAMP}), the level clamp "
          f"binds for a median **{100*summary['median_frac_ticks_level_clamp_binds']:.0f}%** of "
          "readings — i.e. for most people Tier-2 is pulled back to the same ceiling Tier-1 "
          f"would allow; dose change vs profile falls to a median "
          f"**{summary['median_abs_dose_change_pct_clamped']}%**.",
          f"- vs Tier-1: median ISF ratio **{summary['median_t2_over_t1_unclamped']}** unclamped, "
          f"**{summary['median_t2_over_t1_clamped']}** clamped.",
          "\n## Sensitivity ratio vs profile, by glucose band\n",
          "| BG band | Tier-2 no clamp | Tier-2 clamped |", "|---|---|---|"]
    for l in BG_LABELS:
        md.append(f"| {l} | {summary['sensitivity_ratio_vs_profile_by_bg_unclamped'][l]} | "
                  f"{summary['sensitivity_ratio_vs_profile_by_bg_clamped'][l]} |")
    md += ["\n*<1 = more aggressive (more insulin than the current profile). Unclamped Tier-2 is "
           "well below 1 across the range; the clamp lifts the level back toward profile, leaving "
           "the glucose-shape behaviour intact.*\n",
           "## By per-user median TDD\n",
           "| TDD band | users | median measured/profile | median frac clamp binds |",
           "|---|---|---|---|"]
    for row in by_tdd:
        md.append(f"| {row['tdd_band']} | {row['n_users']} | "
                  f"{row['median_measured_over_profile']} | {row['median_frac_clamp_binds']} |")
    md += ["\n![Tier-2 shadow eval](charts/inv008/fig_shadow_eval_tier2.png)\n",
           "**Reading.** Tier-2 is materially more aggressive than the person's current profile — "
           f"about {summary['implied_dose_vs_profile_at_level']}× the correction insulin at the "
           "level, before glucose scaling. The §8.2 level clamp is doing real work here: it binds "
           "for the majority of readings and converts Tier-2 into \"no more than 1.5× stronger "
           "than profile\", which is also the Tier-1 ceiling. So clamped Tier-2 and Tier-1 differ "
           "mainly where measured sensitivity is *weaker* than profile (a minority of users).\n",
           "This is the dosing-magnitude side of the Tier-2 question only. The data-derived study "
           "(Phases 5–6) showed the measured-sensitivity anchor is hypo-biased — it reads most "
           "sensitive for the people who already run low — so even the clamped form needs forward, "
           "outcome-based validation before it doses. It is not a deployable default.\n",
           "*Caveat: counterfactual decision-level replay, not closed-loop outcomes; measured ISF "
           "carries carb/endogenous-glucose confounds and per-user CI; single cohort; median-TDD "
           "anchor stands in for the weekly recalibration.*"]
    (OUT / "phase12_shadow_eval_tier2.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
