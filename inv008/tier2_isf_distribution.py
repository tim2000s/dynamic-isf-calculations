#!/usr/bin/env python3
"""Distribution of the Tier-2 v-next ISF each person's equation produces over observed data.

For every reading in a person's replay we compute the Tier-2 ISF:
    ISF_t2(BG, TDD) = measured_ISF · √(median TDD / TDD) · g(BG)
where measured_ISF is the person's ΔIOB fasting-window sensitivity (the Tier-2 anchor) and
g(BG) is the Diabeloop quartic. This is the raw calculation (no §8.2 clamp); the spread within
a person comes from their glucose range (g(BG)) and their TDD swings (√(median/TDD)).

We report, per person, the distribution of those ISF values (percentiles), and then aggregate
to the population: the distribution of per-person medians, and the pooled distribution of all
per-reading Tier-2 ISF values.

Patients: 114 users with a usable measured-sensitivity fit (r² ≥ 0.10) + replay + profile.
Output: results/tier2_isf_distribution.{json,md},
        charts/inv008/fig_tier2_per_person.png, charts/inv008/fig_tier2_population.png
Run: python -m inv008.tier2_isf_distribution
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
from inv008.dynisf import g_quartic

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
COHORT = config.ROOT / "canonical_cohort.json"
EMPIRICAL = config.ROOT / "empirical_isf_v5.json"
R2_MIN = 0.10
CLAMP = 1.5
PCTS = [5, 10, 25, 50, 75, 90, 95]
BINS = np.logspace(np.log10(5), np.log10(1000), 61)   # ISF mg/dL per U, log-spaced


def to_mgdl(v):
    v = float(v)
    return v * 18.018 if v < 20 else v


def person_dist(args):
    fpath, profile_isf, measured_isf = args
    d = pd.read_parquet(fpath, columns=["bg", "tdd"]).dropna()
    d = d[(d.bg > 0) & (d.tdd > 0)]
    if len(d) < 500:
        return None
    bg, tdd = d.bg.to_numpy(), d.tdd.to_numpy()
    med_tdd = float(np.median(tdd))
    isf_t2 = measured_isf * np.sqrt(med_tdd / tdd) * np.asarray(g_quartic(bg))
    isf_t2 = isf_t2[np.isfinite(isf_t2) & (isf_t2 > 0)]
    if len(isf_t2) < 500:
        return None
    pct = {str(p): float(np.percentile(isf_t2, p)) for p in PCTS}
    counts, _ = np.histogram(isf_t2, bins=BINS)
    return {
        "user": Path(fpath).stem, "n": int(len(isf_t2)),
        "median_tdd": round(med_tdd, 1),
        "profile_isf": round(profile_isf, 1), "measured_isf": round(measured_isf, 1),
        "isf_t2_pct": {k: round(v, 1) for k, v in pct.items()},
        "isf_t2_mean": round(float(isf_t2.mean()), 1),
        "hist": counts.tolist(),
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
    print(f"Tier-2 ISF distribution: {len(jobs)} users on {n_workers} workers")
    with mp.Pool(n_workers) as pool:
        res = [r for r in pool.map(person_dist, jobs) if r]
    res.sort(key=lambda r: r["isf_t2_pct"]["50"])     # sort people by median Tier-2 ISF
    n = len(res)

    med = np.array([r["isf_t2_pct"]["50"] for r in res])
    pooled_hist = np.sum([np.asarray(r["hist"]) for r in res], axis=0)
    centres = np.sqrt(BINS[:-1] * BINS[1:])

    def pooled_pct(p):  # population percentile over the pooled per-reading histogram
        c = np.cumsum(pooled_hist) / pooled_hist.sum()
        return float(centres[np.searchsorted(c, p / 100.0)])

    summary = {
        "n_patients": n, "total_readings": int(sum(r["n"] for r in res)),
        "anchor": "measured_ISF * sqrt(median TDD / TDD) * g(BG)  (unclamped)",
        "per_person_median_isf": {
            "min": round(float(med.min()), 1), "p25": round(float(np.percentile(med, 25)), 1),
            "median": round(float(np.median(med)), 1),
            "p75": round(float(np.percentile(med, 75)), 1), "max": round(float(med.max()), 1)},
        "population_pooled_isf_pct": {str(p): round(pooled_pct(p), 1)
                                      for p in (5, 25, 50, 75, 95)},
        "median_within_person_p10_p90_span": round(float(np.median(
            [r["isf_t2_pct"]["90"] / r["isf_t2_pct"]["10"] for r in res])), 2),
        "per_person": [
            {"user": r["user"], "median_tdd": r["median_tdd"], "n": r["n"],
             "profile_isf": r["profile_isf"], "measured_isf": r["measured_isf"],
             "isf_t2_pct": r["isf_t2_pct"]} for r in res],
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    (OUT / "tier2_isf_distribution.json").write_text(json.dumps(summary, indent=1))

    # ---- figure 1: per-person Tier-2 ISF distribution (sorted by median) ----
    fig, ax = plt.subplots(figsize=(15, 6))
    x = np.arange(n)
    p10 = np.array([r["isf_t2_pct"]["10"] for r in res])
    p25 = np.array([r["isf_t2_pct"]["25"] for r in res])
    p50 = np.array([r["isf_t2_pct"]["50"] for r in res])
    p75 = np.array([r["isf_t2_pct"]["75"] for r in res])
    p90 = np.array([r["isf_t2_pct"]["90"] for r in res])
    meas = np.array([r["measured_isf"] for r in res])
    prof = np.array([r["profile_isf"] for r in res])
    ax.vlines(x, p10, p90, color="#c6c6e8", lw=1.4)                  # p10-p90 whisker
    ax.vlines(x, p25, p75, color="#6a6ad6", lw=3.0)                  # p25-p75 box
    ax.plot(x, p50, ".", color="#1a1a8c", ms=4, label="Tier-2 ISF median (p25–p75, p10–p90)")
    ax.plot(x, meas, "_", color="#d62728", ms=7, mew=1.6, label="measured ISF (anchor)")
    ax.plot(x, prof, "x", color="#2ca02c", ms=4, mew=1.0, label="profile ISF")
    ax.set_yscale("log")
    ax.set_xlabel("person (114, sorted by median Tier-2 ISF)")
    ax.set_ylabel("ISF (mg/dL per U, log scale)")
    ax.set_title("Tier-2 v-next ISF computed over each person's observed data\n"
                 "spread within a person = glucose range × TDD swings; anchor = measured sensitivity")
    ax.legend(fontsize=9, loc="upper left"); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(CHART / "fig_tier2_per_person.png", dpi=150); plt.close(fig)

    # ---- figure 2: population view ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].hist(med, bins=np.logspace(np.log10(med.min()), np.log10(med.max()), 30),
               color="#6a6ad6", alpha=0.85)
    ax[0].axvline(np.median(med), color="k", ls="--", lw=1.5,
                  label=f"median {np.median(med):.0f}")
    ax[0].set_xscale("log"); ax[0].set_xlabel("per-person median Tier-2 ISF (mg/dL per U)")
    ax[0].set_ylabel("people"); ax[0].set_title(f"Per-person median Tier-2 ISF ({n} people)")
    ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3, which="both")
    ax[1].bar(centres, pooled_hist / pooled_hist.sum(), width=np.diff(BINS),
              align="center", color="#1a1a8c", alpha=0.8)
    ax[1].axvline(pooled_pct(50), color="k", ls="--", lw=1.5,
                  label=f"median {pooled_pct(50):.0f}")
    ax[1].set_xscale("log"); ax[1].set_xlabel("Tier-2 ISF per reading (mg/dL per U)")
    ax[1].set_ylabel("fraction of readings")
    ax[1].set_title(f"Pooled distribution of every Tier-2 calculation\n({summary['total_readings']:,} readings)")
    ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(CHART / "fig_tier2_population.png", dpi=150); plt.close(fig)

    md = ["# Tier-2 v-next ISF — distribution over observed data\n",
          f"{n} people, {summary['total_readings']:,} readings. For each reading, "
          "ISF_Tier-2 = measured_ISF · √(median TDD / TDD) · g(BG) (raw, unclamped).\n",
          "## Population\n",
          f"- Per-person median Tier-2 ISF spans {summary['per_person_median_isf']['min']}–"
          f"{summary['per_person_median_isf']['max']} mg/dL per U "
          f"(IQR {summary['per_person_median_isf']['p25']}–{summary['per_person_median_isf']['p75']}, "
          f"median {summary['per_person_median_isf']['median']}).",
          f"- Pooled over every reading: median {summary['population_pooled_isf_pct']['50']}, "
          f"p5–p95 {summary['population_pooled_isf_pct']['5']}–{summary['population_pooled_isf_pct']['95']}.",
          f"- Within a person the Tier-2 ISF spreads a median "
          f"{summary['median_within_person_p10_p90_span']}× from p10 to p90 (glucose range + TDD swings).",
          "\n![Per-person Tier-2 ISF distribution](charts/inv008/fig_tier2_per_person.png)\n",
          "![Population view](charts/inv008/fig_tier2_population.png)\n",
          "## Per-person (sorted by median Tier-2 ISF)\n",
          "| person | TDD | n | profile ISF | measured ISF | Tier-2 p10 | p50 | p90 |",
          "|---|---|---|---|---|---|---|---|"]
    for r in res:
        p = r["isf_t2_pct"]
        md.append(f"| {r['user']} | {r['median_tdd']:.0f} | {r['n']} | {r['profile_isf']:.0f} | "
                  f"{r['measured_isf']:.0f} | {p['10']:.0f} | {p['50']:.0f} | {p['90']:.0f} |")
    md.append("\n*Raw Tier-2 calculation, no §8.2 clamp. The median sits near the person's "
              "measured ISF at typical TDD; readings above target pull ISF down (g(BG)<1, firmer "
              "corrections) and readings below target push it up. Measured ISF is well below "
              "profile ISF for most people, which is why Tier-2 doses more than the current "
              "profile — see the Tier-2 shadow evaluation.*")
    (OUT / "tier2_isf_distribution.md").write_text("\n".join(md))
    print("\n".join(md[:14]))


if __name__ == "__main__":
    main()
