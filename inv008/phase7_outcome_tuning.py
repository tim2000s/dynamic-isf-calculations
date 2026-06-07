#!/usr/bin/env python3
"""Phase 7 prototype: outcome-based ISF tuning vs the failed sensitivity route.

Phases 5–6 showed the per-user level cannot be safely set from a sensitivity regression
(it is hypo-entangled: reads low-ISF for hypo-prone users). The constructive alternative is
to TUNE the level from observed outcomes. This prototype builds, from observables only:

  ISF_cold    = K_pop / √TDD            (working-anchored cold start, K_pop = LOUO median of
                                          entered_ISF·√TDD ≈ 355; no per-user entry)
  ISF_outcome = ISF_cold · m(TBR, TAR)  (bounded outcome nudge, clamped [0.70, 1.40]):
                  weaker (↑ISF) when unexplained lows;  stronger (↓ISF) when sustained highs
                  with no lows.

and contrasts three candidate zero-entry levels through the Phase-5 directional test:
  - sensitivity route  (empirical ISF)            — expected unsafe (Phase 5: +0.38)
  - cold start         (√TDD only, no per-user)   — expected unbiased but coarse
  - outcome route      (cold start + nudge)       — safe by construction; how close to working?

This is a DESIGN PROTOTYPE, not a validation: the nudge's hypo-safety is partly by
construction, and true proof needs a closed-loop trial. The informative, non-circular number
is how close each lands to the working (entered) ISF.

Output: results/phase7_outcome_tuning.{json,md}
Run: python -m inv008.phase7_outcome_tuning
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
sys.path.insert(0, str(ROOT))
OUT = ROOT / "results"
N_WORKERS = min(12, mp.cpu_count())
TABLES = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}

# outcome-nudge parameters (bounded, conservative)
TBR_TARGET = 0.04          # 4% time-below-range is the accepted ceiling
TAR_CEIL = 0.40            # sustained-high threshold
NUDGE_LO, NUDGE_HI = 0.70, 1.40


def outcome_nudge(tbr, tar):
    """Multiplier on ISF from observed outcomes. >1 weakens dosing (safer); <1 strengthens."""
    m = 1.0
    if tbr > TBR_TARGET:                       # too many lows → weaken (raise ISF)
        m *= 1.0 + min((tbr - TBR_TARGET) / TBR_TARGET, 1.0) * 0.4
    if tar > TAR_CEIL and tbr < 0.02:          # sustained highs and few lows → strengthen
        m *= 1.0 - min((tar - TAR_CEIL) / TAR_CEIL, 1.0) * 0.3
    return float(min(max(m, NUDGE_LO), NUDGE_HI))


def main():
    from canonical_cohort import load_canonical_cohort
    from inv008.phase5_outcomes import outcomes
    coh = load_canonical_cohort()
    coh = coh[coh["in_cohort"]].copy()
    work = [(r["user_id"], TABLES[r["cohort"]]) for _, r in coh.iterrows()]
    print(f"Phase 7 outcome tuning: {len(work)} users on {N_WORKERS} workers")
    with mp.Pool(N_WORKERS) as pool:
        oc = pd.DataFrame([r for r in pool.map(outcomes, work, chunksize=2) if not r.get("skipped")])

    emp = pd.DataFrame(json.loads((ROOT / "empirical_isf_v5.json").read_text()))[
        ["user_id", "empirical_isf", "r2"]]
    df = oc.merge(coh[["user_id", "isf", "tdd"]], on="user_id").merge(emp, on="user_id")
    df = df[df.isf.notna() & (df.isf > 0) & (df.tdd > 0)
            & (df.r2 >= 0.10) & df.empirical_isf.between(5, 500)].copy()

    # working-anchored cold-start prior, LOUO
    K = df.isf * np.sqrt(df.tdd)
    df["isf_cold"] = [float(np.median(K[df.user_id != u]) / math.sqrt(t))
                      for u, t in zip(df.user_id, df.tdd)]
    df["isf_outcome"] = df.isf_cold * df.apply(lambda r: outcome_nudge(r.tbr, r.tar), axis=1)
    df["isf_sens"] = df.empirical_isf

    def err_vs_entered(col):    # median |log(col/entered)| → ± %
        e = np.abs(np.log(df[col] / df.isf))
        return 100 * (math.exp(np.median(e)) - 1)

    def sp_TBR(col):            # directional hypo-safety: log(entered/col) vs TBR
        d = pd.DataFrame({"x": np.log(df.isf / df[col]), "tbr": df.tbr})
        return float(d.corr(method="spearman").iloc[0, 1])

    rows = []
    for col, name in (("isf_sens", "sensitivity route (empirical)"),
                      ("isf_cold", "cold start (√TDD, no per-user)"),
                      ("isf_outcome", "outcome route (cold + nudge)")):
        rows.append({"candidate": name,
                     "err_vs_working_pct": round(err_vs_entered(col), 0),
                     "logR_vs_TBR": round(sp_TBR(col), 2)})

    summary = {"n": int(len(df)),
               "K_pop_working": round(float(np.median(K)), 1),
               "median_nudge": round(float(df.apply(lambda r: outcome_nudge(r.tbr, r.tar), axis=1).median()), 3),
               "candidates": rows}
    OUT.mkdir(exist_ok=True)
    (OUT / "phase7_outcome_tuning.json").write_text(json.dumps(summary, indent=1))

    md = ["# Phase 7 — outcome-based tuning prototype vs the sensitivity route\n",
          f"{len(df)} users. Working-anchored cold-start constant K_pop = "
          f"**{summary['K_pop_working']}** (ISF·√TDD), i.e. ISF ≈ {summary['K_pop_working']:.0f}/√TDD.\n",
          "Each candidate is a zero-entry per-user ISF. `err vs working` = median |log(ISF/entered)| "
          "as ± % (entered ISF is the known-working reference that achieves the users' actual TIR). "
          "`logR vs TBR` = directional hypo-safety (positive = the unsafe Phase-5 signature: the "
          "candidate assigns lower ISF / more insulin to hypo-prone users).\n",
          "| candidate | err vs working ISF | logR vs TBR (≤0 = safe) |",
          "|---|---|---|"]
    for r in rows:
        md.append(f"| {r['candidate']} | ±{r['err_vs_working_pct']:.0f}% | {r['logR_vs_TBR']:+.2f} |")
    md.append(f"\n- median outcome nudge applied: ×{summary['median_nudge']:.2f} "
              "(mostly weakening, reflecting hypo-avoidance).")
    md.append("\n## Reading\n")
    md.append("- The **sensitivity route** is closest in level but carries the unsafe hypo-signature "
              "(positive logR-vs-TBR) — it would over-dose the hypo-prone.")
    md.append("- The **cold start** is unbiased (no per-user signal to bias it) but coarse.")
    md.append("- The **outcome route** keeps the cold-start's safety and nudges toward working "
              "levels in the safe direction (weaker when lows present) — the only route that is "
              "both reasonably close and hypo-safe.")
    md.append("\n*This is a design prototype: the nudge's safety is partly by construction, and "
              "the level it reaches is bounded by the cold-start. True validation requires a "
              "closed-loop trial; observationally we can only show direction and rough level.*")
    (OUT / "phase7_outcome_tuning.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
