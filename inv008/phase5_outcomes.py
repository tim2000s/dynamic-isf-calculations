#!/usr/bin/env python3
"""Phase 5: does the derived (sensitivity-anchored) ISF actually beat tuned ISF?

A full closed-loop forward-simulation is too assumption-laden to give a trustworthy
verdict, so this uses a real-outcome, directional test that needs no simulation:

  If the derived ISF (≈ the user's measured empirical sensitivity) is *correct*, then a
  user whose ENTERED ISF is far weaker than it has been under-correcting, and should show
  MORE time-above-range; far stronger → over-correcting → MORE time-below-range.

So we test the sign of the association between the entered-vs-derived ISF gap and realised
glycaemia across users:

  R = entered_ISF / empirical_ISF                  (>1 ⇒ doses weaker than measured sensitivity)
  expect:  log R  ↑  with TAR(+)   and  ↓  with TBR(−)   if empirical/derived ISF is right
           |log R| ↑  ⇒ TIR ↓                            (any mismatch hurts)

A null / wrong-sign result is itself decisive: it means the empirical/derived level is biased
(the 2.4× cohort gap is not real under-dosing) and a sensitivity-anchored ISF would be unsafe —
favouring the profile-anchored (Tier-1) design that preserves the user's working level.

Outcomes are computed from realised CGM. Confound caveat: in a closed loop, basal/SMB/autosens
partially mask a mis-set ISF, attenuating the association.

Output: results/phase5_outcomes.{json,md}
Run: python -m inv008.phase5_outcomes
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
sys.path.insert(0, str(ROOT))
OUT = ROOT / "results"
DSN = "dbname=oref"
N_WORKERS = min(12, mp.cpu_count())
TABLES = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}


def outcomes(args):
    user_id, table = args
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(f"SELECT cgm_mgdl FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL",
                         conn, params=(user_id,))
    finally:
        conn.close()
    g = df["cgm_mgdl"].to_numpy(dtype=float)
    g = g[(g > 10) & (g < 600)]
    if len(g) < 2000:
        return {"user_id": user_id, "skipped": True}
    return {
        "user_id": user_id, "n_cgm": int(len(g)),
        "tir": float(np.mean((g >= 70) & (g <= 180))),
        "tbr": float(np.mean(g < 70)), "tbr54": float(np.mean(g < 54)),
        "tar": float(np.mean(g > 180)), "tar250": float(np.mean(g > 250)),
        "mean_bg": float(np.mean(g)), "cv": float(np.std(g) / np.mean(g)),
    }


def main():
    from canonical_cohort import load_canonical_cohort
    coh = load_canonical_cohort()
    coh = coh[coh["in_cohort"]].copy()
    work = [(r["user_id"], TABLES[r["cohort"]]) for _, r in coh.iterrows()]
    print(f"Phase 5 outcomes: {len(work)} users on {N_WORKERS} workers")
    with mp.Pool(N_WORKERS) as pool:
        oc = pd.DataFrame([r for r in pool.map(outcomes, work, chunksize=2) if not r.get("skipped")])

    emp = pd.DataFrame(json.loads((ROOT / "empirical_isf_v5.json").read_text()))[
        ["user_id", "empirical_isf", "r2"]]
    df = oc.merge(coh[["user_id", "isf", "tdd"]], on="user_id").merge(emp, on="user_id")
    df = df[(df.r2 >= 0.10) & df.empirical_isf.between(5, 500) & df.isf.notna() & (df.isf > 0)]
    df["R"] = df.isf / df.empirical_isf                 # entered / derived(empirical)
    df["logR"] = np.log(df.R)
    df["mismatch"] = df.logR.abs()
    n = len(df)

    def sp(a, b):
        return float(df[[a, b]].corr(method="spearman").iloc[0, 1])

    tests = {
        "logR_vs_TAR (expect +)": sp("logR", "tar"),
        "logR_vs_TBR (expect -)": sp("logR", "tbr"),
        "logR_vs_TIR": sp("logR", "tir"),
        "mismatch_vs_TIR (expect -)": sp("mismatch", "tir"),
        "mismatch_vs_TBR": sp("mismatch", "tbr"),
    }
    # tertiles of R
    df["Rband"] = pd.qcut(df.R, 3, labels=["low (≈ doses to sensitivity)", "mid", "high (doses weak)"])
    band = df.groupby("Rband").agg(
        n=("user_id", "size"), median_R=("R", "median"),
        TIR=("tir", "median"), TAR=("tar", "median"), TBR=("tbr", "median")).reset_index()

    summary = {"n": n, "median_R_entered_over_empirical": round(float(df.R.median()), 2),
               "spearman": {k: round(v, 3) for k, v in tests.items()},
               "tertiles": band.to_dict("records")}
    OUT.mkdir(exist_ok=True)
    (OUT / "phase5_outcomes.json").write_text(json.dumps(summary, indent=1, default=str))

    md = ["# Phase 5 — does the derived (sensitivity-anchored) ISF beat tuned ISF?\n",
          f"{n} users with realised CGM outcomes + valid empirical ISF.",
          f"Median R = entered/empirical = **{summary['median_R_entered_over_empirical']}** "
          "(users dose, on average, well weaker than their measured sensitivity).\n",
          "## Directional test (Spearman across users)\n",
          "If the derived/empirical ISF is *right*, weaker-than-sensitivity dosing (high R) "
          "should mean more time-high and less time-low.\n",
          "| association | ρ | expected sign |",
          "|---|---|---|"]
    exp = {"logR_vs_TAR (expect +)": "+", "logR_vs_TBR (expect -)": "−",
           "mismatch_vs_TIR (expect -)": "−", "logR_vs_TIR": "?", "mismatch_vs_TBR": "?"}
    for k, v in tests.items():
        md.append(f"| {k} | {v:+.2f} | {exp.get(k,'')} |")
    md.append("\n## Outcomes by entered/empirical ratio tertile\n")
    md.append("| R band | n | median R | TIR | TAR | TBR |")
    md.append("|---|---|---|---|---|---|")
    for r in band.to_dict("records"):
        md.append(f"| {r['Rband']} | {r['n']} | {r['median_R']:.1f} | "
                  f"{100*r['TIR']:.0f}% | {100*r['TAR']:.0f}% | {100*r['TBR']:.1f}% |")
    # verdict
    consistent = (tests["logR_vs_TAR (expect +)"] > 0.1 and tests["logR_vs_TBR (expect -)"] < -0.1)
    md.append("\n## Reading\n")
    if consistent:
        md.append("- Associations run in the predicted direction → the gap between entered and "
                  "measured sensitivity tracks real glycaemic cost; moving toward the derived "
                  "ISF would plausibly help. Supports a (cautious) sensitivity-anchored direction.")
    else:
        md.append("- Associations are weak / not in the predicted direction → the large "
                  "entered-vs-empirical gap does **not** translate into the expected glycaemic "
                  "signal. Most likely the empirical level is biased (over-estimates sensitivity), "
                  "so a sensitivity-anchored ISF would dose too strongly. **Favours the "
                  "profile-anchored (Tier-1) design**, which preserves the user's working level "
                  "and only adds the √TDD shape.")
    md.append("\n*Caveat:* closed-loop basal/SMB/autosens partially compensate a mis-set ISF, "
              "attenuating these associations; this is decision-level/observational, single cohort.")
    (OUT / "phase5_outcomes.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
