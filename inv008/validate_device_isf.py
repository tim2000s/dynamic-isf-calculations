#!/usr/bin/env python3
"""Validate the Python V1 ISF replay against device-logged ISF (oref_v5.sug_ISF).

For Trio (v5) DynISF users we have the ISF the device itself computed and logged at
each cycle. This script compares our replayed V1 ISF against that ground truth, after
correcting a data-quality issue: four v5 users switched the *units* of sug_ISF
mid-history (mmol/L-per-U vs mg/dL-per-U). Values < 20 are treated as mmol and scaled
by 18.018.

We expect, per user:
  * a positive log-log correlation (the formula tracks the device curve shape), and
  * a roughly constant multiplicative offset (device applies adjustmentFactor, insulin
    divisor and autosens, which the replay deliberately does not model) — i.e. a tight
    ratio IQR, not a ratio of exactly 1.

Output: results/device_isf_validation.{json,md}
"""
from __future__ import annotations

import json
import warnings
from datetime import datetime
import os
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
REPLAY = ROOT / "inv008_cache" / "replay"
OUT = ROOT / "results"
MMOL_THRESHOLD = 20.0       # sug_ISF below this is logged in mmol/L per U
MMOL_TO_MGDL = 18.018


def correct_units(dev_isf: pd.Series) -> pd.Series:
    """Per-tick unit correction for the four mixed-unit v5 users."""
    return np.where(dev_isf < MMOL_THRESHOLD, dev_isf * MMOL_TO_MGDL, dev_isf)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    coh = {r["user_id"]: r for r in
           json.loads((ROOT / "canonical_cohort.json").read_text())}
    conn = psycopg2.connect(host="localhost", dbname="oref")

    rows = []
    for uid, rec in coh.items():
        if rec.get("src_table") != "oref_v5" or not rec.get("dynisf_user"):
            continue
        pq = REPLAY / f"{uid}.parquet"
        if not pq.exists():
            continue
        rep = pd.read_parquet(pq)
        db = pd.read_sql('SELECT ts_relative_sec, "sug_ISF" AS dev_isf FROM oref_v5 '
                         'WHERE user_id=%s AND "sug_ISF" IS NOT NULL',
                         conn, params=(uid,))
        m = rep.merge(db, on="ts_relative_sec").dropna(subset=["isf_v1", "dev_isf"])
        m = m[(m.dev_isf > 1) & (m.dev_isf < 500) & np.isfinite(m.isf_v1)]
        if len(m) < 1000:
            continue
        m["dev_mgdl"] = correct_units(m["dev_isf"])
        mixed = bool((m.dev_isf < MMOL_THRESHOLD).mean() > 0.05
                     and (m.dev_isf >= MMOL_THRESHOLD).mean() > 0.05)
        corr = float(np.corrcoef(np.log(m.isf_v1), np.log(m.dev_mgdl))[0, 1])
        ratio = m.isf_v1 / m.dev_mgdl
        rows.append({
            "user": uid, "formula": rec.get("formula"), "n_ticks": int(len(m)),
            "mixed_units": mixed,
            "log_corr": round(corr, 3),
            "median_ratio": round(float(ratio.median()), 3),
            "ratio_iqr": round(float(ratio.quantile(.75) - ratio.quantile(.25)), 3),
        })
    conn.close()

    df = pd.DataFrame(rows).sort_values(["formula", "user"])
    summary = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_users": int(len(df)),
        "mixed_unit_users": df.loc[df.mixed_units, "user"].tolist(),
        "log_dynisf": {
            "n": int((df.formula == "log").sum()),
            "median_log_corr": round(float(df[df.formula == "log"].log_corr.median()), 3),
            "median_ratio_iqr": round(float(df[df.formula == "log"].ratio_iqr.median()), 3),
        },
        "per_user": df.to_dict("records"),
    }
    (OUT / "device_isf_validation.json").write_text(json.dumps(summary, indent=1))

    lines = ["# Device-ISF validation (oref_v5 sug_ISF vs replayed V1)",
             f"\n{summary['generated']} · {summary['n_users']} Trio DynISF users",
             f"\nMixed-unit users (sug_ISF switches mmol↔mgdl mid-history, "
             f"corrected <{MMOL_THRESHOLD:.0f}→×{MMOL_TO_MGDL}): "
             f"{', '.join(summary['mixed_unit_users']) or 'none'}",
             f"\nlog-DynISF subgroup: median per-tick log-correlation "
             f"{summary['log_dynisf']['median_log_corr']}, median ratio IQR "
             f"{summary['log_dynisf']['median_ratio_iqr']}\n",
             "| user | formula | n ticks | mixed units | log corr | median ratio | ratio IQR |",
             "|---|---|---|---|---|---|---|"]
    for r in summary["per_user"]:
        lines.append(f"| {r['user']} | {r['formula']} | {r['n_ticks']} | "
                     f"{'yes' if r['mixed_units'] else ''} | {r['log_corr']} | "
                     f"{r['median_ratio']} | {r['ratio_iqr']} |")
    (OUT / "device_isf_validation.md").write_text("\n".join(lines))
    print(df.to_string(index=False))
    print(f"\n→ {OUT/'device_isf_validation.md'}")


if __name__ == "__main__":
    main()
