#!/usr/bin/env python3
"""Canonical cohort and TDD definition for the entire analysis stack.

ONE definition that every downstream artefact (paper, advisor, Phase 4) uses.

TDD = the most inclusive AND defensible per-user estimate:
  • v6 / v7  →  basal_profile + Σ treatments insulin / recording-span days
                (direct from the user's pump-uploaded delivery events)
  The cohort is the OpenAPS Data Commons (v6 and v7) only; the Trio (v5) cohort, from donated
  Nightscout sites, was removed from public analyses on 23 November 2026.

Quality filters (identical across all artefacts):
  • ISF        ∈ [10, 300]   mg/dL
  • CR         ∈ [2, 50]     g/U
  • target_low ∈ [70, 130]   mg/dL
  • TDD        ∈ [5, 200]    U/day
  • n_days     ≥ 14          days of decision history

This module exposes:
  load_canonical_cohort() → DataFrame with one row per user
                            and columns user_id, cohort, isf, cr, target_low,
                            basal, tdd, tdd_method, n_days, dynisf_user,
                            formula, group, in_cohort
"""
from __future__ import annotations

import json
import warnings
import os
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2

warnings.filterwarnings("ignore")

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))

QUALITY = {
    "isf": (10, 300),
    "cr":  (2, 50),
    "target_low": (70, 130),
    "tdd": (5, 200),
    "n_days_min": 14,
}


def _load_profiles() -> pd.DataFrame:
    # OpenAPS Data Commons only (AndroidAPS uploads and oref0). The Trio cohort came from donated
    # Nightscout sites outside the Commons and is not used in public analyses (23 November 2026).
    v6 = pd.DataFrame(json.loads((ROOT / "user_profiles_v6_enriched.json").read_text()).values())
    v6["cohort"] = "v6_aaps_classic"
    v7 = pd.DataFrame(json.loads((ROOT / "user_profiles_v7_enriched.json").read_text()).values())
    v7["cohort"] = "v7_oref0"
    df = pd.concat([v6, v7], ignore_index=True)
    df = df.drop_duplicates("user_id", keep="last").reset_index(drop=True)
    # Keep only the columns we use; some v5 records have stray fields that collide on merge.
    keep_cols = [c for c in (
        "user_id", "cohort", "device", "is_mmol",
        "profile_isf_mean_mgdl", "profile_cr_mean", "profile_target_low_mgdl",
        "profile_total_basal", "profile_max_basal", "profile_mean_basal",
        "profile_n_basal_segments", "inferred_max_iob", "mean_tdd",
    ) if c in df.columns]
    return df[keep_cols].copy()


def _load_hourly_basal_total(uid_to_total: dict) -> dict:
    p = ROOT / "user_basal_profiles.json"
    if not p.exists():
        return uid_to_total
    raw = json.loads(p.read_text())
    for uid, rec in raw.items():
        rates = rec.get("hourly_rates")
        if rates and len(rates) == 24:
            uid_to_total[uid] = float(sum(rates))
    return uid_to_total


def _load_treatments_insulin_per_day(uid_to_pdy: dict) -> dict:
    p = ROOT / "user_treatments_tdd.json"
    if not p.exists():
        return uid_to_pdy
    raw = json.loads(p.read_text())
    for uid, rec in raw.items():
        span = rec.get("span_days", 0)
        total = rec.get("total_insulin_units", 0)
        if span and span > 0 and total > 0:
            uid_to_pdy[uid] = total / span
    return uid_to_pdy


def _per_user_db_metrics(conn) -> pd.DataFrame:
    """Per-user n_days, dynisf_frac, and Sigmoid/Log group seeds across v6/v7."""
    parts = []
    for table in ("oref_v6", "oref_v7"):
        sql = f"""
            SELECT user_id,
                   (MAX(ts_relative_sec)-MIN(ts_relative_sec))/86400.0 AS n_days,
                   AVG(CASE WHEN has_dynisf > 0 THEN 1.0 ELSE 0.0 END) AS dynisf_frac,
                   '{table}' AS src_table
            FROM {table} GROUP BY user_id
        """
        parts.append(pd.read_sql(sql, conn))
    df = pd.concat(parts, ignore_index=True)
    df = df.sort_values("n_days", ascending=False).drop_duplicates("user_id", keep="first")
    df["dynisf_user"] = df["dynisf_frac"] > 0.5
    return df


def load_canonical_cohort() -> pd.DataFrame:
    """Single canonical cohort frame — the one the paper, advisor, and Phase 4 share."""
    profiles = _load_profiles()

    # Resolve total basal per user — prefer extracted hourly sum, fall back to enriched.
    uid_to_basal = {}
    _load_hourly_basal_total(uid_to_basal)
    profiles["basal"] = profiles.apply(
        lambda r: float(uid_to_basal.get(r["user_id"], r.get("profile_total_basal") or 0)),
        axis=1)

    # Resolve treatments-derived insulin per day where available.
    uid_to_tx = {}
    _load_treatments_insulin_per_day(uid_to_tx)
    profiles["tx_per_day"] = profiles["user_id"].map(uid_to_tx)

    # Canonical TDD assembly — most inclusive AND defensible.
    def _canonical_tdd(r):
        # v6 / v7: use basal + treatments/day if treatments available.
        if pd.notna(r["tx_per_day"]) and r["tx_per_day"] > 0 and r["basal"] > 0:
            return r["basal"] + r["tx_per_day"], "treatments_plus_basal"
        # no treatments-derived dose: fall back to a live-extracted mean_tdd where one exists
        m = r.get("mean_tdd")
        if pd.notna(m) and float(m) > 0:
            return float(m), "live_extracted_mean_tdd"
        # last resort (should be empty given filters): no TDD
        return None, None

    out = profiles.apply(_canonical_tdd, axis=1, result_type="expand")
    out.columns = ["tdd", "tdd_method"]
    profiles = pd.concat([profiles, out], axis=1)

    # DB metrics (n_days, DynISF status)
    with psycopg2.connect("dbname=oref") as conn:
        db = _per_user_db_metrics(conn)
    df = profiles.merge(db, on="user_id", how="left")

    # Rename the three core profile fields
    df = df.rename(columns={
        "profile_isf_mean_mgdl": "isf",
        "profile_cr_mean": "cr",
        "profile_target_low_mgdl": "target_low",
    })

    # the per-person formula map covered Trio people only and went with the Trio cohort
    df["formula"] = None
    df["group"] = df["dynisf_user"].map({False: "no_dynisf", True: "dynisf_unknown"})

    # Quality filter
    keep = (
        df["isf"].between(*QUALITY["isf"])
        & df["cr"].between(*QUALITY["cr"])
        & df["target_low"].between(*QUALITY["target_low"])
        & df["tdd"].between(*QUALITY["tdd"])
        & (df["n_days"] >= QUALITY["n_days_min"])
    )
    df["in_cohort"] = keep
    return df


if __name__ == "__main__":
    df = load_canonical_cohort()
    n_total = int(df["in_cohort"].sum())
    print(f"Canonical cohort size: {n_total}")
    print()
    print("By cohort source × TDD method:")
    print(df[df["in_cohort"]].groupby(["cohort", "tdd_method"]).size().to_string())
    print()
    print("By DynISF group:")
    print(df[df["in_cohort"]].groupby("group").size().to_string())
    out = ROOT / "canonical_cohort.json"
    df[df["in_cohort"]].to_json(out, orient="records", indent=2)
    print(f"\nWrote {out}")
