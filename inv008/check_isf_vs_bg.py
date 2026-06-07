#!/usr/bin/env python3
"""Sanity check: does the per-window local ISF actually rise with BG?
Pools overnight clean windows across users, normalises local-ISF to each user's median,
and reports the median normalised local-ISF by BG band."""
from __future__ import annotations
import multiprocessing as mp, os, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd, psycopg2
warnings.filterwarnings("ignore")
ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd())); sys.path.insert(0, str(ROOT))
from inv008.phase1_convergence import _compute_rows, COL_MAP
NIGHT = set(range(0, 6)); TABLES = {"v5_trio": "oref_v5", "v6_aaps_classic": "oref_v6", "v7_oref0": "oref_v7"}


def run_user(args):
    uid, table = args
    cm = COL_MAP[table]
    conn = psycopg2.connect("dbname=oref")
    try:
        df = pd.read_sql(f"SELECT ts_relative_sec,cgm_mgdl,iob_iob,{cm['cob']} AS cob,"
                         f"sug_smb_units,hour FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL "
                         f"AND iob_iob IS NOT NULL ORDER BY ts_relative_sec", conn, params=(uid,))
    finally:
        conn.close()
    if len(df) < 1000: return None
    ts, keep, diob, dbg, trend = _compute_rows(df, table)
    hour = df["hour"].values.astype(int); bg = df["cgm_mgdl"].values.astype(float)
    m = keep & np.isin(hour, list(NIGHT)) & (diob >= 0.25)
    if m.sum() < 50: return None
    X = np.column_stack([np.ones(m.sum()), diob[m], trend[m]])
    beta, *_ = np.linalg.lstsq(X, dbg[m], rcond=None); c = beta[2]
    local = -(dbg[m] - c * trend[m]) / diob[m]; bgm = bg[m]
    ok = np.isfinite(local) & (local >= 5) & (local <= 600) & (bgm > 40) & (bgm < 360)
    local, bgm = local[ok], bgm[ok]
    if len(local) < 50: return None
    return pd.DataFrame({"bg": bgm, "norm_isf": local / np.median(local)})


def main():
    from canonical_cohort import load_canonical_cohort
    coh = load_canonical_cohort(); coh = coh[coh["in_cohort"]]
    work = [(r["user_id"], TABLES[r["cohort"]]) for _, r in coh.iterrows()]
    with mp.Pool(12) as pool:
        parts = [p for p in pool.map(run_user, work, chunksize=2) if p is not None]
    d = pd.concat(parts, ignore_index=True)
    bins = [40, 70, 90, 110, 140, 180, 360]
    d["band"] = pd.cut(d.bg, bins)
    g = d.groupby("band").agg(n=("norm_isf", "size"),
                              median_norm_isf=("norm_isf", "median"),
                              median_bg=("bg", "median"))
    print(f"{len(parts)} users, {len(d)} overnight windows. Local-ISF normalised to each "
          f"user's median (so 1.0 = that user's typical). If DynISF theory held, this would "
          f"FALL across rising BG bands.\n")
    print(f"{'BG band':>12} {'n':>8} {'med BG':>8} {'norm local-ISF':>16}")
    for band, r in g.iterrows():
        print(f"{str(band):>12} {int(r.n):>8} {r.median_bg:>8.0f} {r.median_norm_isf:>16.2f}")


if __name__ == "__main__":
    main()
