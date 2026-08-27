"""Check the reconstructed insulin on board against what REPLACE-BG's pumps recorded.

Every sensitivity estimate in this work rests on insulin action reconstructed from
delivery records, and until now that reconstruction had been validated only in the
Loop cohort, against Loop's own recorded insulin on board, at r = 0.927. REPLACE-BG
carries the entire measured sensitivity result because it is the only cohort with
no algorithm intervening, and its bolus calculator records insulin on board in
232,329 rows from 204 people. That column was previously loaded as null.

Two things are being tested at once. Whether the reconstruction is right, and what
the recorded number means: a 2015 pump wizard conventionally counts bolus insulin
only, where Loop counts delivery net of the programmed basal. Both conventions are
fitted here rather than assumed, across the durations a pump of that era offered.

The recorded value precedes the bolus it was recommending, so that bolus is removed
from the reconstruction before comparing.
"""
from __future__ import annotations

import json
import multiprocessing as mp

import numpy as np
import pandas as pd

from . import config, db, grid as gridmod, insulin_models as M

MIN_RECORDS = 40
# A 2015 pump wizard used a curvilinear decay with a duration the user set,
# typically three to six hours, which is the Walsh family. The oref exponentials
# are included because they are what this work assumes everywhere else.
MODELS = ("walsh_3h", "walsh_4h", "walsh_5h", "walsh_6h",
          "oref_5h75", "oref_6h75", "oref_6h55", "oref_7h75")
CONVENTIONS = ("bolus_only", "net_of_schedule", "total_delivery")


def _series(g: pd.DataFrame, convention: str) -> np.ndarray:
    if convention == "bolus_only":
        return g.bolus_u.to_numpy(float)
    if convention == "net_of_schedule":
        return (g.total_u - g.sched_u.fillna(0.0)).to_numpy(float)
    return g.total_u.to_numpy(float)


def assess(subject_id: str) -> list[dict]:
    streams = db.streams(subject_id)
    wiz = streams["wizard"]
    if wiz.empty or wiz.iob_u.notna().sum() < MIN_RECORDS:
        return []
    g = gridmod.build_grid(streams)
    if g is None or g.empty:
        return []
    ts = g.ts.values
    bolus = g.bolus_u.to_numpy(float)

    w = wiz.dropna(subset=["iob_u"])
    idx = np.searchsorted(ts, w.ts_local.values) - 1
    ok = (idx >= 0) & (idx < len(g))
    idx, recorded = idx[ok], w.iob_u.to_numpy(float)[ok]
    if len(idx) < MIN_RECORDS:
        return []

    rows = []
    for conv in CONVENTIONS:
        u = _series(g, conv)
        if not np.isfinite(u).any():
            continue
        for m in MODELS:
            iob = np.convolve(u, M.kernel(m))[:len(u)]
            mine = iob[idx] - bolus[idx]
            d = mine - recorded
            rows.append(dict(
                subject_id=subject_id, convention=conv, model=m, n=int(len(idx)),
                mad=float(np.median(np.abs(d))), bias=float(np.median(d)),
                corr=float(np.corrcoef(mine, recorded)[0, 1]) if len(idx) > 3 else np.nan,
                recorded_p75=float(np.percentile(np.abs(recorded), 75))))
    return rows


def main() -> int:
    config.ensure_dirs()
    subs = db.subjects("ReplaceBG").subject_id.tolist()
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        out = pool.map(assess, subs)
    r = pd.DataFrame([x for sub in out for x in sub])
    if r.empty:
        print("no subjects with recorded insulin on board")
        return 1
    r.to_parquet(config.RESULTS / "inv009_iob_validation_replacebg.parquet", index=False)

    print(f"REPLACE-BG: reconstructed insulin on board against the pump's own record")
    print(f"{r.subject_id.nunique()} people, {int(r.n.iloc[0]):,}+ comparisons each\n")
    print(f"{'convention':<18s}{'model':<11s}{'median MAD':>11s}{'bias':>8s}{'corr':>7s}")
    g = (r.groupby(["convention", "model"])
           .agg(mad=("mad", "median"), bias=("bias", "median"), corr=("corr", "median"))
           .reset_index().sort_values("mad"))
    for _, x in g.iterrows():
        print(f"{x.convention:<18s}{x.model:<11s}{x.mad:11.3f}{x.bias:8.3f}{x['corr']:7.3f}")
    best = g.iloc[0]
    print(f"\nBest: {best.convention} on {best.model}, median absolute difference "
          f"{best.mad:.3f} U, correlation {best['corr']:.3f}")
    per = r[(r.convention == best.convention) & (r.model == best.model)]
    print(f"Across people: correlation quartiles "
          f"{per['corr'].quantile(.25):.3f} / {per['corr'].median():.3f} / "
          f"{per['corr'].quantile(.75):.3f}")
    (config.RESULTS / "inv009_iob_validation_replacebg.json").write_text(json.dumps(
        dict(n_people=int(r.subject_id.nunique()),
             best_convention=str(best.convention), best_model=str(best.model),
             mad=float(best.mad), corr=float(best["corr"]),
             table=g.to_dict("records")), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
