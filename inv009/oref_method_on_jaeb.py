"""The oref construction, run unchanged on the Jaeb cohorts.

The earlier work returned a ratio of measured to entered sensitivity of 0.41. This
work returns 0.30 to 0.92 with a median of 0.58 on the same question. Cohort is one
candidate explanation and construction is another, and only the second can be
tested by holding it fixed.

The oref construction is:

    observed sensitivity at T = [ BG(T) - BG(T+4h) ] / IOB(T)

for readings starting between 23:00 and 02:55, with a carbohydrate screen that
discards any window in which glucose rises faster than 2 mg/dL per five minutes,
and a floor on IOB(T) so the ratio has something to divide by.

Two choices in it differ from the one used here. The horizon is four hours rather
than six, and insulin on board at T does not finish acting inside four hours, so
the fall observed is smaller than the fall that insulin will produce. And the
denominator is insulin on board at T alone: anything delivered inside the window
lowers glucose and enlarges the numerator without appearing below the line.

Both are held fixed here and only the horizon is varied, so the size of that
choice can be read off.
"""
from __future__ import annotations

import json
import multiprocessing as mp

import numpy as np
import pandas as pd

from . import config, db, grid as gridmod, insulin_models as M
from .forward_isf import schedule_units

HOURS = (23, 0, 1, 2)
HORIZONS = (240, 300, 360, 420)
IOB_MIN = 0.5
RISE_MAX_PER_5MIN = 2.0
MODEL = "oref_6h75"
MIN_POINTS = 30


def analyse(subject_id: str) -> dict | None:
    st = db.streams(subject_id)
    g = gridmod.build_grid(st)
    if g is None or g.empty:
        return None
    n = len(g)
    tot = g.total_u.to_numpy(float)
    sched, _ = schedule_units(g, subject_id, st)
    net = tot - np.nan_to_num(sched)
    iob = np.convolve(net, M.kernel(MODEL))[:n]
    bg = g.cgm.to_numpy(float)
    hour = g.ts.dt.hour.to_numpy()
    tdd = float(np.nansum(tot) / max((n * config.GRID_MIN) / 1440.0, 1e-9))

    # Carbohydrate screen by shape rather than by log, as the oref work did: a
    # window is discarded if glucose ever climbs faster than 2 mg/dL per bin.
    rise = np.full(n, np.nan)
    rise[:-3] = (bg[3:] - bg[:-3]) / 3.0
    out = dict(subject_id=subject_id, study=subject_id.split(":")[0], tdd_u=tdd)
    for hz in HORIZONS:
        h = int(hz / config.GRID_MIN)
        if n < 2 * h:
            continue
        ok = np.zeros(n, dtype=bool)
        cand = np.flatnonzero(np.isin(hour, HOURS))
        cand = cand[(cand >= h) & (cand < n - h)]
        for i in cand:
            if iob[i] < IOB_MIN or not np.isfinite(bg[i]) or not np.isfinite(bg[i + h]):
                continue
            w = rise[i:i + h]
            if np.nanmax(w) > RISE_MAX_PER_5MIN:
                continue
            ok[i] = True
        if ok.sum() < MIN_POINTS:
            continue
        v = (bg[ok] - bg[np.flatnonzero(ok) + h]) / iob[ok]
        out[f"isf_{hz}"] = float(np.median(v))
        out[f"n_{hz}"] = int(ok.sum())
    return out if any(k.startswith("isf_") for k in out) else None


def main() -> int:
    config.ensure_dirs()
    subs = []
    for s in ("ReplaceBG", "Loop", "DCLP3", "DCLP5", "PEDAP", "IOBP2"):
        subs += db.subjects(s).subject_id.tolist()
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        rows = [r for r in pool.map(analyse, subs) if r]
    r = pd.DataFrame(rows)
    if r.empty:
        print("nothing survived")
        return 1
    r.to_parquet(config.RESULTS / "inv009_oref_method.parquet", index=False)

    ent = db.entered_isf().set_index("subject_id").isf
    print("The oref construction on the Jaeb cohorts: fall over the horizon divided")
    print("by insulin on board at the start, overnight, carbohydrate screened by shape\n")
    print(f"{'cohort':<12s}{'people':>7s}" + "".join(f"{hz//60}h".rjust(9) for hz in HORIZONS)
          + f"{'entered':>9s}" + "".join(f"r@{hz//60}h".rjust(8) for hz in HORIZONS))
    out = []
    for s, d in r.groupby("study"):
        e = float(ent.reindex(d.subject_id).dropna().median())
        vals = {hz: float(d[f"isf_{hz}"].median()) if f"isf_{hz}" in d else np.nan
                for hz in HORIZONS}
        row = dict(cohort=s, n_people=int(len(d)), entered=e,
                   isf={str(k): v for k, v in vals.items()},
                   ratio={str(k): (v / e if np.isfinite(e) and e else np.nan)
                          for k, v in vals.items()})
        out.append(row)
        f = lambda v, w=9: (f"%{w}.1f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        es = f"{e:9.1f}" if np.isfinite(e) else f"{'-':>9}"
        print(f"{s:<12s}{len(d):7d}" + "".join(f(vals[hz]) for hz in HORIZONS) + es
              + "".join((f"{vals[hz]/e:8.2f}" if np.isfinite(e) and e else f"{'-':>8}")
                        for hz in HORIZONS))
    (config.RESULTS / "inv009_oref_method.json").write_text(json.dumps(
        dict(horizons=list(HORIZONS), iob_min=IOB_MIN, hours=list(HOURS),
             by_cohort=out), indent=1, default=float))
    print("\nr@Nh is the ratio of measured to entered sensitivity at that horizon.")
    print("The oref work used four hours and reported 0.41.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
