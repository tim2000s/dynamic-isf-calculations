"""Where does glucose actually land after a correction, against where it was aimed?

No sensitivity factor appears in this. The pump recorded the glucose it was given,
the target it was aiming at, and the correction it recommended. If glucose lands
at target then the entered factor is doing its job and any estimator that says
otherwise is wrong. If it lands short, the shortfall is the thing to measure.
"""
from __future__ import annotations

import json
import multiprocessing as mp

import numpy as np
import pandas as pd

from . import db, grid as gm, config

H = int(360 / config.GRID_MIN)
H4 = int(240 / config.GRID_MIN)


def wizard_full(sid):
    """db.streams() does not select the target columns, so fetch the row properly."""
    return db._q("SELECT ts_local, isf_mgdl_per_u, bg_input_mgdl, carb_input_g, "
                 "target_low_mgdl, target_high_mgdl, rec_correction_u "
                 "FROM studies.wizard WHERE subject_id=%s ORDER BY ts_local", (sid,))


_ENTERED = None


def entered_for(sid):
    """Entered sensitivity for a cohort with no bolus calculator.

    DCLP3, DCLP5 and PEDAP record settings on the case report form rather than at
    each dose, so there is no per-dose target either. The achieved sensitivity does
    not need one: it is the fall divided by the units given. Only the comparison
    against what was entered does, and 110 mg/dL is used as the target, which is
    the median the two calculator cohorts actually recorded.
    """
    global _ENTERED
    if _ENTERED is None:
        e = db.entered_isf()
        _ENTERED = dict(zip(e.subject_id, e.isf))
    return _ENTERED.get(sid, np.nan)


def one(sid):
    st = db.streams(sid)
    wiz = wizard_full(sid)
    if wiz.empty:
        isf0 = entered_for(sid)
        if not np.isfinite(isf0):
            return None
        wiz = None
    g = gm.build_grid(st)
    if g is None or len(g) < 2 * H:
        return None
    n = len(g)
    bg = g.cgm.to_numpy(float)
    bol = g.bolus_u.to_numpy(float)
    carbs = g.carbs_g.to_numpy(float)
    ts = g.ts.values

    if wiz is None:
        # No calculator: take every CGM point above threshold as a candidate and
        # let the dose and carbohydrate screens below do the selection.
        step = int(30 / config.GRID_MIN)
        cand = np.arange(H, n - H, step)
        cand = cand[(bg[cand] > 150) & (bg[cand] < 350)]
        if len(cand) < 5:
            return None
        idx = cand
        tgt = np.full(len(idx), 110.0)
        isf = np.full(len(idx), entered_for(sid))
        bgin = bg[idx]
        return _collect(sid, g, bg, bol, carbs, n, idx, tgt, isf, bgin)

    w = wiz.dropna(subset=["bg_input_mgdl", "isf_mgdl_per_u"]).copy()
    w = w[(w.bg_input_mgdl > 150) & (w.bg_input_mgdl < 350) & (w.isf_mgdl_per_u > 0)]
    if "carb_input_g" in w.columns:
        w = w[(w.carb_input_g.isna()) | (w.carb_input_g <= 0)]        # correction only
    if w.empty:
        return None
    tl = w.target_low_mgdl.to_numpy(float)
    th = w.target_high_mgdl.to_numpy(float)
    tgt = np.where(np.isfinite(th), (tl + th) / 2.0, tl)

    idx = np.searchsorted(ts, w.ts_local.values) - 1
    ok = (idx >= 0) & (idx < n - H)
    idx = idx[ok]; tgt = tgt[ok]
    isf = w.isf_mgdl_per_u.to_numpy(float)[ok]
    bgin = w.bg_input_mgdl.to_numpy(float)[ok]
    if len(idx) < 5:
        return None

    return _collect(sid, g, bg, bol, carbs, n, idx, tgt, isf, bgin)


def _collect(sid, g, bg, bol, carbs, n, idx, tgt, isf, bgin):
    ccum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(carbs))])
    bcum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(bol))])
    rows = []
    for i, j in enumerate(idx):
        if not np.isfinite(tgt[i]) or tgt[i] <= 0:
            continue
        if ccum[min(j + H, n)] - ccum[j] > 0:                # no carbs in the window
            continue
        if j >= H and ccum[j] - ccum[j - H] > 0:             # none before either
            continue
        given = bcum[min(j + 6, n)] - bcum[j]                # the dose, within 30 min
        if given < 0.5:
            continue
        extra = bcum[min(j + H, n)] - bcum[min(j + 6, n)]    # nothing else after
        if extra > 0.1:
            continue
        rows.append(dict(sid=sid, bg_start=bg[j], bg_wizard=bgin[i], target=tgt[i],
                         isf=isf[i], given=given,
                         bg_4h=bg[min(j + H4, n - 1)], bg_6h=bg[min(j + H, n - 1)]))
    return rows or None


if __name__ == '__main__':
    subs = []
    for s in ['ReplaceBG', 'Loop', 'DCLP3', 'DCLP5', 'PEDAP', 'IOBP2']:
        subs += db.subjects(s).subject_id.tolist()
    with mp.Pool(7, maxtasksperchild=8) as p:
        out = p.map(one, subs)
    r = pd.DataFrame([x for sub in out if sub for x in sub])
    r = r[np.isfinite(r.bg_6h) & np.isfinite(r.bg_start)]
    r["study"] = r.sid.str.split(":").str[0]
    r["predicted_6h"] = r.bg_start - r.isf * r.given
    r.to_parquet(config.RESULTS / 'inv009_correction_landing.parquet', index=False)
    print("Isolated corrections, no carbohydrate either side. Where did glucose land?\n")
    print(f"{'cohort':<12s}{'n':>6s}{'people':>7s}{'start':>7s}{'target':>8s}{'dose':>7s}"
          f"{'ISF':>6s}{'at 4h':>7s}{'at 6h':>7s}{'predicted':>10s}")
    for s, d in r.groupby("study"):
        print(f"{s:<12s}{len(d):6d}{d.sid.nunique():7d}{d.bg_start.median():7.0f}"
              f"{d.target.median():8.0f}{d.given.median():7.2f}{d.isf.median():6.0f}"
              f"{d.bg_4h.median():7.0f}{d.bg_6h.median():7.0f}{d.predicted_6h.median():10.0f}")
    print("\nPer correction, the fall achieved as a fraction of the fall aimed at:")
    for s, d in r.groupby("study"):
        aim = d.bg_start - d.target
        got = d.bg_start - d.bg_6h
        f = (got / aim).replace([np.inf, -np.inf], np.nan).dropna()
        print(f"  {s:<12s} median {f.median():.2f}   quartiles {f.quantile(.25):.2f} to {f.quantile(.75):.2f}")
    print("\nAnd the sensitivity that would have landed them on target:")
    rows = []
    for s, d in r.groupby("study"):
        need = ((d.bg_start - d.target) / d.given).replace([np.inf, -np.inf], np.nan).dropna()
        got = ((d.bg_start - d.bg_6h) / d.given).replace([np.inf, -np.inf], np.nan).dropna()
        per = d.groupby("sid").apply(
            lambda x: pd.Series({"isf": ((x.bg_start - x.bg_6h) / x.given).median(),
                                 "tdd": np.nan}), include_groups=False)
        rows.append(dict(cohort=s, n=int(len(d)), n_people=int(d.sid.nunique()),
                         entered=float(d.isf.median()), implied=float(need.median()),
                         achieved=float(got.median()),
                         ratio=float((got.median() / d.isf.median()))))
        print(f"  {s:<12s} entered {d.isf.median():.0f}   implied by the dose they gave "
              f"{need.median():.0f}   achieved {got.median():.0f}")
    (config.RESULTS / "inv009_correction_landing.json").write_text(json.dumps(rows, indent=1))
