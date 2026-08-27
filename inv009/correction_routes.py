"""Corrections delivered by any route, not only the ones a person pressed.

Every estimator in this package that isolated a correction did so by looking for a
bolus. Control-IQ, the bionic pancreas and a do-it-yourself loop all correct a
raised glucose by running basal above the programmed rate instead, so a
bolus-only definition sees a fraction of the corrections in those cohorts and
none at all in the bionic pancreas.

The dose here is insulin delivered above the programmed basal within thirty
minutes of T, by whichever route it arrived. The denominator is the action of all
above-schedule insulin across T to T+6h, because a temporary basal delivers its
units over an hour or more and the units delivered up front would understate what
acted.

Reported first is the split by route, which is the size of what a bolus-only
definition was missing.
"""
from __future__ import annotations

import json
import multiprocessing as mp

import numpy as np
import pandas as pd

from . import config, db, grid as gridmod, insulin_models as M
from .forward_isf import schedule_units, _loop_model

HORIZON_MIN = 360.0
DOSE_WINDOW_MIN = 30.0
MIN_DOSE_FRAC_TDD = 0.02       # 0.8 U at a daily dose of 40, 0.27 U at 13.5
BG_LO, BG_HI = 150.0, 350.0
STEP_MIN = 30.0
MIN_EVENTS = 10

STUDIES = ("ReplaceBG", "Loop", "DCLP3", "DCLP5", "PEDAP", "IOBP2")


def analyse(subject_id: str) -> dict | None:
    study = subject_id.split(":")[0]
    model = _loop_model(subject_id) if study == "Loop" else "oref_6h75"
    st = db.streams(subject_id)
    g = gridmod.build_grid(st)
    if g is None or g.empty:
        return None
    n = len(g)
    h = int(HORIZON_MIN / config.GRID_MIN)
    dw = max(1, int(DOSE_WINDOW_MIN / config.GRID_MIN))
    step = max(1, int(STEP_MIN / config.GRID_MIN))
    if n < 2 * h:
        return None

    tot = g.total_u.to_numpy(float)
    bol = np.nan_to_num(g.bolus_u.to_numpy(float))
    sched, _ = schedule_units(g, subject_id, st)
    net = tot - np.nan_to_num(sched)
    tempb = net - bol                       # above-schedule basal, the other route

    tdd = float(np.nansum(tot) / max((n * config.GRID_MIN) / 1440.0, 1e-9))
    min_dose = max(MIN_DOSE_FRAC_TDD * tdd, 0.2)

    iob = np.convolve(net, M.kernel(model))[:n]
    ncum = np.concatenate([[0.0], np.cumsum(net)])
    bcum = np.concatenate([[0.0], np.cumsum(bol)])
    tcum = np.concatenate([[0.0], np.cumsum(np.clip(tempb, 0, None))])
    carbs = np.nan_to_num(g.carbs_g.to_numpy(float))
    if carbs.sum() <= 0:
        carbs = np.where(bol >= config.MEAL_BOLUS_FRAC_TDD * tdd, 1.0, 0.0)
    ccum = np.concatenate([[0.0], np.cumsum(carbs)])
    bg = g.cgm.to_numpy(float)

    rows = []
    for i in range(h, n - h, step):
        if not (BG_LO <= bg[i] <= BG_HI) or not np.isfinite(bg[i + h]):
            continue
        if ccum[i] - ccum[i - h] > 0 or ccum[i + h] - ccum[i] > 0:
            continue
        j = min(i + dw, n)
        given = ncum[j] - ncum[i]
        if given < min_dose:
            continue
        action = iob[i] - iob[i + h] + (ncum[i + h] - ncum[i + 1])
        if action < min_dose:
            continue
        # Whether the system went on pushing insulin after the correcting dose.
        # Episodes where it did are ones where glucose stayed up for a reason the
        # record does not show, and they bias the estimate down; the validated
        # bolus estimator excluded them, so the same condition is applied here on
        # net insulin rather than on boluses alone.
        after = ncum[i + h] - ncum[j]
        settled = after <= 0.5 * given
        rows.append((bg[i], bg[i + h], given, action,
                     bcum[j] - bcum[i], tcum[j] - tcum[i], float(settled), after))
    if len(rows) < MIN_EVENTS:
        return None
    a = np.array(rows, dtype=float)
    fall = a[:, 0] - a[:, 1]
    st_m = a[:, 6] > 0.5
    return dict(subject_id=subject_id, study=study, n_events=len(a), tdd_u=tdd,
                n_settled=int(st_m.sum()),
                isf_settled=float(np.median((fall[st_m] / a[st_m, 3]))) if st_m.sum() >= MIN_EVENTS else np.nan,
                bg_end_settled=float(np.median(a[st_m, 1])) if st_m.sum() >= MIN_EVENTS else np.nan,
                bg_start=float(np.median(a[:, 0])), bg_end=float(np.median(a[:, 1])),
                given=float(np.median(a[:, 2])), action=float(np.median(a[:, 3])),
                bolus_part=float(a[:, 4].sum()), temp_part=float(a[:, 5].sum()),
                isf_given=float(np.median(fall / a[:, 2])),
                isf_action=float(np.median(fall / a[:, 3])))


def main() -> int:
    config.ensure_dirs()
    subs = []
    for s in STUDIES:
        subs += db.subjects(s).subject_id.tolist()
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        rows = [r for r in pool.map(analyse, subs) if r]
    r = pd.DataFrame(rows)
    if r.empty:
        print("no events")
        return 1
    r.to_parquet(config.RESULTS / "inv009_correction_routes.parquet", index=False)

    ent = db.entered_isf().set_index("subject_id").isf
    print("Corrections counted by any route above the programmed basal\n")
    print(f"{'cohort':<12s}{'people':>7s}{'events':>8s}{'temp basal':>11s}"
          f"{'start':>7s}{'ISF all':>8s}{'settled n':>10s}{'ISF settled':>12s}"
          f"{'K settled':>10s}{'entered':>8s}")
    out = []
    for s, d in r.groupby("study"):
        share = d.temp_part.sum() / max(d.temp_part.sum() + d.bolus_part.sum(), 1e-9)
        e = float(ent.reindex(d.subject_id).dropna().median())
        row = dict(cohort=s, n_people=int(len(d)), n_events=int(d.n_events.sum()),
                   temp_share=float(share), bg_start=float(d.bg_start.median()),
                   bg_end=float(d.bg_end.median()),
                   isf=float(d.isf_action.median()), entered=e,
                   k=float((d.isf_action * d.tdd_u).median()))
        out.append(row)
        row["isf_settled"] = float(d.isf_settled.median())
        row["k_settled"] = float((d.isf_settled * d.tdd_u).median())
        row["n_settled"] = int(d.n_settled.sum())
        es = f"{e:8.0f}" if np.isfinite(e) else f"{'-':>8}"
        iss = (f"{row['isf_settled']:12.1f}" if np.isfinite(row['isf_settled'])
               else f"{'-':>12}")
        ks = (f"{row['k_settled']:10.0f}" if np.isfinite(row['k_settled'])
              else f"{'-':>10}")
        print(f"{s:<12s}{row['n_people']:7d}{row['n_events']:8,d}{100*share:10.0f}%"
              f"{row['bg_start']:7.0f}{row['isf']:8.1f}{row['n_settled']:10,d}{iss}{ks}{es}")
    (config.RESULTS / "inv009_correction_routes.json").write_text(json.dumps(out, indent=1))
    print("\n'temp basal' is the share of the correcting dose that arrived as basal")
    print("above the programmed rate rather than as a bolus.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
