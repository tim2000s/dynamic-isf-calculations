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
# How much net insulin may already be on board at T, as a multiple of the
# correcting dose. Without this the denominator is not the correction: in the
# closed-loop cohorts the action over six hours ran 3.4 to 4.2 times the dose
# given, the excess being meal insulin still working from before the window. The
# fall was then divided by all of it. Requiring the slate to be nearly clean makes
# the action and the dose the same quantity, which is the condition under which
# this estimator was validated against a known sensitivity.
MAX_IOB_FRAC_DOSE = 0.5
BG_LO, BG_HI = 150.0, 350.0
STEP_MIN = 30.0
MIN_EVENTS = 10

STUDIES = ("ReplaceBG", "Loop", "DCLP3", "DCLP5", "PEDAP", "IOBP2")

# The action model is the physiological one everywhere, including Loop.
#
# Loop users were previously given the model that best reproduced the insulin on
# board their app displayed, and 120 of 158 matched Walsh at three hours. That is
# a good description of the app's display and a poor one of what the insulin did.
# Computing action on a three-hour curve treats a dose given four hours before T
# as finished when it is still working, which understates IOB(T), understates
# action, and inflates the sensitivity that comes out. The app model is still
# computed alongside so the size of that difference is reported rather than
# assumed.
MODEL_ACTION = "oref_6h75"


def analyse(subject_id: str) -> dict | None:
    study = subject_id.split(":")[0]
    app_model = _loop_model(subject_id) if study == "Loop" else MODEL_ACTION
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

    iob = np.convolve(net, M.kernel(MODEL_ACTION))[:n]
    iob_app = np.convolve(net, M.kernel(app_model))[:n]
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
        action_app = iob_app[i] - iob_app[i + h] + (ncum[i + h] - ncum[i + 1])
        clean = abs(iob[i]) <= MAX_IOB_FRAC_DOSE * given
        # Whether the system went on pushing insulin after the correcting dose.
        # Episodes where it did are ones where glucose stayed up for a reason the
        # record does not show, and they bias the estimate down; the validated
        # bolus estimator excluded them, so the same condition is applied here on
        # net insulin rather than on boluses alone.
        after = ncum[i + h] - ncum[j]
        settled = after <= 0.5 * given
        rows.append((bg[i], bg[i + h], given, action,
                     bcum[j] - bcum[i], tcum[j] - tcum[i], float(settled), after,
                     iob[i], action_app, float(clean)))
    if len(rows) < MIN_EVENTS:
        return None
    a = np.array(rows, dtype=float)
    fall = a[:, 0] - a[:, 1]
    st_m = a[:, 6] > 0.5
    # An episode starting with a net deficit is one where delivery had been held
    # below the programme beforehand. The arithmetic handles it, but those episodes
    # begin with glucose already rising for that reason, so they are reported
    # separately rather than pooled on the assumption they behave the same.
    pos = st_m & (a[:, 8] >= 0)
    neg = st_m & (a[:, 8] < 0)
    cl = st_m & (a[:, 10] > 0.5)
    med = lambda m, col=3: (float(np.median(fall[m] / a[m, col]))
                            if m.sum() >= MIN_EVENTS else np.nan)
    return dict(subject_id=subject_id, study=study, n_events=len(a), tdd_u=tdd,
                n_settled=int(st_m.sum()),
                n_pos=int(pos.sum()), n_neg=int(neg.sum()),
                isf_settled=med(st_m), isf_pos=med(pos), isf_neg=med(neg),
                isf_app=med(st_m, 9), n_clean=int(cl.sum()),
                isf_clean=med(cl), isf_clean_given=med(cl, 2),
                action_over_given=float(np.median(a[st_m, 3] / a[st_m, 2])),
                action_over_given_clean=(float(np.median(a[cl, 3] / a[cl, 2]))
                                         if cl.sum() >= MIN_EVENTS else np.nan),
                iob_start=float(np.median(a[st_m, 8])),
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
    print(f"{'cohort':<12s}{'settled':>9s}{'ISF':>6s}{'K':>6s}"
          f"{'clean n':>9s}{'act/dose':>9s}{'ISF':>6s}{'K':>6s}{'K by dose':>11s}"
          f"{'entered':>8s}")
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
        for tag, col in (("settled", "isf_settled"), ("pos", "isf_pos"),
                         ("neg", "isf_neg"), ("app", "isf_app"),
                         ("clean", "isf_clean"), ("cleang", "isf_clean_given")):
            row[f"isf_{tag}"] = float(d[col].median())
            row[f"k_{tag}"] = float((d[col] * d.tdd_u).median())
        row["n_settled"] = int(d.n_settled.sum())
        row["n_pos"] = int(d.n_pos.sum())
        row["n_neg"] = int(d.n_neg.sum())
        row["iob_start"] = float(d.iob_start.median())
        f1 = lambda v, w=6: (f"%{w}.1f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        f0 = lambda v, w=6: (f"%{w}.0f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        row["n_clean"] = int(d.n_clean.sum())
        row["aog"] = float(d.action_over_given.median())
        row["aog_clean"] = float(d.action_over_given_clean.median())
        print(f"{s:<12s}{row['n_settled']:9,d}"
              f"{f1(row['isf_settled'])}{f0(row['k_settled'])}"
              f"{row['n_clean']:9,d}{f1(row['aog_clean'], 9)}"
              f"{f1(row['isf_clean'])}{f0(row['k_clean'])}"
              f"{f0(row['k_cleang'], 11)}{f0(e, 8)}")
    (config.RESULTS / "inv009_correction_routes.json").write_text(json.dumps(out, indent=1))
    print("\n'temp basal' is the share of the correcting dose that arrived as basal")
    print("above the programmed rate rather than as a bolus.")
    print("'IOB>=0' and 'IOB<0' split on whether the episode opened with a net")
    print("insulin deficit; 'app model' is the same episodes computed on the curve")
    print("each Loop user's app displayed, which is the comparison, not the answer.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
