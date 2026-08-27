"""Run every sensitivity estimator against records whose true sensitivity is known.

Each estimator in this package returns a number. None of them can say whether that
number is the sensitivity that generated the data, because in a real record the
sensitivity is not observed. Here it is: records are simulated at a stated
sensitivity and each estimator is asked to recover it, running through the same
functions the real analyses call rather than through a copy.

Three arms. Open loop, where delivery is a fixed schedule plus the person's own
boluses. Reactive, where a controller adds insulin in proportion to how far glucose
sits above target, which is the confounding-by-indication case. And a harsh arm
with unannounced carbohydrate and error in how much of a dose acted, which is what
a real record looks like.

An estimator that recovers the truth in the open arm and loses it in the reactive
arm is measuring sensitivity and being confounded. One that loses it in the open
arm is not measuring sensitivity at all.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import config, forward_isf as F, synthetic as S

TRUE_ISF = 50.0
TDD = 40.0
N_SUBJECTS = 30
DAYS = 200
HORIZONS = (240, 360, 480)


def landing_estimate(g: pd.DataFrame, horizon_min: float = 360.0,
                     min_dose: float = 0.5, bg_lo: float = 150.0,
                     bg_hi: float = 350.0) -> float:
    """Fall per unit given on isolated corrections, the estimator that survived.

    Mirrors inv009.correction_landing: a dose at raised glucose with no
    carbohydrate either side and nothing further given, then the fall by the
    horizon divided by the units delivered.
    """
    n = len(g)
    h = int(horizon_min / config.GRID_MIN)
    bg = g.cgm.to_numpy(float)
    sched = g.sched_u.to_numpy(float)
    extra = (g.total_u.to_numpy(float) - sched)      # boluses and corrections
    carbs = np.nan_to_num(g.carbs_g.to_numpy(float))
    ccum = np.concatenate([[0.0], np.cumsum(carbs)])
    ecum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(extra))])
    out = []
    for i in range(h, n - h):
        if bg[i] < bg_lo or bg[i] > bg_hi or not np.isfinite(bg[i]):
            continue
        if ccum[i] - ccum[i - h] > 0 or ccum[i + h] - ccum[i] > 0:
            continue
        given = ecum[min(i + 6, n)] - ecum[i]
        if given < min_dose:
            continue
        if ecum[i + h] - ecum[min(i + 6, n)] > 0.1:
            continue
        if not np.isfinite(bg[i + h]):
            continue
        out.append((bg[i] - bg[i + h]) / given)
    return float(np.median(out)) if len(out) >= 10 else np.nan


def cohort(reactive: float = 0.0, unannounced: float = 0.0,
           action_noise: float = 0.0, seed: int = 7) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for k in range(N_SUBJECTS):
        g = S.simulate_subject(rng, tdd=TDD, isf=TRUE_ISF, days=DAYS,
                               reactive=reactive, unannounced=unannounced,
                               action_noise=action_noise)
        r = F.analyse(f"Synth:{k}", grid=g, model="oref_6h75") or {}
        r["landing"] = landing_estimate(g)
        for hz in HORIZONS:
            r[f"landing_{hz}"] = landing_estimate(g, horizon_min=hz)
        r["bg_median_all"] = float(np.nanmedian(g.cgm.to_numpy(float)))
        rows.append(r)
    return pd.DataFrame(rows)


def main() -> int:
    config.ensure_dirs()
    # Every arm carries unannounced carbohydrate. Without it the simulated person
    # sits at a median glucose of 100 and never rises far enough to correct, so the
    # estimator that depends on corrections has nothing to read and the comparison
    # is not a comparison. What varies between arms is the controller and the error
    # in insulin action, which is what the audit is about.
    arms = [("open loop", dict(unannounced=0.45)),
            ("reactive controller", dict(unannounced=0.45, reactive=0.02)),
            ("open loop + action error", dict(unannounced=0.45, action_noise=0.35)),
            ("reactive + action error",
             dict(unannounced=0.45, reactive=0.02, action_noise=0.35))]
    out = []
    print(f"True sensitivity {TRUE_ISF:.0f} mg/dL/U at a daily dose of {TDD:.0f} U, "
          f"{N_SUBJECTS} people, {DAYS} days each\n")
    print(f"{'arm':<36s}{'every point':>13s}{'at 200+':>10s}{'landing 6h':>12s}"
          f"{'median BG':>11s}")
    for label, kw in arms:
        d = cohort(**kw)
        row = dict(arm=label,
                   every_point=float(d.isf_ratio.median()) if "isf_ratio" in d else np.nan,
                   at_200=float(d.isf_head.median()) if "isf_head" in d else np.nan,
                   landing=float(d.landing.median()),
                   bg=float(d.bg_median_all.median()))
        for hz in HORIZONS:
            row[f"landing_{hz}"] = float(d[f"landing_{hz}"].median())
        out.append(row)
        f = lambda v, w=13: (f"%{w}.1f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        print(f"{label:<36s}{f(row['every_point'])}{f(row['at_200'], 10)}"
              f"{f(row['landing'], 12)}{row['bg']:11.0f}")

    print("\nThe landing estimator by horizon, open-loop arm "
          "(true sensitivity is still 50):\n")
    print(f"{'horizon':<12s}{'recovered':>11s}{'recovered / true':>19s}")
    o = out[0]
    for hz in HORIZONS:
        v = o[f"landing_{hz}"]
        print(f"{hz / 60:<12.0f}{v:11.1f}{v / TRUE_ISF:19.2f}")

    (config.RESULTS / "inv009_audit_estimators.json").write_text(json.dumps(
        dict(true_isf=TRUE_ISF, tdd=TDD, n_subjects=N_SUBJECTS, days=DAYS,
             arms=out), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
