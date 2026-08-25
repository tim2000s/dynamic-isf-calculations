"""A non-parametric second opinion on the dose exponent.

Everything else here fits a shape and reports its parameter, which means the
answer is only as good as the shape assumed. A gradient boosted model assumes
none: it is free to find that sensitivity falls with dose, rises with it, does so
in steps, or does not depend on it at all.

The quantity to read is NOT how important daily dose is on its own. Dose has a
large main effect on how far glucose falls overnight that has nothing to do with
sensitivity: people on more insulin are different people. Sensitivity is the
multiplier ON INSULIN, so what carries it is the INTERACTION between insulin
action and dose. SHAP interaction values give exactly that, and the slope of the
implied multiplier against log dose is the same exponent the regressions report,
arrived at without assuming it exists.

    python3 -m inv009.ml_shap
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from . import config, data, stats

FEATURES = ["a_pre", "log_tdd_u", "log_tdd_blend", "bg0", "pre_slope",
            "bg_m60", "bg_m120", "hour", "age", "iob0"]
# Interaction values cost O(features squared) per tree per row, so this is a
# few thousand rather than the whole set. The quantity read off them is a median
# within six dose bands, which a few thousand rows pins down perfectly well.
SAMPLE = 4_000
SEED = 0


def build() -> pd.DataFrame:
    frames = []
    for study in config.COHORTS:
        w = data.load(study)
        if w.empty:
            continue
        keep = ["subject_id", "study", "drop", "a_pre", "tdd_u", "tdd_blend", "bg0",
                "pre_slope", "bg_m60", "bg_m120", "hour", "age", "iob0", "closed_loop"]
        frames.append(w[keep])
    d = pd.concat(frames, ignore_index=True)
    d = d[(d.tdd_u > 0) & (d.tdd_blend > 0)].dropna(subset=["drop", "a_pre", "bg0"])
    d["log_tdd_u"] = np.log(d.tdd_u)
    d["log_tdd_blend"] = np.log(d.tdd_blend)
    return d.dropna(subset=FEATURES)


def main() -> int:
    import lightgbm as lgb
    import shap
    from sklearn.model_selection import GroupKFold

    config.ensure_dirs()
    d = build()
    print(f"  {len(d):,} windows from {d.subject_id.nunique():,} people")

    X, y = d[FEATURES], d["drop"].to_numpy(float)
    groups = d.subject_id.to_numpy()
    gkf = GroupKFold(n_splits=5)
    maes = []
    model = None
    for tr, te in gkf.split(X, y, groups):
        m = lgb.LGBMRegressor(n_estimators=300, learning_rate=0.06, num_leaves=31,
                              min_child_samples=80, subsample=0.8, colsample_bytree=0.8,
                              random_state=SEED, verbose=-1)
        m.fit(X.iloc[tr], y[tr])
        maes.append(float(np.mean(np.abs(y[te] - m.predict(X.iloc[te])))))
        model = m
    print(f"  out-of-person error {np.mean(maes):.2f} mg/dL")

    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(d), size=min(SAMPLE, len(d)), replace=False)
    Xs = X.iloc[idx]
    inter = shap.TreeExplainer(model).shap_interaction_values(Xs)
    ia, it = FEATURES.index("a_pre"), FEATURES.index("log_tdd_u")
    ib = FEATURES.index("bg0")

    # The interaction term is symmetric and split across the pair, so both halves
    # are added back to get the whole of it.
    inter_tdd = inter[:, ia, it] + inter[:, it, ia]
    inter_bg = inter[:, ia, ib] + inter[:, ib, ia]
    a = Xs.a_pre.to_numpy(float)
    lt = Xs.log_tdd_u.to_numpy(float)

    # The interaction is a contribution in mg/dL. Divided by the insulin it acted
    # on, it is a shift in sensitivity, in mg/dL per unit, at that dose.
    ok = np.abs(a) > 0.25
    shift = inter_tdd[ok] / a[ok]
    dose = np.exp(lt[ok])
    band = pd.qcut(dose, 6, duplicates="drop")
    # Named sens_shift, not shift: a column called shift shadows DataFrame.shift.
    prof = pd.DataFrame({"dose": dose, "sens_shift": shift}).groupby(band, observed=True).agg(
        dose=("dose", "median"), sens_shift=("sens_shift", "median"),
        n=("sens_shift", "size"))

    # A sensitivity that scales as dose^b has d(sensitivity)/d(log dose) = b *
    # sensitivity, so regressing the shift on the level gives the exponent back.

    res = dict(n_windows=int(len(d)), n_subjects=int(d.subject_id.nunique()),
               cv_mae=float(np.mean(maes)),
               importance={f: float(v) for f, v in
                           zip(FEATURES, np.abs(inter).sum(axis=2).mean(axis=0))},
               dose_interaction_profile=[dict(dose=float(r.dose),
                                              sens_shift=float(r.sens_shift),
                                              n=int(r.n)) for r in prof.itertuples()],
               mean_abs_interaction_tdd=float(np.mean(np.abs(inter_tdd))),
               mean_abs_interaction_bg=float(np.mean(np.abs(inter_bg))))
    (config.RESULTS / "inv009_ml_shap.json").write_text(json.dumps(res, indent=2, default=float))

    print("\n  Sensitivity shift attributable to daily dose, by dose band")
    print(f"  {'dose (U/day)':>13s} {'shift (mg/dL per U)':>21s} {'n':>7s}")
    for r in prof.itertuples():
        print(f"  {r.dose:13.1f} {r.sens_shift:21.2f} {r.n:7d}")
    lo, hi = float(prof.sens_shift.iloc[0]), float(prof.sens_shift.iloc[-1])
    print(f"\n  Sensitivity falls as dose rises: {'yes' if hi < lo else 'no'}"
          f"  ({lo:+.2f} at the lowest doses to {hi:+.2f} at the highest)")
    print(f"  Interaction with dose is {res['mean_abs_interaction_tdd'] / max(res['mean_abs_interaction_bg'], 1e-9):.2f}"
          f" times the size of the interaction with glucose")
    return 0


if __name__ == "__main__":
    sys.exit(main())
