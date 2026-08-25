"""Put the equations to work and see how well they predict the night.

Everything so far has been about exponents. This is the practical question: if
you had used v1's sensitivity, or v2's, to predict how far glucose would fall
overnight, how wrong would you have been, against the alternatives.

Each candidate gets a per-subject intercept fitted on the training half. That is
deliberate and it favours the equations rather than handicapping them. Overnight,
most insulin is basal and is there to offset endogenous glucose production, so a
prediction of the fall has to carry that offset somewhere; charging it to the
sensitivity factor would be blaming an equation for something it was never
supposed to supply. The intercept absorbs it, and what remains for the
sensitivity to get right is the part of the fall that scales with insulin.

Split is by time, not at random, so the test half is the future of the training
half and overlapping windows cannot leak across it.

    python3 -m inv009.head_to_head
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from inv008 import dynisf

from . import config, data, db, stats

TRAIN_FRAC = 0.7
MIN_TEST = 20


def candidates(d: pd.DataFrame, train: np.ndarray, entered: float | None) -> dict:
    """Sensitivity in mg/dL per unit for each candidate, per window."""
    bg = d.bg0.to_numpy(float)
    tdd_w = d.tdd_blend.to_numpy(float)
    tdd_u = float(d.tdd_u.iloc[0])
    a = d.a_pre.to_numpy(float)
    y = d["drop"].to_numpy(float)

    out: dict[str, np.ndarray] = {}
    # What the equations say, with no freedom at all.
    out["v1"] = dynisf.isf_v1(bg, tdd_w)
    out["v2"] = dynisf.isf_v2(bg, tdd_w)
    # The rules of thumb, static across the night.
    out["rule_1800"] = np.full(len(d), 1800.0 / max(tdd_u, 1e-6))
    out["root_tdd"] = np.full(len(d), 355.0 / np.sqrt(max(tdd_u, 1e-6)))
    if entered is not None and np.isfinite(entered) and entered > 0:
        out["entered"] = np.full(len(d), float(entered))
    # The best a single number could do for this person, learned on the train half.
    denom = float(np.sum(a[train] ** 2))
    fitted = float(np.sum(a[train] * (y[train] - np.mean(y[train]))) / denom) if denom > 0 else np.nan
    out["fitted_flat"] = np.full(len(d), fitted)
    return out


def score_subject(d: pd.DataFrame, entered: float | None) -> dict | None:
    d = d.dropna(subset=["a_pre", "drop", "bg0", "tdd_blend", "tdd_u"]).sort_values("t0")
    if len(d) < config.MIN_WINDOWS:
        return None
    n_tr = int(len(d) * TRAIN_FRAC)
    if len(d) - n_tr < MIN_TEST:
        return None
    train = np.zeros(len(d), bool)
    train[:n_tr] = True
    test = ~train
    a = d.a_pre.to_numpy(float)
    y = d["drop"].to_numpy(float)

    res = dict(subject_id=d.subject_id.iloc[0], study=d.study.iloc[0],
               tdd_u=float(d.tdd_u.iloc[0]), age=float(d.age.iloc[0]),
               n_test=int(test.sum()), bg0=float(np.median(d.bg0)))
    for name, isf in candidates(d, train, entered).items():
        if not np.all(np.isfinite(isf)):
            continue
        pred_no_intercept = isf * a
        # The intercept is the only thing fitted, and only on the train half.
        b = float(np.mean(y[train] - pred_no_intercept[train]))
        err = y[test] - (pred_no_intercept[test] + b)
        res[f"mae_{name}"] = float(np.mean(np.abs(err)))
        res[f"bias_{name}"] = float(np.median(err))
    return res


def main() -> int:
    config.ensure_dirs()
    entered = db.entered_isf().set_index("subject_id").isf.to_dict()
    rows = []
    for study in config.COHORTS:
        w = data.load(study)
        if w.empty:
            continue
        for sid, d in w.groupby("subject_id", sort=False):
            r = score_subject(d, entered.get(sid))
            if r:
                rows.append(r)
        print(f"  {study:10s} scored {sum(1 for r in rows if r['study'] == study):4d}",
              flush=True)
    R = pd.DataFrame(rows)
    R.to_parquet(config.RESULTS / "inv009_head_to_head.parquet", index=False)

    names = [c[4:] for c in R.columns if c.startswith("mae_")]
    res: dict = {"by_study": [], "by_tdd_band": [], "overall": {}}

    def summarise(d: pd.DataFrame) -> dict:
        out = {}
        for n in names:
            col = f"mae_{n}"
            if col in d and d[col].notna().any():
                out[n] = dict(mae=float(d[col].median()),
                              bias=float(d[f"bias_{n}"].median()),
                              n=int(d[col].notna().sum()))
        return out

    res["overall"] = summarise(R)
    print("\nMedian absolute error of the predicted overnight fall, mg/dL "
          "(lower is better)")
    order = sorted(res["overall"], key=lambda n: res["overall"][n]["mae"])
    print(f"{'cohort':26s} " + "  ".join(f"{n:>11s}" for n in order))
    for study, d in R.groupby("study"):
        s = summarise(d)
        res["by_study"].append(dict(label=config.COHORTS[study]["label"], study=study, **s))
        print(f"{config.COHORTS[study]['label']:26s} " +
              "  ".join(f"{s[n]['mae']:11.2f}" if n in s else f"{'-':>11s}" for n in order))
    print(f"{'ALL':26s} " + "  ".join(f"{res['overall'][n]['mae']:11.2f}" for n in order))
    print("\nMedian bias (positive means the fall was bigger than predicted)")
    print(f"{'ALL':26s} " + "  ".join(f"{res['overall'][n]['bias']:11.2f}" for n in order))

    print("\nBy daily dose, median absolute error")
    R["tdd_band"] = stats.band_of(R.tdd_u, [(lo, hi, f"{lo}-{hi if hi < 999 else '+'}")
                                            for lo, hi in config.TDD_BANDS])
    print(f"{'band (U/day)':14s} {'n':>5s} " + "  ".join(f"{n:>11s}" for n in order))
    for band, d in R.groupby("tdd_band"):
        if not band or len(d) < 20:
            continue
        s = summarise(d)
        res["by_tdd_band"].append(dict(band=band, n=int(len(d)), **s))
        print(f"{band:14s} {len(d):5d} " +
              "  ".join(f"{s[n]['mae']:11.2f}" if n in s else f"{'-':>11s}" for n in order))

    # How often each candidate is the best one for a person.
    mae_cols = [f"mae_{n}" for n in order]
    best = R[mae_cols].idxmin(axis=1).str[4:].value_counts()
    res["best_per_subject"] = {k: int(v) for k, v in best.items()}
    print("\nBest candidate per person: " +
          ", ".join(f"{k} {v}" for k, v in best.items()))
    (config.RESULTS / "inv009_head_to_head.json").write_text(json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
