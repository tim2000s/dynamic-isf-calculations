"""Does sensitivity change with the glucose it is acting on?

Both equations say yes, through a logarithmic term that makes a unit of insulin
do less when glucose is high. Testing that has one trap in it, and everything
here is built around avoiding it.

High glucose falls further than low glucose whatever insulin is doing, because
of mass action, renal clearance and simple regression to the mean, and because
near target the body actively defends against a further fall. So a model that
lets glucose explain the size of the fall directly, and then asks whether the
fall per unit of insulin also depends on glucose, is asking two different
questions and will attribute the first answer to the second. The additive
glucose term is therefore mandatory in every fit here. What carries the claim is
the INTERACTION between insulin action and glucose, which is the only thing that
means sensitivity itself changed.

    python3 -m inv009.glucose_axis
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from inv008 import dynisf

from . import config, data, stats

CONTROLS = ["bg0", "pre_slope", "bg_m60", "bg_m120"]
REF_BG = 100.0


def shapes(bg: np.ndarray, tdd: np.ndarray) -> dict[str, np.ndarray]:
    """Candidate multipliers on sensitivity, each normalised to one at 100 mg/dL.

    The two equation shapes are taken from the equations themselves rather than
    approximated, so what is being tested is what the code would actually do.
    """
    ref = np.full_like(bg, REF_BG)
    out = {"flat": np.ones_like(bg)}
    for k in (1.0, 2.0, 3.0):
        out[f"power_k{k:g}"] = (REF_BG / bg) ** k
    with np.errstate(divide="ignore", invalid="ignore"):
        out["v1_log"] = dynisf.isf_v1(bg, tdd) / dynisf.isf_v1(ref, tdd)
        out["v2_log"] = dynisf.isf_v2(bg, tdd) / dynisf.isf_v2(ref, tdd)
    return {k: np.nan_to_num(v, nan=1.0, posinf=1.0, neginf=1.0) for k, v in out.items()}


def _design(d: pd.DataFrame, extra: np.ndarray | None = None):
    cols = [np.ones(len(d)), d.a_pre.to_numpy(float),
            d.bg0.to_numpy() - 100.0, d.pre_slope.to_numpy(float),
            d.bg_m60.to_numpy() - 100.0, d.bg_m120.to_numpy() - 100.0]
    if extra is not None:
        cols.insert(2, extra)
    for hh in sorted(d.hour.unique())[1:]:
        cols.append((d.hour == hh).to_numpy(float))
    return np.column_stack(cols)


def linearised_k(study: str) -> pd.DataFrame:
    """The data's own glucose exponent, per subject.

    Sensitivity written as s * (100/BG)^k is, for modest departures from 100,
    s * (1 + k*log(100/BG)). So the interaction of insulin action with that log
    ratio has coefficient s*k, and the exponent is that divided by s. Same device
    as the within-person dose test, so the two are directly comparable.
    """
    w = data.load(study)
    if w.empty:
        return pd.DataFrame()
    w = w[w.bg0 > 0].copy()
    w["logr"] = np.log(REF_BG / w.bg0)
    rows = []
    for sid, d in w.groupby("subject_id", sort=False):
        d = d.dropna(subset=["a_pre", "drop", "logr"] + CONTROLS)
        if len(d) < config.MIN_WINDOWS or d.a_pre.std() < config.MIN_A_SD:
            continue
        inter = d.a_pre.to_numpy(float) * d.logr.to_numpy(float)
        X = _design(d, extra=inter)
        y = d["drop"].to_numpy(float)
        xtx = np.linalg.pinv(X.T @ X)
        beta = xtx @ (X.T @ y)
        resid = y - X @ beta
        n, kk = X.shape
        meat = (X * (resid ** 2)[:, None]).T @ X
        cov = xtx @ meat @ xtx * (n / max(n - kk, 1))
        s, c = float(beta[3]), float(beta[2])   # a_pre sits after the inserted interaction
        if s <= 0:
            continue
        rows.append(dict(subject_id=sid, study=study, n=int(n), s=s, c=c,
                         k=c / s, se_k=float(np.sqrt(max(cov[2, 2], 0))) / s,
                         tdd_u=d.tdd_u.iloc[0]))
    return pd.DataFrame(rows)


def shape_contest(study: str, folds: int = 5) -> pd.DataFrame:
    """Which glucose shape predicts the overnight fall best, out of sample.

    Blocks are contiguous in time rather than random, because windows overlap and
    a random split would put a window's near-twin in the other fold and make
    every shape look better than it is.
    """
    w = data.load(study)
    if w.empty:
        return pd.DataFrame()
    rows = []
    for sid, d in w.groupby("subject_id", sort=False):
        d = d.dropna(subset=["a_pre", "drop"] + CONTROLS).sort_values("t0")
        if len(d) < max(config.MIN_WINDOWS, folds * 15):
            continue
        sh = shapes(d.bg0.to_numpy(float), d.tdd_blend.to_numpy(float))
        block = np.minimum((np.arange(len(d)) * folds) // len(d), folds - 1)
        res = {"subject_id": sid, "study": study, "n": len(d), "tdd_u": d.tdd_u.iloc[0]}
        # Built once over the whole subject, then split by row. Building it inside
        # the fold would give train and test different hour dummies and different
        # widths whenever a fold happens to miss an hour.
        X_full = _design(d)
        a = d.a_pre.to_numpy(float)
        y = d["drop"].to_numpy(float)
        for name, g in sh.items():
            X = X_full.copy()
            X[:, 1] = a * g
            err = []
            for f in range(folds):
                tr, te = block != f, block == f
                if tr.sum() < 20 or te.sum() < 3:
                    continue
                try:
                    beta = np.linalg.pinv(X[tr].T @ X[tr]) @ (X[tr].T @ y[tr])
                except np.linalg.LinAlgError:
                    continue
                err.append(np.abs(y[te] - X[te] @ beta))
            res[f"mae_{name}"] = float(np.mean(np.concatenate(err))) if err else np.nan
        rows.append(res)
    return pd.DataFrame(rows)


def main() -> int:
    config.ensure_dirs()
    ks, contests = [], []
    for study in config.COHORTS:
        k = linearised_k(study)
        c = shape_contest(study)
        if not k.empty:
            ks.append(k)
        if not c.empty:
            contests.append(c)
        print(f"  {study:10s} k fitted {len(k):4d}   shape contest {len(c):4d}", flush=True)
    K = pd.concat(ks, ignore_index=True)
    C = pd.concat(contests, ignore_index=True)
    K.to_parquet(config.RESULTS / "inv009_glucose_k.parquet", index=False)
    C.to_parquet(config.RESULTS / "inv009_glucose_shapes.parquet", index=False)

    res: dict = {"k_by_study": [], "shapes_by_study": []}
    print("\nGlucose exponent inside each person (positive means less effective when high)")
    for label, d in list(K.groupby("study")) + [("ALL", K)]:
        d = d[np.isfinite(d.k) & np.isfinite(d.se_k) & (d.se_k > 0)]
        if len(d) < 10:
            continue
        p = stats.dersimonian_laird(d.k.to_numpy(), d.se_k.to_numpy()) or {}
        row = dict(label=label, n=int(len(d)), pooled_k=float(p.get("b_re", np.nan)),
                   se=float(p.get("se_re", np.nan)), pval=float(p.get("p", np.nan)),
                   i2=float(p.get("I2_pct", np.nan)), median_k=float(d.k.median()),
                   frac_positive=float((d.k > 0).mean()))
        res["k_by_study"].append(row)
        print(f"  {label:11s} n={row['n']:4d}  pooled k={row['pooled_k']:+.3f} "
              f"(se {row['se']:.3f}, p={row['pval']:.2g})  median {row['median_k']:+.3f}  "
              f"share positive {100 * row['frac_positive']:.0f}%")

    print("\nWhich glucose shape predicts the overnight fall best, out of sample")
    cols = [c for c in C.columns if c.startswith("mae_")]
    for label, d in list(C.groupby("study")) + [("ALL", C)]:
        med = {c[4:]: float(d[c].median()) for c in cols if d[c].notna().any()}
        if not med:
            continue
        best = min(med, key=med.get)
        flat = med.get("flat", np.nan)
        res["shapes_by_study"].append(dict(label=label, n=int(len(d)), best=best,
                                           medians=med,
                                           gain_over_flat=float(flat - med[best])))
        line = "  ".join(f"{k}={v:.2f}" for k, v in sorted(med.items(), key=lambda kv: kv[1]))
        print(f"  {label:11s} n={len(d):4d}  best={best:9s} (beats flat by "
              f"{flat - med[best]:+.2f} mg/dL)   {line}")
    (config.RESULTS / "inv009_glucose_axis.json").write_text(json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
