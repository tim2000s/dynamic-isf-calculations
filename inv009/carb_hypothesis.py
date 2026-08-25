"""Is the dynamic ISF effect physiology, or carbohydrate the model has not caught?

This is the question the rest of the analysis measures around rather than
answers. Both explanations predict the same observation. If sensitivity genuinely
falls when glucose is high and when recent dose is large, the equations describe
something real. If instead carbohydrate from a recent meal is still absorbing,
unaccounted for, then glucose is high for that reason, insulin appears to achieve
less than it should, and recent dose is large because the person ate. The
observation is identical and the implication is not.

The archives can separate them because two cohorts recorded what people ate, so
time since the last meal is known rather than inferred. Three tests follow, each
with a prediction that differs between the two explanations.

Time since the last meal. Under the carbohydrate explanation the apparent glucose
dependence is strongest close to a meal and decays as absorption completes. Under
the physiological explanation it is constant, because a person's sensitivity does
not know when they last ate.

Recent carbohydrate as a control. Under the carbohydrate explanation, holding
recent grams constant removes the dose relationship, since dose is large when
someone has eaten. Under the physiological explanation the dose relationship
survives the control.

Cohorts that logged meals against cohorts that did not. Under the carbohydrate
explanation the effect is larger where meals cannot be screened out, because
unrecorded carbohydrate stays in the windows.

    python3 -m inv009.carb_hypothesis
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from . import config, data, stats, windows as windowmod

CONTROLS = ["bg0", "pre_slope", "bg_m60", "bg_m120"]
REF_BG = 100.0
CARB_BANDS = [(0, 4, "under 4 h"), (4, 6, "4 to 6 h"), (6, 9, "6 to 9 h"),
              (9, 14, "9 to 14 h"), (14, 1e9, "over 14 h")]


def base_screen(w: pd.DataFrame) -> pd.Series:
    """Everything the usual screen asks except the carbohydrate clearance.

    The clearance is what is being varied, so it cannot also be a fixed entry
    condition. Carbohydrate inside the window is still excluded: a window with a
    meal in it measures the meal.
    """
    h = config.HORIZON_MIN / config.GRID_MIN
    return (w.bg0.between(config.BG0_MIN, config.BG0_MAX)
            & (w.n_cgm >= config.MIN_CGM_FRACTION * h)
            & w.bg_end.notna() & w.tdd_blend.notna()
            & (w.carbs_in_g.fillna(1) <= 0)).fillna(False)


def _robust_fit(d: pd.DataFrame, inter: np.ndarray, extra: list[np.ndarray] | None = None):
    """Fit the fall on insulin action plus an interaction, return both and its error."""
    cols = [np.ones(len(d)), d.a_pre.to_numpy(float), inter,
            d.bg0.to_numpy() - 100.0, d.pre_slope.to_numpy(float),
            d.bg_m60.to_numpy() - 100.0, d.bg_m120.to_numpy() - 100.0]
    if extra:
        cols.extend(extra)
    for hh in sorted(d.hour.unique())[1:]:
        cols.append((d.hour == hh).to_numpy(float))
    X = np.column_stack(cols)
    y = d["drop"].to_numpy(float)
    xtx = np.linalg.pinv(X.T @ X)
    beta = xtx @ (X.T @ y)
    resid = y - X @ beta
    n, k = X.shape
    meat = (X * (resid ** 2)[:, None]).T @ X
    cov = xtx @ meat @ xtx * (n / max(n - k, 1))
    return float(beta[1]), float(beta[2]), float(np.sqrt(max(cov[2, 2], 0.0)))


MIN_LOGGED_DAY_G = 50.0


def glucose_k_by_carb_gap(studies=("Loop", "ReplaceBG"),
                          require_logged_day: bool = False,
                          bg_range: tuple[float, float] | None = None) -> pd.DataFrame:
    """The apparent glucose dependence, split by how long since the person ate.

    bg_range narrows the glucose window the exponent is measured over. It
    matters more than expected. Near target the body opposes a further fall, so
    those windows read as high sensitivity at low glucose and drag the exponent
    down. Excluding them raises every band by about 0.8, which is larger than
    the effect of carbohydrate being measured here, so the range has to be
    stated alongside any value.
    """
    rows = []
    for study in studies:
        w = data.load(study, screened=False)
        if w.empty or not bool(w.has_carb_stream.iloc[0]):
            continue
        w = w[base_screen(w) & w.carb_free_h.notna()].copy()
        if bg_range is not None:
            w = w[w.bg0.between(*bg_range)]
        if require_logged_day:
            # A long gap since the last recorded meal means either a genuine
            # fast or a meal that was never written down, and the two point
            # opposite ways. The preceding day tells them apart: windows more
            # than fourteen hours clear of a recorded meal show a median of
            # ZERO grams across the whole previous twenty-four hours, while
            # carrying the second highest insulin on board of any band. Nobody
            # eats nothing for a day and boluses through it. Requiring the day
            # to have been recorded at all removes those nights, and it has to
            # be done per night: filtering on how much a person records on
            # average leaves their unrecorded nights in.
            w = w[w.carbs_prev_24h >= MIN_LOGGED_DAY_G]
            if w.empty:
                continue
        w["logr"] = np.log(REF_BG / w.bg0)
        for lo, hi, label in CARB_BANDS:
            sub = w[(w.carb_free_h >= lo) & (w.carb_free_h < hi)]
            for sid, d in sub.groupby("subject_id", sort=False):
                d = d.dropna(subset=["a_pre", "drop", "logr"] + CONTROLS)
                if len(d) < 30 or d.a_pre.std() < config.MIN_A_SD:
                    continue
                s, c, se = _robust_fit(d, d.a_pre.to_numpy(float) * d.logr.to_numpy(float))
                if s <= 0:
                    continue
                rows.append(dict(subject_id=sid, study=study, band=label, n=len(d),
                                 s=s, k=c / s, se_k=se / s,
                                 median_gap=float(d.carb_free_h.median()),
                                 median_bg0=float(d.bg0.median())))
    return pd.DataFrame(rows)


def tdd_with_carb_control(studies=("Loop", "ReplaceBG")) -> pd.DataFrame:
    """The dose relationship, with and without recent carbohydrate held constant."""
    rows = []
    for study in studies:
        w = data.load(study)
        if w.empty or "carbs_prev_24h" not in w.columns:
            continue
        w = w[(w.tdd_blend > 0) & (w.tdd_u > 0)].copy()
        w["logratio"] = np.log(w.tdd_blend / w.tdd_u)
        for sid, d in w.groupby("subject_id", sort=False):
            d = d.dropna(subset=["a_pre", "drop", "logratio", "carbs_prev_24h",
                                 "carbs_prev_8h"] + CONTROLS)
            if len(d) < config.MIN_WINDOWS or d.a_pre.std() < config.MIN_A_SD:
                continue
            a = d.a_pre.to_numpy(float)
            inter = a * d.logratio.to_numpy(float)
            lr = d.logratio.to_numpy(float)
            c24 = d.carbs_prev_24h.to_numpy(float)
            c8 = d.carbs_prev_8h.to_numpy(float)
            if np.std(c24) < 1.0:
                continue
            s0, b0, se0 = _robust_fit(d, inter, extra=[lr])
            # Recent grams enter directly and as an interaction with insulin, so
            # the control removes both a shift in the fall and a shift in what a
            # unit achieves.
            s1, b1, se1 = _robust_fit(d, inter, extra=[lr, c24, c8, a * c24, a * c8])
            if s0 <= 0 or s1 <= 0:
                continue
            rows.append(dict(subject_id=sid, study=study, n=len(d),
                             b_raw=b0 / s0, se_raw=se0 / s0,
                             b_adj=b1 / s1, se_adj=se1 / s1,
                             carbs24_sd=float(np.std(c24))))
    return pd.DataFrame(rows)


def sensitivity_by_carb_gap(studies=("Loop", "ReplaceBG")) -> pd.DataFrame:
    """Measured sensitivity itself, split by how long since the person ate."""
    rows = []
    for study in studies:
        w = data.load(study, screened=False)
        if w.empty or not bool(w.has_carb_stream.iloc[0]):
            continue
        w = w[base_screen(w) & w.carb_free_h.notna()].copy()
        for lo, hi, label in CARB_BANDS:
            sub = w[(w.carb_free_h >= lo) & (w.carb_free_h < hi)]
            for sid, d in sub.groupby("subject_id", sort=False):
                d = d.dropna(subset=["a_pre", "drop"] + CONTROLS)
                if len(d) < 30 or d.a_pre.std() < config.MIN_A_SD:
                    continue
                s, _, _ = _robust_fit(d, np.zeros(len(d)))
                rows.append(dict(subject_id=sid, study=study, band=label, n=len(d), s=s))
    return pd.DataFrame(rows)


def _pool(d: pd.DataFrame, val: str, se: str) -> dict:
    d = d[np.isfinite(d[val]) & np.isfinite(d[se]) & (d[se] > 0)]
    if len(d) < 10:
        return {}
    p = stats.dersimonian_laird(d[val].to_numpy(), d[se].to_numpy()) or {}
    return dict(n=int(len(d)), pooled=float(p.get("b_re", np.nan)),
                se=float(p.get("se_re", np.nan)), pval=float(p.get("p", np.nan)),
                median=float(d[val].median()))


def main() -> int:
    config.ensure_dirs()
    res: dict = {}

    print("TEST 1  Apparent glucose dependence, by time since the last meal")
    print("  v1's scaler corresponds to k about +0.62 and v2's to about +1.77\n")
    arms = (("every night, glucose 90 to 300", dict(), "glucose_k_by_gap"),
            ("only nights where the day's food was recorded",
             dict(require_logged_day=True), "glucose_k_by_gap_logged"),
            ("glucose 120 to 220, near-target windows excluded",
             dict(bg_range=(120.0, 220.0)), "glucose_k_by_gap_midrange"))
    for tag, kw, key in arms:
        K = glucose_k_by_carb_gap(**kw)
        if K.empty:
            continue
        K.to_parquet(config.RESULTS / f"inv009_{key}.parquet", index=False)
        res[key] = []
        print(f"  {tag}")
        print(f"  {'gap':>12s} {'people':>7s} {'pooled k':>9s} {'se':>6s} {'p':>9s} "
              f"{'median bg0':>11s}")
        for _, _, label in CARB_BANDS:
            d = K[K.band == label]
            p = _pool(d, "k", "se_k")
            if not p:
                continue
            p.update(band=label, median_bg0=float(d.median_bg0.median()),
                     median_gap=float(d.median_gap.median()))
            res[key].append(p)
            print(f"  {label:>12s} {p['n']:7d} {p['pooled']:+9.2f} {p['se']:6.2f} "
                  f"{p['pval']:9.1e} {p['median_bg0']:11.0f}")
        print()

    print("\nTEST 2  Sensitivity itself, by time since the last meal")
    S = sensitivity_by_carb_gap()
    S.to_parquet(config.RESULTS / "inv009_carb_sensitivity.parquet", index=False)
    res["sensitivity_by_gap"] = []
    print(f"  {'gap':>12s} {'people':>7s} {'median sensitivity':>19s}")
    for _, _, label in CARB_BANDS:
        d = S[S.band == label]
        if len(d) < 10:
            continue
        row = dict(band=label, n=int(len(d)), median_s=float(d.s.median()))
        res["sensitivity_by_gap"].append(row)
        print(f"  {label:>12s} {row['n']:7d} {row['median_s']:19.2f}")

    print("\nTEST 3  Dose relationship, before and after holding recent carbohydrate constant")
    T = tdd_with_carb_control()
    T.to_parquet(config.RESULTS / "inv009_carb_tdd.parquet", index=False)
    res["tdd_carb_control"] = []
    for label, d in list(T.groupby("study")) + [("ALL", T)]:
        raw, adj = _pool(d, "b_raw", "se_raw"), _pool(d, "b_adj", "se_adj")
        if not raw or not adj:
            continue
        row = dict(label=label, n=raw["n"], b_raw=raw["pooled"], se_raw=raw["se"],
                   b_adj=adj["pooled"], se_adj=adj["se"],
                   shift=adj["pooled"] - raw["pooled"],
                   retained=float(adj["pooled"] / raw["pooled"]) if raw["pooled"] else np.nan)
        res["tdd_carb_control"].append(row)
        print(f"  {label:11s} n={row['n']:4d}  unadjusted {row['b_raw']:+.3f} "
              f"(se {row['se_raw']:.3f})   adjusted {row['b_adj']:+.3f} "
              f"(se {row['se_adj']:.3f})   {100 * row['retained']:.0f}% retained")

    print("\nTEST 4  Cohorts that logged meals against cohorts that did not")
    res["by_logging"] = []
    glu = json.loads((config.RESULTS / "inv009_glucose_axis.json").read_text())
    for r in glu["k_by_study"]:
        if r["label"] == "ALL":
            continue
        logs = config.COHORTS[r["label"]]["carbs"]
        res["by_logging"].append(dict(study=r["label"], logs_carbs=logs,
                                      pooled_k=r["pooled_k"], n=r["n"]))
    for logs in (True, False):
        sel = [x for x in res["by_logging"] if x["logs_carbs"] is logs]
        if sel:
            print(f"  {'logged meals' if logs else 'did not log meals':>18s}: "
                  f"pooled k by cohort " +
                  ", ".join(f"{x['study']} {x['pooled_k']:+.1f}" for x in sel))

    (config.RESULTS / "inv009_carb_hypothesis.json").write_text(
        json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
