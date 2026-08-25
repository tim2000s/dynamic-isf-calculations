"""Sensitivity against glucose and dose together, and whether the two combine cleanly.

Both equations are separable. v1 multiplies a dose term by a glucose scaler, v2
does the same with different exponents, and in neither does the glucose term
depend on how much insulin a person uses. Written out, both assert

    sensitivity = f(dose) x g(glucose)

with no term in which the two axes meet. That assumption has not been tested, and
it is the part of the question about whether the three quantities hold a clear
relationship rather than two separate ones.

Three measurements follow.

Separability. Fitting glucose and dose terms together with a term in which they
multiply gives an interaction coefficient. Both equations predict zero.

Whether the glucose term depends on the dose. Each person's own glucose exponent
is regressed on their daily dose. Both equations predict a flat line.

The surface itself. Sensitivity is measured in each combination of dose band and
glucose band, expressed relative to each person's own overall sensitivity, and
placed beside the surfaces the two equations predict for the same cells. A
separable law has identical rows across dose bands, because the glucose profile
does not change with dose.

    python3 -m inv009.joint_surface
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
BG_BANDS = [(90, 120, "90-120"), (120, 150, "120-150"), (150, 190, "150-190"),
            (190, 300, "190-300")]
TDD_BANDS = [(0, 25, "under 25"), (25, 40, "25-40"), (40, 60, "40-60"),
             (60, 1000, "over 60")]


def _fit(d: pd.DataFrame, extras: dict[str, np.ndarray]):
    """Fall on insulin action plus named interaction terms. Returns a coefficient map."""
    named = [("const", np.ones(len(d))), ("a_pre", d.a_pre.to_numpy(float))]
    named += list(extras.items())
    named += [("bg0", d.bg0.to_numpy() - 100.0), ("pre_slope", d.pre_slope.to_numpy(float)),
              ("bg_m60", d.bg_m60.to_numpy() - 100.0),
              ("bg_m120", d.bg_m120.to_numpy() - 100.0)]
    for hh in sorted(d.hour.unique())[1:]:
        named.append((f"hour_{hh}", (d.hour == hh).to_numpy(float)))
    X = np.column_stack([v for _, v in named])
    idx = {n: i for i, (n, _) in enumerate(named)}
    y = d["drop"].to_numpy(float)
    xtx = np.linalg.pinv(X.T @ X)
    beta = xtx @ (X.T @ y)
    resid = y - X @ beta
    n, k = X.shape
    meat = (X * (resid ** 2)[:, None]).T @ X
    cov = xtx @ meat @ xtx * (n / max(n - k, 1))
    return ({n: float(beta[i]) for n, i in idx.items()},
            {n: float(np.sqrt(max(cov[i, i], 0.0))) for n, i in idx.items()})


def separability(studies=tuple(config.COHORTS)) -> pd.DataFrame:
    """Per person: a glucose exponent, a dose exponent, and the term where they meet."""
    rows = []
    for study in studies:
        w = data.load(study)
        if w.empty:
            continue
        w = w[(w.tdd_blend > 0) & (w.tdd_u > 0) & (w.bg0 > 0)].copy()
        # Named in full. Short column names collide with DataFrame methods:
        # lt is less-than, and shift and corr have caught this project already.
        w["log_bg_ratio"] = np.log(REF_BG / w.bg0)
        w["log_tdd_ratio"] = np.log(w.tdd_blend / w.tdd_u)
        for sid, d in w.groupby("subject_id", sort=False):
            d = d.dropna(subset=["a_pre", "drop", "log_bg_ratio", "log_tdd_ratio"]
                         + CONTROLS)
            if len(d) < 80 or d.a_pre.std() < config.MIN_A_SD:
                continue
            a = d.a_pre.to_numpy(float)
            lg = d["log_bg_ratio"].to_numpy(float)
            lt = d["log_tdd_ratio"].to_numpy(float)
            b, se = _fit(d, {"x_g": a * lg, "x_t": a * lt, "x_gt": a * lg * lt,
                             "main_bg": lg, "main_tdd": lt})
            s = b["a_pre"]
            if s <= 0:
                continue
            rows.append(dict(subject_id=sid, study=study, n=len(d), s=s,
                             tdd_u=float(d.tdd_u.iloc[0]),
                             k_glucose=b["x_g"] / s, se_glucose=se["x_g"] / s,
                             k_dose=b["x_t"] / s, se_dose=se["x_t"] / s,
                             k_joint=b["x_gt"] / s, se_joint=se["x_gt"] / s))
    return pd.DataFrame(rows)


def band_profile(studies=("Loop", "ReplaceBG", "IOBP2", "DCLP3", "DCLP5", "PEDAP"),
                 ref: str = "120-150") -> pd.DataFrame:
    """Sensitivity by glucose band, imposing no shape, from one fit per person.

    Insulin action is interacted with glucose band indicators instead of with a
    log of glucose. Same rows, same controls and same single regression as the
    power-law version, with the functional form removed.

    This is what settles a disagreement between the two. The power-law exponent
    comes out near zero and changes sign with the glucose range selected, which
    reads as a weak effect. The band profile shows why: the relationship is not
    a power law. Sensitivity is materially lower below 120 mg/dL and close to
    flat from there to 300, so fitting a monotone curve across the whole range
    averages a step into a slope of almost nothing.
    """
    rows = []
    for study in studies:
        w = data.load(study)
        if w.empty:
            continue
        w = w[w.bg0 > 0].copy()
        for sid, d in w.groupby("subject_id", sort=False):
            d = d.dropna(subset=["a_pre", "drop"] + CONTROLS)
            if len(d) < 150 or d.a_pre.std() < config.MIN_A_SD:
                continue
            a = d.a_pre.to_numpy(float)
            extras, ok = {}, True
            for lo, hi, band in BG_BANDS:
                if band == ref:
                    continue
                m = d.bg0.between(lo, hi).to_numpy(float)
                if m.sum() < 25:
                    ok = False
                    break
                extras[f"a_x_{band}"] = a * m
            if not ok:
                continue
            b, _ = _fit(d, extras)
            s = b["a_pre"]
            if s <= 0:
                continue
            r = dict(subject_id=sid, study=study, tdd_u=float(d.tdd_u.iloc[0]), n=len(d))
            for lo, hi, band in BG_BANDS:
                r[band] = 1.0 if band == ref else 1.0 + b[f"a_x_{band}"] / s
            rows.append(r)
    return pd.DataFrame(rows)


def surface(studies=("Loop", "ReplaceBG", "IOBP2", "DCLP3", "DCLP5", "PEDAP")) -> pd.DataFrame:
    """Sensitivity in each glucose band, relative to that person's own overall value."""
    rows = []
    for study in studies:
        w = data.load(study)
        if w.empty:
            continue
        for sid, d in w.groupby("subject_id", sort=False):
            d = d.dropna(subset=["a_pre", "drop"] + CONTROLS)
            if len(d) < 120 or d.a_pre.std() < config.MIN_A_SD:
                continue
            b_all, _ = _fit(d, {})
            s_all = b_all["a_pre"]
            if s_all <= 0:
                continue
            tdd = float(d.tdd_u.iloc[0])
            for lo, hi, band in BG_BANDS:
                sub = d[d.bg0.between(lo, hi)]
                if len(sub) < 30 or sub.a_pre.std() < config.MIN_A_SD:
                    continue
                b, _ = _fit(sub, {})
                rows.append(dict(subject_id=sid, study=study, tdd_u=tdd,
                                 bg_band=band, bg_mid=(lo + hi) / 2,
                                 s_band=b["a_pre"], s_all=s_all,
                                 ratio=b["a_pre"] / s_all, n=len(sub)))
    return pd.DataFrame(rows)


def predicted_surface() -> pd.DataFrame:
    """What each equation says the same cells should hold."""
    rows = []
    for tlo, thi, tband in TDD_BANDS:
        tdd = np.array([min((tlo + thi) / 2, 90.0)])
        for lo, hi, band in BG_BANDS:
            bg = np.array([(lo + hi) / 2])
            ref = np.array([REF_BG])
            rows.append(dict(tdd_band=tband, bg_band=band,
                             v1=float(dynisf.isf_v1(bg, tdd)[0] / dynisf.isf_v1(ref, tdd)[0]),
                             v2=float(dynisf.isf_v2(bg, tdd)[0] / dynisf.isf_v2(ref, tdd)[0])))
    return pd.DataFrame(rows)


def main() -> int:
    config.ensure_dirs()
    res: dict = {}

    print("Is the relationship separable, as both equations assume?")
    S = separability()
    S.to_parquet(config.RESULTS / "inv009_separability.parquet", index=False)
    print(f"  {'cohort':11s} {'people':>6s} {'glucose k':>18s} {'dose k':>18s} "
          f"{'glucose x dose':>18s}")
    res["separability"] = []
    for label, d in list(S.groupby("study")) + [("ALL", S)]:
        out = dict(label=label, n=int(len(d)))
        cells = []
        for name, val, se in (("glucose", "k_glucose", "se_glucose"),
                              ("dose", "k_dose", "se_dose"),
                              ("joint", "k_joint", "se_joint")):
            dd = d[np.isfinite(d[val]) & np.isfinite(d[se]) & (d[se] > 0)]
            p = stats.dersimonian_laird(dd[val].to_numpy(), dd[se].to_numpy()) or {}
            out[name] = float(p.get("b_re", np.nan))
            out[f"{name}_se"] = float(p.get("se_re", np.nan))
            out[f"{name}_p"] = float(p.get("p", np.nan))
            cells.append(f"{out[name]:+7.3f} (p{p.get('p', float('nan')):7.0e})")
        res["separability"].append(out)
        print(f"  {label:11s} {out['n']:6d} " + " ".join(f"{c:>18s}" for c in cells))
    print("\n  both equations predict the joint term is zero")

    print("\nDoes a person's glucose exponent depend on how much insulin they use?")
    d = S[np.isfinite(S.k_glucose) & (S.tdd_u > 0)]
    fit = stats.loglog(d.tdd_u, d.k_glucose - d.k_glucose.min() + 1e-6)
    rho = float(d[["tdd_u", "k_glucose"]].corr(method="spearman").iloc[0, 1])
    lo_t = d[d.tdd_u < d.tdd_u.median()].k_glucose.median()
    hi_t = d[d.tdd_u >= d.tdd_u.median()].k_glucose.median()
    res["glucose_k_vs_dose"] = dict(n=int(len(d)), spearman=rho,
                                    median_k_low_dose=float(lo_t),
                                    median_k_high_dose=float(hi_t))
    print(f"  n={len(d)}  Spearman(dose, glucose exponent) = {rho:+.3f}")
    print(f"  median glucose exponent: {lo_t:+.3f} below the median dose, "
          f"{hi_t:+.3f} above it")

    print("\nSensitivity by glucose band, with no shape imposed (relative to 120-150)")
    B = band_profile()
    B.to_parquet(config.RESULTS / "inv009_band_profile.parquet", index=False)
    res["band_profile"] = dict(n=int(len(B)))
    for lo, hi, band in BG_BANDS:
        res["band_profile"][band] = float(B[band].median())
        print(f"  {band:>9s}  {B[band].median():.2f}")
    print("  v1 predicts " + ", ".join(
        f"{b[2]} {float(dynisf.isf_v1(np.array([(b[0]+b[1])/2]), np.array([40.0]))[0] / dynisf.isf_v1(np.array([135.0]), np.array([40.0]))[0]):.2f}"
        for b in BG_BANDS))
    print("  v2 predicts " + ", ".join(
        f"{b[2]} {float(dynisf.isf_v2(np.array([(b[0]+b[1])/2]), np.array([40.0]))[0] / dynisf.isf_v2(np.array([135.0]), np.array([40.0]))[0]):.2f}"
        for b in BG_BANDS))

    print("\nThe measured surface, as a fraction of each person's own sensitivity")
    F = surface()
    F.to_parquet(config.RESULTS / "inv009_surface.parquet", index=False)
    F["tdd_band"] = stats.band_of(F.tdd_u, TDD_BANDS)
    piv = F[F.tdd_band != ""].groupby(["tdd_band", "bg_band"]).ratio.median().unstack()
    piv = piv.reindex(index=[b[2] for b in TDD_BANDS], columns=[b[2] for b in BG_BANDS])
    cnt = F[F.tdd_band != ""].groupby(["tdd_band", "bg_band"]).ratio.size().unstack()
    cnt = cnt.reindex(index=[b[2] for b in TDD_BANDS], columns=[b[2] for b in BG_BANDS])
    print(f"  {'dose (U/day)':>14s} " + " ".join(f"{b[2]:>10s}" for b in BG_BANDS))
    for band in [b[2] for b in TDD_BANDS]:
        cells = " ".join(f"{piv.loc[band, c]:10.2f}" if pd.notna(piv.loc[band, c])
                         else f"{'-':>10s}" for c in piv.columns)
        print(f"  {band:>14s} {cells}")
    print(f"  {'people per cell':>14s} " +
          " ".join(f"{int(cnt.loc[TDD_BANDS[0][2], c]):10d}" if pd.notna(cnt.loc[TDD_BANDS[0][2], c])
                   else f"{'-':>10s}" for c in cnt.columns))
    res["surface"] = {t: {b: (float(piv.loc[t, b]) if pd.notna(piv.loc[t, b]) else None)
                          for b in piv.columns} for t in piv.index}

    P = predicted_surface()
    res["predicted"] = P.to_dict("records")
    print("\n  what the equations predict for the same cells (identical down every column,")
    print("  because a separable law's glucose profile does not change with dose)")
    for name in ("v1", "v2"):
        pv = P.pivot(index="tdd_band", columns="bg_band", values=name)
        pv = pv.reindex(index=[b[2] for b in TDD_BANDS], columns=[b[2] for b in BG_BANDS])
        print(f"  {name}: " + "  ".join(f"{c} {pv.loc[TDD_BANDS[0][2], c]:.2f}"
                                        for c in pv.columns))

    (config.RESULTS / "inv009_joint_surface.json").write_text(
        json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
