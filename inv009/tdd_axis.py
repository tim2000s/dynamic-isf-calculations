"""Does sensitivity scale with total daily dose, and with what exponent?

This is the question the two equations disagree about. v1 says sensitivity goes
as one over daily dose, v2 says one over its square. Between two people whose
doses differ by a factor of three, v1 predicts a threefold difference in
sensitivity and v2 a ninefold one, so the exponent is not a detail.

It is asked twice, because the two versions are not the same question.

Between people: does someone on eighty units a day have half the sensitivity of
someone on forty? That is what the 1800 rule asserts and what the equations were
fitted to.

Within a person: when this person's own recent dose runs high, does their
sensitivity fall in the same proportion? That is what the equations actually DO,
several times an hour, and nothing guarantees a law that holds across people
holds inside one. This is the version that has never been tested.

    python3 -m inv009.tdd_axis
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from . import config, data, stats

CONTROLS = ["bg0", "pre_slope", "bg_m60", "bg_m120"]
CANDIDATES = {"v1": -1.0, "v2": -2.0, "root_tdd": -0.5, "flat": 0.0}


def power_fit(tdd, s, se=None):
    """Fit s = K * TDD^b directly, without taking logs.

    Taking logs would mean dropping every subject whose estimated sensitivity
    came out negative, and those are not a random subset: they are the ones whose
    controller reacted hardest. Fitting on the natural scale keeps them, at the
    cost of a fit that weights large sensitivities more heavily.
    """
    tdd = np.asarray(tdd, float)
    s = np.asarray(s, float)
    m = np.isfinite(tdd) & np.isfinite(s) & (tdd > 0)
    if se is not None:
        se = np.asarray(se, float)
        m &= np.isfinite(se) & (se > 0)
        wt = 1.0 / se[m]
    else:
        wt = np.ones(int(m.sum()))
    if m.sum() < 20:
        return dict(k=np.nan, b=np.nan, n=int(m.sum()))
    t, y = tdd[m], s[m]
    wt = wt / wt.mean()

    def resid(p):
        return wt * (y - p[0] * np.power(t, p[1]))

    try:
        r = least_squares(resid, x0=[np.median(y) * np.median(t), -1.0],
                          bounds=([-1e6, -4.0], [1e6, 2.0]), max_nfev=2000)
        return dict(k=float(r.x[0]), b=float(r.x[1]), n=int(m.sum()))
    except Exception:
        return dict(k=np.nan, b=np.nan, n=int(m.sum()))


def between_subject(S: pd.DataFrame, col: str = "s_pre") -> list[dict]:
    """The exponent across people, per cohort and pooled."""
    out = []
    for label, d in list(S.groupby("study")) + [("ALL", S)]:
        d = d[d.tdd_u > 0]
        pos = d[d[col] > 0]
        wt = 1.0 / np.square(pos[f"se_{col[2:]}"].replace(0, np.nan))
        r = stats.loglog(pos.tdd_u, pos[col], w=wt)
        ci = stats.boot_slope_ci(pos.tdd_u, pos[col], w=wt)
        ru = stats.loglog(pos.tdd_u, pos[col])
        ciu = stats.boot_slope_ci(pos.tdd_u, pos[col])
        pw = power_fit(d.tdd_u, d[col], d[f"se_{col[2:]}"])
        # No Deming here. These sensitivities are noisy regression slopes and
        # daily dose is a well measured average, so the error is essentially all
        # in y, where it costs precision but no bias.
        out.append(dict(
            label=label, n_all=int(len(d)), n_positive=int(len(pos)),
            frac_positive=float((d[col] > 0).mean()),
            loglog_slope=r["slope"], loglog_se=r["se"], ci_lo=ci[0], ci_hi=ci[1],
            unweighted_slope=ru["slope"], unweighted_ci=[ciu[0], ciu[1]],
            power_b=pw["b"], power_k=pw["k"],
            summary=stats.describe_slope(r, ci)))
    return out


def within_subject(study: str, acol: str = "a_pre") -> pd.DataFrame:
    """The exponent inside each person, from how their own recent dose moves.

    Linearised: sensitivity is written as s * (tdd_blend / tdd_u)^b, and for
    small departures from a person's usual dose that is s * (1 + b*log ratio), so
    the interaction of insulin action with the log dose ratio has coefficient
    s*b. The exponent is that coefficient divided by the sensitivity itself.
    """
    w = data.load(study)
    if w.empty:
        return pd.DataFrame()
    w = w[(w.tdd_blend > 0) & (w.tdd_u > 0)].copy()
    w["logratio"] = np.log(w.tdd_blend / w.tdd_u)
    rows = []
    for sid, d in w.groupby("subject_id", sort=False):
        d = d.dropna(subset=[acol, "drop", "logratio"] + CONTROLS)
        if len(d) < config.MIN_WINDOWS or d[acol].std() < config.MIN_A_SD:
            continue
        a = d[acol].to_numpy(float)
        inter = a * d.logratio.to_numpy(float)
        cols = [np.ones(len(d)), a, inter, d.logratio.to_numpy(float),
                d.bg0.to_numpy() - 100.0, d.pre_slope.to_numpy(),
                d.bg_m60.to_numpy() - 100.0, d.bg_m120.to_numpy() - 100.0]
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
        s, c = float(beta[1]), float(beta[2])
        if s <= 0:
            continue                       # an exponent about a negative sensitivity is meaningless
        rows.append(dict(subject_id=sid, study=study, tdd_u=d.tdd_u.iloc[0],
                         n=int(n), s=s, c=c, se_c=float(np.sqrt(max(cov[2, 2], 0))),
                         b=c / s, se_b=float(np.sqrt(max(cov[2, 2], 0))) / s,
                         logratio_sd=float(d.logratio.std())))
    return pd.DataFrame(rows)


def main() -> int:
    config.ensure_dirs()
    S = pd.read_parquet(config.RESULTS / "inv009_user_isf.parquet")
    S = S[S.n_windows >= config.MIN_WINDOWS]
    res: dict = {}

    print("BETWEEN PEOPLE: exponent of sensitivity against daily dose")
    for col in ("s_pre", "s_net_pre", "s_total"):
        res[f"between_{col}"] = between_subject(S, col)
        print(f"\n  regressor: {col}")
        for r in res[f"between_{col}"]:
            print(f"    {r['label']:11s} n={r['n_positive']:4d}/{r['n_all']:4d} "
                  f"loglog {r['loglog_slope']:+.3f} [{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}]  "
                  f"unweighted {r['unweighted_slope']:+.3f}  "
                  f"power-fit(all) b={r['power_b']:+.3f}")

    print("\nWITHIN A PERSON: does their own recent dose move their sensitivity?")
    within = []
    for study in config.COHORTS:
        d = within_subject(study)
        if not d.empty:
            within.append(d)
    W = pd.concat(within, ignore_index=True) if within else pd.DataFrame()
    if not W.empty:
        W.to_parquet(config.RESULTS / "inv009_within_tdd.parquet", index=False)
        res["within"] = []
        for label, d in list(W.groupby("study")) + [("ALL", W)]:
            d = d[np.isfinite(d.b) & np.isfinite(d.se_b) & (d.se_b > 0)]
            if len(d) < 10:
                continue
            pooled = stats.dersimonian_laird(d.b.to_numpy(), d.se_b.to_numpy()) or {}
            row = dict(label=label, n=int(len(d)),
                       pooled_b=float(pooled.get("b_re", np.nan)),
                       pooled_se=float(pooled.get("se_re", np.nan)),
                       pooled_p=float(pooled.get("p", np.nan)),
                       i2=float(pooled.get("I2_pct", np.nan)),
                       frac_same_sign=float(pooled.get("frac_same_sign", np.nan)),
                       median_b=float(d.b.median()),
                       frac_negative=float((d.b < 0).mean()))
            res["within"].append(row)
            print(f"    {label:11s} n={row['n']:4d}  pooled b={row['pooled_b']:+.3f} "
                  f"(se {row['pooled_se']:.3f}, p={row['pooled_p']:.2g}, I2={row['i2']:.0f}%)  "
                  f"median {row['median_b']:+.3f}   (v1 predicts -1, v2 -2)")

    (config.RESULTS / "inv009_tdd_axis.json").write_text(json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
