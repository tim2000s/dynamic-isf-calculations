"""The sensitivity constant fitted from observed glucose response, by cohort.

A constant is fitted by asking which law of the form ISF = K / TDD^b best predicts
the overnight fall on periods the fit never saw. Every person receives their own
intercept, because overnight most insulin is basal offsetting hepatic glucose
output and no sensitivity factor should be asked to carry that.

The choice of exposure decides what the answer means, and an earlier version of
this work got it wrong. Fitting on `a_pre`, the action within the window from
insulin given before it opened, was chosen because pre-committed insulin cannot
respond to what happens inside the window, which protects against reverse
causation where an algorithm is dosing. It also discards the correction: a dose
given at the start of a window is delivered inside it, so on isolated correction
nights only 22 to 29% of the action is pre-window and the remainder is the basal
tail. The constant that came out of that fit, 880, was the coefficient on the
basal tail, and it was roughly half what a correction is observed to do.

Fitting on `a_tot` uses the whole dose. In the open-loop cohort that is close to
exogenous once each person carries their own intercept, since a correction is
decided from the glucose reading and then nothing further intervenes. Where an
algorithm is running it is not exogenous, and those cohorts are reported for
completeness rather than as estimates of sensitivity.

Scoring is a per-person 70/30 split in time, median absolute error on the held-out
tail, pooled as the median across people.
"""
from __future__ import annotations

import json

import numpy as np

from . import config
from .recommend import _subject_frames

SOURCE = "tdd_7d"
GRID = [round(x, 2) for x in np.arange(0.0, 1.65, 0.05)]
MIN_PEOPLE = 25
N_BOOT = 400
EXPOSURES = ("a_tot", "a_pre")


def fit_constant(frames, source: str, b: float, expo: str) -> float:
    """K in ISF = K / TDD^b, one number for everyone, each person their own intercept."""
    num = den = 0.0
    for _, _, d, n_tr in frames:
        tr = d.iloc[:n_tr]
        a = tr[expo].to_numpy(float)
        y = tr["drop"].to_numpy(float)
        x = a / np.power(tr[source].to_numpy(float), b)
        x = x - x.mean()
        y = y - y.mean()
        num += float(np.sum(x * y))
        den += float(np.sum(x * x))
    return num / den if den > 0 else np.nan


def score(frames, source: str, b: float, k: float, expo: str) -> float:
    errs = []
    for _, _, d, n_tr in frames:
        tr, te = d.iloc[:n_tr], d.iloc[n_tr:]
        isf_tr = k / np.power(tr[source].to_numpy(float), b)
        isf_te = k / np.power(te[source].to_numpy(float), b)
        off = float(np.mean(tr["drop"].to_numpy(float) - isf_tr * tr[expo].to_numpy(float)))
        e = te["drop"].to_numpy(float) - (isf_te * te[expo].to_numpy(float) + off)
        errs.append(np.median(np.abs(e)))
    return float(np.median(errs))


def curve(frames, expo: str, source: str = SOURCE):
    out = []
    for b in GRID:
        k = fit_constant(frames, source, b, expo)
        out.append((float(b), float(k), score(frames, source, b, k, expo)))
    return out


def at(c, b: float):
    row = min(c, key=lambda t: abs(t[0] - b))
    return row[1], row[2]


def main() -> int:
    config.ensure_dirs()
    frames = list(_subject_frames())
    by: dict[str, list] = {}
    for t in frames:
        by.setdefault(t[1], []).append(t)

    res = {"n_people": len(frames), "source": SOURCE,
           "horizon_min": config.HORIZON_MIN, "by_exposure": {}}
    rng = np.random.default_rng(20260826)

    for expo in EXPOSURES:
        pooled = curve(frames, expo)
        b_star, k_star, mae_star = min(pooled, key=lambda t: t[2])
        k1, mae1 = at(pooled, 1.0)
        boot = [fit_constant([frames[i] for i in rng.integers(0, len(frames), len(frames))],
                             SOURCE, 1.0, expo) for _ in range(N_BOOT)]
        ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
        rows = []
        for st, fr in sorted(by.items()):
            if len(fr) < MIN_PEOPLE:
                continue
            c = curve(fr, expo)
            b_c, k_c, m_c = min(c, key=lambda t: t[2])
            k_c1, m_c1 = at(c, 1.0)
            rows.append(dict(cohort=st, n=len(fr), best_b=b_c, k_at_best=k_c,
                             k_at_1=k_c1, mae_at_1=m_c1))
        res["by_exposure"][expo] = dict(
            pooled=dict(best_b=b_star, k_at_best=k_star, mae_at_best=mae_star,
                        k_at_1=k1, mae_at_1=mae1, k_at_1_ci=list(ci)),
            by_cohort=rows)

    (config.RESULTS / "inv009_constant_by_cohort.json").write_text(json.dumps(res, indent=1))

    print(f"Sensitivity constant, {len(frames)} people, "
          f"{config.HORIZON_MIN / 60:.0f}-hour window\n")
    for expo in EXPOSURES:
        e = res["by_exposure"][expo]
        p = e["pooled"]
        what = ("the whole dose" if expo == "a_tot"
                else "insulin committed before the window, mostly the basal tail")
        print(f"--- exposure {expo}: {what} ---")
        print(f"  at Walsh's exponent 1.0 : K = {p['k_at_1']:.0f} "
              f"[{p['k_at_1_ci'][0]:.0f}, {p['k_at_1_ci'][1]:.0f}]   "
              f"MAE {p['mae_at_1']:.2f} mg/dL")
        print(f"{'cohort':<12s}{'people':>7s}{'K at b=1':>10s}{'MAE':>8s}")
        for r in e["by_cohort"]:
            print(f"{r['cohort']:<12s}{r['n']:7d}{r['k_at_1']:10.0f}{r['mae_at_1']:8.2f}")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# ---------------------------------------------------------------------------
# The settings-comparable constant.
#
# Both exposures above mix insulin that lowers glucose with insulin that holds it
# level, and a per-person intercept does not separate them: within a person the
# variation in a_pre is mostly basal tail, and the variation in a_tot is mostly
# basal too. Fitting either produces a coefficient roughly half what a correction
# is observed to do, which is how 880 came about.
#
# What a settings file claims is narrower and can be measured directly: the fall
# after a correction, per unit given, with nothing subtracted. Restricting to
# isolated corrections at a raised glucose and a six-hour window, so the response
# is inside the observation, that is a quantity comparable with Walsh's 1700.

MIN_BOLUS_U = 1.0
BG_LO, BG_HI = 150.0, 300.0
MIN_NIGHTS = 3


def per_unit_given(study: str, rng=None):
    """Median fall per unit given on isolated corrections, and that times TDD."""
    from . import data
    rng = np.random.default_rng(20260826) if rng is None else rng
    w = data.load(study)
    if w.empty:
        return None
    q = w[w.quiet_before_h >= config.EVENT_QUIET_H]
    t = q[(q.bolus_first30_u >= MIN_BOLUS_U)
          & np.isclose(q.bolus_in_u, q.bolus_first30_u)
          & (q.bg0 >= BG_LO) & (q.bg0 < BG_HI)].copy()
    if len(t) < 50:
        return None
    t["isf"] = t["drop"] / t.bolus_first30_u
    g = t.groupby("subject_id").agg(isf=("isf", "median"), tdd=("tdd_u", "median"),
                                    n=("isf", "size"))
    g = g[(g.isf > 0) & (g.tdd > 0) & (g.n >= MIN_NIGHTS)]
    if len(g) < 20:
        return None
    k = (g.isf * g.tdd).to_numpy(float)
    b = [np.median(k[rng.integers(0, len(k), len(k))]) for _ in range(4000)]
    return dict(cohort=study, n_people=int(len(g)), n_nights=int(g.n.sum()),
                isf=float(g.isf.median()), k=float(np.median(k)),
                k_lo=float(np.percentile(b, 2.5)), k_hi=float(np.percentile(b, 97.5)))
