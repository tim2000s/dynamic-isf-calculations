"""The three Walsh constants, recomputed across every cohort we hold.

Walsh's rules of thumb set a sensitivity factor at 1700 divided by total daily
dose, a carb ratio at 500 divided by it, and basal at half of it. A May analysis
tested those against 138 open-source loop users and found all three wrong. This
repeats it on roughly five times the people, adds seven cohorts that were not
available then, and reports each separately rather than pooling.

Two sources. The JAEB public archives supply settings people had entered along
with a daily dose measured from their own pump record. The OpenAPS Commons
extract supplies the original 138, unchanged, so the earlier result can be read
against the new ones rather than replaced by them.

One definitional caution runs through the basal figures. Walsh's half refers to
the basal a profile programmes. Under an automated system what is delivered is
not what is programmed, so both are reported where a study records the schedule,
and they are not interchangeable.

    python3 -m inv009.walsh_constants
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

from . import config

WALSH = {"isf_x_tdd": 1700.0, "cr_x_tdd": 500.0, "basal_frac": 0.50, "slope": -1.0}
MIN_N = 15
N_BOOT = 4000


def boot_ci(x, fn=np.median, n_boot=N_BOOT, seed=0, lo=2.5, hi=97.5):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    vals = np.array([fn(x[rng.integers(0, len(x), len(x))]) for _ in range(n_boot)])
    return tuple(float(v) for v in np.percentile(vals, [lo, hi]))


def loglog_slope(tdd, isf, seed=0, n_boot=1500):
    """Slope of log ISF on log TDD, with a bootstrap interval. Walsh implies -1."""
    t = np.asarray(tdd, float)
    v = np.asarray(isf, float)
    m = np.isfinite(t) & np.isfinite(v) & (t > 0) & (v > 0)
    lt, lv = np.log(t[m]), np.log(v[m])
    if len(lt) < 8:
        return np.nan, (np.nan, np.nan), int(len(lt))
    slope = float(np.polyfit(lt, lv, 1)[0])
    rng = np.random.default_rng(seed)
    bs = np.empty(n_boot)
    for i in range(n_boot):
        j = rng.integers(0, len(lt), len(lt))
        bs[i] = np.polyfit(lt[j], lv[j], 1)[0] if np.std(lt[j]) > 1e-9 else np.nan
    return slope, tuple(float(x) for x in np.nanpercentile(bs, [2.5, 97.5])), int(len(lt))


def load() -> pd.DataFrame:
    frames = []

    j = pd.read_parquet(config.RESULTS / "inv009_entered_isf.parquet")
    j = j[(j.isf > 0) & (j.tdd_u > 0) & (j.n_days >= config.MIN_DAYS_TDD)]
    frames.append(pd.DataFrame(dict(
        user_id=j.subject_id, cohort=j.study, source="JAEB",
        isf=j.isf, cr=j.cr, tdd=j.tdd_u,
        basal_delivered=j.tdd_basal_u, basal_profile=np.nan,
        dynisf=False, n_days=j.n_days)))

    # The May cohort file lives in the older working directory, so look in both.
    candidates = [config.ROOT / "canonical_cohort.json",
                  pathlib.Path.home() / "oref-investigations" / "canonical_cohort.json"]
    p = next((x for x in candidates if x.exists()), candidates[0])
    if p.exists():
        c = pd.DataFrame(json.loads(p.read_text()))
        c = c[c.get("in_cohort", True).astype(bool)] if "in_cohort" in c else c
        c = c[(c.isf > 0) & (c.tdd > 0)]
        frames.append(pd.DataFrame(dict(
            user_id=c.user_id, cohort=c.cohort, source="OpenAPS Commons",
            isf=c.isf, cr=c.cr, tdd=c.tdd,
            basal_delivered=np.nan, basal_profile=c.basal,
            dynisf=(c.dynisf_frac.fillna(0) > 0.1) if "dynisf_frac" in c else False,
            n_days=c.n_days)))

    d = pd.concat(frames, ignore_index=True)
    d["isf_x_tdd"] = d.isf * d.tdd
    d["cr_x_tdd"] = d.cr * d.tdd
    # Delivered where that is what a study records, programmed where it is.
    d["basal_frac"] = np.where(d.basal_profile.notna(),
                               d.basal_profile / d.tdd,
                               d.basal_delivered / d.tdd)
    d["basal_basis"] = np.where(d.basal_profile.notna(), "programmed", "delivered")
    d.loc[~d.basal_frac.between(0.05, 0.95), "basal_frac"] = np.nan
    return d


def summarise(d: pd.DataFrame, label: str, basis: str) -> dict:
    row = dict(cohort=label, n=int(len(d)), basal_basis=basis)
    for key in ("isf_x_tdd", "cr_x_tdd", "basal_frac"):
        v = d[key].dropna()
        row[f"{key}_n"] = int(len(v))
        if len(v) < 5:
            row[f"{key}_median"] = np.nan
            row[f"{key}_ci"] = [np.nan, np.nan]
            row[f"{key}_excludes_walsh"] = None
            continue
        med = float(v.median())
        ci = boot_ci(v)
        row[f"{key}_median"] = med
        row[f"{key}_ci"] = list(ci)
        row[f"{key}_excludes_walsh"] = bool(WALSH[key] < ci[0] or WALSH[key] > ci[1])
    s, ci, n = loglog_slope(d.tdd, d.isf)
    row.update(slope=s, slope_ci=list(ci), slope_n=n,
               slope_excludes_minus_one=bool(np.isfinite(ci[0]) and (-1.0 < ci[0] or -1.0 > ci[1])))
    row["median_tdd"] = float(d.tdd.median())
    row["median_isf"] = float(d.isf.median())
    return row


def main() -> int:
    config.ensure_dirs()
    d = load()
    d.to_parquet(config.RESULTS / "inv009_walsh_cohort.parquet", index=False)
    print(f"{len(d)} people: {int((d.source=='JAEB').sum())} from the JAEB archives, "
          f"{int((d.source=='OpenAPS Commons').sum())} from OpenAPS Commons\n")

    rows = []
    for label, sub in d.groupby("cohort"):
        if len(sub) < MIN_N:
            continue
        rows.append(summarise(sub, label, sub.basal_basis.mode().iloc[0]))
    for label, sub in d.groupby("source"):
        rows.append(summarise(sub, f"ALL {label}", "mixed"))
    rows.append(summarise(d, "EVERYTHING", "mixed"))

    # Age is not a nuisance here, it is the finding. Walsh's constants were set
    # for adults, and pooling children with adults moves both the constant and
    # the slope. The JAEB archives carry age; OpenAPS Commons does not, but that
    # cohort is adults.
    j = pd.read_parquet(config.RESULTS / "inv009_entered_isf.parquet")
    j = j[(j.isf > 0) & (j.tdd_u > 0) & (j.n_days >= config.MIN_DAYS_TDD)]
    j = j.rename(columns={"subject_id": "user_id", "tdd_u": "tdd",
                          "tdd_basal_u": "basal_delivered"})
    j["basal_profile"] = np.nan
    j["isf_x_tdd"] = j.isf * j.tdd
    j["cr_x_tdd"] = j.cr * j.tdd
    j["basal_frac"] = j.basal_delivered / j.tdd
    j.loc[~j.basal_frac.between(0.05, 0.95), "basal_frac"] = np.nan
    age_rows = []
    for band, sub in j.groupby("age_band"):
        if len(sub) < MIN_N or not band:
            continue
        age_rows.append(summarise(sub, f"JAEB {band}", "delivered"))
    adults = j[j.age >= 18]
    age_rows.append(summarise(adults, "JAEB adults (18+)", "delivered"))
    kids = j[j.age < 18]
    age_rows.append(summarise(kids, "JAEB under 18", "delivered"))
    rows.extend(age_rows)
    R = pd.DataFrame(rows)
    R.to_parquet(config.RESULTS / "inv009_walsh_by_cohort.parquet", index=False)

    def fmt(r, key, dp=0):
        if not np.isfinite(r[f"{key}_median"]):
            return f"{'-':>22s}"
        lo, hi = r[f"{key}_ci"]
        star = "*" if r[f"{key}_excludes_walsh"] else " "
        return f"{r[f'{key}_median']:8.{dp}f} [{lo:.{dp}f}, {hi:.{dp}f}]{star}"

    print("ISF x TDD   (Walsh 1700)        CR x TDD   (Walsh 500)      basal share (Walsh 0.50)")
    print(f"{'cohort':>20s} {'n':>5s} {'ISF x TDD':>23s} {'CR x TDD':>21s} {'basal':>20s}")
    for _, r in R.iterrows():
        print(f"{r['cohort']:>20s} {r['n']:5d} {fmt(r,'isf_x_tdd')} "
              f"{fmt(r,'cr_x_tdd')} {fmt(r,'basal_frac',2)}")
    print("\n  * the 95% interval excludes the Walsh value\n")

    print("The parametric form: slope of log ISF on log TDD, which Walsh puts at -1")
    print(f"{'cohort':>20s} {'n':>5s} {'slope':>8s} {'95% interval':>20s} {'excludes -1':>12s}")
    for _, r in R.iterrows():
        if not np.isfinite(r["slope"]):
            continue
        lo, hi = r["slope_ci"]
        print(f"{r['cohort']:>20s} {r['slope_n']:5d} {r['slope']:+8.2f} "
              f"[{lo:+.2f}, {hi:+.2f}]{'':>4s} {'yes' if r['slope_excludes_minus_one'] else 'no':>12s}")

    (config.RESULTS / "inv009_walsh_constants.json").write_text(
        json.dumps(rows, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
