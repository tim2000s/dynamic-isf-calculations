"""The constant as the data gives it, rather than as people entered it.

Every figure in the Walsh work so far multiplies a setting somebody typed into a
pump by a daily dose. That answers what people run, which is a decision. It does
not answer what the constant should be.

This measures the sensitivity from each person's own glucose response and
multiplies that by their own measured daily dose. It needs no settings file, so
it covers DCLP3 and IOBP2 as well, and those two carry 446 people between them
that every earlier table omitted.

THE LEVEL IS NOT RECOVERABLE FROM THIS DATA, and the run below shows why rather
than working around it.

Fitting the overnight fall against insulin already committed gives a slope that
is attenuated: pulled toward zero because insulin action is reconstructed rather
than observed, because residual carbohydrate rides along with meal boluses, and
because a loop delivers most insulin exactly when glucose is refusing to move.
Measured against people who also have an entered sensitivity, the factor is 6.78.

Calibrating on that factor would make the answer depend on entered settings,
which is the thing this module exists to avoid. Worse, the matched-correction
estimator, which needs no basal model and no settings, disagrees: it puts the
factor between 2.01 in the open-loop cohort and 5.31 in the youngest one. Two
estimators differing by nearly two-fold on the level, and an attenuation that
varies by cohort, mean neither the absolute constant nor the comparison between
cohorts can be trusted from these data.

What survives is structure that a common multiplier cannot change, principally
whether the constant drifts with daily dose. That question is answered here. The
constants printed below carry the calibration and are reported for completeness,
not as an answer.

    python3 -m inv009.measured_constant
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from . import config, data, db, stats
from .effective_isf import _fit_subject

MIN_WINDOWS = 60
N_BOOT = 3000
WALSH = 1700.0


def per_subject() -> pd.DataFrame:
    rows = []
    for study in config.COHORTS:
        w = data.load(study)
        if w.empty:
            continue
        for sid, d in w.groupby("subject_id", sort=False):
            if len(d) < MIN_WINDOWS:
                continue
            f = _fit_subject(d, "a_pre")
            if not f:
                continue
            rows.append(dict(subject_id=sid, study=study,
                             isf_measured=f["s"], se=f["se"], n_windows=f["n"],
                             tdd=float(d.tdd_u.iloc[0]),
                             basal_frac=float(d.basal_frac.iloc[0]),
                             age=float(d.age.iloc[0]) if pd.notna(d.age.iloc[0]) else np.nan))
    R = pd.DataFrame(rows)
    R["age_band"] = stats.band_of(R.age, config.AGE_BANDS)
    R["k_raw"] = R.isf_measured * R.tdd
    return R


def boot_median_ci(x, n_boot=N_BOOT, seed=0):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    v = np.array([np.median(x[rng.integers(0, len(x), len(x))]) for _ in range(n_boot)])
    return tuple(float(a) for a in np.percentile(v, [2.5, 97.5]))


def main() -> int:
    config.ensure_dirs()
    R = per_subject()
    print(f"{len(R)} people with a measurable sensitivity, across "
          f"{R.study.nunique()} cohorts")
    pos = R[R.isf_measured > 0]
    print(f"  {len(pos)} of them positive, which is the group a constant can be "
          f"formed from\n")

    # The calibration. Where somebody has both a measured and an entered
    # sensitivity, the ratio between them is what this method costs.
    ent = db.entered_isf().set_index("subject_id").isf
    pos = pos.assign(isf_entered=pos.subject_id.map(ent))
    both = pos.dropna(subset=["isf_entered"])
    both = both[both.isf_entered > 0]
    ratio = (both.isf_entered / both.isf_measured).replace([np.inf, -np.inf], np.nan).dropna()
    cal = float(np.median(ratio))
    cal_ci = boot_median_ci(ratio)
    print(f"CALIBRATION  {len(both)} people have both a measured and an entered value")
    print(f"  entered divided by measured: median {cal:.2f} "
          f"(95% interval {cal_ci[0]:.2f} to {cal_ci[1]:.2f})")
    print(f"  so this method reads about {1/cal:.0%} of a person's working sensitivity,")
    print(f"  and the raw product below is multiplied by {cal:.2f} to give the constant\n")

    R["k_calibrated"] = R.k_raw * cal
    pos = pos.assign(k_calibrated=pos.k_raw * cal)
    R.to_parquet(config.RESULTS / "inv009_measured_constant.parquet", index=False)

    def block(d, label):
        k = d.k_calibrated.dropna()
        if len(k) < 12:
            return None
        lo, hi = boot_median_ci(k)
        return dict(label=label, n=int(len(k)), median=float(k.median()),
                    ci_lo=lo, ci_hi=hi,
                    covers_walsh=bool(lo <= WALSH <= hi),
                    median_tdd=float(d.tdd.median()),
                    median_isf=float(d.isf_measured.median() * cal))

    out = {"calibration": cal, "calibration_ci": list(cal_ci),
           "n_calibration": int(len(both)), "by_cohort": [], "by_age": [], "overall": None}

    print("MEASURED CONSTANT by cohort, calibrated, against Walsh's 1700")
    print(f"  {'cohort':>28s} {'n':>5s} {'constant':>9s} {'95% interval':>18s} "
          f"{'covers 1700':>12s} {'median TDD':>11s}")
    for study, d in pos.groupby("study"):
        b = block(d, config.COHORTS[study]["label"])
        if b:
            out["by_cohort"].append(b)
            print(f"  {b['label']:>28s} {b['n']:5d} {b['median']:9.0f} "
                  f"[{b['ci_lo']:6.0f}, {b['ci_hi']:6.0f}] {'yes' if b['covers_walsh'] else 'no':>12s} "
                  f"{b['median_tdd']:11.1f}")
    b = block(pos, "EVERYONE")
    out["overall"] = b
    print(f"  {b['label']:>28s} {b['n']:5d} {b['median']:9.0f} "
          f"[{b['ci_lo']:6.0f}, {b['ci_hi']:6.0f}] {'yes' if b['covers_walsh'] else 'no':>12s} "
          f"{b['median_tdd']:11.1f}")

    print("\nby age")
    print(f"  {'band':>28s} {'n':>5s} {'constant':>9s} {'95% interval':>18s} {'covers 1700':>12s}")
    for band, d in pos.groupby("age_band"):
        if not band:
            continue
        b = block(d, band)
        if b:
            out["by_age"].append(b)
            print(f"  {b['label']:>28s} {b['n']:5d} {b['median']:9.0f} "
                  f"[{b['ci_lo']:6.0f}, {b['ci_hi']:6.0f}] "
                  f"{'yes' if b['covers_walsh'] else 'no':>12s}")

    # Does the measured constant drift with daily dose? A true constant should not.
    r = stats.loglog(pos.tdd, pos.k_calibrated)
    ci = stats.boot_slope_ci(pos.tdd, pos.k_calibrated)
    out["constant_vs_tdd_slope"] = dict(slope=r["slope"], ci=list(ci), n=r["n"])
    print(f"\nDoes the measured constant hold across doses? A real constant would not drift.")
    print(f"  slope of the constant against daily dose: {r['slope']:+.3f} "
          f"[{ci[0]:+.3f}, {ci[1]:+.3f}]")
    print(f"  a slope of 0 means one constant fits every dose; "
          f"{'it does not' if not (ci[0] <= 0 <= ci[1]) else 'consistent with 0'}")

    (config.RESULTS / "inv009_measured_constant.json").write_text(
        json.dumps(out, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
