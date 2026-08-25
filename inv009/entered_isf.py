"""What people had entered, against how much insulin they use.

This is the cleanest question in the investigation and the narrowest. It asks
whether the sensitivity factors people were actually running scale with total
daily dose the way v1 assumes, the way v2 assumes, or some other way. It says
nothing directly about physiology: an entered setting is a decision, informed by
a clinic, a rule of thumb and a person's own experience of what works.

That is worth stating plainly rather than hedging around, because the 1800 rule
is itself the source of v1's exponent. If entered settings sit on 1800/TDD, that
is partly people following the rule the equation was built from, and the finding
is about the rule's grip rather than about the body. The value of the comparison
is that v2 proposes a squared law, and nothing in practice looks like one.

    python3 -m inv009.entered_isf
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from . import config, data, db, stats

# The constants each equation implies at a normal target, from inv008.dynisf:
# v1 gives 1800/(TDD*ln(NT/75+1)) and v2 gives 2300/(ln(NT/75)*TDD^2*0.02).
V1_CONSTANT = 1800.0 / np.log(config_nt := 99.0 / 75.0 + 1.0)
V2_CONSTANT = 2300.0 / (np.log(99.0 / 75.0) * 0.02)


def subject_table() -> pd.DataFrame:
    """One row per subject: entered settings, daily dose, age, study."""
    ent = db.entered_isf()
    rows = []
    for f in sorted(config.WINDOW_CACHE.glob("*.parquet")):
        d = pd.read_parquet(f, columns=["subject_id", "study", "tdd_u", "tdd_basal_u",
                                        "basal_frac", "n_days", "age"])
        if d.empty:
            continue
        rows.append(d.iloc[[0]])
    tdd = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    tdd = tdd.drop(columns=["study"])          # already on the settings side
    t = ent.merge(tdd, on="subject_id", how="inner")
    t["age_band"] = stats.band_of(t.age, config.AGE_BANDS)
    t["isf_x_tdd"] = t.isf * t.tdd_u
    t["isf_x_tdd2"] = t.isf * t.tdd_u ** 2
    t["cr_x_tdd"] = t.cr * t.tdd_u
    return t


def analyse(t: pd.DataFrame) -> dict:
    out: dict = {}
    use = t.dropna(subset=["isf", "tdd_u"])
    use = use[(use.tdd_u > 0) & (use.isf > 0) & (use.n_days >= config.MIN_DAYS_TDD)]

    def block(d: pd.DataFrame, label: str) -> dict:
        r = stats.loglog(d.tdd_u, d.isf)
        ci = stats.boot_slope_ci(d.tdd_u, d.isf)
        return dict(label=label, n=r["n"], slope=r["slope"], se=r["se"], r2=r["r2"],
                    ci_lo=ci[0], ci_hi=ci[1],
                    deming=stats.deming(d.tdd_u, d.isf),
                    isf_x_tdd_median=float(d.isf_x_tdd.median()),
                    isf_x_tdd_iqr=[float(d.isf_x_tdd.quantile(.25)),
                                   float(d.isf_x_tdd.quantile(.75))],
                    rho_rule_vs_tdd=float(d[["isf_x_tdd", "tdd_u"]].corr(method="spearman")
                                          .iloc[0, 1]),
                    rho_v2_vs_tdd=float(d[["isf_x_tdd2", "tdd_u"]].corr(method="spearman")
                                        .iloc[0, 1]),
                    summary=stats.describe_slope(r, ci))

    out["overall"] = block(use, "all subjects")
    out["by_study"] = [block(d, s) for s, d in use.groupby("study") if len(d) >= 20]
    out["by_age_band"] = [block(d, b) for b, d in use.groupby("age_band") if len(d) >= 20]
    out["by_source"] = [block(d, s) for s, d in use.groupby("source") if len(d) >= 20]

    # Basal-only dose cannot itself have been chosen by the correction factor, so
    # it is the version of the x axis that a sensitivity factor cannot feed back into.
    b = use.dropna(subset=["tdd_basal_u"])
    b = b[b.tdd_basal_u > 0]
    rb = stats.loglog(b.tdd_basal_u, b.isf)
    out["vs_basal_dose"] = dict(n=rb["n"], slope=rb["slope"], se=rb["se"],
                                ci=list(stats.boot_slope_ci(b.tdd_basal_u, b.isf)),
                                note="basal dose as the regressor: a correction factor "
                                     "cannot influence it the way it influences total dose")

    cr = use.dropna(subset=["cr"])
    cr = cr[cr.cr > 0]
    rc = stats.loglog(cr.tdd_u, cr.cr)
    out["carb_ratio"] = dict(n=rc["n"], slope=rc["slope"], se=rc["se"],
                             ci=list(stats.boot_slope_ci(cr.tdd_u, cr.cr)),
                             cr_x_tdd_median=float(cr.cr_x_tdd.median()),
                             summary=stats.describe_slope(rc, stats.boot_slope_ci(cr.tdd_u, cr.cr)))

    out["constants"] = dict(
        v1_implied_isf_x_tdd=float(V1_CONSTANT),
        v2_implied_isf_x_tdd2=float(V2_CONSTANT),
        observed_isf_x_tdd_median=float(use.isf_x_tdd.median()),
        observed_isf_x_tdd2_median=float(use.isf_x_tdd2.median()))
    return out


def main() -> int:
    config.ensure_dirs()
    t = subject_table()
    res = analyse(t)
    t.to_parquet(config.RESULTS / "inv009_entered_isf.parquet", index=False)
    (config.RESULTS / "inv009_entered_isf.json").write_text(json.dumps(res, indent=2))

    o = res["overall"]
    print(f"Entered sensitivity against daily dose, {o['n']} people")
    print(f"  log-log slope {o['summary']}")
    print(f"  Deming slope  {o['deming']:+.3f}   (v1 implies -1, v2 implies -2)")
    print(f"  ISF x TDD median {o['isf_x_tdd_median']:.0f} "
          f"(IQR {o['isf_x_tdd_iqr'][0]:.0f}-{o['isf_x_tdd_iqr'][1]:.0f}); "
          f"v1 implies {V1_CONSTANT:.0f}")
    print(f"  Spearman(ISF x TDD, TDD) = {o['rho_rule_vs_tdd']:+.3f}  "
          f"(zero would mean the 1800 rule holds)")
    print(f"  Spearman(ISF x TDD^2, TDD) = {o['rho_v2_vs_tdd']:+.3f}  "
          f"(zero would mean v2's squared law holds)")
    print("\nBy study")
    for b in res["by_study"]:
        print(f"  {b['label']:12s} n={b['n']:4d}  slope {b['slope']:+.3f} "
              f"[{b['ci_lo']:+.3f},{b['ci_hi']:+.3f}]  Deming {b['deming']:+.3f}")
    print("\nBy age")
    for b in res["by_age_band"]:
        print(f"  {b['label']:8s} n={b['n']:4d}  slope {b['slope']:+.3f} "
              f"[{b['ci_lo']:+.3f},{b['ci_hi']:+.3f}]")
    print(f"\nBasal dose only: slope {res['vs_basal_dose']['slope']:+.3f} "
          f"[{res['vs_basal_dose']['ci'][0]:+.3f},{res['vs_basal_dose']['ci'][1]:+.3f}]")
    print(f"Carb ratio: {res['carb_ratio']['summary']}")
    print(f"  CR x TDD median {res['carb_ratio']['cr_x_tdd_median']:.0f} (the 500 rule)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
