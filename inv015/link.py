#!/usr/bin/env python3
"""Which number scales the published equation for one person: the profile ISF or the daily dose?

The glucose shape of Romero-Ugalde et al. (doi:10.2337/db25-920-P) is taken as given
here. Two links turn it into a sensitivity for one person. The daily-dose link sets the
value at one glucose, the anchor, to a constant divided by the person's total daily dose,
and the profile link sets it to the person's entered sensitivity. They agree for anyone whose entered
sensitivity x daily dose equals the constant, and they differ by however far a person's
entered setting sits from that rule. This script asks which of the two numbers a person's
own corrections follow.

The records are the isolated corrections of REPLACE-BG already collected by
`inv009.correction_landing`: a calculator correction at 150 to 350 mg/dL of at least
0.5 U, no carbohydrate for six hours either side and no further bolus for six hours
after. No algorithm doses in that trial, so the dose is the calculator's and the
programmed basal is the only other insulin. Loop and the Control-IQ cohorts are left out
because their algorithms keep dosing after the correction and the table counts the bolus
alone (REVIEW.md item 1).

What is measured, per person: the median over their corrections of the fall in six hours
divided by the units given, divided again by the published shape at the starting glucose,
which puts every correction on the footing of the anchor glucose. It is an outcome ratio. The fall
includes whatever glucose would have done without the correction, so its level is not
read. What is read is how it varies across people:

    log(achieved) = a + b_dose x log(daily dose) + b_profile x log(entered sensitivity)

The daily-dose link predicts b_dose = -1 and b_profile = 0. The profile link predicts
b_dose = 0 and b_profile = 1.

Expected before the run. The dose was sized by the entered sensitivity, so a person with
a high entered value gave fewer units for the same excess glucose, and any fall that did
not come from the correction is then divided by fewer units. That pushes b_profile up and
b_dose towards 0 whichever link is right, so the test leans towards the profile. The
size of the lean is measured first, by putting constructed falls with a known answer
through the same regression on the real people, doses and starting glucose (REVIEW.md
item 4). A b_dose clearly below 0 with the entered value in the model would say the
daily dose carries information about a person that their profile does not. Entered
sensitivity and daily dose are correlated across these adults (slope -0.72 in
docs/preprint-dosing-constants.md), so the two coefficients are expected to carry wide
intervals. No size was predicted for either.

Added after the first run, which returned b_profile near 1. A fall that follows the
excess over target whatever the dose was would return the profile link by construction,
because the calculator sized the dose from that same excess and the entered value. The
overnight windows of INV-009 put the fall on nights with no bolus at 46% of the fall per
unit on correction nights in this cohort, so the constructed falls were extended with
that share and with 75%. In type 1 diabetes that fall is basal and earlier doses acting.

The run, 11 November 2026 (results/inv015_link.md). On 150 adults and 3,433 corrections
b_dose was -0.07 (-0.43 to 0.26) and b_profile 0.95 (0.55 to 1.32). Constructed falls
under the daily-dose link returned -0.54 and 0.37 with 46% of the fall independent of the
dose, outside both intervals, and -0.30 and 0.60 with 75%, inside both. Under the profile
link they returned between -0.03 and 0.06 and between 0.88 and 0.97. The corrections
followed the entered sensitivity, and the daily-dose link is excluded unless about three
quarters of the six-hour fall was independent of the dose. This is one trial of adults
on pumps with no algorithm. It does not test whether sensitivity moves with recent dosing
within a person, which the profile link ignores and the daily-dose link follows.

Changed on 12 November 2026: the anchor in `inv013.candidate` moved from 180 to 140
mg/dL. Expected before the rerun: the published value at 140 is 1.44 times its value at
180, so every achieved value is 1.44 times larger. That is a constant in logs, so both
coefficients, the constructed falls and the misses should be unchanged, and the last
column of the group table should be 1.44 times its earlier values (0.96, 1.22, 1.04 and
1.02). The rerun agreed. A value of that column above 1 says the six-hour fall per unit
was larger than the rule anchored at 140 mg/dL predicts. The level is still not read as a
sensitivity, since the fall includes basal and earlier doses acting.

Modelling choices: five or more corrections per person. The person's entered value is
the median of the values the calculator recorded at those corrections. The shape is read
at the starting glucose, which is where a controller would read it. In the constructed
falls, the flat fall unrelated to the correction is 0 or 30 mg/dL and the noise has a standard
deviation of 40 mg/dL, and neither figure is a measurement. Intervals are 95%, from
resampling people 2,000 times.

Run as `python3 -m inv015.link` to write results/inv015_link.md.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inv013.candidate import ANCHOR_GLUCOSE, RESULTS, curve_shape

COHORT = "ReplaceBG"
MIN_CORRECTIONS = 5
N_BOOT = 2000
TRUTHS = ("daily dose", "profile")
# 0.46 is the fall on nights with no bolus over the fall on correction nights, per unit, in
# REPLACE-BG (8.12 of 17.65 mg/dL per U, results/inv009_measured_isf_by_study.parquet).
UNRELATED = (("none", 0.0, 0.0), ("30 mg/dL at any glucose", 30.0, 0.0),
             ("46% of the fall, following the excess over target", 0.0, 0.46),
             ("75% of the fall, following the excess over target", 0.0, 0.75))


def people(corr: pd.DataFrame, fall: np.ndarray | None = None) -> pd.DataFrame:
    """One row per person: daily dose, entered sensitivity, achieved sensitivity at the anchor."""
    c = corr.assign(fall=corr.bg_start - corr.bg_6h if fall is None else fall)
    c["achieved"] = c.fall / c.given / curve_shape(c.bg_start)
    p = c.groupby("sid").agg(n=("fall", "size"), tdd=("tdd", "first"), isf=("isf", "median"),
                             achieved=("achieved", "median"), start=("bg_start", "median"),
                             land=("bg_6h", "median"), target=("target", "median"))
    return p[(p.n >= MIN_CORRECTIONS) & (p.achieved > 0) & (p.tdd > 0)]


def fit(p: pd.DataFrame) -> np.ndarray:
    """b_dose and b_profile from least squares in logs."""
    X = np.column_stack([np.ones(len(p)), np.log(p.tdd), np.log(p.isf)])
    return np.linalg.lstsq(X, np.log(p.achieved.to_numpy(float)), rcond=None)[0][1:]


def misses(p: pd.DataFrame) -> np.ndarray:
    """Median miss of each link as a factor, its constant set to the cohort median, and their ratio."""
    out = []
    for pred in (1.0 / p.tdd, p.isf):
        r = np.log(p.achieved / pred)
        out.append(float(np.exp(np.median(np.abs(r - np.median(r))))))
    return np.array(out + [out[0] / out[1]])


def boot(p: pd.DataFrame, stat, n: int = N_BOOT) -> list[tuple[float, float, float]]:
    rng = np.random.default_rng(0)
    draws = np.array([stat(p.iloc[rng.integers(0, len(p), len(p))]) for _ in range(n)])
    return [(float(e), float(lo), float(hi))
            for e, lo, hi in zip(stat(p), *np.percentile(draws, [2.5, 97.5], axis=0))]


def constructed(corr: pd.DataFrame, truth: str, f0: float, noise: float, share: float = 0.0,
                seed: int = 0) -> np.ndarray:
    """Falls for the real people and doses under a known link, level matched to the cohort.

    `f0` is a fall unrelated to the correction that is the same at every glucose. `share` is the
    part of the fall that follows the excess over target whatever the dose was.
    """
    base = 1.0 / corr.tdd if truth == "daily dose" else corr.isf
    real = corr.bg_start - corr.bg_6h
    k = np.median(real / corr.given / curve_shape(corr.bg_start) / base)
    c = np.median(real / (corr.bg_start - corr.target))
    rng = np.random.default_rng(seed)
    return ((1 - share) * k * base * curve_shape(corr.bg_start) * corr.given
            + share * c * (corr.bg_start - corr.target) + f0
            + rng.normal(0.0, noise, len(corr))).to_numpy(float)


def main() -> None:
    corr = pd.read_parquet(RESULTS / "inv009_correction_landing.parquet")
    # the landing table holds repeated rows (86 of 3,583 here), REVIEW.md item 8
    corr = corr[corr.study == COHORT].drop_duplicates().copy()
    ent = pd.read_parquet(RESULTS / "inv009_entered_isf.parquet").set_index("subject_id")
    corr["tdd"] = corr.sid.map(ent.tdd_u)
    assert (corr.sid.map(ent.age).dropna() >= 18).all(), "REPLACE-BG enrolled adults"
    p = people(corr)
    corr = corr[corr.sid.isin(p.index)]

    f = lambda v, d=2: f"{v[0]:.{d}f} ({v[1]:.{d}f} to {v[2]:.{d}f})"
    b, m = boot(p, fit), boot(p, misses)
    flat = boot(people(corr.assign(bg_start=ANCHOR_GLUCOSE, bg_6h=ANCHOR_GLUCOSE - (corr.bg_start - corr.bg_6h))), fit)
    product = p.isf * p.tdd
    q = np.percentile(product, [10, 50, 90])
    L = ["# INV-015: profile ISF or daily dose as the scale of the equation", "",
         "Generated by `inv015/link.py`. Isolated calculator corrections in REPLACE-BG, where no",
         "algorithm doses. The achieved sensitivity is the six-hour fall per unit given, put on the",
         f"footing of {ANCHOR_GLUCOSE:.0f} mg/dL by the published shape. It is an outcome ratio and its level is not",
         "read. Intervals are 95% over people.", "",
         f"{len(p)} adults with {MIN_CORRECTIONS} or more corrections, {len(corr)} corrections. Their entered",
         f"sensitivity x daily dose has a median of {q[1]:.0f}, with the tenth and ninetieth percentiles at",
         f"{q[0]:.0f} and {q[2]:.0f}. Entered sensitivity and daily dose correlate at "
         f"{np.corrcoef(np.log(p.isf), np.log(p.tdd))[0, 1]:.2f} in logs.", "",
         "## What the regression returns when the link is known", "",
         "Constructed falls for the same people, doses and starting glucose. The 30 mg/dL, the 75% and",
         "the noise (standard deviation 40 mg/dL) are modelling choices. The 46% is the fall on nights",
         "with no bolus as a share of the fall per unit on correction nights (INV-009, overnight).", "",
         "| True link | Fall unrelated to the dose | People | b_dose | b_profile |", "|---|---|---:|---:|---:|"]
    for truth in TRUTHS:
        for label, f0, share in UNRELATED:
            s = people(corr, constructed(corr, truth, f0, 40.0, share))
            bb = boot(s, fit)
            L.append(f"| {truth} | {label} | {len(s)} | {f(bb[0])} | {f(bb[1])} |")
    L += ["", "## What the corrections returned", "",
          "| | b_dose | b_profile |", "|---|---:|---:|",
          "| Daily-dose link predicts | -1 | 0 |", "| Profile link predicts | 0 | 1 |",
          f"| Measured | {f(b[0])} | {f(b[1])} |",
          f"| Measured, shape not applied | {f(flat[0])} | {f(flat[1])} |", "",
          "## How far each link sits from a person's achieved value", "",
          "Each link is given one free constant, set to the cohort median, so only the ordering of",
          "people is scored. The miss is the median factor between the link and the achieved value.", "",
          "| Link | Median miss, as a factor |", "|---|---:|",
          f"| Constant / daily dose | {f(m[0])} |", f"| Constant x entered sensitivity | {f(m[1])} |",
          f"| Ratio of the two misses | {f(m[2])} |", "",
          "## Where corrections landed, by how far the entered setting sits from the rule", "",
          "People in four groups of equal size by entered sensitivity x daily dose. The profile link",
          "predicts the last column is level across the groups. The daily-dose link predicts it falls",
          "in proportion to 1 / the product.", "",
          "| Entered sensitivity x daily dose | People | Start, mg/dL | Target | At 6 h | Achieved / entered |",
          "|---|---:|---:|---:|---:|---:|"]
    grp = pd.qcut(product, 4)
    for g, d in p.groupby(grp, observed=True):
        r = boot(d, lambda x: np.array([np.median(x.achieved / x.isf)]))[0]
        L.append(f"| {g.left:.0f}–{g.right:.0f} | {len(d)} | {d.start.median():.0f} | {d.target.median():.0f} | "
                 f"{d.land.median():.0f} | {f(r)} |")
    out = RESULTS / "inv015_link.md"
    out.write_text("\n".join(L) + "\n")
    print(out.read_text())


if __name__ == "__main__":
    main()
