#!/usr/bin/env python3
"""How the published curve attaches to a person's profile ISF: multiplied or shifted?

The glucose curve of Romero-Ugalde et al. (doi:10.2337/db25-920-P), f, is taken as given.
INV-015 found that a person's corrections follow their entered sensitivity, so the curve is
tied to that value at the anchor, 140 mg/dL. Two transformations do that. Multiplying gives
g(x) = f(x) x profile / f(140), which stretches the curve and keeps its proportions.
Shifting gives g(x) = f(x) + profile - f(140), which slides it and keeps its absolute
drop. The published curve is 50.9 mg/dL per U at 140 mg/dL, so the two agree for a
profile near 50 and separate on either side: below it the shift falls faster with glucose
and reaches zero for a profile under 35.6, above it the shift is flatter.

The records are the isolated calculator corrections of REPLACE-BG used by INV-015
(`inv009.correction_landing`): no algorithm doses, no carbohydrate for six hours either
side, no further bolus for six hours after. The fall in six hours is set against what each
form predicts, sensitivity at the starting glucose x units given. The fall includes basal
and earlier doses acting, so its level is not read and every statistic is a ratio in which
a level cancels.

Two statistics, fixed before the run.
  Across people. People are put in four groups by entered sensitivity. In each, the total
  fall is divided by the total a form predicts. The right form gives the same value in
  every group, so the value in the lowest group over the value in the highest is 1.
  Within glucose. The fall over the profile's own prediction is taken above and below 230
  mg/dL (the median starting glucose is 227), and the upper value over the lower is the
  steepness. Multiplying predicts the same steepness whatever the entered sensitivity.
  Shifting predicts a steeper fall with glucose in people with a low entered value, so the
  steepness in the lower half of people over the steepness in the upper half is below 1.

Known truth first (REVIEW.md item 4). Falls are constructed for the same people, doses and
starting glucose under each form, and both statistics are run on them. The calculator sized
the dose from the excess over target and the entered value, so any part of the fall that
follows the excess whatever the dose was is proportional to profile x units. That part
looks like the profile alone and, across people, like the multiplied form, so the first
statistic leans towards multiplying. The constructions carry 0%, 46% and 75% of the fall
in that part, the shares used in INV-015, to measure the lean. The second statistic has
no such lean, since that part has no steepness in either half, but it is diluted by it.

Expected before the run. The predictors alone show that the shift gives a sensitivity of
zero or less at 299 of the 3,433 corrections, in 32 of the 150 people, and INV-015 returned
an exponent of 0.95 (0.55 to 1.32) on the entered value with the multiplied shape applied.
So the multiplied form's ratio across people should be near 1, the shifted form's should
be above 1 (people with a low entered value falling further than the shift predicts), and
the steepness ratio should be nearer 1 than the value the shift constructs. No size was
predicted. A verdict needs the measured values to sit with the constructions of one form
and outside those of the other at 46%.

The run, 12 November 2026 (results/inv017_shift_or_scale.md). Across people, the
multiplied form's ratio was 1.53 (0.81 to 2.89) and the shifted form's 6.88 (3.05 to
17.78). Falls constructed under multiplying gave 1.00, 1.16 and 1.29 for the first and
4.49, 5.22 and 5.81 for the second at 0%, 46% and 75%, all inside the measured intervals.
Falls constructed under shifting gave 0.22, 0.56 and 0.92, and 0.99, 2.52 and 4.15. The
first two of each lie outside the measured intervals and the third inside. The corrections
therefore sit with multiplying, and shifting is excluded unless about three quarters of
the six-hour fall was independent of the dose, the same margin INV-015 left the daily-dose
link. In the 44 people with an entered value of 16 to 35 mg/dL per U the fall was 20.6
mg/dL per unit, where multiplying predicts 14.8 and shifting 4.1. The ratio above 1 for
both forms comes from the highest group, whose fall per unit was lower than the group
below it under every form.

The within-glucose statistic decided nothing. The measured ratio was 0.75 (0.40 to 1.23),
which holds every constructed value from 0.52 to 1.07. The measured steepness itself was
above 1 in both halves (1.41 and 1.88), which no construction reproduced, so the six-hour
fall rises with starting glucose faster than the constructions allow and that statistic is
not read further.

Modelling choices: five or more corrections per person, as in INV-015. A person's group
is set by the median of the values the calculator recorded, and each correction is
predicted from its own recorded value. The split at 230 mg/dL, the noise in the
constructed falls (standard deviation 40 mg/dL) and the 75% share are choices and not
measurements. Intervals are 95%, from resampling people 2,000 times.

Run as `python3 -m inv017.shift_or_scale` to write results/inv017_shift_or_scale.md.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inv008.dynisf import quartic_isf
from inv013.candidate import ANCHOR_GLUCOSE, RESULTS, curve_shape
from inv015.link import COHORT, MIN_CORRECTIONS, N_BOOT

FORMS = ("Profile alone", "Multiplied", "Shifted")
SPLIT = 230.0
SHARES = (0.0, 0.46, 0.75)
STATS = ("Profile alone, lowest over highest group", "Multiplied, lowest over highest group",
         "Shifted, lowest over highest group", "Steepness, lower half of people",
         "Steepness, upper half of people", "Steepness, lower half over upper half")


def sensitivities(isf, bg) -> dict:
    """Sensitivity at a glucose under each form, all equal to the entered value at the anchor."""
    shape = curve_shape(bg)
    return dict(zip(FORMS, (isf + 0.0 * shape, isf * shape,
                            isf + (shape - 1.0) * float(quartic_isf(ANCHOR_GLUCOSE)))))


def person_sums(corr: pd.DataFrame, fall=None) -> pd.DataFrame:
    """One row per person: entered value, and the sums every statistic is a ratio of."""
    c = corr.assign(fall=corr.bg_start - corr.bg_6h if fall is None else fall)
    for form, s in sensitivities(c.isf, c.bg_start).items():
        c[form] = s * c.given
    hi = c.bg_start >= SPLIT
    for name, keep in (("hi", hi), ("lo", ~hi)):
        c["fall_" + name], c["flat_" + name] = c.fall.where(keep, 0.0), c[FORMS[0]].where(keep, 0.0)
    agg = {k: (k, "sum") for k in ("fall", "given", *FORMS, "fall_hi", "flat_hi", "fall_lo", "flat_lo")}
    p = c.groupby("sid").agg(n=("fall", "size"), isf=("isf", "median"), start=("bg_start", "median"), **agg)
    return p[p.n >= MIN_CORRECTIONS]


def stats(p: pd.DataFrame, edges) -> np.ndarray:
    """The six ratios of STATS. `edges` are the quartiles of entered sensitivity, fixed outside."""
    low, high = p[p.isf <= edges[0]], p[p.isf > edges[2]]
    across = [(low.fall.sum() / low[f].sum()) / (high.fall.sum() / high[f].sum()) for f in FORMS]
    steep = [(h.fall_hi.sum() / h.flat_hi.sum()) / (h.fall_lo.sum() / h.flat_lo.sum())
             for h in (p[p.isf <= edges[1]], p[p.isf > edges[1]])]
    return np.array(across + steep + [steep[0] / steep[1]])


def boot(p: pd.DataFrame, edges, n: int = N_BOOT) -> list[tuple[float, float, float]]:
    rng = np.random.default_rng(0)
    draws = np.array([stats(p.iloc[rng.integers(0, len(p), len(p))], edges) for _ in range(n)])
    return [(float(e), float(lo), float(hi))
            for e, lo, hi in zip(stats(p, edges), *np.percentile(draws, [2.5, 97.5], axis=0))]


def constructed(corr: pd.DataFrame, form: str, share: float, noise: float = 40.0, seed: int = 0) -> np.ndarray:
    """Falls under a known form, level matched to the cohort. `share` follows the excess over target."""
    pred = sensitivities(corr.isf, corr.bg_start)[form] * corr.given
    real = corr.bg_start - corr.bg_6h
    k = np.median(real) / np.median(pred)
    c = np.median(real / (corr.bg_start - corr.target))
    rng = np.random.default_rng(seed)
    return ((1 - share) * k * pred + share * c * (corr.bg_start - corr.target)
            + rng.normal(0.0, noise, len(corr))).to_numpy(float)


def main() -> None:
    corr = pd.read_parquet(RESULTS / "inv009_correction_landing.parquet")
    # the landing table holds repeated rows, REVIEW.md item 8
    corr = corr[corr.study == COHORT].drop_duplicates().copy()
    p = person_sums(corr)
    corr = corr[corr.sid.isin(p.index)]
    edges = np.percentile(p.isf, [25, 50, 75])
    shifted = sensitivities(corr.isf, corr.bg_start)["Shifted"]

    f = lambda v: f"{v[0]:.2f} ({v[1]:.2f} to {v[2]:.2f})"
    L = ["# INV-017: the published curve multiplied or shifted onto the profile ISF", "",
         "Generated by `inv017/shift_or_scale.py`, whose docstring holds the modelling choices and",
         "what was expected before the run. Isolated calculator corrections in REPLACE-BG, where no",
         "algorithm doses. The six-hour fall includes basal and earlier doses acting, so its level is",
         "not read. Intervals are 95% over people.", "",
         f"{len(p)} adults with {MIN_CORRECTIONS} or more corrections, {len(corr)} corrections. The shifted form gives a",
         f"sensitivity of zero or less at {(shifted <= 0).sum()} of them, in {corr.sid[shifted <= 0].nunique()} people.", "",
         "## Fall per unit by entered sensitivity", "",
         "Total fall over total units in each group, beside what each form predicts. A form that",
         "attaches the curve correctly has the same achieved-over-predicted value in every group.", "",
         "| Entered sensitivity, mg/dL per U | People | Median start, mg/dL | Fall per unit | " +
         " | ".join(f"{x}, predicted" for x in FORMS) + " | " + " | ".join(f"{x}, achieved / predicted" for x in FORMS) + " |",
         "|---|---:|---:|---:|" + "---:|" * (2 * len(FORMS))]
    bounds = [-np.inf, *edges, np.inf]
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        g = p[(p.isf > lo) & (p.isf <= hi)]
        L.append(f"| {g.isf.min():.0f} to {g.isf.max():.0f} | {len(g)} | {g.start.median():.0f} | {g.fall.sum() / g.given.sum():.1f} | "
                 + " | ".join(f"{g[x].sum() / g.given.sum():.1f}" for x in FORMS) + " | "
                 + " | ".join(f"{g.fall.sum() / g[x].sum():.2f}" for x in FORMS) + " |")
    L += ["", "## The two statistics, measured and under a known form", "",
          "Multiplying predicts 1 in the second and last columns. Shifting predicts 1 in the third.",
          "The constructed rows use the same people, doses and starting glucose, with the stated share",
          "of the fall following the excess over target whatever the dose was.", "",
          "| Falls | " + " | ".join(STATS) + " |", "|---|" + "---:|" * len(STATS),
          "| Measured | " + " | ".join(f(v) for v in boot(p, edges)) + " |"]
    for form in FORMS[1:]:
        for share in SHARES:
            b = boot(person_sums(corr, constructed(corr, form, share)), edges)
            L.append(f"| Constructed, {form.lower()}, {100 * share:.0f}% | " + " | ".join(f(v) for v in b) + " |")
    out = RESULTS / "inv017_shift_or_scale.md"
    out.write_text("\n".join(L) + "\n")
    print(out.read_text())


if __name__ == "__main__":
    main()
