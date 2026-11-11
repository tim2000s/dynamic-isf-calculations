#!/usr/bin/env python3
"""Reference forms for docs/dynamic-isf-8-proposal.md.

The candidate is the glucose curve published by Romero-Ugalde et al. (Diabetes 2025;
74 Suppl 1: 920-P, doi:10.2337/db25-920-P), linked to total daily dose. The published
curve is a sensitivity in mg/dL per unit with no daily dose in it: 79 at 100 mg/dL,
35 at 180, below 20 above 300. The poster's figure fits it in two pieces, a quartic
for glucose in (100, 380] (`inv008.dynisf.quartic_isf`) and a quadratic for [0, 100]
(`_QUADRATIC` below). Both sets of coefficients were checked against the figure on
11 November 2026.

The link to daily dose is this repository's and is a modelling choice. The entered
settings' product for adults (sensitivity x daily dose = 1915, document 6) is taken to
hold at 180 mg/dL, the reference glucose the abstract uses, which gives

    sensitivity = 1915 x quartic(glucose) / quartic(180) / TDD

`curve_profile` is the alternative link, the same shape scaled by the profile value.
Below 100 mg/dL the curve is held at its value at 100 where the published quadratic
falls to 65 near 60 mg/dL, so it gives less insulin there than the published curve;
the March 2026 work declined the quadratic for dosing on the same ground. Above 380
mg/dL it is held at its value at 380, where the published fit ends. Both are
modelling choices. Nothing here was fitted.

The comparators are written from the stock sources read on 10 November 2026:
AndroidAPS upstream master f77abd9b82 (OpenAPSSMBPlugin.calculateRawDynIsf) and
Trio upstream main 0b29cd534 (DynamicISF.calculate). Both were read again on 11
November 2026 at the dev heads (AndroidAPS 0c593a4b09, Trio 902be75e8) and neither
calculation had changed. Trio's logarithmic mode is the AndroidAPS equation with the
profile value cancelling; its defaults (factor 0.8, limits 0.7 and 1.2) are what
separate the two columns. `inv008.dynisf.isf_v1` is
the Boost form with a glucose cap and a velocity term, so the stock equation is
restated here and tested against it with both switched off.

Expected before the run: the curve sits above the stock equation at every glucose from
70 to 300 mg/dL (less insulin), by most near target, and falls more steeply than it.
The reference dose implied by 1915 at 180 mg/dL is expected near 54 U, about twice the
24 U implied by the 1800 rule at 105 mg/dL used in the March 2026 work, and below the
dose implied by this repository's own peak-to-trough levels. The run agreed.

Added 11 November 2026: the measured level in each glucose bin of INV-012 divided by
the published value there, which is the daily dose at which the published curve would
return that measurement. A row that stays level means the published shape fits the
cohort. Expected from the within-person ratios already in results/inv012_peak_trough.md:
level from 100 to 200 mg/dL where no algorithm doses (2.55 against the published 2.61),
a dip at 150 mg/dL under Loop, Trio and oref0 (their 100:150 ratios are 2.6 to 3.1
against the published 1.76), and rising above 200 mg/dL in every cohort that reaches
it (250:200 of 0.84 and 0.90 against 0.80). The run agreed.

`candidate` is the profile-anchored form of this document's first draft (commit
fd3662e), kept as the flat arm of the evaluation.

Run as `python3 -m inv013.candidate` to write results/inv013_comparator.md.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from inv008.dynisf import quartic_isf

DIVISORS = {"rapid, 75 min peak": 55.0, "ultra-rapid, 55 min peak": 65.0, "45 min peak": 75.0}
BG_POINTS = (70.0, 80.0, 100.0, 120.0, 150.0, 180.0, 200.0, 250.0, 300.0, 380.0)
# The worked person: the simulator's, close to the adult entered median of 1915.
PROFILE_ISF, TDD = 50.0, 40.0
# Tenth and ninetieth percentile of blended over own average dose, 960 people (document 6).
TDD_SWING = (0.72, 1.14)
# The abstract's reference glucose, and the range over which the quartic is used.
REF_GLUCOSE, CURVE_RANGE = 180.0, (100.0, 380.0)
# The poster's fit for glucose in [0, 100] mg/dL: constant, linear and squared terms.
_QUADRATIC = (98.03, -1.077, 0.008868)
RESULTS = Path(__file__).resolve().parents[1] / "results"
# The INV-012 cohorts that hold carbohydrate records, by their name in its results file.
OWN_COHORTS = {"Loop": "Loop", "v5_trio": "Trio", "v7_oref0": "oref0", "ReplaceBG": "no algorithm (REPLACE-BG)"}
OWN_BINS = (100, 150, 200, 250, 300)
TARGET = 100.0
# Sensitivity x daily dose in entered settings: adults, and the pooled trial record.
ENTERED = {"adults": 1915.0, "pooled, 51% children": 2067.0}


def published(g: float) -> float:
    """The poster's fit in mg/dL per unit: the quadratic up to 100 mg/dL, the quartic above."""
    return float(quartic_isf(g)) if g > 100.0 else float(np.polyval(_QUADRATIC[::-1], g))


def own_levels(path: Path = RESULTS / "inv012_peak_trough.md") -> dict:
    """Peak-to-trough sensitivity x daily dose on this repository's records: people, then
    (median, lower, upper) by glucose bin, from the positive-insulin table of INV-012."""
    block = path.read_text().split("## Median over events with positive insulin used")[1].split("\n## ")[0]
    out = {}
    for row in block.splitlines():
        cells = [c.strip() for c in row.split("|")]
        if len(cells) > 8 and cells[1] in OWN_COHORTS:
            out[OWN_COHORTS[cells[1]]] = (int(cells[2]), {
                g: tuple(float(x) for x in re.findall(r"-?\d+", c)) for g, c in zip(OWN_BINS, cells[4:9]) if c != "-"})
    return {n: out[n] for n in OWN_COHORTS.values()}


def curve_shape(bg):
    """The published quartic relative to its value at 180 mg/dL, held outside its range."""
    return quartic_isf(np.clip(np.asarray(bg, dtype=float), *CURVE_RANGE)) / quartic_isf(REF_GLUCOSE)


def curve_k(bg, product: float = 1915.0):
    """Sensitivity x daily dose: the entered product at 180 mg/dL, shaped by the quartic."""
    return product * curve_shape(bg)


def curve_profile(bg, profile_isf):
    """The alternative link: the profile value at 180 mg/dL, shaped by the quartic."""
    return profile_isf * curve_shape(bg)


def curve_isf(bg, tdd, factor: float = 1.0):
    """The curve as a sensitivity: the stock equation with its glucose term swapped."""
    return curve_k(bg) / (factor * tdd)


def units_to_target(isf_at, bg: float, read: str = "start", target: float = TARGET) -> float:
    """Units for the distance to target, with the sensitivity read at the starting glucose,
    at the midpoint (how the curve files its events) or along the whole path."""
    if read == "start":
        return (bg - target) / float(isf_at(bg))
    if read == "midpoint":
        return (bg - target) / float(isf_at((bg + target) / 2.0))
    path = np.linspace(target, bg, 2001)
    return float(np.trapezoid(1.0 / isf_at(path), path))


def stock_v1(bg, tdd, divisor: float, factor: float = 1.0):
    """AndroidAPS: 1800 / (TDD x ln(BG / divisor + 1)), the adjustment factor scaling TDD."""
    return 1800.0 / (factor * tdd * np.log(np.asarray(bg, dtype=float) / divisor + 1.0))


def trio_log(bg, tdd, profile_isf, insulin_factor: float = 55.0, factor: float = 0.8,
             lo: float = 0.7, hi: float = 1.2):
    """Trio logarithmic: the same equation as a ratio of the profile value, clamped."""
    ratio = profile_isf * factor * tdd * np.log(np.asarray(bg, dtype=float) / insulin_factor + 1.0) / 1800.0
    return profile_isf / np.clip(ratio, lo, hi)


def trio_sigmoid(bg, profile_isf, target: float = 100.0, factor: float = 0.5,
                 tdd_ratio: float = 1.0, lo: float = 0.7, hi: float = 1.2):
    """Trio sigmoid: ratio 1 at target, bounded by the autosens limits."""
    offset = np.log(1.0 / (hi - 1.0) - lo / (hi - 1.0))
    exponent = (np.asarray(bg, dtype=float) - target) * 0.0555 * factor * tdd_ratio + offset
    return profile_isf / ((hi - lo) / (1.0 + np.exp(-exponent)) + lo)


def candidate(bg, profile_isf, tdd, tdd_ref, b: float = 1.0, pivot: float = 100.0,
              divisor: float = 55.0, lo: float = 0.7, hi: float = 1.2,
              high_from: float | None = None):
    """Profile value over a clamped dose ratio, eased below the pivot only."""
    bg = np.asarray(bg, dtype=float)
    level = profile_isf / np.clip((tdd / tdd_ref) ** b, lo, hi)
    log = lambda g: np.log(g / divisor + 1.0)
    out = level * log(pivot) / log(np.minimum(bg, pivot))
    if high_from is not None:
        out = out * log(high_from) / log(np.maximum(bg, high_from))
    return out


def main() -> None:
    bg = np.array(BG_POINTS)
    lines = [
        "# INV-013: the curve against the stock dynamic ISF",
        "",
        "## Sensitivity x daily dose by glucose",
        "",
        "The stock equation is 1800 / ln(glucose / divisor + 1). The curve is",
        "1915 x quartic(glucose) / quartic(180), held below 100 and above 380 mg/dL.",
        "The published column is the poster's fit in mg/dL per unit: the quadratic up to",
        "100 mg/dL and the quartic above. The curve uses the quartic's value at 100 (80.6).",
        "",
        "| Glucose, mg/dL | Stock, divisor 55 | Stock, divisor 65 | Stock, divisor 75 | Published | Curve | Curve over stock 55 | Curve over stock 75 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for g in bg:
        st = [float(stock_v1(g, 1.0, d)) for d in DIVISORS.values()]
        k = float(curve_k(g))
        q = f"{published(g):.1f}"
        lines.append(f"| {g:.0f} | {st[0]:.0f} | {st[1]:.0f} | {st[2]:.0f} | {q} | {k:.0f} | "
                     f"{k / st[0]:.2f} | {k / st[2]:.2f} |")
    lines += ["", "Ratio of the value at 100 to the value at 200 mg/dL: curve "
              f"{float(curve_k(100.0) / curve_k(200.0)):.2f}, stock "
              f"{float(stock_v1(100.0, 1.0, 55.0) / stock_v1(200.0, 1.0, 55.0)):.2f}. At 300 over 200: curve "
              f"{float(curve_k(300.0) / curve_k(200.0)):.2f}, stock "
              f"{float(stock_v1(300.0, 1.0, 55.0) / stock_v1(200.0, 1.0, 55.0)):.2f}.",
              f"Stock level at 100 mg/dL, divisor 75 over divisor 55: "
              f"{float(stock_v1(100.0, 1.0, 75.0) / stock_v1(100.0, 1.0, 55.0)):.2f}."]

    fine = np.linspace(100.0, 380.0, 2801)
    cross = lambda other: fine[np.argmin(np.abs(curve_k(fine) - other))]
    lines += ["", "Glucose at which the curve equals the stock equation: "
              + ", ".join(f"{cross(stock_v1(fine, 1.0, d)):.0f} mg/dL (divisor {d:.0f})"
                          if np.min(np.abs(curve_k(fine) - stock_v1(fine, 1.0, d))) < 1.0
                          else f"none from 100 to 380 mg/dL (divisor {d:.0f})"
                          for d in DIVISORS.values()) + "."]

    q = lambda g: float(quartic_isf(g))
    lines += ["", "## Reference dose", "",
              "The published curve carries no daily dose. Linking it to one needs the dose at",
              "which the published values apply, which is a level divided by the quartic.", "",
              "| Level taken | Reference dose, U | Curve x daily dose at 200 mg/dL |", "|---|---:|---:|"]
    levels = {"Adult entered product 1915 at 180 mg/dL (this draft)": (ENTERED["adults"], REF_GLUCOSE),
              "Pooled entered product 2067 at 180 mg/dL": (ENTERED["pooled, 51% children"], REF_GLUCOSE),
              "1800 rule at 105 mg/dL (March 2026 work)": (1800.0, 105.0)}
    own = own_levels()
    levels |= {f"Own peak-to-trough level under {n} at 200 mg/dL": (v[200][0], 200.0) for n, (_, v) in own.items()}
    for name, (level, at) in levels.items():
        lines.append(f"| {name} | {level / q(at):.1f} | {level / q(at) * q(200.0):.0f} |")
    lines += ["", "The adult product carries a 95% interval of 1828 to 1986, which at 200 mg/dL is "
              f"{float(curve_k(200.0, 1828.0)):.0f} to {float(curve_k(200.0, 1986.0)):.0f}.",
              f"At 180 mg/dL the TDD link gives {float(curve_isf(180.0, 20.0)):.1f} mg/dL per unit on 20 U a day and "
              f"{float(curve_isf(180.0, 80.0)):.1f} on 80 U a day. The published value there is {q(180.0):.1f} for both."]
    lines += ["", f"Dose multiple against the reference glucose: {1.0 / float(curve_shape(120.0)):.2f} at 120 mg/dL "
              f"and {1.0 / float(curve_shape(340.0)):.2f} at 340 mg/dL. The abstract gives about a half and about double."]

    lines += ["", "## The published curve against the measured levels", "",
              "Measured sensitivity x daily dose in each glucose bin (INV-012, median over events with",
              "positive insulin used, 95% interval over people) divided by the published value at that",
              "glucose. The result is the daily dose, in U, at which the published curve returns the",
              "measurement. A row that stays level means the published shape fits that cohort.", "",
              "| Cohort | People | " + " | ".join(str(g) for g in OWN_BINS) + " |", "|---|---:|" + "---:|" * len(OWN_BINS)]
    for name, (people, v) in own.items():
        lines.append(f"| {name} | {people} | " + " | ".join(
            "{:.0f} ({:.0f} to {:.0f})".format(*(x / published(g) for x in v[g])) if g in v else "-"
            for g in OWN_BINS) + " |")

    cols = {
        "AndroidAPS, divisor 55": stock_v1(bg, TDD, 55.0),
        "Trio logarithmic": trio_log(bg, TDD, PROFILE_ISF),
        "Trio sigmoid": trio_sigmoid(bg, PROFILE_ISF),
        "Curve, TDD link": curve_isf(bg, TDD),
        "Curve, profile link": curve_profile(bg, PROFILE_ISF),
        "Profile, eased below 100": candidate(bg, PROFILE_ISF, TDD, TDD),
    }
    lines += ["", "## One person", "",
              f"Profile sensitivity {PROFILE_ISF:.0f} mg/dL per unit, total daily dose {TDD:.0f} U, recent dose",
              "equal to the person's own average. AndroidAPS at an adjustment factor of 100%, Trio at",
              "its defaults (0.8 logarithmic, 0.5 sigmoid, limits 0.7 and 1.2, insulin factor 55).",
              "Values are mg/dL per unit.", "",
              "| Glucose, mg/dL | " + " | ".join(cols) + " |", "|---:|" + "---:|" * len(cols)]
    for i, g in enumerate(bg):
        lines.append(f"| {g:.0f} | " + " | ".join(f"{v[i]:.1f}" for v in cols.values()) + " |")

    # the worked comparisons below use the profile link, the candidate since INV-015
    link = lambda g: curve_profile(g, PROFILE_ISF)
    meets = lambda other: fine[np.argmin(np.abs(link(fine) - other))]
    for other in (stock_v1(fine, TDD, 55.0), stock_v1(fine, TDD, 75.0),
                  trio_log(fine, TDD, PROFILE_ISF), trio_sigmoid(fine, PROFILE_ISF)):
        assert np.min(np.abs(link(fine) - other)) < 0.05, "no crossing from 100 to 380 mg/dL"
    lines += ["", "The profile link over AndroidAPS at 100% for this person: "
              + ", ".join(f"{float(link(g) / stock_v1(g, TDD, 55.0)):.2f} at {g:.0f}" for g in (100.0, 150.0, 200.0, 300.0))
              + f" mg/dL. The two are equal at {meets(stock_v1(fine, TDD, 55.0)):.0f} mg/dL with a divisor of 55 and at "
              f"{meets(stock_v1(fine, TDD, 75.0)):.0f} with 75. The profile link equals Trio's logarithmic mode at "
              f"{meets(trio_log(fine, TDD, PROFILE_ISF)):.0f} mg/dL and its sigmoid at {meets(trio_sigmoid(fine, PROFILE_ISF)):.0f}."]
    inside = fine[(link(fine) <= PROFILE_ISF / 0.7) & (link(fine) >= PROFILE_ISF / 1.2)]
    pinned = fine[np.isclose(trio_log(fine, TDD, PROFILE_ISF), PROFILE_ISF / 1.2)].min()
    lines += ["", "Trio's logarithmic mode is the AndroidAPS equation divided by its adjustment factor, held inside",
              f"its limits. At the default factor of 0.8 it is {100 * (1 / 0.8 - 1):.0f}% above AndroidAPS at 100% until the upper limit",
              f"binds, which for this person is from {pinned:.0f} mg/dL. With the factor at 1.0 and the limits wide it",
              f"returns the AndroidAPS value: {float(trio_log(200.0, TDD, PROFILE_ISF, factor=1.0, lo=0.01, hi=100.0)):.1f} "
              f"against {float(stock_v1(200.0, TDD, 55.0)):.1f} mg/dL per unit at 200 mg/dL."]
    lines += ["", f"Inside Trio's default limits ({PROFILE_ISF / 1.2:.1f} to {PROFILE_ISF / 0.7:.1f}) the curve on the profile link "
              f"runs only from {inside.min():.0f} to {inside.max():.0f} mg/dL for this person."]

    lines += ["", "## Units to bring that person to 100 mg/dL", "",
              "Distance from 100 mg/dL over the sensitivity, with the curve on the profile link. This is",
              "the arithmetic of the sensitivity alone and leaves out insulin on board, predictions and",
              "delivery limits. An event is taken to be filed under the glucose halfway down its fall,",
              "which the abstract does not state. On that reading, using the curve at the starting glucose asks for more insulin",
              "than its construction implies.", "",
              "| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Profile | Curve, read at start | Curve, read at midpoint | Curve, along the path |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
    stock = lambda g: stock_v1(g, TDD, 55.0)
    trio = lambda g: trio_log(g, TDD, PROFILE_ISF)
    curve = link
    for g in (150.0, 200.0, 250.0, 300.0):
        lines.append(f"| {g:.0f} | {units_to_target(stock, g):.2f} | {units_to_target(trio, g):.2f} | "
                     f"{(g - TARGET) / PROFILE_ISF:.2f} | {units_to_target(curve, g):.2f} | "
                     f"{units_to_target(curve, g, 'midpoint'):.2f} | {units_to_target(curve, g, 'path'):.2f} |")
    extra = [units_to_target(curve, g) / units_to_target(curve, g, "midpoint") - 1.0 for g in (150.0, 200.0, 250.0, 300.0)]
    lines += ["", f"Reading at the start asks for {100 * min(extra):.0f} to {100 * max(extra):.0f}% more than reading at the midpoint."]
    same = fine[np.argmin(np.abs(curve_isf(fine, TDD) - PROFILE_ISF))]
    lines += [f"For this person the TDD link equals the profile value at {same:.0f} mg/dL; the profile link does so at {REF_GLUCOSE:.0f}."]

    lines += ["", "## Movement with recent dose", "",
              "The daily dose in the link is this draft's to define. Dividing by",
              "the existing blend moves the value by the same multiplier as the stock equation:",
              ", ".join(f"{1.0 / r:.2f} at a blend of {r:.2f} of the person's average" for r in TDD_SWING) + "."]

    out = RESULTS / "inv013_comparator.md"
    out.write_text("\n".join(lines) + "\n")
    print(out.read_text())


if __name__ == "__main__":
    main()
