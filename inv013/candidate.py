#!/usr/bin/env python3
"""Reference form of the candidate sensitivity in docs/dynamic-isf-8-proposal.md.

The comparators are written from the stock sources read on 10 November 2026:
AndroidAPS upstream master f77abd9b82 (OpenAPSSMBPlugin.calculateRawDynIsf) and
Trio upstream main 0b29cd534 (DynamicISF.calculate). `inv008.dynisf.isf_v1` is
the Boost form with a glucose cap and a velocity term, so the stock equation is
restated here and tested against it with both switched off.

The candidate is a modelling choice assembled from findings, and no parameter in
it was fitted. Level: the entered profile value, scaled by recent dose against
the person's own average and held inside the autosens limits. Glucose: the stock
logarithm below a pivot, where it raises the value, and nothing between the
pivot and `high_from`. Above `high_from` the stock logarithm is an option, off
unless a glucose is passed.

Run as `python3 -m inv013.candidate` to write results/inv013_comparator.md.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

DIVISORS = {"rapid, 75 min peak": 55.0, "ultra-rapid, 55 min peak": 65.0, "45 min peak": 75.0}
BG_POINTS = (70.0, 80.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0)
# The worked person: the simulator's, close to the adult entered median of 1915.
PROFILE_ISF, TDD = 50.0, 40.0
# Tenth and ninetieth percentile of blended over own average dose, 960 people (document 6).
TDD_SWING = (0.72, 1.14)


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
    cols = {f"AndroidAPS, divisor {int(d)}": stock_v1(bg, TDD, d) for d in DIVISORS.values()}
    cols["Trio logarithmic"] = trio_log(bg, TDD, PROFILE_ISF)
    cols["Trio sigmoid"] = trio_sigmoid(bg, PROFILE_ISF)
    cols["Candidate"] = candidate(bg, PROFILE_ISF, TDD, TDD)
    cols["Candidate, option above 200"] = candidate(bg, PROFILE_ISF, TDD, TDD, high_from=200.0)

    lines = [
        "# INV-013: what each equation calculates for one person",
        "",
        f"Profile sensitivity {PROFILE_ISF:.0f} mg/dL per unit, total daily dose {TDD:.0f} U, recent dose",
        "equal to the person's own average. AndroidAPS at an adjustment factor of 100%, Trio at",
        "its defaults (0.8 logarithmic, 0.5 sigmoid, limits 0.7 and 1.2, insulin factor 55).",
        "Values are mg/dL per unit.",
        "",
        "| Glucose, mg/dL | " + " | ".join(cols) + " |",
        "|---|" + "---:|" * len(cols),
    ]
    for i, g in enumerate(bg):
        lines.append(f"| {g:.0f} | " + " | ".join(f"{v[i]:.1f}" for v in cols.values()) + " |")

    lines += ["", "## Level the stock equation implies at 100 mg/dL", "",
              "| Insulin setting | Divisor | Sensitivity x daily dose |", "|---|---:|---:|"]
    for name, d in DIVISORS.items():
        lines.append(f"| {name} | {d:.0f} | {float(stock_v1(100.0, 1.0, d)):.0f} |")

    lines += ["", "## Correction size against the stock equation", "",
              "Stock AndroidAPS value at divisor 55 divided into the candidate's: above 1 the",
              "candidate asks for less insulin for the same distance from target.", "",
              "| Glucose, mg/dL | Candidate | Candidate, option above 200 |", "|---|---:|---:|"]
    stock = cols["AndroidAPS, divisor 55"]
    for i, g in enumerate(bg):
        lines.append(f"| {g:.0f} | {cols['Candidate'][i] / stock[i]:.2f} | "
                     f"{cols['Candidate, option above 200'][i] / stock[i]:.2f} |")
    lines += ["", f"Stock ratio of the value at 300 to the value at 200 mg/dL: "
              + ", ".join(f"{float(stock_v1(300.0, 1.0, d) / stock_v1(200.0, 1.0, d)):.2f} (divisor {d:.0f})"
                          for d in DIVISORS.values()) + "."]

    lines += ["", "## Movement with recent dose", "",
              "Multiplier on the sensitivity when the blended dose sits at the tenth and ninetieth",
              "percentile of its ratio to the person's own average.", "",
              "| Blended over average | Stock | Candidate, exponent 1 | Candidate, exponent 0.65 |",
              "|---:|---:|---:|---:|"]
    for r in TDD_SWING:
        lines.append(f"| {r:.2f} | {1.0 / r:.2f} | "
                     f"{float(candidate(150.0, 1.0, r, 1.0)):.2f} | "
                     f"{float(candidate(150.0, 1.0, r, 1.0, b=0.65)):.2f} |")

    lines += ["", "## Units to bring the worked person to 100 mg/dL", "",
              "Distance from 100 mg/dL divided by the sensitivity at the starting glucose. This is",
              "the arithmetic of the sensitivity alone and leaves out insulin on board, predictions",
              "and delivery limits.", "",
              "| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Candidate |",
              "|---:|---:|---:|---:|"]
    for i, g in enumerate(bg):
        if g > 100.0:
            lines.append(f"| {g:.0f} | " + " | ".join(
                f"{(g - 100.0) / cols[c][i]:.2f}" for c in ("AndroidAPS, divisor 55", "Trio logarithmic", "Candidate")) + " |")

    lines += ["", "## Shape alone", "",
              "Each value divided by the same equation's value at 100 mg/dL, divisor 55. The last two",
              "columns are the candidate at the other two divisors.", "",
              "| Glucose, mg/dL | Stock | Candidate | Candidate, divisor 65 | Candidate, divisor 75 |",
              "|---:|---:|---:|---:|---:|"]
    for g in bg:
        lines.append(f"| {g:.0f} | {float(stock_v1(g, 1.0, 55.0) / stock_v1(100.0, 1.0, 55.0)):.2f} | "
                     + " | ".join(f"{float(candidate(g, 1.0, 1.0, 1.0, divisor=d)):.2f}" for d in (55.0, 65.0, 75.0)) + " |")

    lines += ["", "## Figures quoted in the proposal", "",
              f"- Level at 100 mg/dL, divisor 75 over divisor 55: "
              f"{float(stock_v1(100.0, 1.0, 75.0) / stock_v1(100.0, 1.0, 55.0)):.2f}.",
              f"- Trio logarithmic over AndroidAPS (divisor 55) at 200 mg/dL: "
              f"{float(trio_log(200.0, TDD, PROFILE_ISF) / stock_v1(200.0, TDD, 55.0)):.2f}.",
              f"- Candidate at 80 mg/dL with recent dose 1.10 of average: "
              f"{float(candidate(80.0, PROFILE_ISF, 1.1 * TDD, TDD)):.1f}; from 100 mg/dL upward: "
              f"{float(candidate(100.0, PROFILE_ISF, 1.1 * TDD, TDD)):.1f}.",
              f"- Spread of the easing at 70 mg/dL across divisors 55 to 75: "
              f"{float(candidate(70.0, 1.0, 1.0, 1.0, divisor=75.0) - candidate(70.0, 1.0, 1.0, 1.0)):.3f}."]

    out = Path(__file__).resolve().parents[1] / "results" / "inv013_comparator.md"
    out.write_text("\n".join(lines) + "\n")
    print(out.read_text())


if __name__ == "__main__":
    main()
