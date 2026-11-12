#!/usr/bin/env python3
"""The poster's quartic refitted to this repository's peak-to-trough records.

Romero-Ugalde et al. (doi:10.2337/db25-920-P) fit a quartic in glucose to their binned
sensitivity for glucose in (100, 380] mg/dL. `inv013.candidate` set each of five measured
bins of INV-012 against that fit. This script does the other half: it fits the same form,
with all five coefficients free, to the INV-012 bins of each cohort, and compares the
result with the published fit held to its own shape at the best single scale.

What is fitted is the cohort median, in each 10 mg/dL bin of midpoint glucose, of each
person's median sensitivity x daily dose over events with positive insulin used. That
quantity is an outcome ratio: an event ends where glucose stopped falling, and simulated
people with a constant sensitivity return a ratio that falls with glucose
(results/inv012_peak_trough.md). The refit describes that ratio on these records. It is
a description of the measurement and makes no claim about physiological sensitivity.

Modelling choices, fixed before the run:
  Bins from 100 to 380 mg/dL, the published range, that hold 10 people (the INV-012
  rule). Trio and oref0 stop near 200 and 240 mg/dL, so their refits cover less.
  Least squares on the bin medians, each weighted by the square root of its people,
  since a median over 12 people is noisier than one over 800. The poster does not
  state its weighting.
  The published fit is given one free number, a scale, by the same weighted least
  squares. That scale is the daily dose, in U, at which the published equation best
  returns the cohort.
  Shape is reported relative to the value at 180 mg/dL, the reference glucose of the
  abstract. Intervals are 95%, from resampling people 2,000 times.

Expected before the run, from the bin medians already tabulated on 11 November 2026.
Where no algorithm doses (REPLACE-BG) the published shape matched three bins from 100 to
200 mg/dL at about 100 U, so the refit relative to 180 mg/dL should sit near the
published 2.28 at 100 and 1.30 at 150, and above the published 0.70 and 0.60 at 250 and
300, where the measured values level off and the published fit keeps falling. (This
paragraph first gave 1.42 for 150 mg/dL, an arithmetic slip found when the run printed
the published row. The directions stated were not changed.) Under Loop
the refit should be above the published shape at 100, below it at 150 and above it from
250. The published fit at its best scale should miss the Loop bins by more than it
misses the REPLACE-BG bins, and the refit should miss both by less, having four more
free numbers. A quartic through a curve that falls and then levels off may bend upward
at the top of its range, which would be the form and not the records.

The run agreed on each point. The upward bend appeared under Loop between 200 and 250
mg/dL, where the bin medians themselves go from 2,182 to 2,025, so the refit's values
above 200 mg/dL are read with the bins beside them.

The published fit is taken as its quartic in every bin, the bin at 100 mg/dL included,
as `inv013.candidate.curve_shape` does. The poster's quadratic gives 79.0 there and the
quartic 80.6.

Reads results/inv012_peak_trough_person_bins.parquet, which is not in git (parquet files
are ignored) and is rewritten by `python3 -m inv012.peak_trough`. The bin medians are
checked against the committed results/inv012_peak_trough.md before anything is fitted.

Run as `python3 -m inv013.refit` to write results/inv013_refit.md.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inv008.dynisf import quartic_isf
from inv013.candidate import OWN_COHORTS, REF_GLUCOSE, RESULTS, own_levels

FIT_RANGE = (100, 380)
MIN_PEOPLE = 10
N_BOOT = 2000
POINTS = (100, 150, 200, 250, 300)


def bin_medians(wide: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Median over people in each bin (columns), and the people behind it."""
    n = np.isfinite(wide).sum(axis=0)
    with np.errstate(all="ignore"):
        return np.nanmedian(wide, axis=0), n


def fit_quartic(g, y, n) -> np.ndarray:
    """Quartic in glucose / 100, highest power first, weighted by root people."""
    return np.polyfit(np.asarray(g, float) / 100.0, y, 4, w=np.sqrt(n))


def fit_scale(g, y, n) -> float:
    """The one number that best lays the published fit over the bins: a daily dose in U."""
    q = quartic_isf(np.asarray(g, float))
    return float(np.sum(n * y * q) / np.sum(n * q * q))


def miss(fitted, y, n) -> float:
    """Root mean square miss as a share of the bin median, weighted by people."""
    return float(np.sqrt(np.sum(n * ((fitted - y) / y) ** 2) / np.sum(n)))


def refit(wide: pd.DataFrame, n_boot: int = N_BOOT) -> dict:
    """Refit one cohort from its people x bins table, with intervals over people."""
    med, n = bin_medians(wide.to_numpy(float))
    keep = (wide.columns >= FIT_RANGE[0]) & (wide.columns <= FIT_RANGE[1]) & (n >= MIN_PEOPLE)
    g, x = wide.columns[keep].to_numpy(float), wide.to_numpy(float)[:, keep]

    def one(sample: np.ndarray) -> np.ndarray:
        y, k = bin_medians(sample)
        ok = np.isfinite(y)
        c = fit_quartic(g[ok], y[ok], k[ok])
        at = lambda v: np.polyval(c, np.asarray(v, float) / 100.0)
        pts = [p for p in POINTS if p <= g.max()]
        return np.concatenate([[at(REF_GLUCOSE), fit_scale(g[ok], y[ok], k[ok])],
                               at(pts), at(pts) / at(REF_GLUCOSE)])

    rng = np.random.default_rng(0)
    draws = np.array([one(x[rng.integers(0, len(x), len(x))]) for _ in range(n_boot)])
    est = one(x)
    val = [(float(e), float(lo), float(hi))
           for e, lo, hi in zip(est, *np.percentile(draws, [2.5, 97.5], axis=0))]
    pts = [p for p in POINTS if p <= g.max()]
    y, k = med[keep], n[keep]
    c = fit_quartic(g, y, k)
    return dict(bins=g, people=int(len(x)), med=y, n=k, coef=c / 100.0 ** np.arange(4, -1, -1),
                at_ref=val[0], scale=val[1],
                level=dict(zip(pts, val[2:2 + len(pts)])), shape=dict(zip(pts, val[2 + len(pts):])),
                miss_refit=miss(np.polyval(c, g / 100.0), y, k),
                miss_published=miss(fit_scale(g, y, k) * quartic_isf(g), y, k))


def main() -> None:
    pb = pd.read_parquet(RESULTS / "inv012_peak_trough_person_bins.parquet")
    pb = pb.replace([np.inf, -np.inf], np.nan)
    own, fits = own_levels(), {}
    for cohort, name in OWN_COHORTS.items():
        wide = pb[pb.cohort == cohort].pivot(index="subject", columns="mid", values="med_pos")
        med, _ = bin_medians(wide.to_numpy(float))
        for b, v in own[name][1].items():      # the parquet and the committed table are one run
            assert round(float(med[list(wide.columns).index(b)])) == v[0], (name, b)
        fits[name] = refit(wide)

    f = lambda v, d=0: f"{v[0]:.{d}f} ({v[1]:.{d}f} to {v[2]:.{d}f})"
    cell = lambda r, key, d: [f(r[key][p], d) if p in r[key] else "-" for p in POINTS]
    head = "| Cohort | People | " + " | ".join(str(p) for p in POINTS) + " |"
    rule = "|---|---:|" + "---:|" * len(POINTS)
    L = ["# INV-013: the published quartic refitted to the INV-012 bins", "",
         "Generated by `inv013/refit.py`. The fitted quantity is the cohort median of each person's",
         "median sensitivity x daily dose (events with positive insulin used) in 10 mg/dL bins of",
         "midpoint glucose. It is an outcome ratio. Intervals are 95% over people.", "",
         "## Fit quality", "",
         "The published fit keeps its shape and is given one scale. The refit has five free",
         "coefficients. The miss is the root mean square difference from the bin medians, as a",
         "percentage of each median, weighted by people.", "",
         "| Cohort | People | Bins, mg/dL | Best scale of the published fit, U | Miss of the published fit | Miss of the refit | Refit at 180 mg/dL |",
         "|---|---:|---|---:|---:|---:|---:|"]
    for name, r in fits.items():
        L.append(f"| {name} | {r['people']} | {r['bins'].min():.0f} to {r['bins'].max():.0f} | {f(r['scale'])} | "
                 f"{100 * r['miss_published']:.0f}% | {100 * r['miss_refit']:.0f}% | {f(r['at_ref'])} |")
    L += ["", "## Shape relative to 180 mg/dL", "", head, rule,
          "| Published fit | | " + " | ".join(f"{float(quartic_isf(p) / quartic_isf(REF_GLUCOSE)):.2f}" for p in POINTS) + " |"]
    L += [f"| {name} | {r['people']} | " + " | ".join(cell(r, "shape", 2)) + " |" for name, r in fits.items()]
    L += ["", "## Refitted sensitivity x daily dose", "", head, rule]
    L += [f"| {name} | {r['people']} | " + " | ".join(cell(r, "level", 0)) + " |" for name, r in fits.items()]
    L += ["", "## Refitted coefficients", "",
          "Sensitivity x daily dose = a + b G + c G^2 + d G^3 + e G^4, G in mg/dL, valid over the bins fitted.", "",
          "| Cohort | a | b | c | d | e |", "|---|---:|---:|---:|---:|---:|"]
    L += [f"| {name} | " + " | ".join(f"{v:.4g}" for v in r["coef"][::-1]) + " |" for name, r in fits.items()]
    L += ["", "## Bin medians relative to the refit at 180 mg/dL", "",
          "Each fitted bin's median over the cohort's refit at 180 mg/dL, with the people behind it in brackets.", "",
          "| Bin, mg/dL | " + " | ".join(fits) + " |", "|---:|" + "---:|" * len(fits)]
    for b in sorted(set().union(*(r["bins"] for r in fits.values()))):
        L.append(f"| {b:.0f} | " + " | ".join(
            next((f"{m / r['at_ref'][0]:.2f} ({k})" for g, m, k in zip(r["bins"], r["med"], r["n"]) if g == b), "-")
            for r in fits.values()) + " |")
    out = RESULTS / "inv013_refit.md"
    out.write_text("\n".join(L) + "\n")
    print(out.read_text())


if __name__ == "__main__":
    main()
