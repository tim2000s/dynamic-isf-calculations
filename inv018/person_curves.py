#!/usr/bin/env python3
"""Each REPLACE-BG person's own curve from the interval technique, and the range across people.

The published curve (Romero-Ugalde et al., doi:10.2337/db25-920-P) is one curve for
everybody. This script applies the technique behind it, as `inv012.peak_trough` implements
it, to every person in REPLACE-BG separately and draws the curves side by side. REPLACE-BG
is the cohort with no algorithm dosing: 196 adults on a pump and a bolus calculator.

What a curve is. An event runs from a glucose maximum away from a meal to the next
minimum. Its sensitivity is the fall divided by the insulin used above scheduled basal.
Events are filed by midpoint glucose in 10 mg/dL bins, and a person's value in a bin is
the median over their events with positive insulin used, in mg/dL per unit. No daily dose
enters: `inv012` multiplies by it and this script divides it out again.

The reading is an outcome ratio. An event ends where glucose stopped falling, so events
are selected on what glucose did. In `results/inv012_peak_trough.md` simulated people
whose sensitivity is constant return a ratio that falls with glucose and sits above the
truth. The curves here describe the measurement on these records. Their height is not a
sensitivity to dose from.

Each curve is drawn twice: as measured, and divided by the person's own reading at
140 mg/dL, the anchor of the proposal. The second drawing shows how far the shape differs
between people once the height is taken out.

Modelling choices, fixed before the run:
  Bins from 100 to 380 mg/dL, the range of the published quartic. A person's bin needs 5
  events (the `inv012` rule) and a percentile over people needs 10 people in the bin.
  A person's reading at 140 is the median of their bins from 130 to 150 mg/dL, and at 250
  the median of their bins from 230 to 270, since one bin of five events is noisy.
  The entered ISF is the person's median over the calculator corrections of INV-015.
  People are split into thirds by it.
  Intervals are 95%, from resampling people 2,000 times.

Expected before the run. An uncommitted look at the `inv012` bins on 12 November 2026
came first, so these expectations are not blind. The readings at 140 mg/dL sit two to
three times above the entered ISF. Across thirds of entered ISF the share of the 140
reading kept at 250 is about equal (0.56 to 0.61), which is what curves that differ by a
multiple give, and the absolute drop is larger in the higher thirds, where curves that
differ by a shift would give equal drops. The spread across people is expected to be
wide as measured and narrower once each curve is divided by its own reading at 140.

The run agreed on the first two and only partly on the last. The reading at 140 was 2.5
to 2.9 times the entered ISF in each third. The share kept was 0.58, 0.56 and 0.61, and
the drop in the highest third was 2.32 (1.23 to 5.51) times that in the lowest. Dividing
by the reading at 140 narrowed the spread at 180 and 200 mg/dL (75th over 25th
percentile 1.80 and 2.02, from 2.59 and 2.95) and did not narrow it at 250 (2.82 from
2.56). Single curves cross one another throughout. A bin holds as few as five events, so
this run cannot say how much of the difference in shape between people would repeat on
a second sample of the same person's events.

    python3 -m inv018.person_curves            # from the repository root, DYNISF_ROOT unset

Writes results/inv018_person_curves.md and charts/inv018/fig_person_curves.png. The
per-person bins go to results/inv018_person_bins.parquet first, which is not in git.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008.dynisf import quartic_isf
from inv009 import config, db
from inv012.peak_trough import MIN_PEOPLE, PROMINENCE, _run
from inv013.candidate import ANCHOR_GLUCOSE
from inv015.link import boot

COHORT = "ReplaceBG"
RANGE = (100, 380)
ANCHOR_BINS, HIGH_BINS, HIGH = (130, 150), (230, 270), 250
POINTS = (100, 140, 180, 200, 250, 300)
QS = (10, 25, 50, 75, 90)
THIRDS = ("Lowest third", "Middle third", "Highest third")
CHART = config.ROOT / "charts" / "inv018" / "fig_person_curves.png"
BLUE, ORANGE, INK, MUTED, SURFACE = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"


def curves(pb: pd.DataFrame) -> pd.DataFrame:
    """People x bins: each person's median fall per unit of insulin used, mg/dL per U."""
    pb = pb[(pb.mid >= RANGE[0]) & (pb.mid <= RANGE[1])]
    isf = (pb.med_pos / pb.tdd).replace([np.inf, -np.inf], np.nan)
    return pb.assign(isf=isf).pivot(index="subject", columns="mid", values="isf")


def band(wide: pd.DataFrame, lo: int, hi: int) -> pd.Series:
    """A person's reading over a band of bins: the median of their bin medians."""
    return wide.loc[:, lo:hi].median(axis=1)


def relative(wide: pd.DataFrame) -> pd.DataFrame:
    """Each curve divided by the person's own reading at the anchor."""
    return wide.div(band(wide, *ANCHOR_BINS), axis=0)


def spread(wide: pd.DataFrame) -> pd.DataFrame:
    """Percentiles over people in every bin that holds enough of them."""
    wide = wide.loc[:, wide.notna().sum() >= MIN_PEOPLE]
    out = pd.DataFrame(np.nanpercentile(wide, QS, axis=0).T, index=wide.columns, columns=QS)
    return out.assign(n=wide.notna().sum())


def by_entered(wide: pd.DataFrame, entered: pd.Series) -> tuple[pd.DataFrame, list]:
    """Readings at the anchor and at 250 mg/dL by thirds of entered ISF.

    A shift between people gives every third the same drop, and a multiple gives every
    third the same share kept. The last two numbers returned are the highest third over
    the lowest for each, so one of them should be 1.
    """
    p = pd.DataFrame({"at_anchor": band(wide, *ANCHOR_BINS), "at_high": band(wide, *HIGH_BINS),
                      "entered": entered}).dropna()
    p["third"] = pd.qcut(p.entered, 3, labels=THIRDS)
    drop = lambda d: (d.at_anchor - d.at_high).median()
    kept = lambda d: (d.at_high / d.at_anchor).median()
    rows = []
    for name, d in p.groupby("third", observed=True):
        iv = boot(d, lambda x: np.array([drop(x), kept(x)]))
        rows.append(dict(third=name, people=len(d), entered=d.entered.median(),
                         at_anchor=d.at_anchor.median(), at_high=d.at_high.median(),
                         over_entered=(d.at_anchor / d.entered).median(), drop=iv[0], kept=iv[1]))
    lo_hi = lambda x: (x[x.third == THIRDS[0]], x[x.third == THIRDS[-1]])
    contrast = boot(p, lambda x: np.array([drop(lo_hi(x)[1]) / drop(lo_hi(x)[0]),
                                           kept(lo_hi(x)[1]) / kept(lo_hi(x)[0])]))
    return pd.DataFrame(rows), contrast


def figure(wide: pd.DataFrame, path) -> None:
    g = np.arange(RANGE[0], RANGE[1] + 1, 5.0)
    panels = ((wide, quartic_isf(g), "Sensitivity, mg/dL per U", "As measured"),
              (relative(wide), quartic_isf(g) / quartic_isf(ANCHOR_GLUCOSE),
               f"Relative to own reading at {ANCHOR_GLUCOSE:.0f} mg/dL",
               "Each curve divided by its own reading at 140 mg/dL"))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor=SURFACE)
    for ax, (w, pub, ylabel, title) in zip(axes, panels):
        s = spread(w)
        for _, row in w.iterrows():
            ax.plot(w.columns, row, color=BLUE, lw=0.5, alpha=0.18)
        ax.fill_between(s.index, s[10], s[90], color=BLUE, alpha=0.12, lw=0, label="10th to 90th percentile of people")
        ax.fill_between(s.index, s[25], s[75], color=BLUE, alpha=0.25, lw=0, label="25th to 75th percentile")
        ax.plot(s.index, s[50], color=BLUE, lw=2, label="Median person")
        ax.plot(g, pub, color=ORANGE, lw=2, ls="--", label="Published curve")
        ax.plot([], [], color=BLUE, lw=0.5, alpha=0.6, label="One person")
        ax.set_ylim(0, 3 * s[50].iloc[0])      # single curves and the sparse top bins run off the scale
        ax.set_xlim(*RANGE)
        ax.set_facecolor(SURFACE)
        ax.set_title(title, loc="left", fontsize=11, color=INK)
        ax.set_xlabel("Midpoint glucose of the event, mg/dL", color=MUTED)
        ax.set_ylabel(ylabel, color=MUTED)
        ax.tick_params(colors=MUTED, length=0)
        ax.grid(color="#e4e3df", lw=0.6)
        for side in ax.spines.values():
            side.set_visible(False)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="lower center", ncol=5, frameon=False,
               fontsize=9, labelcolor=INK)
    fig.suptitle(f"Interval technique applied to each of {len(wide)} adults in REPLACE-BG",
                 x=0.01, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def write_md(pb: pd.DataFrame, wide: pd.DataFrame, thirds: pd.DataFrame, contrast: list, path) -> None:
    f = lambda v, d=1: f"{v[0]:.{d}f} ({v[1]:.{d}f} to {v[2]:.{d}f})"
    people = pb.drop_duplicates("subject")
    none_used = pb.pivot(index="subject", columns="mid", values="frac_nonpos").median()
    rel, q = relative(wide), lambda v: float(quartic_isf(v))
    med = lambda w, p: f(boot(w[[p]].dropna(), lambda x: np.array([x[p].median()]))[0], 1 if w is wide else 2)
    L = ["# INV-018: the interval technique applied to each person in REPLACE-BG", "",
         "Generated by `inv018/person_curves.py`. Each person's curve is their median fall per unit of",
         "insulin used above scheduled basal, over events with positive insulin used, in 10 mg/dL bins",
         "of the event's midpoint glucose. It is an outcome ratio: an event ends where glucose stopped",
         "falling, and simulated people with a constant sensitivity return a falling curve by the same",
         "technique (`results/inv012_peak_trough.md`).", "",
         f"{len(wide)} adults have a bin of 5 events between {RANGE[0]} and {RANGE[1]} mg/dL, out of "
         f"{len(people)} with any bin. They hold {int(people.n_events.sum())} events, a median of "
         f"{people.n_events.median():.0f} each, and a median of {wide.notna().sum(axis=1).median():.0f} "
         f"bins each.", "",
         "![Each person's curve](../charts/inv018/fig_person_curves.png)", "",
         "## Range across people, as measured", "",
         "Sensitivity in mg/dL per U. The median carries a 95% interval over people. The last column is",
         "the median person's share of events in the bin with no net insulin used, which the curve leaves out.", "",
         "| Midpoint glucose, mg/dL | People | 10th | 25th | Median | 75th | 90th | 75th over 25th | Published | No net insulin |",
         "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    s = spread(wide)
    L += [f"| {p} | {s.n[p]} | {s[10][p]:.0f} | {s[25][p]:.0f} | {med(wide, p)} | {s[75][p]:.0f} | {s[90][p]:.0f} | "
          f"{s[75][p] / s[25][p]:.2f} | {q(p):.1f} | {100 * none_used[p]:.0f}% |" for p in POINTS if p in s.index]
    L += ["", f"## Range across people, each curve divided by its own reading at {ANCHOR_GLUCOSE:.0f} mg/dL", "",
          f"The reading at {ANCHOR_GLUCOSE:.0f} is the median of the person's bins from {ANCHOR_BINS[0]} to {ANCHOR_BINS[1]} mg/dL.", "",
          "| Midpoint glucose, mg/dL | People | 10th | 25th | Median | 75th | 90th | 75th over 25th | Published |",
          "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    s = spread(rel)
    L += [f"| {p} | {s.n[p]} | {s[10][p]:.2f} | {s[25][p]:.2f} | {med(rel, p)} | {s[75][p]:.2f} | {s[90][p]:.2f} | "
          f"{s[75][p] / s[25][p]:.2f} | {q(p) / q(ANCHOR_GLUCOSE):.2f} |" for p in POINTS if p in s.index]
    L += ["", "## By thirds of entered ISF", "",
          f"Readings at {ANCHOR_GLUCOSE:.0f} and {HIGH} mg/dL (bins {HIGH_BINS[0]} to {HIGH_BINS[1]}), mg/dL per U. The drop and the share kept are",
          "medians of each person's own, with 95% intervals over people.", "",
          f"| Entered ISF | People | Median entered | Reading at {ANCHOR_GLUCOSE:.0f} | Reading at {HIGH} | Reading at {ANCHOR_GLUCOSE:.0f} over entered | Drop | Share kept |",
          "|---|---:|---:|---:|---:|---:|---:|---:|"]
    L += [f"| {r.third} | {r.people} | {r.entered:.0f} | {r.at_anchor:.0f} | {r.at_high:.0f} | {r.over_entered:.2f} | "
          f"{f(r.drop)} | {f(r.kept, 2)} |" for r in thirds.itertuples()]
    L += ["", f"Highest third over lowest: the drop is {f(contrast[0], 2)} times as large, where curves that differ by",
          f"a shift give 1. The share kept is {f(contrast[1], 2)} times as large, where curves that differ by a multiple",
          f"give 1. The published curve keeps {q(HIGH) / q(ANCHOR_GLUCOSE):.2f} of its value at {ANCHOR_GLUCOSE:.0f} by {HIGH} mg/dL."]
    path.write_text("\n".join(L) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="people, for a pilot; percentiles need 10 in a bin")
    ap.add_argument("--workers", type=int, default=config.WORKERS)
    a = ap.parse_args()
    ids = db.subjects(COHORT).subject_id.tolist()
    jobs = [("jaeb", (sid,), PROMINENCE) for sid in (ids[:a.limit] if a.limit else ids)]
    t0 = time.time()
    with mp.Pool(a.workers, maxtasksperchild=8) as pool:
        res = [r for r in pool.imap_unordered(_run, jobs, chunksize=1) if r is not None and len(r)]
    print(f"{len(jobs)} people, {len(res)} with a usable bin, {time.time() - t0:.0f} s")
    tag = "_limit" if a.limit else ""
    pb = pd.concat(res, ignore_index=True).sort_values(["subject", "mid"])
    pb.to_parquet(config.RESULTS / f"inv018_person_bins{tag}.parquet", index=False)

    wide = curves(pb).dropna(how="all")
    figure(wide, CHART.with_name(CHART.stem + tag + CHART.suffix))
    corr = pd.read_parquet(config.RESULTS / "inv009_correction_landing.parquet")
    entered = corr[corr.study == COHORT].drop_duplicates().groupby("sid").isf.median()
    out = config.RESULTS / f"inv018_person_curves{tag}.md"
    write_md(pb, wide, *by_entered(wide, entered), out)
    print(out.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
