#!/usr/bin/env python3
"""Does a person's own interval curve repeat? Each REPLACE-BG person's events in two halves.

`inv018.person_curves` drew one curve per person and found them far apart, both in
height and in shape. A bin holds as few as five events, so that spread could be
differences between people or sampling noise. This script splits each person's events
into two halves, builds a curve from each half by the same technique, and asks how well
one half's curve predicts the other's.

Two things are compared between halves. The height is the reading at 140 mg/dL. The
shape is the share of that reading kept at 190 and at 250 mg/dL. A per-user height needs
the first to repeat and a per-user shape needs the second.

Two ways of halving. Alternate events (first, third, fifth against second, fourth,
sixth) put both halves in the same months, so they differ by sampling alone. Earlier
against later events adds any drift over the study, and is the case a device would face:
a curve built from past events, used on later ones.

Statistics:
  The rank correlation across people between the two halves. It is close to the share
  of the spread between people, in half a person's events, that repeats.
  The same stepped up to all of a person's events, 2r / (1 + r) (Spearman-Brown), for
  alternate events only, where the halves are exchangeable.
  Two misses on the second half, each the typical factor between a prediction and the
  reading. One predicts from the person's own first half. The other gives everyone the
  median of the first halves. A curve evaluated per person earns its place where the
  first miss is the smaller.

Modelling choices, fixed before the run:
  Bins, bands and the five-event rule are those of `inv018.person_curves`, applied to
  each half, so a bin needs about ten events in the whole record to appear in both.
  The reading at 190 is the median of the person's bins from 180 to 200 mg/dL.
  A quantity is called repeatable when its reliability on all events is 0.7 or more with
  an interval that stays above 0.5, and not repeatable when the interval stays below
  0.5. Anything else is undecided. 0.7 is a convention for a measurement used on one
  person, and it is a choice.
  Intervals are 95%, from resampling people 2,000 times.

Known truth. `inv018/tests/test_split_half.py` builds people who share one shape and
differ only in height, and people whose shapes differ. The shape correlation has to be
near zero in the first and well above it in the second, with the height repeating in
both.

Expected before the run, having seen the curves of `inv018.person_curves` and no half.
The height repeats: people differ nearly three-fold between quartiles and the band at
140 holds dozens of events, so a correlation above 0.7 between alternate halves. The
shape repeats weakly, between 0 and 0.5 at 250 mg/dL, because the thirds of entered ISF
shared one shape and single curves cross one another. Earlier against later should be
lower than alternate for both. A shape correlation near the height's would mean people
have curves of their own and one published shape misses that.

Two changes after the first cohort run, both to statistics and neither to the records.
The correlation was first taken on the logs of the readings. A few people have almost
no insulin counted in a band and read up to 6e13 mg/dL per U, and those points set the
figure (0.36 against a rank correlation of 0.71 for the reading at 250 mg/dL, earlier
against later), so the rank correlation replaced it. The second column of factors first
compared a person with the median person inside one half, which is not a prediction of
anything, and became the miss described above.

The run. The height repeated as expected. The shape repeated weakly, also as expected,
and earlier against later was not lower than alternate for the shape, which was not
expected. On the second half, a person's own shape from the first half missed by more
than the median person's shape at 190 mg/dL and by the same at 250, under both ways of
halving, while the person's own height missed by less than the median height. On about
450 events a person, evaluating the shape per person bought nothing over one shape for
everyone. The two shape rows of the alternate table share 1.57 (1.43 to 1.79) in one
column, which was checked: the inputs differ and the figures part at the third decimal.

    python3 -m inv018.split_half            # from the repository root, DYNISF_ROOT unset

Writes results/inv018_split_half.md and charts/inv018/fig_split_half.png. The per-person
bins of each half go to results/inv018_split_half_bins.parquet first, which is not in git.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd

from inv009 import config, db
from inv012.peak_trough import PROMINENCE, jaeb_events, person_bins
from inv013.candidate import ANCHOR_GLUCOSE
from inv015.link import boot
from inv018.person_curves import (ANCHOR_BINS, BLUE, CHART, COHORT, HIGH, HIGH_BINS, INK, MUTED,
                                  SURFACE, band, curves)

MID, MID_BINS = 190, (180, 200)
SPLITS = ("Alternate events", "Earlier and later events")
HEIGHT, SHAPE = f"Reading at {ANCHOR_GLUCOSE:.0f} mg/dL", f"Share kept at {HIGH} mg/dL"
MEASURES = (HEIGHT, f"Reading at {HIGH} mg/dL", f"Share kept at {MID} mg/dL", SHAPE)
REPEATABLE, NOT_REPEATABLE = (0.7, 0.5), 0.5


def halves(n: int) -> dict:
    """For n events in time order, which fall in the first half under each way of halving."""
    i = np.arange(n)
    return dict(zip(SPLITS, (i % 2 == 0, i < n / 2)))


def _person(sid: str) -> pd.DataFrame | None:
    r = jaeb_events(sid, PROMINENCE)
    if r is None:
        return None
    ev, tdd, _ = r
    parts = [person_bins(ev[m if first else ~m], tdd).assign(split=split, half="ab"[not first])
             for split, m in halves(len(ev)).items() for first in (True, False)]
    return pd.concat(parts, ignore_index=True).assign(subject=sid, tdd=tdd, n_events=len(ev))


def measures(wide: pd.DataFrame) -> pd.DataFrame:
    """Per person: the height at the anchor, the reading at 250, and the shares kept."""
    a, hi = band(wide, *ANCHOR_BINS), band(wide, *HIGH_BINS)
    return pd.DataFrame(dict(zip(MEASURES, (a, hi, band(wide, *MID_BINS) / a, hi / a))))


def agreement(a: pd.Series, b: pd.Series) -> dict:
    """How well the first half predicts the second across people."""
    p = np.log(pd.DataFrame({"a": a, "b": b})).replace([np.inf, -np.inf], np.nan).dropna()
    def stat(x):
        own = np.median(np.abs(x.b - x.a))                    # predicted from the person's own first half
        shared = np.median(np.abs(x.b - x.a.median()))        # predicted from the median first half
        # Spearman by ranks, quicker than Series.corr. The last is above 1 where the own half wins.
        return np.array([np.corrcoef(x.a.rank(), x.b.rank())[0, 1],
                         np.exp(own), np.exp(shared), np.exp(shared - own)])
    r, own, shared, gain = boot(p, stat)
    return dict(people=len(p), r=r, full=tuple(2 * v / (1 + v) for v in r), own=own, shared=shared, gain=gain)


def verdict(full: tuple) -> str:
    if full[0] >= REPEATABLE[0] and full[1] > REPEATABLE[1]:
        return "repeatable"
    return "not repeatable" if full[2] < NOT_REPEATABLE else "undecided"


def compare(pb: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Agreement for every way of halving and every measure, and the pairs behind it."""
    rows, pairs = [], {}
    for split in SPLITS:
        a, b = (measures(curves(pb[(pb.split == split) & (pb.half == h)])) for h in "ab")
        for m in MEASURES:
            pairs[split, m] = pd.DataFrame({"a": a[m], "b": b[m]}).dropna()
            rows.append(dict(split=split, measure=m, **agreement(a[m], b[m])))
    return pd.DataFrame(rows), pairs


def figure(pairs: dict, table: pd.DataFrame, path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4), facecolor=SURFACE)
    for ax, m in zip(axes, (HEIGHT, SHAPE)):
        p = pairs[SPLITS[0], m]
        r = table[(table.split == SPLITS[0]) & (table.measure == m)].r.iloc[0]
        lim = np.percentile(p.to_numpy(), [1, 99]) * [0.8, 1.25]
        ax.plot(lim, lim, color=MUTED, lw=1, ls="--")
        ax.scatter(p.a, p.b, s=22, color=BLUE, alpha=0.6, edgecolor=SURFACE, lw=0.6)
        ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim, aspect="equal")
        ax.set_facecolor(SURFACE)
        ax.set_title(f"{m}, {len(p)} people\nrank correlation {r[0]:.2f} ({r[1]:.2f} to {r[2]:.2f})",
                     loc="left", fontsize=10.5, color=INK)
        ax.set_xlabel("From one half of the person's events", color=MUTED)
        ax.set_ylabel("From the other half", color=MUTED)
        ax.tick_params(colors=MUTED, length=0, which="both")
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter("{x:g}")
            axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(color="#e4e3df", lw=0.6, which="both")
        for side in ax.spines.values():
            side.set_visible(False)
    fig.suptitle("Each REPLACE-BG adult's interval curve from alternate halves of their events. "
                 "Dashed line: the halves agree.", x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


def write_md(pb: pd.DataFrame, table: pd.DataFrame, path) -> None:
    f = lambda v: f"{v[0]:.2f} ({v[1]:.2f} to {v[2]:.2f})"
    n = pb.drop_duplicates("subject").n_events
    L = ["# INV-018: does a person's own interval curve repeat?", "",
         "Generated by `inv018/split_half.py`. Each REPLACE-BG adult's events were split into two halves",
         "and a curve built from each half as in `inv018/person_curves.py`. The readings are outcome",
         "ratios (`results/inv012_peak_trough.md`), so what repeats is the measurement on that person.",
         "Events with no net insulin used are left out of every reading, 14% to 41% of events by glucose",
         "(`results/inv018_person_curves.md`).", "",
         f"{len(n)} adults, a median of {n.median():.0f} events each, so about {n.median() / 2:.0f} in a half.", "",
         "The rank correlation is across people, between the two halves. The reliability on all events",
         "steps it up by 2r / (1 + r) and is given for alternate events only. The last two columns are the",
         "typical factor by which a prediction of the second half missed: from the person's own first half,",
         "and from the median of everyone's first half. The last column divides the second miss by the first,",
         "so the person's own half predicts better where it is above 1.", "",
         "![Half against half](../charts/inv018/fig_split_half.png)", ""]
    for split in SPLITS:
        L += [f"## {split}", "",
              "| Quantity | People | Rank correlation between halves | Reliability on all events | Verdict | Miss from own first half | Miss from the median first half | Median miss over own miss |",
              "|---|---:|---:|---:|---|---:|---:|---:|"]
        for r in table[table.split == split].itertuples():
            alt = split == SPLITS[0]
            L.append(f"| {r.measure} | {r.people} | {f(r.r)} | {f(r.full) if alt else '-'} | "
                     f"{verdict(r.full) if alt else '-'} | {f(r.own)} | {f(r.shared)} | {f(r.gain)} |")
        L.append("")
    L += [f"A quantity is repeatable at a reliability of {REPEATABLE[0]} or more with an interval above {REPEATABLE[1]},",
          f"and not repeatable with an interval below {NOT_REPEATABLE}. The thresholds were set before the run.",
          "The rank correlation replaced a correlation of logs after the run, for the reason in the script."]
    path.write_text("\n".join(L) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="people, for a pilot")
    ap.add_argument("--workers", type=int, default=config.WORKERS)
    a = ap.parse_args()
    ids = db.subjects(COHORT).subject_id.tolist()
    ids = ids[:a.limit] if a.limit else ids
    t0 = time.time()
    with mp.Pool(a.workers, maxtasksperchild=8) as pool:
        res = [r for r in pool.imap_unordered(_person, ids, chunksize=1) if r is not None and len(r)]
    print(f"{len(ids)} people, {len(res)} with a usable bin, {time.time() - t0:.0f} s")
    tag = "_limit" if a.limit else ""
    pb = pd.concat(res, ignore_index=True).sort_values(["subject", "split", "half", "mid"])
    pb.to_parquet(config.RESULTS / f"inv018_split_half_bins{tag}.parquet", index=False)

    table, pairs = compare(pb)
    figure(pairs, table, CHART.with_name(f"fig_split_half{tag}.png"))
    out = config.RESULTS / f"inv018_split_half{tag}.md"
    write_md(pb, table, out)
    print(out.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
