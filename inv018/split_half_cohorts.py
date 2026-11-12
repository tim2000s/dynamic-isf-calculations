#!/usr/bin/env python3
"""Does the split-half result hold outside REPLACE-BG, and does a longer record change it?

`inv018.split_half` found, in 196 REPLACE-BG adults with about 450 events each, that the
height of a person's interval curve repeated between halves of their events and its
shape did not, and that a person's own shape predicted their other half no better than
the median person's shape. It left two things unchecked, and this script checks them
with the same statistics.

Other cohorts. The same comparison is run in each cohort of `inv012.peak_trough`,
reported separately. In all but REPLACE-BG an algorithm doses, and its basal changes
are counted as insulin. What repeats there is the measurement on a person and their
algorithm together.

Longer records. Within a cohort, the people with at least k events are measured on
their first 200, 400, 800 and 1,600 events, k being the longest of those lengths that
30 people reach. The people are the same at every length, so a change along the row is
the length and not the people.

Modelling choices, fixed before the run:
  Everything in `inv018.split_half`: bins, bands, the five-event rule, the rank
  correlation, the two misses, the verdict thresholds, 2,000 resamples of people.
  Halves are alternate events throughout. Earlier against later is added on all events.
  A row is printed where 20 people have both halves.
  In DCLP3, DCLP5, PEDAP and IOBP2 no carbohydrate is recorded and a bolus of 5% of the
  daily dose stands in for a meal, as in INV-012.

Expected before the run.
  The height repeats in every cohort, rank correlation above 0.6.
  The shape repeats more under an algorithm than in REPLACE-BG, because the algorithm
  and the person's settings tie insulin to glucose the same way in both halves. A shape
  that repeats there would not show that people differ in physiology.
  Along a row of lengths the shape correlation rises. From 0.29 at about 450 events in
  REPLACE-BG, 2r / (1 + r) gives about 0.45 at 900 and 0.6 at 1,800 for a real
  difference between people diluted by noise. A flat row near zero would mean there is
  no difference to find.
  A person's own shape misses by less than the median shape only where the correlation
  between halves passes about 0.5. That is expected nowhere at 800 events or fewer.

The second miss divided by the first, with its own interval, was added after the first
cohort run, because two misses with overlapping intervals on the same people do not say
which is smaller. Nothing else changed.

The run. The height repeated in every cohort (0.74 to 0.95 between alternate halves) and
a person's own height predicted their other half better than the median height in every
cohort. The shape repeated more than in REPLACE-BG under Loop only (0.47 at 190 and 0.54
at 250 mg/dL), and there it fell to 0.28 and 0.27 between earlier and later events, so
part of what repeated belonged to the period. The other algorithm cohorts were near
REPLACE-BG. Along the lengths the shape correlation rose in Loop (0.45, 0.51, 0.57 at
190 mg/dL) and DCLP5 and did not in DCLP3 or PEDAP. A person's own shape did not predict
their other half better than the median shape in any cohort at any length: the ratio of
misses had no interval above 1, and the highest was 1.06 (0.93 to 1.20), Loop at 1,600
events. In REPLACE-BG, IOBP2 and PEDAP the own shape predicted worse at 190 mg/dL. The
Trio and oref0 cohorts had too few people with bins at 190 mg/dL in both halves to
enter the shape tables.

    python3 -m inv018.split_half_cohorts     # from the repository root, DYNISF_ROOT unset

Writes results/inv018_split_half_cohorts.md. The per-person bins of each half go to
results/inv018_split_half_cohorts_bins.parquet first, which is not in git.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd

from inv009 import config, db
from inv012.peak_trough import OREF, PROMINENCE, STUDIES, jaeb_events, oref_events, person_bins
from inv018.person_curves import curves
from inv018.split_half import HEIGHT, MEASURES, SHAPE, SPLITS, agreement, halves, measures, verdict

LENGTHS = (200, 400, 800, 1600)
ALL = 0                       # the `events` value for a person's whole record
MIN_PAIRS, MIN_TIER = 20, 30
SHARES = (MEASURES[2], SHAPE)


def _person(job) -> pd.DataFrame | None:
    kind, cohort, sid, tdd = job
    if kind == "jaeb":
        r = jaeb_events(sid, PROMINENCE)
        ev, tdd = (r[0], r[1]) if r else (None, None)
    else:
        ev = oref_events(sid, cohort, PROMINENCE)
    if ev is None or ev.empty:
        return None
    parts = []
    for k in (ALL,) + tuple(k for k in LENGTHS if k <= len(ev)):
        e = ev.iloc[:k] if k else ev
        for split, m in halves(len(e)).items():
            if k == ALL or split == SPLITS[0]:
                parts += [person_bins(e[h], tdd).assign(events=k, split=split, half=name)
                          for h, name in ((m, "a"), (~m, "b"))]
    return pd.concat(parts, ignore_index=True).assign(subject=sid, cohort=cohort, tdd=tdd, n_events=len(ev))


def tier(n_events: pd.Series) -> int:
    """The longest length that 30 people reach, or 0."""
    return max((k for k in LENGTHS if (n_events >= k).sum() >= MIN_TIER), default=0)


def rows(d: pd.DataFrame, split: str = SPLITS[0]) -> dict:
    """Agreement between the halves in `d` for each measure with enough people."""
    a, b = (measures(curves(d[(d.split == split) & (d.half == h)])) for h in "ab")
    both = lambda m: ((a[m] > 0) & (b[m] > 0) & np.isfinite(a[m] * b[m])).sum()   # index-aligned
    return {m: agreement(a[m], b[m]) for m in MEASURES if both(m) >= MIN_PAIRS}


def write_md(pb: pd.DataFrame, path) -> None:
    f = lambda v: f"{v[0]:.2f} ({v[1]:.2f} to {v[2]:.2f})"
    cohorts = list(dict.fromkeys(pb.cohort))
    L = ["# INV-018: the split-half comparison in every cohort, and at longer records", "",
         "Generated by `inv018/split_half_cohorts.py`, with the statistics of `inv018/split_half.py`.",
         "The readings are outcome ratios (`results/inv012_peak_trough.md`), and under an algorithm what",
         "repeats is the measurement on a person and their algorithm together. Events with no net insulin",
         "used are left out of every reading. Cohorts are reported separately. Intervals are 95% over people.",
         "A miss is the typical factor by which a prediction of the second half was wrong: from the person's",
         "own first half, or from the median of everyone's first half. The second divided by the first is",
         "above 1 where the person's own half predicts better. A cohort is absent from a table where",
         f"fewer than {MIN_PAIRS} people have the bins in both halves.", ""]
    whole = pb[pb.events == ALL]
    both = {c: [rows(whole[whole.cohort == c], split) for split in SPLITS] for c in cohorts}
    for m in (HEIGHT,) + SHARES:
        L += [f"## {m}, all of each person's events", "",
              "| Cohort | People | Events per person | Rank correlation, alternate halves | Reliability on all events | Verdict "
              "| Miss from own half | Miss from the median | Median miss over own miss "
              "| Rank correlation, earlier and later | Miss from own half | Miss from the median | Median miss over own miss |",
              "|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
        for c in cohorts:
            d = whole[whole.cohort == c]
            alt, late = (r.get(m) for r in both[c])
            if alt and late:
                L.append(f"| {c} | {alt['people']} | {d.drop_duplicates('subject').n_events.median():.0f} | {f(alt['r'])} | "
                         f"{f(alt['full'])} | {verdict(alt['full'])} | {f(alt['own'])} | {f(alt['shared'])} | {f(alt['gain'])} | "
                         f"{f(late['r'])} | {f(late['own'])} | {f(late['shared'])} | {f(late['gain'])} |")
        L.append("")
    L += ["## The same people at longer records", "",
          "In each cohort, the people whose record reaches the longest length that 30 people reach, measured",
          "on their first 200, 400, 800 and 1,600 events. Alternate halves, so each half holds half the events.",
          "The group is the same along a cohort's rows. Fewer of them have the bins at the shorter lengths.", ""]
    long = {}
    for c in cohorts:
        d = pb[pb.cohort == c]
        top = tier(d.drop_duplicates("subject").n_events)
        long[c] = {k: rows(d[(d.n_events >= top) & (d.events == k)]) for k in LENGTHS if k <= top}
    for m in (HEIGHT,) + SHARES:
        L += [f"### {m}", "",
              "| Cohort | Events used | People | Rank correlation | Miss from own half | Miss from the median | Median miss over own miss |",
              "|---|---:|---:|---:|---:|---:|---:|"]
        L += [f"| {c} | {k} | {r[m]['people']} | {f(r[m]['r'])} | {f(r[m]['own'])} | {f(r[m]['shared'])} | {f(r[m]['gain'])} |"
              for c in cohorts for k, r in long[c].items() if m in r]
        L.append("")
    path.write_text("\n".join(L) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="people per cohort, for a pilot")
    ap.add_argument("--workers", type=int, default=config.WORKERS)
    ap.add_argument("--oref-cohort", type=Path,
                    default=Path.home() / "oref-investigations" / "canonical_cohort.json")
    a = ap.parse_args()
    cut = lambda x: x[:a.limit] if a.limit else x
    jobs = [("jaeb", s, sid, None) for s in STUDIES for sid in cut(db.subjects(s).subject_id.tolist())]
    coh = json.loads(a.oref_cohort.read_text())
    for c in OREF:
        jobs += cut([("oref", c, r["user_id"], float(r["tdd"])) for r in coh
                     if r.get("cohort") == c and r.get("in_cohort") and r.get("tdd")])
    t0 = time.time()
    with mp.Pool(a.workers, maxtasksperchild=8) as pool:
        res = [r for r in pool.imap_unordered(_person, jobs, chunksize=1) if r is not None and len(r)]
    print(f"{len(jobs)} people, {len(res)} with a usable bin, {time.time() - t0:.0f} s")
    tag = "_limit" if a.limit else ""
    order = {c: i for i, c in enumerate(STUDIES + tuple(OREF))}
    pb = pd.concat(res, ignore_index=True)
    pb = pb.sort_values(["cohort", "subject", "events", "split", "half", "mid"],
                        key=lambda s: s.map(order) if s.name == "cohort" else s)
    pb.to_parquet(config.RESULTS / f"inv018_split_half_cohorts_bins{tag}.parquet", index=False)
    out = config.RESULTS / f"inv018_split_half_cohorts{tag}.md"
    write_md(pb, out)
    print(out.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
