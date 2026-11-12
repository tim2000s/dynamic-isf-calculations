#!/usr/bin/env python3
"""Stage one of the proposal's evaluation: each arm replayed on recorded histories.

Every recorded loop decision of a person (a glucose reading and the daily dose the device
held at that moment) is given to each of the eight arms of
docs/dynamic-isf-8-proposal.md, and the insulin each arm asks for is set against the
person's profile sensitivity alone. Nothing is simulated. The histories were produced
under whatever each person ran at the time, so the stage shows what each arm would have
asked for at those readings and says nothing about where glucose would then have gone.
The stage is descriptive and has no pass rule.

What is reported is the insulin asked for per mg/dL of distance from 100 mg/dL, as a
multiple of what the profile sensitivity asks for. Above 100 mg/dL a value of 1.5 is a
correction half as large again as the profile's. Below 100 mg/dL it is the same multiple
of the insulin withheld. The total above 100 mg/dL is the sum of an arm's correction
units over every reading there, over the same sum under the profile.

Inputs. The per-reading glucose and daily dose come from the INV-008 replay cache
(inv008/stage2_replay.py): the device's own figure on Trio, a reconstruction on
AndroidAPS and oref0. Profile sensitivity is one value per person from
canonical_cohort.json, in mg/dL per U as that file holds it. The first run converted
values below 20 as mmol/L, the rule of inv008/phase11_shadow_eval.py. The only three such
values (10, 14 and 14) belong to people averaging 62 to 115 U a day, where the converted
figures of 180 to 252 mg/dL per U are not credible, so the rule was removed and the run
repeated. Carbohydrate on board is read from the decision tables to time meals.

Modelling choices, fixed before the run:
  Target 100 mg/dL for every person, the target of the proposal's worked tables.
  The reference arms read the sensitivity at the current glucose, as the systems do. The
  forms of the equation read it at the midpoint between the reading and target, except
  the arm named as read at the start.
  The logarithmic model as shipped is the AndroidAPS equation with divisor 55 for the
  AndroidAPS and oref0 cohorts (the insulin each person used is not recorded) and Trio's
  form at its defaults (factor 0.8, limits 0.7 and 1.2) for the Trio cohort. The sigmoid
  is Trio's at its defaults for every cohort.
  The average daily dose of a person is the mean of the per-reading daily dose. It serves
  the average-dose arm and is the reference dose of the first draft.
  A meal entry is a rise in carbohydrate on board of 5 g or more between consecutive
  decisions. A person is counted as announcing meals at one entry or more per recorded
  day, a recorded day being 288 decisions.
  Repeated decisions at one timestamp are collapsed to one. Readings from 40 to 400 mg/dL
  with a positive daily dose are kept. A person enters a row with 20 readings in it.
  Cohorts are reported separately. A cell is the median over people, with a 95% interval
  from resampling people 2,000 times.

Expected before the run, from the curve's shape (2.29, 1.30, 0.87, 0.70 and 0.60 of its
value at 180 mg/dL, at 100, 150, 200, 250 and 300 mg/dL):
  Profile link read at the midpoint: below 1 for readings up to 260 mg/dL, where the
  midpoint reaches 180, and above 1 beyond. About 0.5 from 100 to 150 mg/dL.
  Profile link read at the start: below 1 up to 180 mg/dL and above it beyond, so higher
  than the midpoint reading in every row above 100 mg/dL.
  Both links to the daily dose: the profile link's row multiplied by a person's profile
  sensitivity x daily dose over 1915, so the spread across people is wider and the level
  follows each cohort's product.
  Logarithmic as shipped: rises with glucose in every cohort, held between 0.7 and 1.2 on
  Trio.
  Sigmoid: 1 at 100 mg/dL, rising to 1.2 and falling to 0.7.
  First draft: near 1 above 100 mg/dL and below 1 under it.
  Time since a meal enters no arm. Within a glucose band it should change nothing, and
  the totals should be highest in the two hours after an entry, where glucose is highest.

Added on 12 November 2026, after the first run: a ninth arm, the profile link read at the
start with the shape held at its value at 200 mg/dL (0.87) for every reading above 200.
On the REPLACE-BG records the refitted shape follows the published one from 100 to 200
mg/dL and stops falling above (results/inv013_refit.md: 0.81 at 250 and 0.90 at 300
against the published 0.70 and 0.60). Expected before its run: equal to the start-read
arm up to 200 mg/dL, 1.15 in the two rows above, and a total between the midpoint arm's
0.70 and the start-read arm's 0.97 to 1.00, nearer the second. It returned 1.14 in those rows and
totals of 0.90, 0.90 and 0.93 in the three cohorts.

Added on 12 November 2026, after the ninth arm: two arms that take the shape from the
REPLACE-BG records over the whole range they cover, in place of the published shape cut at
200 mg/dL. The shape is the refit of results/inv013_refit.md (the published form with its
five coefficients free, fitted to the bins from 100 to 370 mg/dL), relative to its value at
180 mg/dL and read at the current glucose. The refit is lowest at 244 mg/dL, at 0.81, and
rises to 1.15 at 370. The bin medians tabulated in that file do not show the rise: they sit
between 0.81 and 0.88 from 210 to 290 mg/dL and scatter above 300, where 60 people or fewer
remain. The first arm therefore holds the shape at its lowest value for every reading above
244 mg/dL, which is a modelling choice. The second uses the refit as fitted up to 370 mg/dL.
Below 100 mg/dL both hold the value at 100 (2.16). Expected before the run: both within
0.05 of the start-read arm up to 200 mg/dL and about 1.2 from 200 to 260. Above 260 the
first should give 1.24 and the second less, near 1.1. The totals should lie between the
ninth arm's 0.90 to 0.93 and the start-read arm's 0.96 to 1.00, within 0.02 of each other.
The refitted quantity is an outcome ratio (inv013/refit.py), so only its shape is used and
its level is not read. The run agreed. From 100 to 150 mg/dL both gave 0.59 to 0.60 against
the start-read arm's 0.57 to 0.58, from 150 to 200 they gave 0.93 against 0.92, and from
200 to 260 they gave 1.19. Above 260 the first gave 1.24 and the second 1.13 to 1.15.
Totals were 0.93, 0.92 and 0.95 for the first and 0.92, 0.92 and 0.95 for the second.

The run agreed with each expectation but one: the totals were as high from two to four
hours after an entry as in the first two. Only 0, 2 and 5 people in the three cohorts fell
below one entry a day, so the stage describes people who announce meals.

Limits. One profile sensitivity per person stands in for a schedule that may vary by
hour. Every decision is counted, including those at which the loop would not have dosed.
The daily dose on AndroidAPS and oref0 is a reconstruction, and the two arms linked to it
and the logarithmic arm inherit its errors there. Nobody in the AndroidAPS cohort ran a
dynamic sensitivity (no decision of the 1,307,076 in `oref_v6` is flagged for it), so the
logarithmic arm there is the equation at 100% on profiles never set up for it, and the
stage holds no record of the model as AndroidAPS users run it. On Trio most people ran a dynamic sensitivity, so their glucose record already reflects
one of the reference arms.

Run as
  python3 -m inv016.replay --replay <dir of INV-008 replay parquets> --cohort <canonical_cohort.json>
to write results/inv016_replay.md. `--pilot DIR` runs two people over two days into DIR.
"""
from __future__ import annotations

import argparse
import glob
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd

from inv008.config import DB
from inv013.candidate import (RESULTS, TARGET, candidate, curve_isf, curve_shape,
                              stock_v1, trio_log, trio_sigmoid)

ARMS = ("Profile link, midpoint", "Profile link, start", "Profile link, start, flat above 200",
        "Profile link, start, REPLACE-BG shape, level above 244",
        "Profile link, start, REPLACE-BG shape as fitted",
        "TDD link", "Average-dose link", "Logarithmic as shipped", "Sigmoid", "First draft")
G_EDGES = (40, 80, 100, 150, 200, 260, 401)
G_LABELS = ("below 80", "80 to 100", "100 to 150", "150 to 200", "200 to 260", "above 260")
M_EDGES = (0, 2, 4, 6, np.inf)
M_LABELS = ("0 to 2 h", "2 to 4 h", "4 to 6 h", "over 6 h")
MEAL_RISE = 5.0
MIN_READINGS = 20
N_BOOT = 2000
WORKERS = 7
COHORTS = {"v5_trio": "Trio", "v6_aaps_classic": "AndroidAPS", "v7_oref0": "oref0"}
# the REPLACE-BG refit of results/inv013_refit.md, highest power first, as printed there
RBG_COEF = (1.161e-06, -0.001719, 0.9756, -236.4, 2.341e4)
RBG_MIN, RBG_TOP = 244.0, 370.0     # where the refit is lowest, and its last fitted bin


def rbg_shape(bg, top: float = RBG_MIN):
    """REPLACE-BG refit over its value at 180 mg/dL, held at its ends outside 100 to `top`."""
    return np.polyval(RBG_COEF, np.clip(bg, 100.0, top)) / np.polyval(RBG_COEF, 180.0)


def arms(bg, tdd, pisf: float, mean_tdd: float, trio: bool) -> dict:
    """Insulin asked for per mg/dL from target under each arm, over the profile's."""
    mid = (bg + TARGET) / 2.0
    log = trio_log(bg, tdd, pisf) if trio else stock_v1(bg, tdd, 55.0)
    return dict(zip(ARMS, (
        1.0 / curve_shape(mid), 1.0 / curve_shape(bg), 1.0 / curve_shape(np.minimum(bg, 200.0)),
        1.0 / rbg_shape(bg), 1.0 / rbg_shape(bg, RBG_TOP),
        pisf / curve_isf(mid, tdd), pisf / curve_isf(mid, mean_tdd),
        pisf / log, pisf / trio_sigmoid(bg, pisf), pisf / candidate(bg, pisf, tdd, mean_tdd))))


def hours_since_meal(ts, cob) -> tuple[np.ndarray, int]:
    """Hours since the last rise in carbohydrate on board, and the number of rises."""
    c = pd.Series(cob, dtype=float).ffill().to_numpy()
    rise = np.r_[False, np.diff(c) >= MEAL_RISE]
    # shortcut: an entry made during a gap in the record is timed at the first decision after it
    last = pd.Series(np.where(rise, ts, np.nan)).ffill().to_numpy()
    return (ts - last) / 3600.0, int(rise.sum())


def person(job: tuple) -> pd.DataFrame | None:
    """One person's sums by glucose band and time since a meal. Reads the database."""
    import psycopg2
    user, cohort, path, table, pisf, days = job
    d = (pd.read_parquet(path, columns=["ts_relative_sec", "bg", "tdd"])
         .drop_duplicates("ts_relative_sec").sort_values("ts_relative_sec"))
    col = '"sug_COB"' if table == "oref_v5" else "sug_cob"
    with psycopg2.connect(**DB) as conn, conn.cursor() as cur:
        cur.execute(f"SELECT ts_relative_sec, max({col}) FROM {table} WHERE user_id = %s GROUP BY 1", (user,))
        cob = pd.DataFrame(cur.fetchall(), columns=["ts_relative_sec", "cob"])
    d = d.merge(cob, on="ts_relative_sec", how="left")
    if days:
        d = d[d.ts_relative_sec < d.ts_relative_sec.iloc[0] + days * 86400]
    hrs, meals = hours_since_meal(d.ts_relative_sec.to_numpy(float), d.cob.to_numpy(float))
    ok = (d.bg.between(40, 400) & (d.tdd > 0)).to_numpy()
    if ok.sum() < MIN_READINGS:
        return None
    bg, tdd, hrs = d.bg.to_numpy(float)[ok], d.tdd.to_numpy(float)[ok], hrs[ok]
    mean_tdd, dist = float(tdd.mean()), np.abs(bg - TARGET)
    cells = pd.DataFrame({"g": np.digitize(bg, G_EDGES) - 1,
                          "m": np.where(np.isnan(hrs), -1, np.digitize(hrs, M_EDGES) - 1),
                          "n": 1, "dist": dist})
    for name, r in arms(bg, tdd, pisf, mean_tdd, cohort == "v5_trio").items():
        cells["r_" + name], cells["rd_" + name] = r, r * dist
    out = cells.groupby(["g", "m"], as_index=False).sum()
    return out.assign(user=user, cohort=cohort, product=pisf * mean_tdd,
                      meals_per_day=meals / (len(d) / 288.0))


def med_ci(x, rng) -> tuple[float, float, float]:
    """Median over people with a 95% interval from resampling people."""
    x = np.asarray(x, float)
    draws = np.median(x[rng.integers(0, len(x), (N_BOOT, len(x)))], axis=1)
    return float(np.median(x)), *np.percentile(draws, [2.5, 97.5])


def row(sub: pd.DataFrame, num: str, den: str, rng) -> tuple[int, list[str]]:
    """Per person, sum(num) / sum(den) for each arm; then the median over people."""
    p = sub.groupby("user").sum(numeric_only=True)
    p = p[p.n >= MIN_READINGS]
    if len(p) < 5:
        return len(p), ["-"] * len(ARMS)
    return len(p), ["{:.2f} ({:.2f} to {:.2f})".format(*med_ci(p[num + a] / p[den], rng)) for a in ARMS]


def report(cells: pd.DataFrame) -> str:
    rng = np.random.default_rng(0)
    head = "| People | " + " | ".join(ARMS) + " |"
    rule = "|---:|" + "---:|" * len(ARMS)
    L = ["# INV-016: the proposal's arms replayed on recorded histories", "",
         "Generated by `inv016/replay.py`, whose docstring holds the modelling choices and what was",
         "expected before the run. Each cell is the insulin an arm asks for per mg/dL of distance",
         "from 100 mg/dL, as a multiple of what the person's profile sensitivity asks for: the",
         "median over people with a 95% interval. The profile alone is 1.00 in every cell.",
         "The stage is descriptive and has no pass rule.", "",
         "## People", "",
         "| Cohort | People | Decisions | Announce meals | Profile sensitivity x average daily dose, median (10th to 90th centile) | Readings above 100 mg/dL |",
         "|---|---:|---:|---:|---:|---:|"]
    per = cells.groupby(["cohort", "user"]).agg(n=("n", "sum"), product=("product", "first"),
                                                meals=("meals_per_day", "first"))
    above = cells[cells.g >= 2].groupby(["cohort", "user"]).n.sum()
    for c, name in COHORTS.items():
        p = per.loc[c]
        L.append(f"| {name} | {len(p)} | {p.n.sum():,} | {(p.meals >= 1).sum()} | "
                 f"{p['product'].median():.0f} ({p['product'].quantile(.1):.0f} to {p['product'].quantile(.9):.0f}) | "
                 f"{100 * above.loc[c].sum() / p.n.sum():.0f}% |")
    cells = cells.assign(announces=cells.meals_per_day >= 1)
    for c, name in COHORTS.items():
        sub = cells[cells.cohort == c]
        L += ["", f"## {name}", "", "By glucose at the reading, mg/dL.", "", "| Glucose " + head, "|---" + rule]
        for i, lab in enumerate(G_LABELS):
            n, v = row(sub[sub.g == i], "r_", "n", rng)
            L.append(f"| {lab} | {n} | " + " | ".join(v) + " |")
        hi = sub[sub.g >= 2]
        L += ["", "Total correction units over every reading above 100 mg/dL.", "",
              "| Readings " + head, "|---" + rule]
        groups = [("All people", hi), ("Announce meals", hi[hi.announces]),
                  ("Do not announce meals", hi[~hi.announces])]
        groups += [(f"Announce meals, {lab} after an entry", hi[hi.announces & (hi.m == i)])
                   for i, lab in enumerate(M_LABELS)]
        for lab, g in groups:
            n, v = row(g, "rd_", "dist", rng)
            L.append(f"| {lab} | {n} | " + " | ".join(v) + " |")
        hi = hi.groupby("user").sum(numeric_only=True)
        L += ["", "The same total across people: 10th and 90th centile.", "",
              "| " + " | ".join(ARMS) + " |", "|" + "---:|" * len(ARMS),
              "| " + " | ".join("{:.2f} to {:.2f}".format(*(hi["rd_" + a] / hi.dist).quantile([.1, .9]))
                                for a in ARMS) + " |"]
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replay", type=Path, required=True)
    ap.add_argument("--cohort", type=Path, required=True)
    ap.add_argument("--pilot", type=Path, help="two people, two days, written to this directory")
    a = ap.parse_args()
    cohort = {r["user_id"]: r for r in json.loads(a.cohort.read_text())}
    jobs = []
    for f in sorted(glob.glob(str(a.replay / "*.parquet"))):
        r = cohort.get(Path(f).stem)
        if r and r["cohort"] in COHORTS:
            jobs.append((r["user_id"], r["cohort"], f, r["src_table"], float(r["isf"]),
                         2 if a.pilot else 0))
    out = RESULTS
    if a.pilot:
        jobs, out = [jobs[0], jobs[-1]], a.pilot
    t = time.time()
    first = person(jobs[0])       # one worker timed before the pool starts
    print(f"{len(jobs)} people; the first took {time.time() - t:.1f} s")
    with mp.Pool(WORKERS, maxtasksperchild=4) as pool:
        rest = pool.map(person, jobs[1:], chunksize=1)
    cells = pd.concat([x for x in [first, *rest] if x is not None], ignore_index=True)
    out.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(out / "inv016_replay_cells.parquet", index=False)   # raw, before reduction
    if a.pilot:
        print(cells.groupby("user").sum(numeric_only=True).T)
        return
    (out / "inv016_replay.md").write_text(report(cells))
    print((out / "inv016_replay.md").read_text())


if __name__ == "__main__":
    main()
