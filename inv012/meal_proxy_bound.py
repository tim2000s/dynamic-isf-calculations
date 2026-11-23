#!/usr/bin/env python3
"""How much the meal-by-bolus substitution moves the peak-to-trough ratio.

In `inv012/peak_trough.py` an event needs no carbohydrate entry in the three
hours before its peak or during its fall. DCLP3, DCLP5, PEDAP and IOBP2 hold no
carbohydrate records, so there a bolus of at least 5% of the person's daily dose
stands in for a carbohydrate entry. Section 4.9 of the paper inferred that the
substitution probably lowers the within-person ratio of the value at 100 mg/dL to
the value at 200; the reviewer (round 4, point M3) thought it would raise it.
Loop and REPLACE-BG hold carbohydrate records, so in those two cohorts both rules
can be applied to the same people. This script runs the committed estimator there
twice: once as committed, with the carbohydrate records, and once with the
carbohydrate column emptied, which makes `grid_events` apply the bolus rule it
applies to the other four cohorts. Nothing else changes: the same grid, insulin
model, reference basal, daily dose and prominence.

Reported for each cohort and each rule: the pre-declared statistic (median over
every event) and the post hoc one (median over events with positive insulin used)
as the within-person 100-over-200 ratio, the share of events at 100 mg/dL with no
net insulin used, and the event counts. The effect of the substitution is the
bolus-rule cohort median minus the carbohydrate-rule cohort median, both taken
over the people who have the ratio under both rules, with a 95% percentile
interval from 2,000 bootstrap draws that resample those people. The median of each
person's own difference is given beside it.

Expected sign, written before the run. Two things change. The bolus rule removes
events holding a bolus of at least 5% of the daily dose, which includes large
corrections; those are commoner at high glucose and carry a large insulin used, so
removing them raises the value at 200 mg/dL and lowers the ratio. The bolus rule
also stops removing events near a carbohydrate entry that had no large bolus,
such as treatment of a low with carbohydrate alone and small snacks; such events
open near target with little insulin used and so have large or sign-unstable
values at 100 mg/dL, which raises the ratio under the pre-declared statistic and,
where insulin used is positive but small, under the post hoc one. The first effect
is expected to dominate under the post hoc statistic, so the ratio is expected to
fall with the bolus rule there, in line with the paper's inference. Under the
pre-declared statistic the direction is not predicted, because the share of
events with no net insulin used at 100 mg/dL is expected to rise and the median
over every event moves with that share.

The run of 18 November 2026 did not bear out the last expectation: under the bolus rule the
share of events at 100 mg/dL with no net insulin used fell, from 41% to 36% in REPLACE-BG and
from 40% to 33% in Loop (results/inv012_meal_proxy_bound.md).

The value is an outcome ratio, since an event ends where glucose stops falling.

    python3 -m inv012.meal_proxy_bound        # DYNISF_ROOT left unset

Outputs results/inv012_meal_proxy_bound.{json,md} (aggregates only).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd

from inv009 import db, grid as gridmod
from inv012.peak_trough import (MIN_PEOPLE, MODEL_DEFAULT, PROMINENCE, _loop_model, boot,
                                grid_events, person_bins, summarise)

STUDIES = ("ReplaceBG", "Loop")
RULES = ("carbohydrate records", "bolus of 5% of daily dose")
N_BOOT, SEED = 2000, 0
OUT = Path(__file__).resolve().parent.parent / "results"


def both_rules(job):
    subject_id, prom = job
    st = db.streams(subject_id)
    g = gridmod.build_grid(st)
    if g is None or g.empty:
        return None
    model = _loop_model(subject_id) if subject_id.startswith("Loop:") else MODEL_DEFAULT
    out = []
    for rule, grid in ((RULES[0], g), (RULES[1], g.assign(carbs_g=0.0))):
        ev, tdd, proxy = grid_events(grid, subject_id, st, model, prom)
        # A person with no carbohydrate entry at all already gets the bolus rule in the
        # committed run; they are kept, as committed, and counted in the output.
        if rule == RULES[1] and not proxy:
            raise RuntimeError(f"{subject_id}: the bolus rule did not reach grid_events")
        pb = person_bins(ev, tdd)
        if len(pb):
            out.append(pb.assign(subject=subject_id, cohort=f"{subject_id.split(':')[0]} | {rule}",
                                 tdd=tdd, n_events=len(ev), meal_proxy=proxy))
    return pd.concat(out, ignore_index=True) if out else None


def ratio_by_person(pb: pd.DataFrame, col: str) -> pd.Series:
    wide = pb.pivot(index="subject", columns="mid", values=col).replace([np.inf, -np.inf], np.nan)
    if 100 not in wide or 200 not in wide:
        return pd.Series(dtype=float)
    r = (wide[100] / wide[200]).replace([np.inf, -np.inf], np.nan).dropna()
    return r


def paired_effect(r_carb: pd.Series, r_bolus: pd.Series, n_boot: int = N_BOOT,
                  seed: int = SEED) -> dict | None:
    both = r_carb.index.intersection(r_bolus.index)
    if len(both) < MIN_PEOPLE:
        return None
    c, b = r_carb[both].to_numpy(float), r_bolus[both].to_numpy(float)
    idx = np.random.default_rng(seed).integers(0, len(both), (n_boot, len(both)))
    d = np.median(b[idx], axis=1) - np.median(c[idx], axis=1)
    pd_ = np.median((b - c)[idx], axis=1)
    lo, hi = np.percentile(d, [2.5, 97.5])
    plo, phi = np.percentile(pd_, [2.5, 97.5])
    verdict = ("interval above zero" if lo > 0 else "interval below zero" if hi < 0
               else "interval includes zero")
    return dict(n_people=int(len(both)), carb_median=float(np.median(c)),
                bolus_median=float(np.median(b)), diff=float(np.median(b) - np.median(c)),
                lo=float(lo), hi=float(hi), verdict=verdict,
                median_person_diff=float(np.median(b - c)), person_lo=float(plo),
                person_hi=float(phi))


def self_check() -> None:
    r1 = pd.Series([2.0, 3.0, 4.0] * 5, index=[f"p{i}" for i in range(15)])
    e = paired_effect(r1, r1 - 1.0, n_boot=100)
    assert abs(e["diff"] + 1.0) < 1e-12 and e["verdict"] == "interval below zero"


def markdown(res: dict) -> str:
    L = ["# Meal-by-bolus substitution in the peak-to-trough ratio, bounded in Loop and REPLACE-BG", "",
         "Generated by `inv012/meal_proxy_bound.py`. The committed estimator of `inv012/peak_trough.py` "
         "run on the same people with the carbohydrate records and with the bolus rule used in the four "
         "cohorts without carbohydrate records. Within-person ratio of the value at 100 mg/dL to the value "
         f"at 200, cohort median with a 95% interval from {N_BOOT} draws resampling people. Outcome "
         "ratios, not measured ISF.", ""]
    for s, r in res["cohorts"].items():
        L += [f"## {s}", "", "| Rule | People | Events | Pre-declared, all events: 100 over 200 (95% interval), people "
              "| Post hoc, positive insulin used: 100 over 200 (95% interval), people "
              "| Share of events at 100 mg/dL with no net insulin used |", "|---|---:|---:|---|---|---|"]
        for rule in RULES:
            x = r["rules"][rule]
            f = lambda k: ("%.2f (%.2f to %.2f), %d" % tuple(x[k])) if x.get(k) else "-"
            L.append(f"| {rule} | {x['n_people']} | {x['n_events']:,} | {f('med_all_100_over_200')} | "
                     f"{f('med_pos_100_over_200')} | {100 * x['frac_nonpos_100']:.0f}% |")
        L += ["", "Effect of the substitution (bolus rule minus carbohydrate rule), people with the ratio under both.", "",
              "| Statistic | People | Carbohydrate rule | Bolus rule | Difference of medians (95% interval) | Verdict "
              "| Median within-person difference (95% interval) |", "|---|---:|---:|---:|---|---|---|"]
        for k, lab in (("med_all", "Pre-declared, all events"), ("med_pos", "Post hoc, positive insulin used")):
            e = r["effect"].get(k)
            if e:
                L.append(f"| {lab} | {e['n_people']} | {e['carb_median']:.2f} | {e['bolus_median']:.2f} | "
                         f"{e['diff']:+.2f} ({e['lo']:+.2f} to {e['hi']:+.2f}) | {e['verdict']} | "
                         f"{e['median_person_diff']:+.2f} ({e['person_lo']:+.2f} to {e['person_hi']:+.2f}) |")
        L += ["", "Bin medians: fall (mg/dL) and insulin used (U) at 100 and 200 mg/dL.", "",
              "| Rule | Fall 100 | Fall 200 | Used 100 | Used 200 | Events per person, median |",
              "|---|---:|---:|---:|---:|---:|"]
        for rule in RULES:
            x = r["rules"][rule]
            L.append(f"| {rule} | {x['dg_100']:.0f} | {x['dg_200']:.0f} | {x['used_100']:.2f} | "
                     f"{x['used_200']:.2f} | {x['events_per_person']:.0f} |")
        L.append("")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--raw", type=Path, default=None, help="person-bin table, kept out of git")
    a = ap.parse_args()
    self_check()
    jobs = []
    for s in STUDIES:
        ids = db.subjects(s).subject_id.tolist()
        jobs += [(i, PROMINENCE) for i in (ids[:a.limit] if a.limit else ids)]
    t0 = time.time()
    with mp.Pool(min(a.workers, 7), maxtasksperchild=8) as pool:
        parts = [p for p in pool.imap_unordered(both_rules, jobs, chunksize=1) if p is not None]
    print(f"{len(jobs)} people, {len(parts)} with a usable bin, {time.time() - t0:.0f} s")
    pb = pd.concat(parts, ignore_index=True)
    if a.raw:
        pb.to_parquet(a.raw, index=False)
    summ = {s["cohort"]: s for s in summarise(pb)}
    res = {"prominence": PROMINENCE, "cohorts": {}}
    for s in STUDIES:
        r = {"rules": {}, "effect": {}}
        ratios = {}
        for rule in RULES:
            key = f"{s} | {rule}"
            x = summ[key]
            d = pb[pb.cohort == key]
            b100, b200 = x["bins"].get(100, {}), x["bins"].get(200, {})
            r["rules"][rule] = dict(
                n_people=x["n_people"], n_events=x["n_events"],
                med_all_100_over_200=x.get("med_all_100_over_200"),
                med_pos_100_over_200=x.get("med_pos_100_over_200"),
                frac_nonpos_100=b100.get("frac_nonpos", [np.nan])[0],
                dg_100=b100.get("dg", [np.nan])[0], dg_200=b200.get("dg", [np.nan])[0],
                used_100=b100.get("used", [np.nan])[0], used_200=b200.get("used", [np.nan])[0],
                events_per_person=float(d.drop_duplicates("subject").n_events.median()),
                people_already_on_bolus_rule=int(d.drop_duplicates("subject").meal_proxy.sum()))
            ratios[rule] = {c: ratio_by_person(d, c) for c in ("med_all", "med_pos")}
        for c in ("med_all", "med_pos"):
            r["effect"][c] = paired_effect(ratios[RULES[0]][c], ratios[RULES[1]][c])
        res["cohorts"][s] = r
    assert any(res["cohorts"][s]["rules"][RULES[0]]["n_events"] != res["cohorts"][s]["rules"][RULES[1]]["n_events"]
               for s in STUDIES), "the two rules selected identical events"
    a.out.mkdir(exist_ok=True)
    (a.out / "inv012_meal_proxy_bound.json").write_text(json.dumps(res, indent=1, default=float))
    md = markdown(res)
    (a.out / "inv012_meal_proxy_bound.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
