#!/usr/bin/env python3
"""Peak-to-trough fall per unit of insulin used, binned by midpoint glucose.

This runs the peak-to-trough estimator behind the published curve of sensitivity
against glucose (Romero-Ugalde et al., doi:10.2337/db25-920-P) on the records where
our own analyses found no glucose-dependent sensitivity: the JAEB cohorts of INV-009 and
the Trio and oref0 cohorts of INV-008. It also runs it on simulated records whose
sensitivity is constant by construction.

The estimator. An event starts at a glucose maximum away from a meal
and ends at the next minimum. The fall is G_max - G_min. The insulin used is
IOB(t1) - IOB(t2) plus delivery between the two, net of reference basal. Their
ratio is the event's ISF. Events are binned by midpoint glucose rounded to
10 mg/dL, and each person contributes the median of ISF x TDD in each bin.

Modelling choices, fixed before any run and before the abstract was read. Two differ
from the published rules, which ask for an hour clear of a meal with under 2.5 g of
carbohydrate on board and a minimum more than 20 minutes after the peak:
  A maximum or minimum is a turning point with a prominence of 15 mg/dL, since
  five-minute sensor data turns on noise far more often than on physiology.
  "Away from a meal" means no carbohydrate entry in the 3 h before the peak or
  during the fall.
  DCLP3, DCLP5, PEDAP and IOBP2 hold no carbohydrate records, so a bolus of 5%
  of the daily dose stands in for a meal there. That also removes events holding
  a large correction, so those four cohorts do not apply the meal rule and
  the output marks them.
  A person's bin needs 5 events.
  The primary value is the median over every event, with no event excluded.
  The median over events with positive insulin used, and the ratio of
  summed fall to summed insulin, are reported beside it.

Expected before running. The simulated arms have a true ISF x TDD of 2000 at
every glucose. An event is selected on glucose having fallen, so the numerator
is positive by construction while the denominator is small near target, and the
curve is expected to fall with glucose in the simulated arms as well as the real
ones. The simulator's pull toward 110 mg/dL adds fall at high glucose that no
insulin caused, which works in the opposite direction. A flat simulated curve at
2000 would mean the estimator recovers a constant and the real curves can be
read as sensitivity.

Changed after the first run, on 10 November 2026, and so not pre-specified.
The median over every event was declared primary and is still tabulated, but
the comparison of shape no longer rests on it. Near target the insulin
used is zero or negative in 9% to 58% of events, the ratio changes sign with
it, and the median lands on either side of zero (DCLP3 is negative at
100 mg/dL). The comparison uses the median over events with positive insulin
used. A line is also fitted through each person's events in a bin, fall against
insulin used, where the bin holds 10 events. Its intercept is the fall at no net
insulin and its slope is the fall per further unit.

Expected before the line was fitted. The intercept is positive, since a
qualifying fall is at least the prominence whatever the insulin. The slope is
positive. A sensitivity that falls with glucose as the published fit does puts the
slope x TDD at 100 mg/dL well above its value at 200. A constant sensitivity
puts the two about equal and leaves the intercept to produce the falling ratio.
The simulated arms should return a slope x TDD near 2000. Three things bias the
slope and none is corrected: error in the insulin estimate pulls it toward zero,
dosing that reacts to a stubborn glucose does the same, and a longer event has
both more fall and more insulin acted.

The simulated arms, corrected later the same day. Until then the simulator
wrote every unlogged meal into the carbohydrate record and absorbed it at twice
its size (inv009/tests/test_synthetic.py now holds this), so the meal rule
removed every excursion and the arms reached no bin above 130 mg/dL. With
unlogged meals out of the record the arms reach 250 mg/dL. Two gentler
controllers were added, at a tenth and a fifth of the original gain, because
the original spends a tenth of its time at the 40 mg/dL floor. The two gains
were chosen after looking at 0.002, 0.004 and 0.008 on 30 simulated people:
they bracket the insulin used per event at 200 mg/dL in the real cohorts. The
250:200 ratio was added because no simulated arm fills the 300 bin. Expected
before these arms were run through the full script: a ratio above 2000 that
falls with glucose, and a slope x TDD near 2000 at 150 to 250 mg/dL. A slope
that falls from 200 to 250 under this constant sensitivity would mean the
decline seen in REPLACE-BG and Loop needs no change in sensitivity to explain it.

    python3 -m inv012.peak_trough            # from the repository root

DYNISF_ROOT is left unset. It moves the INV-009 results directory as well, and
the per-person Loop insulin model is read from there.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from inv009 import config, db, grid as gridmod, insulin_models as M, synthetic as S
from inv009.forward_isf import MODEL_DEFAULT, _loop_model, schedule_units

PROMINENCE = 15.0          # mg/dL
CARB_CLEAR_MIN = 180.0
MAX_GAP_MIN = 10.0         # a longer step between readings ends the trace
BG_LO, BG_HI = 40.0, 400.0
# A rise in carbs on board this large is a carb entry. Fat and protein entered as
# carbohydrate equivalents (U001 in Trio) arrive as many small entries, some under it.
COB_JUMP_G = 2.0
MIN_PER_BIN = 5
MIN_PEOPLE = 10
FIT_MIN = 10               # events in a person's bin before a line is fitted
STUDIES = ("ReplaceBG", "Loop", "DCLP3", "DCLP5", "PEDAP", "IOBP2")
OREF = {"v5_trio": "oref_v5", "v7_oref0": "oref_v7"}
SYN_ISF, SYN_TDD, SYN_N, SYN_DAYS = 50.0, 40.0, 30, 200
SYN_ARMS = [("open loop, meals logged", dict()),
            ("open loop, 45% unlogged", dict(unannounced=0.45)),
            ("reactive, 45% unlogged", dict(unannounced=0.45, reactive=0.02)),
            ("reactive at a tenth the gain, 45% unlogged", dict(unannounced=0.45, reactive=0.002)),
            ("reactive at a fifth the gain, 45% unlogged", dict(unannounced=0.45, reactive=0.004))]
VALUES = ("med_all", "med_pos", "pooled", "slope_tdd")
# 100:150 because few oref people reach 200, 250:200 because no simulated arm reaches 300
RATIOS = ((100, 150), (100, 200), (250, 200), (300, 200))


def acted_cum(net: np.ndarray, kern: np.ndarray) -> np.ndarray:
    """Running insulin acted, so that acted over (i, j] is the difference.

    IOB(i) - IOB(j) + delivered in (i, j] rearranges to A[j] - A[i] with
    A = cumulative delivery - IOB. A dose in bin i is inside IOB(i) already.
    """
    return np.cumsum(net) - np.convolve(net, kern)[:len(net)]


def find_events(t, bg, acted, carb, prominence: float = PROMINENCE) -> pd.DataFrame:
    """Every maximum-to-next-minimum fall with no carbohydrate near it.

    `t` is minutes, `bg` is NaN where unusable, `acted` is the running insulin
    acted and `carb` is positive on rows where carbohydrate was entered.
    """
    ok = np.isfinite(bg) & np.isfinite(acted)
    brk = ~ok | (np.diff(t, prepend=t[0]) > MAX_GAP_MIN)
    ccum = np.concatenate([[0.0], np.cumsum(np.asarray(carb) > 0)])
    rows = []
    for idx in np.split(np.arange(len(t)), np.flatnonzero(brk)):
        idx = idx[ok[idx]]
        if len(idx) < 3:
            continue
        pk = find_peaks(bg[idx], prominence=prominence)[0]
        tr = find_peaks(-bg[idx], prominence=prominence)[0]
        for k, p in enumerate(pk):
            j = np.searchsorted(tr, p)
            if j == len(tr) or (k + 1 < len(pk) and tr[j] > pk[k + 1]):
                continue
            i1, i2 = idx[p], idx[tr[j]]
            if ccum[i2 + 1] - ccum[np.searchsorted(t, t[i1] - CARB_CLEAR_MIN)] > 0:
                continue
            rows.append((t[i2] - t[i1], bg[i1], bg[i2], acted[i2] - acted[i1]))
    return pd.DataFrame(rows, columns=["dur", "g_max", "g_min", "used"])


def _line(d: pd.DataFrame) -> pd.Series:
    """Fall = intercept + slope x insulin used, by least squares, in one bin."""
    if len(d) < FIT_MIN or np.ptp(d.used) == 0:
        return pd.Series({"icpt": np.nan, "slope": np.nan})
    slope, icpt = np.polyfit(d.used, d.dg, 1)
    return pd.Series({"icpt": icpt, "slope": slope})


def person_bins(ev: pd.DataFrame, tdd: float) -> pd.DataFrame:
    """One person's ISF x TDD by midpoint glucose, as three ratios and a fitted line."""
    if ev.empty:
        return ev
    ev = ev.assign(mid=(np.round((ev.g_max + ev.g_min) / 20.0) * 10).astype(int),
                   dg=ev.g_max - ev.g_min)
    with np.errstate(divide="ignore"):
        ev["isf_tdd"] = ev.dg / ev.used * tdd
    g = ev.groupby("mid")
    fit = g[["used", "dg"]].apply(_line)
    out = pd.DataFrame({
        "icpt": fit.icpt, "slope_tdd": fit.slope * tdd,
        "n": g.size(),
        "med_all": g.isf_tdd.median(),
        "med_pos": ev[ev.used > 0].groupby("mid").isf_tdd.median(),
        "pooled": g.dg.sum() / g.used.sum() * tdd,
        "frac_nonpos": g.used.apply(lambda u: float((u <= 0).mean())),
        "dg": g.dg.median(), "used": g.used.median(), "dur": g.dur.median(),
        "dg_nonpos": ev[ev.used <= 0].groupby("mid").dg.median()})
    return out[out.n >= MIN_PER_BIN].reset_index()


def _from_grid(g: pd.DataFrame, subject: str, cohort: str, streams: dict,
               model: str, prom: float, tdd: float | None = None) -> pd.DataFrame:
    n = len(g)
    tot = np.nan_to_num(g.total_u.to_numpy(float))
    sched, _ = schedule_units(g, subject, streams)
    if tdd is None:
        tdd = float(tot.sum() / max(n * config.GRID_MIN / 1440.0, 1e-9))
    carbs = np.nan_to_num(g.carbs_g.to_numpy(float))
    proxy = carbs.sum() <= 0
    if proxy:                 # no carbohydrate on record: a meal-sized bolus stands in
        carbs = np.nan_to_num(g.bolus_u.to_numpy(float)) >= config.MEAL_BOLUS_FRAC_TDD * tdd
    bg = g.cgm.to_numpy(float).copy()
    bg[(bg < BG_LO) | (bg > BG_HI)] = np.nan
    ev = find_events(np.arange(n) * config.GRID_MIN, bg,
                     acted_cum(tot - np.nan_to_num(sched), M.kernel(model)), carbs, prom)
    return person_bins(ev, tdd).assign(subject=subject, cohort=cohort, tdd=tdd,
                                       n_events=len(ev), meal_proxy=proxy)


def _jaeb(subject_id: str, prom: float) -> pd.DataFrame | None:
    st = db.streams(subject_id)
    g = gridmod.build_grid(st)
    if g is None or g.empty:
        return None
    study = subject_id.split(":")[0]
    return _from_grid(g, subject_id, study, st,
                      _loop_model(subject_id) if study == "Loop" else MODEL_DEFAULT, prom)


def _synth(seed: int, arm: str, kw: dict, prom: float) -> pd.DataFrame:
    g = S.simulate_subject(np.random.default_rng(seed), tdd=SYN_TDD, isf=SYN_ISF,
                           days=SYN_DAYS, **kw)
    # The nominal dose scales every bin alike, so the truth is 50 x 40 = 2000.
    return _from_grid(g, f"Synth:{arm}:{seed}", f"synthetic: {arm}", {},
                      MODEL_DEFAULT, prom, tdd=SYN_TDD)


def _oref(user_id: str, cohort: str, tdd: float, prom: float) -> pd.DataFrame | None:
    """The loop's own insulin activity, already net of profile basal, integrated.

    Activity is used in place of IOB difference plus delivery because the tables
    hold suggested rates and microboluses, which are not proof of delivery, and
    no manual boluses at all. The activity integral is the same quantity from the
    loop's own record of what it believed was on board.
    """
    table = OREF[cohort]
    cob = '"sug_COB"' if table == "oref_v5" else "sug_cob"
    d = db._q(f"SELECT ts_relative_sec, cgm_mgdl, iob_activity, {cob} AS cob FROM {table} "
              "WHERE user_id=%s ORDER BY ts_relative_sec", (user_id,))
    d = d.drop_duplicates("ts_relative_sec").apply(pd.to_numeric, errors="coerce")
    if len(d) < 1000:
        return None
    t = d.ts_relative_sec.to_numpy(float) / 60.0
    bg = d.cgm_mgdl.to_numpy(float).copy()
    if np.nanmedian(bg) < 40:
        bg = bg * 18.018
    bg[(bg < BG_LO) | (bg > BG_HI)] = np.nan
    act = d.iob_activity.to_numpy(float)              # units per minute
    a = np.nan_to_num(act)
    dt = np.clip(np.diff(t, prepend=t[0]), 0, MAX_GAP_MIN)
    acted = np.cumsum(0.5 * (a + np.roll(a, 1)) * dt)
    acted[~np.isfinite(act)] = np.nan
    carb = np.diff(d.cob.ffill().fillna(0).to_numpy(float), prepend=0.0) >= COB_JUMP_G
    ev = find_events(t, bg, acted, carb, prom)
    return person_bins(ev, tdd).assign(subject=user_id, cohort=cohort, tdd=tdd,
                                       n_events=len(ev), meal_proxy=False)


def _run(job):
    kind, args, prom = job
    return {"jaeb": _jaeb, "synth": _synth, "oref": _oref}[kind](*args, prom)


def boot(v, n: int = 2000) -> tuple[float, float, float]:
    """Median with a 95% interval from resampling people."""
    v = np.asarray(v, float)
    m = np.median(np.random.default_rng(0).choice(v, (n, len(v))), axis=1)
    return float(np.median(v)), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def summarise(pb: pd.DataFrame) -> list[dict]:
    out = []
    for cohort, d in pb.groupby("cohort", sort=False):
        row = dict(cohort=cohort, n_people=int(d.subject.nunique()),
                   n_events=int(d.drop_duplicates("subject").n_events.sum()),
                   meal_proxy=float(d.drop_duplicates("subject").meal_proxy.mean()), bins={})
        for col in VALUES + ("frac_nonpos", "dg", "dg_nonpos", "used", "dur", "icpt"):
            wide = d.pivot(index="subject", columns="mid", values=col)
            wide = wide.replace([np.inf, -np.inf], np.nan)
            for mid in wide.columns:
                v = wide[mid].dropna()
                if len(v) >= MIN_PEOPLE:
                    row["bins"].setdefault(int(mid), dict(n=int(len(v))))[col] = boot(v)
            if col in VALUES:      # each person's own ratio, so people are their own baseline
                for a, b in RATIOS:
                    if a in wide and b in wide:
                        r = (wide[a] / wide[b]).dropna()
                        r = r[np.isfinite(r)]
                        if len(r) >= MIN_PEOPLE:
                            row[f"{col}_{a}_over_{b}"] = boot(r) + (int(len(r)),)
        out.append(row)
    return out


def write_md(summary: list[dict], path) -> None:
    f = lambda b, c: "%.0f [%.0f, %.0f]" % b[c] if b and c in b else "-"
    r = lambda s, k: "%.2f [%.2f, %.2f] (n=%d)" % s[k] if k in s else "-"
    L = ["# Peak-to-trough ISF x TDD by midpoint glucose", "",
         "Generated by `inv012/peak_trough.py`. Cohort median of per-person bin medians,",
         "95% bootstrap interval over people.",
         "\"Meal by bolus\" is the share of people with no carbohydrate records, for whom a "
         "bolus of 5% of the daily dose stood in for a meal.", "",
         "The median over every event was declared primary before the first run. It is "
         "unstable near target, where the insulin used is zero or negative in a large share "
         "of events and the ratio changes sign, so the comparison with the published shape uses the "
         "median over events with positive insulin used. That change was made after seeing "
         "the first run. The fitted line was added at the same time.", ""]
    for col, title in zip(VALUES, ("Median over every event (declared primary before the first run)",
                                   "Median over events with positive insulin used "
                                   "(used for the comparison)",
                                   "Summed fall over summed insulin used",
                                   "Slope of fall on insulin used within a bin, x TDD")):
        L += [f"## {title}", "",
              "| Cohort | People | Meal by bolus | 100 | 150 | 200 | 250 | 300 "
              "| " + " | ".join(f"{a}:{b} within person" for a, b in RATIOS) + " |",
              "|" + "---|" * (8 + len(RATIOS))]
        for s in summary:
            L.append(f"| {s['cohort']} | {s['n_people']} | {100 * s['meal_proxy']:.0f}% | "
                     + " | ".join(f(s["bins"].get(m), col) for m in (100, 150, 200, 250, 300))
                     + " | " + " | ".join(r(s, f"{col}_{a}_over_{b}") for a, b in RATIOS) + " |")
        L.append("")
    L += ["## Fall at no net insulin used", "",
          f"Intercept (mg/dL) of the line whose slope is tabulated above, fitted per person in "
          f"bins holding {FIT_MIN} events.", "",
          "| Cohort | 100 | 150 | 200 | 250 | 300 |", "|---|---|---|---|---|---|"]
    for s in summary:
        L.append(f"| {s['cohort']} | "
                 + " | ".join(f(s["bins"].get(m), "icpt") for m in (100, 150, 200, 250, 300)) + " |")
    L += ["", "## What goes into the ratio", "",
          "Median fall (mg/dL), median insulin used (U), the share of events with no "
          "net insulin used, and the median fall in those events, at 100, 200 and 300 mg/dL.", "",
          "| Cohort | Fall 100 | Fall 200 | Fall 300 | Used 100 | Used 200 | Used 300 "
          "| Share <=0 at 100 | at 200 | at 300 | Fall with none used 100 | 200 |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    p = lambda b, c, fmt: fmt % b[c][0] if b and c in b else "-"
    for s in summary:
        b = [s["bins"].get(m) for m in (100, 200, 300)]
        L.append(f"| {s['cohort']} | " + " | ".join(
            [p(x, "dg", "%.0f") for x in b] + [p(x, "used", "%.2f") for x in b]
            + [p(x, "frac_nonpos", "%.2f") for x in b]
            + [p(x, "dg_nonpos", "%.0f") for x in b[:2]]) + " |")
    path.write_text("\n".join(L) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="people per cohort, for timing")
    ap.add_argument("--workers", type=int, default=config.WORKERS)
    ap.add_argument("--prominence", type=float, default=PROMINENCE)
    ap.add_argument("--oref-cohort", type=Path,
                    default=Path.home() / "oref-investigations" / "canonical_cohort.json")
    a = ap.parse_args()
    config.ensure_dirs()

    jobs = [("synth", (1000 * i + k, arm, kw))
            for i, (arm, kw) in enumerate(SYN_ARMS) for k in range(a.limit or SYN_N)]
    for s in STUDIES:
        ids = db.subjects(s).subject_id.tolist()
        jobs += [("jaeb", (sid,)) for sid in (ids[:a.limit] if a.limit else ids)]
    coh = json.loads(a.oref_cohort.read_text())
    for c in OREF:
        rs = [r for r in coh if r.get("cohort") == c and r.get("in_cohort") and r.get("tdd")]
        jobs += [("oref", (r["user_id"], c, float(r["tdd"])))
                 for r in (rs[:a.limit] if a.limit else rs)]

    jobs = [j + (a.prominence,) for j in jobs]
    t0 = time.time()
    with mp.Pool(a.workers, maxtasksperchild=8) as pool:
        res = [r for r in pool.imap_unordered(_run, jobs, chunksize=1)
               if r is not None and len(r)]
    print(f"{len(jobs)} people, {len(res)} with a usable bin, {time.time() - t0:.0f} s")
    pb = pd.concat(res, ignore_index=True)
    order = {c: i for i, c in enumerate([f"synthetic: {x}" for x, _ in SYN_ARMS]
                                        + list(STUDIES) + list(OREF))}
    pb = pb.sort_values(["cohort", "subject", "mid"],
                        key=lambda s: s.map(order) if s.name == "cohort" else s)
    tag = ("" if a.prominence == PROMINENCE else f"_prom{a.prominence:.0f}") \
        + ("_limit" if a.limit else "")
    pb.to_parquet(config.RESULTS / f"inv012_peak_trough_person_bins{tag}.parquet", index=False)
    summary = summarise(pb)
    (config.RESULTS / f"inv012_peak_trough{tag}.json").write_text(json.dumps(
        dict(prominence=a.prominence, carb_clear_min=CARB_CLEAR_MIN, min_per_bin=MIN_PER_BIN,
             synthetic_truth=SYN_ISF * SYN_TDD, cohorts=summary), indent=1))
    write_md(summary, config.RESULTS / f"inv012_peak_trough{tag}.md")
    print((config.RESULTS / f"inv012_peak_trough{tag}.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
