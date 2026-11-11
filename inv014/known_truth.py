#!/usr/bin/env python3
"""The peak-to-trough method on simulated people whose sensitivity was measured first.

INV-012 ran the peak-to-trough ratio on a small simulator with a constant sensitivity
(inv009/synthetic.py). Its closed-loop arms returned a ratio that fell from 150 to 250
mg/dL. Its open-loop arm did not fall (250:200 of 1.31, 1.00 to 1.94), and no arm reached
100 or 300 mg/dL. This script repeats the check on TimSim, which runs the UVA/Padova
S2008 model through simglucose with layers for behaviour, sensor and infusion site.
TimSim's own audit (docs/SPEC_VS_OUTPUT.md) records that its validation figures predate
a rebuild, so none is quoted here and the simulator is not described as validated.

UVA/Padova has no single sensitivity: utilisation saturates in glucose and production
responds to both glucose and insulin. The true sensitivity at each glucose is therefore
measured before the method is read.

Truth, per person. Insulin is held at a constant rate for 40 h with no meals and every
layer off, which settles glucose at a level set by the rate. Two runs are made at each
rate. One continues unchanged. The other gives a bolus of 30 minutes of scheduled basal
(about 1% of the daily dose) at 40 h. The true sensitivity at that glucose is the largest
difference between the two glucose traces in the following 8 h, divided by the bolus
delivered. Twelve rates from 0.3 to 1.75 of scheduled basal give twelve levels, and the
sensitivity at 100, 150, 200, 250 and 300 mg/dL is interpolated between them with no
extrapolation.

Method, per person. The same people live 28 days under two regimes: therapy with no
algorithm (TimSim's "basal-bolus": scheduled basal, a bolus for announced carbohydrate,
a correction only with a meal) and oref0 with SMB. Each regime is run with the infusion
site layer on and off. `inv012.peak_trough.find_events` and `person_bins` are used
unchanged, and the comparison is on the median over events with positive insulin used,
as in INV-012. (Four things in this paragraph changed after a pilot and are listed below.)

Modelling choices, fixed before the run:
  TimSim at the commit recorded in the results file, its 40 adults with titrated
  settings (`examples.panel.dosable`), the bench seed.
  The truth is a fasting steady state. A free-living high usually has carbohydrate still
  absorbing and more insulin on board, so the truth describes the body at that glucose
  and not every state the body passes through at it. The physiology layer varies
  sensitivity by time of day and from day to day about a mean of one, so the truth is
  the centre of what the free-living runs had.
  A probe is dropped where the unchanged run moved by more than 5 mg/dL over the 6 h
  before the bolus, or where the run ended early.
  The record given to the method is what a real record holds: sensor glucose, pump
  insulin (including what a failed site did not absorb) net of scheduled basal, and
  carbohydrate as the person announced it. Unannounced meals and rescue carbohydrate
  are absent from it.
  Insulin on board uses the INV-009 default action curve, as in INV-012. The
  simulator's own absorption differs, as a real person's does.
  The daily dose is the insulin delivered per day in that run.
  The site-off arms were declared in advance because a failed site puts insulin into
  the record that never acted, and that lowers a fall-per-unit ratio at high glucose
  for a reason that is separate from the method.
  Arms are compared on the people who completed all four. Intervals are 95%, from
  resampling people 2,000 times. A bin needs 5 events in a person and 10 people.

Expected before the run. The probe was piloted on one person on 11 November 2026 and
returned 9.1, 13.8, 15.5 and 16.5 mg/dL per U at 93, 162, 226 and 287 mg/dL, so the true
sensitivity is expected to rise with glucose across the panel: below its 200 mg/dL
value at 100 and at or above it at 250 and 300. The method is expected to return more
than the truth x daily dose at 100 and 150 mg/dL, and a 100:200 ratio within person
above 1 where the truth's is below 1. From 200 to 300 mg/dL the method is expected to
stay level or fall while the truth rises. A method ratio that follows the truth's
within its interval at every pair of bins would count against the claim that the
method's shape is produced by the method.

Changed on 11 November 2026 after the four arms were piloted on one person, and before
any cohort run:
  Bins are 50 mg/dL wide and centred on 100, 150, 200, 250 and 300 mg/dL, where INV-012
  and the published method round the midpoint to 10 mg/dL. The pilot found about 120
  events in 28 days. The 10 mg/dL bins at 250 and 300 mg/dL held fewer than 5 of them
  under both regimes, and the bin at 200 mg/dL held none under oref0. The truth is read
  at the centre of each bin. Runs are 56 days.
  A regime of corrections by hand was added, with the site layer on and off: scheduled
  basal, the same meal bolus, and a full correction towards 140 mg/dL net of insulin on
  board whenever the sensor reads above 180 mg/dL, at most once in 3 h. The built-in
  therapy corrects only with a meal, so in the pilot more than half of its events had
  no net insulin used in nine of eleven bins, and the remaining ratios rested on doses
  near zero. The
  threshold, the target and the interval are choices and were not fitted to anything.
  The built-in therapy is kept with the site layer on, as the arm in which nobody doses
  between meals, and its site-off arm was dropped.
  Scheduled basal is rounded to the pump's increment before it is subtracted. The
  pilot person's schedule of 1.111 U/h is delivered as 1.10 U/h, which had put a
  standing negative dose into the record.
The expectations above were not changed, and they apply to the corrections-by-hand arms
as written.

    python -m inv014.known_truth --timsim ~/StudioProjects/TimSim

from the repository root, with an interpreter that has simglucose, and DYNISF_ROOT
unset. Writes results/inv014_known_truth.md. The per-person tables are parquet and stay
out of git, since TimSim does not redistribute anything derived from the published
parameter tables.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from inv009 import insulin_models as M
from inv009.forward_isf import MODEL_DEFAULT
from inv012.peak_trough import (BG_HI, BG_LO, MIN_PEOPLE, RATIOS, acted_cum, boot,
                                find_events, person_bins)

RESULTS = Path(__file__).resolve().parent.parent / "results"
POINTS = (100, 150, 200, 250, 300)
RATES = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.15, 1.3, 1.5, 1.75)   # of scheduled basal
HOLD_H, WATCH_H, SETTLE_H = 40, 8, 6
BOLUS_MIN_OF_BASAL = 30.0
MAX_DRIFT = 5.0            # mg/dL over SETTLE_H
PANEL, DAYS, WORKERS = 40, 56, 7
BIN_WIDTH = 50.0
HAND = dict(target=140.0, correct_above=180.0, min_interval_min=180.0)
ARMS = (("oref0 SMB", "oref", True), ("oref0 SMB, site off", "oref", False),
        ("corrections by hand", "hand", True), ("corrections by hand, site off", "hand", False),
        ("meal boluses only", "basal-bolus", True))


def _sim(timsim: str):
    """TimSim is found by path, in the worker, because macOS workers start empty."""
    if timsim not in sys.path:
        sys.path.insert(0, timsim)
    from simglucose.controller.base import Action
    from timsim.simulation.run import DEFAULT_LAYERS, simulate
    from timsim.therapy import titrate
    return Action, simulate, titrate, DEFAULT_LAYERS


class Hold:
    """A constant insulin rate that ignores the sensor, with one bolus at a set minute."""
    name = "hold"

    def __init__(self, action, rate_u_per_min, bolus_u, at_min):
        self.action, self.rate, self.bolus_u, self.at = action, rate_u_per_min, bolus_u, at_min

    def policy(self, observation, reward, done, **info):
        due = abs(info["clock_minutes"] - self.at) < 1e-6
        return self.action(basal=self.rate, bolus=self.bolus_u / info["sample_time"] if due else 0.0)


def probe(control: np.ndarray, dosed: np.ndarray, k: int, units: float, step: float) -> dict:
    """True sensitivity from two traces that differ by one bolus at sample `k`."""
    gap = control[k:] - dosed[k:]
    back = int(SETTLE_H * 60 / step)
    return dict(level=float(control[k]), sens=float(gap.max() / units),
                drift=float(control[k] - control[k - back]),
                at_h=float(np.argmax(gap) * step / 60.0))


def truth_at(level, sens, points=POINTS) -> np.ndarray:
    """Sensitivity at each glucose, interpolated inside the levels reached and NaN outside."""
    o = np.argsort(level)
    return np.interp(points, np.asarray(level, float)[o], np.asarray(sens, float)[o],
                     left=np.nan, right=np.nan)


def truth_job(job) -> dict:
    timsim, name, rate, seed = job
    Action, simulate, titrate, _ = _sim(timsim)
    prof = titrate.profile_for(name)
    basal, at = prof.basal_at(0), HOLD_H * 60.0
    out = dict(subject=name, rate=rate, error="")
    try:
        runs = [simulate(name, days=(HOLD_H + WATCH_H) / 24.0, seed=seed, layers=(), meals=[],
                         controller=Hold(Action, rate * basal, u, at), profile=prof)
                for u in (0.0, BOLUS_MIN_OF_BASAL * basal)]
    except Exception as e:      # a run that raised is reported apart from one that stopped short
        return {**out, "error": repr(e)[:200]}
    ctl, dosed = runs
    out["termination"] = " ".join(str(r.termination) for r in runs if r.termination)
    out["meals"] = len(ctl.meal_log) + len(dosed.meal_log) + ctl.rescue_events + dosed.rescue_events
    if out["termination"]:
        return out
    step, units = float(ctl.sample_time), float(dosed.bolus_series.sum())
    k = int(at / step)
    # the two runs received the same insulin up to the bolus, and the unchanged run no bolus
    assert ctl.bolus_series.sum() == 0 and units > 0, (name, rate, units)
    assert np.allclose(ctl.insulin_series[:k - 1], dosed.insulin_series[:k - 1]), (name, rate)
    return {**out, **probe(ctl.bg, dosed.bg, k, units, step), "units": units,
            "rate_u_h": float(ctl.insulin_series[:k - 1].mean() / step * 60.0)}


def arm_job(job):
    timsim, name, arm, driver, site, days, seed = job
    Action, simulate, titrate, default_layers = _sim(timsim)
    layers = tuple(l for l in default_layers if site or l != "site")
    prof = titrate.profile_for(name)
    ctrl = driver
    if driver == "oref":
        from timsim.controller import oref
        ctrl = oref.OrefController(name, profile=prof)
    if driver == "hand":       # the built-in therapy with a correction rule and basal left at schedule
        from timsim.controller.loop import LoopSettings
        inf = float("inf")
        ctrl = LoopSettings(**HAND, correction_fraction=1.0, max_correction_u=prof.limits.max_bolus_u,
                            max_iob_u=prof.limits.max_iob_u, predict_min=0.0, trend_window_min=0.0,
                            suspend_below=-inf, reduce_below=-inf, boost_above=inf,
                            max_basal_above=inf, max_basal_scale=1.0)
    row = dict(subject=name, arm=arm, error="", termination="")
    try:
        r = simulate(name, days=days, seed=seed, layers=layers, controller=ctrl, profile=prof)
    except Exception as e:
        return {**row, "error": repr(e)[:200]}, None
    finally:
        if hasattr(ctrl, "close"):
            ctrl.close()
    step, n = float(r.sample_time), len(r.cgm)
    cgm = np.asarray(r.cgm, float)
    bg = cgm.copy()
    bg[(bg < BG_LO) | (bg > BG_HI)] = np.nan
    from timsim.therapy.limits import _to_nearest
    inc = prof.limits.basal_increment_u_per_hour
    sched = np.array([_to_nearest(prof.basal_at(m - step) * 60.0, inc) for m in r.minutes]) / 60.0 * step
    carbs = np.zeros(n)
    for minute, _, announced in r.announced_log:
        i = int(np.searchsorted(r.minutes, minute))
        if announced > 0 and i < n:
            carbs[i] += announced
    tdd = float(r.insulin_units / (n * step / 1440.0))
    ev = find_events(r.minutes, bg, acted_cum(r.insulin_series - sched, M.kernel(MODEL_DEFAULT)), carbs)
    pct = np.percentile(cgm, [5, 25, 50, 75, 95])
    end = r.termination or {}      # a dict when the run ended early, which parquet will not take
    row.update(termination=str(end.get("direction", "")) + (", site not absorbing" if end.get("site_failing_at_end") else ""),
               completed=bool(r.completed), controller=r.controller,
               layers=",".join(r.layers), tdd=tdd, n_events=len(ev), days=n * step / 1440.0,
               announced=int((carbs > 0).sum()), unannounced=int(r.meals_unannounced),
               rescues=int(r.rescue_events),
               basal_net=float((r.insulin_series - sched - r.bolus_series).sum() / (n * step / 1440.0)),
               bolus=float(r.bolus_series.sum() / (n * step / 1440.0)), tar180=float(100 * np.mean(cgm > 180)),
               tar250=float(100 * np.mean(cgm > 250)), tbr70=float(100 * np.mean(cgm < 70)),
               **{f"p{q}": float(v) for q, v in zip((5, 25, 50, 75, 95), pct)})
    return row, person_bins(ev, tdd, BIN_WIDTH).assign(subject=name, arm=arm, tdd=tdd)


def _pool(fn, jobs, workers):
    t0 = time.time()
    with mp.Pool(workers) as pool:
        out = pool.map(fn, jobs, chunksize=1)
    print(f"{fn.__name__}: {len(jobs)} jobs in {time.time() - t0:.0f} s on {workers} workers", flush=True)
    return out


def _b(v, d=0) -> str:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return "-" if len(v) < MIN_PEOPLE else ("%.{0}f [%.{0}f, %.{0}f] (n=%d)".format(d) % (*boot(v), len(v)))


def report(truth: pd.DataFrame, runs: pd.DataFrame, pb: pd.DataFrame, stamp: str) -> str:
    L = ["# INV-014: the peak-to-trough method against a measured truth in TimSim", "",
         f"Generated by `inv014/known_truth.py`. {stamp}", "",
         "Medians over people with a 95% interval from resampling people. A dash marks fewer "
         f"than {MIN_PEOPLE} people. Bins are {BIN_WIDTH:.0f} mg/dL wide.", "", "## The true sensitivity by glucose", ""]
    bad = truth[(truth.error != "") | (truth.termination.fillna("") != "")]
    good = truth.drop(bad.index)
    settled = good[good.drift.abs() <= MAX_DRIFT]
    L += [f"{len(truth)} probes (two runs each): {int((truth.error != '').sum())} raised, "
          f"{len(bad) - int((truth.error != '').sum())} ended early, "
          f"{len(good) - len(settled)} had not settled, {len(settled)} used. "
          f"Meals or rescue carbohydrate in a probe run: {int(truth.meals.fillna(0).sum())}. "
          f"The largest difference came {settled.at_h.median():.1f} h after the bolus at the median, "
          f"and at the end of the {WATCH_H} h window in "
          f"{int((settled.at_h >= WATCH_H - 0.2).sum())} probes.", ""]
    tw = pd.DataFrame({s: truth_at(d.level, d.sens) for s, d in settled.groupby("subject")},
                      index=list(POINTS)).T
    head = "| | " + " | ".join(str(p) for p in POINTS) + " |"
    rule = "|---|" + "---:|" * len(POINTS)
    L += [head, rule,
          "| True sensitivity, mg/dL per U | " + " | ".join(_b(tw[p], 1) for p in POINTS) + " |",
          "| Relative to the person's own value at 200 mg/dL | "
          + " | ".join(_b(tw[p] / tw[200], 2) for p in POINTS) + " |", ""]

    L += ["## The runs", ""]
    ok = runs[(runs.error == "") & (runs.termination == "")]
    per = ok.groupby("subject").arm.nunique()
    panel = per[per == len(ARMS)].index
    L += [f"{len(panel)} people completed all {len(ARMS)} arms and are the panel for every table below.", "",
          "| Arm | Runs | Raised | Ended early | Controller and layers as run | Daily dose, U | "
          "Bolus, U/day | Basal away from schedule, U/day | Events | "
          "Sensor glucose 5th, 50th, 95th percentile, mg/dL | Above 180 | Above 250 | Below 70 |",
          "|---|---:|---:|---:|---|---:|---:|---:|---:|---|---:|---:|---:|"]
    for arm, d in runs.groupby("arm", sort=False):
        p = d[d.subject.isin(panel)]
        ends = d[d.termination != ""].termination.value_counts()
        L.append(f"| {arm} | {len(d)} | {int((d.error != '').sum())} | "
                 f"{', '.join(f'{k} {v}' for k, v in ends.items()) or 0} | "
                 f"{', '.join(sorted(set(p.controller)))}: {', '.join(sorted(set(p.layers)))} | "
                 f"{p.tdd.median():.1f} | {p.bolus.median():.1f} | {p.basal_net.median():.2f} | {int(p.n_events.sum())} | "
                 f"{p.p5.median():.0f}, {p.p50.median():.0f}, {p.p95.median():.0f} | "
                 f"{p.tar180.median():.1f}% | {p.tar250.median():.1f}% | {p.tbr70.median():.1f}% |")
    L.append("")

    long = tw.stack().rename("sens").rename_axis(["subject", "mid"]).reset_index()
    pb = pb[pb.subject.isin(panel) & pb.mid.isin(POINTS)].merge(long, on=["subject", "mid"], how="left")
    pb = pb.replace([np.inf, -np.inf], np.nan)
    pb["truth_tdd"] = pb.sens * pb.tdd
    pb["rel"] = pb.med_pos / pb.truth_tdd
    for arm, d in pb.groupby("arm", sort=False):
        w = {c: d.pivot(index="subject", columns="mid", values=c).reindex(columns=list(POINTS))
             for c in ("med_pos", "truth_tdd", "rel", "frac_nonpos", "dg", "used")}
        both = w["med_pos"].where(w["truth_tdd"].notna())
        L += [f"## {arm}", "", head, rule,
              "| Method: fall per unit x daily dose | " + " | ".join(_b(both[p]) for p in POINTS) + " |",
              "| Truth x daily dose, same people | "
              + " | ".join(_b(w["truth_tdd"][p].where(both[p].notna())) for p in POINTS) + " |",
              "| Method over truth, within person | " + " | ".join(_b(w["rel"][p], 2) for p in POINTS) + " |",
              "| Share of events with no net insulin used | "
              + " | ".join(_b(w["frac_nonpos"][p], 2) for p in POINTS) + " |",
              "| Fall, mg/dL | " + " | ".join(_b(w["dg"][p]) for p in POINTS) + " |",
              "| Net insulin used, U | " + " | ".join(_b(w["used"][p], 2) for p in POINTS) + " |", "",
              "| Within person | Method | Truth | Method over truth |", "|---|---:|---:|---:|"]
        for a, b in RATIOS:
            m, t = both[a] / both[b], w["truth_tdd"][a] / w["truth_tdd"][b]
            t = t.where(m.notna())
            L.append(f"| {a}:{b} | {_b(m, 2)} | {_b(t, 2)} | {_b(m / t, 2)} |")
        L.append("")
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--timsim", required=True, help="path to the TimSim checkout")
    ap.add_argument("--stage", choices=("truth", "arms", "report", "all"), default="all")
    ap.add_argument("--people", type=int, default=PANEL)
    ap.add_argument("--days", type=int, default=DAYS)
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--tag", default="", help="suffix for a pilot, so it cannot overwrite the run")
    a = ap.parse_args()
    timsim = str(Path(a.timsim).expanduser().resolve())
    _sim(timsim)
    from examples.panel import dosable
    from timsim import provenance, seeds
    names, seed = dosable(a.people), seeds.BENCH
    f = lambda kind: RESULTS / f"inv014_{kind}{a.tag}.parquet"

    if a.stage in ("truth", "all"):
        rows = _pool(truth_job, [(timsim, n, r, seed) for n in names for r in RATES], a.workers)
        pd.DataFrame(rows).to_parquet(f("truth"))
    if a.stage in ("arms", "all"):
        jobs = [(timsim, n, arm, drv, site, a.days, seed) for arm, drv, site in ARMS for n in names]
        out = _pool(arm_job, jobs, a.workers)
        pd.to_pickle(out, RESULTS / f"inv014_arms_raw{a.tag}.pkl")     # before anything that can fail
        pd.DataFrame([r for r, _ in out]).to_parquet(f("runs"))
        pd.concat([b for _, b in out if b is not None], ignore_index=True).to_parquet(f("person_bins"))
    if a.stage in ("report", "all"):
        s = provenance.stamp()
        dirty = subprocess.run(["git", "-C", timsim, "status", "--porcelain"], capture_output=True, text=True).stdout
        stamp = (f"TimSim commit {str(s['commit'])[:7]}{' (uncommitted changes)' if dirty.strip() else ''}, "
                 f"population {s['population']}, settings {s['settings']}, simglucose {s['simglucose']}, "
                 f"seed {seed}, {a.people} adults, {a.days} days.")
        text = report(pd.read_parquet(f("truth")), pd.read_parquet(f("runs")),
                      pd.read_parquet(f("person_bins")), stamp)
        (RESULTS / f"inv014_known_truth{a.tag}.md").write_text(text)
        print(text)


if __name__ == "__main__":
    main()
