#!/usr/bin/env python3
"""INV-019: glucose outcomes under each sensitivity rule, in simulation, with the harm rule.

This is the third stage of the evaluation in docs/dynamic-isf-8-proposal.md. Each arm is the
shipped oref0 algorithm run by TimSim's adapter, with the sensitivity it is handed each
five-minute cycle replaced by the arm's rule read at the current sensor glucose. The person,
the meals, the sensor noise and the titrated settings are the same in every arm (one seed),
so a difference between two arms is the rule. Everything below the line "Results" was
written after the run. Everything above it was committed before the run. A pilot on two
people for nine days ran before that commit, to exercise the code from the first cycle to
the report, and no arm, margin or expectation was changed after it.

Arms. Every arm with a rule runs autosens pinned at 1.0, which is what AndroidAPS does in
dynamic mode with its default settings (`ApsDynIsfAdjustSensitivity` off, read in
nightscout/AndroidAPS dev f66398c9d0).

  profile, pinned            the profile ISF unchanged. The reference for every arm.
  placebo                    the profile ISF x 1.01. The null arm: it changes nothing of
                             size, so what it shows is how far two runs drift apart alone.
  profile, autosens          oref0 as published, autosens detecting between 0.7 and 1.2.
  logarithmic, AndroidAPS    1800 / (TDD x ln(G / 55 + 1)), adjustment factor 100%, no limits.
  logarithmic, matched       the same with a factor set per person so that it returns the
                             profile ISF at 140 mg/dL at the person's titrated daily dose.
  logarithmic, Trio          factor 0.8, the ratio to the profile held in 0.7 to 1.2.
  Trio sigmoid               Trio's defaults, the profile ISF at the middle of the target.
  candidate                  profile ISF x Q(G) / Q(140), Q the published quartic.
  candidate, Trio limits     the candidate with its ratio to the profile held in 0.7 to 1.2.
  candidate, anchor 180      profile ISF x Q(G) / Q(180).
  equation on TDD            1915 x Q(G) / Q(140) / TDD.
  flat above 100             the proposal's first draft: the profile above 100, eased below.
  profile x 0.7, x 1.4       the profile ISF mis-set, pinned.
  candidate on x 0.7, x 1.4  the candidate on the mis-set profile.

Modelling choices, fixed before the run.
  TDD is AndroidAPS's blend, ((1.4 x last 4 h + 0.6 x the 4 h before) x 3) x 0.33 + seven-day
  daily mean x 0.34 + last 24 h x 0.33, from insulin the simulated pump delivered. AndroidAPS
  takes the previous calendar day where this takes the last 24 h. Trio's own blend is not
  reproduced, so the two logarithmic arms differ by factor and limits only.
  A rule switches on when seven days of delivery are held, as AndroidAPS requires, and this
  is applied to every rule, the ones that need no TDD included. Week one is therefore the
  same run in every pinned arm on the unscaled profile, and the report checks that it is.
  Outcomes are read from day 8 to day 35 on true glucose, which a compression artefact in
  the sensor cannot move. People are the 40 adults of INV-014, settings titrated by TimSim.
  A person counts in a comparison when both runs completed. Runs that ended early are
  counted by direction and judged separately.
  Intervals are 95%, from resampling people 2,000 times, on the mean paired difference.

The harm rule, with two sides. It is applied to the candidate against two comparators:
`logarithmic, matched`, which has the candidate's level at 140 mg/dL and differs in shape,
and `profile, autosens`, which is what a person has with dynamic sensitivity switched off.
  Low side: the difference in percentage of time below 54 mg/dL, candidate less comparator.
  High side: the same for time above 250 mg/dL.
  Margins: 0.25 points below 54 and 1.25 points above 250, a quarter of the consensus limits
  of 1% and 5% (Battelino et al., Diabetes Care 2019;42:1593). The fraction is a choice.
  A side reads harm when the interval lies wholly above zero and the mean is above the
  margin. It reads cleared when the upper bound is below the margin. Anything else is
  unresolved.
  Early ends: a side also reads harm when the people whose candidate run ended early in that
  direction while the comparator's completed outnumber the reverse by three or more.
  Null check: each criterion is computed first for the placebo against `profile, pinned`.
  The rule is used only where the placebo reads cleared. Where it does not, the margin sits
  inside the drift of the design, the side is reported as undecidable at this panel size,
  and the margins are left where they are.
  What follows from a verdict. A harm reading on either side is reported with the arm's
  glucose percentiles and looked into, and it does not by itself stop the candidate. Tim
  decided this on 12 November 2026, after the two-person pilot and before the full run,
  because the simulated body departs from the people the rule is written for in the ways
  the next paragraph sets out. The first version of this paragraph, committed with the
  plan, made harm a stop before any code that doses. Cleared on both sides against both
  comparators supports starting the fourth stage. Unresolved leads to a larger panel or a
  second seed. Code that calculates the value in shadow and doses nothing is unaffected
  by any reading.
  `logarithmic, AndroidAPS` is reported beside these and is outside the rule, for the reason
  in the next paragraph.

What this panel is, from results/inv019_panel.md (written by `--stage panel`, which runs no
simulation) and results/inv014_known_truth.md. The titrated profile ISF sits near the true
sensitivity TimSim's body shows at 150 mg/dL. Profile ISF x daily dose is about half of the
1,915 that adults enter, so the logarithmic model at a factor of 100% doses these people
more weakly than their profile at every glucose for the median person. In the worked person
of the proposal it doses more strongly than the profile below 211 mg/dL. The matched arm
removes that difference of level. The body differs from the equation's premise as well: its
true sensitivity rises from 100 to 300 mg/dL, and the equation's falls.

Expected before the run, each against `profile, pinned` unless stated.
  placebo: every interval includes zero.
  profile, autosens: time below 70 mg/dL lower or unchanged.
  candidate: it asks for more insulin above 140 mg/dL and less below, in a body where a unit
  does more at high glucose. Time above 250 lower. Time below 54 higher, and the low
  side of the harm rule is expected to read harm or unresolved, with oref0's ceilings
  (insulin on board, microbolus size) deciding which. The sign of the change in mean glucose is not
  predicted.
  candidate, Trio limits: the same signs, smaller, since the ratio stops at 1.2.
  candidate, anchor 180: sensitivity 1.44 times the candidate's at every glucose, so less
  insulin than the candidate throughout: time below 54 lower than the candidate's and time
  above 250 higher.
  candidate against logarithmic, matched: equal at 140, the candidate steeper on both
  sides. Time above 250 lower and time below 54 higher.
  logarithmic, AndroidAPS: less insulin than the profile for the median person. Mean glucose
  and time above 250 higher, time below 54 lower.
  profile x 0.7 doses 1.43 times the insulin per correction and x 1.4 doses 0.71 times.
  Time below 54 rises under x 0.7 and time above 250 under x 1.4, and the candidate carries
  each error through unchanged, so its differences from its own mis-set profile keep the
  signs given for the candidate.

Limits. The vendored oref0 takes one sensitivity per cycle. AndroidAPS in dynamic mode also
uses the value along its forecast, and that is not exercised. A harm verdict here describes
the equation in a body built the other way round and bounded by oref0's ceilings. TimSim
orders options and finds failures, and its own documents say an absence of events on a
bench of this size does not show safety. The unannounced-meal arm, a mis-set adjustment
factor and the plain average-dose arm of the proposal are left for a second run.

    python -m inv019.outcome_test --timsim ~/StudioProjects/TimSim

from the repository root, with an interpreter that has simglucose, and DYNISF_ROOT unset.
Each run is written to results/inv019_raw/ before anything is reduced, and a run whose file
is already there is not repeated, so an interrupted run resumes and a fresh one starts from
an empty directory. That directory and
every per-person table stay out of git, since TimSim does not redistribute anything derived
from the published parameter tables. Writes results/inv019_outcome_test.md.

Results
(none yet)
"""
from __future__ import annotations

import argparse
import subprocess
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from inv008.dynisf import quartic_isf
from inv012.peak_trough import MIN_PEOPLE
from inv013.candidate import (ANCHOR_GLUCOSE, REF_GLUCOSE, candidate, curve_isf, curve_profile,
                              curve_shape, stock_v1, trio_log, trio_sigmoid)
from inv014.known_truth import _pool, _sim

RESULTS = Path(__file__).resolve().parent.parent / "results"
PANEL, DAYS, WARM_DAYS, WORKERS = 40, 35, 7, 7
DIVISOR = 55.0                       # AndroidAPS, insulin peaking at 75 min
LIMITS = (0.7, 1.2)                  # Trio's default autosens limits
MARGIN = {"tbr54": 0.25, "tar250": 1.25}        # points of time: a quarter of 1% and of 5%
EARLY_EXCESS = 3
ENDS = {"tbr54": "hypo", "tar250": "hyper"}
N_BOOT = 2000

# rule(glucose, profile ISF, blended TDD, matched factor, target) -> mg/dL per U
RULES = {
    "profile": lambda g, isf, tdd, f, t: isf,
    "log": lambda g, isf, tdd, f, t: stock_v1(g, tdd, DIVISOR),
    "log_matched": lambda g, isf, tdd, f, t: stock_v1(g, tdd, DIVISOR, f),
    "log_trio": lambda g, isf, tdd, f, t: trio_log(g, tdd, isf, DIVISOR, 0.8, *LIMITS),
    "sigmoid": lambda g, isf, tdd, f, t: trio_sigmoid(g, isf, t),
    "candidate": lambda g, isf, tdd, f, t: curve_profile(g, isf),
    "candidate_limits": lambda g, isf, tdd, f, t: isf / np.clip(1.0 / curve_shape(g), *LIMITS),
    "candidate_180": lambda g, isf, tdd, f, t: curve_profile(g, isf) * quartic_isf(ANCHOR_GLUCOSE) / quartic_isf(REF_GLUCOSE),
    "curve_tdd": lambda g, isf, tdd, f, t: curve_isf(g, tdd),
    "flat": lambda g, isf, tdd, f, t: candidate(g, isf, 1.0, 1.0),
}
# arm: (rule, multiplier on the profile ISF, autosens: a number pins it and None detects)
ARMS = {
    "profile, pinned": ("profile", 1.0, 1.0),
    "placebo": ("profile", 1.01, 1.0),
    "profile, autosens": ("profile", 1.0, None),
    "logarithmic, AndroidAPS": ("log", 1.0, 1.0),
    "logarithmic, matched": ("log_matched", 1.0, 1.0),
    "logarithmic, Trio": ("log_trio", 1.0, 1.0),
    "Trio sigmoid": ("sigmoid", 1.0, 1.0),
    "candidate": ("candidate", 1.0, 1.0),
    "candidate, Trio limits": ("candidate_limits", 1.0, 1.0),
    "candidate, anchor 180": ("candidate_180", 1.0, 1.0),
    "equation on TDD": ("curve_tdd", 1.0, 1.0),
    "flat above 100": ("flat", 1.0, 1.0),
    "profile x 0.7": ("profile", 0.7, 1.0),
    "candidate on x 0.7": ("candidate", 0.7, 1.0),
    "profile x 1.4": ("profile", 1.4, 1.0),
    "candidate on x 1.4": ("candidate", 1.4, 1.0),
}
NULL, REFERENCE = "placebo", "profile, pinned"
REFERENCES = {"candidate on x 0.7": "profile x 0.7", "candidate on x 1.4": "profile x 1.4"}
COMPARATORS = ("logarithmic, matched", "profile, autosens")
METRICS = ("tbr54", "tar250", "tbr70", "tar180", "tir", "mean", "dose")


def blended_tdd(cum, step: float = 5.0):
    """AndroidAPS's blended daily dose from cumulative units per step. None before seven days."""
    if len(cum) <= 7 * 1440 / step:
        return None
    back = lambda hours: cum[-1] - cum[-1 - int(hours * 60 / step)]
    return (1.4 * back(4) + 0.6 * (back(8) - back(4))) * 3 * 0.33 + back(168) / 7 * 0.34 + back(24) * 0.33


def matched_factor(isf: float, tdd: float) -> float:
    """The adjustment factor at which the logarithmic model returns `isf` at the anchor."""
    return 1800.0 / (tdd * np.log(ANCHOR_GLUCOSE / DIVISOR + 1.0) * isf)


def outcomes(v) -> dict:
    """Percentages of time and the mean, for one glucose series in mg/dL."""
    v = np.asarray(v, float)
    p5, p50, p95 = (float(x) for x in np.percentile(v, [5, 50, 95]))
    return dict(p5=p5, p50=p50, p95=p95, mean=float(v.mean()), tir=float(100 * np.mean((v >= 70) & (v <= 180))),
                tbr70=float(100 * np.mean(v < 70)), tbr54=float(100 * np.mean(v < 54)),
                tar180=float(100 * np.mean(v > 180)), tar250=float(100 * np.mean(v > 250)))


def mean_ci(d) -> tuple[float, float, float]:
    """Mean with a 95% interval from resampling people."""
    d = np.asarray(d, float)
    m = np.random.default_rng(0).choice(d, (N_BOOT, len(d))).mean(axis=1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def verdict(d, margin: float) -> str:
    """Harm, cleared or unresolved for paired differences, by the rule in the docstring."""
    m, lo, hi = mean_ci(d)
    return "harm" if lo > 0 and m > margin else "cleared" if hi < margin else "unresolved"


def controller(timsim: str, name: str, base, rule: str, autosens):
    """oref0 on `base`, a profile, with each cycle's sensitivity set by `rule`."""
    _sim(timsim)
    from timsim.controller.oref import OrefController
    factor = matched_factor(float(base.isf[0]), float(base.total_daily_insulin))

    class Arm(OrefController):
        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            self.cum, self.sens_log = [0.0], []

        def policy(self, observation, reward, done, **info):
            given = sum(info.get(k) or 0.0 for k in ("delivered_basal_u", "delivered_bolus_u"))
            self.cum.append(self.cum[-1] + given)
            minute = float(info.get("clock_minutes", 0.0))
            isf, tdd = base.isf_at(minute), blended_tdd(self.cum, float(info.get("sample_time", 5.0)))
            sens = isf if tdd is None else float(RULES[rule](float(observation.CGM), isf, tdd, factor, base.target))
            assert np.isfinite(sens) and sens > 0, (name, rule, sens)
            self.profile = replace(base, isf=np.full(24, sens))
            self.sens_log.append(sens)
            return super().policy(observation, reward, done, **info)

    return Arm(name, profile=base, autosens=autosens)


def arm_job(job) -> dict:
    timsim, name, arm, days, seed, raw = job
    _, simulate, titrate, _ = _sim(timsim)
    rule, scale, autosens = ARMS[arm]
    prof = titrate.profile_for(name)
    prof = replace(prof, isf=prof.isf * scale)
    ctrl = controller(timsim, name, prof, rule, autosens)
    row, out = dict(subject=name, arm=arm, error=""), Path(raw) / f"{arm}__{name}.npz"
    if out.exists():
        ctrl.close()
        return row
    try:
        r = simulate(name, days=days, seed=seed, controller=ctrl, profile=prof)
    except Exception as e:      # a run that raised is reported apart from one that stopped short
        return {**row, "error": repr(e)[:200]}
    finally:
        ctrl.close()
    end = r.termination or {}      # a dict when the run ended early
    np.savez_compressed(
        out, subject=name, arm=arm, days=days, seed=seed,
        completed=bool(r.completed), termination=str(end.get("direction", "")),
        rescues=int(r.rescue_events), layers=",".join(r.layers), step=float(r.sample_time),
        minutes=np.asarray(r.minutes, float), bg=np.asarray(r.bg, np.float32),
        cgm=np.asarray(r.cgm, np.float32), insulin=np.asarray(r.insulin_series, np.float32),
        sens=np.asarray(ctrl.sens_log, np.float32), isf=float(prof.isf[0]))
    return row


def reduce(path: Path) -> dict:
    """One run's raw file to one row: week one for the checks, the rest for the outcomes."""
    z = np.load(path)
    warm = z["minutes"] >= WARM_DAYS * 1440.0
    per_day = lambda x: float(x.sum() / (len(x) * float(z["step"]) / 1440.0)) if len(x) else np.nan
    row = {k: z[k].item() for k in ("subject", "arm", "days", "seed", "completed", "termination", "rescues", "layers")}
    live = z["sens"][int(WARM_DAYS * 1440 / float(z["step"])) + 1:]
    row.update(week1_bg=float(z["bg"][~warm].mean()), week1_dose=per_day(z["insulin"][~warm]),
               dose=per_day(z["insulin"][warm]),
               sens_ratio=float(np.median(live) / z["isf"]) if len(live) else np.nan,
               sens_moved=float(np.mean(~np.isclose(live, z["isf"], rtol=1e-4))) if len(live) else np.nan)
    if warm.any():
        row.update(outcomes(z["bg"][warm]), **{f"cgm_{k}": v for k, v in outcomes(z["cgm"][warm]).items()})
    return row


def _f(v, d=2) -> str:
    return f"{v[0]:+.{d}f} ({v[1]:+.{d}f} to {v[2]:+.{d}f})"


def paired(runs: pd.DataFrame, arm: str, ref: str, metric: str) -> np.ndarray:
    """Arm less reference for the people who completed both."""
    w = runs[runs.completed].pivot(index="subject", columns="arm", values=metric)
    return (w[arm] - w[ref]).dropna().to_numpy() if {arm, ref} <= set(w.columns) else np.array([])


def early_excess(runs: pd.DataFrame, arm: str, ref: str, direction: str) -> int:
    """People whose `arm` run ended early in `direction` while `ref` completed, less the reverse."""
    ended = runs.assign(e=runs.termination.str.startswith(direction)).pivot(index="subject", columns="arm", values="e")
    done = runs.pivot(index="subject", columns="arm", values="completed")
    return int((ended[arm] & done[ref]).sum() - (ended[ref] & done[arm]).sum())


def harm_rows(runs: pd.DataFrame, arm: str, ref: str) -> list[str]:
    out = []
    for metric, margin in MARGIN.items():
        d, extra = paired(runs, arm, ref, metric), early_excess(runs, arm, ref, ENDS[metric])
        word = "-" if len(d) < MIN_PEOPLE else "harm" if extra >= EARLY_EXCESS else verdict(d, margin)
        out.append(f"| {arm} | {ref} | {'below 54' if metric == 'tbr54' else 'above 250'} | {len(d)} | "
                   f"{_f(mean_ci(d)) if len(d) >= MIN_PEOPLE else '-'} | {margin} | {extra:+d} | {word} |")
    return out


def report(runs: pd.DataFrame, raised: int, stamp: str) -> str:
    assert runs.days.nunique() == 1 and runs.seed.nunique() == 1, "rows disagree on the design"
    L = ["# INV-019: glucose outcomes under each sensitivity rule, in TimSim", "",
         f"Generated by `inv019/outcome_test.py`. {stamp}", "",
         f"True glucose from day {WARM_DAYS + 1}. Differences are means over people who completed both "
         f"runs, with a 95% interval from resampling people. A dash marks fewer than {MIN_PEOPLE} people. "
         f"{raised} runs raised.", "", "## Checks on the runs", "",
         "Week one is the same run in every pinned arm on the unscaled profile, so its largest "
         "difference from the reference should be zero. The applied sensitivity over the profile "
         "ISF shows that each rule reached its runs.", "",
         "| Arm | Runs | Largest week-one difference in mean glucose, mg/dL | Median applied sensitivity over profile ISF | Cycles where it differed from the profile |",
         "|---|---:|---:|---:|---:|"]
    w1 = runs.pivot(index="subject", columns="arm", values="week1_bg")
    for arm, d in runs.groupby("arm", sort=False):
        same = ARMS[arm][1:] == (1.0, 1.0)
        L.append(f"| {arm} | {len(d)} | {(w1[arm] - w1[REFERENCE]).abs().max():.3f}{'' if same else ' (differs by design)'} | "
                 f"{d.sens_ratio.median():.3f} | {100 * d.sens_moved.median():.0f}% |")
    L += ["", "## The arms", "",
          "Medians over the people who completed the arm.", "",
          "| Arm | Completed | Ended early, hypoglycaemic | Ended early, hyperglycaemic | Daily dose, U | Glucose 5th / 50th / 95th percentile, mg/dL | Mean glucose, mg/dL | 70 to 180 | Below 70 | Below 54 | Above 180 | Above 250 |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for arm, d in runs.groupby("arm", sort=False):
        c = d[d.completed]
        L.append(f"| {arm} | {len(c)} | {int(d.termination.str.startswith('hypo').sum())} | "
                 f"{int(d.termination.str.startswith('hyper').sum())} | {c.dose.median():.1f} | {c.p5.median():.0f} / {c.p50.median():.0f} / {c.p95.median():.0f} | {c['mean'].median():.0f} | "
                 + " | ".join(f"{c[m].median():.2f}%" for m in ("tir", "tbr70", "tbr54", "tar180", "tar250")) + " |")
    L += ["", "## Each arm against its reference", "",
          "| Arm | Reference | People | " + " | ".join(METRICS) + " |", "|---|---|---:|" + "---:|" * len(METRICS)]
    for arm in runs.arm.unique():
        ref = REFERENCES.get(arm, REFERENCE)
        if arm == ref:
            continue
        ds = [paired(runs, arm, ref, m) for m in METRICS]
        L.append(f"| {arm} | {ref} | {len(ds[0])} | "
                 + " | ".join(_f(mean_ci(d)) if len(d) >= MIN_PEOPLE else "-" for d in ds) + " |")
    L += ["", "## The harm rule", "",
          "The first two rows are the null check. A side is decidable only where the placebo reads cleared.", "",
          "| Arm | Against | Side | People | Difference, points of time | Margin | Excess early ends | Reading |",
          "|---|---|---|---:|---:|---:|---:|---|"]
    pairs = [(NULL, REFERENCE)] + [(a, c) for a in ("candidate",) for c in COMPARATORS] \
        + [("candidate", "logarithmic, AndroidAPS")]
    for arm, ref in pairs:
        if {arm, ref} <= set(runs.arm):
            L += harm_rows(runs, arm, ref)
    worst = [f"{paired(runs, 'candidate', c, 'tbr54').max():+.2f} points against {c}"
             for c in COMPARATORS if len(paired(runs, "candidate", c, "tbr54"))]
    L += ["", "The last pair is outside the rule and is shown for comparison."]
    if worst:
        L += ["", "Largest rise in one person's time below 54 mg/dL, candidate less comparator: " + ", ".join(worst) + "."]
    return "\n".join(L) + "\n"


def panel_note(timsim: str, names: list[str]) -> str:
    """The panel's settings in aggregate, from the titrated profiles and no simulation."""
    _, _, titrate, _ = _sim(timsim)
    p = [titrate.profile_for(n) for n in names]
    isf, tdd = np.array([float(x.isf.mean()) for x in p]), np.array([float(x.total_daily_insulin) for x in p])
    q = lambda v, d=1: "%.{0}f (%.{0}f to %.{0}f)".format(d) % tuple(np.percentile(v, [50, 10, 90]))
    ratio = lambda g: curve_profile(g, isf) / stock_v1(g, tdd, DIVISOR)
    return "\n".join([
        "# INV-019: the simulated panel's settings", "",
        f"Generated by `inv019/outcome_test.py --stage panel`. {len(names)} adults, settings titrated by TimSim. "
        "Median with the 10th to 90th percentile over people.", "",
        "| Quantity | Value |", "|---|---:|",
        f"| Profile ISF, mg/dL per U | {q(isf)} |",
        f"| People whose ISF schedule is flat across the day | {sum(len(set(np.round(x.isf, 6))) == 1 for x in p)} |",
        f"| Titrated daily dose, U | {q(tdd)} |",
        f"| Profile ISF x daily dose | {q(isf * tdd, 0)} |",
        f"| Factor at which the logarithmic model returns the profile ISF at {ANCHOR_GLUCOSE:.0f} mg/dL | {q(matched_factor(isf, tdd), 2)} |",
        *[f"| Candidate sensitivity over logarithmic at 100%, at {g} mg/dL | {q(ratio(float(g)), 2)} |" for g in (100, 140, 180, 250, 300)],
    ]) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--timsim", required=True, help="path to the TimSim checkout")
    ap.add_argument("--stage", choices=("panel", "arms", "report", "all"), default="all")
    ap.add_argument("--people", type=int, default=PANEL)
    ap.add_argument("--days", type=int, default=DAYS)
    ap.add_argument("--workers", type=int, default=WORKERS)
    ap.add_argument("--arms", nargs="*", default=list(ARMS), help="a subset, for a pilot")
    ap.add_argument("--tag", default="", help="suffix for a pilot, so it cannot overwrite the run")
    a = ap.parse_args()
    timsim = str(Path(a.timsim).expanduser().resolve())
    _sim(timsim)
    from examples.panel import dosable
    from timsim import provenance, seeds
    names, seed = dosable(a.people), seeds.BENCH
    raw = RESULTS / f"inv019_raw{a.tag}"

    if a.stage == "panel":
        (RESULTS / "inv019_panel.md").write_text(panel_note(timsim, names))
        print((RESULTS / "inv019_panel.md").read_text())
        return
    if a.stage in ("arms", "all"):
        raw.mkdir(exist_ok=True)
        rows = _pool(arm_job, [(timsim, n, arm, a.days, seed, str(raw)) for arm in a.arms for n in names], a.workers)
        pd.DataFrame(rows).to_csv(raw / "status.csv", index=False)
    if a.stage in ("report", "all"):
        status = pd.read_csv(raw / "status.csv", keep_default_na=False)
        runs = pd.DataFrame([reduce(raw / f"{r.arm}__{r.subject}.npz") for r in status[status.error == ""].itertuples()])
        runs["arm"] = pd.Categorical(runs.arm, [x for x in ARMS if x in set(runs.arm)])
        runs = runs.sort_values(["arm", "subject"]).astype({"arm": str})
        s = provenance.stamp()
        dirty = subprocess.run(["git", "-C", timsim, "status", "--porcelain"], capture_output=True, text=True).stdout
        stamp = (f"TimSim commit {str(s['commit'])[:7]}{' (uncommitted changes)' if dirty.strip() else ''}, "
                 f"population {s['population']}, settings {s['settings']}, simglucose {s['simglucose']}, "
                 f"seed {seed}, {a.people} adults, {a.days} days.")
        text = report(runs, int((status.error != "").sum()), stamp)
        (RESULTS / f"inv019_outcome_test{a.tag}.md").write_text(text)
        print(text)


if __name__ == "__main__":
    main()
