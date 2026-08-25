"""Effective sensitivity: what a unit of insulin actually did, per person.

Two estimators, because one of them is honest about a problem the other cannot
escape.

The regression fits, per subject, the fall in glucose against the insulin that
acted during a window, holding starting glucose, its recent trend, the hour and
where glucose had been. Its slope is milligrams per decilitre per unit, which is
what a sensitivity factor is. Its weakness is that under a closed loop insulin is
close to a function of recent glucose, so conditioning on recent glucose leaves
very little independent variation in insulin to read. Where that bites, the
fitted slope goes negative, because the controller gives most insulin exactly
when glucose is refusing to come down. The fraction of subjects with a negative
slope is therefore reported for every cohort as a diagnostic rather than hidden.

The event estimator avoids that by comparing like with like: a night with an
isolated correction bolus against the same person's other nights that started at
the same glucose and the same hour with no bolus at all. The difference in fall,
divided by the extra insulin, is a sensitivity factor that needs no model of
basal need and no assumption about what the controller was thinking. Corrections
overnight are rare, so this is pooled across a cohort rather than fitted per
person.

    python3 -m inv009.effective_isf
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from . import config, data, stats

SPECS = {
    "pre":     "a_pre",       # insulin already committed when the window opened
    "total":   "a_tot",       # everything that acted, reactive insulin included
    "net_pre": "a_net_pre",   # committed insulin, as a deviation from routine
}
CONTROLS = ["bg0", "pre_slope", "bg_m60", "bg_m120"]


def _fit_subject(d: pd.DataFrame, acol: str) -> dict:
    """Ordinary least squares with a robust standard error, per subject."""
    d = d.dropna(subset=[acol, "drop"] + CONTROLS)
    if len(d) < config.MIN_WINDOWS:
        return {}
    a = d[acol].to_numpy(dtype=float)
    if np.std(a) < config.MIN_A_SD:
        return {}
    cols = [np.ones(len(d)), a, d.bg0.to_numpy() - 100.0, d.pre_slope.to_numpy(),
            d.bg_m60.to_numpy() - 100.0, d.bg_m120.to_numpy() - 100.0]
    hours = sorted(d.hour.unique())[1:]
    for hh in hours:
        cols.append((d.hour == hh).to_numpy(dtype=float))
    X = np.column_stack(cols)
    y = d["drop"].to_numpy(dtype=float)
    xtx = np.linalg.pinv(X.T @ X)
    beta = xtx @ (X.T @ y)
    resid = y - X @ beta
    # HC1 sandwich: windows overlap and their residuals are not identically
    # scattered, so a textbook standard error would be too confident.
    n, k = X.shape
    meat = (X * (resid ** 2)[:, None]).T @ X
    cov = xtx @ meat @ xtx * (n / max(n - k, 1))
    return dict(s=float(beta[1]), se=float(np.sqrt(max(cov[1, 1], 0.0))),
                n=int(n), sd_a=float(np.std(a)),
                rmse=float(np.sqrt(np.mean(resid ** 2))))


def per_subject(study: str) -> pd.DataFrame:
    w = data.load(study)
    if w.empty:
        return pd.DataFrame()
    rows = []
    for sid, d in w.groupby("subject_id", sort=False):
        base = dict(subject_id=sid, study=study, tdd_u=d.tdd_u.iloc[0],
                    tdd_basal_u=d.tdd_basal_u.iloc[0], basal_frac=d.basal_frac.iloc[0],
                    age=d.age.iloc[0], n_windows=len(d),
                    closed_loop=bool(d.closed_loop.iloc[0]))
        for name, acol in SPECS.items():
            f = _fit_subject(d, acol)
            base[f"s_{name}"] = f.get("s", np.nan)
            base[f"se_{name}"] = f.get("se", np.nan)
            base[f"n_{name}"] = f.get("n", 0)
            base[f"sd_a_{name}"] = f.get("sd_a", np.nan)
        rows.append(base)
    return pd.DataFrame(rows)


def event_isf(study: str, bg_tol: float = None, min_u: float = None) -> pd.DataFrame:
    """One sensitivity estimate per isolated overnight correction bolus.

    Treated nights carry a single correction near the start of the window and
    nothing else. Control nights are the same person, the same hour, and a
    starting glucose within tolerance, with no bolus at all. Matching within
    subject means any residual carbohydrate tail, any error in the person's basal
    need and any quirk of their insulin model is present in both arms and cancels.
    """
    bg_tol = config.EVENT_BG_MATCH if bg_tol is None else bg_tol
    min_u = config.EVENT_MIN_U if min_u is None else min_u
    w = data.load(study)
    if w.empty:
        return pd.DataFrame()
    out = []
    for sid, d in w.groupby("subject_id", sort=False):
        quiet = d.quiet_before_h >= config.EVENT_QUIET_H
        treated = d[quiet & (d.bolus_first30_u >= min_u)
                    & np.isclose(d.bolus_in_u, d.bolus_first30_u)]
        control = d[quiet & (d.bolus_in_u <= 1e-9)]
        if treated.empty or len(control) < 3:
            continue
        for t in treated.itertuples():
            c = control[(control.hour == t.hour) & ((control.bg0 - t.bg0).abs() <= bg_tol)]
            if len(c) < 3:
                continue
            extra = t.a_tot - c.a_tot.mean()
            if extra <= 0.05:
                continue
            out.append(dict(subject_id=sid, study=study, t0=t.t0, hour=t.hour,
                            bg0=t.bg0, tdd_u=t.tdd_u, age=t.age,
                            bolus_u=t.bolus_first30_u, extra_action_u=extra,
                            n_control=len(c),
                            isf_event=float((t.drop - c["drop"].mean()) / extra)))
    return pd.DataFrame(out)


def main() -> int:
    config.ensure_dirs()
    subj, events = [], []
    for study in config.COHORTS:
        s = per_subject(study)
        e = event_isf(study)
        if not s.empty:
            subj.append(s)
        if not e.empty:
            events.append(e)
        print(f"  {study:10s} subjects fitted {len(s):4d}  events {len(e):5d}", flush=True)
    S = pd.concat(subj, ignore_index=True)
    E = pd.concat(events, ignore_index=True) if events else pd.DataFrame()
    S.to_parquet(config.RESULTS / "inv009_user_isf.parquet", index=False)
    if not E.empty:
        E.to_parquet(config.RESULTS / "inv009_events.parquet", index=False)

    res: dict = {"specs": list(SPECS), "by_study": []}
    for study, d in S.groupby("study"):
        row = dict(study=study, label=config.COHORTS[study]["label"],
                   closed_loop=config.COHORTS[study]["closed_loop"], n=len(d))
        for name in SPECS:
            v = d[f"s_{name}"].dropna()
            row[f"median_{name}"] = float(v.median()) if len(v) else np.nan
            row[f"frac_positive_{name}"] = float((v > 0).mean()) if len(v) else np.nan
            row[f"n_{name}"] = int(len(v))
        if not E.empty and (E.study == study).any():
            ev = E[E.study == study]
            per = ev.groupby("subject_id").isf_event.median()
            row["event_n"] = int(len(ev))
            row["event_subjects"] = int(len(per))
            row["event_isf_median"] = float(np.median(ev.isf_event))
            row["event_isf_ci"] = [float(x) for x in
                                   np.percentile(ev.isf_event, [25, 75])]
        res["by_study"].append(row)
    (config.RESULTS / "inv009_effective_isf.json").write_text(json.dumps(res, indent=2))

    print("\nPer-subject regression slopes (mg/dL per acting unit)")
    print(f"{'cohort':26s} {'loop':>6s} {'n':>5s}  " +
          "  ".join(f"{k:>16s}" for k in SPECS))
    for r in res["by_study"]:
        cells = "  ".join(f"{r['median_' + k]:7.1f} ({100 * r['frac_positive_' + k]:3.0f}%+)"
                          for k in SPECS)
        print(f"{r['label']:26s} {'closed' if r['closed_loop'] else 'OPEN':>6s} "
              f"{r['n']:5d}  {cells}")
    print("\n  the percentage is the share of subjects with a positive slope; under a")
    print("  reactive controller it falls, which is the confound made visible")
    print("\nMatched correction-bolus estimator")
    for r in res["by_study"]:
        if "event_isf_median" in r:
            print(f"  {r['label']:26s} {r['event_n']:5d} events / {r['event_subjects']:3d} "
                  f"people   ISF {r['event_isf_median']:6.1f} "
                  f"(IQR {r['event_isf_ci'][0]:.1f} to {r['event_isf_ci'][1]:.1f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
