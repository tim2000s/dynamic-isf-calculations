"""Sensitivity from every unit the system delivered, not only the ones a person pressed.

An algorithm that raises temporary basal above the scheduled profile has given a
correction, and one that suspends has given a negative one. Treating those nights
as untreated controls, which the bolus-matched estimator did, discards most of the
dosing variation in every closed-loop study and leaves the contrast to a handful
of manual boluses. Here the exposure is a_net, the action within the window from
insulin delivered above or below the person's own scheduled profile, whatever
route it arrived by.

The identification problem is that the algorithm decides from the CGM trace, so
its dose is a response to glucose as well as a cause of it. The design answers it
by stratifying finely on the trace the algorithm itself saw: the same person, the
same hour, the same starting glucose, the same slope into the window and the same
level an hour earlier. Conditional on all of that the algorithm's intent is close
to fixed, and what varies within a stratum is insulin on board, the shape of the
person's own basal profile at that hour and the granularity of temporary basal.
The estimate is a within-stratum regression of the fall on that residual dose,
pooled across strata, with intervals bootstrapped over people rather than windows
because windows from one person are not independent.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import config, data

STUDIES = ["ReplaceBG", "Loop", "DCLP3", "DCLP5", "PEDAP", "IOBP2"]
LABEL = {"ReplaceBG": "REPLACE-BG (open loop)", "Loop": "Loop (DIY closed loop)",
         "DCLP3": "DCLP3 (Control-IQ, adult)", "DCLP5": "DCLP5 (Control-IQ, 6-13y)",
         "PEDAP": "PEDAP (Control-IQ, 2-5y)", "IOBP2": "IOBP2 (bionic pancreas)"}

BG_BIN, SLOPE_BIN, M60_BIN = 15.0, 1.0, 25.0
MIN_PER_STRATUM = 4
MIN_SPREAD_U = 0.4


def strata(w: pd.DataFrame, use_m60: bool = True) -> pd.Series:
    """Label each window by the trace the algorithm saw before it dosed."""
    k = [w.subject_id.astype(str), w.hour.astype(int).astype(str),
         np.round(w.bg0 / BG_BIN).astype(int).astype(str),
         np.round(w.pre_slope / SLOPE_BIN).astype(int).astype(str)]
    if use_m60:
        k.append(np.round(w.bg_m60 / M60_BIN).astype(int).astype(str))
    return pd.Series(["|".join(t) for t in zip(*k)], index=w.index)


def within(w: pd.DataFrame, expo: str, use_m60: bool = True) -> pd.DataFrame:
    """Demean the fall and the dose inside each stratum, keeping usable strata."""
    d = w[["subject_id", "hour", "bg0", "pre_slope", "bg_m60", "drop", expo]].dropna()
    if d.empty:
        return d
    d = d.assign(k=strata(d, use_m60))
    g = d.groupby("k")[expo]
    keep = (g.transform("size") >= MIN_PER_STRATUM) & \
           (g.transform(lambda x: x.max() - x.min()) >= MIN_SPREAD_U)
    d = d[keep]
    if d.empty:
        return d
    return d.assign(dy=d["drop"] - d.groupby("k")["drop"].transform("mean"),
                    dx=d[expo] - d.groupby("k")[expo].transform("mean"))


def fit(d: pd.DataFrame, n_boot: int = 400, seed: int = 20260826) -> dict:
    """Pooled within-stratum slope, with people resampled rather than windows."""
    if len(d) < 200 or (d.dx ** 2).sum() <= 0:
        return dict(isf=np.nan, lo=np.nan, hi=np.nan, n=int(len(d)),
                    n_subj=int(d.subject_id.nunique()) if len(d) else 0, n_strata=0)
    slope = float((d.dx * d.dy).sum() / (d.dx ** 2).sum())
    subs = d.subject_id.unique()
    idx = {s: g for s, g in d.groupby("subject_id", sort=False)}
    rng = np.random.default_rng(seed)
    b = []
    for _ in range(n_boot):
        pick = pd.concat([idx[s] for s in rng.choice(subs, len(subs), replace=True)])
        den = (pick.dx ** 2).sum()
        if den > 0:
            b.append((pick.dx * pick.dy).sum() / den)
    return dict(isf=slope, lo=float(np.percentile(b, 2.5)),
                hi=float(np.percentile(b, 97.5)), n=int(len(d)),
                n_subj=int(len(subs)), n_strata=int(d.k.nunique()))


def endogeneity(w: pd.DataFrame) -> dict:
    """How tightly the system's later dosing tracks what it already committed.

    Within a stratum the algorithm has seen the same trace, so any remaining
    association between insulin committed before the window and insulin given
    during it is the algorithm continuing to respond to something the trace does
    not show. A human bolusing once and going to sleep produces almost none of
    this; a controller sampling every five minutes produces a great deal, and
    that is what makes its dose endogenous to the outcome being measured.
    """
    d = w[["subject_id", "hour", "bg0", "pre_slope", "bg_m60",
           "a_net_pre", "a_net_in"]].dropna()
    if len(d) < 500:
        return dict(n=int(len(d)), corr=np.nan)
    d = d.assign(k=strata(d, True))
    d = d[d.groupby("k").a_net_pre.transform("size") >= MIN_PER_STRATUM]
    if len(d) < 500:
        return dict(n=int(len(d)), corr=np.nan)
    x = d.a_net_pre - d.groupby("k").a_net_pre.transform("mean")
    y = d.a_net_in - d.groupby("k").a_net_in.transform("mean")
    if x.std() <= 0 or y.std() <= 0:
        return dict(n=int(len(d)), corr=np.nan)
    return dict(n=int(len(d)), corr=float(np.corrcoef(x, y)[0, 1]))


def by_dose_size(ev: pd.DataFrame) -> list[dict]:
    """Matched treated-minus-control fall, split by how much was given.

    Reported because the pooled ratio hides two things: very small doses where
    noise dominates, and very large ones which are given by people whose daily
    dose is high and whose sensitivity is therefore genuinely lower.
    """
    rows = []
    d = ev[(ev.bg0 >= 150) & (ev.bg0 < 260)]
    for lo, hi in [(0.3, 1.0), (1.0, 2.0), (2.0, 3.5), (3.5, 20.0)]:
        b = d[(d.bolus_u >= lo) & (d.bolus_u < hi)]
        if len(b) < 25:
            continue
        dt, dc, ex = b.drop_t.mean(), b.drop_c.mean(), b.extra_action_u.mean()
        rows.append(dict(lo=lo, hi=hi, n=int(len(b)), treated=float(dt),
                         control=float(dc), acting_u=float(ex),
                         isf=float((dt - dc) / ex) if ex > 0 else np.nan))
    return rows


def main() -> int:
    config.ensure_dirs()
    ui = pd.read_parquet(config.RESULTS / "inv009_user_isf.parquet")
    ent = pd.read_parquet(config.RESULTS / "inv009_entered_isf.parquet")
    ent = ent[ent.isf > 0].set_index("subject_id").isf
    ev_all = pd.read_parquet(config.RESULTS / "inv009_events_full.parquet")

    rows = []
    for s in STUDIES:
        w = data.load(s)
        if w.empty:
            continue
        u = ui[ui.study == s]
        tdd = float(u.tdd_u.median()) if len(u) else np.nan
        e = float(ent.reindex(u.subject_id).dropna().median()) if len(u) else np.nan
        r = dict(study=s, label=LABEL[s], tdd_median=tdd, isf_entered=e,
                 closed_loop=bool(u.closed_loop.iloc[0]) if len(u) else None)
        for tag, expo, m60 in [("net", "a_net", True), ("netpre", "a_net_pre", True),
                               ("net_nom60", "a_net", False)]:
            f = fit(within(w, expo, m60))
            r.update({f"{tag}_{k}": v for k, v in f.items()})
        r["k_net"] = r["net_isf"] * tdd
        r["k_netpre"] = r["netpre_isf"] * tdd
        r["k_entered"] = e * tdd if np.isfinite(e) else np.nan
        r["net_over_entered"] = r["net_isf"] / e if e else np.nan
        r["endog_corr"] = endogeneity(w)["corr"]
        r["by_dose"] = by_dose_size(ev_all[ev_all.study == s])
        rows.append(r)

    t = pd.DataFrame([{k: v for k, v in r.items() if k != "by_dose"} for r in rows])
    t.to_parquet(config.RESULTS / "inv009_net_dose_isf.parquet", index=False)
    (config.RESULTS / "inv009_net_dose_isf.json").write_text(
        json.dumps(rows, indent=1, default=float))

    print("Sensitivity from insulin above or below the scheduled profile, whatever route")
    print("it arrived by, within person / hour / glucose / slope / level-an-hour-back\n")
    print(f"{'study':<26s}{'windows':>9s}{'people':>7s}{'strata':>7s}"
          f"{'ISF':>7s}{'95% CI':>16s}{'entered':>8s}{'ratio':>7s}")
    for r in rows:
        n = lambda v, w=7, d=1: (f"%{w}.{d}f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        ci = (f"{r['net_lo']:6.1f} to{r['net_hi']:6.1f}"
              if np.isfinite(r.get('net_lo', np.nan)) else f"{'-':>16}")
        print(f"{r['label']:<26s}{r['net_n']:9d}{r['net_n_subj']:7d}{r['net_n_strata']:7d}"
              f"{n(r['net_isf'])}{ci}{n(r['isf_entered'], 8, 0)}"
              f"{n(r.get('net_over_entered', np.nan), 7, 2)}")
    print("\nUsing only insulin committed before the window opened (a_net_pre),")
    print("which the algorithm could not have chosen in response to the window\n")
    print(f"{'study':<26s}{'windows':>9s}{'ISF':>7s}{'95% CI':>16s}")
    for r in rows:
        n = lambda v: f"{v:7.1f}" if np.isfinite(v) else f"{'-':>7}"
        ci = (f"{r['netpre_lo']:6.1f} to{r['netpre_hi']:6.1f}"
              if np.isfinite(r.get('netpre_lo', np.nan)) else f"{'-':>16}")
        print(f"{r['label']:<26s}{r['netpre_n']:9d}{n(r['netpre_isf'])}{ci}")
    print("\nAs a Walsh-style constant (sensitivity x median daily dose)\n")
    print(f"{'study':<26s}{'net dose':>10s}{'committed':>11s}{'entered':>9s}")
    for r in rows:
        n = lambda v, w=10: (f"%{w}.0f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        print(f"{r['label']:<26s}{n(r['k_net'])}"
              f"{n(r['netpre_isf'] * r['tdd_median'], 11)}{n(r['k_entered'], 9)}")
    print("\nHow endogenous the dose is: correlation between insulin committed before")
    print("the window and insulin given during it, within an identical trace\n")
    print(f"{'study':<26s}{'corr':>8s}")
    for r in rows:
        c = r.get("endog_corr", np.nan)
        print(f"{r['label']:<26s}" + (f"{c:8.3f}" if np.isfinite(c) else f"{'-':>8}"))
    print("\nMatched fall by dose given, 150-260 mg/dL (mg/dL, and ISF per acting unit)\n")
    for r in rows:
        if not r["by_dose"]:
            continue
        print(f"  {r['label']}")
        print(f"    {'dose U':>10s}{'n':>6s}{'treated':>9s}{'control':>9s}"
              f"{'diff':>7s}{'acting':>8s}{'ISF':>7s}")
        for b in r["by_dose"]:
            print(f"    {b['lo']:5.1f}-{b['hi']:<4.1f}{b['n']:6d}{b['treated']:9.0f}"
                  f"{b['control']:9.0f}{b['treated'] - b['control']:7.0f}"
                  f"{b['acting_u']:8.2f}{b['isf']:7.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
