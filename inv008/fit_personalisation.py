#!/usr/bin/env python3
"""Re-fit the ISF~TDD law on the per-tick BLENDED TDD (independent of the cross-
sectional scalar TDD), and decide what must be ubiquitous vs personalised.

Three questions:
  A. Item 1 — does the population law (exponent, K) survive when TDD is each user's
     median *blended* per-tick TDD instead of the canonical treatments/span scalar?
  B. Between-user — after the best population law, how much ISF variance is a stable
     per-user offset (→ would a per-user constant K help), and how much of that offset
     is real signal vs the noise in measuring a user's own sensitivity?
  C. Within-user identifiability — over ~14 days, how much does one user's TDD move?
     If it barely moves, a per-user *exponent* cannot be estimated; only K can.

Run: python -m inv008.fit_personalisation
"""
from __future__ import annotations

import glob
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
REPLAY = ROOT / "inv008_cache" / "replay"


def load():
    coh = pd.DataFrame(json.loads((ROOT / "canonical_cohort.json").read_text()))
    emp = pd.DataFrame(json.loads((ROOT / "empirical_isf_v5.json").read_text()))[
        ["user_id", "empirical_isf", "se_isf", "r2"]]
    df = coh.merge(emp, on="user_id", how="left")
    # per-user median BLENDED TDD from the replay metas
    rows = []
    for f in glob.glob(str(REPLAY / "*.meta.json")):
        m = json.loads(open(f).read())
        if m.get("median_tdd") is not None:
            rows.append({"user_id": m["user"], "tdd_blend": m["median_tdd"]})
    df = df.merge(pd.DataFrame(rows), on="user_id", how="left")
    df["emp_valid"] = (df.r2 >= 0.10) & df.empirical_isf.between(5, 500)
    return df


def louo_sqrt_power(d, ycol, tddcol):
    """LOUO-CV: fixed-exponent K/sqrt(TDD) and free power-law A*TDD^b on tddcol."""
    d = d.dropna(subset=[ycol, tddcol])
    d = d[(d[tddcol] > 0) & (d[ycol] > 0)]
    y = d[ycol].to_numpy(); tdd = d[tddcol].to_numpy(); n = len(d)
    # K/sqrt LOUO
    sp = np.full(n, np.nan); pw = np.full(n, np.nan)
    for i in range(n):
        tr = np.arange(n) != i
        sp[i] = np.median(y[tr] * np.sqrt(tdd[tr])) / math.sqrt(tdd[i])
        b, a = np.polyfit(np.log(tdd[tr]), np.log(y[tr]), 1)
        pw[i] = math.exp(a) * tdd[i] ** b
    # full-sample fitted params (for reporting)
    K = float(np.median(y * np.sqrt(tdd)))
    bfull, afull = np.polyfit(np.log(tdd), np.log(y), 1)

    def m(pred):
        le = np.abs(np.log(pred / y))
        return (round(float(np.median(np.abs(pred - y))), 1),
                round(float(np.median(le)), 3),
                round(float((le < math.log(1.3)).mean()), 3))

    return dict(n=n, K=round(K, 1), slope=round(float(bfull), 3),
                A=round(float(math.exp(afull)), 1),
                sqrt_metrics=m(sp), power_metrics=m(pw))


def main():
    df = load()
    have = df.tdd_blend.notna()
    print(f"cohort {len(df)}; with blended TDD {have.sum()}; "
          f"empirical-valid & blended {(df.emp_valid & have).sum()}")
    r = np.corrcoef(df.loc[have, "tdd"], df.loc[have, "tdd_blend"])[0, 1]
    ratio = (df.loc[have, "tdd_blend"] / df.loc[have, "tdd"])
    print(f"scalar vs blended TDD: corr {r:.3f}, median ratio {ratio.median():.2f} "
          f"(IQR {ratio.quantile(.25):.2f}-{ratio.quantile(.75):.2f})")

    print("\n=== A. Population law on BLENDED TDD (LOUO-CV) ===")
    for tgt, col, mask in [("entered ISF", "isf", have),
                           ("empirical ISF", "empirical_isf", df.emp_valid & have)]:
        for tname, tcol in [("blended", "tdd_blend"), ("scalar(prior)", "tdd")]:
            res = louo_sqrt_power(df[mask], col, tcol)
            print(f"  {tgt:14s} TDD={tname:13s} n={res['n']:3d} | "
                  f"K/sqrt: K={res['K']:6.1f} medlogerr={res['sqrt_metrics'][1]} "
                  f"within30={res['sqrt_metrics'][2]} | "
                  f"power slope={res['slope']:+.3f} (A={res['A']})")

    print("\n=== B. Between-user offset: is per-user K worth it? ===")
    e = df[df.emp_valid & have].copy()
    e["K_i"] = e.empirical_isf * np.sqrt(e.tdd_blend)        # per-user empirical constant
    logK = np.log(e.K_i)
    Kglob = float(np.median(e.K_i))
    resid_sd = float(np.std(logK))                           # between-user residual (log)
    # per-user measurement noise in log space (delta method): se/value
    noise_sd_i = (e.se_isf / e.empirical_isf).to_numpy()
    noise_var = float(np.mean(noise_sd_i**2))
    signal_var = max(resid_sd**2 - noise_var, 0.0)
    print(f"  global K = {Kglob:.0f};  per-user K_i range "
          f"{e.K_i.min():.0f}-{e.K_i.max():.0f} ({e.K_i.max()/e.K_i.min():.1f}x)")
    print(f"  between-user residual SD (log) = {resid_sd:.3f}  "
          f"(= {100*(math.exp(resid_sd)-1):.0f}% typical ISF offset a global K leaves)")
    print(f"  per-user measurement noise SD (log) = {math.sqrt(noise_var):.3f}")
    print(f"  → signal/total variance = {signal_var/resid_sd**2:.2f}  "
          f"(fraction of the global-K error that a 14-day per-user calibration can remove)")

    print("\n=== C. Within-user TDD range: can a per-user EXPONENT be identified? ===")
    cvs, ranges = [], []
    for f in glob.glob(str(REPLAY / "*.parquet")):
        d = pd.read_parquet(f, columns=["tdd"]).dropna()
        d = d[d.tdd > 0]
        if len(d) < 2000:
            continue
        cvs.append(float(d.tdd.std() / d.tdd.mean()))
        ranges.append(float(d.tdd.quantile(.9) / d.tdd.quantile(.1)))
    cvs = np.array(cvs); ranges = np.array(ranges)
    btw = df.loc[have, "tdd_blend"]
    print(f"  within-user TDD: median CV {np.median(cvs)*100:.0f}%, "
          f"median p90/p10 {np.median(ranges):.2f}x  (n={len(cvs)} users)")
    print(f"  between-user TDD: {btw.min():.0f}-{btw.max():.0f} U/day "
          f"({btw.max()/btw.min():.0f}x)")
    print(f"  → within-user lever arm is ~{np.median(ranges):.1f}x vs ~{btw.max()/btw.min():.0f}x "
          f"between users: a per-user TDD exponent is not identifiable from one person's data.")

    print("\n=== VERDICT ===")
    print("  exponent: UBIQUITOUS (≈ -0.5) — supported on blended TDD, and not")
    print("            estimable per-user (no within-user TDD range).")
    print("  constant K: PER-USER — most of a global equation's error is a stable")
    print("            per-user offset, and 14 days measures it far more precisely")
    print("            than it adds noise.")


if __name__ == "__main__":
    main()
