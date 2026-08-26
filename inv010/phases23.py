#!/usr/bin/env python3
"""Phases two and three, on recorded device fields only.

Phase one showed the INV-009 reconstruction does not track this device, so
anything leaning on it is set aside. What follows uses only values the device
wrote down, so it stands whatever the port does.

Phase two asks how often oref autosens is actually doing anything.
Phase three asks what Boost's dose-derived ratio and oref's autosens do to each
other, since both are adjusting the same quantity from different evidence.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import psycopg2
from scipy import stats as sps

DSN = "dbname=oref"
LO, HI = 0.7, 1.2


def q(sql, p=()):
    with psycopg2.connect(DSN) as c, c.cursor() as cur:
        cur.execute(sql, p)
        return pd.DataFrame(cur.fetchall(), columns=[d[0] for d in cur.description])


def main() -> int:
    d = q("""
        SELECT user_id, ts_utc,
               NULLIF(openaps->'suggested'->>'boostAutosens_orefRatio','')::float AS oref,
               (openaps->'suggested'->>'sensitivityRatio')::float                 AS applied,
               NULLIF(openaps->'suggested'->>'variable_sens','')::float           AS vsens,
               NULLIF(openaps->'suggested'->>'TDD','')::float                     AS tdd,
               (openaps->'suggested'->>'bg')::float                               AS bg
        FROM boost_devicestatus_raw
        WHERE openaps->'suggested' ? 'sensitivityRatio'""")
    for c in ("oref", "applied", "vsens", "tdd", "bg"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["ts_utc"] = pd.to_datetime(d.ts_utc, utc=True)

    print("PHASE TWO  how often is oref autosens doing anything?\n")
    print(f"  {'user':>6s} {'n':>7s} {'exactly 1':>10s} {'at ceiling':>11s} "
          f"{'at floor':>9s} {'active':>7s} {'median bg':>10s} {'% below 70':>11s}")
    rows = []
    for u, g in d.dropna(subset=["oref"]).groupby("user_id"):
        v = g.oref.to_numpy(float)
        neutral = float(np.mean(np.isclose(v, 1.0, atol=1e-6)))
        ceil = float(np.mean(v >= HI - 1e-6))
        floor = float(np.mean(v <= LO + 1e-6))
        bg = g.bg.dropna()
        rows.append(dict(user=u, n=int(len(g)), neutral=neutral, ceiling=ceil,
                         floor=floor, active=1 - neutral,
                         median_bg=float(bg.median()) if len(bg) else np.nan,
                         pct_below_70=float((bg < 70).mean()) if len(bg) else np.nan))
        print(f"  {u:>6s} {len(g):7d} {100*neutral:9.1f}% {100*ceil:10.1f}% "
              f"{100*floor:8.1f}% {100*(1-neutral):6.1f}% {rows[-1]['median_bg']:10.0f} "
              f"{100*rows[-1]['pct_below_70']:10.1f}%")
    P2 = pd.DataFrame(rows)
    print(f"\n  across {len(P2)} users: autosens sits exactly at 1.000 for a median of "
          f"{100*P2.neutral.median():.0f}% of decisions")
    print(f"  it never moves at all for {int((P2.active < 0.01).sum())} of them")

    print("\nPHASE THREE  what the two mechanisms do to each other\n")
    print("  Boost applies a dose-derived ratio; oref autosens computes its own from")
    print("  recent deviations. Both scale the same sensitivity.\n")
    print(f"  {'user':>6s} {'applied':>8s} {'oref':>7s} {'corr':>8s} {'p':>9s} "
          f"{'both move':>10s} {'opposed':>9s}")
    rows3 = []
    for u, g in d.dropna(subset=["oref", "applied"]).groupby("user_id"):
        a, o = g.applied.to_numpy(float), g.oref.to_numpy(float)
        if np.std(a) < 1e-6 or np.std(o) < 1e-6:
            print(f"  {u:>6s} {np.mean(a):8.3f} {np.mean(o):7.3f} "
                  f"{'':>8s} {'':>9s} {'one is flat':>10s}")
            rows3.append(dict(user=u, applied=float(np.mean(a)), oref=float(np.mean(o)),
                              flat=True))
            continue
        r = sps.spearmanr(a, o)
        both = (np.abs(a - 1) > 0.02) & (np.abs(o - 1) > 0.02)
        opposed = float(np.mean(np.sign(a[both] - 1) != np.sign(o[both] - 1))) if both.sum() > 50 else np.nan
        print(f"  {u:>6s} {np.mean(a):8.3f} {np.mean(o):7.3f} {r.statistic:+8.3f} "
              f"{r.pvalue:9.1e} {100*both.mean():9.1f}% {100*opposed:8.0f}%")
        rows3.append(dict(user=u, applied=float(np.mean(a)), oref=float(np.mean(o)),
                          spearman=float(r.statistic), p=float(r.pvalue),
                          both_move=float(both.mean()), opposed=opposed, flat=False))
    P3 = pd.DataFrame(rows3)
    mv = P3[~P3.flat.astype(bool)]
    if len(mv):
        print(f"\n  of {len(mv)} users where both move: median correlation "
              f"{mv.spearman.median():+.3f}")
        opp = mv.opposed.dropna()
        if len(opp):
            print(f"  when both are away from neutral they point opposite ways "
                  f"{100*opp.median():.0f}% of the time")

    json.dump(dict(phase2=rows, phase3=rows3), open("results/inv010_phases23.json", "w"),
              indent=2, default=float)
    return 0


if __name__ == "__main__":
    sys.exit(main())
