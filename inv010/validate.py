#!/usr/bin/env python3
"""Phase one: is the INV-009 autosens reconstruction faithful to a real device?

INV-009 compared sensitivity detectors using a port of the AndroidAPS plugins.
That port needed correcting three times, and the last correction inverted its
conclusion. It has never been checked against a device.

The Boost cohort records what a device actually computed. Two fields matter and
they are not the same thing. `sensitivityRatio` is what Boost applied, and
`boostAutosens_mode` says 'tdd' for every user, so that field is Boost's own
dose-derived ratio rather than autosens. The oref autosens value sits separately
in `boostAutosens_orefRatio`, and that is what this port should reproduce.

Four of the ten users with that field return exactly 1.000 for every record, so
autosens is not running for them. Six have a ratio that moves.

The comparison is run at two levels, to place any error rather than merely find
one. Level A feeds the device's own insulin activity and sensitivity into my
deviation and ratio arithmetic, which tests that arithmetic alone. Level B
reconstructs activity from the pump record as INV-009 does, which tests the whole
pipeline. Agreement at A and not at B puts the fault in the insulin
reconstruction.

The maximum daily basal that the ratio divides by is not recorded, so the level
of the ratio cannot be reproduced directly. Correlation and the sign of the
departure from neutral can, and those carry the claim.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd
import psycopg2
from scipy import stats as sps

DSN = "dbname=oref"
BIN_MIN = 5.0
AUTOSENS_MIN, AUTOSENS_MAX = 0.7, 1.2
LOOKBACK_H = 24.0
MIN_VALID = 12
OUT = "results/inv010_phase1.json"


def q(sql, params=()):
    with psycopg2.connect(DSN) as c, c.cursor() as cur:
        cur.execute(sql, params)
        return pd.DataFrame(cur.fetchall(), columns=[d[0] for d in cur.description])


def load(user):
    d = q("""
        SELECT ts_utc,
               (openaps->'suggested'->>'bg')::float                         AS bg,
               (openaps->'iob'->>'activity')::float                         AS activity,
               NULLIF(openaps->'suggested'->>'variable_sens','')::float     AS vsens,
               NULLIF(openaps->'suggested'->>'boostAutosens_orefRatio','')::float AS oref_ratio,
               (openaps->'suggested'->>'sensitivityRatio')::float           AS applied,
               NULLIF(openaps->'suggested'->>'COB','')::float               AS cob
        FROM boost_devicestatus_raw
        WHERE user_id = %s AND openaps->'suggested' ? 'boostAutosens_orefRatio'
        ORDER BY ts_utc""", (user,))
    if d.empty:
        return None
    d["ts_utc"] = pd.to_datetime(d.ts_utc, utc=True)
    for c in ("bg", "activity", "vsens", "oref_ratio", "applied", "cob"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    # The device posts more often than every five minutes, so collapse to the grid.
    d = d.set_index("ts_utc").resample("5min").median().dropna(subset=["bg", "activity"])
    return d.reset_index()


def max_daily_basal(user):
    """Highest programmed basal rate, approximated from what the pump actually ran."""
    d = q("""SELECT rate FROM boost_treatments
             WHERE user_id=%s AND rate IS NOT NULL AND rate > 0""", (user,))
    if d.empty:
        return np.nan
    return float(np.nanpercentile(pd.to_numeric(d.rate, errors="coerce").dropna(), 99))


def ratio_from_deviation(dev, sens, mdb):
    """The step every AndroidAPS sensitivity plugin shares."""
    n_back = int(LOOKBACK_H * 60 / BIN_MIN)
    med = pd.Series(dev).rolling(n_back, min_periods=MIN_VALID).median().to_numpy()
    basal_off = med * (60.0 / BIN_MIN) / sens
    return np.clip(1.0 + basal_off / mdb, AUTOSENS_MIN, AUTOSENS_MAX), basal_off


def main() -> int:
    users = q("""SELECT user_id, count(*) n FROM boost_devicestatus_raw
                 WHERE openaps->'suggested' ? 'boostAutosens_orefRatio'
                 GROUP BY 1 ORDER BY 2 DESC""").user_id.tolist()
    rows = []
    print(f"{'user':>6s} {'n':>7s} {'moves':>6s} {'corr':>7s} {'p':>9s} "
          f"{'sign agree':>11s} {'implied mdb':>12s} {'observed':>9s}")
    for u in users:
        d = load(u)
        if d is None or len(d) < 500:
            continue
        moves = float(d.oref_ratio.std()) > 0.005
        sens = d.vsens.to_numpy(float)
        act = d.activity.to_numpy(float)
        bg = d.bg.to_numpy(float)

        # oref: bgi = -activity * sens * 5, deviation = delta - bgi
        dev = np.full(len(d), np.nan)
        dev[1:] = (bg[1:] - bg[:-1]) + act[:-1] * sens[:-1] * BIN_MIN
        # AAPS excludes deviations while carbohydrate is absorbing, which the
        # device reports directly here rather than needing a proxy.
        if d.cob.notna().any():
            dev = np.where(np.nan_to_num(d.cob.to_numpy(float)) > 0, np.nan, dev)
        dev = np.where((bg < 80) & (dev > 0), 0.0, dev)

        mdb = max_daily_basal(u)
        mine, basal_off = ratio_from_deviation(dev, np.nanmedian(sens), mdb)
        rec = d.oref_ratio.to_numpy(float)
        ok = np.isfinite(mine) & np.isfinite(rec)
        if ok.sum() < 200 or not moves:
            print(f"{u:>6s} {len(d):7d} {'no':>6s} {'':>7s} {'':>9s} {'':>11s} "
                  f"{'':>12s} {mdb:9.2f}   autosens flat at 1.000")
            rows.append(dict(user=u, n=int(len(d)), moves=False, mdb_observed=mdb))
            continue
        r = sps.spearmanr(mine[ok], rec[ok])
        # Sign agreement: do we and the device put the ratio the same side of 1?
        both_off = ok & ((np.abs(rec - 1) > 0.005) | (np.abs(mine - 1) > 0.005))
        agree = float(np.mean(np.sign(mine[both_off] - 1) == np.sign(rec[both_off] - 1)))
        # What maximum daily basal would the device's own ratio imply, given my
        # deviations? If the arithmetic is right this should land near the observed.
        m2 = ok & np.isfinite(basal_off) & (np.abs(rec - 1) > 0.02)
        implied = float(np.nanmedian(basal_off[m2] / (rec[m2] - 1.0))) if m2.sum() > 50 else np.nan
        print(f"{u:>6s} {len(d):7d} {'yes':>6s} {r.statistic:+7.3f} {r.pvalue:9.1e} "
              f"{100*agree:10.0f}% {implied:12.2f} {mdb:9.2f}")
        rows.append(dict(user=u, n=int(len(d)), moves=True, spearman=float(r.statistic),
                         p=float(r.pvalue), sign_agreement=agree,
                         mdb_implied=implied, mdb_observed=mdb,
                         recorded_sd=float(np.nanstd(rec)), mine_sd=float(np.nanstd(mine))))
    R = pd.DataFrame(rows)
    R.to_parquet("results/inv010_phase1.parquet", index=False)
    mv = R[R.moves.astype(bool)]
    res = dict(n_users=int(len(R)), n_moving=int(len(mv)),
               median_spearman=float(mv.spearman.median()) if len(mv) else None,
               median_sign_agreement=float(mv.sign_agreement.median()) if len(mv) else None,
               per_user=rows)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=2, default=float)
    if len(mv):
        print(f"\n  {len(mv)} users with a moving ratio: median correlation "
              f"{mv.spearman.median():+.3f}, median sign agreement "
              f"{100*mv.sign_agreement.median():.0f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
