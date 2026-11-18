#!/usr/bin/env python3
"""Sensitivity variants of the four-hour overnight outcome ratio of section 4.5.

`inv008/overnight_sensitivity.py` divides the glucose fall over four hours by the
modelled IOB at a start between 23:00 and 02:55, keeps starts at or above target
with at least 0.30 U of IOB, and rejects a window in which any forward 15-minute
rise exceeds 2 mg/dL per five minutes. Its pooled median by starting-glucose band,
normalised to the lowest band, rose from 1.00 near 106 mg/dL to 2.26 near
198 mg/dL. Review 1, point M2, asked whether that rise survives other rise
thresholds and other choices of the denominator. This script reruns the same
estimator, from the same database rows and the same user list, under:

  the rise screen at 1, 2 (the committed value), 3 and 5 mg/dL per five minutes,
  and with no screen at all, each at the committed IOB minimum of 0.30 U;
  the IOB minimum at 0.10, 0.30, 0.50 and 1.00 U, each with the committed screen.

For each variant it reports the pooled median ratio by band, the curve normalised
to the lowest band, and a 95% percentile interval for every normalised value from
2,000 bootstrap draws that resample people. Every other rule of the committed
script is kept: the plausibility range of -50 to 400 mg/dL per U for the pooled
view, at least 30 windows per person, at least 200 windows per band.

Repeated records. The committed script read the oref tables without collapsing
repeated (user, timestamp) rows (REVIEW.md item 8). The baseline is computed both
ways: as committed, to show that this script reproduces 81,814 windows and the
published curve, and with rows collapsed to one per timestamp, which is the set
every variant is reported on.

Expected sign, written before the run. The committed curve rises with starting
glucose. Loosening the screen admits windows in which glucose rose later, which
lowers the ratio in every band; such windows are expected to be commoner where the
start is near target, because a fall from near target is limited by the loop
suspending, so loosening is expected to steepen the rise rather than remove it.
Tightening it removes windows in the other direction and is expected to flatten
the curve somewhat. Raising the IOB minimum removes windows with a small
denominator, which are commoner near target where the loop has withheld insulin.
The normalised value in the top band is expected to stay above 1 under every
variant; a variant in which its interval includes 1 would mean the rise does not
survive that change.

Added after the first run, on 18 November 2026. With no rise screen the pooled
median in the lowest band was -0.7 mg/dL per U, so the normalised curve divides
by a value near zero and its intervals are meaningless. The pooled median of every
band is therefore also given with its own interval, with the top band minus the
lowest band as the test of the rise, which needs no reference band away from zero.
The first run also showed that the committed estimator, run unchanged on the
database as it now stands, returns 81,769 windows where the committed output has
81,814: three Trio people return different counts, so the oref_v5 rows for them
have changed since that output was written. This script reproduces 81,769.

The ratio is an outcome ratio and is not read as a measured ISF.

    python3 -m inv008.overnight_variants --users ~/oref-investigations/empirical_isf_v5.json

Outputs results/overnight_variants.{json,md} (aggregates only).
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psycopg2

from inv008.overnight_sensitivity import (BG_CENTRES, BG_EDGES, HORIZON_S, SENS_RANGE,
                                          START_HOURS, TARGET, TOL_S)

DSN = "dbname=oref"
RISE_VARIANTS = (1.0, 2.0, 3.0, 5.0, None)      # None = no screen
IOB_VARIANTS = (0.10, 0.30, 0.50, 1.00)
BASE_RISE, BASE_IOB = 2.0, 0.30
MIN_PERSON, MIN_BAND = 30, 200
N_BOOT, SEED = 2000, 0
OUT = Path(__file__).resolve().parent.parent / "results"
PLATFORM = {"oref_v5": "Trio", "oref_v6": "AAPS", "oref_v7": "OpenAPS"}


def candidate_windows(df: pd.DataFrame) -> pd.DataFrame:
    """Every overnight start at or above target with a 4 h endpoint, before any screen.

    Carries the maximum forward 15-minute slope inside the window so that any rise
    threshold, and the IOB minimum, can be applied afterwards.
    """
    n = len(df)
    ts = df.ts_relative_sec.to_numpy(float)
    hour = df.hour.to_numpy(float)
    bg = df.cgm_mgdl.to_numpy(float)
    iob = df.iob_iob.to_numpy(float)
    end4 = np.searchsorted(ts, ts + HORIZON_S)
    end4 = np.where(end4 < n, end4, n - 1)
    ok_end = np.abs(ts[end4] - (ts + HORIZON_S)) <= TOL_S
    p15 = np.searchsorted(ts, ts + 900)
    p15 = np.where(p15 < n, p15, n - 1)
    ok15 = np.abs(ts[p15] - (ts + 900)) <= TOL_S
    slope15 = np.where(ok15, (bg[p15] - bg) / 3.0, 0.0)
    rows = []
    start = np.isin(hour, list(START_HOURS)) & ok_end & (iob >= min(IOB_VARIANTS)) & (bg >= TARGET)
    for i in np.flatnonzero(start):
        j = end4[i]
        seg = slope15[i:j + 1]
        if seg.size == 0:
            continue
        rows.append((bg[i], iob[i], bg[i] - bg[j], float(np.nanmax(seg))))
    return pd.DataFrame(rows, columns=["bg0", "iob0", "drop", "maxrise"])


def analyse(job):
    user_id, table = job
    conn = psycopg2.connect(DSN)
    try:
        df = pd.read_sql(f"SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob FROM {table} "
                         "WHERE user_id=%s AND cgm_mgdl IS NOT NULL AND iob_iob IS NOT NULL "
                         "ORDER BY ts_relative_sec", conn, params=(user_id,))
    finally:
        conn.close()
    if len(df) < 500:
        return None
    raw = candidate_windows(df)
    dd = df.drop_duplicates("ts_relative_sec").reset_index(drop=True)
    col = candidate_windows(dd) if len(dd) >= 500 else raw.iloc[:0]
    out = []
    for tag, w in (("raw", raw), ("collapsed", col)):
        if len(w):
            out.append(w.assign(user=user_id, platform=PLATFORM[table], rows=tag))
    return {"user": user_id, "n_rows": len(df), "n_distinct": len(dd),
            "windows": pd.concat(out, ignore_index=True) if out else None}


def select(W: pd.DataFrame, rise: float | None, iob_min: float) -> pd.DataFrame:
    """Apply one variant's screen, then the committed per-person minimum."""
    m = W.iob0 >= iob_min
    if rise is not None:
        m &= W.maxrise <= rise
    s = W[m]
    keep = s.groupby("user").bg0.transform("size") >= MIN_PERSON
    return s[keep]


def curve(S: pd.DataFrame, n_boot: int = N_BOOT, seed: int = SEED) -> dict:
    """Pooled median ratio by band, normalised to the lowest band, people resampled."""
    sens = (S["drop"] / S.iob0).to_numpy(float)
    ok = (sens > SENS_RANGE[0]) & (sens < SENS_RANGE[1])
    band = np.digitize(S.bg0.to_numpy(float), BG_EDGES) - 1
    users = S.user.to_numpy()
    sens, band, users = sens[ok], band[ok], users[ok]
    codes, uniq = pd.factorize(users)
    nb = len(BG_CENTRES)

    def bands(sv, bv):
        med = np.full(nb, np.nan)
        for k in range(nb):
            m = bv == k
            if m.sum() >= MIN_BAND:
                med[k] = np.median(sv[m])
        return med

    point = bands(sens, band)
    order = np.argsort(codes, kind="stable")
    groups = np.split(order, np.flatnonzero(np.diff(codes[order])) + 1)
    rng = np.random.default_rng(seed)
    absd = np.full((n_boot, nb), np.nan)
    for b in range(n_boot):
        idx = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
        absd[b] = bands(sens[idx], band[idx])
    top = nb - 1
    with np.errstate(all="ignore"):
        draws = absd / absd[:, :1]
        norm = point / point[0]
        lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0)
        alo, ahi = np.nanpercentile(absd, [2.5, 97.5], axis=0)
        diff = absd[:, top] - absd[:, 0]
        dlo, dhi = np.nanpercentile(diff, [2.5, 97.5])
    per_person = S.assign(s=S["drop"] / S.iob0).groupby("user").s.median()
    # A normalised curve means nothing when its reference band can be zero.
    stable = bool(alo[0] > 0)
    return {
        "n_people": int(S.user.nunique()), "n_windows": int(len(S)),
        "per_person_median_ratio": float(per_person.median()),
        "bands": {f"{c:.0f}": dict(pooled_median=None if np.isnan(point[k]) else float(point[k]),
                                   pooled_lo=float(alo[k]), pooled_hi=float(ahi[k]),
                                   normalised=None if np.isnan(norm[k]) else float(norm[k]),
                                   lo=float(lo[k]), hi=float(hi[k]),
                                   n_windows=int((band == k).sum()))
                  for k, c in enumerate(BG_CENTRES)},
        "normalisation_stable": stable,
        "top_minus_lowest": dict(diff=float(point[top] - point[0]), lo=float(dlo), hi=float(dhi),
                                 verdict=("interval above zero" if dlo > 0 else
                                          "interval below zero" if dhi < 0 else
                                          "interval includes zero")),
        "top_band_verdict": ("reference band interval includes zero" if not stable else
                             "interval above 1" if lo[top] > 1 else
                             "interval below 1" if hi[top] < 1 else "interval includes 1"),
    }


def self_check() -> None:
    """The screen and the IOB minimum must each change the selected windows."""
    W = pd.DataFrame({"bg0": [120.0] * 4 * MIN_PERSON, "iob0": [0.2, 0.4, 0.6, 1.2] * MIN_PERSON,
                      "drop": [10.0] * 4 * MIN_PERSON, "maxrise": [0.5, 1.5, 2.5, 6.0] * MIN_PERSON,
                      "user": ["a"] * 4 * MIN_PERSON})
    assert len(select(W, None, 0.1)) == 4 * MIN_PERSON
    assert len(select(W, 2.0, 0.1)) == 2 * MIN_PERSON
    assert len(select(W, None, 0.5)) == 2 * MIN_PERSON


def variants() -> list[tuple[str, float | None, float]]:
    v = [(f"rise screen {'none' if r is None else f'{r:g}'}, IOB >= {BASE_IOB:.2f} U", r, BASE_IOB)
         for r in RISE_VARIANTS]
    v += [(f"rise screen {BASE_RISE:g}, IOB >= {i:.2f} U", BASE_RISE, i)
          for i in IOB_VARIANTS if i != BASE_IOB]
    return v


def markdown(res: dict) -> str:
    cs = [f"{c:.0f}" for c in BG_CENTRES]
    L = ["# Four-hour overnight outcome ratio under other screens", "",
         "Generated by `inv008/overnight_variants.py`. Outcome ratio = glucose fall over four hours "
         "divided by modelled IOB at the start, mg/dL per U; it is not a measured ISF. Pooled median "
         "by starting-glucose band, normalised to the lowest band, with a 95% percentile interval "
         f"from {N_BOOT} bootstrap draws resampling people.", "",
         f"Rows read: {res['rows']['raw']:,}; distinct (user, timestamp): {res['rows']['distinct']:,}.", "",
         "## Reproduction of the committed run (repeated rows not collapsed)", ""]
    b = res["baseline_as_committed"]
    L.append(f"{b['n_people']} people, {b['n_windows']:,} windows. Pooled medians: "
             + ", ".join(f"{c} mg/dL {b['bands'][c]['pooled_median']:.1f}" for c in cs) + ".")
    L += ["", "## Variants, repeated rows collapsed", "",
          "| Variant | People | Windows | " + " | ".join(f"{c} mg/dL" for c in cs) + " | Top band verdict |",
          "|---|---:|---:|" + "---|" * len(cs) + "---|"]
    for name, r in res["variants"].items():
        cells = [f"{r['bands'][c]['normalised']:.2f} ({r['bands'][c]['lo']:.2f} to {r['bands'][c]['hi']:.2f})"
                 if r['bands'][c]['normalised'] is not None else "-" for c in cs]
        L.append(f"| {name} | {r['n_people']} | {r['n_windows']:,} | " + " | ".join(cells)
                 + f" | {r['top_band_verdict']} |")
    L += ["", "Pooled median ratio in the lowest band (the normalising value), mg/dL per U.", "",
          "| Variant | " + " | ".join(f"{c} mg/dL" for c in cs) + " |", "|---|" + "---:|" * len(cs)]
    for name, r in res["variants"].items():
        L.append(f"| {name} | " + " | ".join("-" if r['bands'][c]['pooled_median'] is None
                                              else f"{r['bands'][c]['pooled_median']:.1f}" for c in cs) + " |")
    L += ["", "Pooled median ratio in the top band minus that in the lowest band, mg/dL per U "
          "(95% interval). This does not depend on a reference band away from zero.", "",
          "| Variant | Lowest band (95% interval) | Top band (95% interval) | Top minus lowest | Verdict |",
          "|---|---|---|---|---|"]
    for name, r in list(res["variants"].items()) + list(res["by_platform"].items()):
        b0, bt, d = r["bands"][cs[0]], r["bands"][cs[-1]], r["top_minus_lowest"]
        L.append(f"| {name.replace(' | ', ', ')} | {b0['pooled_median']:.1f} ({b0['pooled_lo']:.1f} to {b0['pooled_hi']:.1f}) | "
                 f"{bt['pooled_median']:.1f} ({bt['pooled_lo']:.1f} to {bt['pooled_hi']:.1f}) | "
                 f"{d['diff']:+.1f} ({d['lo']:+.1f} to {d['hi']:+.1f}) | {d['verdict']} |")
    L += ["", "## By platform, committed screen and no screen, repeated rows collapsed", "",
          "| Platform | Variant | People | Windows | " + " | ".join(f"{c} mg/dL" for c in cs) + " | Top band verdict |",
          "|---|---|---:|---:|" + "---|" * len(cs) + "---|"]
    for key, r in res["by_platform"].items():
        p, name = key.split(" | ")
        cells = [f"{r['bands'][c]['normalised']:.2f} ({r['bands'][c]['lo']:.2f} to {r['bands'][c]['hi']:.2f})"
                 if r['bands'][c]['normalised'] is not None else "-" for c in cs]
        L.append(f"| {p} | {name} | {r['n_people']} | {r['n_windows']:,} | " + " | ".join(cells)
                 + f" | {r['top_band_verdict']} |")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--users", type=Path, required=True, help="empirical_isf_v5.json")
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    self_check()
    jobs = [(e["user_id"], e["table"]) for e in json.loads(a.users.read_text())]
    if a.limit:
        jobs = jobs[:a.limit]
    t0 = time.time()
    with mp.Pool(min(a.workers, 7), maxtasksperchild=8) as pool:
        res = [r for r in pool.imap_unordered(analyse, jobs) if r]
    print(f"{len(jobs)} users, {len(res)} read, {time.time() - t0:.0f} s")
    W = pd.concat([r["windows"] for r in res if r["windows"] is not None], ignore_index=True)
    raw, col = W[W.rows == "raw"], W[W.rows == "collapsed"]
    out = {"rows": {"raw": int(sum(r["n_rows"] for r in res)),
                    "distinct": int(sum(r["n_distinct"] for r in res))},
           "baseline_as_committed": curve(select(raw, BASE_RISE, BASE_IOB)),
           "variants": {}, "by_platform": {}}
    for name, rise, iob in variants():
        S = select(col, rise, iob)
        out["variants"][name] = curve(S)
        print(f"  {name:40s} people {S.user.nunique():4d} windows {len(S):7d}", flush=True)
    for p in ("Trio", "AAPS", "OpenAPS"):
        for name, rise in ((f"rise screen {BASE_RISE:g}", BASE_RISE), ("no rise screen", None)):
            S = select(col[col.platform == p], rise, BASE_IOB)
            if S.user.nunique() >= 5:
                out["by_platform"][f"{p} | {name}"] = curve(S)
    a.out.mkdir(exist_ok=True)
    (a.out / "overnight_variants.json").write_text(json.dumps(out, indent=1))
    md = markdown(out)
    (a.out / "overnight_variants.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
