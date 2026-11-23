#!/usr/bin/env python3
"""How far the UTC day boundary moves the blended TDD and the V1 and V2 ISF.

The reconstructed delivery of `inv008/stage1_tdd.py` forms calendar-day totals
on UTC midnight (`tdd_windows.windowed_tdd`), while AndroidAPS forms them on local
midnight. The day totals feed two terms of the blend, yesterday's total and the
seven-day average, so a person far from UTC has delivery moved between adjacent
days. Section 3.1 of the paper states that the size of this was not measured.

This script rebuilds each reconstructed person's delivery grid exactly as stage 1
does and recomputes the five TDD windows and the blend with the day boundary moved
to midnight at UTC + offset, for offsets of -12, -9, -6, -3, +3, +6, +9 and
+12 hours, and for the person's own recorded offset where the raw records carry
one. It compares each shifted blend with the UTC blend at every five-minute bin
where both are defined.

The recorded offset. OpenAPS records (oref_v7, Nightscout documents) carry
timestamps with an explicit offset such as -05:00, and a `utcOffset` field in
minutes; AndroidAPS exports (oref_v6) carry `utcOffset` in milliseconds on their
profile-switch and temporary-basal records. The person's offset is the commonest
value among explicit-offset timestamp strings and non-zero `utcOffset` fields
(a zero `utcOffset` is a default in some uploaders and is not counted unless no
other value exists). shortcut: one offset per person, so a daylight-saving change
moves the boundary by up to an hour for part of the year; a per-record offset would
remove that.

Trio (oref_v5) is not recomputed: its TDD is the device's own logged value, which
the replay does not form from calendar days.

What is reported, per cohort: the median and 95th percentile of the absolute
relative change in blended TDD, and in each equation's ISF. At the same glucose
the glucose term cancels, so V1's ISF changes by TDD_utc / TDD_shifted - 1 and
V2's by (TDD_utc / TDD_shifted)^2 - 1; both are computed through `inv008/dynisf.py`
at 120 mg/dL as a check. Each person contributes an equal number of bins (a random
2,000, or all if fewer), so a long record does not dominate, and the 95% interval
on each figure is from 2,000 bootstrap draws that resample people. The share of
bins whose blend becomes defined or undefined with the shift is given beside it.

Expected sign, written before the run. A shifted boundary moves delivery between
neighbouring days without changing the total, so the seven-day average changes
little and yesterday's total changes most where delivery is lumpy around midnight,
which overnight boluses and temporary basal are. The blend gives yesterday a weight
of 0.33, so a typical change of a few per cent in yesterday's total is expected to
move the blend by about 1% at the median, rising with the size of the shift, with a
95th percentile of a few per cent. V2's change is about twice V1's because it squares
TDD. Changes larger than 10% at the median would mean the UTC boundary matters for
the replayed ISF.

    python3 -m inv008.day_boundary --root ~/oref-investigations

Outputs results/day_boundary.{json,md} (aggregates only). Per-person offsets are
not written.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import re
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

OFFSETS_H = (-12, -9, -6, -3, 3, 6, 9, 12)
SAMPLE_BINS = 2000
N_BOOT, SEED = 2000, 0
OUT = Path(__file__).resolve().parent.parent / "results"
COHORT = {"v6": "AAPS", "v7": "OpenAPS"}
_TZ = re.compile(r"T\d\d:\d\d(?::\d\d(?:\.\d+)?)?([+-])(\d\d):?(\d\d)$")


def offset_from_strings(values) -> Counter:
    c = Counter()
    for v in values:
        if isinstance(v, str):
            m = _TZ.search(v.strip())
            if m:
                sign = 1 if m.group(1) == "+" else -1
                c[sign * (int(m.group(2)) * 60 + int(m.group(3)))] += 1
    return c


def recorded_offset_min(platform: str, raw_dir: Path) -> int | None:
    """Commonest recorded offset in minutes east of UTC, or None."""
    strings, fields = Counter(), Counter()
    if platform == "v7":
        for f in sorted((raw_dir / "direct-sharing-31").glob("*treatments*.json")):
            try:
                recs = json.loads(f.read_text())
            except Exception:
                continue
            recs = recs if isinstance(recs, list) else [recs]
            recs = [r for r in recs if isinstance(r, dict)]
            strings += offset_from_strings(r.get("timestamp") or r.get("created_at") for r in recs)
            fields.update(int(r["utcOffset"]) for r in recs
                          if isinstance(r.get("utcOffset"), (int, float)))
    else:
        base = raw_dir / "direct-sharing-396"
        for name in ("ProfileSwitches.json", "TemporaryBasals.json"):
            for f in sorted(base.glob(f"upload-*/{name}")):
                try:
                    recs = json.loads(f.read_text())
                except Exception:
                    continue
                if isinstance(recs, list):
                    fields.update(int(round(r["utcOffset"] / 60000.0)) for r in recs
                                  if isinstance(r, dict) and isinstance(r.get("utcOffset"), (int, float)))
            if fields:
                break
    nonzero = Counter({k: v for k, v in fields.items() if k != 0})
    pool = strings + nonzero
    if not pool and fields:
        pool = fields
    return pool.most_common(1)[0][0] if pool else None


def rel_changes(base: np.ndarray, shifted: np.ndarray) -> dict[str, np.ndarray]:
    ok = np.isfinite(base) & np.isfinite(shifted) & (base > 0) & (shifted > 0)
    r = base[ok] / shifted[ok]
    return {"tdd": np.abs(shifted[ok] / base[ok] - 1.0), "v1": np.abs(r - 1.0),
            "v2": np.abs(r ** 2 - 1.0)}


def analyse(job):
    """One person: grid as stage 1 builds it, blends at every offset, sampled changes."""
    user_id, platform, raw_dir, root, seed = job
    import os
    os.environ["DYNISF_ROOT"] = root          # inv008.config reads it at import, in the worker
    from inv008 import dynisf, sources
    from inv008.tdd_windows import build_delivery_grid, windowed_tdd
    hourly = sources.load_hourly_basal(user_id)
    if hourly is None:
        return None
    if platform == "v7":
        bolus_ts, bolus_u, temps = sources.load_v7_events(raw_dir)
    else:
        bolus_ts, bolus_u, temps = sources.load_v6_events(raw_dir)
    if len(bolus_ts) < 50:
        return None
    anchor = sources.recover_anchor(platform, raw_dir)
    t0, t1 = float(bolus_ts.min()), float(bolus_ts.max())
    if temps:
        t0, t1 = min(t0, temps[0].ts), max(t1, temps[-1].ts)
    if anchor:
        t0 = min(t0, anchor)
    if (t1 - t0) / 86400.0 < 9:
        return None
    grid = build_delivery_grid(int(t0), int(t1), bolus_ts, bolus_u, temps, hourly)

    def blend(off_s):
        w = windowed_tdd(grid, day_offset_sec=off_s)
        return dynisf.blend_tdd(w.tdd_4h, w.tdd_8to4h, w.tdd_1d, w.tdd_7d)

    base = blend(0)
    if not np.isfinite(base).any():
        return None
    rec = recorded_offset_min(platform, Path(root) / "240 NS samples" / raw_dir)
    arms = {f"{h:+d} h": h * 3600 for h in OFFSETS_H}
    if rec is not None and rec != 0:
        arms["recorded offset"] = rec * 60
    rng = np.random.default_rng(seed)
    out = {"user": user_id, "cohort": COHORT[platform], "recorded_offset_min": rec, "arms": {}}
    for name, off in arms.items():
        sh = blend(off)
        ch = rel_changes(base, sh)
        n = len(ch["tdd"])
        if n == 0:
            continue
        take = rng.choice(n, min(n, SAMPLE_BINS), replace=False)
        fb, fs = np.isfinite(base), np.isfinite(sh)
        out["arms"][name] = {"sample": {k: v[take] for k, v in ch.items()},
                             "n_bins": int(n),
                             "frac_validity_changed": float(np.mean(fb != fs)),
                             "offset_s": int(off)}
    # check the closed forms against the equations at one glucose
    b, s = base[np.isfinite(base)][:5], blend(3 * 3600)[np.isfinite(base)][:5]
    ok = np.isfinite(s) & (s > 0) & (b > 0)
    if ok.any():
        r = dynisf.isf_v2(120.0, s[ok]) / dynisf.isf_v2(120.0, b[ok])
        assert np.allclose(r, (b[ok] / s[ok]) ** 2)
    return out


def summarise(rows: list[dict], arm: str, n_boot: int = N_BOOT, seed: int = SEED) -> dict | None:
    per = [r["arms"][arm] for r in rows if arm in r["arms"]]
    if len(per) < 2:
        return None
    rng = np.random.default_rng(seed)
    out = {"n_people": len(per), "n_bins_compared": int(sum(p["n_bins"] for p in per)),
           "frac_validity_changed_median": float(np.median([p["frac_validity_changed"] for p in per]))}
    for k in ("tdd", "v1", "v2"):
        arrs = [p["sample"][k] for p in per]
        allv = np.concatenate(arrs)
        med, p95 = np.median(allv), np.percentile(allv, 95)
        bm, bp = np.empty(n_boot), np.empty(n_boot)
        for i in range(n_boot):
            v = np.concatenate([arrs[j] for j in rng.integers(0, len(arrs), len(arrs))])
            bm[i], bp[i] = np.median(v), np.percentile(v, 95)
        out[k] = {"median": float(med), "median_lo": float(np.percentile(bm, 2.5)),
                  "median_hi": float(np.percentile(bm, 97.5)), "p95": float(p95),
                  "p95_lo": float(np.percentile(bp, 2.5)), "p95_hi": float(np.percentile(bp, 97.5))}
    return out


def self_check() -> None:
    b = np.array([40.0, 50.0, np.nan])
    s = np.array([44.0, 50.0, 30.0])
    c = rel_changes(b, s)
    assert np.allclose(c["tdd"], [0.1, 0.0]) and np.allclose(c["v1"], [abs(40 / 44 - 1), 0.0])
    assert np.allclose(c["v2"], [abs((40 / 44) ** 2 - 1), 0.0])
    assert offset_from_strings(["2019-01-01T10:00:00-05:00", "2019-01-01T10:00:00Z",
                                "2019-01-01T10:00:00.000+0530"]) == Counter({-300: 1, 330: 1})


def markdown(res: dict) -> str:
    L = ["# The UTC day boundary and the blended TDD", "",
         "Generated by `inv008/day_boundary.py`. Calendar-day totals recomputed with the day boundary at "
         "midnight UTC + offset, compared with the committed UTC boundary at every five-minute bin where "
         "both blends are defined. Absolute relative change, as a percentage; each person contributes an "
         f"equal number of bins; 95% intervals from {N_BOOT} bootstrap draws resampling people.", ""]
    for c, r in res["cohorts"].items():
        L += [f"## {c}", "", f"{r['n_people']} people recomputed. Recorded offset found for "
              f"{r['n_recorded']} ({r['n_recorded_nonzero']} non-zero)." if r['n_recorded'] else
              f"{r['n_people']} people recomputed. No recorded offset found.", "",
              "| Shift | People | Blended TDD median | TDD 95th pct | V1 ISF median | V1 95th pct "
              "| V2 ISF median | V2 95th pct | Bins whose blend validity changes (median person) |",
              "|---|---:|---|---|---|---|---|---|---:|"]
        for arm, s in r["arms"].items():
            if s is None:
                continue
            f = lambda k, q: f"{100 * s[k][q]:.2f}% ({100 * s[k][q + '_lo']:.2f} to {100 * s[k][q + '_hi']:.2f})"
            L.append(f"| {arm} | {s['n_people']} | {f('tdd', 'median')} | {f('tdd', 'p95')} | "
                     f"{f('v1', 'median')} | {f('v1', 'p95')} | {f('v2', 'median')} | {f('v2', 'p95')} | "
                     f"{100 * s['frac_validity_changed_median']:.2f}% |")
        L.append("")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, required=True,
                    help="the INV-008 input root (holds '240 NS samples' and the user mappings)")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--raw", type=Path, default=None,
                    help="pickle of the per-person samples, written before any reduction; kept out of git")
    a = ap.parse_args()
    self_check()
    jobs = []
    for plat, fname in (("v6", "user_mapping_v6.json"), ("v7", "user_mapping_v7.json")):
        m = json.loads((a.root / fname).read_text())
        items = sorted(m.items())[: a.limit or None]
        jobs += [(u, plat, raw, str(a.root), i) for i, (u, raw) in enumerate(items)]
    t0 = time.time()
    rows = []
    with mp.Pool(min(a.workers, 7), maxtasksperchild=4) as pool:
        for i, r in enumerate(pool.imap_unordered(analyse, jobs), 1):
            if r is not None:
                rows.append(r)
            if i % 20 == 0:
                print(f"  {i}/{len(jobs)} {time.time() - t0:.0f} s", flush=True)
    rows.sort(key=lambda r: r["user"])  # the pool returns people in any order (REVIEW.md item 21)
    if a.raw:
        import pickle
        a.raw.write_bytes(pickle.dumps(rows))     # raw output first (REVIEW.md item 16)
    res = {"offsets_h": list(OFFSETS_H), "cohorts": {}}
    for c in ("AAPS", "OpenAPS"):
        rc = [r for r in rows if r["cohort"] == c]
        recs = [r["recorded_offset_min"] for r in rc if r["recorded_offset_min"] is not None]
        arms = [f"{h:+d} h" for h in OFFSETS_H] + ["recorded offset"]
        res["cohorts"][c] = {
            "n_people": len(rc), "n_recorded": len(recs),
            "n_recorded_nonzero": sum(1 for x in recs if x != 0),
            "arms": {arm: summarise(rc, arm) for arm in arms}}
        print(f"  {c}: {len(rc)} people", flush=True)
    a.out.mkdir(exist_ok=True)
    (a.out / "day_boundary.json").write_text(json.dumps(res, indent=1))
    md = markdown(res)
    (a.out / "day_boundary.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
