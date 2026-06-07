#!/usr/bin/env python3
"""Extract per-user treatments-derived TDD from raw NS / AAPS dumps.

Replaces the hybrid TDD `max(2 × basal, basal + SMB/day, sug_tdd)` with
ground-truth `basal_profile + sum(treatments insulin) / span_days`.

Sources:
  • v7 (oref0):        /240 NS samples/<dir>/direct-sharing-31/<id>_treatments_*.json
  • v6 (AAPS classic): /240 NS samples/<dir>/direct-sharing-396/upload-*/Treatments.json
  • v5 (Trio):         /240 NS samples/<dir>/direct-sharing-31/...  (same NS format as v7
                       where present)

Output: user_treatments_tdd.json keyed by user_id with fields:
  total_insulin_units, span_days, n_events, source
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
import os
from pathlib import Path

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
NS = ROOT / "240 NS samples"
OUT = ROOT / "user_treatments_tdd.json"

MV6 = json.loads((ROOT / "user_mapping_v6.json").read_text())
MV7 = json.loads((ROOT / "user_mapping_v7.json").read_text())


def parse_ts(s):
    """Lenient ISO-8601 parser. Returns datetime in UTC or None."""
    if s is None:
        return None
    if isinstance(s, (int, float)):
        # Some treatments use epoch ms
        try:
            ts = float(s)
            if ts > 1e12:  # ms
                ts /= 1000.0
            return datetime.fromtimestamp(ts, tz=timezone.utc)
        except Exception:
            return None
    s = str(s)
    try:
        # Strip trailing 'Z', accept optional fractional / TZ
        if s.endswith("Z"):
            return datetime.fromisoformat(s[:-1] + "+00:00").astimezone(timezone.utc)
        d = datetime.fromisoformat(s[:25])
        if d.tzinfo is None:
            d = d.replace(tzinfo=timezone.utc)
        return d.astimezone(timezone.utc)
    except Exception:
        return None


def sum_treatments(records: list) -> tuple[float, int, datetime, datetime]:
    """Sum insulin>0 events; return (total_units, n_events, t_min, t_max)."""
    total = 0.0
    n = 0
    t_min = None
    t_max = None
    for d in records:
        if not isinstance(d, dict):
            continue
        ins = d.get("insulin")
        if ins is None:
            continue
        try:
            ins = float(ins)
        except (TypeError, ValueError):
            continue
        if ins <= 0:
            continue
        total += ins
        n += 1
        ts = parse_ts(d.get("timestamp") or d.get("created_at") or d.get("date"))
        if ts is None:
            continue
        if t_min is None or ts < t_min:
            t_min = ts
        if t_max is None or ts > t_max:
            t_max = ts
    return total, n, t_min, t_max


def extract_v7(user_dir: Path) -> dict | None:
    """oref0 / v7: <id>_treatments_*.json under direct-sharing-31."""
    ds31 = user_dir / "direct-sharing-31"
    if not ds31.is_dir():
        return None
    files = sorted(ds31.glob("*treatments*.json"))
    if not files:
        return None
    total = 0.0
    n = 0
    t_min = t_max = None
    for f in files:
        try:
            data = json.loads(f.read_text())
        except Exception:
            continue
        if not isinstance(data, list):
            data = [data]
        s, k, a, b = sum_treatments(data)
        total += s
        n += k
        if a and (t_min is None or a < t_min):
            t_min = a
        if b and (t_max is None or b > t_max):
            t_max = b
    if n == 0 or t_min is None or t_max is None:
        return None
    span = (t_max - t_min).total_seconds() / 86400.0
    if span <= 0:
        return None
    return {
        "total_insulin_units": round(total, 2),
        "span_days": round(span, 2),
        "n_events": n,
        "source": "v7/" + user_dir.name,
    }


def extract_v6(user_dir: Path) -> dict | None:
    """AAPS / v6: Treatments.json under direct-sharing-396/upload-*/."""
    ds396 = user_dir / "direct-sharing-396"
    if not ds396.is_dir():
        return None
    candidates = list(ds396.glob("upload-*/Treatments.json"))
    if not candidates:
        return None
    total = 0.0
    n = 0
    t_min = t_max = None
    for f in candidates:
        try:
            data = json.loads(f.read_text())
        except Exception:
            continue
        if not isinstance(data, list):
            data = [data]
        s, k, a, b = sum_treatments(data)
        total += s
        n += k
        if a and (t_min is None or a < t_min):
            t_min = a
        if b and (t_max is None or b > t_max):
            t_max = b
    if n == 0 or t_min is None or t_max is None:
        return None
    span = (t_max - t_min).total_seconds() / 86400.0
    if span <= 0:
        return None
    return {
        "total_insulin_units": round(total, 2),
        "span_days": round(span, 2),
        "n_events": n,
        "source": "v6/" + user_dir.name,
    }


def main():
    out = {}
    for uid, raw in MV7.items():
        ud = NS / raw
        if not ud.is_dir():
            continue
        rec = extract_v7(ud)
        if rec:
            out[uid] = rec
    for uid, raw in MV6.items():
        ud = NS / raw
        if not ud.is_dir():
            continue
        rec = extract_v6(ud)
        if rec:
            out[uid] = rec
    OUT.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT} with {len(out)} users")
    # Distribution sanity
    spans = sorted(r["span_days"] for r in out.values())
    totals = sorted(r["total_insulin_units"] / r["span_days"] for r in out.values())
    if totals:
        print(f"  per-day insulin (treatments / span): min={totals[0]:.1f}, "
              f"med={totals[len(totals)//2]:.1f}, max={totals[-1]:.1f} U/day")
        print(f"  span_days: min={spans[0]:.0f}, med={spans[len(spans)//2]:.0f}, max={spans[-1]:.0f}")


if __name__ == "__main__":
    main()
