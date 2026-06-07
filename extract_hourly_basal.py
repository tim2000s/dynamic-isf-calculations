#!/usr/bin/env python3
"""Extract a 24-element hourly basal-rate vector per user from raw NS / AAPS dumps.

Output: user_basal_profiles.json  ->  {user_id: {"hourly_rates": [b0..b23], "source": "...", "total": float, "n_segments": int}}

Sources (in priority order):
  1. v7 (oref0 / OpenAPS): /240 NS samples/<dir>/direct-sharing-31/profile_*.json
     Nightscout profile docs with `store/<defaultProfile>/basal` segments.
  2. v6 (AAPS classic):   /240 NS samples/<dir>/direct-sharing-396/upload-*/ProfileSwitches.json
     AAPS records with `profile.basal` segments.

v3 / v5 cohorts: if their NS source dirs aren't present, they're skipped here and
will be omitted from the hourly model (still contribute to scalar models in train_profile_advisor.py).
"""
from __future__ import annotations

import json
import re
from collections import Counter
import os
from pathlib import Path

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
NS_SAMPLES = ROOT / "240 NS samples"
OUT_PATH = ROOT / "user_basal_profiles.json"

MAPPING_V6 = json.loads((ROOT / "user_mapping_v6.json").read_text())
MAPPING_V7 = json.loads((ROOT / "user_mapping_v7.json").read_text())


def segments_to_hourly(segments: list[dict]) -> list[float] | None:
    """Convert NS basal segment list (sorted by timeAsSeconds) to a 24-hour vector.

    Each segment carries a `value` (U/h) effective from `timeAsSeconds` until the next
    segment's timeAsSeconds (wrapping at 86400). We sample the rate at the top of each hour.
    """
    if not segments:
        return None
    norm = []
    for seg in segments:
        try:
            t = int(float(seg.get("timeAsSeconds", 0)))
            v = float(seg.get("value"))
        except (TypeError, ValueError):
            continue
        norm.append((t, v))
    if not norm:
        return None
    norm.sort()
    # Ensure first segment starts at t=0 (NS convention guarantees this; defensively pad).
    if norm[0][0] != 0:
        norm.insert(0, (0, norm[0][1]))
    out = []
    for hour in range(24):
        sec = hour * 3600
        rate = norm[0][1]
        for t, v in norm:
            if t <= sec:
                rate = v
            else:
                break
        out.append(round(rate, 4))
    return out


def pick_v7_profile(profile_doc: dict) -> dict | None:
    """Pick the active profile from a Nightscout `profile` document."""
    store = profile_doc.get("store") or {}
    if not store:
        return None
    default_name = profile_doc.get("defaultProfile")
    if default_name and default_name in store:
        return store[default_name]
    # Fallback: first store entry.
    return next(iter(store.values()))


def extract_v7_user(user_dir: Path) -> tuple[list[float], str] | None:
    """Look in direct-sharing-31 for profile JSONs; return hourly basal + source filename."""
    ds31 = user_dir / "direct-sharing-31"
    if not ds31.is_dir():
        return None
    # Filenames vary: "profile_*.json", "<rawid>_profile_*.json", "<rawid>_profile.json".
    # Match any JSON whose name contains "profile" but excludes other doc types.
    profile_files = sorted(p for p in ds31.glob("*profile*.json")
                           if not any(t in p.name for t in
                                      ("treatment", "entries", "devicestatus")))
    if not profile_files:
        return None
    # Walk newest-first; pick the first profile that yields a parseable basal vector.
    for pf in reversed(profile_files):
        try:
            data = json.loads(pf.read_text())
        except Exception:
            continue
        records = data if isinstance(data, list) else [data]
        # NS returns newest-first within the array; iterate as-is.
        for record in records:
            ps = pick_v7_profile(record)
            if not ps:
                continue
            hourly = segments_to_hourly(ps.get("basal") or [])
            if hourly:
                return hourly, f"{pf.name}#{record.get('defaultProfile', '?')}"
    return None


def extract_v6_user(user_dir: Path) -> tuple[list[float], str] | None:
    """Look in direct-sharing-396/upload-*/ProfileSwitches.json; return hourly basal."""
    ds396 = user_dir / "direct-sharing-396"
    if not ds396.is_dir():
        return None
    upload_dirs = sorted(ds396.glob("upload-*"))
    if not upload_dirs:
        return None
    # Aggregate: most-frequently-used profile across all switches wins, since one-off
    # profile switches (eg activity) shouldn't define the user's baseline.
    rate_counter: Counter[tuple] = Counter()
    sample_source = None
    for ud in upload_dirs:
        ps_file = ud / "ProfileSwitches.json"
        if not ps_file.is_file():
            continue
        try:
            switches = json.loads(ps_file.read_text())
        except Exception:
            continue
        if not isinstance(switches, list):
            continue
        for sw in switches:
            if sw.get("isDeletion"):
                continue
            if sw.get("isValid") is False:
                continue
            prof = sw.get("profile") or {}
            hourly = segments_to_hourly(prof.get("basal") or [])
            if hourly:
                rate_counter[tuple(hourly)] += 1
                if sample_source is None:
                    sample_source = f"{ud.name}/ProfileSwitches.json"
    if not rate_counter:
        return None
    most_common, _ = rate_counter.most_common(1)[0]
    return list(most_common), sample_source or "ProfileSwitches.json"


def main() -> None:
    out: dict[str, dict] = {}
    skipped: dict[str, list[str]] = {"v6_no_dir": [], "v6_no_basal": [], "v7_no_dir": [], "v7_no_basal": []}

    for uid, raw_id in MAPPING_V7.items():
        user_dir = NS_SAMPLES / raw_id
        if not user_dir.is_dir():
            skipped["v7_no_dir"].append(uid)
            continue
        result = extract_v7_user(user_dir)
        if result is None:
            skipped["v7_no_basal"].append(uid)
            continue
        hourly, source = result
        out[uid] = {
            "hourly_rates": hourly,
            "source": f"v7/{raw_id}/{source}",
            "total_24h": round(sum(hourly), 3),
            "n_distinct_rates": len(set(hourly)),
        }

    for uid, raw_id in MAPPING_V6.items():
        user_dir = NS_SAMPLES / raw_id
        if not user_dir.is_dir():
            skipped["v6_no_dir"].append(uid)
            continue
        result = extract_v6_user(user_dir)
        if result is None:
            skipped["v6_no_basal"].append(uid)
            continue
        hourly, source = result
        out[uid] = {
            "hourly_rates": hourly,
            "source": f"v6/{raw_id}/{source}",
            "total_24h": round(sum(hourly), 3),
            "n_distinct_rates": len(set(hourly)),
        }

    OUT_PATH.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT_PATH} with {len(out)} users")
    for k, v in skipped.items():
        if v:
            print(f"  skipped {k}: {len(v)} -> {v[:5]}{'...' if len(v) > 5 else ''}")


if __name__ == "__main__":
    main()
