"""Raw-data adapters for OREF-INV-008: per-user treatments, temp basals, basal profiles
and absolute-time anchors.

Formats (verified against the archive 2026-06-07):
  v7 (direct-sharing-31): NS-style docs. Boluses = any record with insulin > 0.
      Temp basals = eventType "Temp Basal" with duration + absolute (U/h) or percent;
      records with neither are cancellations. isFakedTempBasal → extended-bolus
      emulation, treated like any other temp segment.
  v6 (direct-sharing-396/upload-*): AAPS export. Treatments.json = flat bolus records
      (date ms, insulin, isSMB, isValid, isDeletion). NO temp basal records exist in
      these uploads → basal contribution approximated by the profile schedule
      (flagged basal_source="profile_only").
  v5 (Trio): no raw files needed — oref_v5.sug_TDD is the device-logged TDD.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from inv008 import config
from inv008.tdd_windows import TempBasalEvent


def parse_ts(s) -> float | None:
    """Lenient timestamp → epoch seconds (same rules as extract_treatments_tdd.py)."""
    if s is None:
        return None
    if isinstance(s, (int, float)):
        ts = float(s)
        if ts > 1e12:
            ts /= 1000.0
        return ts if ts > 0 else None
    s = str(s)
    try:
        if s.endswith("Z"):
            return datetime.fromisoformat(s[:-1] + "+00:00").timestamp()
        d = datetime.fromisoformat(s[:25])
        if d.tzinfo is None:
            d = d.replace(tzinfo=timezone.utc)
        return d.timestamp()
    except Exception:
        return None


def load_mappings() -> dict[str, tuple[str, str]]:
    """user_id → (platform, raw directory name)."""
    out = {}
    for uid, raw in json.loads(config.USER_MAPPING_V7.read_text()).items():
        out[uid] = ("v7", raw)
    for uid, raw in json.loads(config.USER_MAPPING_V6.read_text()).items():
        out[uid] = ("v6", raw)
    return out


def load_hourly_basal(user_id: str) -> np.ndarray | None:
    profiles = json.loads(config.BASAL_PROFILES.read_text())
    rec = profiles.get(user_id)
    if rec:
        return np.asarray(rec["hourly_rates"], dtype=float)
    # v6 fallback: enriched profile scalars. The users missing from
    # user_basal_profiles.json all have profile_n_basal_segments == 1, so a flat
    # mean-basal vector reproduces their schedule exactly.
    enriched = config.ROOT / "user_profiles_v6_enriched.json"
    if enriched.exists():
        rec = json.loads(enriched.read_text()).get(user_id)
        if rec and rec.get("profile_mean_basal"):
            return np.full(24, float(rec["profile_mean_basal"]))
    return None


# ---------------------------------------------------------------- v7 (NS docs)

def load_v7_events(raw_dir: str):
    """→ (bolus_ts, bolus_units, [TempBasalEvent]) from all *treatments*.json files."""
    ds31 = config.NS_SAMPLES / raw_dir / "direct-sharing-31"
    bolus_ts, bolus_u, temps = [], [], []
    for f in sorted(ds31.glob("*treatments*.json")):
        try:
            recs = json.loads(f.read_text())
        except Exception:
            continue
        if not isinstance(recs, list):
            recs = [recs]
        for r in recs:
            if not isinstance(r, dict):
                continue
            ts = parse_ts(r.get("timestamp") or r.get("created_at") or r.get("date"))
            if ts is None:
                continue
            ins = r.get("insulin")
            if ins is not None:
                try:
                    ins = float(ins)
                except (TypeError, ValueError):
                    ins = None
                if ins and ins > 0:
                    bolus_ts.append(ts)
                    bolus_u.append(ins)
            if r.get("eventType") == "Temp Basal":
                dur = r.get("duration")
                try:
                    dur = float(dur) if dur is not None else 0.0
                except (TypeError, ValueError):
                    dur = 0.0
                rate = r.get("absolute", r.get("rate"))
                pct = r.get("percent")
                try:
                    rate = float(rate) if rate is not None else None
                except (TypeError, ValueError):
                    rate = None
                try:
                    pct = 100.0 + float(pct) if pct is not None else None
                    # NS `percent` is relative to profile: -100..+ (0 = profile rate)
                except (TypeError, ValueError):
                    pct = None
                temps.append(TempBasalEvent(ts=ts, duration_min=dur,
                                            rate_u_h=rate, percent=pct))
    order = np.argsort(bolus_ts) if bolus_ts else []
    return (np.asarray(bolus_ts, dtype=float)[order] if len(bolus_ts) else np.array([]),
            np.asarray(bolus_u, dtype=float)[order] if len(bolus_u) else np.array([]),
            sorted(temps, key=lambda e: e.ts))


# ------------------------------------------------------------ v6 (AAPS export)

def load_v6_events(raw_dir: str):
    """→ (bolus_ts, bolus_units, []) from upload-*/Treatments.json. No temp basals."""
    ds396 = config.NS_SAMPLES / raw_dir / "direct-sharing-396"
    bolus_ts, bolus_u = [], []
    for f in sorted(ds396.glob("upload-*/Treatments.json")):
        try:
            recs = json.loads(f.read_text())
        except Exception:
            continue
        if not isinstance(recs, list):
            continue
        for r in recs:
            if not isinstance(r, dict) or r.get("isDeletion") or r.get("isValid") is False:
                continue
            ins = r.get("insulin")
            try:
                ins = float(ins) if ins is not None else 0.0
            except (TypeError, ValueError):
                continue
            if ins <= 0:
                continue
            ts = parse_ts(r.get("date") or r.get("timestamp"))
            if ts is None:
                continue
            bolus_ts.append(ts)
            bolus_u.append(ins)
    order = np.argsort(bolus_ts) if bolus_ts else []
    return (np.asarray(bolus_ts, dtype=float)[order] if len(bolus_ts) else np.array([]),
            np.asarray(bolus_u, dtype=float)[order] if len(bolus_u) else np.array([]),
            [])


# --------------------------------------------------------------- time anchors

def recover_anchor(platform: str, raw_dir: str) -> float | None:
    """Earliest decision timestamp in the raw archive ≈ the extractors' min(ts),
    which defined ts_relative_sec = ts - min(ts). Validated downstream against the
    DB `hour` column (see stage2.validate_anchor)."""
    if platform == "v7":
        ds = config.NS_SAMPLES / raw_dir / "direct-sharing-31"
        best = None
        for f in sorted(ds.glob("*devicestatus*.json")):
            try:
                recs = json.loads(f.read_text())
            except Exception:
                continue
            if not isinstance(recs, list):
                recs = [recs]
            for r in recs:
                if not isinstance(r, dict):
                    continue
                oap = r.get("openaps")
                if not isinstance(oap, dict):
                    continue
                sug = oap.get("suggested")
                if not isinstance(sug, dict):
                    continue
                ts = parse_ts(sug.get("timestamp") or r.get("created_at"))
                if ts and (best is None or ts < best):
                    best = ts
            # files are date-ranged; the first file with any decision bounds the min
            if best is not None:
                break
        return best
    if platform == "v6":
        ds = config.NS_SAMPLES / raw_dir / "direct-sharing-396"
        best = None
        for f in sorted(ds.glob("upload-*/BgReadings.json")):
            try:
                recs = json.loads(f.read_text())
            except Exception:
                continue
            for r in recs if isinstance(recs, list) else []:
                ts = parse_ts(r.get("date"))
                if ts and (best is None or ts < best):
                    best = ts
            if best is not None:
                break
        return best
    return None
