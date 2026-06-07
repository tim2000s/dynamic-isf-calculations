"""Stage 1 — per-user windowed-TDD reconstruction (v6/v7 only; v5 uses device sug_TDD).

Worker function `run_user` is process-safe: it does its own file I/O and returns a small
status dict. Output: inv008_cache/tdd/<user>.parquet + <user>.meta.json
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from inv008 import config, sources
from inv008.tdd_windows import build_delivery_grid, windowed_tdd
from inv008.dynisf import blend_tdd


def run_user(args: tuple[str, str, str]) -> dict:
    user_id, platform, raw_dir = args
    out_pq = config.TDD_DIR / f"{user_id}.parquet"
    out_meta = config.TDD_DIR / f"{user_id}.meta.json"

    hourly = sources.load_hourly_basal(user_id)
    if hourly is None:
        return {"user": user_id, "status": "skip", "reason": "no basal profile"}

    if platform == "v7":
        bolus_ts, bolus_u, temps = sources.load_v7_events(raw_dir)
    elif platform == "v6":
        bolus_ts, bolus_u, temps = sources.load_v6_events(raw_dir)
    else:
        return {"user": user_id, "status": "skip", "reason": f"stage1 n/a for {platform}"}

    if len(bolus_ts) < 50:
        return {"user": user_id, "status": "skip", "reason": f"only {len(bolus_ts)} boluses"}

    anchor = sources.recover_anchor(platform, raw_dir)

    t0 = float(bolus_ts.min())
    t1 = float(bolus_ts.max())
    if temps:
        t0 = min(t0, temps[0].ts)
        t1 = max(t1, temps[-1].ts)
    if anchor:
        t0 = min(t0, anchor)
    span_days = (t1 - t0) / 86400.0
    if span_days < 9:  # need ≥7 prior days for the blend to ever be valid
        return {"user": user_id, "status": "skip", "reason": f"span {span_days:.1f}d too short"}

    grid = build_delivery_grid(int(t0), int(t1), bolus_ts, bolus_u, temps, hourly)
    win = windowed_tdd(grid)
    win["tdd_blend"] = blend_tdd(win["tdd_4h"], win["tdd_8to4h"], win["tdd_1d"], win["tdd_7d"])

    tmp = out_pq.with_suffix(".parquet.tmp")
    win.to_parquet(tmp, index=False)
    tmp.rename(out_pq)

    valid_frac = float(np.isfinite(win["tdd_blend"]).mean())
    meta = {
        "user": user_id,
        "platform": platform,
        "raw_dir": raw_dir,
        "anchor_epoch_sec": anchor,
        "basal_source": "profile+temp" if temps else "profile_only",
        "n_boluses": int(len(bolus_ts)),
        "n_temp_events": int(len(temps)),
        "span_days": round(span_days, 1),
        "blend_valid_frac": round(valid_frac, 4),
        "tdd_blend_median": (round(float(np.nanmedian(win["tdd_blend"])), 2)
                             if valid_frac > 0 else None),
        "source_commit": config.SOURCE_COMMIT,
    }
    out_meta.write_text(json.dumps(meta, indent=1))
    return {"user": user_id, "status": "ok", **{k: meta[k] for k in
            ("basal_source", "span_days", "blend_valid_frac", "tdd_blend_median")}}


def user_list(platforms: tuple[str, ...]) -> list[tuple[str, str, str]]:
    mapping = sources.load_mappings()
    return [(uid, plat, raw) for uid, (plat, raw) in sorted(mapping.items())
            if plat in platforms]
