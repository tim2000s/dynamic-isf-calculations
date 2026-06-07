"""Stage 2 — per-user ISF replay: every CGM tick → variableSens under V1 and V2.

v5:    TDD = device-logged sug_TDD (forward-filled ≤ 60 min).
v6/v7: TDD = stage-1 blended windows, joined via the recovered absolute-time anchor
       (validated against the DB `hour` column; corrected mod-24h if needed).

Each row also carries a flat-TDD arm (treatments total / span) so the report can show
how much the W8H blend itself matters.

Output: inv008_cache/replay/<user>.parquet + <user>.meta.json
"""
from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message="pandas only supports SQLAlchemy")

from inv008 import config
from inv008.dynisf import (isf_v1, isf_v2, sens_normal_target_v1,
                              sens_normal_target_v2)

FLAT_TDD = json.loads((config.ROOT / "user_treatments_tdd.json").read_text()) \
    if (config.ROOT / "user_treatments_tdd.json").exists() else {}


def _db():
    import psycopg2
    return psycopg2.connect(**config.DB)


def _load_ticks(user_id: str, table: str, with_tdd: bool) -> pd.DataFrame:
    cols = 'ts_relative_sec, cgm_mgdl, hour' + (', "sug_TDD" AS sug_tdd' if with_tdd else '')
    if table in ("oref_v6", "oref_v7"):
        cols = "ts_relative_sec, cgm_mgdl, hour"
    with _db() as conn:
        df = pd.read_sql(
            f"SELECT {cols} FROM {table} WHERE user_id = %s "
            f"AND cgm_mgdl IS NOT NULL ORDER BY ts_relative_sec",
            conn, params=(user_id,))
    return df


def _hour_mismatch(sample: pd.DataFrame, anchor: float) -> float:
    implied = ((anchor + sample["ts_relative_sec"].to_numpy()) % 86400) // 3600
    return float((implied != sample["hour"].to_numpy()).mean())


def refine_anchor(df: pd.DataFrame, anchor: float) -> tuple[float, float]:
    """Refine an anchor against the DB hour column: first the mod-24h shift
    (circular mode of the hour offset), then a ±30 min sub-hour sweep."""
    sample = df.iloc[:: max(1, len(df) // 2000)]
    implied = ((anchor + sample["ts_relative_sec"].to_numpy()) % 86400) // 3600
    diff = (sample["hour"].to_numpy() - implied) % 24
    vals, counts = np.unique(diff, return_counts=True)
    anchor += int(vals[np.argmax(counts)]) * 3600
    best = (anchor, _hour_mismatch(sample, anchor))
    for delta in range(-1800, 1801, 60):
        m = _hour_mismatch(sample, anchor + delta)
        if m < best[1]:
            best = (anchor + delta, m)
    return best


def _join_tdd(df: pd.DataFrame, win: pd.DataFrame, anchor: float) -> pd.DataFrame:
    out = df.copy()
    out["abs_ts"] = (anchor + out["ts_relative_sec"]).astype("int64")
    return pd.merge_asof(out.sort_values("abs_ts"),
                         win.rename(columns={"ts": "abs_ts"}),
                         on="abs_ts", direction="backward",
                         tolerance=2 * config.GRID_SEC)


def choose_anchor(df: pd.DataFrame, win: pd.DataFrame,
                  anchor_start: float | None) -> tuple[pd.DataFrame, dict]:
    """Try candidate anchors (archive start; TDD-span end aligned to DB end),
    refine each against the hour column, keep the one with the best TDD-join
    coverage (ties → lower hour mismatch)."""
    rel_max = float(df["ts_relative_sec"].max())
    candidates = []
    if anchor_start is not None:
        candidates.append(("archive_start", float(anchor_start)))
    candidates.append(("span_end", float(win["ts"].max()) - rel_max))

    best = None
    for name, raw in candidates:
        a, mism = refine_anchor(df, raw)
        joined = _join_tdd(df, win, a)
        cov = float(joined["tdd_blend"].notna().mean())
        cand = (cov, -mism, name, a, joined)
        if best is None or cand[:2] > best[:2]:
            best = cand
    cov, neg_mism, name, a, joined = best
    return joined, {"anchor_epoch_sec": a, "anchor_candidate": name,
                    "anchor_hour_mismatch": round(-neg_mism, 4),
                    "anchor_uncertain": -neg_mism > config.ANCHOR_HOUR_MISMATCH_MAX}


def run_user(args: tuple[str, str]) -> dict:
    user_id, platform = args
    out_pq = config.REPLAY_DIR / f"{user_id}.parquet"
    out_meta = config.REPLAY_DIR / f"{user_id}.meta.json"
    table = {"v5": "oref_v5", "v6": "oref_v6", "v7": "oref_v7"}[platform]

    df = _load_ticks(user_id, table, with_tdd=(platform == "v5"))
    if len(df) < 1000:
        return {"user": user_id, "status": "skip", "reason": f"only {len(df)} ticks"}

    meta: dict = {"user": user_id, "platform": platform, "n_ticks": int(len(df)),
                  "source_commit": config.SOURCE_COMMIT}

    if platform == "v5":
        tdd = df["sug_tdd"].ffill(limit=12)  # ≤ 60 min carry-forward
        meta["tdd_source"] = "device sug_TDD"
        meta["tdd_coverage"] = round(float(tdd.notna().mean()), 4)
    else:
        meta_path = config.TDD_DIR / f"{user_id}.meta.json"
        pq_path = config.TDD_DIR / f"{user_id}.parquet"
        if not (meta_path.exists() and pq_path.exists()):
            return {"user": user_id, "status": "skip", "reason": "no stage-1 output"}
        s1 = json.loads(meta_path.read_text())
        win = pd.read_parquet(pq_path)
        df, anchor_meta = choose_anchor(df, win, s1.get("anchor_epoch_sec"))
        meta.update(anchor_meta)
        meta["basal_source"] = s1["basal_source"]
        tdd = df["tdd_blend"]
        meta["tdd_source"] = f"reconstructed ({s1['basal_source']})"
        meta["tdd_coverage"] = round(float(tdd.notna().mean()), 4)

    bg = df["cgm_mgdl"].to_numpy()
    tdd_arr = tdd.to_numpy(dtype=float)

    out = pd.DataFrame({
        "ts_relative_sec": df["ts_relative_sec"].to_numpy(),
        "bg": bg,
        "tdd": tdd_arr,
        "sens_nt_v1": sens_normal_target_v1(tdd_arr),
        "sens_nt_v2": sens_normal_target_v2(tdd_arr),
        "isf_v1": isf_v1(bg, tdd_arr),
        "isf_v2": isf_v2(bg, tdd_arr),
    })

    flat = FLAT_TDD.get(user_id)
    if flat and flat.get("span_days", 0) > 0:
        flat_tdd = flat["total_insulin_units"] / flat["span_days"]
        # v6 flat TDD from boluses-only sources still includes basal via profile in
        # stage 1; the flat arm here is treatments-total/span as recorded.
        out["isf_v1_flat"] = isf_v1(bg, flat_tdd)
        out["isf_v2_flat"] = isf_v2(bg, flat_tdd)
        meta["flat_tdd"] = round(flat_tdd, 2)

    tmp = out_pq.with_suffix(".parquet.tmp")
    out.to_parquet(tmp, index=False)
    tmp.rename(out_pq)

    valid = np.isfinite(out["isf_v1"]) & np.isfinite(out["isf_v2"])
    meta["valid_frac"] = round(float(valid.mean()), 4)
    if valid.any():
        meta["median_isf_v1"] = round(float(np.nanmedian(out.loc[valid, "isf_v1"])), 1)
        meta["median_isf_v2"] = round(float(np.nanmedian(out.loc[valid, "isf_v2"])), 1)
        meta["median_tdd"] = round(float(np.nanmedian(out.loc[valid, "tdd"])), 1)
    out_meta.write_text(json.dumps(meta, indent=1))
    return {"user": user_id, "status": "ok",
            **{k: meta.get(k) for k in ("tdd_source", "tdd_coverage", "valid_frac",
                                        "median_tdd", "median_isf_v1", "median_isf_v2",
                                        "anchor_hour_mismatch")}}


def user_list(platforms: tuple[str, ...]) -> list[tuple[str, str]]:
    out = []
    if "v5" in platforms:
        import psycopg2
        with psycopg2.connect(**config.DB) as conn, conn.cursor() as cur:
            cur.execute('SELECT DISTINCT user_id FROM oref_v5 WHERE "sug_TDD" IS NOT NULL')
            out += [(r[0], "v5") for r in cur.fetchall()]
    from inv008 import sources
    mapping = sources.load_mappings()
    out += [(uid, plat) for uid, (plat, _) in sorted(mapping.items()) if plat in platforms]
    return out
