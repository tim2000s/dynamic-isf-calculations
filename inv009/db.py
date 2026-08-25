"""Reading one subject's record out of the studies schema.

Connections are opened inside the worker that uses them. A connection handed
across a fork is a connection two processes then write to.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import psycopg2

from . import config


def _q(sql: str, params: tuple) -> pd.DataFrame:
    with psycopg2.connect(config.DSN) as c, c.cursor() as cur:
        cur.execute(sql, params)
        cols = [d[0] for d in cur.description]
        return pd.DataFrame(cur.fetchall(), columns=cols)


def subjects(study: str | None = None) -> pd.DataFrame:
    """Subjects with both glucose and insulin, which is the minimum to be usable."""
    sql = """
        SELECT s.subject_id, s.study_name, s.age_years
        FROM studies.subject s
        WHERE EXISTS (SELECT 1 FROM studies.cgm c WHERE c.subject_id = s.subject_id)
          AND (EXISTS (SELECT 1 FROM studies.basal b WHERE b.subject_id = s.subject_id)
               OR EXISTS (SELECT 1 FROM studies.bolus o WHERE o.subject_id = s.subject_id))
          AND (%s IS NULL OR s.study_name = %s)
        ORDER BY s.subject_id
    """
    return _q(sql, (study, study))


def streams(subject_id: str) -> dict[str, pd.DataFrame]:
    """Every stream this subject has, each sorted by time."""
    out = {
        "cgm": _q("SELECT ts_local, cgm_mgdl FROM studies.cgm WHERE subject_id=%s "
                  "ORDER BY ts_local", (subject_id,)),
        "basal": _q("SELECT ts_local, rate_u_hr FROM studies.basal WHERE subject_id=%s "
                    "ORDER BY ts_local", (subject_id,)),
        "bolus": _q("SELECT ts_local, bolus_u, delivery_duration_s FROM studies.bolus "
                    "WHERE subject_id=%s ORDER BY ts_local", (subject_id,)),
        "carbs": _q("SELECT ts_local, carbs_g FROM studies.carbs WHERE subject_id=%s "
                    "ORDER BY ts_local", (subject_id,)),
        "sched": _q("SELECT ts_local, sched_rate_u_hr FROM studies.basal_sched "
                    "WHERE subject_id=%s ORDER BY ts_local", (subject_id,)),
        "wizard": _q("SELECT ts_local, iob_u, isf_mgdl_per_u, cr_g_per_u, bg_input_mgdl, "
                     "carb_input_g FROM studies.wizard WHERE subject_id=%s ORDER BY ts_local",
                     (subject_id,)),
    }
    for k, df in out.items():
        if not df.empty:
            df["ts_local"] = pd.to_datetime(df["ts_local"])
            for c in df.columns:
                if c != "ts_local":
                    df[c] = pd.to_numeric(df[c], errors="coerce")
    return out


def basal_schedule(subject_id: str) -> np.ndarray | None:
    """The 48 half-hourly programmed basal rates, where a study recorded them.

    DCLP5 and PEDAP ship the pump programme at each visit. Loop's schedule comes
    from the rate each temp basal suppressed and arrives through streams().
    Everything else has no schedule on record.
    """
    d = _q("SELECT basal_hh FROM studies.pump_settings WHERE subject_id=%s "
           "AND basal_hh IS NOT NULL ORDER BY visit_seq", (subject_id,))
    if d.empty:
        return None
    arrs = [np.array([np.nan if v is None else float(v) for v in row], dtype=float)
            for row in d.basal_hh if row is not None and len(row) == 48]
    if not arrs:
        return None
    return np.nanmedian(np.vstack(arrs), axis=0)


def entered_isf() -> pd.DataFrame:
    """One entered sensitivity factor and carb ratio per subject, from either source.

    The bolus calculator gives a value per dose, so the median over a subject's
    record is what they were mostly running. The case report forms give a value
    per visit, already time weighted across the day's segments, so the median
    over visits is the same idea.
    """
    wiz = _q("""SELECT subject_id,
                       percentile_cont(0.5) WITHIN GROUP (ORDER BY isf_mgdl_per_u) AS isf,
                       percentile_cont(0.5) WITHIN GROUP (ORDER BY cr_g_per_u)     AS cr,
                       percentile_cont(0.5) WITHIN GROUP (ORDER BY target_low_mgdl) AS target,
                       count(*) AS n_records
                FROM studies.wizard GROUP BY 1""", ())
    ps = _q("""SELECT subject_id,
                      percentile_cont(0.5) WITHIN GROUP (ORDER BY isf_tw_mgdl_per_u) AS isf,
                      percentile_cont(0.5) WITHIN GROUP (ORDER BY cr_tw_g_per_u)     AS cr,
                      NULL::float AS target,
                      count(*) AS n_records
               FROM studies.pump_settings GROUP BY 1""", ())
    wiz["source"], ps["source"] = "calculator", "pump_settings"
    out = pd.concat([wiz, ps], ignore_index=True)
    for c in ("isf", "cr", "target", "n_records"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out["study"] = out.subject_id.str.split(":").str[0]
    return out
