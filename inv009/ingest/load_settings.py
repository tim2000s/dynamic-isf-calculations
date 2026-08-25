#!/usr/bin/env python3
"""Load pump settings and bolus calculator records that BabelBetes does not extract.

BabelBetes emits five event streams and stops there. That is the right scope for
it, but it leaves out the number a study of insulin sensitivity most wants: the
sensitivity factor the person actually had entered. Those numbers sit unloaded
in the raw archives, in two different shapes.

    Loop        LOOPDeviceWizard.txt              -> studies.wizard
    REPLACE-BG  HDeviceWizard.txt                 -> studies.wizard
    DCLP5       InsulinPumpSettings.txt           -> studies.pump_settings
    PEDAP       PEDAPInsulinDeliveryDetails.txt   -> studies.pump_settings
    Loop        LOOPDeviceBasal{1,2,3}.txt        -> studies.basal_sched

Three things about this data are worth knowing before using it.

Units. Loop and REPLACE-BG store the calculator fields in mmol/L whatever the
person saw on screen. The conversion is decided per subject on the median, not
per row, so one stray value cannot flip a subject's units.

Timestamps. The wizard files are put on the same clock as the event streams, by
the same rule BabelBetes used, so a wizard row lines up with studies.cgm. Rows
outside a subject's loaded glucose record are dropped, which is what enforces
the study windows without re-deriving each study's exclusion rules.

Scheduled basal. Loop computes insulin on board net of scheduled basal, so
reproducing what the app believed needs the schedule and not just what was
delivered. The only record of it is the rate each temp basal suppressed. It is
observed only while a temp was running and it is a step function, so it is
stored as changes and should be forward filled.

Usage:
    python3 inv009/ingest/load_settings.py               # everything
    python3 inv009/ingest/load_settings.py --source loop_wizard
    python3 inv009/ingest/load_settings.py --dry-run     # parse, report, load nothing
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path.home() / "babelbetes/data/raw"
PSQL = "/opt/homebrew/opt/postgresql@17/bin/psql"
DB = os.environ.get("BABELBETES_DB", "oref")
VENV_PY = Path.home() / ".venvs/babelbetes/bin/python"
MMOL_TO_MGDL = 18.018
# A subject whose median entered sensitivity is below this is working in mmol/L.
# Real mg/dL sensitivities start around 10 and mmol ones top out around 10, but
# the gap in practice is wide: the cohort medians cluster near 2 and near 45.
MMOL_ISF_CUTOFF = 20.0

WIZARD_COLS = ("subject_id, ts_local, bolus_rec_id, bg_input_mgdl, carb_input_g, iob_u, "
               "cr_g_per_u, isf_mgdl_per_u, target_low_mgdl, target_high_mgdl, "
               "rec_correction_u, rec_net_u")
SETTINGS_COLS = ("subject_id, visit_seq, visit_label, therapy_dt, tdd_reported_u, "
                 "basal_reported_u, basal_hh, basal_total_u, isf_segments, "
                 "cr_segments, isf_tw_mgdl_per_u, cr_tw_g_per_u")
SCHED_COLS = "subject_id, ts_local, sched_rate_u_hr"


def psql(sql: str) -> str:
    return subprocess.run([PSQL, "-d", DB, "-v", "ON_ERROR_STOP=1", "-Atc", sql],
                          check=True, text=True, capture_output=True).stdout.strip()


def _num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _convert_units(df: pd.DataFrame) -> pd.DataFrame:
    """Scale the glucose-dimensioned columns for subjects recorded in mmol/L.

    Decided per subject on the median sensitivity. Carb ratio is grams per unit
    in both unit systems and the recommendation columns are units of insulin, so
    neither is touched.
    """
    # A zero or negative sensitivity or carb ratio is a skipped field, not a
    # setting, and one of them would drag a subject's median across the units
    # cutoff. Clear them first.
    for col in ("isf_mgdl_per_u", "cr_g_per_u", "target_low_mgdl", "target_high_mgdl"):
        df.loc[df[col] <= 0, col] = np.nan
    med = df.groupby("subject_id")["isf_mgdl_per_u"].transform("median")
    is_mmol = med < MMOL_ISF_CUTOFF
    for col in ("isf_mgdl_per_u", "bg_input_mgdl", "target_low_mgdl", "target_high_mgdl"):
        df.loc[is_mmol, col] = df.loc[is_mmol, col] * MMOL_TO_MGDL
    n_sub = df.loc[is_mmol, "subject_id"].nunique()
    print(f"    mmol/L subjects converted: {n_sub} of {df.subject_id.nunique()}")
    return df


def loop_wizard() -> pd.DataFrame:
    """Loop's bolus calculator. ts_local = UTC plus the roster's per-subject offset."""
    path = RAW / "Loop/Data Tables/LOOPDeviceWizard.txt"
    df = pd.read_csv(path, sep="|", low_memory=False,
                     usecols=["PtID", "UTCDtTm", "RecommendedCorrection", "RecommendedNet",
                              "BgInput", "CarbInput", "InsulinOnBoard", "InsulinCarbRatio",
                              "InsulinSensitivity", "BgTargetLow", "BgTargetHigh",
                              "BgTargetTarget", "BolusRecID"])
    roster = pd.read_csv(RAW / "Loop/Data Tables/PtRoster.txt", sep="|",
                         usecols=["PtID", "PtTimezoneOffset"])
    df = df.merge(roster, on="PtID", how="left")
    df = df.dropna(subset=["UTCDtTm", "PtTimezoneOffset"])
    ts = (pd.to_datetime(df.UTCDtTm, format="%Y-%m-%d %H:%M:%S", errors="coerce")
          + pd.to_timedelta(df.PtTimezoneOffset, unit="hour"))

    # Loop records a target range in either of two ways depending on the version
    # that uploaded, so fall back to the single target when the pair is absent.
    lo, hi = _num(df.BgTargetLow), _num(df.BgTargetHigh)
    tgt = _num(df.BgTargetTarget)
    out = pd.DataFrame({
        "subject_id": "Loop:" + df.PtID.astype(str),
        "ts_local": ts,
        "bolus_rec_id": df.BolusRecID.fillna("").astype(str).str.replace(r"\.0$", "", regex=True),
        "bg_input_mgdl": _num(df.BgInput),
        "carb_input_g": _num(df.CarbInput),
        "iob_u": _num(df.InsulinOnBoard),
        "cr_g_per_u": _num(df.InsulinCarbRatio),
        "isf_mgdl_per_u": _num(df.InsulinSensitivity),
        "target_low_mgdl": lo.fillna(tgt),
        "target_high_mgdl": hi.fillna(tgt),
        "rec_correction_u": _num(df.RecommendedCorrection),
        "rec_net_u": _num(df.RecommendedNet),
    })
    return _convert_units(out.dropna(subset=["ts_local"]))


def replacebg_wizard() -> pd.DataFrame:
    """REPLACE-BG's bolus calculator. Days are counted from a 2015-01-01 epoch."""
    path = RAW / "REPLACE-BG/Data Tables/HDeviceWizard.txt"
    df = pd.read_csv(path, sep="|", low_memory=False,
                     usecols=["PtId", "DeviceDtTmDaysFromEnroll", "DeviceTm",
                              "RecommendedCorrection", "RecommendedNet", "BgInput",
                              "CarbInput", "InsulinCarbRatio", "InsulinSensitivity",
                              "BgTargetLow", "BgTargetHigh", "BgTargetTarget",
                              "ParentHDeviceBolusID"])
    days = _num(df.DeviceDtTmDaysFromEnroll)
    df = df[days >= 0].copy()          # pre-enrolment uploads are outside the loaded record
    days = days[days >= 0]
    ts = (pd.Timestamp("2015-01-01")
          + pd.to_timedelta(days, unit="D")
          + pd.to_timedelta(df.DeviceTm.astype(str), errors="coerce"))

    lo, hi = _num(df.BgTargetLow), _num(df.BgTargetHigh)
    tgt = _num(df.BgTargetTarget)
    out = pd.DataFrame({
        "subject_id": "ReplaceBG:" + df.PtId.astype(str),
        "ts_local": ts,
        "bolus_rec_id": df.ParentHDeviceBolusID.fillna("").astype(str).str.replace(r"\.0$", "", regex=True),
        "bg_input_mgdl": _num(df.BgInput),
        "carb_input_g": _num(df.CarbInput),
        "iob_u": np.nan,               # the column exists but this study did not populate it usefully
        "cr_g_per_u": _num(df.InsulinCarbRatio),
        "isf_mgdl_per_u": _num(df.InsulinSensitivity),
        "target_low_mgdl": lo.fillna(tgt),
        "target_high_mgdl": hi.fillna(tgt),
        "rec_correction_u": _num(df.RecommendedCorrection),
        "rec_net_u": _num(df.RecommendedNet),
    })
    # A keyed glucose of zero means the field was skipped, not a glucose of zero.
    out.loc[out.bg_input_mgdl <= 0, "bg_input_mgdl"] = np.nan
    return _convert_units(out.dropna(subset=["ts_local"]))


def _parse_settings(path: Path, study: str, pid_col: str, visit_col: str | None,
                    tdd_col: str, basal7d_col: str, dt_col: str) -> pd.DataFrame:
    """Shared parser for the two CRF pump-settings files.

    The basal schedule is sparse by convention: a half hour is written only when
    the rate changes, and holds until the next one written. Reading the blanks as
    zero would halve everybody's basal, so they are forward filled, wrapping from
    the end of the day when midnight itself is blank.
    """
    df = pd.read_csv(path, sep="|", low_memory=False, dtype=str)
    hh_cols = [f"InsBasal{h:02d}{m:02d}" for h in range(24) for m in (0, 30)]
    basal = df[hh_cols].apply(_num)
    filled = basal.ffill(axis=1)
    # Midnight blank means the rate carried over from the end of the previous day.
    last = filled.iloc[:, -1]
    for c in hh_cols:
        filled[c] = filled[c].fillna(last)

    rows = []
    seq = df.groupby(pid_col).cumcount() + 1
    for i, (_, r) in enumerate(df.iterrows()):
        isf_seg, cr_seg = [], []
        for k in range(1, 11):
            sh, sm = _num(pd.Series([r.get(f"InsBolusStart{k}Hr"), r.get(f"InsBolusStart{k}Min")]))
            eh, em = _num(pd.Series([r.get(f"InsBolusEnd{k}Hr"), r.get(f"InsBolusEnd{k}Min")]))
            cf = pd.to_numeric(r.get(f"CorrFactor{k}"), errors="coerce")
            cho = pd.to_numeric(r.get(f"CHORatio{k}"), errors="coerce")
            if np.isnan(sh) or np.isnan(eh):
                continue
            start = int(sh) * 60 + int(0 if np.isnan(sm) else sm)
            end = int(eh) * 60 + int(0 if np.isnan(em) else em)
            if end <= start:
                end += 24 * 60          # a segment that wraps past midnight
            if not np.isnan(cf):
                isf_seg.append({"start_min": start, "end_min": end, "value": float(cf)})
            if not np.isnan(cho):
                cr_seg.append({"start_min": start, "end_min": end, "value": float(cho)})

        def tw(segs):
            if not segs:
                return np.nan
            w = sum(s["end_min"] - s["start_min"] for s in segs)
            return sum(s["value"] * (s["end_min"] - s["start_min"]) for s in segs) / w if w else np.nan

        hh = filled.iloc[i].to_numpy(dtype=float)
        rows.append({
            "subject_id": f"{study}:{r[pid_col]}",
            "visit_seq": int(seq.iloc[i]),
            "visit_label": (r.get(visit_col) if visit_col else None),
            "therapy_dt": pd.to_datetime(r.get(dt_col), errors="coerce"),
            "tdd_reported_u": pd.to_numeric(r.get(tdd_col), errors="coerce"),
            "basal_reported_u": pd.to_numeric(r.get(basal7d_col), errors="coerce"),
            "basal_hh": hh,
            "basal_total_u": float(np.nansum(hh) / 2.0) if np.isfinite(hh).any() else np.nan,
            "isf_segments": isf_seg,
            "cr_segments": cr_seg,
            "isf_tw_mgdl_per_u": tw(isf_seg),
            "cr_tw_g_per_u": tw(cr_seg),
        })
    return pd.DataFrame(rows)


def dclp5_settings() -> pd.DataFrame:
    return _parse_settings(RAW / "DCLP5/InsulinPumpSettings.txt", "DCLP5",
                           pid_col="PtID", visit_col=None, tdd_col="CurrTotInsDaily",
                           basal7d_col="TotBasalPreced7Days", dt_col="InsTherapyDt")


def pedap_settings() -> pd.DataFrame:
    return _parse_settings(RAW / "PEDAP/Data Files/PEDAPInsulinDeliveryDetails.txt", "PEDAP",
                           pid_col="PtID", visit_col="Visit", tdd_col="PumpTotInsDaily",
                           basal7d_col="PumpTotBasalPreced7Days", dt_col="PumpInsTherapyDt")


def loop_basal_sched(out_csv: Path, files: list[str] | None = None) -> int:
    """Loop's scheduled basal, from the rate each temp basal suppressed.

    Six and a half gigabytes of raw text describing a step function that changes
    a handful of times a day, so DuckDB reads it under a memory ceiling, spills
    to scratch rather than growing without bound, and emits only the rows where
    the scheduled rate actually changed.
    """
    if files is None:
        files = [str(RAW / f"Loop/Data Tables/LOOPDeviceBasal{i}.txt") for i in (1, 2, 3)]
    spill = Path(tempfile.gettempdir()) / "duckdb_spill"
    spill.mkdir(exist_ok=True)
    # The roster is 919 rows and trips DuckDB's dialect sniffer, so it is read
    # here and inlined rather than being a second thing that can fail.
    roster = pd.read_csv(RAW / "Loop/Data Tables/PtRoster.txt", sep="|",
                         usecols=["PtID", "PtTimezoneOffset"]).dropna()
    tz_values = ", ".join(f"('{int(r.PtID)}', {float(r.PtTimezoneOffset)})"
                          for r in roster.itertuples())
    sql = f"""
    SET memory_limit='8GB';
    SET threads TO 4;
    SET temp_directory='{spill}';
    COPY (
      WITH raw AS (
        SELECT CAST(PtID AS VARCHAR) AS pid,
               strptime(UTCDtTm, '%Y-%m-%d %H:%M:%S') AS utc_ts,
               CAST(SuprRate AS DOUBLE) AS rate
        FROM read_csv({files!r}, delim='|', header=true, all_varchar=true,
                      ignore_errors=true)
        WHERE SuprBasalType = 'scheduled' AND SuprRate IS NOT NULL AND UTCDtTm IS NOT NULL
      ),
      roster(pid, tz) AS (VALUES {tz_values}),
      shifted AS (
        SELECT 'Loop:' || raw.pid AS subject_id,
               raw.utc_ts + INTERVAL (CAST(roster.tz * 60 AS BIGINT)) MINUTE AS ts_local,
               raw.rate AS sched_rate_u_hr
        FROM raw JOIN roster USING (pid)
        WHERE raw.utc_ts IS NOT NULL
      ),
      changes AS (
        SELECT *, lag(sched_rate_u_hr) OVER (PARTITION BY subject_id ORDER BY ts_local) AS prev
        FROM shifted
      )
      SELECT subject_id, ts_local, sched_rate_u_hr FROM changes
      WHERE prev IS NULL OR prev IS DISTINCT FROM sched_rate_u_hr
      ORDER BY subject_id, ts_local
    ) TO '{out_csv}' (FORMAT csv, HEADER false);
    """
    script = out_csv.with_suffix(".sql")
    script.write_text(sql)
    subprocess.run([str(VENV_PY), "-c",
                    "import duckdb,sys;duckdb.connect().execute(open(sys.argv[1]).read())",
                    str(script)], check=True)
    return sum(1 for _ in open(out_csv))


def copy_into(table: str, cols: str, csv_path: Path, conflict_cols: str,
              force_not_null: str = "") -> None:
    """Stage the CSV then upsert, so a re-run replaces rather than collides.

    force_not_null is for the text key columns: COPY reads an empty field as
    NULL, and a record with no bolus link still needs a key.
    """
    stage = table.split(".")[-1] + "_stage"
    psql(f"DROP TABLE IF EXISTS studies.{stage}")
    psql(f"CREATE UNLOGGED TABLE studies.{stage} (LIKE {table} INCLUDING DEFAULTS)")
    opts = "FORMAT csv" + (f", FORCE_NOT_NULL ({force_not_null})" if force_not_null else "")
    subprocess.run([PSQL, "-d", DB, "-v", "ON_ERROR_STOP=1", "-q",
                    "-c", f"\\copy studies.{stage} ({cols}) FROM '{csv_path}' ({opts})"],
                   check=True)
    updates = ", ".join(f"{c.strip()} = EXCLUDED.{c.strip()}"
                        for c in cols.split(",") if c.strip() not in conflict_cols)
    psql(f"""INSERT INTO {table} ({cols})
             SELECT DISTINCT ON ({conflict_cols}) {cols} FROM studies.{stage}
             ORDER BY {conflict_cols}
             ON CONFLICT ({conflict_cols}) DO UPDATE SET {updates}""")
    psql(f"DROP TABLE studies.{stage}")
    psql(f"ANALYZE {table}")


def write_csv(df: pd.DataFrame, path: Path) -> None:
    out = df.copy()
    for c in out.columns:
        if out[c].dtype == object and len(out) and isinstance(out[c].dropna().iloc[0] if out[c].notna().any() else None, (list, np.ndarray)):
            if c == "basal_hh":
                # Postgres array literals need the word NULL for a missing
                # element, and a row whose whole schedule is blank (the CRF's
                # "settings unchanged" case) is no schedule at all.
                out[c] = out[c].map(
                    lambda a: "" if not np.isfinite(np.asarray(a, dtype=float)).any()
                    else "{" + ",".join("NULL" if not np.isfinite(v) else f"{v:g}" for v in a) + "}")
            else:
                out[c] = out[c].map(json.dumps)
    out.to_csv(path, index=False, header=False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default="all",
                    choices=["all", "loop_wizard", "replacebg_wizard", "dclp5", "pedap",
                             "loop_basal_sched"])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    want = (lambda s: args.source in ("all", s))
    tmpdir = Path(tempfile.mkdtemp(prefix="inv009_settings_"))

    if want("loop_wizard") or want("replacebg_wizard"):
        frames = []
        if want("loop_wizard"):
            print("  Loop bolus calculator")
            frames.append(loop_wizard())
        if want("replacebg_wizard"):
            print("  REPLACE-BG bolus calculator")
            frames.append(replacebg_wizard())
        wiz = pd.concat(frames, ignore_index=True)
        print(f"    {len(wiz):,} rows, {wiz.subject_id.nunique()} subjects, "
              f"{wiz.iob_u.notna().sum():,} with insulin on board")
        if not args.dry_run:
            p = tmpdir / "wizard.csv"
            write_csv(wiz[[c.strip() for c in WIZARD_COLS.split(",")]], p)
            copy_into("studies.wizard", WIZARD_COLS, p, "subject_id, ts_local, bolus_rec_id",
                      force_not_null="bolus_rec_id")

    if want("dclp5") or want("pedap"):
        frames = []
        if want("dclp5"):
            print("  DCLP5 pump settings")
            frames.append(dclp5_settings())
        if want("pedap"):
            print("  PEDAP pump settings")
            frames.append(pedap_settings())
        st = pd.concat(frames, ignore_index=True)
        print(f"    {len(st):,} visit records, {st.subject_id.nunique()} subjects, "
              f"{st.isf_tw_mgdl_per_u.notna().sum():,} with a correction factor")
        if not args.dry_run:
            p = tmpdir / "settings.csv"
            write_csv(st[[c.strip() for c in SETTINGS_COLS.split(",")]], p)
            copy_into("studies.pump_settings", SETTINGS_COLS, p, "subject_id, visit_seq")

    if want("loop_basal_sched"):
        print("  Loop scheduled basal (this reads 6.7 GB, allow a few minutes)")
        if not args.dry_run:
            p = tmpdir / "sched.csv"
            n = loop_basal_sched(p)
            print(f"    {n:,} rate changes")
            psql("TRUNCATE studies.basal_sched")
            subprocess.run([PSQL, "-d", DB, "-v", "ON_ERROR_STOP=1", "-q", "-c",
                            f"\\copy studies.basal_sched ({SCHED_COLS}) FROM '{p}' (FORMAT csv)"],
                           check=True)
            psql("ANALYZE studies.basal_sched")

    if not args.dry_run and want("loop_wizard") and want("replacebg_wizard"):
        # Keep the calculator records consistent with the event streams: a row
        # for a subject BabelBetes excluded, or from outside the glucose record
        # it loaded, has nothing to be analysed against. This enforces each
        # study's own inclusion window without re-deriving its exclusion rules.
        n0 = int(psql("SELECT count(*) FROM studies.wizard"))
        psql("DELETE FROM studies.wizard w WHERE NOT EXISTS "
             "(SELECT 1 FROM studies.subject s WHERE s.subject_id = w.subject_id)")
        psql("""DELETE FROM studies.wizard w USING (
                    SELECT subject_id, min(ts_local) t0, max(ts_local) t1
                    FROM studies.cgm GROUP BY 1) s
                WHERE w.subject_id = s.subject_id
                  AND (w.ts_local < s.t0 OR w.ts_local > s.t1)""")
        n1 = int(psql("SELECT count(*) FROM studies.wizard"))
        print(f"    trimmed to the loaded record: {n0 - n1:,} rows dropped")

    if not args.dry_run:
        for t in ("wizard", "pump_settings", "basal_sched"):
            n = psql(f"SELECT count(*) FROM studies.{t}")
            s = psql(f"SELECT count(DISTINCT subject_id) FROM studies.{t}")
            print(f"  studies.{t}: {int(n):,} rows, {s} subjects")
    return 0


if __name__ == "__main__":
    sys.exit(main())
