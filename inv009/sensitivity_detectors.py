"""How well do sensitivity detectors work, against the ones AndroidAPS ships?

The previous step showed that adjusting a static sensitivity from recent error
is worth more than the choice of dose exponent. It did not say which detector to
use. This ports the three AndroidAPS ships, puts alternatives beside them, and
scores them all the same way.

The three shipped detectors share a final step and differ only in the statistic
they take over recent deviations:

    basalOff = statistic * (60/5) / sens
    ratio    = 1 + basalOff / maxDailyBasal
    ratio    = clamp(ratio, autosens_min, autosens_max)   defaults 0.7 and 1.2

A deviation is what glucose did over five minutes minus what the insulin acting
in those five minutes should have done. Positive means resistance, and the ratio
divides the sensitivity factor, so a ratio above one gives larger corrections.

    SensitivityAAPSPlugin            median over AutosensPeriod hours, default 24
    SensitivityOref1Plugin           median over 8 and 24 hours, positive
                                     deviations zeroed below 80 mg/dL, and short
                                     histories padded with up to 18 zeros
    SensitivityWeightedAveragePlugin recency-weighted mean over the same period

Against those: a nightly exponentially weighted ratio, a variance-weighted
update in the spirit of a Kalman filter, and a trimmed mean. All are scored on
the overnight fall, with only what a pump could fit, on nights the fit never saw.

    python3 -m inv009.sensitivity_detectors
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys

import numpy as np
import pandas as pd

from . import action, config, db, grid as gridmod, insulin_models as M, tdd as tddmod

MODEL = "oref_6h75"
BIN = config.GRID_MIN
AUTOSENS_MIN, AUTOSENS_MAX = 0.7, 1.2      # DoubleKey.AutosensMin / AutosensMax
MIN_HOURS, MIN_HOURS_FULL = 1.0, 4.0       # Sensitivity.kt
CARB_ABSORPTION_H = 5.0                    # how long a meal invalidates deviations
PAD_MAX = 18                               # oref1's zero padding cap
BASE_CONSTANT = 1800.0                     # ISF = 1800 / seven day average


def _ramp(n_valid: np.ndarray) -> np.ndarray:
    """AbstractSensitivityPlugin.fillResult: fade the ratio in over 1 to 4 hours of data."""
    hours = np.clip(n_valid / 12.0, MIN_HOURS, MIN_HOURS_FULL)
    return (hours - MIN_HOURS) / (MIN_HOURS_FULL - MIN_HOURS)


def _apply(stat: np.ndarray, n_valid: np.ndarray, sens: float, max_daily_basal: float,
           lo: float = AUTOSENS_MIN, hi: float = AUTOSENS_MAX) -> np.ndarray:
    """The step all three shipped detectors share."""
    basal_off = stat * (60.0 / BIN) / sens
    ratio = 1.0 + basal_off / max_daily_basal
    ratio = np.clip(ratio, lo, hi)
    ratio = _ramp(n_valid) * (ratio - 1.0) + 1.0
    return np.where(np.isfinite(ratio), ratio, 1.0)


def _windows_matrix(x: np.ndarray, starts: np.ndarray, back: int) -> np.ndarray:
    """Rows of the `back` values preceding each start, padded with NaN at the edges."""
    idx = starts[:, None] - np.arange(back, 0, -1)[None, :]
    out = np.where(idx >= 0, x[np.clip(idx, 0, len(x) - 1)], np.nan)
    return out


def analyse(job):
    subject_id, study = job
    try:
        streams = db.streams(subject_id)
        g = gridmod.build_grid(streams)
        if g is None:
            return None
        subj = tddmod.subject_level(g)
        if not np.isfinite(subj.get("tdd_u", np.nan)):
            return None
        tw = tddmod.windowed(g)
        n = len(g)
        cgm = g.cgm.to_numpy(float)
        total = g.total_u.to_numpy(float)
        carbs = g.carbs_g.to_numpy(float)
        tdd7 = tw["tdd_7d"].to_numpy(float)

        # The sensitivity the detector is adjusting, and the one it divides by.
        base_isf = BASE_CONSTANT / np.where(tdd7 > 0, tdd7, np.nan)
        sens_ref = float(np.nanmedian(base_isf))
        if not np.isfinite(sens_ref) or sens_ref <= 0:
            return None

        # Profile.getMaxDailyBasal(). No profile is stored here, so the highest
        # hourly basal rate the person actually ran stands in for it, taken high
        # in the distribution rather than at the maximum so one temp basal spike
        # cannot set it.
        hourly = pd.Series(g.basal_u.to_numpy(float)).rolling(12, min_periods=12).sum()
        max_daily_basal = float(np.nanpercentile(hourly.to_numpy(), 99))
        if not np.isfinite(max_daily_basal) or max_daily_basal <= 0.05:
            return None

        # Units acting in each bin, from every dose contributing to it.
        kern = M.kernel(MODEL)
        act_kern = -np.diff(np.append(kern, 0.0))
        acting = np.convolve(total, act_kern)[:n]

        # deviation = what glucose did, minus what the insulin should have done.
        dev = np.full(n, np.nan)
        dev[:-1] = (cgm[1:] - cgm[:-1]) + acting[:-1] * sens_ref

        # AAPS marks a deviation invalid while carbohydrate is still absorbing.
        if (carbs > 0).any():
            absorbing = pd.Series((carbs > 0).astype(float)).rolling(
                int(CARB_ABSORPTION_H * 60 / BIN), min_periods=1).max().to_numpy() > 0
            dev = np.where(absorbing, np.nan, dev)

        # oref1 only: a positive deviation below 80 mg/dL is not treated as resistance.
        dev_oref1 = np.where((cgm < 80) & (dev > 0), 0.0, dev)

        hour = g.ts.dt.hour.to_numpy()
        minute = g.ts.dt.minute.to_numpy()
        ok = np.isin(hour, config.START_HOURS) & (minute == 0)
        h = int(config.HORIZON_MIN / BIN)
        ok &= (np.arange(n) + h) < n
        starts = np.flatnonzero(ok)
        if len(starts) < 60:
            return None

        a_pre, a_in = action.window_action(total, kern, h)
        med_bins = max(int(config.ENDPOINT_MEDIAN_MIN / BIN), 1)
        bg_s = pd.Series(cgm).rolling(med_bins, center=True, min_periods=1).median().to_numpy()

        out = {"subject_id": subject_id, "study": study,
               "t0": g.ts.to_numpy()[starts],
               "bg0": bg_s[starts], "bg_end": bg_s[starts + h],
               "a_pre": a_pre[starts],
               "tdd_7d": tdd7[starts], "tdd_u": subj["tdd_u"],
               "base_isf": base_isf[starts],
               "sens_ref": sens_ref, "max_daily_basal": max_daily_basal}
        out["drop"] = out["bg0"] - out["bg_end"]

        for tag, series, back in (("aaps24", dev, 288), ("aaps8", dev, 96),
                                  ("oref1_8", dev_oref1, 96), ("oref1_24", dev_oref1, 288),
                                  ("wavg24", dev, 288)):
            W = _windows_matrix(series, starts, back)
            cnt = np.sum(np.isfinite(W), axis=1)
            if tag.startswith("oref1"):
                # Pad short histories with zeros, as the plugin does, by adding
                # that many zeros before taking the median.
                pad = np.rint((1.0 - np.clip(cnt / back, 0, 1)) * PAD_MAX).astype(int)
                stat = np.array([
                    np.median(np.concatenate([row[np.isfinite(row)], np.zeros(p)]))
                    if (np.isfinite(row).sum() + p) > 0 else np.nan
                    for row, p in zip(W, pad)])
            elif tag == "wavg24":
                # Weight rises linearly with recency, as the plugin's key
                # difference does, normalised over the valid entries.
                wt = np.arange(1, back + 1, dtype=float)[None, :]
                m = np.isfinite(W)
                den = np.sum(wt * m, axis=1)
                stat = np.where(den > 0, np.nansum(np.where(m, W, 0.0) * wt, axis=1) / den, np.nan)
            else:
                stat = np.nanmedian(W, axis=1)
            # Both clamps are applied to the statistic itself, which is what
            # changing the preference does. Rescaling an already-clamped ratio
            # would not be the same operation and would flatter the wider one.
            out[f"ratio_{tag}"] = _apply(stat, cnt, sens_ref, max_daily_basal)
            out[f"ratio_{tag}_wide"] = _apply(stat, cnt, sens_ref, max_daily_basal,
                                              lo=0.5, hi=1.5)
            out[f"nvalid_{tag}"] = cnt
        return pd.DataFrame(out)
    except Exception:
        return None


def sequential_ratios(d: pd.DataFrame, kind: str, half_life: float,
                      lo: float, hi: float) -> np.ndarray:
    """Detectors that update once per night rather than from a rolling window."""
    a = d.a_pre.to_numpy(float)
    y = d["drop"].to_numpy(float)
    isf = d.base_isf.to_numpy(float)
    off = d.attrs.get("offset", 0.0)
    ratio, var = 1.0, 0.25
    decay = 0.5 ** (1.0 / half_life)
    out = np.empty(len(d))
    for i in range(len(d)):
        out[i] = ratio
        if not (a[i] > 0.3 and np.isfinite(isf[i]) and np.isfinite(y[i])):
            continue
        # The ratio that would have made tonight right. Sensitivity divides, so a
        # night that fell short of prediction wants a ratio above one.
        pred_drop = isf[i] * a[i]
        if pred_drop <= 0:
            continue
        want = pred_drop / max(y[i] - off, 1e-6) if (y[i] - off) > 0 else hi
        want = float(np.clip(want, lo, hi))
        if kind == "ewma":
            ratio = decay * ratio + (1 - decay) * want
        elif kind == "kalman":
            # Weight each night by how much insulin was acting: a night with
            # little insulin says little about sensitivity.
            obs_var = 1.0 / max(a[i], 0.1)
            gain = var / (var + obs_var)
            ratio = ratio + gain * (want - ratio)
            var = (1 - gain) * var + (1 - decay) * 0.05
        ratio = float(np.clip(ratio, lo, hi))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    config.ensure_dirs()
    jobs = [(r.subject_id, r.study_name)
            for study in config.COHORTS
            for r in db.subjects(study).itertuples()]
    if args.limit:
        jobs = jobs[:args.limit]
    print(f"{len(jobs)} people, {config.WORKERS} workers")
    frames = []
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        for i, out in enumerate(pool.imap_unordered(analyse, jobs, chunksize=4), 1):
            if out is not None:
                frames.append(out)
            if i % 300 == 0:
                print(f"  {i}/{len(jobs)}", flush=True)
    D = pd.concat(frames, ignore_index=True)
    D = D.dropna(subset=["drop", "a_pre", "base_isf"])
    D.to_parquet(config.RESULTS / "inv009_detectors.parquet", index=False)
    print(f"{len(D):,} nights from {D.subject_id.nunique():,} people")

    shipped = ["aaps24", "aaps8", "oref1_8", "oref1_24", "wavg24"]
    alts = [("ewma", 1.0), ("ewma", 3.0), ("ewma", 7.0), ("kalman", 3.0)]
    rows = []
    for sid, d in D.groupby("subject_id", sort=False):
        d = d.sort_values("t0")
        if len(d) < 80:
            continue
        n_tr = int(len(d) * 0.7)
        if len(d) - n_tr < 20:
            continue
        tr, te = d.iloc[:n_tr], d.iloc[n_tr:]
        rec = {"subject_id": sid, "study": d.study.iloc[0], "n": len(d)}

        def score(ratio_all):
            r_tr = ratio_all[:n_tr]
            r_te = ratio_all[n_tr:]
            p_tr = (tr.base_isf.to_numpy(float) / r_tr) * tr.a_pre.to_numpy(float)
            off = float(np.mean(tr["drop"].to_numpy(float) - p_tr))
            p_te = (te.base_isf.to_numpy(float) / r_te) * te.a_pre.to_numpy(float)
            return float(np.median(np.abs(te["drop"].to_numpy(float) - (p_te + off))))

        rec["static"] = score(np.ones(len(d)))
        for tag in shipped:
            rec[tag] = score(d[f"ratio_{tag}"].to_numpy(float))
            rec[f"{tag}_wide"] = score(d[f"ratio_{tag}_wide"].to_numpy(float))
        # Sequential detectors need the offset from a static pass first.
        p_tr0 = tr.base_isf.to_numpy(float) * tr.a_pre.to_numpy(float)
        d.attrs["offset"] = float(np.mean(tr["drop"].to_numpy(float) - p_tr0))
        for kind, hl in alts:
            for lo, hi in ((0.7, 1.2), (0.5, 1.5)):
                key = f"{kind}{hl:g}_{lo}-{hi}"
                rec[key] = score(sequential_ratios(d, kind, hl, lo, hi))
        rows.append(rec)
    R = pd.DataFrame(rows)
    R.to_parquet(config.RESULTS / "inv009_detector_scores.parquet", index=False)

    cols = [c for c in R.columns if c not in ("subject_id", "study", "n")]
    med = {c: float(R[c].median()) for c in cols}
    base = med["static"]
    print(f"\n{len(R)} people scored. Error in the predicted overnight fall, mg/dL.\n")
    print(f"  {'detector':>26s} {'error':>8s} {'vs static':>10s} {'best for':>9s}")
    best_counts = R[cols].idxmin(axis=1).value_counts().to_dict()
    for c, v in sorted(med.items(), key=lambda kv: kv[1]):
        print(f"  {c:>26s} {v:8.2f} {v - base:>+10.2f} {best_counts.get(c, 0):9d}")
    res = dict(n_people=int(len(R)), n_nights=int(len(D)), medians=med,
               best_counts={k: int(v) for k, v in best_counts.items()},
               shipped=shipped)
    (config.RESULTS / "inv009_detectors.json").write_text(json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
