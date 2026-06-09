#!/usr/bin/env python3
"""Test #2: model-INDEPENDENT marginal dose-response of correction insulin vs glucose.

The realised_isf used elsewhere is derived from the loop's own IOB prediction (reason_iobpredbg), so
it inherits any bias in the loop's insulin model. This replicates the Diabeloop method instead: take
the CORRECTION insulin the algorithm actually delivered (sug_smb_units, summed over the window) and
ask how much glucose actually dropped per unit — by glucose level — WITHOUT using the IOB prediction.

Per overnight (23:00-02:00), carb-screened, 4-hour window:
    drop            = cgm_start − cgm_end
    correction      = Σ sug_smb_units over the window (delivered correction insulin, model-free)
    iob_start       = iob_iob at the decision (background insulin already on board)
    pre_slope       = 30-min backward slope (entry trajectory; mean-reversion control)

Within each glucose band we regress  drop ~ b0 + b_corr·correction + b_iob·iob_start. The coefficient
**b_corr is the effective ISF** — extra mg/dL drop per unit of delivered correction, holding starting
IOB fixed. Insulin-INDEPENDENT clearance (renal, glucose effectiveness) scales with BG, not with the
correction dose, so it lands in b0, partly shielding b_corr from that confound.

Verdict: if b_corr FALLS with glucose ⇒ a unit of insulin does less at high BG ⇒ genuine resistance ⇒
the power law / Tim is right at the top end. If b_corr is flat/rising ⇒ no acute resistance overnight.

Single-process, one user at a time (no mp.Pool — the crash hazard does not apply). v6 excluded.
Output: results/dose_response_db.{json,md}, charts/inv008/fig_dose_response_db.png
Run: python -m inv008.dose_response_db
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psycopg2

from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
TBL = {"v5_trio": "oref_v5", "v7_oref0": "oref_v7"}
SMBC = {"oref_v5": "sug_smb_units", "oref_v7": "sug_smb_units"}
START_HOURS = {23, 0, 1, 2}
HZ, TOL, RISE_MAX = 4 * 3600, 300, 2.0
BANDS = [(100, 120), (120, 145), (145, 175), (175, 205), (205, 260)]
LBL = [f"{a}-{b}" for a, b in BANDS]
CTR = [(a + b) / 2 for a, b in BANDS]


def user_windows(user_id, table):
    conn = psycopg2.connect("dbname=oref")
    try:
        d = pd.read_sql(f"""SELECT ts_relative_sec, hour, cgm_mgdl, iob_iob,
                                   {SMBC[table]} AS smb
                            FROM {table} WHERE user_id=%s AND cgm_mgdl IS NOT NULL
                            ORDER BY ts_relative_sec""", conn, params=(user_id,))
    finally:
        conn.close()
    if len(d) < 500:
        return None
    ts = d.ts_relative_sec.values.astype(float)
    bg = d.cgm_mgdl.values.astype(float)
    hr = d.hour.values.astype(float)
    iob = pd.to_numeric(d.iob_iob, errors="coerce").values.astype(float)
    smb = pd.to_numeric(d.smb, errors="coerce").fillna(0).values.astype(float)
    n = len(d)
    csum = np.concatenate([[0.0], np.cumsum(smb)])               # for window sums
    end4 = np.searchsorted(ts, ts + HZ); end4 = np.where(end4 < n, end4, n - 1)
    oke = np.abs(ts[end4] - (ts + HZ)) <= TOL
    p15 = np.searchsorted(ts, ts + 900); p15 = np.where(p15 < n, p15, n - 1)
    sl = np.where(np.abs(ts[p15] - (ts + 900)) <= TOL, (bg[p15] - bg) / 3.0, 0.0)
    pre = np.searchsorted(ts, ts - 1800); pre = np.clip(pre, 0, n - 1)
    pre_slope = np.where(np.abs(ts[pre] - (ts - 1800)) <= TOL, (bg - bg[pre]) / 6.0, np.nan)
    rows = []
    for i in range(n):
        if hr[i] not in START_HOURS or not oke[i] or bg[i] < 80 or bg[i] > 260:
            continue
        j = end4[i]
        if np.nanmax(sl[i:j + 1]) > RISE_MAX:                    # carb screen
            continue
        corr = csum[j + 1] - csum[i]                              # delivered correction over window
        rows.append((bg[i], bg[i] - bg[j], corr, iob[i], pre_slope[i]))
    if len(rows) < 40:
        return None
    df = pd.DataFrame(rows, columns=["bg", "drop", "corr", "iob", "pre_slope"])
    df["user"] = user_id
    return df


def band_isf(d):
    """Per glucose band: b_corr (effective ISF) from drop ~ corr + iob, with bootstrap CI."""
    out = []
    for a, b in BANDS:
        m = (d.bg >= a) & (d.bg < b) & d["corr"].notna() & d["iob"].notna()
        g = d[m]
        if len(g) < 200:
            out.append({"band": f"{a}-{b}", "n": int(len(g)), "isf": None, "ci": None})
            continue
        X = np.column_stack([np.ones(len(g)), g["corr"].values, g["iob"].values])
        fit = ec.ols(g["drop"].values, X)
        if fit is None:
            out.append({"band": f"{a}-{b}", "n": int(len(g)), "isf": None, "ci": None})
            continue
        beta, se, _ = fit
        out.append({"band": f"{a}-{b}", "n": int(len(g)), "isf": round(float(beta[1]), 1),
                    "ci": [round(float(beta[1] - 1.96 * se[1]), 1), round(float(beta[1] + 1.96 * se[1]), 1)],
                    "median_corr": round(float(g["corr"].median()), 2)})
    return out


def main():
    coh = {r["user_id"]: r for r in json.load(open(config.ROOT / "canonical_cohort.json"))}
    jobs = [(u, TBL[r["cohort"]]) for u, r in coh.items()
            if r.get("cohort") in TBL and r.get("isf")]
    print(f"dose-response: {len(jobs)} v5/v7 users, single-process")
    parts = []
    for k, (u, table) in enumerate(jobs):
        df = user_windows(u, table)
        if df is not None:
            parts.append(df)
        if k % 25 == 0:
            print(f"  {k}/{len(jobs)} users, {sum(len(p) for p in parts):,} windows")
    D = pd.concat(parts, ignore_index=True)
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)

    allb = band_isf(D)
    clean = band_isf(D[D.pre_slope > -1.0])                       # non-falling entry (mean-rev removed)
    falling = band_isf(D[D.pre_slope <= -1.0])

    def slope_k(curve):
        pts = [(CTR[i], curve[i]["isf"]) for i in range(len(curve)) if curve[i]["isf"] and CTR[i] >= 120]
        if len(pts) < 3:
            return None
        x = np.log(np.array([p[0] for p in pts]) / 100.0)
        y = np.log(np.clip([p[1] for p in pts], 1, None))
        b = np.polyfit(x, y, 1)[0]
        return round(-float(b), 2)                                # k>0 ⇒ ISF falls with BG = resistance

    summary = {
        "n_users": int(D.user.nunique()), "n_windows": int(len(D)),
        "method": "model-FREE: effective ISF = slope of drop on delivered correction (sug_smb_units), "
                  "controlling starting IOB, per glucose band. Renal/mass-action clearance ∝ BG lands "
                  "in the intercept, not this slope.",
        "effective_isf_by_bg_all": allb,
        "effective_isf_by_bg_clean_nonfalling": clean,
        "effective_isf_by_bg_falling": falling,
        "implied_k": {"all": slope_k(allb), "clean": slope_k(clean), "falling": slope_k(falling),
                      "note": "k>0 = effective ISF FALLS with glucose = resistance = power-law direction"},
        "verdict": None,
    }
    kc = slope_k(clean)
    # identification check: a real ISF is ~20-100 mg/dL/U. Negative/tiny slopes ⇒ the dose is
    # endogenous (loop doses MORE when glucose responds LESS) → confounding-by-indication → the
    # dose-response slope does NOT identify ISF. Detect and refuse to interpret.
    isf_vals = [c["isf"] for c in clean if c["isf"] is not None]
    plausible = isf_vals and np.median(isf_vals) >= 5.0
    summary["identification_valid"] = bool(plausible)
    summary["reverse_causation_evidence"] = ("effective-ISF slopes are negative/implausible "
        f"(median {np.median(isf_vals):.1f} mg/dL/U); the loop's correction dose is reactive, so the "
        "dose→drop relationship is confounded by indication and cannot identify ISF") if not plausible else None
    summary["verdict"] = (("RESISTANCE confirmed" if (kc and kc > 0.3) else
                           "no acute resistance — effective ISF flat/rising with glucose") + " (model-free)"
                          if plausible else
                          "INCONCLUSIVE — dose-response confounded by reactive dosing (negative/implausible "
                          "ISF). Closed-loop doses are endogenous; this method cannot identify ISF.")
    (OUT / "dose_response_db.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for curve, col, lab in [(allb, "#999", "all"), (clean, "#2ca02c", "clean (non-falling)"),
                            (falling, "#1f77b4", "falling")]:
        xs = [CTR[i] for i in range(len(curve)) if curve[i]["isf"]]
        ys = [curve[i]["isf"] for i in range(len(curve)) if curve[i]["isf"]]
        ax[0].plot(xs, ys, "o-", color=col, lw=2, label=lab)
    ax[0].set_xlabel("glucose (mg/dL)"); ax[0].set_ylabel("effective ISF = mg/dL drop per U correction")
    ax[0].set_title("Model-free dose-response: effective ISF vs glucose\n(falls ⇒ resistance / power law)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    # clean with CIs
    xs = [CTR[i] for i in range(len(clean)) if clean[i]["isf"]]
    ys = [clean[i]["isf"] for i in range(len(clean)) if clean[i]["isf"]]
    lo = [clean[i]["ci"][0] for i in range(len(clean)) if clean[i]["isf"]]
    hi = [clean[i]["ci"][1] for i in range(len(clean)) if clean[i]["isf"]]
    ax[1].plot(xs, ys, "o-", color="#2ca02c", lw=2.5, label=f"clean (k={kc})")
    ax[1].fill_between(xs, lo, hi, color="#2ca02c", alpha=0.2, label="95% CI")
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("effective ISF (mg/dL per U)")
    ax[1].set_title("Clean (mean-reversion removed), with CI"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_dose_response_db.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Test #2: model-free correction-dose response vs glucose\n",
          f"{summary['n_users']} users, {len(D):,} windows. Effective ISF = slope of (drop) on "
          "(delivered correction insulin), controlling starting IOB, per glucose band. **Model-free** "
          "(does not use the loop's IOB prediction). Renal/mass-action clearance ∝ BG → intercept, not "
          "this slope.\n",
          "## Effective ISF (mg/dL per U of correction) by glucose\n",
          "| BG band | all | clean (non-falling) | falling | n(clean) |", "|---|---|---|---|---|"]
    for i, (a, b) in enumerate(BANDS):
        md.append(f"| {a}-{b} | {allb[i]['isf']} | {clean[i]['isf']} | {falling[i]['isf']} | {clean[i]['n']:,} |")
    md += [f"\n**Implied k** (k>0 = ISF falls with BG = resistance): all {summary['implied_k']['all']}, "
           f"**clean {kc}**, falling {summary['implied_k']['falling']}.\n",
           f"**Verdict: {summary['verdict']}.**\n",
           "![dose response](charts/inv008/fig_dose_response_db.png)\n",
           "*Caveat: corrections are larger at high BG (selection); if the loop compensates for poor "
           "response with more insulin, b_corr could be biased. Manual boluses and temp-basal delivery "
           "are not counted (SMB corrections only).*"]
    (OUT / "dose_response_db.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
