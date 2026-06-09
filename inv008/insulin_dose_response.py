#!/usr/bin/env python3
"""Insulin dose-response vs glucose, with mean-reversion removed (physiological, not statistical).

Glucose in a fasting T1D window falls because INSULIN drives it down (+ insulin-independent renal/
mass-action clearance at high BG). So the realised-ISF-vs-BG curve is contaminated when a high BG was
*already falling on the way in* — that drop is partly the prior trajectory resolving, not the marginal
insulin dose-response at that glucose. To see whether insulin is genuinely less effective at high
glucose (glucotoxic resistance → the Diabeloop power law) we must isolate windows where the drop had
to be INITIATED by insulin: glucose flat or rising at entry, not already mean-reverting.

ratio = actual_drop / predicted_drop_static = realised_ISF / profile_ISF   (<1 ⇒ insulin did less
than the static profile expected). Stratify by entry trajectory (pre_slope, mg/dL per 5 min):
    falling  (already coming down — mean-reverting, contaminated)
    flat     (stable set-point — insulin must drive the drop: the CLEAN dose-response)
    rising   (being pushed up — insulin must overcome it: also clean for insulin attribution)

Hypothesis under test (Tim / Diabeloop): in the CLEAN (non-falling) strata, realised ISF should
FALL with glucose (resistance) even though the pooled curve rises — i.e. our rising curve was a
mean-reversion artifact and the power law is right for genuine hyperglycaemia.
Renal/mass-action clearance (BG≳180) is insulin-independent and flagged, not removed here.
Parquet-only. Output: results/insulin_dose_response.{json,md}, charts/inv008/fig_insulin_dose_response.png
Run: python -m inv008.insulin_dose_response
"""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config, err_common as ec

OUT = config.ROOT / "results"
CHART = config.ROOT / "charts" / "inv008"
BANDS = [(100, 120), (120, 145), (145, 175), (175, 205), (205, 260)]
LBL = [f"{a}-{b}" for a, b in BANDS]
CTR = [(a + b) / 2 for a, b in BANDS]
MIN_CELL = 30          # min windows for a (stratum, band) median


def main():
    d = ec.load_windows().copy()
    d["actual_drop"] = d.bg - d.bg_end
    d["pred_drop"] = d.err_static + d.actual_drop
    d = d[(d.pred_drop > 0) & d.pre_slope.notna()].reset_index(drop=True)
    d["ratio"] = (d.actual_drop / d.pred_drop).clip(-1.0, 3.0)   # realised ÷ profile ISF

    # entry strata (pre_slope is mg/dL per 5 min)
    d["entry"] = np.where(d.pre_slope < -1.0, "falling",
                  np.where(d.pre_slope > 1.0, "rising", "flat"))
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)

    def curve(sub):
        out = []
        for a, b in BANDS:
            m = (sub.bg >= a) & (sub.bg < b)
            out.append({"band": f"{a}-{b}", "n": int(m.sum()),
                        "ratio": round(float(sub.ratio[m].median()), 2) if m.sum() >= MIN_CELL else None,
                        "iob": round(float(sub.iob[m].median()), 1) if m.sum() >= MIN_CELL else None})
        return out

    strata = {s: curve(d[d.entry == s]) for s in ["falling", "flat", "rising"]}
    strata["all"] = curve(d)
    clean = curve(d[d.entry != "falling"])   # flat+rising = insulin-must-drive-it
    strata["clean_nonfalling"] = clean

    # implied power-law k from the clean curve, fit on the actionable range (>=120)
    def fit_k(c):
        pts = [(CTR[i], c[i]["ratio"]) for i in range(len(c)) if c[i]["ratio"] and CTR[i] >= 120]
        if len(pts) < 3:
            return None
        bg = np.array([p[0] for p in pts]); r = np.array([p[1] for p in pts])
        r = np.clip(r, 0.05, None)
        # ratio ≈ (BG/100)^(-k_eff)·scale  →  log r = log scale − k·log(BG/100)
        A = np.column_stack([np.ones(len(bg)), np.log(bg / 100.0)])
        beta, *_ = np.linalg.lstsq(A, np.log(r), rcond=None)
        return round(-float(beta[1]), 2)   # k>0 ⇒ resistance (ratio falls with BG) = power-law direction

    k_clean = fit_k(clean)
    k_falling = fit_k(strata["falling"])
    k_all = fit_k(strata["all"])

    # high-BG read: ratio at 175-260 vs 120-145, per stratum
    def hi_vs_mid(c):
        d_ = {x["band"]: x["ratio"] for x in c}
        hi = [d_[b] for b in ["175-205", "205-260"] if d_.get(b)]
        mid = [d_[b] for b in ["120-145", "145-175"] if d_.get(b)]
        return (round(np.mean(hi) - np.mean(mid), 2) if hi and mid else None)

    summary = {
        "n_windows": int(len(d)),
        "entry_counts": {s: int((d.entry == s).sum()) for s in ["falling", "flat", "rising"]},
        "ratio_by_bg_and_entry": strata,
        "implied_k": {"clean_nonfalling": k_clean, "falling": k_falling, "all": k_all,
                      "note": "k>0 = realised ISF FALLS with glucose = resistance = power-law direction; "
                              "k≤0 = flat/rising = insulin still effective at high BG."},
        "hi_minus_mid_ratio": {s: hi_vs_mid(strata[s]) for s in ["falling", "flat", "rising", "clean_nonfalling"]},
        "verdict": ("POWER-LAW vindicated in clean windows" if (k_clean and k_clean > 0.3)
                    else "still flat/rising once mean-reversion removed — power law NOT recovered"),
        "caveat": "BG≳180 carries insulin-INDEPENDENT renal/mass-action clearance that inflates the "
                  "ratio at the top band regardless of entry trajectory; a true insulin-only resistance "
                  "test would also subtract that (needs dose-level data).",
    }
    (OUT / "insulin_dose_response.json").write_text(json.dumps(summary, indent=1))

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    cols = {"all": "#999", "falling": "#1f77b4", "flat": "#2ca02c", "rising": "#d62728"}
    for s in ["all", "falling", "flat", "rising"]:
        ys = [x["ratio"] for x in strata[s]]
        ax[0].plot(CTR, ys, "o-", color=cols[s], lw=2, label=f"{s} (n={summary['entry_counts'].get(s, len(d))})")
    ax[0].axhline(1.0, color="k", ls="--", lw=1)
    ax[0].set_xlabel("glucose (mg/dL)"); ax[0].set_ylabel("realised ÷ profile ISF")
    ax[0].set_title("Insulin dose-response vs glucose, by entry trajectory\n(falling = mean-reverting; flat/rising = insulin-driven)")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    # clean vs falling, with the power-law expectation overlaid
    ax[1].plot(CTR, [x["ratio"] for x in clean], "o-", color="#2ca02c", lw=2.5, label=f"clean (k={k_clean})")
    ax[1].plot(CTR, [x["ratio"] for x in strata["falling"]], "s--", color="#1f77b4", lw=2, label=f"falling (k={k_falling})")
    pl = [(120.0 / b) ** 1.0 for b in CTR]   # a gentle power-law reference (k=1), normalised at 120
    pl = [p / pl[1] * clean[1]["ratio"] for p in pl] if clean[1]["ratio"] else pl
    ax[1].plot(CTR, pl, ":", color="#d62728", lw=2, label="power-law k=1 (resistance ref)")
    ax[1].axhline(1.0, color="k", ls="--", lw=1)
    ax[1].set_xlabel("glucose (mg/dL)"); ax[1].set_ylabel("realised ÷ profile ISF")
    ax[1].set_title("Mean-reversion removed: does resistance appear?"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_insulin_dose_response.png", dpi=150); plt.close(fig)

    # ---- markdown ----
    md = ["# Insulin dose-response vs glucose, mean-reversion removed\n",
          f"{len(d):,} windows. ratio = realised ÷ profile ISF (<1 ⇒ insulin did less than profile "
          "expected). Stratified by entry trajectory: **falling** = already mean-reverting; "
          "**flat/rising** = insulin must drive the drop (clean dose-response).\n",
          f"Entry counts: falling {summary['entry_counts']['falling']:,}, flat "
          f"{summary['entry_counts']['flat']:,}, rising {summary['entry_counts']['rising']:,}.\n",
          "## realised ÷ profile ISF by glucose and entry\n",
          "| BG band | all | falling | flat | rising | clean(flat+rising) |", "|---|---|---|---|---|---|"]
    for i, (a, b) in enumerate(BANDS):
        row = [f"{a}-{b}"]
        for s in ["all", "falling", "flat", "rising", "clean_nonfalling"]:
            v = strata[s][i]["ratio"]; row.append(str(v) if v is not None else "–")
        md.append("| " + " | ".join(row) + " |")
    md += [f"\n**Implied k** (k>0 = resistance / power-law direction): clean **{k_clean}**, "
           f"falling {k_falling}, all {k_all}.\n",
           f"High−mid ratio shift: clean {summary['hi_minus_mid_ratio']['clean_nonfalling']}, "
           f"falling {summary['hi_minus_mid_ratio']['falling']}.\n",
           f"**Verdict: {summary['verdict']}.**\n",
           "![insulin dose-response](charts/inv008/fig_insulin_dose_response.png)\n",
           "*" + summary["caveat"] + "*"]
    (OUT / "insulin_dose_response.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
