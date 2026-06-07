#!/usr/bin/env python3
"""Phase 10: multi-patient prediction backtest to set a DEFAULT glucose exponent k.

Per-site fitting (Phase 9c) left k weakly identified because each site's MAE-vs-k curve is
flat. Here we instead share ONE global k across patients and fit each patient's level (the
scale α) individually — the right structure for "universal exponent + per-patient level".
Pooling many patients' prediction errors at a shared k aggregates the weak per-site signal,
and a leave-one-patient-out (LOPO) loop gives an honest, out-of-sample default k.

Model:   ISF_i(BG) = α_i · (target/BG)^k     (α_i per patient, k global)
Score:   prediction-error at end-of-insulin-action — scale the loop's predicted drop by
         ISF_cand/ISF_loop, compare to actual_bg_end, MAE.

Patients: the 12-site multisite 4h cache + the single large boost 4h cache (13 total),
all at the end-of-insulin-action horizon. Parallel over patients.

Output: results/phase10_multipatient_k.{json,md}, charts/inv008/fig_multipatient_k.png
Run: python -m inv008.phase10_multipatient_k
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DRIVE = Path("/Users/timstreet/Library/CloudStorage/GoogleDrive-tim.street@liveintheirshoes.com/My Drive/Dynamic ISF data")
MULTI = Path(os.environ.get("MULTISITE_4H_CACHE", DRIVE / "multisite_4h_sample_cache.pkl"))
BOOST = Path(os.environ.get("BOOST_4H_CACHE", DRIVE / "boost_4h_cache.pkl"))
OUT = Path(os.environ.get("DYNISF_ROOT", Path.cwd())) / "results"
CHART = Path(os.environ.get("DYNISF_ROOT", Path.cwd())) / "charts" / "inv008"
TARGET, D = 99.0, 75.0
KGRID = np.round(np.arange(0.5, 5.01, 0.1), 2)
N_WORKERS = min(12, mp.cpu_count())


def load_patients():
    pts = []
    for s in pickle.load(open(MULTI, "rb")):
        pts.append(dict(name=s["name"], tdd=float(s["tdd_median"]),
                        bg=np.asarray(s["bg"], float), isf_loop=np.asarray(s["isf_actual"], float),
                        pred_drop=np.asarray(s["pred_drop"], float),
                        actual_end=np.asarray(s["actual_bg_end"], float)))
    b = pickle.load(open(BOOST, "rb"))["strict"]
    pts.append(dict(name="boost(N=1)", tdd=float(b["tdd_7day"].median()),
                    bg=b["bg"].to_numpy(float), isf_loop=b["variable_sens"].to_numpy(float),
                    pred_drop=(b["bg"] - b["pred_iob_final"]).to_numpy(float),
                    actual_end=b["actual_bg_end"].to_numpy(float)))
    return pts


def patient_curve(p):
    bg, isf_loop, pd_, ae = p["bg"], p["isf_loop"], p["pred_drop"], p["actual_end"]
    m = (np.isfinite(bg) & np.isfinite(isf_loop) & np.isfinite(pd_) & np.isfinite(ae)
         & (isf_loop > 0) & (bg > 0) & (np.abs(pd_) >= 3))
    bg, isf_loop, pd_, ae = bg[m], isf_loop[m], pd_[m], ae[m]
    if len(bg) < 30:
        return None

    def mae(isf):
        return float(np.abs(ae - (bg - pd_ * (isf / isf_loop))).mean())

    def best_scale(base):
        a0 = np.median(isf_loop) / np.median(base)
        return min(mae(a * base) for a in a0 * np.linspace(0.3, 2.5, 121))

    mae_pl = [best_scale((TARGET / bg) ** k) for k in KGRID]
    mae_log = best_scale(math.log(TARGET / D + 1.0) / np.log(bg / D + 1.0))
    return {"name": p["name"], "tdd": p["tdd"], "n": int(len(bg)),
            "mae_pl": mae_pl, "mae_log": round(mae_log, 2),
            "mae_loop": round(mae(isf_loop), 2),
            "best_k": float(KGRID[int(np.argmin(mae_pl))])}


def main():
    pts = load_patients()
    print(f"Phase 10: {len(pts)} patients on {N_WORKERS} workers")
    with mp.Pool(N_WORKERS) as pool:
        res = [r for r in pool.map(patient_curve, pts) if r]
    M = np.array([r["mae_pl"] for r in res])              # patients × k
    names = [r["name"] for r in res]
    ns = np.array([r["n"] for r in res])

    mean_curve = M.mean(axis=0)                            # equal patient weight
    pooled_curve = np.average(M, axis=0, weights=ns)       # window-weighted
    k_eq = float(KGRID[int(np.argmin(mean_curve))])
    k_pool = float(KGRID[int(np.argmin(pooled_curve))])

    # LOPO: choose k on the other patients (equal weight), evaluate held-out at that k
    lopo_k, lopo_mae = [], []
    for i in range(len(res)):
        others = np.delete(M, i, axis=0).mean(axis=0)
        ki = int(np.argmin(others))
        lopo_k.append(float(KGRID[ki]))
        lopo_mae.append(float(M[i, ki]))

    # flatness: mean per-patient MAE at representative k
    def mean_at(k):
        return round(float(M[:, int(np.argmin(np.abs(KGRID - k)))].mean()), 2)
    flat = {str(k): mean_at(k) for k in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0)}

    mae_pl_best = round(float(mean_curve.min()), 2)
    mae_log_mean = round(float(np.mean([r["mae_log"] for r in res])), 2)
    mae_loop_mean = round(float(np.mean([r["mae_loop"] for r in res])), 2)

    summary = {
        "n_patients": len(res), "total_windows": int(ns.sum()),
        "k_default_equal_weight": k_eq, "k_default_pooled": k_pool,
        "lopo_k_median": round(float(np.median(lopo_k)), 2),
        "lopo_k_range": [round(min(lopo_k), 2), round(max(lopo_k), 2)],
        "lopo_mae_median": round(float(np.median(lopo_mae)), 2),
        "mean_mae_at_best_k": mae_pl_best, "mean_mae_log": mae_log_mean,
        "mean_mae_loop": mae_loop_mean,
        "flatness_mean_mae_by_k": flat,
        "per_patient_best_k": {r["name"]: r["best_k"] for r in res},
    }
    OUT.mkdir(exist_ok=True); CHART.mkdir(parents=True, exist_ok=True)
    (OUT / "phase10_multipatient_k.json").write_text(json.dumps(summary, indent=1))

    # figure: per-patient MAE-vs-k (faint) + mean curve + the default
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for r in res:
        ax.plot(KGRID, np.array(r["mae_pl"]) - min(r["mae_pl"]), color="#bbb", lw=0.8, alpha=0.6)
    ax.plot(KGRID, mean_curve - mean_curve.min(), "b-", lw=2.5,
            label=f"mean across patients (min at k={k_eq})")
    ax.axvline(k_eq, color="b", ls="--", lw=1)
    ax.set_xlabel("glucose exponent k"); ax.set_ylabel("MAE − per-curve min (mg/dL)")
    ax.set_title(f"Multi-patient prediction backtest: default k\n"
                 f"{len(res)} patients, end-of-insulin-action; LOPO k median {summary['lopo_k_median']}")
    ax.set_ylim(0, 8); ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(CHART / "fig_multipatient_k.png", dpi=150); plt.close(fig)

    md = ["# Phase 10 — multi-patient prediction backtest for the default glucose exponent k\n",
          f"{len(res)} patients, {int(ns.sum()):,} windows, end-of-insulin-action horizon. "
          "Global k shared, per-patient level α fit individually; prediction-error scored.\n",
          f"- **default k (equal patient weight): {k_eq}**; window-weighted: {k_pool}",
          f"- **leave-one-patient-out k: median {summary['lopo_k_median']}** "
          f"[range {summary['lopo_k_range'][0]}–{summary['lopo_k_range'][1]}]",
          f"- mean per-patient MAE: power-law@best-k **{mae_pl_best}** vs log {mae_log_mean} "
          f"vs loop {mae_loop_mean}",
          "\n## How much does k matter? (mean per-patient MAE by k)\n",
          "| k | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 |", "|---|---|---|---|---|---|---|",
          "| MAE | " + " | ".join(str(flat[str(k)]) for k in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0)) + " |",
          "\n![Multi-patient MAE vs k](charts/inv008/fig_multipatient_k.png)\n",
          "## Per-patient best k\n", "| patient | TDD | n | best k |", "|---|---|---|---|"]
    for r in sorted(res, key=lambda x: x["tdd"]):
        md.append(f"| {r['name'][:18]} | {r['tdd']:.0f} | {r['n']} | {r['best_k']:.1f} |")
    rng = max(flat.values()) - min(flat.values())
    md.append(f"\n**Reading:** the population objective is minimised at k≈{k_eq} (LOPO median "
              f"{summary['lopo_k_median']}); but MAE varies only ~{rng:.1f} mg/dL across k 1–4, "
              "so the exponent is a weak lever — any moderate k in ~1.5–3 is near-optimal. "
              "Power-law beats both log and the loop on the mean. Use the LOPO k as the default; "
              "treat the curve shape as more important than the precise exponent.")
    (OUT / "phase10_multipatient_k.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
