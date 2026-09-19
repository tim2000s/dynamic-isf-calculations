#!/usr/bin/env python3
"""Run the peak-70, divisor-55 Lyumjev model-mismatch analysis.

The principal comparison uses a six-hour oref exponential action model. Four
arms separate the 45 to 70 minute peak change from the 75 to 55 divisor change.
Independent grid cells are evaluated with seven spawned workers.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing as mp
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv011.model import (
    SCENARIOS,
    Scenario,
    activity_exponential,
    dynamic_isf,
    first_crossing_minutes,
    remaining_exponential,
)


DIA_MIN = 360.0
REFERENCE = SCENARIOS[0]
BG_POINTS = (80.0, 99.0, 120.0, 150.0, 180.0, 210.0)
TDD_POINTS = (20.0, 40.0, 60.0, 80.0)
REPORT_TIMES = (15.0, 30.0, 45.0, 60.0, 70.0, 90.0, 120.0, 180.0, 240.0, 300.0)


def _cell(job: tuple[str, float, float, Scenario, float]) -> dict:
    equation, bg, tdd, scenario, dia_min = job
    time = np.arange(0.0, dia_min + 0.001, 1.0)
    ref_remaining = remaining_exponential(time, dia_min, REFERENCE.peak_min)
    arm_remaining = remaining_exponential(time, dia_min, scenario.peak_min)
    ref_activity = activity_exponential(time, dia_min, REFERENCE.peak_min)
    arm_activity = activity_exponential(time, dia_min, scenario.peak_min)
    ref_isf = float(dynamic_isf(equation, bg, tdd, REFERENCE.divisor))
    arm_isf = float(dynamic_isf(equation, bg, tdd, scenario.divisor))

    with np.errstate(divide="ignore", invalid="ignore"):
        remaining_ratio = arm_remaining * arm_isf / (ref_remaining * ref_isf)
        delivered_ratio = (1.0 - arm_remaining) * arm_isf / (
            (1.0 - ref_remaining) * ref_isf
        )
        instantaneous_ratio = arm_activity * arm_isf / (ref_activity * ref_isf)
    remaining_ratio[(ref_remaining <= 1e-12) | ~np.isfinite(remaining_ratio)] = np.nan
    delivered_ratio[(1.0 - ref_remaining <= 1e-12) | ~np.isfinite(delivered_ratio)] = np.nan
    instantaneous_ratio[(ref_activity <= 1e-12) | ~np.isfinite(instantaneous_ratio)] = np.nan

    at_times = {}
    for minute in REPORT_TIMES:
        i = int(minute)
        at_times[str(int(minute))] = {
            "reference_fraction_remaining": float(ref_remaining[i]),
            "scenario_fraction_remaining": float(arm_remaining[i]),
            "remaining_effect_ratio": float(remaining_ratio[i]),
            "delivered_effect_ratio": float(delivered_ratio[i]),
            "instantaneous_effect_ratio": float(instantaneous_ratio[i]),
        }

    return {
        "equation": equation,
        "bg_mgdl": bg,
        "tdd_u_day": tdd,
        "scenario": scenario.name,
        "peak_min": scenario.peak_min,
        "divisor": scenario.divisor,
        "reference_isf_mgdl_u": ref_isf,
        "scenario_isf_mgdl_u": arm_isf,
        "total_effect_ratio": arm_isf / ref_isf,
        "correction_requirement_ratio": ref_isf / arm_isf,
        "remaining_effect_crossing_min": first_crossing_minutes(time, remaining_ratio),
        "delivered_effect_crossing_min": first_crossing_minutes(time, delivered_ratio),
        "at_times": at_times,
    }


def principal_run(workers: int = 7) -> list[dict]:
    jobs = [
        (equation, bg, tdd, scenario, DIA_MIN)
        for equation in ("v1", "v2")
        for bg in BG_POINTS
        for tdd in TDD_POINTS
        for scenario in SCENARIOS
    ]
    context = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        return list(pool.map(_cell, jobs, chunksize=4))


def _sensitivity_cell(job: tuple[str, float, float, float, float, float]) -> dict:
    equation, bg, dia, reference_peak, assumed_peak, minute = job
    ref_isf = float(dynamic_isf(equation, bg, 40.0, 75.0))
    wrong_isf = float(dynamic_isf(equation, bg, 40.0, 55.0))
    ref_r = float(remaining_exponential(minute, dia, reference_peak))
    wrong_r = float(remaining_exponential(minute, dia, assumed_peak))
    return {
        "equation": equation,
        "bg_mgdl": bg,
        "dia_min": dia,
        "reference_peak_min": reference_peak,
        "assumed_peak_min": assumed_peak,
        "time_min": minute,
        "remaining_effect_ratio": wrong_r * wrong_isf / (ref_r * ref_isf),
        "correction_requirement_ratio": ref_isf / wrong_isf,
    }


def sensitivity_run(workers: int = 7) -> list[dict]:
    jobs = [
        (equation, bg, dia, reference_peak, assumed_peak, minute)
        for equation in ("v1", "v2")
        for bg in BG_POINTS
        for dia in (300.0, 360.0, 420.0)
        for reference_peak in (35.0, 40.0, 45.0, 50.0)
        for assumed_peak in (65.0, 70.0, 75.0)
        for minute in (30.0, 60.0, 90.0, 120.0, 180.0, 240.0)
    ]
    context = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        return list(pool.map(_sensitivity_cell, jobs, chunksize=24))


def summarise_sensitivity(rows: list[dict]) -> list[dict]:
    d = pd.DataFrame(rows)
    out = []
    for (equation, bg, minute), g in d.groupby(["equation", "bg_mgdl", "time_min"]):
        x = g.remaining_effect_ratio.to_numpy(float)
        out.append({
            "equation": equation,
            "bg_mgdl": float(bg),
            "time_min": float(minute),
            "n_parameter_combinations": int(len(x)),
            "remaining_effect_ratio_p05": float(np.quantile(x, 0.05)),
            "remaining_effect_ratio_median": float(np.median(x)),
            "remaining_effect_ratio_p95": float(np.quantile(x, 0.95)),
        })
    return out


def temp_basal_example() -> dict:
    """Two hours at 1 U/h above scheduled basal, on a five-minute grid."""
    step = 5.0
    time = np.arange(0.0, 480.0 + step, step)
    delivery = np.zeros_like(time)
    delivery[(time >= 0.0) & (time < 120.0)] = step / 60.0
    ref_kernel = remaining_exponential(np.arange(0.0, DIA_MIN + step, step),
                                       DIA_MIN, REFERENCE.peak_min)
    wrong_kernel = remaining_exponential(np.arange(0.0, DIA_MIN + step, step),
                                         DIA_MIN, 70.0)
    ref_iob = np.convolve(delivery, ref_kernel, mode="full")[:len(time)]
    wrong_iob = np.convolve(delivery, wrong_kernel, mode="full")[:len(time)]
    out = {
        "description": "1 U/h above scheduled basal for 120 minutes, total 2 U",
        "bg_mgdl": 150.0,
        "tdd_u_day": 40.0,
        "time_min": time.tolist(),
        "delivery_u": delivery.tolist(),
        "reference_iob_u": ref_iob.tolist(),
        "assumed_iob_u": wrong_iob.tolist(),
        "by_equation": {},
    }
    for equation in ("v1", "v2"):
        ref_isf = float(dynamic_isf(equation, 150.0, 40.0, 75.0))
        wrong_isf = float(dynamic_isf(equation, 150.0, 40.0, 55.0))
        ref_effect = ref_iob * ref_isf
        wrong_effect = wrong_iob * wrong_isf
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = wrong_effect / ref_effect
        out["by_equation"][equation] = {
            "reference_effect_to_come_mgdl": ref_effect.tolist(),
            "assumed_effect_to_come_mgdl": wrong_effect.tolist(),
            "ratio": [float(x) if np.isfinite(x) else None for x in ratio],
            "selected_times": {
                str(minute): {
                    "reference_iob_u": float(ref_iob[int(minute / step)]),
                    "assumed_iob_u": float(wrong_iob[int(minute / step)]),
                    "remaining_effect_ratio": float(ratio[int(minute / step)]),
                }
                for minute in (60, 120, 180, 240, 300, 360)
            },
        }
    return out


def make_figures(rows: list[dict], temp_basal: dict, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    time = np.arange(0.0, DIA_MIN + 0.001, 1.0)
    r45 = remaining_exponential(time, DIA_MIN, 45.0)
    r70 = remaining_exponential(time, DIA_MIN, 70.0)
    a45 = activity_exponential(time, DIA_MIN, 45.0)
    a70 = activity_exponential(time, DIA_MIN, 70.0)

    fig, ax = plt.subplots(1, 2, figsize=(11.8, 4.6))
    ax[0].plot(time, a45 * 60.0, lw=2.3, label="45-minute peak")
    ax[0].plot(time, a70 * 60.0, lw=2.3, label="70-minute peak")
    ax[0].set(xlabel="minutes after delivery", ylabel="fraction acting per hour",
              title="The assumed action peak shifts activity later")
    ax[1].plot(time, r45, lw=2.3, label="45-minute peak")
    ax[1].plot(time, r70, lw=2.3, label="70-minute peak")
    ax[1].set(xlabel="minutes after delivery", ylabel="fraction recorded as IOB",
              title="The later curve retains more modelled IOB")
    for a in ax:
        a.grid(alpha=0.25)
        a.legend()
    fig.tight_layout()
    fig.savefig(outdir / "fig_inv011_action_timing.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    d = pd.DataFrame(rows)
    combined = d[(d.scenario == "combined_70_55") & (d.tdd_u_day == 40.0)]
    fig, ax = plt.subplots(1, 2, figsize=(11.8, 4.7))
    for equation, colour in (("v1", "#e66101"), ("v2", "#1b9e77")):
        g = combined[combined.equation == equation].sort_values("bg_mgdl")
        ax[0].plot(g.bg_mgdl, g.correction_requirement_ratio, marker="o", lw=2.3,
                   color=colour, label=equation.upper())
    ax[0].axhline(1.0, color="#555", lw=1)
    ax[0].set(xlabel="glucose (mg/dL)", ylabel="inverse-ISF multiplier",
              title="Divisor 55 assigns more insulin per mg/dL of modelled effect")
    ax[0].grid(alpha=0.25)
    ax[0].legend()

    for equation, ls in (("v1", "-"), ("v2", "--")):
        row = next(r for r in rows if r["scenario"] == "combined_70_55"
                   and r["equation"] == equation and r["bg_mgdl"] == 99.0
                   and r["tdd_u_day"] == 40.0)
        vals = [row["at_times"][str(int(t))]["remaining_effect_ratio"] for t in REPORT_TIMES]
        ax[1].plot(REPORT_TIMES, vals, marker="o", lw=2.3, ls=ls,
                   label=f"{equation.upper()}, glucose 99")
    ax[1].axhline(1.0, color="#555", lw=1)
    ax[1].set(xlabel="minutes after delivery",
              ylabel="modelled effect still to come, wrong / reference",
              title="The direction of the timing error reverses")
    ax[1].grid(alpha=0.25)
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(outdir / "fig_inv011_combined_mismatch.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11.8, 4.7), sharey=True)
    scenarios = [s.name for s in SCENARIOS]
    labels = ["45 / 75\nmatched", "45 / 55\ndivisor", "70 / 75\npeak", "70 / 55\nboth"]
    for j, equation in enumerate(("v1", "v2")):
        vals = []
        for scenario in scenarios:
            row = next(r for r in rows if r["scenario"] == scenario
                       and r["equation"] == equation and r["bg_mgdl"] == 99.0
                       and r["tdd_u_day"] == 40.0)
            vals.append(row["at_times"]["60"]["remaining_effect_ratio"])
        ax[j].bar(np.arange(4), vals, color=("#4c78a8", "#f58518", "#72b7b2", "#e45756"))
        ax[j].axhline(1.0, color="#333", ls="--")
        ax[j].set_xticks(np.arange(4), labels)
        ax[j].set_title(f"{equation.upper()} at glucose 99 mg/dL")
        ax[j].set_xlabel("peak (minutes) / divisor")
        ax[j].grid(alpha=0.2, axis="y")
    ax[0].set_ylabel("modelled effect still to come at 60 minutes\nrelative to 45 / 75")
    fig.tight_layout()
    fig.savefig(outdir / "fig_inv011_four_arms.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    t = np.asarray(temp_basal["time_min"], dtype=float)
    fig, ax = plt.subplots(1, 2, figsize=(11.8, 4.7))
    ax[0].plot(t, temp_basal["reference_iob_u"], lw=2.3, label="45-minute peak")
    ax[0].plot(t, temp_basal["assumed_iob_u"], lw=2.3, label="70-minute peak")
    ax[0].axvspan(0, 120, color="#999", alpha=0.12, label="+1 U/h temp basal")
    ax[0].set(xlabel="minutes from start", ylabel="IOB attributed to the extra basal (U)",
              title="A temporary-basal increase compounds the timing difference")
    ax[0].grid(alpha=0.25)
    ax[0].legend()
    for equation, colour in (("v1", "#e66101"), ("v2", "#1b9e77")):
        ratio = np.asarray(temp_basal["by_equation"][equation]["ratio"], dtype=float)
        ax[1].plot(t, ratio, lw=2.3, color=colour, label=equation.upper())
    ax[1].axhline(1.0, color="#333", ls="--")
    ax[1].set(xlabel="minutes from start",
              ylabel="modelled effect still to come, wrong / reference",
              title="The same reversal appears with temporary basal")
    ax[1].set_ylim(0, 5)
    ax[1].grid(alpha=0.25)
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(outdir / "fig_inv011_temp_basal.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def result_markdown(result: dict) -> str:
    rows = result["principal"]
    lines = [
        "# INV-011: Peak and divisor mismatch",
        "",
        "This is a software-model sensitivity analysis. It does not estimate clinical outcomes.",
        "",
        "## Inverse-ISF effect of divisor 55",
        "",
        "| Glucose | V1, 55 versus 75 | V2, 55 versus 75 |",
        "|---:|---:|---:|",
    ]
    for bg in BG_POINTS:
        vals = {}
        for equation in ("v1", "v2"):
            row = next(r for r in rows if r["scenario"] == "combined_70_55"
                       and r["equation"] == equation and r["bg_mgdl"] == bg
                       and r["tdd_u_day"] == 40.0)
            vals[equation] = row["correction_requirement_ratio"]
        lines.append(f"| {bg:.0f} mg/dL | {vals['v1']:.2f} times | {vals['v2']:.2f} times |")

    lines += [
        "",
        "The multiplier is the insulin assigned per mg/dL of modelled effect. It does not mean "
        "that a positive correction is delivered at or below target.",
        "",
        "## Effect still attributed to one unit at glucose 99 mg/dL",
        "",
        "| Time | Reference 45-minute IOB | Assumed 70-minute IOB | V1 combined ratio | V2 combined ratio |",
        "|---:|---:|---:|---:|---:|",
    ]
    v1 = next(r for r in rows if r["scenario"] == "combined_70_55" and r["equation"] == "v1"
              and r["bg_mgdl"] == 99.0 and r["tdd_u_day"] == 40.0)
    v2 = next(r for r in rows if r["scenario"] == "combined_70_55" and r["equation"] == "v2"
              and r["bg_mgdl"] == 99.0 and r["tdd_u_day"] == 40.0)
    for minute in REPORT_TIMES:
        key = str(int(minute))
        x = v1["at_times"][key]
        lines.append(
            f"| {minute:.0f} min | {x['reference_fraction_remaining']:.3f} | "
            f"{x['scenario_fraction_remaining']:.3f} | {x['remaining_effect_ratio']:.2f} | "
            f"{v2['at_times'][key]['remaining_effect_ratio']:.2f} |"
        )
    lines += [
        "",
        f"V1 crosses from underestimating to overestimating effect still to come at "
        f"{v1['remaining_effect_crossing_min']:.0f} minutes. V2 crosses at "
        f"{v2['remaining_effect_crossing_min']:.0f} minutes in this reference case.",
        "",
        "The mismatch increases calculated correction at delivery, understates early action, "
        "and retains too much modelled effect later. A full controller can moderate or amplify "
        "those errors through IOB limits and low-glucose prediction. This analysis does not "
        "predict the net delivered insulin or clinical outcome of a live controller.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=7)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    principal = principal_run(args.workers)
    sensitivity_raw = sensitivity_run(args.workers)
    temp_basal = temp_basal_example()
    result = {
        "analysis": "INV-011 peak and divisor mismatch",
        "reference": {"peak_min": 45, "divisor": 75, "dia_min": DIA_MIN},
        "implementation": {"peak_min": 70, "divisor": 55, "dia_min": DIA_MIN},
        "interpretation_limit": (
            "Software-model sensitivity only. The 45-minute reference is the AAPS Lyumjev "
            "model parameter and is not an individual pharmacodynamic measurement."
        ),
        "principal": principal,
        "sensitivity_summary": summarise_sensitivity(sensitivity_raw),
        "temp_basal_example": temp_basal,
        "sensitivity_grid": {
            "dia_min": [300, 360, 420],
            "reference_peak_min": [35, 40, 45, 50],
            "assumed_peak_min": [65, 70, 75],
            "glucose_mgdl": list(BG_POINTS),
            "time_min": [30, 60, 90, 120, 180, 240],
        },
    }
    results_dir = args.root / "results"
    charts_dir = args.root / "charts" / "inv011"
    results_dir.mkdir(parents=True, exist_ok=True)
    charts_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "inv011_peak_divisor_mismatch.json").write_text(json.dumps(result, indent=2))
    md = result_markdown(result)
    (results_dir / "inv011_peak_divisor_mismatch.md").write_text(md)
    make_figures(principal, temp_basal, charts_dir)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
