#!/usr/bin/env python3
"""Sensitivity of the Dynamic ISF comparisons to insulinDivisor.

AndroidAPS maps the selected insulin peak to an insulin divisor:

    peak > 65 min  -> divisor 55   (rapid acting, nominal peak 75)
    peak > 50 min  -> divisor 65   (ultra rapid, nominal peak 55)
    otherwise      -> divisor 75   (Lyumjev-style, nominal peak 45)

The principal INV-008 replay fixed divisor 75 because insulin type was not
available consistently. This module recalculates both equations at divisors 55,
65 and 75 without changing glucose, blended TDD, observation windows or the
loop-derived activity integral. It therefore isolates divisor sensitivity.

Outputs:
    results/divisor_sensitivity.json
    results/divisor_sensitivity.md
    charts/inv008/fig_divisor_sensitivity.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config
from inv008.dynisf import isf_v1, isf_v2


DIVISORS = (55.0, 65.0, 75.0)
BG_POINTS = (80.0, 99.0, 120.0, 150.0, 180.0, 210.0)
AAPS_MAPPING = {
    "rapid_acting_peak_75": 55,
    "ultra_rapid_peak_55": 65,
    "lyumjev_style_peak_45": 75,
}


def equation_effects() -> list[dict]:
    rows = []
    for bg in BG_POINTS:
        ref1 = float(isf_v1(bg, 40.0, divisor=75.0))
        ref2 = float(isf_v2(bg, 40.0, divisor=75.0))
        row: dict = {"bg": bg}
        for divisor in DIVISORS:
            v1 = float(isf_v1(bg, 40.0, divisor=divisor))
            v2 = float(isf_v2(bg, 40.0, divisor=divisor))
            key = str(int(divisor))
            row[f"v1_isf_ratio_d{key}_vs_d75"] = v1 / ref1
            row[f"v1_correction_ratio_d{key}_vs_d75"] = ref1 / v1
            row[f"v2_isf_ratio_d{key}_vs_d75"] = v2 / ref2
            row[f"v2_correction_ratio_d{key}_vs_d75"] = ref2 / v2
        rows.append(row)
    return rows


def same_window(path: Path) -> dict:
    d = pd.read_parquet(path).reset_index(drop=True)
    actual_drop = d.bg.to_numpy(float) - d.bg_end.to_numpy(float)
    # Existing V1 error = predicted drop - actual drop. This recovers the
    # common loop-derived activity integral without rerunning window selection.
    activity = (d.err_v1.to_numpy(float) + actual_drop) / d.isf_v1.to_numpy(float)

    errors: dict[str, np.ndarray] = {
        "static": d.err_static.to_numpy(float),
        "loop": d.err_loop.to_numpy(float),
    }
    summary: dict = {
        "n_windows": int(len(d)),
        "n_people": int(d.user.nunique()),
        "static": {
            "mae": float(np.nanmedian(np.abs(errors["static"]))),
            "bias": float(np.nanmedian(errors["static"])),
        },
        "loop": {
            "mae": float(np.nanmedian(np.abs(errors["loop"]))),
            "bias": float(np.nanmedian(errors["loop"])),
        },
        "by_divisor": {},
    }
    bg = d.bg.to_numpy(float)
    tdd = d.tdd.to_numpy(float)
    for divisor in DIVISORS:
        key = str(int(divisor))
        v1 = isf_v1(bg, tdd, divisor=divisor)
        v2 = isf_v2(bg, tdd, divisor=divisor)
        e1 = activity * v1 - actual_drop
        e2 = activity * v2 - actual_drop
        errors[f"v1_d{key}"] = e1
        errors[f"v2_d{key}"] = e2
        summary["by_divisor"][key] = {
            "v1": {"mae": float(np.nanmedian(np.abs(e1))),
                   "bias": float(np.nanmedian(e1))},
            "v2": {"mae": float(np.nanmedian(np.abs(e2))),
                   "bias": float(np.nanmedian(e2))},
            "fraction_v2_isf_above_v1": float(np.nanmean(v2 > v1)),
        }

        wins = {name: 0 for name in ("static", "loop", "v1", "v2")}
        scored = 0
        for _, idx in d.groupby("user").groups.items():
            if len(idx) < 60:
                continue
            ii = np.asarray(idx, dtype=int)
            vals = {
                "static": np.nanmedian(np.abs(errors["static"][ii])),
                "loop": np.nanmedian(np.abs(errors["loop"][ii])),
                "v1": np.nanmedian(np.abs(e1[ii])),
                "v2": np.nanmedian(np.abs(e2[ii])),
            }
            wins[min(vals, key=vals.get)] += 1
            scored += 1
        summary["by_divisor"][key]["per_person_best"] = wins
        summary["by_divisor"][key]["n_people_scored"] = scored
    return summary


def pointwise(path: Path) -> dict:
    d = pd.read_parquet(path, columns=["subject_id", "bg", "tdd_blend", "isf_eff"])
    out: dict = {"n_points": int(len(d)), "n_people": int(d.subject_id.nunique()),
                 "by_divisor": {}}
    bg = d.bg.to_numpy(float)
    tdd = d.tdd_blend.to_numpy(float)
    proxy = d.isf_eff.to_numpy(float)
    for divisor in DIVISORS:
        key = str(int(divisor))
        out["by_divisor"][key] = {}
        for name, fn in (("v1", isf_v1), ("v2", isf_v2)):
            calc = fn(bg, tdd, divisor=divisor)
            pp = pd.DataFrame({"subject_id": d.subject_id, "value": calc}) \
                .groupby("subject_id").value.median()
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = calc / proxy
            ratio = ratio[np.isfinite(ratio) & (ratio > 0) & (ratio < 100)]
            out["by_divisor"][key][name] = {
                "per_person_median_calculated_isf": float(pp.median()),
                "median_positive_ratio_to_proxy": float(np.median(ratio)),
                "fraction_within_30_percent": float(np.mean((ratio > 0.7) & (ratio < 1.3))),
            }
    return out


def make_figure(eq: list[dict], sw: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    bg = np.array([r["bg"] for r in eq])
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    colours = {"v1": "#e66101", "v2": "#1b9e77"}
    for form in ("v1", "v2"):
        for divisor, ls in ((55, "-"), (65, "--")):
            vals = [r[f"{form}_correction_ratio_d{divisor}_vs_d75"] for r in eq]
            ax[0].plot(bg, vals, marker="o", ls=ls, lw=2,
                       color=colours[form], label=f"{form.upper()}, divisor {divisor} vs 75")
    ax[0].axhline(1.0, color="#555", lw=1)
    ax[0].set_xlabel("glucose (mg/dL)")
    ax[0].set_ylabel("correction multiplier relative to divisor 75")
    ax[0].set_title("Divisor mismatch changes correction strength")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.25)

    x = np.arange(len(DIVISORS)); width = 0.34
    v1 = [sw["by_divisor"][str(int(d))]["v1"]["mae"] for d in DIVISORS]
    v2 = [sw["by_divisor"][str(int(d))]["v2"]["mae"] for d in DIVISORS]
    ax[1].bar(x - width / 2, v1, width, color=colours["v1"], label="V1")
    ax[1].bar(x + width / 2, v2, width, color=colours["v2"], label="V2")
    ax[1].axhline(sw["static"]["mae"], color="#333", ls="--",
                  label=f"static ({sw['static']['mae']:.1f})")
    ax[1].set_xticks(x, [str(int(d)) for d in DIVISORS])
    ax[1].set_xlabel("insulin divisor")
    ax[1].set_ylabel("median absolute endpoint error (mg/dL)")
    ax[1].set_title("Same windows, each divisor")
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def markdown(res: dict) -> str:
    lines = [
        "# Insulin-divisor sensitivity",
        "",
        "AndroidAPS maps a nominal 75-minute rapid-insulin peak to divisor 55, "
        "a 55-minute ultra-rapid peak to divisor 65, and a 45-minute Lyumjev-style "
        "peak to divisor 75. The principal audit replay fixed divisor 75 because "
        "insulin type was not consistently available.",
        "",
        "## Same-window endpoint prediction",
        "",
        "| divisor | V1 MAE | V1 bias | V2 MAE | V2 bias |",
        "|---:|---:|---:|---:|---:|",
    ]
    for divisor in DIVISORS:
        d = res["same_window"]["by_divisor"][str(int(divisor))]
        lines.append(f"| {int(divisor)} | {d['v1']['mae']:.1f} | {d['v1']['bias']:+.1f} | "
                     f"{d['v2']['mae']:.1f} | {d['v2']['bias']:+.1f} |")
    lines += [
        "",
        f"Static ISF MAE was {res['same_window']['static']['mae']:.1f} mg/dL on the same windows.",
        "",
        "## Six-hour action-balance comparison",
        "",
        "| divisor | V1 / proxy | V1 within 30% | V2 / proxy | V2 within 30% |",
        "|---:|---:|---:|---:|---:|",
    ]
    for divisor in DIVISORS:
        d = res["pointwise"]["by_divisor"][str(int(divisor))]
        lines.append(f"| {int(divisor)} | {d['v1']['median_positive_ratio_to_proxy']:.2f} | "
                     f"{100*d['v1']['fraction_within_30_percent']:.1f}% | "
                     f"{d['v2']['median_positive_ratio_to_proxy']:.2f} | "
                     f"{100*d['v2']['fraction_within_30_percent']:.1f}% |")
    lines += [
        "",
        "Divisor 55 lowers both calculated ISFs and therefore makes correction stronger. "
        "The effect is modest for V1 and large for V2 near its logarithmic floor. V2 improves "
        "when divisor 55 is used, but remains worse than static ISF and V1 in the same-window test.",
        "",
        "![Divisor sensitivity](charts/inv008/fig_divisor_sensitivity.png)",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--windows", type=Path, required=True,
                    help="head_to_head_windows.parquet from INV-008")
    ap.add_argument("--pointwise", type=Path,
                    default=config.ROOT / "results" / "inv009_pointwise.parquet")
    args = ap.parse_args()
    config.ensure_dirs()
    result = {
        "aaps_peak_to_divisor": AAPS_MAPPING,
        "equation_effects": equation_effects(),
        "same_window": same_window(args.windows),
        "pointwise": pointwise(args.pointwise),
        "interpretation": (
            "Divisor mismatch materially changes absolute ISF, especially V2 near its floor. "
            "Across divisors 55, 65 and 75, neither equation beats static ISF on the selected "
            "same windows; V2 remains the least accurate dynamic candidate."
        ),
    }
    out = config.ROOT / "results"
    out.mkdir(exist_ok=True)
    (out / "divisor_sensitivity.json").write_text(json.dumps(result, indent=2))
    (out / "divisor_sensitivity.md").write_text(markdown(result))
    make_figure(result["equation_effects"], result["same_window"],
                config.ROOT / "charts" / "inv008" / "fig_divisor_sensitivity.png")
    print(markdown(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
