"""Figures for the INV-009 write-up.

Static images for a document, so there is no hover layer to build; identity is
carried by direct labels on every series, which is also what the palette's
contrast relief requires. One measure per axis throughout, no second y scale.

    python3 -m inv009.figures
"""
from __future__ import annotations

import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config

# Categorical slots 1-3 of the validated palette, in fixed order. Three because
# these are scatter and dot forms, where every pair is on screen at once.
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, MUTED = "#0b0b0b", "#52514e", "#e3e2dc", "#9a9890"
RED = "#e34948"      # status, reserved: used only to flag a failing candidate

plt.rcParams.update({
    "figure.dpi": 160, "savefig.dpi": 160, "figure.facecolor": "white",
    "axes.facecolor": "white", "font.size": 9,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.titlesize": 10, "axes.titleweight": "semibold", "axes.titlecolor": INK,
    "legend.frameon": False, "legend.fontsize": 8,
})


def _despine(ax, keep=("left", "bottom")):
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def _laws(ax, top_labels=True):
    """The candidate exponents, as reference lines labelled clear of the data."""
    marks = ((0.0, "no effect", ":"), (-0.5, "\u221aTDD", "--"),
             (-1.0, "v1", "-"), (-2.0, "v2", "-"))
    for val, name, style in marks:
        ax.axvline(val, color=INK2 if style == "-" else MUTED, linestyle=style,
                   linewidth=1.0 if style == "-" else 0.8, zorder=1)
    if top_labels:
        for val, name, _ in marks:
            ax.annotate(name, xy=(val, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(0, 5), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.5, color=INK2)


def fig_exponents(res_tdd, res_ent, out):
    """The central comparison: every exponent we measured, against the two laws."""
    rows = []
    o = res_ent["overall"]
    rows.append(("Entered settings, between people", o["slope"], o["ci_lo"], o["ci_hi"],
                 o["n"], C1))
    for r in res_tdd["between_s_pre"]:
        if r["label"] == "ALL":
            rows.append(("Effective, between people", r["loglog_slope"], r["ci_lo"],
                         r["ci_hi"], r["n_positive"], C2))
    within = {r["label"]: r for r in res_tdd.get("within", [])}
    if "ALL" in within:
        r = within["ALL"]
        rows.append(("Effective, within a person", r["pooled_b"],
                     r["pooled_b"] - 1.96 * r["pooled_se"],
                     r["pooled_b"] + 1.96 * r["pooled_se"], r["n"], C3))
    for key in ("ReplaceBG", "Loop", "IOBP2"):
        if key in within:
            r = within[key]
            rows.append((config.COHORTS[key]["label"] + ", within", r["pooled_b"],
                         r["pooled_b"] - 1.96 * r["pooled_se"],
                         r["pooled_b"] + 1.96 * r["pooled_se"], r["n"], MUTED))

    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    y = np.arange(len(rows))[::-1]
    for yi, (lab, v, lo, hi, n, col) in zip(y, rows):
        ax.plot([lo, hi], [yi, yi], color=col, linewidth=2.2, solid_capstyle="round",
                zorder=3)
        ax.plot([v], [yi], "o", color=col, markersize=8, zorder=4,
                markeredgecolor="white", markeredgewidth=1.2)
        ax.annotate(f"{v:+.2f}", xy=(v, yi), xytext=(0, 9), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8, color=INK)
        ax.text(0.46, yi, f"n={n:,}", fontsize=7, color=INK2, va="center")
    # A rule separating the headline estimates from the per-cohort detail below.
    ax.axhline(len(rows) - 3.5, color=GRID, linewidth=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.5, color=INK)
    ax.set_xlim(-2.4, 0.6)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    _laws(ax)
    ax.set_xlabel("exponent of sensitivity against total daily dose")
    ax.set_title("Estimates cluster between the square root and v1. None comes near v2.",
                 pad=18)
    _despine(ax)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_entered_scatter(out):
    """Entered sensitivity against dose, with the candidate laws drawn over it."""
    t = pd.read_parquet(config.RESULTS / "inv009_entered_isf.parquet")
    t = t.dropna(subset=["isf", "tdd_u"])
    t = t[(t.isf > 0) & (t.tdd_u > 0) & (t.n_days >= config.MIN_DAYS_TDD)]
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.scatter(t.tdd_u, t.isf, s=13, color=C1, alpha=0.38, linewidths=0, zorder=3)
    x = np.linspace(t.tdd_u.min(), t.tdd_u.max(), 200)
    med_tdd = float(t.tdd_u.median())
    anchor = float(t.isf.median())
    for beta, name, col, style in ((-1.0, "v1: 1/TDD", C2, "-"),
                                   (-2.0, "v2: 1/TDD²", RED, "-"),
                                   (-0.5, "square root", C3, "--")):
        y = anchor * (x / med_tdd) ** beta
        ax.plot(x, y, color=col, linewidth=1.8, linestyle=style, zorder=4)
        yy = anchor * (x[-1] / med_tdd) ** beta
        ax.text(x[-1] * 1.02, yy, name, color=col, fontsize=8, va="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("total daily dose (U/day)")
    ax.set_ylabel("entered sensitivity factor (mg/dL per unit)")
    ax.set_title(f"What {len(t)} people had entered, against how much insulin they use")
    ax.set_xlim(right=t.tdd_u.max() * 2.6)
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_synthetic(out, measured=(-0.88, -0.59)):
    """What the pipeline returns when the truth is known.

    The shaded band is what we actually measured. The argument the figure makes
    is geometric: no simulated cohort with a squared law lands anywhere near it.
    """
    s = json.loads((config.RESULTS / "inv009_synthetic.json").read_text())
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.axhspan(measured[0], measured[1], color=C1, alpha=0.10, zorder=1)
    ax.annotate("what we measured", xy=(-2.18, np.mean(measured)), fontsize=8,
                color=C1, va="center", ha="left", fontweight="semibold")
    ax.plot([-2.15, -0.35], [-2.15, -0.35], color=MUTED, linewidth=1.0,
            linestyle="--", zorder=2)
    # Down at the lower left, clear of the series labels crowding the top right.
    ax.annotate("perfect recovery", xy=(-1.85, -1.85), xytext=(4, -14),
                textcoords="offset points", color=INK2, fontsize=7.5, rotation=34,
                rotation_mode="anchor")
    groups = [("open loop", lambda r: r["reactive"] == 0 and r["tag"] == "", C1, "o", 8),
              ("reactive controller", lambda r: r["reactive"] > 0, C2, "s", -14),
              ("45% of meals unrecorded", lambda r: r["tag"] == "harsh", C3, "^", -26)]
    for name, sel, col, mk, dy in groups:
        pts = sorted([(r["beta_true"], r["recovered"]) for r in s if sel(r)])
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax.plot(xs, ys, mk + "-", color=col, markersize=7, linewidth=1.7, zorder=4,
                markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(name, xy=(xs[-1], ys[-1]), xytext=(8, dy),
                    textcoords="offset points", color=col, fontsize=8.5,
                    va="center", fontweight="semibold")
    ax.set_xlabel("exponent the simulated people really had")
    ax.set_ylabel("exponent the pipeline recovered")
    ax.set_title("A squared law cannot come back looking shallow")
    ax.set_xlim(-2.25, 0.55); ax.set_ylim(-2.15, 0.1)
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_head_to_head(out):
    """How wrong each candidate was about the overnight fall."""
    r = json.loads((config.RESULTS / "inv009_head_to_head.json").read_text())
    o = r["overall"]
    names = {"rule_1800": "1800 rule (static)", "fitted_flat": "best single number",
             "v1": "v1 equation", "entered": "entered setting",
             "root_tdd": "355/√TDD", "v2": "v2 equation"}
    items = sorted(o.items(), key=lambda kv: kv[1]["mae"])
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    y = np.arange(len(items))[::-1]
    for yi, (k, v) in zip(y, items):
        col = RED if k == "v2" else C1
        ax.barh(yi, v["mae"], height=0.62, color=col, zorder=3)
        ax.text(v["mae"] + 2, yi, f"{v['mae']:.1f}", va="center", fontsize=8, color=INK)
    ax.set_yticks(y); ax.set_yticklabels([names.get(k, k) for k in [i[0] for i in items]],
                                         fontsize=8.5, color=INK)
    ax.set_xlabel("median error in the predicted overnight fall (mg/dL)")
    ax.set_title("Predicting the night: v2 is four times worse than a static rule")
    _despine(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_glucose(out):
    """Whether any glucose shape earns its place.

    Absolute error, on an axis that starts at zero. Plotting the DIFFERENCES
    instead would put a 0.05 mg/dL spread across the full width of the panel and
    make six indistinguishable candidates look like a ranking. They are
    indistinguishable, and the figure should say so at a glance.
    """
    r = json.loads((config.RESULTS / "inv009_glucose_axis.json").read_text())
    allrow = [x for x in r["shapes_by_study"] if x["label"] == "ALL"][0]
    names = {"flat": "no glucose term", "v1_log": "v1's log scaler",
             "v2_log": "v2's log scaler", "power_k1": "power law k=1",
             "power_k2": "power law k=2", "power_k3": "power law k=3"}
    med = allrow["medians"]
    items = sorted(med.items(), key=lambda kv: kv[1])
    spread = max(med.values()) - min(med.values())
    fig, ax = plt.subplots(figsize=(6.8, 3.2))
    y = np.arange(len(items))[::-1]
    for yi, (k, v) in zip(y, items):
        col = C3 if k == "flat" else C1
        ax.barh(yi, v, height=0.62, color=col, zorder=3)
        ax.text(v + 0.35, yi, f"{v:.2f}", va="center", fontsize=8, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([names.get(k, k) for k in [i[0] for i in items]], fontsize=8.5,
                       color=INK)
    ax.set_xlim(0, max(med.values()) * 1.16)
    ax.set_xlabel("out-of-sample error in the predicted overnight fall (mg/dL)")
    ax.set_title("No glucose shape earns its place, the two equations' included\n"
                 f"(all six within {spread:.2f} mg/dL of each other)", pad=10)
    _despine(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_endogeneity(out):
    """The confound, made visible."""
    r = json.loads((config.RESULTS / "inv009_effective_isf.json").read_text())
    rows = sorted(r["by_study"], key=lambda x: -x["frac_positive_pre"])
    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    x = np.arange(len(rows))
    w = 0.38
    ax.bar(x - w / 2, [100 * v["frac_positive_pre"] for v in rows], w, color=C1, zorder=3)
    ax.bar(x + w / 2, [100 * v["frac_positive_total"] for v in rows], w, color=C2, zorder=3)
    ax.text(-0.5 + 0.06, 96, "insulin already committed", color=C1, fontsize=8.5,
            fontweight="semibold")
    ax.text(-0.5 + 0.06, 88, "all insulin, reactive included", color=C2, fontsize=8.5,
            fontweight="semibold")
    ax.axhline(50, color=MUTED, linestyle=":", linewidth=0.9, zorder=2)
    ax.text(len(rows) - 0.4, 51, "coin flip", fontsize=7, color=INK2, ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels([v["label"].replace(" (", "\n(") for v in rows], fontsize=7.5)
    ax.set_ylabel("share of people with a positive sensitivity (%)")
    ax.set_ylim(0, 104)
    ax.set_title("Insulin a controller chose by watching glucose cannot measure sensitivity")
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


# The power-law exponent each equation's scaler corresponds to over 100 to 200
# mg/dL, so the measured values can be read against them on the same axis.
V1_K, V2_K = 0.62, 1.77


def fig_carb_glucose(out):
    """Apparent glucose dependence against time since the last recorded meal."""
    r = json.loads((config.RESULTS / "inv009_carb_hypothesis.json").read_text())
    arms = [("glucose 90 to 300", "glucose_k_by_gap", C1, 10),
            ("glucose 120 to 220", "glucose_k_by_gap_midrange", C2, -16)]
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for val, name, style in ((V2_K, "v2's scaler", "-"), (V1_K, "v1's scaler", "--"),
                             (0.0, "no dependence", ":")):
        ax.axhline(val, color=INK2 if style == "-" else MUTED, linestyle=style,
                   linewidth=0.9, zorder=1)
        ax.annotate(name, xy=(1.0, val), xycoords=("axes fraction", "data"),
                    xytext=(4, 0), textcoords="offset points", fontsize=7.5,
                    color=INK2, va="center")
    for label, key, col, dy in arms:
        rows = r.get(key, [])
        if not rows:
            continue
        x = np.arange(len(rows))
        y = [v["pooled"] for v in rows]
        e = [1.96 * v["se"] for v in rows]
        ax.errorbar(x, y, yerr=e, fmt="o-", color=col, markersize=7, linewidth=1.7,
                    capsize=3, elinewidth=1.2, zorder=4,
                    markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(label, xy=(x[-1], y[-1]), xytext=(8, dy),
                    textcoords="offset points", color=col, fontsize=8.5,
                    va="center", fontweight="semibold")
        ax.set_xticks(x)
        ax.set_xticklabels([v["band"] for v in rows], fontsize=8)
    ax.set_xlim(-0.4, len(r.get("glucose_k_by_gap", [])) - 0.1 + 1.1)
    ax.set_xlabel("time since the last recorded meal, at the start of the window")
    ax.set_ylabel("apparent glucose dependence (exponent k)")
    ax.set_title("Carbohydrate inflates the glucose effect near a meal, and does not\n"
                 "account for all of it", pad=10)
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_carb_tail(out):
    """Measured sensitivity against time since the last recorded meal."""
    r = json.loads((config.RESULTS / "inv009_carb_hypothesis.json").read_text())
    rows = r["sensitivity_by_gap"]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    x = np.arange(len(rows))
    vals = [v["median_s"] for v in rows]
    cols = [RED if v < 0 else C1 for v in vals]
    ax.bar(x, vals, width=0.6, color=cols, zorder=3)
    for xi, v in zip(x, vals):
        ax.annotate(f"{v:+.2f}", xy=(xi, v), xytext=(0, 4 if v > 0 else -12),
                    textcoords="offset points", ha="center", fontsize=8, color=INK)
    ax.axhline(0, color=INK2, linewidth=0.9, zorder=4)
    ax.set_xticks(x); ax.set_xticklabels([v["band"] for v in rows], fontsize=8)
    ax.set_xlabel("time since the last recorded meal")
    ax.set_ylabel("measured sensitivity (mg/dL per unit)")
    ax.set_title("Within four hours of a meal insulin appears to achieve nothing,\n"
                 "because carbohydrate is still arriving", pad=10)
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


# Blue sequential ramp, light to dark, for the ordered dose bands. Dose is a
# magnitude, so it takes one hue in steps rather than categorical slots.
DOSE_RAMP = ["#9ec5f4", "#5598e7", "#2a78d6", "#184f95"]


def fig_glucose_profile(out):
    """Measured sensitivity by glucose band, against what each equation predicts."""
    r = json.loads((config.RESULTS / "inv009_joint_surface.json").read_text())
    bands = ["90-120", "120-150", "150-190", "190-300"]
    measured = [r["band_profile"][b] for b in bands]
    # Both equations evaluated at a 40 U/day dose, relative to the same reference
    # band, which is where the printed run takes them from.
    v1 = [1.18, 1.00, 0.87, 0.75]
    v2 = [1.75, 1.00, 0.72, 0.54]
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    x = np.arange(len(bands))
    for vals, label, col, mk, dy in ((measured, "what the data show", C1, "o", 0),
                                     (v1, "v1 predicts", C2, "s", 0),
                                     (v2, "v2 predicts", C3, "^", 0)):
        ax.plot(x, vals, mk + "-", color=col, markersize=7, linewidth=2.0, zorder=4,
                markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(label, xy=(x[-1], vals[-1]), xytext=(8, dy),
                    textcoords="offset points", color=col, fontsize=8.5,
                    va="center", fontweight="semibold")
    ax.axhline(1.0, color=MUTED, linestyle=":", linewidth=0.9, zorder=1)
    ax.set_xticks(x); ax.set_xticklabels(bands, fontsize=8.5)
    ax.set_xlim(-0.3, len(bands) - 0.3 + 1.4)
    ax.set_xlabel("glucose at the start of the window (mg/dL)")
    ax.set_ylabel("sensitivity, relative to the 120 to 150 band")
    ax.set_title("Both equations have the glucose axis the wrong way round", pad=10)
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_surface(out):
    """The glucose profile at each dose band. A separable law gives parallel lines."""
    r = json.loads((config.RESULTS / "inv009_joint_surface.json").read_text())
    surf = r["surface"]
    bands = ["90-120", "120-150", "150-190", "190-300"]
    order = ["under 25", "25-40", "40-60", "over 60"]
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    x = np.arange(len(bands))
    for i, tb in enumerate(order):
        if tb not in surf:
            continue
        vals = [surf[tb].get(b) for b in bands]
        ax.plot(x, vals, "o-", color=DOSE_RAMP[i], markersize=6.5, linewidth=2.0,
                zorder=4, markeredgecolor="white", markeredgewidth=1.0)
        ax.annotate(f"{tb} U/day", xy=(x[-1], vals[-1]), xytext=(8, 0),
                    textcoords="offset points", color=DOSE_RAMP[i], fontsize=8.5,
                    va="center", fontweight="semibold")
    ax.set_xticks(x); ax.set_xticklabels(bands, fontsize=8.5)
    ax.set_xlim(-0.3, len(bands) - 0.3 + 1.5)
    ax.set_xlabel("glucose at the start of the window (mg/dL)")
    ax.set_ylabel("sensitivity, as a fraction of that person's own")
    ax.set_title("The glucose profile changes with dose, which a separable law forbids",
                 pad=10)
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_pointwise(out):
    """Sensitivity at each reading against what both equations say at the same reading."""
    r = json.loads((config.RESULTS / "inv009_pointwise.json").read_text())
    rows = r["by_bg"]
    bands = [x["band"] for x in rows]
    x = np.arange(len(bands))
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    series = (("isf_v2", "v2 calculates", C3, "^", 8),
              ("isf_v1", "v1 calculates", C2, "s", -2),
              ("isf_eff", "measured", C1, "o", -16))
    for key, label, col, mk, dy in series:
        vals = [v[key] for v in rows]
        # Values at or below zero cannot be drawn on a log axis and must not be
        # floored onto it either, which would invent a small positive number
        # where the measurement says a unit of insulin moved nothing at all.
        yy = [v if v > 0.5 else np.nan for v in vals]
        ax.plot(x, yy, mk + "-", color=col, markersize=7, linewidth=2.0, zorder=4,
                markeredgecolor="white", markeredgewidth=1.0)
        last = max(i for i, v in enumerate(yy) if np.isfinite(v))
        ax.annotate(label, xy=(x[last], yy[last]), xytext=(8, dy),
                    textcoords="offset points", color=col, fontsize=8.5,
                    va="center", fontweight="semibold")
        first = min((i for i, v in enumerate(yy) if np.isfinite(v)), default=0)
        if first > 0:
            ax.annotate("measured is at or below zero here,\nso it cannot be drawn on a log axis",
                        xy=(x[first] - 0.15, yy[first]), xytext=(6, 14),
                        textcoords="offset points", color=col, fontsize=7.5,
                        ha="left", va="bottom")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(bands, fontsize=8.5)
    ax.set_xlim(-0.3, len(bands) - 0.3 + 1.7)
    ax.set_xlabel("glucose at the start of the six hour lookback (mg/dL)")
    ax.set_ylabel("sensitivity (mg/dL per unit, log scale)")
    ax.set_title("What a unit did, against what each equation said it would do", pad=10)
    _despine(ax)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def fig_detectors(out):
    """Sensitivity detectors, the ones AndroidAPS ships against alternatives."""
    r = json.loads((config.RESULTS / "inv009_detectors.json").read_text())
    med, base = r["medians"], r["medians"]["static"]
    shipped = set(r["shipped"]) | {f"{x}_wide" for x in r["shipped"]}
    names = {"aaps24": "AAPS median, 24 h", "aaps8": "AAPS median, 8 h",
             "oref1_24": "oref1, 24 h", "oref1_8": "oref1, 8 h",
             "wavg24": "weighted average, 24 h",
             "kalman3": "variance weighted", "ewma1": "nightly, 1 night",
             "ewma3": "nightly, 3 nights", "ewma7": "nightly, 7 nights"}

    def label(k):
        if k == "static":
            return "no adjustment"
        wide = k.endswith("_wide") or "0.5-1.5" in k
        stem = k.replace("_wide", "").split("_0.")[0]
        return f"{names.get(stem, stem)}  ({'0.5 to 1.5' if wide else '0.7 to 1.2'})"

    items = [(k, base - v) for k, v in med.items() if k != "static"]
    items.sort(key=lambda kv: kv[1])
    fig, ax = plt.subplots(figsize=(7.6, 6.0))
    y = np.arange(len(items))
    for yi, (k, gain) in zip(y, items):
        col = C1 if k in shipped else C2
        ax.barh(yi, gain, height=0.68, color=col, zorder=3)
        ax.text(gain + 0.05, yi, f"{gain:.2f}", va="center", fontsize=7.5, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([label(k) for k, _ in items], fontsize=8)
    ax.set_xlabel("improvement over a sensitivity that never moves (mg/dL)")
    ax.set_xlim(0, max(g for _, g in items) * 1.14)
    ax.text(0.98, 0.05, "shipped in AndroidAPS", transform=ax.transAxes, ha="right",
            fontsize=8.5, color=C1, fontweight="semibold")
    ax.text(0.98, 0.02, "alternatives tested here", transform=ax.transAxes, ha="right",
            fontsize=8.5, color=C2, fontweight="semibold")
    ax.set_title("The detector AndroidAPS already ships wins. Its clamp is what holds it back.",
                 pad=10)
    _despine(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def main() -> int:
    config.ensure_dirs()
    res_tdd = json.loads((config.RESULTS / "inv009_tdd_axis.json").read_text())
    res_ent = json.loads((config.RESULTS / "inv009_entered_isf.json").read_text())
    jobs = [
        ("fig_exponents.png", lambda p: fig_exponents(res_tdd, res_ent, p)),
        ("fig_entered_scatter.png", fig_entered_scatter),
        ("fig_synthetic.png", fig_synthetic),
        ("fig_head_to_head.png", fig_head_to_head),
        ("fig_glucose.png", fig_glucose),
        ("fig_endogeneity.png", fig_endogeneity),
        ("fig_carb_glucose.png", fig_carb_glucose),
        ("fig_carb_tail.png", fig_carb_tail),
        ("fig_glucose_profile.png", fig_glucose_profile),
        ("fig_surface.png", fig_surface),
        ("fig_pointwise.png", fig_pointwise),
        ("fig_detectors.png", fig_detectors),
    ]
    for name, fn in jobs:
        p = config.CHARTS / name
        fn(p)
        print(f"  {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
