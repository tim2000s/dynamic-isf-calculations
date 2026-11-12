#!/usr/bin/env python3
"""Standalone figures for the proposal (document 8).

Each figure is written to charts/proposal8/ as a PNG that carries its own title, a line
saying what is plotted and a source line, so it can be read without the document.

Figure 1 redraws the results figure of the published poster (Romero-Ugalde et al.,
doi:10.2337/db25-920-P). The two curves are the published equations. The poster gives its
bin medians and quartiles only as a drawing, so they are measured here from the copy of
that drawing held in glucose-isf/: the gridlines calibrate both axes, the top of each bar
gives the median and the caps of each whisker give the quartiles. One pixel is 0.35 mg/dL
per unit. The measurement is checked against the published equations, which were fitted
to those medians, before anything is drawn, and results/inv013_poster_figure.md holds the
values. At 60 and 80 mg/dL both quartiles lie beyond the poster's axis, and at 100 mg/dL
the upper one does; those are left blank.

Figures 2 to 4 draw the proposal's rule and the evidence for it. They compute nothing
new: the models come from `inv013.candidate`, and every measured number is read from the
committed results file of the investigation that produced it (INV-015, INV-017), with an
assertion that the row is there.

Figure 5 draws the outcome test in simulation (INV-019) from results/inv019_outcome_test.md:
each rule's mean paired difference from the same simulated people on their unchanged
profile ISF, with its 95% interval. The harm-rule readings in its heading are read from the
same file. TimSim's results are relative, and the figure says so in its source line.

Modelling choices shown, each made in the proposal and none in the publication: the
profile ISF taken to hold at 140 mg/dL, the worked person (profile ISF 50 mg/dL per unit,
40 U a day), AndroidAPS at an adjustment factor of 100% with a divisor of 55, and Trio's
sigmoid at its defaults.

Run as `python3 -m inv013.figures`.
"""
from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

from inv008.dynisf import quartic_isf
from inv013.candidate import (ANCHOR_GLUCOSE, PROFILE_ISF, RESULTS, TDD, curve_profile, published,
                              stock_v1, trio_sigmoid)

ROOT = Path(__file__).resolve().parents[1]
POSTER = ROOT / "glucose-isf" / "Variable ISF.png"
CHARTS = ROOT / "charts" / "proposal8"
BLUE, ORANGE, AQUA, INK, MUTED, SURFACE, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
BINS = np.arange(60, 381, 20)
CITE = "Romero-Ugalde et al., Diabetes 2025;74(Suppl 1):920-P, doi:10.2337/db25-920-P"
TRIO_LIMITS = (0.7, 1.2)


def _grid_lines(mask: np.ndarray) -> np.ndarray:
    """Centres of the runs of True in a 1-D mask."""
    i = np.where(mask)[0]
    return np.array([r.mean() for r in np.split(i, np.where(np.diff(i) > 1)[0] + 1)])


def digitise_poster(path: Path = POSTER) -> pd.DataFrame:
    """Bin medians and quartiles measured from the poster's results figure, mg/dL per unit."""
    im = np.asarray(Image.open(path).convert("RGB")).astype(int)
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    grey = (abs(r - g) < 8) & (abs(g - b) < 8) & (r > 150) & (r < 225)
    bar = (abs(r - 88) < 25) & (abs(g - 128) < 25) & (abs(b - 182) < 25)
    # whiskers: the dark core, which survives inside a bar, or a thin cap blended towards white
    navy = ((r < 90) & (g < 90) & (b > 90) & (b - r > 40)) | ((abs(r - g) < 15) & (b - r > 30) & (r < 170))
    # gridlines, read in a strip clear of the bars: left of the first bar, and below zero
    # The frame is caught as well: its bottom edge is the last row, its sides the first and last columns.
    rows = _grid_lines(grey[:, 166:184].all(axis=1))[:-1]
    cols = _grid_lines(grey[406:428, :].all(axis=0))[1:-1]
    value = 20 * np.round((rows[-1] - rows) / (rows[-1] - rows[-2]))     # the lowest two lines are 0 and 20
    assert len(rows) >= 5 and value.max() == 120 and len(cols) == 8, (rows, cols)   # 50 up to 400
    fy = np.polyfit(rows, value, 1)
    fx = np.polyfit(np.arange(50, 401, 50), cols, 1)
    assert np.abs(np.polyval(fy, rows) - value).max() < 0.5
    legend_bottom = 140                                         # the legend holds a swatch in the bar colour
    out = []
    for glucose in BINS:
        c = int(round(np.polyval(fx, glucose)))
        tops = [np.where(bar[legend_bottom:, c + d])[0].min() + legend_bottom for d in (-7, -6, 6, 7)]
        top = float(np.median(tops))
        wide = np.where(navy[:, c - 9:c + 10].sum(axis=1) >= 8)[0]   # a cap is wider than the marker
        caps = [k for k in _grid_lines(np.isin(np.arange(len(navy)), wide)) if abs(k - top) > 6] if len(wide) else []
        above, below = [k for k in caps if k < top], [k for k in caps if k > top]
        assert len(above) <= 1 and len(below) <= 1, (glucose, caps)
        val = lambda k: float(np.polyval(fy, k))
        out.append(dict(glucose=int(glucose), median=val(top),
                        q25=val(below[0]) if below else np.nan, q75=val(above[0]) if above else np.nan))
    d = pd.DataFrame(out)
    d["published_fit"] = [published(float(x)) for x in d.glucose]
    return d


def _style(ax, xlabel: str, ylabel: str) -> None:
    ax.set_facecolor(SURFACE)
    ax.set_xlabel(xlabel, color=MUTED)
    ax.set_ylabel(ylabel, color=MUTED)
    ax.tick_params(colors=MUTED, length=0)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_visible(False)


def _finish(fig, title: str, what: str, source: str, name: str) -> Path:
    fig.text(0.012, 0.975, title, fontsize=13, color=INK, va="top", weight="medium")
    fig.text(0.012, 0.925, what, fontsize=9.5, color=MUTED, va="top", linespacing=1.4)
    fig.text(0.012, 0.012, source, fontsize=7.5, color=MUTED, va="bottom", linespacing=1.4)
    CHARTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(CHARTS / name, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return CHARTS / name


def _end_labels(ax, x: float, items: list, gap: float) -> None:
    """Label line ends at x, pushed apart by at least `gap` so that neighbouring labels stay readable."""
    at = -np.inf
    for y, text in sorted(items):
        at = max(y, at + gap)
        ax.annotate(text, (x, at), xytext=(5, 0), textcoords="offset points", va="center", fontsize=9, color=INK)


def fig_published(d: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(9, 5.6), facecolor=SURFACE)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.8, bottom=0.17)
    lo, hi = -10, 120                                           # the poster's own axis
    ax.bar(d.glucose, d["median"], width=9, color=BLUE, alpha=0.75, lw=0, label="Median of the bin")
    q25, q75 = d.q25.fillna(lo), d.q75.fillna(hi)
    ax.vlines(d.glucose, q25, q75, color=INK, lw=1)
    for q in (d.q25, d.q75):
        ax.hlines(q.dropna(), d.glucose[q.notna()] - 3, d.glucose[q.notna()] + 3, color=INK, lw=1)
    ax.plot([], [], color=INK, lw=1, marker="_", ms=7, label="25th to 75th percentile")
    low, high = np.linspace(50, 100, 60), np.linspace(100, 380, 300)
    ax.plot(low, [published(x) for x in low], color=ORANGE, lw=2, ls=":", label="Published fit up to 100 mg/dL (quadratic)")
    ax.plot(high, quartic_isf(high), color=ORANGE, lw=2, ls="--",
            label="Published fit from 100 to 380 mg/dL (quartic)")
    for x in (100, 180, 300):
        m = float(d.loc[d.glucose == x, "median"].iloc[0])
        ax.annotate(f"{m:.0f}", (x, m), xytext=(7, 9), textcoords="offset points", fontsize=9, color=INK)
    ax.set_ylim(lo, hi)
    ax.set_xlim(45, 395)
    _style(ax, "Glucose, mg/dL", "Sensitivity, mg/dL per U")
    ax.legend(frameon=False, fontsize=9, labelcolor=MUTED, loc="upper right")
    return _finish(
        fig, "The published curve: insulin sensitivity against glucose in 11,201 adults",
        "Redrawn from the poster's results figure. Bars are the median sensitivity in each 20 mg/dL bin of glucose, whiskers the\n"
        "quartiles across events. Labelled bars are the medians at 100, 180 and 300 mg/dL.",
        f"Source: {CITE}. Medians and quartiles measured from the published figure (inv013/figures.py);\n"
        "curves are the published equations. Quartiles at 60 and 80 mg/dL, and the upper quartile at 100, lie beyond the poster's axis.",
        "fig1_published_curve.png")


def fig_scaled() -> Path:
    g = np.linspace(70, 380, 311)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5.8), facecolor=SURFACE)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.78, bottom=0.2, wspace=0.16)
    for profile in (30.0, 50.0, 80.0):
        a.plot(g, curve_profile(g, profile), color=BLUE, lw=2)
        a.plot([ANCHOR_GLUCOSE], [profile], "o", color=BLUE, ms=8, mec=SURFACE, mew=2)
    _end_labels(a, 380, [(float(curve_profile(380.0, p)), f"profile ISF {p:.0f}") for p in (30.0, 50.0, 80.0)], 7)
    a.set_title("Three people: the same curve through each profile ISF at 140 mg/dL", loc="left", fontsize=10.5, color=INK)
    lo, hi = PROFILE_ISF / TRIO_LIMITS[1], PROFILE_ISF / TRIO_LIMITS[0]
    b.axhspan(lo, hi, color=MUTED, alpha=0.1, lw=0)
    b.annotate("range Trio's default limits allow (0.7 to 1.2)", (378, hi), xytext=(0, -11), textcoords="offset points",
               ha="right", fontsize=8.5, color=MUTED)
    series = [("Proposed: profile ISF x published curve", "Proposed", curve_profile(g, PROFILE_ISF), BLUE, "-"),
              ("AndroidAPS logarithmic, factor 100%", "AndroidAPS", stock_v1(g, TDD, 55.0), ORANGE, "-"),
              ("Trio sigmoid, defaults", "Trio sigmoid", trio_sigmoid(g, PROFILE_ISF), AQUA, "-"),
              ("Profile ISF, unchanged", "Profile ISF", np.full_like(g, PROFILE_ISF), MUTED, "--")]
    for label, _, y, colour, ls in series:
        b.plot(g, y, color=colour, lw=2, ls=ls, label=label)
    _end_labels(b, 380, [(float(y[-1]), short) for _, short, y, _, _ in series], 5)
    b.set_title(f"One person (profile ISF {PROFILE_ISF:.0f}, {TDD:.0f} U a day) under each model", loc="left", fontsize=10.5, color=INK)
    for ax in (a, b):
        _style(ax, "Current glucose, mg/dL", "Sensitivity used for dosing, mg/dL per U")
        ax.set_xlim(70, 440)                                    # room for the labels at the line ends
        ax.set_xticks(range(100, 351, 50))
        ax.set_ylim(0, 135 if ax is a else 90)
    b.set_ylabel("")
    fig.legend(*b.get_legend_handles_labels(), frameon=False, fontsize=9, labelcolor=MUTED, ncol=4,
               loc="lower right", bbox_to_anchor=(0.985, 0.075))
    return _finish(
        fig, "The proposed rule: the published curve scaled to each person's profile ISF",
        f"Sensitivity = profile ISF x Q(glucose) / Q({ANCHOR_GLUCOSE:.0f}), where Q is the published curve. The profile ISF is taken to hold at\n"
        f"{ANCHOR_GLUCOSE:.0f} mg/dL, a modelling choice of the proposal. Below 100 mg/dL the value at 100 is held.",
        f"Published curve: {CITE}. The link to the profile, the anchor and the worked person are this proposal's choices.\n"
        "AndroidAPS: 1800 / (TDD x ln(glucose / 55 + 1)). Drawn by inv013/figures.py from inv013/candidate.py.",
        "fig2_scaled_to_profile.png")


def _pair(cell: str) -> tuple[float, float, float]:
    v = [float(x) for x in re.findall(r"-?\d+\.\d+", cell)]
    assert len(v) == 3, cell
    return v[0], v[1], v[2]


def link_rows(path: Path = RESULTS / "inv015_link.md") -> tuple[list, tuple]:
    """Constructed rows (true link, dose-unrelated fall, b_dose, b_profile) and the measured row of INV-015."""
    built, measured = [], None
    for row in path.read_text().splitlines():
        c = [x.strip() for x in row.split("|")]
        if len(c) == 7 and c[1] in ("daily dose", "profile"):
            built.append((c[1], c[2], _pair(c[4]), _pair(c[5])))
        if len(c) == 5 and c[1] == "Measured":
            measured = (_pair(c[2]), _pair(c[3]))
    assert len(built) == 8 and measured, (len(built), measured)
    return built, measured


def fig_link() -> Path:
    built, (md, mp) = link_rows()
    people, corrections = re.search(r"(\d+) adults with 5 or more corrections, (\d+) corrections",
                                    (RESULTS / "inv015_link.md").read_text()).groups()
    fig, ax = plt.subplots(figsize=(8.4, 6.6), facecolor=SURFACE)
    fig.subplots_adjust(left=0.09, right=0.97, top=0.8, bottom=0.17)
    leader = dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=2, shrinkB=9)
    for (x, y, label, colour, at) in ((-1, 0, "The daily-dose link predicts (-1, 0)", ORANGE, (-1.0, -0.27)),
                                      (0, 1, "The profile link\npredicts (0, 1)", BLUE, (0.3, 1.3))):
        ax.plot([x], [y], "o", ms=17, mfc="none", mec=colour, mew=2)
        ax.annotate(label, (x, y), xytext=at, ha="center", va="center", fontsize=9, color=INK, arrowprops=leader)
    for true, _, d, p in built:
        colour = ORANGE if true == "daily dose" else BLUE
        ax.plot([d[1], d[2]], [p[0], p[0]], color=colour, lw=0.8, alpha=0.7)
        ax.plot([d[0], d[0]], [p[1], p[2]], color=colour, lw=0.8, alpha=0.7)
        ax.plot([d[0]], [p[0]], "o", ms=6.5, color=colour, mec=SURFACE, mew=1.5)
    ax.plot([], [], "o", ms=6.5, color=ORANGE, label="Constructed falls, true link is daily dose")
    ax.plot([], [], "o", ms=6.5, color=BLUE, label="Constructed falls, true link is profile ISF")
    ax.plot([md[1], md[2]], [mp[0], mp[0]], color=INK, lw=1.6)
    ax.plot([md[0], md[0]], [mp[1], mp[2]], color=INK, lw=1.6)
    ax.plot([md[0]], [mp[0]], "D", ms=9, color=INK, mec=SURFACE, mew=2, label="Measured, with 95% interval over people")
    ax.annotate(f"Measured: {md[0]:.2f} on daily dose,\n{mp[0]:.2f} on entered ISF", (md[0], mp[0]), xytext=(-0.72, 1.22),
                ha="center", va="center", fontsize=9, color=INK, arrowprops=leader)
    ax.set_xlim(-1.25, 0.5)
    ax.set_ylim(-0.35, 1.45)
    _style(ax, "Coefficient on daily dose (the daily-dose link predicts -1)", "Coefficient on entered ISF (the profile link predicts 1)")
    ax.legend(frameon=False, fontsize=8.5, labelcolor=MUTED, loc="lower right")
    return _finish(
        fig, "Which scale corrections followed: the entered ISF or the daily dose",
        f"{int(corrections):,} calculator corrections by {people} adults in REPLACE-BG, where no algorithm doses. The fall per unit is regressed on\n"
        "both scales at once. Small points show what the same regression returns when the true link is known, with a growing share\n"
        "of the fall unrelated to the dose (none, 30 mg/dL, 46%, 75%).",
        "Source: results/inv015_link.md (inv015/link.py). The fall per unit is an outcome ratio, put on the footing of 140 mg/dL by the\n"
        "published shape. Its level is not read. Intervals are 95% over people.",
        "fig3_profile_or_daily_dose.png")


def fig_multiply_or_shift() -> Path:
    text = (RESULTS / "inv017_shift_or_scale.md").read_text()
    people, corrections = re.search(r"(\d+) adults with 5 or more corrections, (\d+) corrections", text).groups()
    zero, zero_people = re.search(r"zero or less at (\d+) of them, in (\d+) people", text).groups()
    g = np.linspace(100, 380, 281)
    q, q140 = quartic_isf(g), float(quartic_isf(ANCHOR_GLUCOSE))
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5.6), facecolor=SURFACE, sharey=True)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.78, bottom=0.17, wspace=0.06)
    ends: dict = {}
    for profile in (20.0, 35.0, 50.0, 80.0):
        for ax, y in ((a, profile * q / q140), (b, q + profile - q140)):
            ax.plot(g, y, color=BLUE, lw=2)
            ax.plot([ANCHOR_GLUCOSE], [profile], "o", color=BLUE, ms=7, mec=SURFACE, mew=2)
            ends.setdefault(ax, []).append((float(y[-1]), f"{profile:.0f}"))
    b.axhspan(-40, 0, color=ORANGE, alpha=0.12, lw=0)
    b.annotate("zero or below: no dose can be calculated", (372, -28), ha="right", va="center", fontsize=9, color=INK)
    for ax in (a, b):
        _end_labels(ax, 380, ends[ax], 8)
    a.set_title("Multiplied: profile ISF x Q(glucose) / Q(140)", loc="left", fontsize=10.5, color=INK)
    b.set_title("Shifted: Q(glucose) + profile ISF - Q(140)", loc="left", fontsize=10.5, color=INK)
    for ax in (a, b):
        _style(ax, "Current glucose, mg/dL", "Sensitivity, mg/dL per U")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xlim(100, 400)                                   # room for the labels at the line ends
        ax.set_xticks(range(100, 351, 50))
        ax.set_ylim(-40, 130)
    b.set_ylabel("")
    return _finish(
        fig, "Two ways to attach the published curve to a profile ISF",
        "Each line is one profile ISF (20, 35, 50 and 80 mg/dL per U, labelled at the right), taken to hold at 140 mg/dL. Multiplying keeps\n"
        f"every person above zero. Shifting gave zero or less at {zero} of {int(corrections):,} corrections, in {zero_people} of {people} adults in REPLACE-BG.",
        f"Published curve Q: {CITE}. Count of corrections from results/inv017_shift_or_scale.md\n(inv017/shift_or_scale.py). Drawn by inv013/figures.py.",
        "fig4_multiply_or_shift.png")


def outcome_rows(path: Path = RESULTS / "inv019_outcome_test.md") -> dict:
    """INV-019, each arm against `profile, pinned`: arm -> {column: (mean, low, high)}."""
    lines = path.read_text().splitlines()
    head = next(i for i, row in enumerate(lines) if row.startswith("| Arm | Reference |"))
    cols = [c.strip() for c in lines[head].split("|")][4:-1]
    out = {}
    for row in lines[head + 2:]:
        c = [x.strip() for x in row.split("|")]
        if len(c) != len(cols) + 5:
            break
        if c[2] == "profile, pinned":                           # the mis-set candidates have their own reference
            out[c[1]] = {k: _pair(v) for k, v in zip(cols, c[4:-1])}
    assert len(out) == 13 and "candidate" in out, list(out)
    return out


def harm_rows(path: Path = RESULTS / "inv019_outcome_test.md") -> dict:
    """INV-019's harm rule: (arm, against, side) -> (difference with interval, reading)."""
    out = {}
    for row in path.read_text().splitlines():
        c = [x.strip() for x in row.split("|")]
        if len(c) == 10 and c[3] in ("below 54", "above 250"):
            out[(c[1], c[2], c[3])] = (_pair(c[5]), c[8])
    assert len(out) == 8 and all(r in ("harm", "cleared", "unresolved") for _, r in out.values()), out
    return out


OUTCOME_ARMS = [   # (arm in the report, label, group), drawn top to bottom
    ("candidate", "Proposed rule (profile ISF at 140 mg/dL)", 0),
    ("candidate, Trio limits", "Proposed rule inside Trio's limits, 0.7 to 1.2", 0),
    ("candidate, anchor 180", "Proposed rule, profile ISF at 180 mg/dL", 0),
    ("equation on TDD", "Published curve scaled by daily dose", 0),
    ("flat above 100", "Profile ISF above 100 mg/dL, curve below", 0),
    ("profile, autosens", "oref0 with autosens, no dynamic rule", 1),
    ("logarithmic, AndroidAPS", "Logarithmic, AndroidAPS at 100%", 1),
    ("logarithmic, matched", "Logarithmic, matched to the profile at 140", 1),
    ("logarithmic, Trio", "Logarithmic, Trio (factor 0.8, limits)", 1),
    ("Trio sigmoid", "Trio sigmoid, defaults", 1),
    ("placebo", "Placebo: profile ISF x 1.01", 2),
    ("profile x 0.7", "Profile ISF entered 30% too low", 2),
    ("profile x 1.4", "Profile ISF entered 40% too high", 2),
]
OUTCOME_GROUPS = [("The published curve, in the forms tested", BLUE), ("Models in use today", ORANGE),
                  ("Checks: a null arm and a mis-set profile", MUTED)]


def fig_outcomes() -> Path:
    rows, harm = outcome_rows(), harm_rows()
    text = (RESULTS / "inv019_outcome_test.md").read_text()
    people, days = re.search(r"(\d+) adults, (\d+) days", text).groups()
    panels = [("tbr54", "Time below 54 mg/dL"), ("tar250", "Time above 250 mg/dL"), ("tir", "Time in range, 70 to 180 mg/dL")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 6.9), facecolor=SURFACE, sharey=True)
    fig.subplots_adjust(left=0.255, right=0.985, top=0.79, bottom=0.2, wspace=0.13)
    at = [-(i + 0.8 * group) for i, (_, _, group) in enumerate(OUTCOME_ARMS)]   # a gap between the groups
    for ax, (col, title) in zip(axes, panels):
        for y, (arm, _, group) in zip(at, OUTCOME_ARMS):
            mean, low, high = rows[arm][col]
            colour = OUTCOME_GROUPS[group][1]
            ax.plot([low, high], [y, y], color=colour, lw=2, solid_capstyle="round")
            ax.plot([mean], [y], "o", ms=8, color=colour, mec=SURFACE, mew=2)
            if arm == "candidate":
                ax.annotate(f"{mean:+.2f}", (mean, y), xytext=(0, 9), textcoords="offset points", ha="center", fontsize=9, color=INK)
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_title(title, loc="left", fontsize=10.5, color=INK)
        _style(ax, "Change, percentage points of time", "")
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(at, [label for _, label, _ in OUTCOME_ARMS])
    axes[0].tick_params(axis="y", labelcolor=INK, labelsize=9)
    axes[0].set_ylim(min(at) - 0.7, 0.9)
    fig.legend([plt.Line2D([], [], marker="o", color=c, lw=2, ms=8, mec=SURFACE, mew=2) for _, c in OUTCOME_GROUPS],
               [name for name, _ in OUTCOME_GROUPS], frameon=False, fontsize=9, labelcolor=MUTED, ncol=3,
               loc="lower right", bbox_to_anchor=(0.985, 0.075))
    (auto, auto_reading), (matched, matched_reading) = (harm[("candidate", who, "below 54")] for who in ("profile, autosens", "logarithmic, matched"))
    high = {harm[("candidate", who, "above 250")][1] for who in ("profile, autosens", "logarithmic, matched")}
    assert len(high) == 1, high                                 # the heading says one reading for both comparators
    return _finish(
        fig, "Glucose outcomes under each sensitivity rule, in simulation",
        f"{people} simulated adults for {days} days on oref0, with only the sensitivity rule changed. Points are the mean difference from the same people on their\n"
        "unchanged profile ISF (days 8 to 35), lines the 95% interval over people. Harm rule for the proposed rule, time below 54 mg/dL: "
        f"{auto_reading} against\noref0 with autosens ({auto[0]:+.2f} points, {auto[1]:.2f} to {auto[2]:.2f}), {matched_reading} against the matched logarithmic model "
        f"({matched[0]:+.2f}, {matched[1]:.2f} to {matched[2]:.2f}). Time above 250 mg/dL: {high.pop()} against both.",
        "Source: results/inv019_outcome_test.md (inv019/outcome_test.py), TimSim on simglucose, one seed. Results are relative. In this simulator the true sensitivity\n"
        "rises with glucose where the published curve's falls, and the titrated profile ISF x daily dose is about half of what adults enter.",
        "fig5_outcomes_in_simulation.png")


def main() -> None:
    d = digitise_poster()
    miss = (d["median"] - d.published_fit).abs()
    L = ["# INV-013: the poster's results figure, measured from its image", "",
         "Generated by `inv013/figures.py` from `glucose-isf/Variable ISF.png`. Sensitivity in mg/dL per unit, by 20 mg/dL",
         "bin of glucose. One pixel of the image is 0.35 mg/dL per unit. A blank quartile lies beyond the poster's axis.", "",
         f"The measured medians sit within {miss.max():.1f} mg/dL per unit of the published equations (median difference {miss.median():.1f}).", "",
         "| Glucose, mg/dL | Median | 25th percentile | 75th percentile | Published equation |", "|---:|---:|---:|---:|---:|"]
    f = lambda v: "" if np.isnan(v) else f"{v:.1f}"
    L += [f"| {r.glucose} | {r.median:.1f} | {f(r.q25)} | {f(r.q75)} | {r.published_fit:.1f} |" for r in d.itertuples()]
    (RESULTS / "inv013_poster_figure.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    for path in (fig_published(d), fig_scaled(), fig_link(), fig_multiply_or_shift(), fig_outcomes()):
        print(path)


if __name__ == "__main__":
    main()
