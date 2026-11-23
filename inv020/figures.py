"""INV-020 figure, from the per-band file in the git-ignored inv020_cache/.

fig1: the curve's sensitivity over the logarithmic model's by glucose band, median over people
per cohort, with the middle half shaded. Only aggregates are drawn. Palette as inv013/figures.py.

    python3 -m inv020.figures
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv020.replay import CHARTS, COHORTS, G_LABELS, PRIVATE

BLUE, ORANGE, INK, MUTED, SURFACE, GRID = (
    "#2a78d6",
    "#eb6834",
    "#0b0b0b",
    "#52514e",
    "#fcfcfb",
    "#e4e3df",
)
STYLE = {"v6": (BLUE, "o"), "v7": (ORANGE, "s")}


def axes(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)


def isf_ratio_by_glucose(bands: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8), facecolor=SURFACE)
    axes(ax)
    x = np.arange(len(G_LABELS))
    for c, name in COHORTS.items():
        colour, marker = STYLE[c]
        b = bands[(bands.cohort == c) & (bands.n >= 20)].groupby("g").ratio
        med = b.median().reindex(range(len(G_LABELS)))
        q1, q3 = (
            b.quantile(0.25).reindex(med.index),
            b.quantile(0.75).reindex(med.index),
        )
        ax.fill_between(x, q1, q3, color=colour, alpha=0.12, linewidth=0)
        ax.plot(
            x,
            med,
            color=colour,
            linewidth=2,
            marker=marker,
            markersize=8,
            markeredgecolor=SURFACE,
            markeredgewidth=2,
            label=name,
        )
    ax.axhline(1, color=MUTED, linewidth=1, linestyle=(0, (4, 3)))
    ax.text(
        -0.3,
        1.02,
        "same as the logarithmic model",
        color=MUTED,
        fontsize=8,
        va="bottom",
    )
    ax.set_yscale("log")
    ax.set_yticks([0.5, 0.75, 1, 1.5, 2, 3, 4])
    ax.get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")
    )
    ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks(x, G_LABELS)
    ax.set_xlim(-0.4, len(G_LABELS) - 0.6)
    ax.set_xlabel("Glucose at the decision, mg/dL", color=INK)
    ax.set_ylabel("Curve ISF / logarithmic ISF", color=INK)
    ax.set_title(
        "Above 1 the curve reads a higher ISF than the logarithmic model\n"
        "(median over people, shading the middle half)",
        color=INK,
        fontsize=11,
        loc="left",
    )
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(CHARTS / "fig1_isf_ratio_by_glucose.png", dpi=160, facecolor=SURFACE)
    plt.close(fig)


def main() -> None:
    CHARTS.mkdir(parents=True, exist_ok=True)
    isf_ratio_by_glucose(pd.read_csv(PRIVATE / "inv020_bands.csv"))


if __name__ == "__main__":
    main()
