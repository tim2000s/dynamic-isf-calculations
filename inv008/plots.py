"""Per-person values drawn as an aggregate, since per-person results are never published (REVIEW.md
item 22): the people are ordered along x and split into bins of at least five, and each bin is drawn
as its median x and median y with the interquartile range of y."""

from __future__ import annotations

import numpy as np

MIN_PEOPLE = 5


def bins(x, y, n_bins: int = 8, min_people: int = MIN_PEOPLE):
    """(median x, median y, 25th and 75th percentile of y, people) per bin, or [] if too few people."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    k = min(n_bins, len(x) // min_people)
    if k < 1:
        return []
    groups = np.array_split(np.argsort(x, kind="stable"), k)
    return [
        (
            float(np.median(x[g])),
            float(np.median(y[g])),
            float(np.percentile(y[g], 25)),
            float(np.percentile(y[g], 75)),
            len(g),
        )
        for g in groups
    ]


def binned(ax, x, y, color, label, n_bins: int = 8):
    """Draw bins() on ax as points with interquartile bars; returns the number of people drawn."""
    b = bins(x, y, n_bins)
    if b:
        mx, my, lo, hi, n = map(np.array, zip(*b))
        ax.errorbar(
            mx,
            my,
            yerr=[my - lo, hi - my],
            fmt="o-",
            color=color,
            capsize=3,
            lw=1.4,
            label=f"{label} (median and IQR in {len(b)} bins of {int(n.min())} or more people)",
        )
        return int(n.sum())
    return 0


def _nice(w: float) -> float:
    """The 1, 2 or 5 times a power of ten at or above w."""
    e = 10.0 ** np.floor(np.log10(w))
    return float(next(m * e for m in (1, 2, 5, 10) if m * e >= w))


def hist_edges(values, bins=20, min_people: int = MIN_PEOPLE):
    """Histogram edges with every bar holding at least min_people: each bin is widened into the next until it does.

    The outer edges are moved out to a round-number grid, since an edge placed at the data's minimum
    or maximum would show that one person's value.
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if len(v) < min_people:
        return np.array([])
    if np.ndim(bins) == 0:
        w = _nice((v.max() - v.min()) / bins) if v.max() > v.min() else 1.0
        edges = (
            np.arange(np.ceil(v.min() / w) - 1, np.floor(v.max() / w) + 2) * w
        )  # strictly outside both extremes
    else:
        edges = np.asarray(bins, float).copy()
        lo, hi = _nice(edges[1] - edges[0]), _nice(edges[-1] - edges[-2])
        if edges[0] >= v.min():
            edges[0] = (np.ceil(v.min() / lo) - 1) * lo
        if edges[-1] <= v.max():
            edges[-1] = (np.floor(v.max() / hi) + 1) * hi
    out = [edges[0]]
    for e in edges[1:]:
        if ((v >= out[-1]) & (v < e)).sum() >= min_people:
            out.append(e)
    if len(out) == 1:
        return np.array([])
    out[-1] = edges[-1]  # a short last run joins the bar before it
    return np.array(out)


def hist(ax, values, bins=20, **kw):
    """ax.hist on hist_edges(), so no bar shows fewer than five people's values or an extreme."""
    edges = hist_edges(values, bins)
    if len(edges):
        ax.hist(np.asarray(values, float), bins=edges, **kw)
