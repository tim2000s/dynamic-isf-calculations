"""Per-person values reach a chart only as bins of at least five people."""

import numpy as np

from inv008.plots import bins


def test_bins_hold_at_least_five_people_and_never_one():
    x = np.arange(23.0)
    b = bins(x, x * 2, n_bins=8)
    assert len(b) == 4 and min(r[4] for r in b) >= 5 and sum(r[4] for r in b) == 23
    assert bins(np.arange(4.0), np.arange(4.0)) == []


def test_histogram_bars_hold_at_least_five_people():
    from inv008.plots import hist_edges

    v = np.r_[np.zeros(12), np.arange(1.0, 9.0), [50.0]]  # sparse tail and one outlier
    e = hist_edges(v, bins=25)
    counts, _ = np.histogram(v, e)
    assert counts.sum() == len(v) and counts.min() >= 5
    assert len(hist_edges(np.arange(4.0))) == 0


def test_histogram_never_shows_four_people_or_an_extreme():
    from inv008.plots import hist_edges

    v = np.r_[
        np.full(10, 1.03), [7.13, 7.2, 7.3, 7.4], np.nan
    ]  # four in the tail must merge, NaN ignored
    e = hist_edges(v, bins=10)
    counts, _ = np.histogram(v[np.isfinite(v)], e)
    assert counts.min() >= 5 and counts.sum() == 14
    assert e[0] < 1.03 and e[-1] > 7.4
    e = hist_edges(np.linspace(3.17, 9.83, 30), bins=np.linspace(3.17, 9.83, 8))
    assert e[0] < 3.17 and e[-1] > 9.83
    v = np.r_[
        20.0, np.arange(25.0, 201.0, 5.0)
    ]  # whole numbers: the minimum sits on the grid
    e = hist_edges(v, bins=20)
    assert e[0] < 20.0 and e[-1] > 200.0
