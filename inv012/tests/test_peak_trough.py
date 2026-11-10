"""The event finder has to return a known sensitivity from a trace built on one."""
from __future__ import annotations

import numpy as np
import pytest

from inv009 import insulin_models as M
from inv012.peak_trough import acted_cum, find_events


def _trace():
    # Glucose climbs to 250, one unit lands at the peak and lowers it by 50 per
    # unit acted, then it climbs again once the dose is spent.
    kern = M.kernel("oref_6h75")
    n, p = 200, 40
    net = np.zeros(n)
    net[p] = 1.0
    acted = acted_cum(net, kern)
    bg = np.where(np.arange(n) < p, 250.0 - 2.0 * (p - np.arange(n)), 250.0 - 50.0 * acted)
    late = np.arange(n) > p + len(kern)
    bg[late] += np.arange(n)[late] - (p + len(kern))
    return np.arange(n) * 5.0, bg, acted, np.zeros(n), p


def test_one_unit_recovers_fifty():
    t, bg, acted, carb, _ = _trace()
    ev = find_events(t, bg, acted, carb)
    assert len(ev) == 1
    assert ev.used[0] == pytest.approx(1.0, abs=0.01)
    assert (ev.g_max[0] - ev.g_min[0]) / ev.used[0] == pytest.approx(50.0, abs=0.5)


def test_carbohydrate_before_the_peak_removes_the_event():
    t, bg, acted, carb, p = _trace()
    carb[p - 12] = 30.0                      # an hour before the peak
    assert find_events(t, bg, acted, carb).empty


def test_a_gap_in_the_fall_removes_the_event():
    t, bg, acted, carb, p = _trace()
    bg[p + 10] = np.nan
    assert find_events(t, bg, acted, carb).empty
