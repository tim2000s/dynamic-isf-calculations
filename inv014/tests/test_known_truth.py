"""The probe has to return a known drop per unit and report a baseline that was moving."""
from __future__ import annotations

import numpy as np
import pytest

from inv014.known_truth import probe, truth_at


def test_probe_returns_drop_per_unit():
    control = np.full(200, 200.0)
    dosed = control.copy()
    dosed[110:] -= 10.0 * np.minimum(1.0, np.arange(90) / 30.0)      # 10 mg/dL lower after 150 min
    p = probe(control, dosed, 100, 0.5, 5.0)
    assert p["sens"] == pytest.approx(20.0)
    assert p["level"] == 200.0 and p["drift"] == 0.0
    assert p["at_h"] == pytest.approx(40 * 5 / 60)


def test_probe_reports_a_moving_baseline():
    control = 200.0 + 0.2 * np.arange(200)
    assert probe(control, control - 1.0, 100, 1.0, 5.0)["drift"] == pytest.approx(0.2 * 72)


def test_truth_is_not_extrapolated():
    out = truth_at([300.0, 120.0, 200.0], [18.0, 10.0, 14.0], (100, 160, 250, 320))
    assert np.isnan(out[0]) and np.isnan(out[3])
    assert out[1] == pytest.approx(12.0) and out[2] == pytest.approx(16.0)
