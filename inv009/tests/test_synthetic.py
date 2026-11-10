"""A meal simulated as unlogged has to be missing from the record, or no screen is tested by it."""
from __future__ import annotations

import numpy as np

from inv009 import synthetic as S


def test_unlogged_meals_leave_no_carbohydrate_entry():
    # Until 10 November 2026 every unlogged meal was also written to carbs_g and
    # absorbed at twice its size, so the arms built on it hid nothing.
    g = S.simulate_subject(np.random.default_rng(1), tdd=40.0, isf=50.0, days=20, unannounced=1.0)
    assert g.carbs_g.sum() == 0.0
    assert g.cgm.max() > 200.0                 # eaten all the same
    logged = S.simulate_subject(np.random.default_rng(1), tdd=40.0, isf=50.0, days=20)
    assert logged.carbs_g.sum() > 0.0
