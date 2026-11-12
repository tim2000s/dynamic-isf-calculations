"""A shared shape must not look repeatable, and shapes that differ between people must."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from inv008.dynisf import quartic_isf
from inv018.split_half import HEIGHT, SHAPE, agreement, halves, measures

MIDS = np.arange(100, 310, 10)


def _half(height, power, seed) -> pd.DataFrame:
    """People x bins: height x (published shape) ** power, with 15% noise in every bin."""
    shape = np.asarray(quartic_isf(MIDS) / quartic_isf(140), float)
    noise = np.random.default_rng(seed).lognormal(0.0, 0.15, (len(height), len(MIDS)))
    return pd.DataFrame(height[:, None] * shape[None, :] ** power[:, None] * noise, columns=MIDS)


@pytest.mark.parametrize("spread,shape_repeats", [(0.0, False), (0.6, True)])
def test_a_shape_repeats_only_where_people_differ_in_it(spread, shape_repeats):
    rng = np.random.default_rng(0)
    height, power = rng.uniform(40, 300, 200), 1 + spread * rng.uniform(-1, 1, 200)
    a, b = measures(_half(height, power, 1)), measures(_half(height, power, 2))
    assert agreement(a[HEIGHT], b[HEIGHT])["r"][0] > 0.8
    r = agreement(a[SHAPE], b[SHAPE])["r"][0]
    assert r > 0.5 if shape_repeats else abs(r) < 0.2


def test_each_way_of_halving_puts_every_event_in_one_half():
    for first in halves(11).values():
        assert first.sum() in (5, 6) and len(first) == 11


def test_the_length_tier_is_the_longest_that_thirty_people_reach():
    from inv018.split_half_cohorts import tier
    assert tier(pd.Series([900] * 40 + [1700] * 29)) == 800
    assert tier(pd.Series([150] * 100)) == 0
