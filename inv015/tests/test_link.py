"""The regression has to return each link when the falls were built from it."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from inv015.link import constructed, fit, people


def _corrections(n_people: int = 80, per: int = 8) -> pd.DataFrame:
    rng = np.random.default_rng(1)
    tdd = np.exp(rng.normal(np.log(40), 0.4, n_people))
    isf = 1915 / tdd * np.exp(rng.normal(0, 0.3, n_people))
    sid = np.repeat(np.arange(n_people), per)
    start = rng.uniform(150, 350, len(sid))
    given = (start - 110) / isf[sid]
    return pd.DataFrame(dict(sid=sid, tdd=tdd[sid], isf=isf[sid], bg_start=start, given=given,
                             bg_6h=start - 0.7 * isf[sid] * given, target=110.0))


@pytest.mark.parametrize("truth,expected", [("daily dose", (-1.0, 0.0)), ("profile", (0.0, 1.0))])
def test_a_known_link_is_recovered_without_noise(truth, expected):
    corr = _corrections()
    b = fit(people(corr, constructed(corr, truth, f0=0.0, noise=0.0)))
    np.testing.assert_allclose(b, expected, atol=1e-6)


def test_an_unrelated_fall_leans_towards_the_profile():
    corr = _corrections()
    b = fit(people(corr, constructed(corr, "daily dose", f0=30.0, noise=0.0)))
    assert -1.0 < b[0] < 0.0 < b[1]
