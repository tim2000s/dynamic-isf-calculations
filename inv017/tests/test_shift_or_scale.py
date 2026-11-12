"""Each statistic has to return 1 for the form the falls were built from, and not for the other."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from inv017.shift_or_scale import STATS, constructed, person_sums, sensitivities, stats


def _corrections(n_people: int = 120, per: int = 10) -> pd.DataFrame:
    rng = np.random.default_rng(1)
    isf = rng.uniform(40.0, 100.0, n_people)        # high enough that the shift stays positive
    sid = np.repeat(np.arange(n_people), per)
    start = rng.uniform(150, 350, len(sid))
    return pd.DataFrame(dict(sid=sid, isf=isf[sid], bg_start=start, given=(start - 110) / isf[sid],
                             bg_6h=110.0, target=110.0))


def test_forms_meet_at_the_anchor_and_separate_above_it():
    s = sensitivities(np.array([30.0, 100.0]), np.array([140.0, 140.0]))
    assert all(np.allclose(v, [30.0, 100.0]) for v in s.values())
    s = sensitivities(np.array([30.0, 100.0]), np.array([250.0, 250.0]))
    assert s["Shifted"][0] < s["Multiplied"][0] and s["Shifted"][1] > s["Multiplied"][1]


@pytest.mark.parametrize("form,own,other", [("Multiplied", 1, 2), ("Shifted", 2, 1)])
def test_a_known_form_is_recovered_without_noise(form, own, other):
    corr = _corrections()
    p = person_sums(corr, constructed(corr, form, share=0.0, noise=0.0))
    r = dict(zip(range(len(STATS)), stats(p, np.percentile(p.isf, [25, 50, 75]))))
    assert r[own] == pytest.approx(1.0, abs=1e-9) and abs(r[other] - 1.0) > 0.1
    # multiplying has one steepness for everyone; shifting is steeper where the entered value is low
    assert (r[5] == pytest.approx(1.0, abs=0.05)) if form == "Multiplied" else (r[5] < 0.9)
