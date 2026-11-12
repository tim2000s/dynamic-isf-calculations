"""People built as multiples of one curve share a shape, and people built as shifts do not."""
from __future__ import annotations

import numpy as np
import pandas as pd

from inv008.dynisf import quartic_isf
from inv018.person_curves import by_entered, curves, relative, spread

MIDS = np.arange(100, 310, 10)
ENTERED = pd.Series(np.linspace(40.0, 120.0, 60), index=[f"p{i}" for i in range(60)])


def _bins(form) -> pd.DataFrame:
    """inv012's layout: `med_pos` is the reading x daily dose, so a dose of 2 U halves it back."""
    return pd.DataFrame([dict(subject=s, mid=m, tdd=2.0, med_pos=2.0 * form(e, float(quartic_isf(m))))
                         for s, e in ENTERED.items() for m in MIDS])


def test_multiples_of_one_curve_share_a_shape():
    wide = curves(_bins(lambda e, q: e * q / float(quartic_isf(140))))
    assert np.allclose(wide.loc["p0", 140], 40.0)                  # the daily dose is divided out
    s = spread(relative(wide))
    assert np.allclose(s[10], s[90])
    _, (drop, kept) = by_entered(wide, ENTERED)
    assert abs(kept[0] - 1) < 1e-9 and drop[0] > 1.5


def test_shifts_of_one_curve_keep_the_drop_and_change_the_shape():
    wide = curves(_bins(lambda e, q: q + e - float(quartic_isf(140))))
    s = spread(relative(wide))
    assert s[90][250] - s[10][250] > 0.1
    _, (drop, kept) = by_entered(wide, ENTERED)
    assert abs(drop[0] - 1) < 1e-9 and kept[0] > 1.2
