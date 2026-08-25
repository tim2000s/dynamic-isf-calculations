"""Total daily dose, as a subject-level property and as the loop would see it.

Two different quantities go by that name and the difference is the whole point
of the exercise. The subject-level dose is how much insulin this person uses,
and it is what a between-person law like the 1800 rule is about. The windowed
dose is what a controller would have computed at that moment from recent
delivery, and it is what a dynamic equation would actually be fed. A law that
holds between people need not hold within one, so both are carried.

The windowed side reuses INV-008's calculator unchanged so the two
investigations can be read against each other.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inv008 import dynisf, tdd_windows

from . import config, grid as gridmod


def windowed(grid: pd.DataFrame) -> pd.DataFrame:
    """Per-bin trailing dose components and the blend the equations use."""
    g = pd.DataFrame({"ts": gridmod.epoch_s(grid.ts),
                      "total_u": grid.total_u.to_numpy(),
                      "bolus_u": grid.bolus_u.to_numpy()})
    w = tdd_windows.windowed_tdd(g, bin_sec=int(config.GRID_MIN * 60))
    w["tdd_blend"] = dynisf.blend_tdd(w.tdd_4h, w.tdd_8to4h, w.tdd_1d, w.tdd_7d)
    return w


def subject_level(grid: pd.DataFrame) -> dict:
    """This person's usual daily dose, over days the pump was actually recording.

    A day missing bins is a day of upload gap, and counting its partial total
    would bias the average down for exactly the people whose data is patchiest.
    """
    d = gridmod.daily_totals(grid)
    ok = d[d.complete]
    if len(ok) < config.MIN_DAYS_TDD:
        return dict(tdd_u=np.nan, tdd_basal_u=np.nan, basal_frac=np.nan, n_days=len(ok))
    return dict(tdd_u=float(ok.total_u.mean()),
                tdd_basal_u=float(ok.basal_u.mean()),
                basal_frac=float(ok.basal_u.sum() / max(ok.total_u.sum(), 1e-9)),
                n_days=int(len(ok)))
