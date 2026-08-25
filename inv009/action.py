"""How much insulin acted during a window, and how much of that was decided in advance.

This is the estimator's engine. There is no insulin on board column in this data
and no loop prediction to lean on, so everything is reconstructed from delivered
basal and boluses by convolving them with an insulin model.

The split matters more than the total. Insulin already in the body when a window
opens was decided before anything in that window happened, so it cannot have been
a reaction to the glucose the window goes on to measure. Insulin delivered during
the window can be, and under a closed loop usually is. Keeping the two apart is
what lets a closed-loop cohort be read at all: the pre-window term is
predetermined, the in-window term is not.

All quantities are UNITS OF INSULIN THAT ACTED inside the window, so a regression
of glucose fall on them has a slope in mg/dL per unit, which is what a sensitivity
factor is.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config


def iob_series(u: np.ndarray, kern: np.ndarray) -> np.ndarray:
    """Insulin on board at each bin, from a delivery series and a model kernel."""
    return np.convolve(u, kern)[:len(u)]


def window_action(u: np.ndarray, kern: np.ndarray, horizon_bins: int):
    """Units acting inside each window, split by whether they were given before it.

    Returns (a_pre, a_in), both indexed by the window's opening bin.

    `a_pre` is what insulin given strictly before the window contributed while
    the window ran. `a_in` is what insulin given during the window contributed
    before the window closed. Together they are all the insulin action the
    window saw.

    The identities used, with K the fraction-remaining kernel and IOB its
    convolution with delivery:

        a_pre[i] = sum_{j<i}      u[j] (K[i-j] - K[i+H-j])
                 = IOB[i] - u[i] - IOB[i+H] + C[i]
        a_in[i]  = sum_{i<=j<=i+H} u[j] (1 - K[i+H-j])
                 = U_in[i] - C[i]

    where C[i] = sum_{i<=j<=i+H} u[j] K[i+H-j], which is the convolution of
    delivery with the kernel truncated to the horizon, read at bin i+H. Writing
    it this way keeps the whole thing to three convolutions over the record
    rather than one pass per window.
    """
    n = len(u)
    h = int(horizon_bins)
    iob = np.convolve(u, kern)[:n]
    c_full = np.convolve(u, kern[:h + 1])[:n]

    def shift(x):
        """x[i+H], zero past the end."""
        out = np.zeros(n)
        if h < n:
            out[:n - h] = x[h:]
        return out

    csum = np.concatenate([[0.0], np.cumsum(u)])
    hi = np.minimum(np.arange(n) + h + 1, n)
    u_in = csum[hi] - csum[np.arange(n)]

    c = shift(c_full)
    a_in = u_in - c
    a_pre = iob - u - shift(iob) + c
    return a_pre, a_in


def reference_profile(ts: pd.Series, u: np.ndarray, days: int = 30) -> np.ndarray:
    """The delivery a subject usually gets at this time of day, in units per bin.

    Basal need is never observed, and overnight it is most of the insulin, so a
    regression on total delivery is asking a question the data cannot answer: the
    dose barely varies from night to night and the intercept absorbs the level.
    Subtracting the subject's own usual pattern turns the regressor into the
    deviation from routine, which is the part that does vary.

    Half-hour of day, median over a trailing window of days, so a profile change
    partway through a record is followed rather than averaged across.
    """
    half_hour = (ts.dt.hour * 2 + ts.dt.minute // 30).to_numpy()
    day = ts.dt.normalize()
    day_idx = ((day - day.iloc[0]).dt.days).to_numpy()
    block = day_idx // days                       # trailing blocks, not a true rolling median
    ref = np.zeros(len(u), dtype=float)
    frame = pd.DataFrame({"block": block, "hh": half_hour, "u": u})
    med = frame.groupby(["block", "hh"], sort=False)["u"].transform("median")
    ref[:] = med.to_numpy()
    return ref


def net_delivery(u: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Delivery as a deviation from the subject's usual pattern."""
    return u - ref
