"""Insulin models, as the apps that made this data implemented them.

Two families are needed. The Loop cohort ran Loop, so it gets LoopKit's models
and not an approximation of them: the exponential presets and the Walsh curves
its era also offered. Everything else is oref-family or a commercial controller,
and gets the oref exponential.

The exponential is the same algebra in both, which is not a coincidence: oref
took it from LoopKit. The difference that matters is the ten minute delay LoopKit
applies before any insulin acts at all, and which oref does not.

Everything here returns the FRACTION REMAINING at a time after the dose, so 1 at
delivery and 0 at the end of the action. Activity is its complement, and the
grid-level kernels below are what the action module convolves with.
"""
from __future__ import annotations

import numpy as np

from . import config

# LoopKit ExponentialInsulinModelPreset: (action duration, peak) in minutes.
# lyumjev and afrezza are omitted deliberately: both postdate this data.
LOOPKIT_PRESETS = {
    "loop_adult": (360.0, 75.0),      # rapidActingAdult, LoopKit's default
    "loop_child": (360.0, 65.0),      # rapidActingChild
    "loop_fiasp": (360.0, 55.0),      # fiasp
}
LOOPKIT_DELAY = 10.0

# oref/AAPS: the same exponential, no delay. The peaks are the ones the
# platform offers for rapid, ultra-rapid and Lyumjev respectively.
OREF_PRESETS = {
    "oref_6h75": (360.0, 75.0),
    "oref_6h65": (360.0, 65.0),
    "oref_6h55": (360.0, 55.0),
    "oref_5h75": (300.0, 75.0),
    "oref_7h75": (420.0, 75.0),
}

# LoopKit WalshInsulinModel, fourth order fits by action duration. Coefficients
# run highest power first and t is in minutes.
WALSH_COEFFS = {
    180.0: (-3.2030e-9, 1.354e-6, -1.759e-4, 9.255e-4, 0.99951),
    240.0: (-3.310e-10, 2.530e-7, -5.510e-5, -9.086e-4, 0.99950),
    300.0: (-2.950e-10, 2.320e-7, -5.550e-5, 4.490e-4, 0.99300),
    360.0: (-1.493e-10, 1.413e-7, -4.095e-5, 6.365e-4, 0.99700),
}


def remaining_exponential(t_min, dia: float, peak: float, delay: float = 0.0):
    """Fraction of a dose still to act, `t_min` after it was given.

    This is LoopKit's ExponentialInsulinModel, which oref also uses. The three
    precomputed terms are its tau, a and S.
    """
    t = np.asarray(t_min, dtype=float) - delay
    tau = peak * (1.0 - peak / dia) / (1.0 - 2.0 * peak / dia)
    a = 2.0 * tau / dia
    s = 1.0 / (1.0 - a + (1.0 + a) * np.exp(-dia / tau))
    with np.errstate(over="ignore", invalid="ignore"):
        r = 1.0 - s * (1.0 - a) * (
            ((t ** 2) / (tau * dia * (1.0 - a)) - t / tau - 1.0) * np.exp(-t / tau) + 1.0)
    r = np.where(t <= 0.0, 1.0, r)
    r = np.where(t >= dia, 0.0, r)
    return np.clip(r, 0.0, 1.0)


def remaining_walsh(t_min, dia: float, delay: float = LOOPKIT_DELAY):
    """Fraction remaining under LoopKit's Walsh model.

    LoopKit rounds the action duration to the nearest whole hour it has a fit
    for, and clamps outside three to six hours, so that is done here too.
    """
    dia = float(min(max(round(dia / 60.0) * 60.0, 180.0), 360.0))
    c4, c3, c2, c1, c0 = WALSH_COEFFS[dia]
    t = np.asarray(t_min, dtype=float) - delay
    r = ((c4 * t + c3) * t + c2) * t * t + c1 * t + c0
    r = np.where(t <= 0.0, 1.0, r)
    r = np.where(t >= dia, 0.0, r)
    return np.clip(r, 0.0, 1.0)


def remaining(model: str, t_min):
    """Fraction remaining under a named model from the registry."""
    if model in LOOPKIT_PRESETS:
        dia, peak = LOOPKIT_PRESETS[model]
        return remaining_exponential(t_min, dia, peak, delay=LOOPKIT_DELAY)
    if model in OREF_PRESETS:
        dia, peak = OREF_PRESETS[model]
        return remaining_exponential(t_min, dia, peak, delay=0.0)
    if model.startswith("walsh_"):
        return remaining_walsh(t_min, float(model.split("_")[1].rstrip("h")) * 60.0)
    raise KeyError(f"unknown insulin model: {model!r}")


def duration(model: str) -> float:
    """Total minutes from a dose until none of it is left, delay included."""
    if model in LOOPKIT_PRESETS:
        return LOOPKIT_PRESETS[model][0] + LOOPKIT_DELAY
    if model in OREF_PRESETS:
        return OREF_PRESETS[model][0]
    if model.startswith("walsh_"):
        dia = min(max(round(float(model.split("_")[1].rstrip("h")) * 60.0 / 60.0) * 60.0, 180.0), 360.0)
        return dia + LOOPKIT_DELAY
    raise KeyError(f"unknown insulin model: {model!r}")


def kernel(model: str, bin_min: float = config.GRID_MIN) -> np.ndarray:
    """Fraction remaining at each whole bin after a dose, for convolution.

    Element m is what is left `m` bins later, so element 0 is 1 and the last
    element is 0. Convolving a delivery series with this gives insulin on board.
    """
    n = int(np.ceil(duration(model) / bin_min)) + 1
    k = remaining(model, np.arange(n) * bin_min)
    # The Walsh fits are fourth order polynomials and wiggle very slightly near
    # the origin (0.0002 at fifteen minutes in the six hour fit). Insulin on
    # board cannot rise, so the kernel is made monotone. The exponentials are
    # already monotone and this leaves them untouched.
    return np.minimum.accumulate(k)


MODELS = tuple(LOOPKIT_PRESETS) + tuple(OREF_PRESETS) + ("walsh_3h", "walsh_4h", "walsh_5h", "walsh_6h")
