"""Pure functions for the INV-011 implementation-mismatch analysis.

The comparison is between two software configurations. The reference uses the
fixed parameters in the AAPS Lyumjev plugin, a 45 minute oref exponential peak
and divisor 75. The implementation under test uses a 70 minute peak and divisor
55. Neither configuration is asserted to be an individual's pharmacology.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from inv008.dynisf import isf_v1, isf_v2


@dataclass(frozen=True)
class Scenario:
    name: str
    peak_min: float
    divisor: float


SCENARIOS = (
    Scenario("matched_45_75", 45.0, 75.0),
    Scenario("divisor_only_45_55", 45.0, 55.0),
    Scenario("peak_only_70_75", 70.0, 75.0),
    Scenario("combined_70_55", 70.0, 55.0),
)


def exponential_terms(dia_min: float, peak_min: float) -> tuple[float, float, float]:
    """Return the oref exponential model's tau, a and scale terms."""
    if not 0.0 < peak_min < dia_min / 2.0:
        raise ValueError("peak must be greater than zero and less than half the DIA")
    tau = peak_min * (1.0 - peak_min / dia_min) / (1.0 - 2.0 * peak_min / dia_min)
    a = 2.0 * tau / dia_min
    scale = 1.0 / (1.0 - a + (1.0 + a) * np.exp(-dia_min / tau))
    return float(tau), float(a), float(scale)


def remaining_exponential(t_min, dia_min: float, peak_min: float):
    """Fraction of a delivery still to act under the oref exponential model."""
    t = np.asarray(t_min, dtype=float)
    tau, a, scale = exponential_terms(dia_min, peak_min)
    with np.errstate(over="ignore", invalid="ignore"):
        remaining = 1.0 - scale * (1.0 - a) * (
            ((t**2) / (tau * dia_min * (1.0 - a)) - t / tau - 1.0)
            * np.exp(-t / tau)
            + 1.0
        )
    remaining = np.where(t <= 0.0, 1.0, remaining)
    remaining = np.where(t >= dia_min, 0.0, remaining)
    return np.clip(remaining, 0.0, 1.0)


def activity_exponential(t_min, dia_min: float, peak_min: float):
    """Fraction of a one-unit delivery acting per minute."""
    t = np.asarray(t_min, dtype=float)
    tau, _, scale = exponential_terms(dia_min, peak_min)
    activity = (scale / tau**2) * t * (1.0 - t / dia_min) * np.exp(-t / tau)
    activity = np.where((t <= 0.0) | (t >= dia_min), 0.0, activity)
    return np.maximum(activity, 0.0)


def dynamic_isf(equation: str, bg, tdd, divisor: float):
    """Evaluate V1 or the author-confirmed no-plus-one V2 equation."""
    if equation == "v1":
        return isf_v1(bg, tdd, divisor=divisor)
    if equation == "v2":
        return isf_v2(bg, tdd, divisor=divisor)
    raise ValueError(f"unknown equation: {equation!r}")


def first_crossing_minutes(time_min, values, level: float = 1.0) -> float | None:
    """First upward crossing of level, linearly interpolated."""
    t = np.asarray(time_min, dtype=float)
    y = np.asarray(values, dtype=float)
    finite = np.isfinite(y)
    for i in range(1, len(y)):
        if finite[i - 1] and finite[i] and y[i - 1] < level <= y[i]:
            if y[i] == y[i - 1]:
                return float(t[i])
            frac = (level - y[i - 1]) / (y[i] - y[i - 1])
            return float(t[i - 1] + frac * (t[i] - t[i - 1]))
    return None

