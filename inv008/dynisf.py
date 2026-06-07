"""Dynamic ISF equations: v1 (Chris Wilson, original) and v2 (revised maths).

v1:
    sensNormalTarget = 1800 / (TDD * ln(NT/divisor + 1))            # ISF proportional to 1/TDD
    ISF              = sensNT * scaler

v2:
    sensNormalTarget = 2300 / (ln(NT/divisor + 1) * TDD^2 * 0.02)   # ISF proportional to 1/TDD^2
    ISF              = sensNT * scaler

(The optional velocity damping of v1's glucose response is held at its default of 1.0,
i.e. the full scaler, throughout — so both equations share the identical glucose term.)

Shared:
    bgAdj  = cap + (bg - cap)/3            if bg > cap
    scaler = ln(NT/divisor + 1) / ln(bgAdj/divisor + 1)
    TDD blend: W8H = (1.4*TDD_4h + 0.6*TDD_8to4h) * 3
        if W8H < 0.75*TDD_7d:  adj7 = W8H + (W8H/TDD_7d)*(TDD_7d - W8H)
                               TDD  = 0.34*adj7 + 0.33*TDD_1d + 0.33*W8H
        else:                  TDD  = 0.33*W8H + 0.34*TDD_7d + 0.33*TDD_1d

All functions are vectorised: bg / tdd components may be numpy arrays or scalars.
"""
from __future__ import annotations

import numpy as np

from inv008 import config


def cap_bg(bg, cap: float = config.BG_CAP):
    """BG compression above the cap: excess contributes at 1/3 weight."""
    bg = np.asarray(bg, dtype=float)
    return np.where(bg > cap, cap + (bg - cap) / 3.0, bg)


def scaler(bg, normal_target: float = config.NORMAL_TARGET,
           divisor: float = config.INSULIN_DIVISOR, cap: float = config.BG_CAP):
    """ln(NT/div + 1) / ln(bgAdj/div + 1); 1.0 where the denominator is invalid."""
    bg_adj = cap_bg(bg, cap)
    sbg = np.log(bg_adj / divisor + 1.0)
    log_term = np.log(normal_target / divisor + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(sbg > 0, log_term / sbg, 1.0)
    return s


def blend_tdd(tdd_4h, tdd_8to4h, tdd_1d, tdd_7d,
              adjust_factor: float = config.ADJUST_FACTOR):
    """Master's TddCalculator blend. Returns NaN where any component is missing
    or tdd_7d <= 0 (the plugins fall back to profile ISF in that case)."""
    t4 = np.asarray(tdd_4h, dtype=float)
    t84 = np.asarray(tdd_8to4h, dtype=float)
    t1 = np.asarray(tdd_1d, dtype=float)
    t7 = np.asarray(tdd_7d, dtype=float)

    w8h = (1.4 * t4 + 0.6 * t84) * 3.0
    adj7 = w8h + (w8h / np.where(t7 > 0, t7, np.nan)) * (t7 - w8h)
    low = w8h < 0.75 * t7
    tdd = np.where(low,
                   0.34 * adj7 + 0.33 * t1 + 0.33 * w8h,
                   0.33 * w8h + 0.34 * t7 + 0.33 * t1)
    tdd = tdd * adjust_factor
    valid = np.isfinite(t4) & np.isfinite(t84) & np.isfinite(t1) & np.isfinite(t7) & (t7 > 0)
    return np.where(valid, tdd, np.nan)


def sens_normal_target_v1(tdd, normal_target: float = config.NORMAL_TARGET,
                          divisor: float = config.INSULIN_DIVISOR):
    """v1: 1800 rule — anchor at normal target, ISF proportional to 1/TDD."""
    tdd = np.asarray(tdd, dtype=float)
    log_term = np.log(normal_target / divisor + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = 1800.0 / (tdd * log_term)
    return np.where((tdd > 0) & (log_term > 0), s, np.nan)


def sens_normal_target_v2(tdd, normal_target: float = config.NORMAL_TARGET,
                          divisor: float = config.INSULIN_DIVISOR):
    """v2: 2300 / (logTerm * TDD^2 * 0.02) — anchor proportional to 1/TDD^2."""
    tdd = np.asarray(tdd, dtype=float)
    log_term = np.log(normal_target / divisor + 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = 2300.0 / (log_term * tdd * tdd * 0.02)
    return np.where((tdd > 0) & (log_term > 0), s, np.nan)


def isf_v1(bg, tdd, normal_target: float = config.NORMAL_TARGET,
           divisor: float = config.INSULIN_DIVISOR, velocity: float = config.VELOCITY,
           cap: float = config.BG_CAP):
    """V1 variableSens at the given BG and (blended) TDD."""
    sens_nt = sens_normal_target_v1(tdd, normal_target, divisor)
    s = scaler(bg, normal_target, divisor, cap)
    return sens_nt * (1.0 - (1.0 - s) * velocity)


def isf_v2(bg, tdd, normal_target: float = config.NORMAL_TARGET,
           divisor: float = config.INSULIN_DIVISOR, cap: float = config.BG_CAP):
    """V2 variableSens at the given BG and (blended) TDD."""
    sens_nt = sens_normal_target_v2(tdd, normal_target, divisor)
    s = scaler(bg, normal_target, divisor, cap)
    return sens_nt * s


def v2_over_v1_ratio(tdd):
    """Closed form: ISF_v2/ISF_v1 = 2300/(0.02*1800*TDD) = 63.888…/TDD
    (exact at velocity=1.0, any BG)."""
    tdd = np.asarray(tdd, dtype=float)
    return 2300.0 / (0.02 * 1800.0 * tdd)


# ---------------------------------------------------------------------------
# v2 (updated) — anchor and glucose term drop the "+1"; BG floored at divisor+1
#   sensNormalTarget = 2300 / (ln(target/divisor) · TDD² · 0.02)
#   scaler           = ln(target/divisor) / ln(BG_floored/divisor)
#   ISF(BG)          = 115000 / (TDD² · ln(BG_floored/divisor))
# BG_floored = clamp(BG): high cap at 210 (excess/3), then floor at divisor+1
# (so ln(BG/divisor) stays positive and finite).
# ---------------------------------------------------------------------------

def _bg_floored_v2u(bg, divisor: float = config.INSULIN_DIVISOR, cap: float = config.BG_CAP):
    bg_adj = cap_bg(bg, cap)                       # existing high cap
    return np.maximum(bg_adj, divisor + 1.0)       # NEW low floor at divisor+1


def sens_normal_target_v2_updated(tdd, normal_target: float = config.NORMAL_TARGET,
                                  divisor: float = config.INSULIN_DIVISOR):
    """v2 updated anchor: 2300 / (ln(target/divisor) · TDD² · 0.02) — no +1."""
    tdd = np.asarray(tdd, dtype=float)
    log_term = np.log(normal_target / divisor)     # no +1
    with np.errstate(divide="ignore", invalid="ignore"):
        s = 2300.0 / (log_term * tdd * tdd * 0.02)
    return np.where((tdd > 0) & (log_term > 0), s, np.nan)


def isf_v2_updated(bg, tdd, normal_target: float = config.NORMAL_TARGET,
                   divisor: float = config.INSULIN_DIVISOR, cap: float = config.BG_CAP):
    """Updated v2 ISF: collapses to 115000 / (TDD² · ln(BG_floored/divisor))."""
    tdd = np.asarray(tdd, dtype=float)
    bg_f = _bg_floored_v2u(bg, divisor, cap)
    log_bg = np.log(bg_f / divisor)                # no +1; BG floored so > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        isf = 2300.0 / (0.02 * tdd * tdd * log_bg)
    return np.where((tdd > 0) & (log_bg > 0), isf, np.nan)


def v2updated_over_v1_ratio(bg, tdd, normal_target: float = config.NORMAL_TARGET,
                            divisor: float = config.INSULIN_DIVISOR, cap: float = config.BG_CAP):
    """ISF_v2updated / ISF_v1 = (63.888…/TDD) · ln(BG/divisor+1)/ln(BG_floored/divisor).
    Now BG-dependent (the glucose terms no longer cancel)."""
    return isf_v2_updated(bg, tdd, normal_target, divisor, cap) / isf_v1(
        bg, tdd, normal_target, divisor, cap=cap)
