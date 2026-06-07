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


# ===========================================================================
# v-next:  ISF(BG) = (K_user / √TDD) · g(BG)
#
#   level:    K_user / √TDD      — universal −½ exponent, per-user constant K_user
#   glucose:  g(BG)              — normalised to 1.0 at normal target, so it
#                                  composes with the level without disturbing the
#                                  per-user anchor
#
# Two glucose curves are provided:
#   • g_quartic  — the Diabeloop clinical-population ISF shape (the default; its
#                  exponent comes from controlled clinical data, not device logs)
#   • g_powerlaw — (target/BG)^k, a one-parameter alternative of the same family
#
# Both equal 1.0 at the normal target. BG is high-capped at config.BG_CAP
# (excess/3, as v1/v2) and low-floored at BG_FLOOR_VNEXT to stay inside the range
# the clinical curve was fit over (no extrapolation past the physiological band).
# ===========================================================================

BG_FLOOR_VNEXT = 54.0      # ~3.0 mmol/L; low end of the clinical fit range

# Diabeloop population ISF–glucose quartic (mg/dL per U vs glucose mg/dL).
_QUARTIC = (272.0, -3.121, 0.01511, -3.305e-5, 2.69e-8)


def quartic_isf(bg):
    """Raw Diabeloop population quartic (un-normalised), vectorised."""
    bg = np.asarray(bg, dtype=float)
    c0, c1, c2, c3, c4 = _QUARTIC
    return c0 + c1 * bg + c2 * bg**2 + c3 * bg**3 + c4 * bg**4


def _bg_clamped_vnext(bg, cap: float = config.BG_CAP, floor: float = BG_FLOOR_VNEXT):
    """High-cap (excess/3 above cap) then low-floor — keeps BG in the clinical range."""
    return np.maximum(cap_bg(bg, cap), floor)


def g_quartic(bg, normal_target: float = config.NORMAL_TARGET,
              cap: float = config.BG_CAP, floor: float = BG_FLOOR_VNEXT):
    """Diabeloop quartic glucose curve, normalised to 1.0 at the normal target."""
    bg_c = _bg_clamped_vnext(bg, cap, floor)
    return quartic_isf(bg_c) / quartic_isf(normal_target)


def g_powerlaw(bg, k: float = 1.3, normal_target: float = config.NORMAL_TARGET,
               cap: float = config.BG_CAP, floor: float = BG_FLOOR_VNEXT):
    """Power-law glucose curve (target/BG)^k, normalised to 1.0 at the normal target.
    Default k≈1.3 matches the Diabeloop quartic's slope over BG 70–250."""
    bg_c = _bg_clamped_vnext(bg, cap, floor)
    return (normal_target / bg_c) ** k


def isf_vnext(bg, tdd, k_user, curve: str = "quartic", k: float = 1.3,
              normal_target: float = config.NORMAL_TARGET,
              cap: float = config.BG_CAP, floor: float = BG_FLOOR_VNEXT):
    """v-next ISF: (K_user/√TDD) · g(BG).

    k_user : per-user constant K_user (scalar or array). Tier-1 anchor is
             profile_ISF · √(median TDD); Tier-2 swaps measured_ISF for profile_ISF.
    curve  : 'quartic' (Diabeloop default) or 'powerlaw'.
    """
    tdd = np.asarray(tdd, dtype=float)
    g = g_quartic(bg, normal_target, cap, floor) if curve == "quartic" \
        else g_powerlaw(bg, k, normal_target, cap, floor)
    with np.errstate(divide="ignore", invalid="ignore"):
        level = np.asarray(k_user, dtype=float) / np.sqrt(tdd)
    isf = level * g
    return np.where(tdd > 0, isf, np.nan)


def k_user_tier1(profile_isf, median_tdd):
    """Tier-1 (profile-anchored) per-user constant: K = profile_ISF · √(median TDD).
    Makes isf_vnext return profile_ISF at the user's median TDD and normal target."""
    return np.asarray(profile_isf, dtype=float) * np.sqrt(np.asarray(median_tdd, dtype=float))
