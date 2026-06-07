"""Hand-computed fixtures for the V1/V2 ISF formulas (master @ 2c3e3276)."""
import math

import numpy as np
import pytest

from inv008.dynisf import (blend_tdd, cap_bg, isf_v1, isf_v2, scaler,
                              sens_normal_target_v1, sens_normal_target_v2,
                              v2_over_v1_ratio)

NT, DIV = 99.0, 75
LOG_TERM = math.log(NT / DIV + 1.0)  # ln(2.32) = 0.841567...


def test_sens_normal_target_v1():
    # 1800 / (50 * ln(2.32)) = 42.78
    assert sens_normal_target_v1(50.0) == pytest.approx(1800.0 / (50.0 * LOG_TERM))
    assert sens_normal_target_v1(50.0) == pytest.approx(42.78, abs=0.01)


def test_sens_normal_target_v2():
    # 2300 / (ln(2.32) * 2500 * 0.02) = 54.66
    assert sens_normal_target_v2(50.0) == pytest.approx(
        2300.0 / (LOG_TERM * 50.0 ** 2 * 0.02))
    assert sens_normal_target_v2(50.0) == pytest.approx(54.66, abs=0.01)


def test_closed_form_ratio():
    # V2/V1 = 63.888../TDD at any BG (velocity = 1)
    for tdd in (20.0, 63.888888, 100.0):
        bgs = np.array([80.0, 120.0, 250.0])
        ratio = isf_v2(bgs, tdd) / isf_v1(bgs, tdd)
        assert ratio == pytest.approx(v2_over_v1_ratio(tdd), rel=1e-9)
    # crossover: identical ISF at TDD = 2300/(0.02*1800)
    assert v2_over_v1_ratio(2300.0 / 36.0) == pytest.approx(1.0)


def test_at_normal_target_scaler_is_one():
    assert scaler(NT) == pytest.approx(1.0)
    assert isf_v1(NT, 40.0) == pytest.approx(sens_normal_target_v1(40.0))
    assert isf_v2(NT, 40.0) == pytest.approx(sens_normal_target_v2(40.0))


def test_bg_cap_compression():
    assert cap_bg(200.0) == 200.0
    assert cap_bg(210.0) == 210.0
    # 270 → 210 + 60/3 = 230
    assert cap_bg(270.0) == pytest.approx(230.0)
    # capped ISF: ISF(270) must equal ISF computed at bgAdj=230 without cap
    assert isf_v1(270.0, 50.0) == pytest.approx(
        float(sens_normal_target_v1(50.0)) * (math.log(NT / DIV + 1) / math.log(230.0 / DIV + 1)))


def test_velocity_dampening_v1():
    # velocity = 0 → variableSens == sensNormalTarget regardless of BG
    assert isf_v1(180.0, 50.0, velocity=0.0) == pytest.approx(
        float(sens_normal_target_v1(50.0)))
    # velocity = 0.5 → halfway between sensNT and full-scaler value
    full = isf_v1(180.0, 50.0, velocity=1.0)
    nt = float(sens_normal_target_v1(50.0))
    assert isf_v1(180.0, 50.0, velocity=0.5) == pytest.approx((full + nt) / 2)


def test_blend_tdd_low_w8h_branch():
    # t4=2, t84=3 → W8H = (1.4*2 + 0.6*3)*3 = 13.8 < 0.75*40
    # adj7 = 13.8 + (13.8/40)*(40-13.8) = 22.839
    # tdd  = 0.34*22.839 + 0.33*38 + 0.33*13.8 = 24.859
    assert blend_tdd(2.0, 3.0, 38.0, 40.0) == pytest.approx(24.859, abs=0.001)


def test_blend_tdd_normal_branch():
    # t4=8, t84=7 → W8H = (11.2+4.2)*3 = 46.2 ≥ 0.75*40
    # tdd = 0.33*46.2 + 0.34*40 + 0.33*38 = 41.386
    assert blend_tdd(8.0, 7.0, 38.0, 40.0) == pytest.approx(41.386, abs=0.001)


def test_blend_tdd_gates():
    assert np.isnan(blend_tdd(np.nan, 3.0, 38.0, 40.0))
    assert np.isnan(blend_tdd(2.0, 3.0, 38.0, 0.0))    # tdd_7d must be > 0
    assert np.isnan(blend_tdd(2.0, 3.0, np.nan, 40.0))


def test_vectorised():
    bg = np.linspace(70, 300, 47)
    out = isf_v2(bg, 45.0)
    assert out.shape == bg.shape
    assert np.all(np.diff(out) <= 1e-12)  # ISF monotonically falls as BG rises


def test_v2_updated_collapse_and_floor():
    import math
    from inv008.dynisf import isf_v2_updated, sens_normal_target_v2_updated
    T, div = 99.0, 75
    # ISF at target == anchor (scaler 1)
    assert isf_v2_updated(T, 50.0) == pytest.approx(float(sens_normal_target_v2_updated(50.0)))
    # collapse: 115000 / (TDD^2 * ln(BG/div)), no +1
    assert float(isf_v2_updated(140.0, 50.0)) == pytest.approx(115000.0 / (2500 * math.log(140.0 / div)))
    # BG floored at divisor+1: BG<=76 clamps to 76
    assert float(isf_v2_updated(70.0, 50.0)) == pytest.approx(float(isf_v2_updated(76.0, 50.0)))
    # anchor dropped the +1 → ~3x the old v2 anchor
    from inv008.dynisf import sens_normal_target_v2
    assert float(sens_normal_target_v2_updated(50.0)) > 2.5 * float(sens_normal_target_v2(50.0))


def test_v2_updated_ratio_is_bg_dependent():
    from inv008.dynisf import v2updated_over_v1_ratio
    # the glucose terms no longer cancel → ratio varies with BG (unlike old 63.9/TDD)
    r_low = float(v2updated_over_v1_ratio(90.0, 50.0))
    r_high = float(v2updated_over_v1_ratio(250.0, 50.0))
    assert r_low > r_high > 1.0


# --- v-next: (K_user/√TDD) · g(BG) -----------------------------------------

def test_g_curves_unity_at_target():
    from inv008.dynisf import g_quartic, g_powerlaw
    # both glucose curves are normalised to 1.0 at the normal target
    assert float(g_quartic(NT)) == pytest.approx(1.0)
    assert float(g_powerlaw(NT)) == pytest.approx(1.0)


def test_g_curves_fall_with_glucose():
    from inv008.dynisf import g_quartic, g_powerlaw
    bg = np.linspace(60, 300, 49)
    for g in (np.asarray(g_quartic(bg)), np.asarray(g_powerlaw(bg))):
        assert np.all(np.diff(g) < 0)          # ISF falls as BG rises
    # hypo-protective below target, more aggressive above
    assert float(g_quartic(70.0)) > 1.0 > float(g_quartic(160.0))


def test_g_quartic_matches_diabeloop_coeffs():
    from inv008.dynisf import g_quartic, quartic_isf
    def q(g): return 272 - 3.121*g + 0.01511*g**2 - 3.305e-5*g**3 + 2.69e-8*g**4
    assert float(quartic_isf(140.0)) == pytest.approx(q(140.0))
    assert float(g_quartic(140.0)) == pytest.approx(q(140.0) / q(NT))


def test_g_bg_clamps():
    from inv008.dynisf import g_quartic, BG_FLOOR_VNEXT
    # low floor: BG below the clinical floor clamps to the floor value
    assert float(g_quartic(40.0)) == pytest.approx(float(g_quartic(BG_FLOOR_VNEXT)))
    # high cap (excess/3 above 210): 270 → 230, evaluated as the raw quartic at 230
    from inv008.dynisf import quartic_isf
    assert float(g_quartic(270.0)) == pytest.approx(quartic_isf(230.0) / quartic_isf(NT))


def test_isf_vnext_anchor_and_tdd_law():
    from inv008.dynisf import isf_vnext, k_user_tier1
    profile_isf, med_tdd = 40.0, 36.0
    K = k_user_tier1(profile_isf, med_tdd)
    # at the user's median TDD and normal target, Tier-1 returns profile ISF exactly
    assert float(isf_vnext(NT, med_tdd, K)) == pytest.approx(profile_isf)
    # √TDD law: doubling TDD multiplies the at-target ISF by 1/√2
    assert float(isf_vnext(NT, 2 * med_tdd, K)) == pytest.approx(profile_isf / math.sqrt(2))
    # glucose curve composes multiplicatively
    from inv008.dynisf import g_quartic
    assert float(isf_vnext(160.0, med_tdd, K)) == pytest.approx(profile_isf * float(g_quartic(160.0)))


def test_isf_vnext_powerlaw_curve_option():
    from inv008.dynisf import isf_vnext, g_powerlaw, k_user_tier1
    K = k_user_tier1(40.0, 36.0)
    assert float(isf_vnext(160.0, 36.0, K, curve="powerlaw", k=1.3)) == pytest.approx(
        40.0 * float(g_powerlaw(160.0, k=1.3)))
