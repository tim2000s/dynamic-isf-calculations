from inv008.divisor_sensitivity import equation_effects


def test_divisor_55_makes_both_equations_more_aggressive_than_75():
    for row in equation_effects():
        assert row["v1_correction_ratio_d55_vs_d75"] > 1.0
        assert row["v2_correction_ratio_d55_vs_d75"] > 1.0


def test_v2_divisor_effect_is_large_near_target():
    row = next(r for r in equation_effects() if r["bg"] == 99.0)
    assert row["v2_correction_ratio_d55_vs_d75"] > 2.0
    assert row["v1_correction_ratio_d55_vs_d75"] < 1.3
