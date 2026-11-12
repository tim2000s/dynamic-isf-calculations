"""The poster digitiser: medians agree with the equations the authors fitted to them."""
import numpy as np

from inv013.figures import OUTCOME_ARMS, digitise_poster, harm_rows, link_rows, outcome_rows


def test_digitised_medians_match_the_published_equations():
    d = digitise_poster()
    assert list(d.glucose) == list(range(60, 381, 20))
    assert np.abs(d["median"] - d.published_fit).max() < 3.0
    at = d.set_index("glucose")
    assert abs(at.loc[100, "median"] - 79) < 1.5 and abs(at.loc[180, "median"] - 35) < 2.5


def test_quartiles_bracket_the_median_and_only_the_off_axis_ones_are_blank():
    d = digitise_poster().set_index("glucose")
    assert list(d.index[d.q25.isna()]) == [60, 80] and list(d.index[d.q75.isna()]) == [60, 80, 100]
    seen = d.dropna()
    assert (seen.q25 < seen["median"]).all() and (seen["median"] < seen.q75).all()
    assert seen.q25.is_monotonic_decreasing and seen.q75.is_monotonic_decreasing


def test_link_rows_reads_the_measured_coefficients():
    built, (dose, profile) = link_rows()
    assert len(built) == 8 and dose[1] < dose[0] < dose[2] and profile[1] < profile[0] < profile[2]


def test_outcome_rows_cover_every_arm_drawn_and_each_interval_holds_its_mean():
    rows = outcome_rows()
    assert {arm for arm, _, _ in OUTCOME_ARMS} == set(rows)
    assert all(low <= mean <= high for cells in rows.values() for mean, low, high in cells.values())
    assert rows["logarithmic, Trio"]["tbr54"][0] < 0 < rows["candidate"]["tbr54"][0]      # a minus sign survives the parse
    assert len(harm_rows()) == 8
