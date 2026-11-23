"""inv020: the per-decision algebra, the percentage switches and the crossover."""

import numpy as np
import pandas as pd

from inv020.replay import changes, cross_glucose, one_per_timestamp

A = np.array


def test_logged_sensitivity_returns_the_logged_prediction_and_request():
    target, ev, mp = A([100.0, 100.0]), A([150.0, 140.0]), A([160.0, 120.0])
    s, prof, iob, act = (
        A([40.0, 55.0]),
        A([45.0, 50.0]),
        A([2.0, -0.5]),
        A([0.02, -0.004]),
    )
    out = changes(target, ev, mp, s, prof, iob, act, isf_log=s, isf_curve=s * 1.5)
    assert np.allclose(out["p_log"], np.fmin(ev, mp))
    assert np.allclose(out["r_log"], (np.fmin(ev, mp) - target) / s)
    # positive IOB above 30 x activity: a higher ISF predicts a lower eventual glucose
    assert out["p_curve"][0] < out["p_log"][0]


def test_negative_iob_is_capped_at_the_profile_and_low_branch_doubles():
    # IOB -1: the IOB term uses min(ISF, profile) = 50 for both arms, so only the deviation term moves
    out = changes(
        A([100.0]),
        A([90.0]),
        A([95.0]),
        A([40.0]),
        A([50.0]),
        A([-1.0]),
        A([0.0]),
        isf_log=A([40.0]),
        isf_curve=A([80.0]),
    )
    assert np.allclose(out["p_curve"], 90.0 + 10.0)  # -(-1) x (50 - 40)
    # eventual glucose below target: oref's 2 x min(0, (eventual - target) / ISF), on eventual only
    assert np.allclose(out["r_log"], 2 * (90.0 - 100.0) / 40.0)
    assert np.allclose(out["r_curve"], 0.0)


def test_crossover_and_duplicates():
    bg = np.repeat([100.0, 120.0, 140.0, 160.0], 25)
    d = pd.DataFrame(
        {"bg": bg, "isf_log": 50.0, "isf_curve": np.where(bg >= 140, 40.0, 60.0)}
    )
    assert cross_glucose(d) == 140.0
    d.loc[d.bg == 160, "isf_curve"] = 60.0  # the top band reverts, so no band qualifies
    assert np.isnan(cross_glucose(d))
    r = pd.DataFrame({"ts": [1, 1, 2, 2, 3], "x": [5, 5, 6, 7, 8]})
    assert one_per_timestamp(r, "ts").ts.tolist() == [1, 3]


def test_logged_prediction_is_recovered_with_the_isf_it_was_made_with():
    # the log arm's ISF (50) differs from the one the logged prediction used (40): its prediction moves
    target, ev, mp = A([100.0]), A([150.0]), A([150.0])
    iob, act = A([2.0]), A([0.0])
    out = changes(
        target,
        ev,
        mp,
        A([40.0]),
        A([45.0]),
        iob,
        act,
        isf_log=A([50.0]),
        isf_curve=A([40.0]),
    )
    assert np.allclose(
        out["p_curve"], 150.0
    )  # the arm with the logged ISF returns the logged prediction
    assert np.allclose(out["p_log"], 150.0 - 2.0 * (50.0 - 40.0))
