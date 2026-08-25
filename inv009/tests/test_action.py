"""The action split has to be right or nothing downstream means anything."""
from __future__ import annotations

import numpy as np
import pytest

from inv009 import insulin_models as M
from inv009.action import iob_series, window_action


@pytest.mark.parametrize("model", M.MODELS)
def test_kernel_shape(model):
    k = M.kernel(model)
    assert k[0] == pytest.approx(1.0, abs=1e-6)     # nothing has acted at delivery
    assert k[-1] == pytest.approx(0.0, abs=1e-6)    # all of it has acted by the end
    assert np.all(np.diff(k) <= 1e-12)              # insulin on board cannot rise


def test_loopkit_adult_matches_published_curve():
    # LoopKit's rapidActingAdult, 360 minutes with a 75 minute peak and a ten
    # minute delay: about a quarter of a dose is left at three hours.
    assert M.remaining("loop_adult", 180) == pytest.approx(0.241, abs=0.005)
    assert M.remaining("loop_adult", 0) == pytest.approx(1.0)
    assert M.remaining("loop_adult", 5) == pytest.approx(1.0)   # inside the delay
    assert M.remaining("loop_adult", 400) == pytest.approx(0.0)


def test_oref_has_no_delay_and_loop_does():
    assert M.remaining("oref_6h75", 5) < 1.0
    assert M.remaining("loop_adult", 5) == pytest.approx(1.0)


def test_single_dose_before_window_is_iob_decay():
    """With no delivery inside it, action is just the fall in insulin on board."""
    k = M.kernel("oref_6h75")
    u = np.zeros(200)
    u[10] = 3.0
    h = 48
    a_pre, a_in = window_action(u, k, h)
    iob = iob_series(u, k)
    i = 60
    assert a_in[i] == pytest.approx(0.0, abs=1e-12)
    assert a_pre[i] == pytest.approx(iob[i] - iob[i + h], abs=1e-9)


def test_dose_inside_window_is_counted_as_in_not_pre():
    k = M.kernel("oref_6h75")
    u = np.zeros(200)
    u[70] = 2.0
    h = 48
    a_pre, a_in = window_action(u, k, h)
    i = 60
    assert a_pre[i] == pytest.approx(0.0, abs=1e-12)
    # Delivered at bin 70, so 38 bins of a 48 bin window remain for it to act in.
    assert a_in[i] == pytest.approx(2.0 * (1.0 - M.remaining("oref_6h75", 38 * 5)), abs=1e-9)


def test_total_action_never_exceeds_insulin_present():
    k = M.kernel("oref_6h75")
    rng = np.random.default_rng(0)
    u = rng.gamma(2.0, 0.05, size=500)
    h = 48
    a_pre, a_in = window_action(u, k, h)
    assert np.all(a_pre[:-h] >= -1e-9)
    assert np.all(a_in[:-h] >= -1e-9)
    iob = iob_series(u, k)
    csum = np.concatenate([[0.0], np.cumsum(u)])
    for i in (100, 200, 300):
        delivered = csum[i + h + 1] - csum[i]
        assert a_pre[i] + a_in[i] <= iob[i] + delivered + 1e-9


def test_action_recovers_a_known_sensitivity():
    """Simulate glucose falling at a known mg/dL per unit and read it back."""
    k = M.kernel("oref_6h75")
    rng = np.random.default_rng(1)
    n, h, isf = 4000, 48, 45.0
    u = np.zeros(n)
    u[::12] = 0.1                                  # a steady basal-like drip
    u[rng.choice(n, 60, replace=False)] += 2.0     # scattered boluses
    a_pre, a_in = window_action(u, k, h)
    starts = np.arange(200, n - h - 1, 7)
    drop = isf * (a_pre + a_in)[starts]
    fit = np.linalg.lstsq(np.c_[np.ones(len(starts)), (a_pre + a_in)[starts]],
                          drop, rcond=None)[0]
    assert fit[1] == pytest.approx(isf, rel=1e-6)


def test_design_matrix_indices_are_named_not_positional():
    """The interaction must not displace the coefficient it is divided by.

    Reading the insulin coefficient from a fixed position, after inserting an
    interaction column ahead of it, returned the starting-glucose coefficient
    instead and flipped the sign of the reported exponent.
    """
    import pandas as pd
    from inv009.glucose_axis import _design

    d = pd.DataFrame({"a_pre": [1.0, 2.0, 3.0, 4.0], "bg0": [100.0, 150.0, 120.0, 180.0],
                      "pre_slope": [0.0, 0.1, -0.1, 0.2], "bg_m60": [100.0] * 4,
                      "bg_m120": [100.0] * 4, "hour": [0, 0, 1, 1]})
    X_plain, idx_plain = _design(d)
    X_int, idx_int = _design(d, extra=d.a_pre.to_numpy() * 0.5)

    assert np.allclose(X_plain[:, idx_plain["a_pre"]], d.a_pre)
    assert np.allclose(X_int[:, idx_int["a_pre"]], d.a_pre)
    assert np.allclose(X_int[:, idx_int["interaction"]], d.a_pre * 0.5)
    # The point of the map: adding the interaction moves bg0, and a positional
    # read of index 3 would silently pick it up as the insulin term.
    assert idx_int["bg0"] != idx_plain["bg0"]
    assert np.allclose(X_int[:, idx_int["bg0"]], d.bg0 - 100.0)
