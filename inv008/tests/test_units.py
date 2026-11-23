"""ISF units in the head-to-head replay."""

import numpy as np

from inv008.config import MMOL_TO_MGDL
from inv008.head_to_head import sug_isf_mgdl


def test_trio_converted_by_reading_oref0_left_alone():
    v = np.array([5.0, 90.0, 14.0])
    assert np.allclose(
        sug_isf_mgdl(v, "oref_v5"), [5.0 * MMOL_TO_MGDL, 90.0, 14.0 * MMOL_TO_MGDL]
    )
    assert np.allclose(sug_isf_mgdl(v, "oref_v7"), v)
