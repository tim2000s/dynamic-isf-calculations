"""Shared constants for INV-009.

Window choices are inherited from INV-008 where there is no reason to differ, so
the two bodies of work can be read against each other. Where they do differ it is
because this data has no loop predictions in it: everything is reconstructed from
delivered insulin, which changes what can be screened and what has to be modelled.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("DYNISF_ROOT", Path(__file__).resolve().parent.parent))
CACHE = ROOT / "inv009_cache"
GRID_CACHE = CACHE / "grid"
WINDOW_CACHE = CACHE / "windows"
RESULTS = ROOT / "results"
CHARTS = ROOT / "charts" / "inv009"

DSN = os.environ.get("INV009_DSN", "dbname=oref")
WORKERS = int(os.environ.get("INV009_WORKERS", "7"))

GRID_MIN = 5.0                  # minutes per bin, which is what CGM reports on
MMOL_TO_MGDL = 18.018

# ---------------------------------------------------------------- windows
# Overnight starts. Fasting is the only regime where a glucose move can be
# attributed to insulin without also modelling a meal, and these hours are the
# ones people are least likely to be eating in.
START_HOURS = (23, 0, 1, 2, 3)
HORIZON_MIN = 240               # 4 h, matching INV-008
ENDPOINT_TOL_MIN = 5.0          # how far an endpoint reading may sit from the mark
ENDPOINT_MEDIAN_MIN = 15.0      # endpoints are a median over this much glucose
MIN_CGM_FRACTION = 0.8          # of the horizon's bins that must carry a reading

BG0_MIN, BG0_MAX = 90.0, 300.0  # below 90 the body is defending, above 300 is rare and odd
CARB_FREE_H = 6.0               # required clearance from the last logged carbohydrate
MEAL_FREE_H = 4.0               # same idea where carbohydrate is not logged
# A meal-sized bolus, as a fraction of the person's own daily dose. Scale free on
# purpose: a unit is a meal to a three year old and a rounding error to an adult,
# so a fixed threshold in units would screen the paediatric cohorts to nothing.
MEAL_BOLUS_FRAC_TDD = 0.05
# Both of these condition on the outcome, so neither is in the primary screen.
# They are the sensitivity analysis: excluding windows where glucose rose drops
# exactly the windows insulin worked least well in, which flatters sensitivity,
# and excluding windows that went low drops the ones it worked best in.
RISE_MAX_30 = 20.0              # mg/dL in 30 min
BG_FLOOR_IN_WINDOW = 70.0

MIN_WINDOWS = 40                # per subject, to enter per-subject statistics
MIN_DAYS_TDD = 30               # per subject, to have a trustworthy daily dose
MIN_A_SD = 0.3                  # within-subject spread of acting insulin, in units

# Isolated correction boluses, for the difference-in-differences estimator.
EVENT_MIN_U = 0.3               # smaller than this moves nothing measurable
EVENT_QUIET_H = 3.0             # no other bolus either side
EVENT_BG_MATCH = 15.0           # mg/dL, how close a control night must start

# ---------------------------------------------------------------- cohorts
# closed_loop marks where the controller reacts to glucose inside the window,
# which is the confound the open-loop arm exists to escape. carbs marks where a
# meal can be screened out rather than inferred.
COHORTS = {
    "Loop":      dict(closed_loop=True,  carbs=True,  model="infer",    label="Loop (DIY)"),
    "ReplaceBG": dict(closed_loop=False, carbs=True,  model="oref_6h75", label="REPLACE-BG (open loop)"),
    "DCLP3":     dict(closed_loop=True,  carbs=False, model="oref_6h75", label="DCLP3 (Control-IQ)"),
    "DCLP5":     dict(closed_loop=True,  carbs=False, model="oref_6h75", label="DCLP5 (Control-IQ, 6-13y)"),
    "PEDAP":     dict(closed_loop=True,  carbs=False, model="oref_6h75", label="PEDAP (Control-IQ, 2-5y)"),
    "IOBP2":     dict(closed_loop=True,  carbs=False, model="oref_6h75", label="IOBP2 (bionic pancreas)"),
}
PRIMARY = ("Loop", "ReplaceBG")     # carbohydrate is logged, so fasting is screened not inferred

# Loop's own presets, plus the Walsh models its era also offered.
LOOP_MODELS = ("loop_adult", "loop_child", "loop_fiasp",
               "walsh_3h", "walsh_4h", "walsh_5h", "walsh_6h")
DEFAULT_LOOP_MODEL = "loop_adult"
# Stored on every window so the model choice is a later select, not a rebuild.
CACHED_MODELS = LOOP_MODELS + ("oref_6h75",)

TDD_BANDS = [(0, 20), (20, 40), (40, 64), (64, 1000)]   # 64 U/day is where v1 and v2 cross
AGE_BANDS = [(0, 6, "2-5"), (6, 13, "6-12"), (13, 18, "13-17"),
             (18, 45, "18-44"), (45, 200, "45+")]


def ensure_dirs() -> None:
    for d in (CACHE, GRID_CACHE, WINDOW_CACHE, RESULTS, CHARTS):
        d.mkdir(parents=True, exist_ok=True)
