"""Shared configuration for the dynamic-ISF equation analysis."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
NS_SAMPLES = ROOT / "240 NS samples"
CACHE = ROOT / "inv008_cache"
TDD_DIR = CACHE / "tdd"
REPLAY_DIR = CACHE / "replay"
LOG_DIR = CACHE / "logs"

DB = dict(host="localhost", dbname="oref")

USER_MAPPING_V6 = ROOT / "user_mapping_v6.json"
USER_MAPPING_V7 = ROOT / "user_mapping_v7.json"
BASAL_PROFILES = ROOT / "user_basal_profiles.json"

# Tag recorded in per-user outputs to identify the equation definitions used.
SOURCE_COMMIT = "dynisf-v1-v2-2026-06"

# ---- formula defaults ----
NORMAL_TARGET = 99.0      # normal target (mg/dL)
BG_CAP = 210.0            # glucose cap (mg/dL); excess above cap at 1/3 weight
VELOCITY = 1.0            # v1 glucose-response damping, held at default (full scaler)
ADJUST_FACTOR = 1.0       # TDD adjustment factor (100%)
INSULIN_DIVISOR = 75      # Lyumjev (peak 45): (90-45)+30. Fiasp=65, rapid=55.

# ---- TDD reconstruction ----
GRID_SEC = 300                    # 5-min delivery grid
MIN_DAYS_FOR_7D = 3               # min valid calendar days to form the 7d average
ANCHOR_HOUR_MISMATCH_MAX = 0.01   # max fraction of sampled ticks whose hour-of-day
                                  # disagrees with DB `hour` before anchor is rejected

# ---- parallelism (Mac mini M4 Pro: 10P + 4E cores, 64 GB) ----
DEFAULT_WORKERS = int(os.environ.get("INV008_WORKERS", "12"))
MAXTASKSPERCHILD = 4              # recycle workers to bound pandas memory creep

PLATFORMS = ("v5", "v6", "v7")


def ensure_dirs() -> None:
    for d in (TDD_DIR, REPLAY_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)
