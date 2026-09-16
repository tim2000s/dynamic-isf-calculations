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
# V2 is the author-confirmed no-+1 form. The AndroidAPS commit cited in the
# original investigation plan contained the earlier +1 implementation and is
# therefore retained as historical provenance only, not as the formula source.
SOURCE_COMMIT = "dynisf-v1-v2-author-confirmed-2026-09"

# ---- formula defaults ----
NORMAL_TARGET = 99.0      # normal target (mg/dL)
BG_CAP = 210.0            # glucose cap (mg/dL); excess above cap at 1/3 weight
VELOCITY = 1.0            # v1 glucose-response damping, held at default (full scaler)
ADJUST_FACTOR = 1.0       # TDD adjustment factor (100%)
# Standard rapid-acting analogue configuration, including NovoRapid. V2 floors
# glucose at divisor + 1, so this configuration becomes 75/76. Other insulin
# configurations must pass their own divisor explicitly.
INSULIN_DIVISOR = 75

# ---- TDD reconstruction ----
GRID_SEC = 300                    # 5-min delivery grid
MIN_DAYS_FOR_7D = 3               # min valid calendar days to form the 7d average
ANCHOR_HOUR_MISMATCH_MAX = 0.01   # max fraction of sampled ticks whose hour-of-day
                                  # disagrees with DB `hour` before anchor is rejected

# ---- parallelism (Mac mini M4 Pro: 10P + 4E cores, 64 GB) ----
DEFAULT_WORKERS = int(os.environ.get("INV008_WORKERS", "7"))
MAXTASKSPERCHILD = 4              # recycle workers to bound pandas memory creep

PLATFORMS = ("v5", "v6", "v7")


def ensure_dirs() -> None:
    for d in (TDD_DIR, REPLAY_DIR, LOG_DIR):
        d.mkdir(parents=True, exist_ok=True)
