#!/usr/bin/env python3
"""Extract per-user treatments-derived TDD from raw NS / AAPS dumps.

Replaces the hybrid TDD `max(2 × basal, basal + SMB/day, sug_tdd)` with
ground-truth `basal_profile + sum(treatments insulin) / span_days`.

Sources:
  • v7 (oref0):        /240 NS samples/<dir>/direct-sharing-31/<id>_treatments_*.json
  • v6 (AAPS classic): /240 NS samples/<dir>/direct-sharing-396/upload-*/Treatments.json
  • v5 (Trio):         /240 NS samples/<dir>/direct-sharing-31/...  (same NS format as v7
                       where present)

Boluses are read with inv008/sources.py. The basal term added downstream (canonical_cohort.py)
is the profile schedule, so temporary basals are not counted here; a modelling choice.

Output: user_treatments_tdd.json keyed by user_id with fields:
  total_insulin_units, span_days, n_events, source
"""
from __future__ import annotations

import json
import os
from pathlib import Path

from inv008 import sources

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
NS = ROOT / "240 NS samples"
OUT = ROOT / "user_treatments_tdd.json"

MV6 = json.loads((ROOT / "user_mapping_v6.json").read_text())
MV7 = json.loads((ROOT / "user_mapping_v7.json").read_text())


def extract(raw: str, loader, platform: str) -> dict | None:
    """Total bolus insulin and span from the INV-008 readers, which count a bolus repeated
    across overlapping files once and skip deleted or invalid AndroidAPS records. Until
    23 November 2026 this script summed every copy itself (REVIEW.md item 20)."""
    ts, units, _ = loader(raw)
    if not len(ts):
        return None
    span = float(ts.max() - ts.min()) / 86400.0
    if span <= 0:
        return None
    return {
        "total_insulin_units": round(float(units.sum()), 2),
        "span_days": round(span, 2),
        "n_events": len(ts),
        "source": f"{platform}/{raw}",
    }


def main():
    out = {}
    for mapping, loader, platform in ((MV7, sources.load_v7_events, "v7"), (MV6, sources.load_v6_events, "v6")):
        for uid, raw in mapping.items():
            if (NS / raw).is_dir() and (rec := extract(raw, loader, platform)):
                out[uid] = rec
    OUT.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT} with {len(out)} users")
    # Distribution sanity
    spans = sorted(r["span_days"] for r in out.values())
    totals = sorted(r["total_insulin_units"] / r["span_days"] for r in out.values())
    if totals:
        print(f"  per-day insulin (treatments / span): min={totals[0]:.1f}, "
              f"med={totals[len(totals)//2]:.1f}, max={totals[-1]:.1f} U/day")
        print(f"  span_days: min={spans[0]:.0f}, med={spans[len(spans)//2]:.0f}, max={spans[-1]:.0f}")


if __name__ == "__main__":
    main()
