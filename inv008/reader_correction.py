"""What the reader correction of 23 November 2026 changed (REVIEW.md item 20).

Two comparisons, each against the files as they stood before the correction, which were kept
under $DYNISF_ROOT with the suffix uncorrected_20261123:
1. Bolus insulin per day for each AAPS and OpenAPS person, as the old extract_treatments_tdd.py
   read it (every copy of a repeated bolus, and deleted AAPS records) against the corrected read.
2. The median per-decision daily dose of each person in the INV-008 replay cache, after over
   before, which is what V1 and V2 were given. For AAPS this combines the two corrections, each
   bolus counted once and temp basals read where the profile basal stood before.
The comparison holds the source tables as they now are, so it also carries any change in them
since the first run (REVIEW.md item 19); that touches the Trio results, which neither part reads.
The timestamp parser was corrected the same day (inv008/sources.py): times written with
milliseconds and an offset had been dropped. Part 3 counts, for each oref0 person, the bolus
and temp-basal records whose time the old parser (the string cut at 25 characters) failed on
and the corrected one reads, as totals.

    DYNISF_ROOT=... python3 -m inv008.reader_correction     writes results/reader_correction.md
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np

from inv008 import config, sources

BEFORE = "uncorrected_20261123"
OUT = Path(__file__).resolve().parent.parent / "results" / "reader_correction.md"
NAMES = {"v6": "AAPS", "v7": "OpenAPS"}


def old_parse(s: str):
    """The parser before 23 November 2026, kept only to count what it dropped."""
    try:
        return datetime.fromisoformat(s if s.endswith("Z") else s[:25])
    except ValueError:
        return None


def bolus_section() -> list[str]:
    old = json.loads(
        (config.ROOT / f"backup_{BEFORE}" / "user_treatments_tdd.json").read_text()
    )
    new = json.loads((config.ROOT / "user_treatments_tdd.json").read_text())
    rows = [
        (
            u,
            n["source"].split("/")[0],
            *[r["total_insulin_units"] / r["span_days"] for r in (old[u], n)],
        )
        for u, n in new.items()
        if u in old
    ]
    L = [
        "## Bolus insulin per day, as read before over the corrected read",
        "",
        "Before counts every copy of a repeated bolus and deleted AAPS records; the corrected read counts",
        "each once and also parses times written with milliseconds and an offset, which were dropped.",
        "",
        "| Cohort | People | Median ratio | People above 1.1 | Highest ratio |",
        "|---|---:|---:|---:|---:|",
    ]
    for c, name in NAMES.items():
        r = np.array([a / b for _, p, a, b in rows if p == c and b > 0])
        L.append(
            f"| {name} | {len(r)} | {np.median(r):.2f} | {(r > 1.1).sum()} | {r.max():.1f} |"
        )
    return L


def dose_section() -> list[str]:
    L = [
        "## Daily dose given to V1 and V2, after over before (AAPS: bolus counted once and temp basals read)",
        "",
        "The source tables are as they now are, so Trio rows moved by drift since the first run are not part",
        "of this comparison, which covers AAPS and OpenAPS only (REVIEW.md item 19).",
        "",
        "Per person, the median per-decision daily dose in the INV-008 replay cache after the correction",
        "divided by the one before it; then the median and range over people.",
        "",
        "| Cohort | People | Median | Lowest | Highest |",
        "|---|---:|---:|---:|---:|",
    ]
    by = {c: [] for c in NAMES}
    for f in sorted((config.CACHE / "replay").glob("*.meta.json")):
        n = json.loads(f.read_text())
        o_path = config.ROOT / f"inv008_cache_{BEFORE}" / "replay" / f.name
        o = json.loads(o_path.read_text()) if o_path.exists() else {}
        if n["platform"] in by and n.get("median_tdd") and o.get("median_tdd"):
            by[n["platform"]].append(n["median_tdd"] / o["median_tdd"])
    for c, name in NAMES.items():
        r = np.array(by[c])
        L.append(
            f"| {name} | {len(r)} | {np.median(r):.2f} | {r.min():.2f} | {r.max():.2f} |"
        )
    return L


def parser_counts(raw: str) -> dict:
    """[dropped by the old parser, read by the new] for boluses (b) and temp basals (t)."""
    n = {"b": [0, 0], "t": [0, 0]}
    for f in sorted(
        (config.NS_SAMPLES / raw / "direct-sharing-31").glob("*treatments*.json")
    ):
        try:
            recs = json.loads(f.read_text())
        except (OSError, ValueError):  # unreadable or malformed dump file
            continue
        for r in recs if isinstance(recs, list) else [recs]:
            if not isinstance(r, dict):
                continue
            t = r.get("timestamp") or r.get("created_at") or r.get("date")
            kind = (
                "t"
                if r.get("eventType") == "Temp Basal"
                else "b"
                if (r.get("insulin") or 0) > 0
                else None
            )
            if kind and sources.parse_ts(t) is not None:
                n[kind][1] += 1
                n[kind][0] += isinstance(t, str) and old_parse(t) is None
    return n


def parser_section() -> list[str]:
    """Totals only: per-person results are not published."""
    people, dropped, read = 0, {"b": 0, "t": 0}, {"b": 0, "t": 0}
    for _, (platform, raw) in sorted(sources.load_mappings().items()):
        n = parser_counts(raw) if platform == "v7" else None
        if n and (n["b"][0] or n["t"][0]):
            people += 1
            for k in dropped:
                dropped[k] += n[k][0]
                read[k] += n[k][1]
    return [
        "## OpenAPS records the old timestamp parser dropped",
        "",
        f"{people} people had records the old parser dropped: {dropped['b']:,} of their {read['b']:,} bolus",
        f"records and {dropped['t']:,} of their {read['t']:,} temp-basal records. Counts include repeats",
        "across files.",
    ]


def main() -> None:
    L = [
        "# The reader correction of 23 November 2026",
        "",
        "Generated by `inv008/reader_correction.py`. See REVIEW.md item 20.",
        "",
        *bolus_section(),
        "",
        *dose_section(),
        "",
        *parser_section(),
    ]
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
