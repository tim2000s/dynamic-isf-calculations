"""Copy result files from the private data folder into the public repository without per-person results.

Per-person results are never published (decided 23 November 2026). The scripts write their full
output to $DYNISF_ROOT/results, which stays private; this copies the named files into results/ of
the repository with every per-person part removed: a JSON key holding per-person records
(`per_person`, `per_user`) or a list of participant IDs is dropped, a JSON file that is nothing but per-person records is not
published, and a Markdown table row that opens with a participant ID is removed, with its table
when no row is left. After copying, every published file is checked for participant IDs and the
run fails if one remains.

    DYNISF_ROOT=... python3 -m inv008.publish NAME [NAME ...]     e.g. overnight_sensitivity.json
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from inv008 import config

REPO = Path(__file__).resolve().parent.parent / "results"
PERSON = re.compile(r"\bU\d{3}\b")
PERSON_KEYS = ("per_person", "per_user")


def is_id_list(v) -> bool:
    return (
        isinstance(v, list)
        and bool(v)
        and all(isinstance(x, str) and PERSON.fullmatch(x) for x in v)
    )


def strip_json(d):
    """The JSON without its per-person records, or None if it holds nothing else."""
    if isinstance(d, list):
        return None
    return {k: v for k, v in d.items() if k not in PERSON_KEYS and not is_id_list(v)}


def strip_md(text: str) -> str:
    """Markdown without table rows that open with a participant ID, or tables left empty."""
    out, table = [], []
    for line in [*text.split("\n"), ""]:
        if line.startswith("|"):
            table.append(line)
            continue
        rows = [r for r in table[2:] if not re.match(r"\|\s*U\d{3}\b", r)]
        if table and (rows or len(table) <= 2):
            out += table[:2] + rows
        table = []
        out.append(line)
    return "\n".join(out[:-1])


def publish(name: str) -> str:
    src, dst = config.ROOT / "results" / name, REPO / name
    if name.endswith(".json"):
        d = strip_json(json.loads(src.read_text()))
        if d is None:
            dst.unlink(missing_ok=True)
            return f"{name}: per-person only, not published"
        text = json.dumps(d, indent=1) + "\n"
    else:
        text = strip_md(src.read_text())
    if PERSON.search(text):
        raise SystemExit(
            f"publish: {name} still holds a participant ID after stripping"
        )
    dst.write_text(text)
    return f"{name}: published"


if __name__ == "__main__":
    for n in sys.argv[1:]:
        print(publish(n))
