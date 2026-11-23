"""No script imports from the data folder, which holds stale copies of repository modules."""

import re
from pathlib import Path

DATA_ON_PATH = re.compile(r"sys\.path\.insert\(\s*0\s*,\s*str\(\s*(config\.)?ROOT\s*\)")
# a value compared with 20 and multiplied by the mmol/L factor on one line: every ISF in the Commons is mg/dL
BELOW_20 = re.compile(
    r"^[^#\n]*(<\s*20\b[^\n]*(18\.0|MMOL_TO_MGDL)|(18\.0|MMOL_TO_MGDL)[^\n]*<\s*20\b)",
    re.M,
)


def test_no_script_puts_the_data_folder_on_the_import_path():
    repo = Path(__file__).resolve().parents[2]
    bad = [
        str(p.relative_to(repo))
        for p in repo.rglob("*.py")
        if p.name != "test_imports.py"
        and DATA_ON_PATH.search(p.read_text(errors="ignore"))
    ]
    assert bad == []


def test_no_isf_is_converted_by_a_below_20_rule():
    repo = Path(__file__).resolve().parents[2]
    bad = [
        str(p.relative_to(repo))
        for p in repo.rglob("*.py")
        if p.name != "test_imports.py" and BELOW_20.search(p.read_text(errors="ignore"))
    ]
    assert bad == []
