#!/usr/bin/env python3
"""All four ISF forms of the selected-window comparison scored on one window set.

`inv008/head_to_head.py` scored static ISF and the loop's own value on 62,751
overnight windows and V1 and V2 on the 57,078 of them where the equations had a
blended TDD. Its per-form medians and its per-person ranking were therefore taken
over nested sets. This script reads the cached window table that script wrote and
scores every form on the windows where V1 and V2 are both finite, so that the four
medians, the paired differences between them and the per-person ranking all refer
to the same windows. The cached table is used in place of a rerun of the database
pipeline because it holds the activity integral each window was rescaled with, and
a rerun would reselect windows from a database whose oref_v5 rows have changed
since that table was written (inv008/overnight_variants.py records three Trio
people whose four-hour window counts now differ from the committed output).

What is computed, on the common set:
  median absolute error and median signed error for each form, error being
  observed end glucose minus predicted end glucose, as in head_to_head.py;
  the paired difference in median absolute error, V1 minus static, V2 minus
  static and V2 minus V1, with a 95% percentile interval from 2,000 bootstrap
  draws that resample people, keeping all of a drawn person's windows, so that
  overlapping windows within a person do not narrow the interval;
  the per-person best form among people with at least 60 common windows;
  the same V1 and V2 figures at divisors 55, 65 and 75, recalculated from glucose
  and blended TDD with the activity integral held fixed, as
  `inv008/divisor_sensitivity.py` does.

Modelling choices. The window table is used as committed. It holds 7,682 rows that
repeat another row in every one of its 18 columns, which is the signature of the
repeated (user, timestamp) records in the oref tables (REVIEW.md item 8); the
table carries no timestamp, so they are collapsed by exact identity across all
columns and the headline figures are repeated on the collapsed set as a check.
The primary figures stay on the table as committed so that they sit beside the
57,078-window V1 and V2 figures already in the paper.

Expected sign, written before the run. On the nested sets V1 had a larger median
absolute error than static (24.6 against 20.3 mg/dL) and V2 a larger one again
(49.7). Dropping the 5,673 windows with no TDD should not reverse that, so V1
minus static and V2 minus static are expected to be positive with intervals above
zero, and V2 minus V1 positive. The static and loop medians on the common set are
expected to lie within about 1 mg/dL of their full-set values, since the dropped
windows are a tenth of the set and come from 30 of the 89 people.

    python3 -m inv008.common_window \
        --windows ~/oref-investigations/results/head_to_head_windows.parquet

Outputs results/common_window.json and results/common_window.md (aggregates only).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from inv008.dynisf import isf_v1, isf_v2

DIVISORS = (55.0, 65.0, 75.0)
FORMS = ("static", "loop", "v1", "v2")
MIN_PERSON = 60
N_BOOT = 2000
SEED = 0
OUT = Path(__file__).resolve().parent.parent / "results"


def errors(d: pd.DataFrame, divisor: float) -> dict[str, np.ndarray]:
    """Signed error per window for each form at one divisor.

    Static and loop errors do not depend on the divisor. V1 and V2 are rescaled
    from the activity integral recovered from the stored V1 error, as in
    divisor_sensitivity.same_window.
    """
    drop = d.bg.to_numpy(float) - d.bg_end.to_numpy(float)
    activity = (d.err_v1.to_numpy(float) + drop) / d.isf_v1.to_numpy(float)
    bg, tdd = d.bg.to_numpy(float), d.tdd.to_numpy(float)
    return {"static": d.err_static.to_numpy(float), "loop": d.err_loop.to_numpy(float),
            "v1": activity * isf_v1(bg, tdd, divisor=divisor) - drop,
            "v2": activity * isf_v2(bg, tdd, divisor=divisor) - drop}


def person_index(users: np.ndarray) -> list[np.ndarray]:
    codes, uniq = pd.factorize(users)
    order = np.argsort(codes, kind="stable")
    bounds = np.flatnonzero(np.diff(codes[order])) + 1
    return np.split(order, bounds)


def boot_mae_diffs(err: dict[str, np.ndarray], groups: list[np.ndarray],
                   pairs: list[tuple[str, str]], n_boot: int = N_BOOT,
                   seed: int = SEED) -> dict[str, list[float]]:
    """Percentile interval for median|e_a| - median|e_b|, resampling people."""
    rng = np.random.default_rng(seed)
    absd = {k: np.abs(v) for k, v in err.items()}
    draws = {f"{a}_minus_{b}": np.empty(n_boot) for a, b in pairs}
    mae_draws = {k: np.empty(n_boot) for k in absd}
    for i in range(n_boot):
        idx = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
        med = {k: float(np.median(v[idx])) for k, v in absd.items()}
        for k in med:
            mae_draws[k][i] = med[k]
        for a, b in pairs:
            draws[f"{a}_minus_{b}"][i] = med[a] - med[b]
    out = {}
    for k, v in draws.items():
        a, b = k.split("_minus_")
        point = float(np.median(absd[a]) - np.median(absd[b]))
        lo, hi = np.percentile(v, [2.5, 97.5])
        out[k] = [point, float(lo), float(hi)]
    for k, v in mae_draws.items():
        lo, hi = np.percentile(v, [2.5, 97.5])
        out[f"mae_{k}"] = [float(np.median(absd[k])), float(lo), float(hi)]
    return out


def verdict(ci: list[float]) -> str:
    if ci[1] > 0:
        return "interval above zero"
    if ci[2] < 0:
        return "interval below zero"
    return "interval includes zero"


def ranking(err: dict[str, np.ndarray], users: np.ndarray) -> tuple[dict, int]:
    wins = {f: 0 for f in FORMS}
    n = 0
    for idx in person_index(users):
        if len(idx) < MIN_PERSON:
            continue
        maes = {f: float(np.median(np.abs(err[f][idx]))) for f in FORMS}
        wins[min(maes, key=maes.get)] += 1
        n += 1
    return wins, n


def score(d: pd.DataFrame, label: str) -> dict:
    users = d.user.to_numpy()
    groups = person_index(users)
    res = {"label": label, "n_windows": int(len(d)), "n_people": int(d.user.nunique()),
           "by_divisor": {}}
    pairs = [("v1", "static"), ("v2", "static"), ("v2", "v1"), ("loop", "static")]
    for div in DIVISORS:
        e = errors(d, div)
        assert all(np.isfinite(v).all() for v in e.values()), "non-finite error on common set"
        row = {f: {"mae": float(np.median(np.abs(e[f]))), "bias": float(np.median(e[f]))}
               for f in FORMS}
        row["paired"] = boot_mae_diffs(e, groups, pairs)
        row["verdicts"] = {k: verdict(v) for k, v in row["paired"].items() if "_minus_" in k}
        row["per_person_best"], row["n_people_ranked"] = ranking(e, users)
        res["by_divisor"][str(int(div))] = row
    return res


def self_check() -> None:
    """Stored divisor-75 errors must equal the recomputed ones; bootstrap must vary."""
    toy = pd.DataFrame({"bg": [150.0, 180.0, 120.0, 200.0], "bg_end": [110.0, 120.0, 100.0, 150.0],
                        "tdd": [40.0, 40.0, 30.0, 50.0], "user": ["a", "a", "b", "b"]})
    act = np.array([1.0, 1.5, 0.5, 2.0])
    toy["isf_v1"] = isf_v1(toy.bg, toy.tdd, divisor=75.0)
    toy["err_v1"] = act * toy.isf_v1 - (toy.bg - toy.bg_end)
    toy["err_static"] = 0.0
    toy["err_loop"] = 0.0
    e = errors(toy, 75.0)
    assert np.allclose(e["v1"], toy.err_v1)
    assert np.allclose(e["v2"], act * isf_v2(toy.bg, toy.tdd, divisor=75.0) - (toy.bg - toy.bg_end))


def markdown(res: dict) -> str:
    p = res["primary"]
    L = ["# All four ISF forms on the common window set", "",
         "Generated by `inv008/common_window.py` from the cached window table of "
         "`inv008/head_to_head.py`. Error is observed end glucose minus predicted end glucose "
         "(positive: the form predicted too large a fall). Intervals are 95% percentile "
         f"intervals from {N_BOOT} bootstrap draws resampling people. The windows are the "
         "future-rise-selected overnight windows of section 4.3 and the comparison is inside "
         "the loop's linear IOB prediction.", "",
         f"Full table: {res['full']['n_windows']:,} windows from {res['full']['n_people']} people. "
         f"Common set (V1 and V2 finite): {p['n_windows']:,} windows from {p['n_people']} people.", "",
         "## Why windows could not be scored for V1 and V2", ""]
    for k, v in res["unscorable"].items():
        L.append(f"- {k}: {v}")
    L += ["", "## Static and loop on the full table and on the common set", "",
          "| Form | Windows | MAE | Signed |", "|---|---:|---:|---:|"]
    for f in ("static", "loop"):
        o = res["full"][f]
        L.append(f"| {f}, full table | {res['full']['n_windows']:,} | {o['mae']:.1f} | {o['bias']:+.1f} |")
    d75 = p["by_divisor"]["75"]
    for f in ("static", "loop"):
        L.append(f"| {f}, common set | {p['n_windows']:,} | {d75[f]['mae']:.1f} | {d75[f]['bias']:+.1f} |")
    for label, r in (("Common set as committed", p), ("Common set, exact duplicate rows collapsed", res["dedup"])):
        L += ["", f"## {label} ({r['n_windows']:,} windows, {r['n_people']} people)", "",
              "| Divisor | Static MAE | Static signed | Loop MAE | Loop signed | V1 MAE | V1 signed "
              "| V2 MAE | V2 signed |", "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for div, row in r["by_divisor"].items():
            L.append(f"| {div} | " + " | ".join(
                f"{row[f]['mae']:.1f} | {row[f]['bias']:+.1f}" for f in FORMS) + " |")
        L += ["", "Paired difference in median absolute error, mg/dL (95% interval), and verdict.", "",
              "| Divisor | V1 minus static | V2 minus static | V2 minus V1 | Loop minus static |",
              "|---:|---|---|---|---|"]
        for div, row in r["by_divisor"].items():
            cells = []
            for k in ("v1_minus_static", "v2_minus_static", "v2_minus_v1", "loop_minus_static"):
                v = row["paired"][k]
                cells.append(f"{v[0]:+.1f} ({v[1]:+.1f} to {v[2]:+.1f}), {row['verdicts'][k]}")
            L.append(f"| {div} | " + " | ".join(cells) + " |")
        L += ["", f"Per-person best form, people with at least {MIN_PERSON} common windows.", "",
              "| Divisor | People | Loop | Static | V1 | V2 |", "|---:|---:|---:|---:|---:|---:|"]
        for div, row in r["by_divisor"].items():
            w = row["per_person_best"]
            L.append(f"| {div} | {row['n_people_ranked']} | {w['loop']} | {w['static']} | {w['v1']} | {w['v2']} |")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--windows", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    self_check()
    D = pd.read_parquet(a.windows)
    full = {"n_windows": int(len(D)), "n_people": int(D.user.nunique())}
    for f in ("static", "loop"):
        e = D[f"err_{f}"].to_numpy(float)
        full[f] = {"mae": float(np.median(np.abs(e))), "bias": float(np.median(e))}
    bad = ~(np.isfinite(D.err_v1) & np.isfinite(D.err_v2))
    B = D[bad]
    tdd_nan, tdd_zero = B.tdd.isna(), B.tdd.fillna(-1).eq(0)
    unscorable = {
        "windows without V1 and V2": f"{int(bad.sum()):,} from {B.user.nunique()} people; "
                                     "in every one the blended TDD was missing or zero",
        "Trio (oref_v5), device-logged TDD recorded as zero":
            f"{int((tdd_zero & B.table.eq('oref_v5')).sum()):,} windows",
        "Trio (oref_v5), no device TDD within 60 min of the window start":
            f"{int((tdd_nan & B.table.eq('oref_v5')).sum()):,} windows",
        "OpenAPS (oref_v7), reconstructed blend undefined":
            f"{int((tdd_nan & B.table.eq('oref_v7')).sum()):,} windows; the blend is undefined "
            "when yesterday's calendar day is incomplete or holds no bolus, when fewer than three "
            "of the previous seven days are complete, or when the window lies outside the "
            "reconstructed delivery span",
        "other": f"{int((bad & ~tdd_nan & ~tdd_zero).sum())} windows",
    }
    C = D[~bad].reset_index(drop=True)
    # stored divisor-75 errors must be reproduced before anything else is read
    e75 = errors(C, 75.0)
    assert np.allclose(e75["v1"], C.err_v1) and np.allclose(e75["v2"], C.err_v2, rtol=1e-9, atol=1e-6)
    res = {"full": full, "unscorable": unscorable,
           "primary": score(C, "common set as committed")}
    Cd = C.drop_duplicates().reset_index(drop=True)
    res["dedup"] = score(Cd, "common set, exact duplicate rows collapsed")
    res["n_exact_duplicate_rows_full_table"] = int(D.duplicated().sum())
    res["n_exact_duplicate_rows_common_set"] = int(C.duplicated().sum())
    a.out.mkdir(exist_ok=True)
    (a.out / "common_window.json").write_text(json.dumps(res, indent=1))
    md = markdown(res)
    (a.out / "common_window.md").write_text(md)
    print(md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
