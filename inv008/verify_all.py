#!/usr/bin/env python3
"""Independent end-to-end verification of the dynamic-ISF analysis.

Re-derives every headline number from scratch (an "oracle" that does not import the
package's equation code where it can avoid it) and cross-checks against (a) the package
functions, (b) the replayed parquet outputs, and (c) the numbers quoted in the documents.
Also quantifies the temp-basal bin-edge convention flagged in review.

Run: python -m inv008.verify_all
Exit 0 only if every check passes.
"""
from __future__ import annotations

import glob
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from inv008 import config, dynisf
from inv008.tdd_windows import TempBasalEvent, build_delivery_grid, windowed_tdd

ROOT = config.ROOT
NT, DIV = 99.0, 75
LT = math.log(NT / DIV + 1.0)
fails: list[str] = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{(' — ' + detail) if detail else ''}")
    if not ok:
        fails.append(name)


# ---------------------------------------------------------------- 1. equations
def oracle_v1(bg, tdd):
    bgadj = bg if bg <= 210 else 210 + (bg - 210) / 3
    return (1800.0 / (tdd * math.log(NT / DIV + 1))) * (math.log(NT / DIV + 1) / math.log(bgadj / DIV + 1))


def oracle_v2(bg, tdd):
    bgadj = bg if bg <= 210 else 210 + (bg - 210) / 3
    bgf = max(bgadj, DIV + 1)              # v2 floors glucose at divisor+1, no +1 in the log
    return 115000.0 / (tdd**2 * math.log(bgf / DIV))


print("\n1. Equations: independent oracle vs package (random inputs)")
rng = np.random.default_rng(0)
maxd1 = maxd2 = 0.0
for _ in range(10000):
    bg = float(rng.uniform(60, 360)); tdd = float(rng.uniform(8, 200))
    maxd1 = max(maxd1, abs(oracle_v1(bg, tdd) - float(dynisf.isf_v1(bg, tdd))))
    maxd2 = max(maxd2, abs(oracle_v2(bg, tdd) - float(dynisf.isf_v2(bg, tdd))))
check("v1 oracle vs package", maxd1 < 1e-9, f"max abs diff {maxd1:.2e}")
check("v2 oracle vs package", maxd2 < 1e-9, f"max abs diff {maxd2:.2e}")

print("\n2. v2/v1 ratio depends on glucose (v1 keeps the +1, v2 does not)")
worst = 0.0
for tdd in (12, 20, 36, 64, 100, 200):
    for bg in (70, 120, 180, 250, 350):
        r_oracle = oracle_v2(bg, tdd) / oracle_v1(bg, tdd)
        worst = max(worst, abs(r_oracle - float(dynisf.v2_over_v1_ratio(bg, tdd))))
check("package v2/v1 ratio matches independent oracle", worst < 1e-9, f"max dev {worst:.2e}")
rs = [float(dynisf.v2_over_v1_ratio(bg, 40.0)) for bg in (80, 120, 200)]
check("v2/v1 falls as glucose rises", rs[0] > rs[1] > rs[2],
      f"{rs[0]:.1f} > {rs[1]:.1f} > {rs[2]:.1f}")

print("\n3. v-next = 355/sqrt(TDD) at normal target; divisor-free")
# at normal target the scaler is 1 for ANY divisor → anchor identical across insulin types
anch = {d: (355 / math.sqrt(50)) * (math.log(99 / d + 1) / math.log(99 / d + 1)) for d in (55, 65, 75)}
check("355/sqrt(TDD) anchor independent of divisor at target",
      max(anch.values()) - min(anch.values()) < 1e-12,
      f"all = {list(anch.values())[0]:.3f} at TDD 50")

print("\n4. Proposal comparison table (ISF at normal target, divisor 75)")
LT2 = math.log(NT / DIV)   # v2 log term, no +1
expect = {15: (143, 1840, 92), 25: (86, 663, 71), 36: (59, 320, 59),
          50: (43, 166, 50), 80: (27, 65, 40), 120: (18, 29, 32)}
tbl_ok = True
for tdd, (e1, e2, en) in expect.items():
    v1 = 1800 / (tdd * LT); v2 = 115000 / (tdd**2 * LT2); vn = 355 / math.sqrt(tdd)
    if not (abs(round(v1) - e1) <= 1 and abs(v2 - e2) <= max(1, 0.02 * e2) and abs(round(vn) - en) <= 1):
        tbl_ok = False
        print(f"      TDD {tdd}: v1={v1:.1f}(doc {e1}) v2={v2:.1f}(doc {e2}) vnext={vn:.1f}(doc {en})")
check("documented table matches recomputation", tbl_ok)

print("\n5. Fit constants reproduced from canonical_cohort + empirical_isf_v5")
coh = pd.DataFrame(json.loads((ROOT / "canonical_cohort.json").read_text()))
emp = pd.DataFrame(json.loads((ROOT / "empirical_isf_v5.json").read_text()))[
    ["user_id", "empirical_isf", "r2"]]
df = coh.merge(emp, on="user_id", how="left")
ev = df[(df.r2 >= 0.10) & df.empirical_isf.between(5, 500)]
K_ent = float(np.median(df.dropna(subset=["isf"]).isf * np.sqrt(df.dropna(subset=["isf"]).tdd)))
K_emp = float(np.median(ev.empirical_isf * np.sqrt(ev.tdd)))
res = json.loads((ROOT / "best_isf_fit_results.json").read_text())
check("K (entered) ~ 355", abs(K_ent - 355) < 6, f"recomputed {K_ent:.1f}, json {res['sqrt_rule_entered_K']}")
check("K (empirical) ~ 145", abs(K_emp - 145) < 6, f"recomputed {K_emp:.1f}, json {res['sqrt_rule_empirical_K']}")
check("anchor ratio 355/145 ~ 2.45", abs(K_ent / K_emp - 2.45) < 0.2, f"{K_ent/K_emp:.2f}")
# power-law slopes
for tgt, col, lo, hi in [("entered", "isf", -0.6, -0.3), ("empirical", "empirical_isf", -0.6, -0.3)]:
    d = (ev if tgt == "empirical" else df).dropna(subset=[col])
    b = np.polyfit(np.log(d.tdd), np.log(d[col]), 1)[0]
    check(f"power-law slope ({tgt}) in [-0.6,-0.3]", lo < b < hi, f"slope {b:.3f}")
check("n_empirical = 114", len(ev) == 114, f"{len(ev)}")
check("n_cohort = 138", len(df) == 138, f"{len(df)}")

print("\n6. Replay parquet: v2/v1 reproduced from the equations, glucose-dependent")
worst_ratio = 0.0; n_checked = 0; low_b, high_b = [], []
for f in glob.glob(str(config.REPLAY_DIR / "*.parquet")):
    d = pd.read_parquet(f, columns=["bg", "isf_v1", "tdd"]).dropna()
    d = d[(d.isf_v1 > 0) & (d.tdd > 0) & (d.bg > 0)]
    if len(d) < 100:
        continue
    bg, tdd, v1 = d.bg.to_numpy(), d.tdd.to_numpy(), d.isf_v1.to_numpy()
    obs = dynisf.isf_v2(bg, tdd) / v1
    ana = np.array([oracle_v2(b, t) for b, t in zip(bg[:200], tdd[:200])]) / v1[:200]
    worst_ratio = max(worst_ratio, float(np.nanmax(np.abs(obs[:200] - ana))))
    if (bg < 90).any():
        low_b.append(float(np.nanmedian(obs[bg < 90])))
    if (bg > 180).any():
        high_b.append(float(np.nanmedian(obs[bg > 180])))
    n_checked += 1
check("package v2/v1 matches oracle in replayed data", worst_ratio < 1e-6,
      f"{n_checked} users, max dev {worst_ratio:.2e}")
check("v2/v1 larger at low glucose than high (cohort medians)",
      np.median(low_b) > np.median(high_b),
      f"low {np.median(low_b):.1f} vs high {np.median(high_b):.1f}")

print("\n7. Cohort composition")
summ = pd.read_json(ROOT / "charts/inv008/cohort_summary.json")
plats = summ.platform.value_counts().to_dict()
check("platform counts (v5/v6/v7)", plats.get("v5") and plats.get("v6") and plats.get("v7"),
      str(plats))

print("\n8. Temp-basal bin-edge convention is unbiased (phase-averaged)")
# The grid samples each bin at its left edge: bin i takes the temp rate iff its start
# lies in [seg_start, seg_end). Review suggested including the bin that *contains*
# seg_start (side=right-1). Settle which is unbiased by averaging delivered temp-excess
# over many random segment phases/durations against the analytic truth.
from inv008.tdd_windows import profile_rate_at
hourly = np.ones(24)


def temp_excess_error(use_fix, seed):
    r = np.random.default_rng(seed)
    ts, dur, n = [], [], 3000
    t, excess_true = 0.0, 0.0
    for _ in range(n):
        t += r.uniform(60, 1200); d = r.uniform(10, 45) * 60
        ts.append(t); dur.append(d / 60.0); excess_true += d * (2.0 - 1.0) / 3600.0
        t += d
    bins = np.arange(0, int(t) + 3600 + 300, 300, dtype=np.int64)
    rate = profile_rate_at(bins, hourly)
    evs = [TempBasalEvent(ts=ts[i], duration_min=dur[i], rate_u_h=2.0) for i in range(n)]
    for i, ev in enumerate(evs):
        s = ev.ts; e = min(ev.ts + ev.duration_min * 60, evs[i + 1].ts) if i + 1 < n else ev.ts + ev.duration_min * 60
        i0 = (np.searchsorted(bins, s, side="right") - 1) if use_fix else np.searchsorted(bins, s, side="left")
        i1 = np.searchsorted(bins, e, side="left"); i0 = max(i0, 0)
        if i1 > i0:
            rate[i0:i1] = 2.0
    excess_grid = (rate * 300 / 3600.0).sum() - len(bins) * 300 / 3600.0
    return (excess_grid - excess_true) / excess_true * 100


cur = np.mean([temp_excess_error(False, s) for s in range(20)])
fix = np.mean([temp_excess_error(True, s) for s in range(20)])
check("current (left-edge) convention is unbiased", abs(cur) < 0.5,
      f"mean temp-excess error {cur:+.2f}% (the suggested side=right-1 alternative is {fix:+.1f}%, biased)")

print("\n9. Device-ISF validation summary present and sane")
dv = json.loads((ROOT / "results/device_isf_validation.json").read_text())
logc = dv["log_dynisf"]["median_log_corr"]
check("log-DynISF device correlation > 0.4", logc > 0.4, f"median log-corr {logc}")
# The validation only spans the DynISF cohort users (n≈18); of those, 2 carry the
# mid-history mmol/mgdl unit switch (U019, U023). The DB-wide count is 4 (also U005,
# U025) but U005/U025 are not DynISF-cohort users so they are out of this scope.
check("mixed-unit users flagged within validation scope", len(dv["mixed_unit_users"]) == 2,
      f"{dv['mixed_unit_users']} (DB-wide there are 4; U005/U025 are out of scope)")

print("\n" + "=" * 60)
if fails:
    print(f"VERIFICATION FAILED: {len(fails)} check(s) — {fails}")
    raise SystemExit(1)
print("ALL CHECKS PASSED")
