#!/usr/bin/env python3
"""
Deep investigation of multi-site backtest results.
Splits sigmoid vs log, examines <105 vs >105, checks scaling factors,
and investigates why Full Diabeloop outperforms hybrid on Trio sites.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

OUT_DIR = Path.home() / 'Downloads'

# Load raw results
with open(OUT_DIR / 'multisite_backtest_results.json') as f:
    results = json.load(f)

print("=" * 80)
print("DEEP INVESTIGATION — SIGMOID vs LOG, <105 vs >105, SCALING")
print("=" * 80)

# ══════════════════════════════════════════════════════════════════════════════
# 1. SPLIT BY MODEL TYPE
# ══════════════════════════════════════════════════════════════════════════════

sigmoid_sites = [r for r in results if r['model'] == 'sigmoid']
log_sites = [r for r in results if r['model'] == 'log']

def weighted_avg(sites, formula_key, metric):
    total_n = sum(r['n_samples'] for r in sites)
    if total_n == 0:
        return 0
    return sum(r['overall'][formula_key][metric] * r['n_samples'] for r in sites) / total_n

print(f"\n{'═'*80}")
print("1. SIGMOID vs LOG — OVERALL MAE")
print(f"{'═'*80}")

for label, sites in [("SIGMOID", sigmoid_sites), ("LOG", log_sites)]:
    total_n = sum(r['n_samples'] for r in sites)
    print(f"\n  {label} ({len(sites)} sites, {total_n:,} samples)")
    print(f"  {'Site':<14s} {'N':>5s} {'TDD':>5s}  {'Loop':>6s} {'Quart':>6s} {'FullDB':>6s} {'Hybrid':>6s}")
    print(f"  {'─'*14} {'─'*5} {'─'*5}  {'─'*6} {'─'*6} {'─'*6} {'─'*6}")
    for r in sites:
        o = r['overall']
        best_val = min(o['loop']['mae'], o['quartic']['mae'], o['full_db']['mae'], o['hybrid']['mae'])
        def mark(v):
            return f"*{v:5.1f}" if abs(v - best_val) < 0.05 else f" {v:5.1f}"
        print(f"  {r['name']:<14s} {r['n_samples']:5d} {r['tdd_median']:5.1f} "
              f"{mark(o['loop']['mae'])} {mark(o['quartic']['mae'])} "
              f"{mark(o['full_db']['mae'])} {mark(o['hybrid']['mae'])}")
    print(f"  {'Weighted avg':<14s} {total_n:5d} {'':5s} "
          f" {weighted_avg(sites, 'loop', 'mae'):5.1f}  {weighted_avg(sites, 'quartic', 'mae'):5.1f} "
          f" {weighted_avg(sites, 'full_db', 'mae'):5.1f}  {weighted_avg(sites, 'hybrid', 'mae'):5.1f}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. <105 vs ≥105 SPLIT (from band data)
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'═'*80}")
print("2. PERFORMANCE SPLIT: <105 vs ≥105")
print(f"{'═'*80}")

for label, sites in [("SIGMOID", sigmoid_sites), ("LOG", log_sites)]:
    print(f"\n  ── {label} ──")

    # Aggregate <105 bands (<90 + 90-105)
    for zone_label, zone_bands in [("<105", ['<90', '90-105']), ("≥105", ['105-120', '120-150', '150+'])]:
        total_n = 0
        sums = {'loop': 0, 'quartic': 0, 'full_db': 0, 'hybrid': 0}
        bias_sums = {'loop': 0, 'quartic': 0, 'full_db': 0, 'hybrid': 0}
        for r in sites:
            for bname in zone_bands:
                bl = r['bands']['loop'].get(bname, {})
                if not bl:
                    continue
                n = bl['n']
                total_n += n
                for fkey in sums:
                    b = r['bands'][fkey].get(bname, {})
                    if b:
                        sums[fkey] += b['mae'] * n
                        bias_sums[fkey] += b['bias'] * n
        if total_n > 0:
            print(f"  {zone_label} (n={total_n:,}):")
            print(f"    {'Formula':<15s}  {'MAE':>6s}  {'Bias':>7s}")
            print(f"    {'─'*15}  {'─'*6}  {'─'*7}")
            for fkey, fname in [('loop', 'Loop actual'), ('quartic', 'Quartic'),
                                ('full_db', 'Full Diabeloop'), ('hybrid', 'Hybrid')]:
                mae = sums[fkey] / total_n
                bias = bias_sums[fkey] / total_n
                print(f"    {fname:<15s}  {mae:6.1f}  {bias:+7.1f}")

# ══════════════════════════════════════════════════════════════════════════════
# 3. PER-SITE <105 vs ≥105 — SIGMOID SITES
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'═'*80}")
print("3. PER-SITE SPLIT: <105 vs ≥105 — SIGMOID")
print(f"{'═'*80}")

for r in sigmoid_sites:
    print(f"\n  {r['name']} (TDD={r['tdd_median']:.1f}, S_quartic={r['S_quartic']:.3f}, S_hybrid={r['S_hybrid']:.3f})")

    for zone_label, zone_bands in [("<105", ['<90', '90-105']), ("≥105", ['105-120', '120-150', '150+'])]:
        total_n = 0
        sums = {'loop': 0, 'quartic': 0, 'full_db': 0, 'hybrid': 0}
        bias_sums = {'loop': 0, 'quartic': 0, 'full_db': 0, 'hybrid': 0}
        for bname in zone_bands:
            bl = r['bands']['loop'].get(bname, {})
            if not bl:
                continue
            n = bl['n']
            total_n += n
            for fkey in sums:
                b = r['bands'][fkey].get(bname, {})
                if b:
                    sums[fkey] += b['mae'] * n
                    bias_sums[fkey] += b['bias'] * n
        if total_n > 0:
            print(f"    {zone_label} (n={total_n}):")
            for fkey, fname in [('loop', 'Loop'), ('quartic', 'Quartic'),
                                ('full_db', 'FullDB'), ('hybrid', 'Hybrid')]:
                mae = sums[fkey] / total_n
                bias = bias_sums[fkey] / total_n
                print(f"      {fname:<8s} MAE={mae:5.1f}  bias={bias:+6.1f}")

# ══════════════════════════════════════════════════════════════════════════════
# 4. SCALING FACTOR ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'═'*80}")
print("4. SCALING FACTOR ANALYSIS")
print(f"{'═'*80}")

print(f"\n  The scaling factor S = (1800 / TDD_median) / anchor")
print(f"  Quartic anchor = {81.6:.1f} (quartic value at BG=99)")
print(f"  Hybrid anchor  = {75.8:.1f} (quartic value at BG=105)")
print(f"  Full DB anchor = {81.6:.1f} (same as quartic — both use quartic at 99)")
print()

print(f"  {'Site':<14s} {'Model':<8s} {'TDD':>5s}  {'1800/TDD':>8s}  "
      f"{'S_quart':>7s}  {'S_hybrid':>8s}  "
      f"{'ISF@80_Q':>8s} {'ISF@80_H':>8s}  "
      f"{'ISF@99_Q':>8s} {'ISF@99_H':>8s}  "
      f"{'ISF@120_Q':>9s} {'ISF@120_H':>9s}")
print(f"  {'─'*14} {'─'*8} {'─'*5}  {'─'*8}  "
      f"{'─'*7}  {'─'*8}  "
      f"{'─'*8} {'─'*8}  "
      f"{'─'*8} {'─'*8}  "
      f"{'─'*9} {'─'*9}")

# ISF formula values
def isf_quartic(bg):
    G = np.asarray(bg, dtype=float)
    return 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4

def isf_full_diabeloop(bg):
    G = np.asarray(bg, dtype=float)
    quad = 98.03 - 1.077*G + 0.008868*G**2
    quart = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    return np.where(G <= 100, quad, quart)

def isf_hybrid(bg):
    G = np.asarray(bg, dtype=float)
    quart = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    power = 75.8 * (105.0 / G) ** 3.5
    return np.where(G >= 105, quart, power)

for r in results:
    tdd = r['tdd_median']
    isf_1800 = 1800.0 / tdd
    S_q = isf_1800 / 81.6
    S_h = isf_1800 / 75.8

    # Scaled ISF values at key glucose points
    isf_q_80 = float(isf_quartic(80)) * S_q
    isf_h_80 = float(isf_hybrid(80)) * S_h
    isf_q_99 = float(isf_quartic(99)) * S_q
    isf_h_99 = float(isf_hybrid(99)) * S_h
    isf_q_120 = float(isf_quartic(120)) * S_q
    isf_h_120 = float(isf_hybrid(120)) * S_h

    print(f"  {r['name']:<14s} {r['model']:<8s} {tdd:5.1f}  {isf_1800:8.1f}  "
          f"{S_q:7.3f}  {S_h:8.3f}  "
          f"{isf_q_80:8.1f} {isf_h_80:8.1f}  "
          f"{isf_q_99:8.1f} {isf_h_99:8.1f}  "
          f"{isf_q_120:9.1f} {isf_h_120:9.1f}")

# ══════════════════════════════════════════════════════════════════════════════
# 5. ISF RATIO ANALYSIS — how do scaled formulas compare to loop's actual ISF?
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'═'*80}")
print("5. ISF RATIO: scaled formula ISF vs loop's actual ISF")
print("     ratio > 1 = formula gives HIGHER ISF = less aggressive")
print("     ratio < 1 = formula gives LOWER ISF = more aggressive")
print(f"{'═'*80}")

print(f"\n  {'Site':<14s} {'Model':<8s} {'TDD':>5s}  "
      f"{'Loop ISF':>8s}  "
      f"{'Q_ratio':>7s} {'DB_ratio':>8s} {'H_ratio':>7s}  "
      f"{'Q_r<105':>7s} {'DB_r<105':>8s} {'H_r<105':>7s}")
print(f"  {'─'*14} {'─'*8} {'─'*5}  {'─'*8}  "
      f"{'─'*7} {'─'*8} {'─'*7}  {'─'*7} {'─'*8} {'─'*7}")

for r in results:
    # We need to compute these from the actual data
    # The ISF ratio is ISF_formula / ISF_actual for the median sample
    # We can estimate from the median BG and median ISF
    tdd = r['tdd_median']
    S_q = (1800.0 / tdd) / 81.6
    S_h = (1800.0 / tdd) / 75.8

    bg_mean = r['bg_mean']
    isf_median = r['isf_actual_median']

    # Overall ratios at mean BG
    ratio_q = float(isf_quartic(bg_mean)) * S_q / isf_median
    ratio_db = float(isf_full_diabeloop(bg_mean)) * S_q / isf_median  # same S as quartic
    ratio_h = float(isf_hybrid(bg_mean)) * S_h / isf_median

    # At BG=85 (typical sub-105)
    bg_low = 85.0
    ratio_q_low = float(isf_quartic(bg_low)) * S_q / isf_median
    ratio_db_low = float(isf_full_diabeloop(bg_low)) * S_q / isf_median
    ratio_h_low = float(isf_hybrid(bg_low)) * S_h / isf_median

    print(f"  {r['name']:<14s} {r['model']:<8s} {tdd:5.1f}  "
          f"{isf_median:8.1f}  "
          f"{ratio_q:7.2f} {ratio_db:8.2f} {ratio_h:7.2f}  "
          f"{ratio_q_low:7.2f} {ratio_db_low:8.2f} {ratio_h_low:7.2f}")


# ══════════════════════════════════════════════════════════════════════════════
# 6. FULL DIABELOOP — WHY DOES IT WIN ON SOME SITES?
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'═'*80}")
print("6. FULL DIABELOOP — WHY DOES IT WIN?")
print(f"{'═'*80}")

# Sites where Full DB was the best formula
fulldb_wins = [r for r in results if r['overall']['full_db']['mae'] ==
               min(r['overall']['loop']['mae'], r['overall']['quartic']['mae'],
                   r['overall']['full_db']['mae'], r['overall']['hybrid']['mae'])]

print(f"\n  Full Diabeloop wins on {len(fulldb_wins)} sites: "
      f"{', '.join(r['name'] for r in fulldb_wins)}")

for r in fulldb_wins:
    tdd = r['tdd_median']
    S_q = (1800.0 / tdd) / 81.6
    isf_med = r['isf_actual_median']

    print(f"\n  {r['name']} ({r['model']}, TDD={tdd:.1f}, median ISF={isf_med:.1f})")
    print(f"    S_quartic={S_q:.3f}, 1800/TDD={1800/tdd:.1f}")

    # Compare Full DB vs Quartic in the <100 zone (where they differ)
    # Full DB uses quadratic ≤100, quartic uses quartic extended
    for bg_val in [70, 80, 85, 90, 95, 100, 105, 120]:
        q = float(isf_quartic(bg_val)) * S_q
        db_val = float(isf_full_diabeloop(bg_val)) * S_q
        h = float(isf_hybrid(bg_val)) * ((1800.0/tdd)/75.8)
        print(f"    BG={bg_val:3d}: Quartic={q:6.1f}  FullDB={db_val:6.1f}  Hybrid={h:6.1f}  LoopISF≈{isf_med:5.1f}  "
              f"Q/Loop={q/isf_med:.2f}  DB/Loop={db_val/isf_med:.2f}  H/Loop={h/isf_med:.2f}")

    # Band detail
    print(f"    Band breakdown:")
    for bname in ['<90', '90-105', '105-120', '120-150', '150+']:
        bl = r['bands']['loop'].get(bname, {})
        bq = r['bands']['quartic'].get(bname, {})
        bd = r['bands']['full_db'].get(bname, {})
        bh = r['bands']['hybrid'].get(bname, {})
        if bl:
            print(f"      {bname:<8s} n={bl['n']:4d}  Loop={bl['mae']:5.1f}({bl['bias']:+5.1f})  "
                  f"Quart={bq.get('mae',0):5.1f}({bq.get('bias',0):+5.1f})  "
                  f"FullDB={bd.get('mae',0):5.1f}({bd.get('bias',0):+5.1f})  "
                  f"Hybrid={bh.get('mae',0):5.1f}({bh.get('bias',0):+5.1f})")


# ══════════════════════════════════════════════════════════════════════════════
# 7. BELOW-90 PREDICTION DETAIL
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'═'*80}")
print("7. BELOW 90 — PREDICTED vs ACTUAL 2h BG")
print("   Higher predicted = formula predicts larger rise = more conservative")
print(f"{'═'*80}")

sites_with_sub90 = [r for r in results if r['sub90']['n'] >= 5]
print(f"\n  {'Site':<14s} {'Model':<8s} {'N':>4s}  {'Actual':>6s}  "
      f"{'Loop':>6s} {'Quart':>6s} {'FullDB':>6s} {'Hybrid':>6s}  "
      f"{'Q-Act':>5s} {'DB-Act':>6s} {'H-Act':>5s}")
print(f"  {'─'*14} {'─'*8} {'─'*4}  {'─'*6}  "
      f"{'─'*6} {'─'*6} {'─'*6} {'─'*6}  "
      f"{'─'*5} {'─'*6} {'─'*5}")

for r in sites_with_sub90:
    s = r['sub90']
    actual = s['actual_mean']
    print(f"  {r['name']:<14s} {r['model']:<8s} {s['n']:4d}  {actual:6.0f}  "
          f"{s['pred_loop']:6.0f} {s['pred_quartic']:6.0f} "
          f"{s['pred_full_db']:6.0f} {s['pred_hybrid']:6.0f}  "
          f"{s['pred_quartic']-actual:+5.0f} {s['pred_full_db']-actual:+6.0f} "
          f"{s['pred_hybrid']-actual:+5.0f}")


# ══════════════════════════════════════════════════════════════════════════════
# 8. KEY QUESTION: Is the hybrid's power-law tail giving ISF values
#    that are TOO HIGH relative to what Trio's sigmoid/log actually uses?
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'═'*80}")
print("8. HYBRID POWER-LAW TAIL vs LOOP ACTUAL ISF — below 105")
print("   If hybrid ISF >> loop ISF, the formula predicts a much larger")
print("   glucose rise than the loop actually predicted, potentially")
print("   making the counterfactual prediction too optimistic.")
print(f"{'═'*80}")

print(f"\n  At BG=80, the unscaled formulas give:")
print(f"    Quartic:  {float(isf_quartic(80)):6.1f}")
print(f"    Quadratic:{float(isf_full_diabeloop(80)):6.1f}")
print(f"    Hybrid:   {float(isf_hybrid(80)):6.1f}")
print(f"    (Hybrid power-law is {float(isf_hybrid(80))/float(isf_quartic(80)):.1f}× the quartic)")

print(f"\n  After TDD scaling, at BG=80:")
print(f"  {'Site':<14s} {'TDD':>5s}  {'Loop_ISF':>8s}  {'Q*S':>6s} {'DB*S':>6s} {'H*S':>6s}  "
      f"{'H/Loop':>6s} {'Q/Loop':>6s} {'DB/Loop':>7s}")
print(f"  {'─'*14} {'─'*5}  {'─'*8}  {'─'*6} {'─'*6} {'─'*6}  "
      f"{'─'*6} {'─'*6} {'─'*7}")

for r in results:
    tdd = r['tdd_median']
    S_q = (1800.0 / tdd) / 81.6
    S_h = (1800.0 / tdd) / 75.8
    isf_med = r['isf_actual_median']

    q80 = float(isf_quartic(80)) * S_q
    db80 = float(isf_full_diabeloop(80)) * S_q
    h80 = float(isf_hybrid(80)) * S_h

    print(f"  {r['name']:<14s} {tdd:5.1f}  {isf_med:8.1f}  "
          f"{q80:6.1f} {db80:6.1f} {h80:6.1f}  "
          f"{h80/isf_med:6.2f} {q80/isf_med:6.2f} {db80/isf_med:7.2f}")


print(f"\n{'═'*80}")
print("DONE")
print(f"{'═'*80}")
