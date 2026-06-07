#!/usr/bin/env python3
"""
DynamicISF v1 vs v2 comparison analysis

v1 formula (BoostPlugin):   ISF = 1800 / (BlendedTDD × 0.70 × ln(BG/82 + 1))
v2 formula (BoostV2Plugin): ISF = 2300 / (BlendedTDD² × 0.02 × ln(BG/82 + 1))

Key insight: ISF_v1 / ISF_v2 = BlendedTDD × (1800×0.02) / (2300×0.70)
                              = BlendedTDD × 0.02226
Crossover (equal ISF) at BlendedTDD ≈ 44.9 U/day.
Below that, v1 is more aggressive (lower ISF); above it, v2 is.
"""

import re
import math
import zipfile
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from datetime import datetime, timedelta

# ── Formula constants ─────────────────────────────────────────────────────────
C_V1      = 1800.0   # empirically confirmed from log back-calculation
C_V2      = 2300.0   # from formula string in logs
DIVISOR   = 82.0     # insulinDivisor (mg/dL)
ADJ_V1    = 0.70
SCALE_V2  = 0.02
TARGET_BG = 99.0     # normalTarget = 5.5 mmol/L


def isf_v1_formula(blended_tdd, bg):
    return C_V1 / (blended_tdd * ADJ_V1 * math.log(bg / DIVISOR + 1))


def isf_v2_formula(blended_tdd, bg):
    return C_V2 / (blended_tdd ** 2 * SCALE_V2 * math.log(bg / DIVISOR + 1))


# ── Log parsing ───────────────────────────────────────────────────────────────

def parse_folder(folder):
    records = []
    seen = set()

    for filename in sorted(os.listdir(folder)):
        if not filename.endswith('.zip'):
            continue
        m = re.search(r'AndroidAPS\._(\d{4}-\d{2}-\d{2})_', filename)
        if not m:
            continue
        file_date_str = m.group(1)

        zip_path = os.path.join(folder, filename)
        with zipfile.ZipFile(zip_path, 'r') as zf:
            for name in zf.namelist():
                with zf.open(name) as f:
                    content = f.read().decode('utf-8', errors='replace')
                parse_blocks(content, file_date_str, records, seen)

    return records


def parse_blocks(content, file_date_str, records, seen):
    lines = content.split('\n')
    prev_secs = None
    date_offset = 0

    for i, line in enumerate(lines):
        # Only match plain-text timestamped calculateBoostIsf lines (skip JSON re-embeds)
        m = re.match(
            r'^(\d{2}:\d{2}:\d{2})\.\d{3} .*calculateBoostIsf.*Boost (V2 )?ISF: TDD data:',
            line
        )
        if not m:
            continue

        time_str   = m.group(1)
        is_v2_block = m.group(2) is not None

        h, mi, s = int(time_str[0:2]), int(time_str[3:5]), int(time_str[6:8])
        cur_secs = h * 3600 + mi * 60 + s

        # Detect midnight rollover
        if prev_secs is not None and cur_secs < prev_secs - 3600:
            date_offset += 1
        prev_secs = cur_secs

        base_date = datetime.strptime(file_date_str, '%Y-%m-%d') + timedelta(days=date_offset)
        ts = base_date.replace(hour=h, minute=mi, second=s)

        block = {
            'timestamp':    ts,
            'formula_type': 'v2' if is_v2_block else 'v1',
        }

        for j in range(i + 1, min(i + 14, len(lines))):
            l = lines[j].strip()
            if not l or l.startswith('"') or l.startswith(','):
                continue  # skip JSON re-embeds

            bm = re.match(r'^Blended TDD=([\d.]+)$', l)
            if bm:
                block['blended_tdd'] = float(bm.group(1))

            fm = re.match(r'^Final TDD=([\d.]+) \(adj factor (\d+)%\)', l)
            if fm:
                block['final_tdd']  = float(fm.group(1))
                block['adj_factor'] = int(fm.group(2))

            im = re.match(r'^(?:TDD|V2) ISF at target: ([\d.]+) mg/dl/U \(profile was ([\d.]+)\)', l)
            if im:
                block['isf_at_target'] = float(im.group(1))
                block['profile_isf']   = float(im.group(2))

            vm = re.match(r'^Variable ISF at BG ([\d.]+): ([\d.]+)', l)
            if vm:
                block['bg']           = float(vm.group(1))
                block['variable_sens'] = float(vm.group(2))
                break  # end of block

        required = ('blended_tdd', 'final_tdd', 'isf_at_target', 'bg', 'variable_sens')
        if not all(k in block for k in required):
            continue
        # Skip records with no CGM reading (BG=0 → variable_sens = Long.MAX_VALUE)
        if block['bg'] <= 0 or block['variable_sens'] > 1e10:
            continue

        # Deduplicate on (minute-precision timestamp, blended_tdd, isf_at_target)
        key = (ts.strftime('%Y-%m-%d %H:%M'), block['blended_tdd'], block['isf_at_target'])
        if key in seen:
            continue
        seen.add(key)

        records.append(block)


# ── Main ──────────────────────────────────────────────────────────────────────

print("Parsing dynisfv1...")
v1_records = parse_folder('dynisfv1')
print(f"  {len(v1_records)} unique records")

print("Parsing dynisfv2...")
v2_records = parse_folder('dynisfv2')
print(f"  {len(v2_records)} unique records (includes transition period v1 blocks)")

df_v1 = pd.DataFrame(v1_records)
df_v2 = pd.DataFrame(v2_records)
df_v1['dataset'] = 'dynisfv1'
df_v2['dataset'] = 'dynisfv2'

df = pd.concat([df_v1, df_v2], ignore_index=True).sort_values('timestamp').reset_index(drop=True)
print(f"\nCombined: {len(df)} records  ({len(df[df.formula_type=='v1'])} v1-formula, {len(df[df.formula_type=='v2'])} v2-formula)")

# ── Compute normalised ISF values ─────────────────────────────────────────────
# For every data point (regardless of which formula generated it), compute what
# both formulas would give at that BlendedTDD and BG.
df['isf_v1_at_bg']     = df.apply(lambda r: isf_v1_formula(r.blended_tdd, r.bg), axis=1)
df['isf_v2_at_bg']     = df.apply(lambda r: isf_v2_formula(r.blended_tdd, r.bg), axis=1)
df['isf_v1_at_target'] = df.apply(lambda r: isf_v1_formula(r.blended_tdd, TARGET_BG), axis=1)
df['isf_v2_at_target'] = df.apply(lambda r: isf_v2_formula(r.blended_tdd, TARGET_BG), axis=1)

# Ratio is purely a function of BlendedTDD: ISF_v1/ISF_v2 = TDD × (C_V1×SCALE_V2)/(C_V2×ADJ_V1)
k = (C_V1 * SCALE_V2) / (C_V2 * ADJ_V1)
df['isf_ratio_v1_over_v2'] = df['blended_tdd'] * k
crossover_tdd = 1.0 / k
print(f"\nISF ratio (v1/v2) = BlendedTDD × {k:.5f}")
print(f"Crossover TDD (equal ISF): {crossover_tdd:.1f} U/day")
ratio_at_27 = 27 * k
print(f"At BlendedTDD=27: ISF_v1/ISF_v2 = {ratio_at_27:.2f}  → v1 ISF is {ratio_at_27*100:.0f}% of v2 ISF")
print(f"  → v1 is {(1/ratio_at_27 - 1)*100:.0f}% MORE aggressive than v2 at this TDD\n")

df.to_csv('dynisf_analysis.csv', index=False)
print("Saved: dynisf_analysis.csv")

# ── Separate actual v1 and v2 formula records for plotting ────────────────────
actual_v1 = df[df['formula_type'] == 'v1'].copy()
actual_v2 = df[df['formula_type'] == 'v2'].copy()

# ── Plots ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(16, 11))
fig.suptitle('DynamicISF v1 vs v2: Formula Comparison', fontsize=14, fontweight='bold')

BLUE  = '#1f77b4'
RED   = '#d62728'
ALPHA = 0.6

# ── Plot 1: ISF at target over time (both formulas, computed from actual data) ─
ax = axes[0, 0]
ax.plot(actual_v1['timestamp'], actual_v1['isf_v1_at_target'],
        color=BLUE, alpha=ALPHA, linewidth=0.8, label='v1 formula (actual)')
ax.plot(actual_v1['timestamp'], actual_v1['isf_v2_at_target'],
        color=RED, alpha=ALPHA, linewidth=0.8, linestyle='--', label='v2 formula (normalised onto v1 data)')

ax.plot(actual_v2['timestamp'], actual_v2['isf_v2_at_target'],
        color=RED, alpha=ALPHA, linewidth=0.8, label='v2 formula (actual)')
ax.plot(actual_v2['timestamp'], actual_v2['isf_v1_at_target'],
        color=BLUE, alpha=ALPHA, linewidth=0.8, linestyle='--', label='v1 formula (normalised onto v2 data)')

ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha='right', fontsize=7)
ax.set_ylabel('ISF at target (mg/dL/U)')
ax.set_title('ISF at Target BG (99 mg/dL) Over Time')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# ── Plot 2: ISF ratio vs BlendedTDD (theoretical line + actual scatter) ────────
ax = axes[0, 1]
tdd_range = np.linspace(df['blended_tdd'].min() * 0.9, max(df['blended_tdd'].max() * 1.1, 50), 200)
ratio_line = tdd_range * k
ax.plot(tdd_range, ratio_line, 'k-', linewidth=2, label=f'ISF_v1/ISF_v2 = TDD × {k:.4f}')
ax.axhline(1.0, color='gray', linestyle=':', linewidth=1)
ax.axvline(crossover_tdd, color='gray', linestyle=':', linewidth=1,
           label=f'Crossover = {crossover_tdd:.1f} U/day')

ax.scatter(actual_v1['blended_tdd'], actual_v1['isf_ratio_v1_over_v2'],
           color=BLUE, alpha=0.4, s=15, label='v1 dataset points')
ax.scatter(actual_v2['blended_tdd'], actual_v2['isf_ratio_v1_over_v2'],
           color=RED, alpha=0.4, s=15, label='v2 dataset points')

ax.set_xlabel('Blended TDD (U/day)')
ax.set_ylabel('ISF_v1 / ISF_v2')
ax.set_title('ISF Ratio vs Blended TDD\n(ratio depends only on TDD, not BG)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.text(0.05, 0.92, '← v2 more aggressive', transform=ax.transAxes, fontsize=8, color='darkred')
ax.text(0.05, 0.08, '← v1 more aggressive', transform=ax.transAxes, fontsize=8, color='steelblue')

# ── Plot 3: ISF vs BG scatter (both formulas, over full BG range) ─────────────
ax = axes[1, 0]
bg_range = np.linspace(40, 200, 300)

# Use median BlendedTDD from each dataset for representative curves
med_tdd_v1 = actual_v1['blended_tdd'].median() if len(actual_v1) else df['blended_tdd'].median()
med_tdd_v2 = actual_v2['blended_tdd'].median() if len(actual_v2) else df['blended_tdd'].median()

curve_v1 = [isf_v1_formula(med_tdd_v1, bg) for bg in bg_range]
curve_v2 = [isf_v2_formula(med_tdd_v2, bg) for bg in bg_range]

ax.plot(bg_range, curve_v1, color=BLUE, linewidth=2,
        label=f'v1 formula (TDD={med_tdd_v1:.1f})')
ax.plot(bg_range, curve_v2, color=RED, linewidth=2,
        label=f'v2 formula (TDD={med_tdd_v2:.1f})')

ax.scatter(actual_v1['bg'], actual_v1['isf_v1_at_bg'],
           color=BLUE, alpha=0.3, s=10)
ax.scatter(actual_v2['bg'], actual_v2['isf_v2_at_bg'],
           color=RED, alpha=0.3, s=10)

ax.axvline(TARGET_BG, color='gray', linestyle=':', linewidth=1, label=f'Target={TARGET_BG:.0f}')
ax.set_xlabel('BG (mg/dL)')
ax.set_ylabel('ISF (mg/dL/U)')
ax.set_title('ISF vs BG — Formula Shape Comparison\n(at median BlendedTDD for each dataset)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(40, 210)

# ── Plot 4: Distribution of ISF at target ─────────────────────────────────────
ax = axes[1, 1]
bins = np.linspace(
    df[['isf_v1_at_target', 'isf_v2_at_target']].min().min() * 0.9,
    df[['isf_v1_at_target', 'isf_v2_at_target']].max().max() * 1.05,
    50
)

ax.hist(df['isf_v1_at_target'], bins=bins, color=BLUE, alpha=0.6, label='v1 formula ISF at target')
ax.hist(df['isf_v2_at_target'], bins=bins, color=RED, alpha=0.6, label='v2 formula ISF at target')

v1_med = df['isf_v1_at_target'].median()
v2_med = df['isf_v2_at_target'].median()
ax.axvline(v1_med, color=BLUE, linewidth=2, linestyle='--', label=f'v1 median = {v1_med:.0f}')
ax.axvline(v2_med, color=RED,  linewidth=2, linestyle='--', label=f'v2 median = {v2_med:.0f}')

ax.set_xlabel('ISF at target (mg/dL/U)')
ax.set_ylabel('Count')
ax.set_title('Distribution of ISF at Target\n(computed for ALL data points using both formulas)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('dynisf_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: dynisf_comparison.png")

# ── Summary stats ──────────────────────────────────────────────────────────────
print("\n── Summary ─────────────────────────────────────────────────────────────")
print(f"{'':30s} {'v1 formula':>12} {'v2 formula':>12}")
print(f"{'ISF at target — median':30s} {df['isf_v1_at_target'].median():12.1f} {df['isf_v2_at_target'].median():12.1f}")
print(f"{'ISF at target — mean':30s} {df['isf_v1_at_target'].mean():12.1f} {df['isf_v2_at_target'].mean():12.1f}")
print(f"{'BlendedTDD — median (v1 data)':30s} {actual_v1['blended_tdd'].median():12.1f} {'':>12}")
print(f"{'BlendedTDD — median (v2 data)':30s} {'':>12} {actual_v2['blended_tdd'].median():12.1f}")

print(f"\nAt median TDD of each dataset:")
if len(actual_v1):
    t = actual_v1['blended_tdd'].median()
    print(f"  v1 data (TDD={t:.1f}): v1 ISF={isf_v1_formula(t, TARGET_BG):.1f}, "
          f"v2 would give={isf_v2_formula(t, TARGET_BG):.1f} "
          f"(v2 is {((isf_v2_formula(t, TARGET_BG)/isf_v1_formula(t, TARGET_BG))-1)*100:+.0f}% vs v1)")
if len(actual_v2):
    t = actual_v2['blended_tdd'].median()
    print(f"  v2 data (TDD={t:.1f}): v2 ISF={isf_v2_formula(t, TARGET_BG):.1f}, "
          f"v1 would give={isf_v1_formula(t, TARGET_BG):.1f} "
          f"(v2 is {((isf_v2_formula(t, TARGET_BG)/isf_v1_formula(t, TARGET_BG))-1)*100:+.0f}% vs v1)")
