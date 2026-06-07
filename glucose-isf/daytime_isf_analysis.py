#!/usr/bin/env python3
"""
Daytime Fasting ISF Analysis
=============================
Compares overnight vs daytime vs all-day ISF model performance
across 13 subjects (12 Trio + 1 AAPS), using only fasting periods
(COB=0, no recent bolus).

Produces charts and summary tables for the paper.
"""

import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

OUT_DIR = Path(__file__).parent
TRIO_CACHE = OUT_DIR / 'multisite_allday_cache.pkl'
BOOST_CACHE = OUT_DIR / 'boost_allday_cache.pkl'


def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

Q_REF = quartic(100)

RATIOS_POP = {76: 1.15, 100: 1.00, 130: 0.80, 170: 0.70}


def ratio_fn(g, ratios=None):
    if ratios is None:
        ratios = RATIOS_POP
    points = sorted(ratios.items())
    gs = [p[0] for p in points]; rs = [p[1] for p in points]
    if g <= gs[0]: return rs[0]
    if g >= gs[-1]: return rs[-1]
    for i in range(len(gs) - 1):
        if gs[i] <= g <= gs[i + 1]:
            t = (g - gs[i]) / (gs[i + 1] - gs[i])
            return rs[i] + t * (rs[i + 1] - rs[i])
    return 1.0


import math
def sigmoid_ratio(g, target=100):
    ln_ref = math.log(target / 120 + 1)
    ln_g = math.log(max(g, 40) / 120 + 1)
    return ln_ref / ln_g if ln_g > 0 else 1.0


# ── Load data ──────────────────────────────────────────────────────────────

with open(TRIO_CACHE, 'rb') as f:
    trio_sites = pickle.load(f)

with open(BOOST_CACHE, 'rb') as f:
    boost_cache = pickle.load(f)


def build_sites(trio_sites, boost_cache, period='allday'):
    sites = []
    for s in trio_sites:
        data = s.get(period)
        if data is None or data['n'] < 10:
            continue
        bg = data['bg']; isf = data['isf_actual']
        m100 = (bg >= 96) & (bg < 104)
        isf100 = np.median(isf[m100]) if m100.sum() >= 5 else np.nan
        sites.append({
            'name': s['name'], 'model': s['model'],
            'tdd': s['tdd_median'], 'n': data['n'],
            'bg': bg, 'isf_actual': isf,
            'pred_drop': data['pred_drop'],
            'actual_bg_end': data['actual_bg_end'],
            'pred_loop': data['pred_loop'],
            'hour': data.get('hour'),
            'isf_true': isf100,
            'isf_tdd': 1800 / s['tdd_median'],
        })

    boost_df = boost_cache.get(period)
    if boost_df is not None and len(boost_df) >= 10:
        bb = boost_df['bg'].values.astype(float)
        bi = boost_df['variable_sens'].values.astype(float)
        m100 = (bb >= 96) & (bb < 104)
        isf100 = np.median(bi[m100]) if m100.sum() >= 5 else np.nan
        tdd = boost_df['tdd_7day'].median()
        sites.append({
            'name': 'User-M', 'model': 'AAPS',
            'tdd': tdd, 'n': len(bb),
            'bg': bb, 'isf_actual': bi,
            'pred_drop': boost_df['pred_drop'].values.astype(float),
            'actual_bg_end': boost_df['actual_bg_end'].values.astype(float),
            'pred_loop': bb - boost_df['pred_drop'].values.astype(float),
            'hour': boost_df['hour'].values.astype(int),
            'isf_true': isf100, 'isf_tdd': 1800 / tdd,
        })
    return sites


def counterfactual_mae(bg, pred_drop, actual_end, isf_actual, isf_model):
    pred = bg - pred_drop * (isf_model / isf_actual)
    return np.mean(np.abs(pred - actual_end))


def counterfactual_bias(bg, pred_drop, actual_end, isf_actual, isf_model):
    pred = bg - pred_drop * (isf_model / isf_actual)
    return np.mean(pred - actual_end)


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: Sample counts and demographics
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 90)
print("ANALYSIS 1: Dataset overview — overnight vs daytime fasting samples")
print("=" * 90)

for period in ['allday', 'overnight', 'daytime']:
    sites = build_sites(trio_sites, boost_cache, period)
    total = sum(s['n'] for s in sites)
    print(f"\n  {period.upper()}: {len(sites)} sites, {total:,} samples")
    for s in sites:
        isf_str = f"{s['isf_true']:.0f}" if not np.isnan(s['isf_true']) else 'N/A'
        print(f"    {s['name']:8s}  n={s['n']:5d}  {s['model']:7s}  "
              f"TDD={s['tdd']:5.1f}  ISF@100={isf_str:>5s}")


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: ISF profile stability — does ISF@100 change between periods?
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("ANALYSIS 2: ISF at target glucose — overnight vs daytime")
print("=" * 90)

on_sites = build_sites(trio_sites, boost_cache, 'overnight')
dt_sites = build_sites(trio_sites, boost_cache, 'daytime')

on_isf = {s['name']: s['isf_true'] for s in on_sites}
dt_isf = {s['name']: s['isf_true'] for s in dt_sites}

print(f"\n  {'Site':8s} {'ON ISF@100':>10s} {'DT ISF@100':>10s} {'Change':>8s} {'%':>6s}")
print("  " + "-" * 50)
isf_changes = []
for name in sorted(set(on_isf.keys()) & set(dt_isf.keys())):
    on_v = on_isf[name]; dt_v = dt_isf[name]
    if np.isnan(on_v) or np.isnan(dt_v): continue
    diff = dt_v - on_v
    pct = (diff / on_v) * 100
    isf_changes.append(pct)
    print(f"  {name:8s} {on_v:10.0f} {dt_v:10.0f} {diff:+8.0f} {pct:+5.0f}%")

if isf_changes:
    print(f"\n  Mean ISF change overnight→daytime: {np.mean(isf_changes):+.1f}%")
    print(f"  Range: {min(isf_changes):+.0f}% to {max(isf_changes):+.0f}%")


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: Full model comparison — overnight vs daytime vs all-day
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("ANALYSIS 3: Model comparison across time periods")
print("=" * 90)

MODEL_DEFS = {
    'Loop (tuned)': lambda s: s['isf_actual'],  # identity — uses loop's own ISF
    'Quartic + 1800/TDD': lambda s: np.array([s['isf_tdd'] * quartic(g) / Q_REF for g in s['bg']]),
    'Profile + pop ratios': lambda s: np.array([s['isf_true'] * ratio_fn(g) for g in s['bg']]),
    'Profile + quartic': lambda s: np.array([s['isf_true'] * quartic(g) / Q_REF for g in s['bg']]),
    'Profile + sigmoid': lambda s: np.array([s['isf_true'] * sigmoid_ratio(g) for g in s['bg']]),
    'Profile + flat': lambda s: np.full(len(s['bg']), s['isf_true']),
    '1800/TDD + flat': lambda s: np.full(len(s['bg']), s['isf_tdd']),
}

all_period_results = {}

for period in ['overnight', 'daytime', 'allday']:
    sites = build_sites(trio_sites, boost_cache, period)
    valid = [s for s in sites if not np.isnan(s['isf_true'])]
    weights = np.array([s['n'] for s in valid])

    period_results = {}
    for model_name, model_fn in MODEL_DEFS.items():
        site_maes = []
        site_biases = []
        for s in valid:
            if model_name == 'Loop (tuned)':
                mae = np.mean(np.abs(s['pred_loop'] - s['actual_bg_end']))
                bias = np.mean(s['pred_loop'] - s['actual_bg_end'])
            else:
                isf_model = model_fn(s)
                mae = counterfactual_mae(s['bg'], s['pred_drop'], s['actual_bg_end'],
                                         s['isf_actual'], isf_model)
                bias = counterfactual_bias(s['bg'], s['pred_drop'], s['actual_bg_end'],
                                           s['isf_actual'], isf_model)
            site_maes.append(mae)
            site_biases.append(bias)

        wmae = np.average(site_maes, weights=weights)
        wbias = np.average(site_biases, weights=weights)
        period_results[model_name] = {
            'wmae': wmae, 'wbias': wbias,
            'per_site_mae': {s['name']: m for s, m in zip(valid, site_maes)},
            'per_site_bias': {s['name']: b for s, b in zip(valid, site_biases)},
        }

    all_period_results[period] = period_results

# Print comparison table
print(f"\n  {'Model':25s} {'Overnight':>10s} {'Daytime':>10s} {'All-Day':>10s} {'Δ Day-Night':>12s}")
print("  " + "-" * 70)
for model_name in MODEL_DEFS:
    on_mae = all_period_results['overnight'][model_name]['wmae']
    dt_mae = all_period_results['daytime'][model_name]['wmae']
    ad_mae = all_period_results['allday'][model_name]['wmae']
    delta = dt_mae - on_mae
    print(f"  {model_name:25s} {on_mae:10.1f} {dt_mae:10.1f} {ad_mae:10.1f} {delta:+12.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 4: Per-site overnight vs daytime — Loop MAE
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("ANALYSIS 4: Per-site Loop MAE — overnight vs daytime")
print("=" * 90)

on_loop = all_period_results['overnight']['Loop (tuned)']['per_site_mae']
dt_loop = all_period_results['daytime']['Loop (tuned)']['per_site_mae']

print(f"\n  {'Site':8s} {'ON MAE':>8s} {'DT MAE':>8s} {'Δ':>8s} {'ON n':>6s} {'DT n':>6s}")
print("  " + "-" * 45)
on_sites_list = build_sites(trio_sites, boost_cache, 'overnight')
dt_sites_list = build_sites(trio_sites, boost_cache, 'daytime')
on_n = {s['name']: s['n'] for s in on_sites_list}
dt_n = {s['name']: s['n'] for s in dt_sites_list}

for name in sorted(set(on_loop.keys()) & set(dt_loop.keys())):
    delta = dt_loop[name] - on_loop[name]
    print(f"  {name:8s} {on_loop[name]:8.1f} {dt_loop[name]:8.1f} {delta:+8.1f} "
          f"{on_n.get(name, 0):6d} {dt_n.get(name, 0):6d}")


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 5: Which model is most ROBUST across time periods?
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("ANALYSIS 5: Model robustness — smallest daytime degradation")
print("=" * 90)

print(f"\n  {'Model':25s} {'ON MAE':>8s} {'DT MAE':>8s} {'Δ':>8s} {'Δ%':>8s}")
print("  " + "-" * 60)
for model_name in MODEL_DEFS:
    on = all_period_results['overnight'][model_name]['wmae']
    dt = all_period_results['daytime'][model_name]['wmae']
    delta = dt - on
    pct = (delta / on) * 100
    print(f"  {model_name:25s} {on:8.1f} {dt:8.1f} {delta:+8.1f} {pct:+7.1f}%")

# Rank by smallest % degradation
ranked = sorted(MODEL_DEFS.keys(),
                key=lambda m: (all_period_results['daytime'][m]['wmae'] -
                               all_period_results['overnight'][m]['wmae']) /
                              all_period_results['overnight'][m]['wmae'])
print(f"\n  Most robust (smallest % degradation):")
for i, m in enumerate(ranked):
    on = all_period_results['overnight'][m]['wmae']
    dt = all_period_results['daytime'][m]['wmae']
    pct = ((dt - on) / on) * 100
    print(f"    {i+1}. {m} ({pct:+.1f}%)")


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 6: Bias direction — do models systematically over/under predict?
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("ANALYSIS 6: Prediction bias (positive = over-predict BG = under-dose)")
print("=" * 90)

print(f"\n  {'Model':25s} {'ON bias':>8s} {'DT bias':>8s} {'Δ':>8s}")
print("  " + "-" * 55)
for model_name in MODEL_DEFS:
    on = all_period_results['overnight'][model_name]['wbias']
    dt = all_period_results['daytime'][model_name]['wbias']
    print(f"  {model_name:25s} {on:+8.1f} {dt:+8.1f} {dt-on:+8.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# CHARTS
# ══════════════════════════════════════════════════════════════════════════════

fig = plt.figure(figsize=(16, 14))
fig.suptitle("Daytime vs Overnight Fasting: Dynamic ISF Model Performance\n"
             "(All samples: COB = 0, no recent bolus, BG 72–200 mg/dL)",
             fontsize=13, fontweight='bold')

# Chart 1: Model comparison bars — overnight vs daytime
ax1 = fig.add_subplot(2, 2, 1)
models_short = ['Loop', 'Q+TDD', 'Prof+pop', 'Prof+Q', 'Prof+Sig', 'Prof+Flat', 'TDD+Flat']
model_keys = list(MODEL_DEFS.keys())
x = np.arange(len(models_short))
width = 0.35
on_vals = [all_period_results['overnight'][m]['wmae'] for m in model_keys]
dt_vals = [all_period_results['daytime'][m]['wmae'] for m in model_keys]
ax1.bar(x - width/2, on_vals, width, label='Overnight', color='tab:blue', alpha=0.8)
ax1.bar(x + width/2, dt_vals, width, label='Daytime', color='tab:orange', alpha=0.8)
ax1.set_xticks(x)
ax1.set_xticklabels(models_short, rotation=35, ha='right', fontsize=8)
ax1.set_ylabel('Weighted Mean MAE (mg/dL)')
ax1.set_title('A. Model Performance: Overnight vs Daytime')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3, axis='y')
ax1.set_ylim(10, 35)

# Chart 2: Per-site Loop MAE
ax2 = fig.add_subplot(2, 2, 2)
common = sorted(set(on_loop.keys()) & set(dt_loop.keys()))
x2 = np.arange(len(common))
on_v2 = [on_loop[n] for n in common]
dt_v2 = [dt_loop[n] for n in common]
ax2.bar(x2 - 0.15, on_v2, 0.3, label='Overnight', color='tab:blue', alpha=0.8)
ax2.bar(x2 + 0.15, dt_v2, 0.3, label='Daytime', color='tab:orange', alpha=0.8)
ax2.set_xticks(x2)
ax2.set_xticklabels(common, rotation=45, ha='right', fontsize=7)
ax2.set_ylabel('Loop MAE (mg/dL)')
ax2.set_title('B. Per-Site Loop MAE: Overnight vs Daytime')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3, axis='y')

# Chart 3: Hour distribution
ax3 = fig.add_subplot(2, 2, 3)
allday_sites = build_sites(trio_sites, boost_cache, 'allday')
all_hours = np.concatenate([s['hour'] for s in allday_sites if s.get('hour') is not None])
hour_counts = np.bincount(all_hours.astype(int), minlength=24)
colors_h = ['tab:blue' if h < 8 else 'tab:orange' for h in range(24)]
ax3.bar(range(24), hour_counts, color=colors_h, alpha=0.7)
ax3.axvline(7.5, color='red', linestyle='--', alpha=0.7, linewidth=1.5, label='Overnight/daytime split')
ax3.set_xlabel('Hour of day')
ax3.set_ylabel('Fasting samples (COB=0)')
ax3.set_title('C. Fasting Sample Distribution by Hour')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3, axis='y')

# Chart 4: Robustness — % degradation overnight → daytime
ax4 = fig.add_subplot(2, 2, 4)
pct_degs = []
for m in model_keys:
    on = all_period_results['overnight'][m]['wmae']
    dt = all_period_results['daytime'][m]['wmae']
    pct_degs.append(((dt - on) / on) * 100)
bar_colors = ['tab:green' if p < 15 else ('tab:orange' if p < 25 else 'tab:red') for p in pct_degs]
bars4 = ax4.bar(range(len(models_short)), pct_degs, color=bar_colors, alpha=0.7)
ax4.set_xticks(range(len(models_short)))
ax4.set_xticklabels(models_short, rotation=35, ha='right', fontsize=8)
ax4.set_ylabel('MAE degradation overnight → daytime (%)')
ax4.set_title('D. Model Robustness: Daytime Degradation')
for bar, val in zip(bars4, pct_degs):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{val:.0f}%', ha='center', fontsize=8)
ax4.grid(True, alpha=0.3, axis='y')
ax4.axhline(0, color='black', linewidth=0.5)

plt.tight_layout()
plt.savefig(OUT_DIR / 'daytime_isf_analysis.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR / 'daytime_isf_analysis.png'}")


# ── Second figure: ISF@100 shift and per-site detail ──────────────────────

fig2, axes2 = plt.subplots(1, 2, figsize=(14, 5))
fig2.suptitle("ISF Profile Stability: Overnight vs Daytime", fontsize=12, fontweight='bold')

# ISF@100 comparison
ax = axes2[0]
common_isf = sorted(set(on_isf.keys()) & set(dt_isf.keys()))
common_isf = [n for n in common_isf if not np.isnan(on_isf[n]) and not np.isnan(dt_isf[n])]
x5 = np.arange(len(common_isf))
on_isf_vals = [on_isf[n] for n in common_isf]
dt_isf_vals = [dt_isf[n] for n in common_isf]
ax.bar(x5 - 0.15, on_isf_vals, 0.3, label='Overnight', color='tab:blue', alpha=0.8)
ax.bar(x5 + 0.15, dt_isf_vals, 0.3, label='Daytime', color='tab:orange', alpha=0.8)
ax.set_xticks(x5)
ax.set_xticklabels(common_isf, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('ISF at 100 mg/dL (mg/dL/U)')
ax.set_title('ISF@100: Overnight vs Daytime')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

# Best model per period
ax = axes2[1]
best_models_on = []
best_models_dt = []
for name in common:
    best_on = min(MODEL_DEFS.keys(),
                  key=lambda m: all_period_results['overnight'][m]['per_site_mae'].get(name, 999))
    best_dt = min(MODEL_DEFS.keys(),
                  key=lambda m: all_period_results['daytime'][m]['per_site_mae'].get(name, 999))
    best_on_mae = all_period_results['overnight'][best_on]['per_site_mae'].get(name, np.nan)
    best_dt_mae = all_period_results['daytime'][best_dt]['per_site_mae'].get(name, np.nan)
    best_models_on.append((name, best_on, best_on_mae))
    best_models_dt.append((name, best_dt, best_dt_mae))

# Show as text table in the chart
text_lines = [f"{'Site':8s}  {'Best overnight':22s} {'MAE':>5s}  {'Best daytime':22s} {'MAE':>5s}"]
text_lines.append("-" * 70)
for (n1, m1, mae1), (n2, m2, mae2) in zip(best_models_on, best_models_dt):
    m1_short = m1[:20]
    m2_short = m2[:20]
    text_lines.append(f"{n1:8s}  {m1_short:22s} {mae1:5.1f}  {m2_short:22s} {mae2:5.1f}")
ax.text(0.02, 0.95, '\n'.join(text_lines), transform=ax.transAxes,
        fontsize=7, verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax.set_title('Best Model Per Site by Period')
ax.axis('off')

plt.tight_layout()
plt.savefig(OUT_DIR / 'daytime_isf_profile.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'daytime_isf_profile.png'}")

print("\nDONE")
