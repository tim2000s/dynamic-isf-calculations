"""
BG Scaling Analysis for the 7D-TDD DynamicISF Formula
======================================================
The TDD_effective analysis found that the sensitivity ratio is
BG-dependent (Kruskal-Wallis p<0.0001):
  - <90 mg/dL:   ratio ≈ 0.95 (formula nearly correct)
  - 90-120:      ratio ≈ 0.78 (over-aggressive)
  - 120-150:     ratio ≈ 0.75 (more over-aggressive)

This means the ln(BG/D+1) scaling drops ISF too steeply as BG rises.
The formula delivers too much insulin at elevated BG.

This script:
  1. Plots the empirical ISF_eff vs BG relationship
  2. Fits several candidate BG scaling functions
  3. Compares counterfactual prediction accuracy for each
  4. Determines optimal BG scaling parameters

Candidate scaling functions (all anchored at target BG = 1700/TDD):
  current:     ISF = K / ln(BG/D+1)                   D=82
  power-log:   ISF = K / ln(BG/D+1)^p                 fit p
  alt-D:       ISF = K' / ln(BG/D'+1)                 fit D'
  power-law:   ISF = ISF_t × (target/BG)^k            fit k
  linear:      ISF = ISF_t × max(0.3, 1 - m×(BG-target))  fit m
  poly-log:    ISF = K / (a × ln(BG/D+1) + b × ln(BG/D+1)²)  fit a,b

Outputs
-------
  ns_bg_scaling_results.png
  ns_bg_scaling_summary.txt
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import minimize_scalar, minimize
from scipy import stats as sp_stats

warnings.filterwarnings('ignore')

D      = 82.0
TARGET = 99.0
LN_TARGET = np.log(TARGET / D + 1)

MIN_BG_DROP = 3.0
BG_LO, BG_HI = 72, 200

# ── 1. Load and filter ───────────────────────────────────────────────────────
print("Loading ns_backtest_overnight.csv ...")
df = pd.read_csv('ns_backtest_overnight.csv')
df['ts'] = pd.to_datetime(df['ts'], format='ISO8601')
df['date'] = pd.to_datetime(df['date']).dt.date
df['hour'] = df['ts'].dt.hour
print(f"  Total overnight rows: {len(df):,}")

v = df.dropna(subset=['pred_iob_24', 'actual_bg_2h', 'isf_v1', 'tdd_7day']).copy()
v = v[(v['bg'] >= BG_LO) & (v['bg'] <= BG_HI)]

v['bg_drop_pred']  = v['bg'] - v['pred_iob_24']
v['bg_drop_actual'] = v['bg'] - v['actual_bg_2h']

v = v[v['bg_drop_pred'].abs() >= MIN_BG_DROP]
v['ratio'] = v['bg_drop_actual'] / v['bg_drop_pred']
v = v[(v['ratio'] > 0) & (v['ratio'] < 5)]
v = v[~((v['bg_drop_pred'] > 0) & (v['bg_drop_actual'] < -9))]

# ISF that would have zeroed the prediction error
v['isf_eff'] = v['isf_v1'] * v['ratio']
v = v[(v['isf_eff'] > 10) & (v['isf_eff'] < 800)]

print(f"  Valid samples: {len(v):,}")

# ── 2. Empirical ISF vs BG ───────────────────────────────────────────────────
print("\n" + "═" * 72)
print("EMPIRICAL ISF vs BG RELATIONSHIP")
print("═" * 72)

# Bin by BG and compute median ISF_eff per bin
v['bg_bin'] = pd.cut(v['bg'], bins=np.arange(70, 205, 5))
bg_profile = v.groupby('bg_bin', observed=True).agg(
    bg_mid=('bg', 'median'),
    n=('isf_eff', 'count'),
    isf_eff_med=('isf_eff', 'median'),
    isf_eff_q25=('isf_eff', lambda x: x.quantile(0.25)),
    isf_eff_q75=('isf_eff', lambda x: x.quantile(0.75)),
    isf_7dtdd_med=('isf_7dtdd', 'median'),
    isf_v1_med=('isf_v1', 'median'),
    tdd_7day_med=('tdd_7day', 'median'),
    ratio_med=('ratio', 'median'),
).dropna()

print(f"\n  {'BG':>5s}  {'n':>4s}  {'ISF_eff':>7s}  {'ISF_7dtdd':>9s}  {'ISF_v1':>6s}  "
      f"{'ratio':>6s}  {'err%':>6s}")
print("  " + "─" * 55)
for _, row in bg_profile.iterrows():
    err_pct = (row['isf_7dtdd_med'] - row['isf_eff_med']) / row['isf_eff_med'] * 100
    print(f"  {row['bg_mid']:5.0f}  {row['n']:4.0f}  {row['isf_eff_med']:7.1f}  "
          f"{row['isf_7dtdd_med']:9.1f}  {row['isf_v1_med']:6.1f}  "
          f"{row['ratio_med']:6.3f}  {err_pct:+6.1f}%")

# ── 3. Define candidate scaling functions ─────────────────────────────────────

# All functions return ISF given BG and parameters
# They are anchored so that at TARGET, ISF = 1700/TDD_7day

def isf_current(bg, tdd):
    """Current: K / ln(BG/D+1), D=82"""
    return (1700.0 / tdd) * LN_TARGET / np.log(bg / D + 1)

def isf_power_log(bg, tdd, p):
    """Power-log: K / ln(BG/D+1)^p"""
    ln_bg = np.log(bg / D + 1)
    return (1700.0 / tdd) * (LN_TARGET ** p) / (ln_bg ** p)

def isf_alt_D(bg, tdd, d_new):
    """Alternative D: K' / ln(BG/D'+1)"""
    ln_t = np.log(TARGET / d_new + 1)
    ln_b = np.log(bg / d_new + 1)
    return (1700.0 / tdd) * ln_t / ln_b

def isf_power_law(bg, tdd, k):
    """Power-law: ISF_target × (target/BG)^k"""
    isf_t = 1700.0 / tdd
    return isf_t * (TARGET / bg) ** k

def isf_linear(bg, tdd, m):
    """Linear: ISF_target × max(0.3, 1 - m×(BG-target))"""
    isf_t = 1700.0 / tdd
    return isf_t * np.maximum(0.3, 1.0 - m * (bg - TARGET))

def isf_sqrt_log(bg, tdd):
    """Sqrt-log: K / sqrt(ln(BG/D+1))  — equivalent to power-log p=0.5"""
    ln_bg = np.log(bg / D + 1)
    return (1700.0 / tdd) * np.sqrt(LN_TARGET) / np.sqrt(ln_bg)

# ── 4. Fit parameters by minimising MAE of counterfactual predictions ─────────

bg_vals  = v['bg'].values
tdd_vals = v['tdd_7day'].values
isf_v1   = v['isf_v1'].values
bg_drop  = v['bg_drop_pred'].values
actual_2h = v['actual_bg_2h'].values

def pred_error_mae(isf_candidate):
    """MAE of counterfactual prediction using candidate ISF"""
    pred = bg_vals - bg_drop * (isf_candidate / isf_v1)
    err = actual_2h - pred
    return np.abs(err).mean()

def pred_error_bias(isf_candidate):
    """Bias of counterfactual prediction"""
    pred = bg_vals - bg_drop * (isf_candidate / isf_v1)
    return (actual_2h - pred).mean()

def pred_errors(isf_candidate):
    """Full error array"""
    pred = bg_vals - bg_drop * (isf_candidate / isf_v1)
    return actual_2h - pred

print("\n" + "═" * 72)
print("FITTING BG SCALING FUNCTIONS")
print("═" * 72)

results = {}

# Current (no fitting)
isf_curr = isf_current(bg_vals, tdd_vals)
e = pred_errors(isf_curr)
results['current'] = {
    'label': f'Current: ln(BG/{D:.0f}+1)',
    'params': f'D={D:.0f}',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd: isf_current(bg, tdd),
}
print(f"\n  Current (D={D:.0f}):  MAE={results['current']['mae']:.2f}  bias={results['current']['bias']:+.2f}")

# Power-log: fit p
print("  Fitting power-log p ...")
def obj_p(p):
    isf = isf_power_log(bg_vals, tdd_vals, p)
    return pred_error_mae(isf)

res_p = minimize_scalar(obj_p, bounds=(0.1, 2.0), method='bounded')
p_best = res_p.x
isf_pl = isf_power_log(bg_vals, tdd_vals, p_best)
e = pred_errors(isf_pl)
results['power_log'] = {
    'label': f'Power-log: ln(BG/D+1)^{p_best:.3f}',
    'params': f'p={p_best:.3f}',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd, p=p_best: isf_power_log(bg, tdd, p),
}
print(f"    p={p_best:.3f}  MAE={results['power_log']['mae']:.2f}  bias={results['power_log']['bias']:+.2f}")

# Alt-D: fit D'
print("  Fitting alternative D ...")
def obj_d(d_new):
    isf = isf_alt_D(bg_vals, tdd_vals, d_new)
    return pred_error_mae(isf)

res_d = minimize_scalar(obj_d, bounds=(20, 500), method='bounded')
d_best = res_d.x
isf_ad = isf_alt_D(bg_vals, tdd_vals, d_best)
e = pred_errors(isf_ad)
results['alt_D'] = {
    'label': f'Alt divisor: ln(BG/{d_best:.0f}+1)',
    'params': f'D={d_best:.1f}',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd, d=d_best: isf_alt_D(bg, tdd, d),
}
print(f"    D'={d_best:.1f}  MAE={results['alt_D']['mae']:.2f}  bias={results['alt_D']['bias']:+.2f}")

# Power-law: fit k
print("  Fitting power-law k ...")
def obj_k(k):
    isf = isf_power_law(bg_vals, tdd_vals, k)
    return pred_error_mae(isf)

res_k = minimize_scalar(obj_k, bounds=(0.01, 3.0), method='bounded')
k_best = res_k.x
isf_pk = isf_power_law(bg_vals, tdd_vals, k_best)
e = pred_errors(isf_pk)
results['power_law'] = {
    'label': f'Power-law: (target/BG)^{k_best:.3f}',
    'params': f'k={k_best:.3f}',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd, k=k_best: isf_power_law(bg, tdd, k),
}
print(f"    k={k_best:.3f}  MAE={results['power_law']['mae']:.2f}  bias={results['power_law']['bias']:+.2f}")

# Linear: fit m
print("  Fitting linear slope m ...")
def obj_m(m):
    isf = isf_linear(bg_vals, tdd_vals, m)
    return pred_error_mae(isf)

res_m = minimize_scalar(obj_m, bounds=(0.0001, 0.02), method='bounded')
m_best = res_m.x
isf_lin = isf_linear(bg_vals, tdd_vals, m_best)
e = pred_errors(isf_lin)
results['linear'] = {
    'label': f'Linear: 1 - {m_best:.4f}×(BG-target)',
    'params': f'm={m_best:.5f}',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd, m=m_best: isf_linear(bg, tdd, m),
}
print(f"    m={m_best:.5f}  MAE={results['linear']['mae']:.2f}  bias={results['linear']['bias']:+.2f}")

# Sqrt-log (power-log with p=0.5)
isf_sq = isf_sqrt_log(bg_vals, tdd_vals)
e = pred_errors(isf_sq)
results['sqrt_log'] = {
    'label': 'Sqrt-log: sqrt(ln(BG/D+1))',
    'params': 'p=0.500',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd: isf_sqrt_log(bg, tdd),
}
print(f"  Sqrt-log (p=0.5):  MAE={results['sqrt_log']['mae']:.2f}  bias={results['sqrt_log']['bias']:+.2f}")

# Flat (no BG scaling — baseline)
isf_flat_vals = 1700.0 / tdd_vals
e = pred_errors(isf_flat_vals)
results['flat'] = {
    'label': 'Flat: 1700/TDD (no BG scaling)',
    'params': 'none',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd: np.full_like(bg, 1.0) * 1700.0 / tdd,
}
print(f"  Flat (no scaling): MAE={results['flat']['mae']:.2f}  bias={results['flat']['bias']:+.2f}")

# ── 5. Joint fit: p and D simultaneously ─────────────────────────────────────
print("\n  Fitting joint (p, D) ...")
def obj_pd(params):
    p, d_val = params
    ln_t = np.log(TARGET / d_val + 1)
    ln_b = np.log(bg_vals / d_val + 1)
    isf = (1700.0 / tdd_vals) * (ln_t ** p) / (ln_b ** p)
    pred = bg_vals - bg_drop * (isf / isf_v1)
    return np.abs(actual_2h - pred).mean()

res_pd = minimize(obj_pd, x0=[0.5, 82], bounds=[(0.05, 2.0), (20, 500)], method='L-BFGS-B')
p_j, d_j = res_pd.x
ln_t_j = np.log(TARGET / d_j + 1)
ln_b_j = np.log(bg_vals / d_j + 1)
isf_joint = (1700.0 / tdd_vals) * (ln_t_j ** p_j) / (ln_b_j ** p_j)
e = pred_errors(isf_joint)
results['joint'] = {
    'label': f'Joint: ln(BG/{d_j:.0f}+1)^{p_j:.3f}',
    'params': f'p={p_j:.3f}, D={d_j:.1f}',
    'mae': np.abs(e).mean(), 'bias': e.mean(), 'rmse': np.sqrt((e**2).mean()),
    'w18': (np.abs(e) <= 18).mean() * 100,
    'isf_fn': lambda bg, tdd, p=p_j, d=d_j: (1700.0/tdd) * (np.log(TARGET/d+1)**p) / (np.log(bg/d+1)**p),
}
print(f"    p={p_j:.3f}, D={d_j:.1f}  MAE={results['joint']['mae']:.2f}  bias={results['joint']['bias']:+.2f}")

# ── 6. Summary table ─────────────────────────────────────────────────────────
print("\n" + "═" * 72)
print("BG SCALING COMPARISON (+2h)")
print("═" * 72)

# Sort by MAE
ranked = sorted(results.items(), key=lambda x: x[1]['mae'])

print(f"\n  {'#':>2s}  {'Scaling':.<45s}  {'MAE':>6s}  {'Bias':>6s}  {'RMSE':>6s}  "
      f"{'±1mmol':>7s}  {'Params':>20s}")
print("  " + "─" * 100)

for i, (name, r) in enumerate(ranked):
    marker = " ←" if name == 'current' else ""
    print(f"  {i+1:2d}  {r['label']:.<45s}  {r['mae']:6.2f}  {r['bias']:+6.2f}  "
          f"{r['rmse']:6.2f}  {r['w18']:6.1f}%  {r['params']:>20s}{marker}")

# ── 7. BG-band breakdown for top formulas ─────────────────────────────────────
print(f"\n── BG-Band MAE Breakdown " + "─" * 50)

v['bg_band'] = pd.cut(v['bg'], bins=[0, 90, 105, 120, 150, 300],
                       labels=['<90', '90-105', '105-120', '120-150', '>150'])

top_formulas = ['current', 'power_log', 'alt_D', 'joint', 'power_law']

# Pre-compute ISF and errors for top formulas
for fname in top_formulas:
    r = results[fname]
    isf_vals = r['isf_fn'](bg_vals, tdd_vals)
    pred = bg_vals - bg_drop * (isf_vals / isf_v1)
    v[f'err_{fname}'] = actual_2h - pred

print(f"\n  {'BG Band':>10s}  {'n':>5s}  │", end='')
for fname in top_formulas:
    print(f"  {fname:>10s}", end='')
print()
print("  " + "─" * (22 + 12 * len(top_formulas)))

for band in ['<90', '90-105', '105-120', '120-150', '>150']:
    mask = v['bg_band'] == band
    n = mask.sum()
    if n < 5:
        continue
    print(f"  {band:>10s}  {n:5d}  │", end='')
    for fname in top_formulas:
        mae = v.loc[mask, f'err_{fname}'].abs().mean()
        print(f"  {mae:10.1f}", end='')
    print()

# Bias by band
print(f"\n  {'BG Band':>10s}  {'n':>5s}  │", end='')
for fname in top_formulas:
    print(f"  {fname:>10s}", end='')
print("   (bias)")
print("  " + "─" * (22 + 12 * len(top_formulas)))

for band in ['<90', '90-105', '105-120', '120-150', '>150']:
    mask = v['bg_band'] == band
    n = mask.sum()
    if n < 5:
        continue
    print(f"  {band:>10s}  {n:5d}  │", end='')
    for fname in top_formulas:
        bias = v.loc[mask, f'err_{fname}'].mean()
        print(f"  {bias:+10.1f}", end='')
    print()

# ── 8. What does the optimal scaling look like? ──────────────────────────────
print(f"\n── ISF Curves at Median TDD = {tdd_vals.mean():.1f} " + "─" * 35)

bg_range = np.arange(70, 201, 5)
tdd_med = np.median(tdd_vals)

print(f"\n  {'BG':>5s}", end='')
for fname in ['current', 'power_log', 'alt_D', 'joint', 'power_law', 'flat']:
    print(f"  {fname:>10s}", end='')
print()
print("  " + "─" * (7 + 12 * 6))

for bg in [75, 85, 95, 99, 105, 115, 125, 140, 160, 180, 200]:
    print(f"  {bg:5d}", end='')
    for fname in ['current', 'power_log', 'alt_D', 'joint', 'power_law', 'flat']:
        isf = results[fname]['isf_fn'](np.array([float(bg)]), np.array([tdd_med]))[0]
        print(f"  {isf:10.1f}", end='')
    print()

# ── 9. Ratio by BG — does the best formula fix BG-dependence? ────────────────
print(f"\n── Ratio by BG Band: Current vs Best " + "─" * 35)

best_name = ranked[0][0]
best_r = ranked[0][1]
isf_best = best_r['isf_fn'](bg_vals, tdd_vals)
v['ratio_best'] = (v['bg_drop_actual'] / v['bg_drop_pred']) * (isf_v1 / isf_best)

print(f"\n  Best formula: {best_r['label']}")
print(f"\n  {'BG Band':>10s}  {'n':>5s}  {'Current ratio':>13s}  {'Best ratio':>10s}  {'Δ':>6s}")
print("  " + "─" * 52)

for band in ['<90', '90-105', '105-120', '120-150', '>150']:
    mask = v['bg_band'] == band
    n = mask.sum()
    if n < 5:
        continue
    r_curr = v.loc[mask, 'ratio'].median()
    r_best = v.loc[mask, 'ratio_best'].median()
    delta = r_best - r_curr
    print(f"  {band:>10s}  {n:5d}  {r_curr:13.3f}  {r_best:10.3f}  {delta:+6.3f}")

# ── 10. Figure ────────────────────────────────────────────────────────────────
print("\nGenerating figure...")

BG_C = '#0f0f0f'; PANEL = '#1a1a2e'; GRID = '#2a2a4a'; TXT = '#e0e0ff'

COLORS = {
    'current': '#4fc3f7',
    'power_log': '#66bb6a',
    'alt_D': '#ffb74d',
    'joint': '#f44336',
    'power_law': '#ce93d8',
    'linear': '#80deea',
    'sqrt_log': '#fff176',
    'flat': '#a5d6a7',
}

def style(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=TXT, labelsize=7)
    ax.set_title(title, color=TXT, fontsize=9, fontweight='bold')
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.grid(True, color=GRID, lw=0.5, ls='--', alpha=0.7)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.xaxis.label.set_fontsize(8); ax.yaxis.label.set_fontsize(8)

fig = plt.figure(figsize=(22, 22))
fig.patch.set_facecolor(BG_C)
gs = gridspec.GridSpec(4, 3, figure=fig, hspace=0.50, wspace=0.35)

bg_range_fine = np.linspace(70, 200, 300)

# P1: ISF curves — all formulas vs empirical (full width)
ax1 = fig.add_subplot(gs[0, :])
style(ax1, f'ISF vs BG: Candidate Scaling Functions (TDD = {tdd_med:.1f} U/day)')

# Empirical points
ax1.errorbar(bg_profile['bg_mid'], bg_profile['isf_eff_med'],
             yerr=[bg_profile['isf_eff_med'] - bg_profile['isf_eff_q25'],
                   bg_profile['isf_eff_q75'] - bg_profile['isf_eff_med']],
             fmt='o', color='white', ms=5, lw=1, capsize=3, capthick=1,
             label='Empirical ISF_eff (median ± IQR)', zorder=10)

# Formula curves
for fname in ['current', 'power_log', 'alt_D', 'joint', 'power_law', 'sqrt_log', 'flat']:
    r = results[fname]
    isf_curve = r['isf_fn'](bg_range_fine, np.full_like(bg_range_fine, tdd_med))
    lw = 2.5 if fname in ['current', 'joint'] else 1.5
    ls = '--' if fname == 'flat' else '-'
    ax1.plot(bg_range_fine, isf_curve, color=COLORS[fname], lw=lw, ls=ls,
             label=r['label'], alpha=0.9)

ax1.axvline(TARGET, color='white', lw=0.8, ls=':', alpha=0.5, label=f'Target {TARGET:.0f}')
ax1.set_xlabel('BG (mg/dL)'); ax1.set_ylabel('ISF (mg/dL per U)')
ax1.set_xlim(70, 200); ax1.set_ylim(0, 400)
ax1.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL, loc='upper right', ncol=2)

# P2: Empirical ISF_eff scatter with density
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'ISF_eff vs BG (individual samples)')
ax2.scatter(v['bg'], v['isf_eff'], alpha=0.08, s=4, color='#4fc3f7', edgecolors='none')
# Overlay empirical median
ax2.plot(bg_profile['bg_mid'], bg_profile['isf_eff_med'], 'o-', color='white', ms=4, lw=2, zorder=10)
# Best fit curve
isf_best_curve = results[best_name]['isf_fn'](bg_range_fine, np.full_like(bg_range_fine, tdd_med))
ax2.plot(bg_range_fine, isf_best_curve, '-', color=COLORS[best_name], lw=2, label=f'Best: {best_name}')
ax2.set_xlabel('BG (mg/dL)'); ax2.set_ylabel('ISF_eff (mg/dL/U)')
ax2.set_xlim(70, 200); ax2.set_ylim(0, 500)
ax2.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P3: MAE bar chart
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'MAE by Scaling Function (+2h)')
names_sorted = [n for n, _ in ranked]
maes_sorted = [r['mae'] for _, r in ranked]
cols_sorted = [COLORS.get(n, '#888888') for n in names_sorted]
labels_sorted = [r['label'].split(':')[0] for _, r in ranked]
bars = ax3.barh(range(len(ranked)), maes_sorted, color=cols_sorted, alpha=0.85)
ax3.set_yticks(range(len(ranked)))
ax3.set_yticklabels(labels_sorted, fontsize=7)
for bar, val in zip(bars, maes_sorted):
    ax3.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
             f'{val:.2f}', va='center', fontsize=7, color=TXT)
ax3.set_xlabel('MAE (mg/dL)')
ax3.invert_yaxis()

# P4: Bias bar chart
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, 'Prediction Bias by Scaling Function (+2h)')
biases_sorted = [r['bias'] for _, r in ranked]
bars = ax4.barh(range(len(ranked)), biases_sorted, color=cols_sorted, alpha=0.85)
ax4.set_yticks(range(len(ranked)))
ax4.set_yticklabels(labels_sorted, fontsize=7)
for bar, val in zip(bars, biases_sorted):
    x = bar.get_width() + 0.2 if val >= 0 else bar.get_width() - 0.8
    ax4.text(x, bar.get_y() + bar.get_height()/2,
             f'{val:+.2f}', va='center', fontsize=7, color=TXT)
ax4.axvline(0, color='white', lw=0.8, ls='--')
ax4.set_xlabel('Mean Error (mg/dL)')
ax4.invert_yaxis()

# P5: BG-band MAE heatmap for top formulas
ax5 = fig.add_subplot(gs[2, 0:2])
style(ax5, 'MAE by BG Band: Current vs Best Alternatives')
bands = ['<90', '90-105', '105-120', '120-150']
top_show = ['current', 'power_log', 'alt_D', 'joint', 'power_law']
x = np.arange(len(bands))
w = 0.15
for i, fname in enumerate(top_show):
    band_maes = []
    for band in bands:
        mask = v['bg_band'] == band
        if mask.sum() >= 5:
            band_maes.append(v.loc[mask, f'err_{fname}'].abs().mean())
        else:
            band_maes.append(0)
    ax5.bar(x + (i - 2) * w, band_maes, w, color=COLORS[fname], alpha=0.85,
            label=results[fname]['label'].split(':')[0])
ax5.set_xticks(x); ax5.set_xticklabels(bands, fontsize=8)
ax5.set_xlabel('BG Band (mg/dL)'); ax5.set_ylabel('MAE (mg/dL)')
ax5.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL, ncol=2)

# P6: Ratio by BG — current vs best
ax6 = fig.add_subplot(gs[2, 2])
style(ax6, 'Sensitivity Ratio by BG: Current vs Best')
bg_bins_r = np.arange(72, 165, 8)
v['bg_bin_r'] = pd.cut(v['bg'], bins=bg_bins_r)
ratio_by_bg = v.groupby('bg_bin_r', observed=True).agg(
    bg_mid=('bg', 'median'),
    ratio_curr=('ratio', 'median'),
    ratio_best=('ratio_best', 'median'),
)
ax6.plot(ratio_by_bg['bg_mid'], ratio_by_bg['ratio_curr'], 'o-', color=COLORS['current'],
         lw=2, ms=5, label='Current')
ax6.plot(ratio_by_bg['bg_mid'], ratio_by_bg['ratio_best'], 's-', color=COLORS[best_name],
         lw=2, ms=5, label=f'Best ({best_name})')
ax6.axhline(1.0, color='white', lw=0.8, ls='--', alpha=0.5)
ax6.set_xlabel('BG (mg/dL)'); ax6.set_ylabel('Ratio (actual/predicted)')
ax6.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P7: Error distribution — current vs best
ax7 = fig.add_subplot(gs[3, 0])
style(ax7, f'Error Distribution: Current vs {best_name}')
ax7.hist(v[f'err_current'].clip(-80, 80), bins=40, alpha=0.5, color=COLORS['current'],
         density=True, label='Current')
ax7.hist(v[f'err_{best_name}'].clip(-80, 80), bins=40, alpha=0.5, color=COLORS[best_name],
         density=True, label=best_name)
ax7.axvline(0, color='white', lw=0.8, ls='--')
ax7.set_xlabel('Pred Error (mg/dL)'); ax7.set_ylabel('Density')
ax7.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P8: ISF ratio (formula / empirical) by BG — shows where each formula diverges
ax8 = fig.add_subplot(gs[3, 1])
style(ax8, 'ISF Formula / ISF_eff by BG')
for fname in ['current', 'power_log', 'joint', 'power_law']:
    isf_f = results[fname]['isf_fn'](bg_profile['bg_mid'].values,
                                      np.full(len(bg_profile), tdd_med))
    ratio_f = isf_f / bg_profile['isf_eff_med'].values
    ax8.plot(bg_profile['bg_mid'], ratio_f, 'o-', color=COLORS[fname], lw=1.5, ms=4,
             label=results[fname]['label'].split(':')[0])
ax8.axhline(1.0, color='white', lw=1, ls='--', alpha=0.5, label='Perfect calibration')
ax8.set_xlabel('BG (mg/dL)'); ax8.set_ylabel('ISF_formula / ISF_empirical')
ax8.set_ylim(0.3, 2.0)
ax8.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P9: Residual BG-dependence check for best formula
ax9 = fig.add_subplot(gs[3, 2])
style(ax9, f'Residual BG-Dependence: {best_name}')
bg_bins_err = np.arange(72, 165, 8)
v['bg_bin_err'] = pd.cut(v['bg'], bins=bg_bins_err)
err_by_bg = v.groupby('bg_bin_err', observed=True).agg(
    bg_mid=('bg', 'median'),
    err_curr=(f'err_current', 'mean'),
    err_best=(f'err_{best_name}', 'mean'),
)
ax9.plot(err_by_bg['bg_mid'], err_by_bg['err_curr'], 'o-', color=COLORS['current'],
         lw=2, ms=5, label='Current')
ax9.plot(err_by_bg['bg_mid'], err_by_bg['err_best'], 's-', color=COLORS[best_name],
         lw=2, ms=5, label=best_name)
ax9.axhline(0, color='white', lw=0.8, ls='--')
ax9.set_xlabel('BG (mg/dL)'); ax9.set_ylabel('Mean Error (mg/dL)')
ax9.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

fig.suptitle('BG Scaling Analysis: Finding the Optimal ISF vs BG Relationship\n'
             f'Overnight 00:00–07:00  |  n={len(v):,} samples  |  +2h horizon  |  '
             f'Median TDD={tdd_med:.1f} U/day',
             color=TXT, fontsize=12, fontweight='bold', y=0.995)

plt.savefig('ns_bg_scaling_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_bg_scaling_results.png")

# ── 11. Summary text ──────────────────────────────────────────────────────────
lines = [
    "BG SCALING ANALYSIS — OVERNIGHT 00:00–07:00",
    "=" * 62,
    f"Samples:  {len(v):,}",
    f"BG range: {BG_LO}–{BG_HI} mg/dL",
    f"Horizon:  +2h",
    f"Anchor:   ISF = 1700/TDD at target {TARGET:.0f} mg/dL",
    "",
    "Question: The current ln(BG/D+1) scaling drops ISF too steeply",
    "as BG rises. What functional form best matches observed overnight",
    "insulin sensitivity?",
    "",
    "Ranking by MAE:",
]
for i, (name, r) in enumerate(ranked):
    lines.append(f"  {i+1}. {r['label']:45s}  MAE={r['mae']:.2f}  bias={r['bias']:+.2f}  "
                 f"±1mmol={r['w18']:.1f}%  [{r['params']}]")

lines += [
    "",
    f"Best: {ranked[0][1]['label']}",
    f"Improvement over current: MAE {results['current']['mae']:.2f} → {ranked[0][1]['mae']:.2f} "
    f"({(results['current']['mae'] - ranked[0][1]['mae']):.2f} mg/dL)",
    "",
    "BG-band MAE (current → best):",
]
for band in ['<90', '90-105', '105-120', '120-150']:
    mask = v['bg_band'] == band
    if mask.sum() < 5:
        continue
    mae_curr = v.loc[mask, 'err_current'].abs().mean()
    mae_best = v.loc[mask, f'err_{best_name}'].abs().mean()
    lines.append(f"  {band:>10s}: {mae_curr:.1f} → {mae_best:.1f}  ({mae_best-mae_curr:+.1f})")

lines += [
    "",
    f"Key finding: Optimal power exponent p = {p_best:.3f} (current uses p = 1.0).",
    f"The current formula's ln scaling is too steep — ISF drops too quickly",
    f"as BG rises. A flatter curve (p < 1) reduces over-aggression at high BG",
    f"without sacrificing accuracy at low BG.",
]

with open('ns_bg_scaling_summary.txt', 'w') as f:
    f.write('\n'.join(lines))

print('\n' + '\n'.join(lines))
print("\nSaved: ns_bg_scaling_summary.txt")
