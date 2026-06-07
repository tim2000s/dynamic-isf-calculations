"""
TDD_effective v2 — Refined Overnight Sensitivity Regression
============================================================
Refinements over v1:
  1. Confidence-weighted nightly aggregation (weight by predicted BG drop)
  2. Dawn Phenomenon stratification (00:00–03:30 vs 03:30–07:00)
  3. Adaptive learning rate (scales α by nightly consistency)
  4. BG-band ratio analysis (tests whether correction is BG-dependent)
  5. Direct ratio approach (TDD_eff = TDD_7day / ratio, no v1 routing)

Outputs
-------
  ns_tdd_effective_v2_results.png
  ns_tdd_effective_v2_summary.txt
  ns_tdd_effective_v2_nightly.csv
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats as sp_stats

warnings.filterwarnings('ignore')

D      = 82.0
TARGET = 99.0
LN_TARGET = np.log(TARGET / D + 1)

# Filters
MIN_BG_DROP    = 3.0
BG_LO, BG_HI  = 72, 200
MIN_SAMPLES    = 5

# Learning rate
ALPHA_BASE     = 0.15

# Dawn split
DAWN_HOUR      = 3.5   # 03:30

# ── 1. Load ───────────────────────────────────────────────────────────────────
print("Loading ns_backtest_overnight.csv ...")
df = pd.read_csv('ns_backtest_overnight.csv')
df['ts'] = pd.to_datetime(df['ts'], format='ISO8601')
df['date'] = pd.to_datetime(df['date']).dt.date
df['hour'] = df['ts'].dt.hour
df['hour_frac'] = df['ts'].dt.hour + df['ts'].dt.minute / 60.0
print(f"  Total overnight rows: {len(df):,}")

# ── 2. Filter ─────────────────────────────────────────────────────────────────
v = df.dropna(subset=['pred_iob_24', 'actual_bg_2h', 'isf_v1', 'tdd_7day']).copy()
v = v[(v['bg'] >= BG_LO) & (v['bg'] <= BG_HI)]

v['bg_drop_pred']   = v['bg'] - v['pred_iob_24']
v['bg_drop_actual']  = v['bg'] - v['actual_bg_2h']

v = v[v['bg_drop_pred'].abs() >= MIN_BG_DROP]
v['ratio'] = v['bg_drop_actual'] / v['bg_drop_pred']
v = v[(v['ratio'] > 0) & (v['ratio'] < 5)]
v = v[~((v['bg_drop_pred'] > 0) & (v['bg_drop_actual'] < -9))]

# Back-calculate TDD_eff per sample (via ISF_v1 route)
v['isf_eff'] = v['isf_v1'] * v['ratio']
v['tdd_eff_sample'] = (1700.0 * LN_TARGET) / (v['isf_eff'] * v['ln_bg'])
v = v[(v['tdd_eff_sample'] > 3) & (v['tdd_eff_sample'] < 120)]

# Confidence weight: magnitude of predicted BG drop
v['weight'] = v['bg_drop_pred'].abs()

# Dawn stratification
v['period'] = np.where(v['hour_frac'] < DAWN_HOUR, 'deep_night', 'pre_dawn')

# BG bands
v['bg_band'] = pd.cut(v['bg'], bins=[0, 90, 120, 150, 300],
                       labels=['<90', '90-120', '120-150', '>150'])

print(f"  Valid 2h samples: {len(v):,}")
print(f"    Deep night (00:00–03:30): {(v['period']=='deep_night').sum():,}")
print(f"    Pre-dawn   (03:30–07:00): {(v['period']=='pre_dawn').sum():,}")

# ── 3. BG-Band Ratio Analysis ────────────────────────────────────────────────
print("\n" + "═" * 72)
print("BG-BAND SENSITIVITY RATIO ANALYSIS")
print("═" * 72)
print(f"\n  {'BG Band':>10s}  {'n':>5s}  {'ratio_med':>9s}  {'ratio_IQR':>16s}  {'TDD_eff_med':>11s}")
print("  " + "─" * 60)

bg_band_ratios = {}
for band, grp in v.groupby('bg_band', observed=True):
    r = grp['ratio']
    t = grp['tdd_eff_sample']
    print(f"  {band:>10s}  {len(grp):5d}  {r.median():9.3f}  "
          f"[{r.quantile(.25):.3f}, {r.quantile(.75):.3f}]  {t.median():11.1f}")
    bg_band_ratios[band] = r.median()

# Kruskal-Wallis test: is ratio significantly different across BG bands?
groups = [grp['ratio'].values for _, grp in v.groupby('bg_band', observed=True)]
if len(groups) >= 2 and all(len(g) >= 5 for g in groups):
    stat, pval = sp_stats.kruskal(*groups)
    print(f"\n  Kruskal-Wallis test: H={stat:.2f}, p={pval:.4f}")
    if pval < 0.05:
        print("  → Ratio varies significantly by BG band. The correction IS BG-dependent.")
        print("    A single TDD_effective may not fully capture this variation.")
    else:
        print("  → No significant BG-band dependence. A scalar TDD correction is appropriate.")

# ── 4. Dawn Phenomenon Analysis ───────────────────────────────────────────────
print("\n" + "═" * 72)
print("DAWN PHENOMENON STRATIFICATION")
print("═" * 72)

for period in ['deep_night', 'pre_dawn']:
    grp = v[v['period'] == period]
    r = grp['ratio']
    t = grp['tdd_eff_sample']
    label = "Deep Night (00:00–03:30)" if period == 'deep_night' else "Pre-Dawn   (03:30–07:00)"
    print(f"\n  {label}  n={len(grp):,}")
    print(f"    ratio:      median={r.median():.3f}  IQR=[{r.quantile(.25):.3f}, {r.quantile(.75):.3f}]")
    print(f"    TDD_eff:    median={t.median():.1f}  IQR=[{t.quantile(.25):.1f}, {t.quantile(.75):.1f}]")
    print(f"    ISF_eff:    median={grp['isf_eff'].median():.1f}")

# Test dawn vs deep night
dn = v[v['period'] == 'deep_night']['ratio']
pd_r = v[v['period'] == 'pre_dawn']['ratio']
if len(dn) >= 10 and len(pd_r) >= 10:
    stat, pval = sp_stats.mannwhitneyu(dn, pd_r, alternative='two-sided')
    print(f"\n  Mann-Whitney U test (deep night vs pre-dawn): U={stat:.0f}, p={pval:.4f}")
    dawn_effect = pd_r.median() - dn.median()
    print(f"  Dawn effect on ratio: {dawn_effect:+.3f} "
          f"({'less sensitive pre-dawn' if dawn_effect < 0 else 'more sensitive pre-dawn'})")

# ── 5. Nightly Aggregation — Three Methods ────────────────────────────────────

def weighted_median(values, weights):
    """Weighted median: sort by value, find the weight that crosses 50%."""
    s = np.argsort(values)
    v_s, w_s = values[s], weights[s]
    cumw = np.cumsum(w_s)
    cutoff = cumw[-1] * 0.5
    idx = np.searchsorted(cumw, cutoff)
    return v_s[min(idx, len(v_s)-1)]

nightly_rows = []
for date, grp in v.groupby('date'):
    if len(grp) < MIN_SAMPLES:
        continue

    row = {'date': date, 'n_samples': len(grp)}

    # Method A: Simple median (v1 approach)
    row['tdd_impl_simple'] = grp['tdd_eff_sample'].median()
    row['ratio_simple']    = grp['ratio'].median()

    # Method B: Confidence-weighted median
    vals = grp['tdd_eff_sample'].values
    wts  = grp['weight'].values
    row['tdd_impl_weighted'] = weighted_median(vals, wts)
    row['ratio_weighted']    = weighted_median(grp['ratio'].values, wts)

    # Method C: Direct ratio approach (TDD_eff = TDD_7day / ratio)
    row['ratio_direct']    = grp['ratio'].median()
    row['tdd_7day_med']    = grp['tdd_7day'].median()
    row['tdd_impl_direct'] = row['tdd_7day_med'] / row['ratio_direct'] if row['ratio_direct'] > 0 else np.nan

    # Deep-night only
    dn_grp = grp[grp['period'] == 'deep_night']
    if len(dn_grp) >= 3:
        row['tdd_impl_deepnight'] = dn_grp['tdd_eff_sample'].median()
        row['ratio_deepnight']    = dn_grp['ratio'].median()
        row['n_deepnight']        = len(dn_grp)
    else:
        row['tdd_impl_deepnight'] = row['tdd_impl_simple']  # fallback
        row['ratio_deepnight']    = row['ratio_simple']
        row['n_deepnight']        = len(dn_grp)

    # IQR for adaptive α
    iqr = grp['tdd_eff_sample'].quantile(0.75) - grp['tdd_eff_sample'].quantile(0.25)
    med = grp['tdd_eff_sample'].median()
    row['cv_night'] = (iqr / med) if med > 0 else 99.0

    row['bg_med']        = grp['bg'].median()
    row['isf_eff_med']   = grp['isf_eff'].median()
    row['isf_7dtdd_med'] = grp['isf_7dtdd'].median()

    nightly_rows.append(row)

nightly = pd.DataFrame(nightly_rows).sort_values('date').reset_index(drop=True)
print(f"\n  Nights with ≥{MIN_SAMPLES} valid samples: {len(nightly)}")

# ── 6. Rolling TDD_effective — All Methods ────────────────────────────────────

methods = [
    ('simple',     'tdd_impl_simple',    'Fixed α=0.15'),
    ('weighted',   'tdd_impl_weighted',  'Weighted median, fixed α=0.15'),
    ('deepnight',  'tdd_impl_deepnight', 'Deep-night only, fixed α=0.15'),
    ('adaptive',   'tdd_impl_weighted',  'Weighted median, adaptive α'),
    ('direct',     'tdd_impl_direct',    'Direct ratio, fixed α=0.15'),
]

for mname, impl_col, mlabel in methods:
    col_eff  = f'tdd_eff_{mname}'
    col_next = f'tdd_eff_next_{mname}'
    col_alpha = f'alpha_{mname}'

    vals = []
    alphas_used = []
    tdd_eff_prev = None

    for i, row in nightly.iterrows():
        if tdd_eff_prev is None:
            tdd_eff_prev = row['tdd_7day_med']

        # Adaptive α: scale by nightly consistency
        if mname == 'adaptive':
            cv = row['cv_night']
            alpha = ALPHA_BASE * min(1.0, 1.0 / (1.0 + cv))
            alpha = max(0.05, min(0.30, alpha))  # floor/ceiling
        else:
            alpha = ALPHA_BASE

        impl = row[impl_col]
        if pd.isna(impl):
            vals.append(tdd_eff_prev)
            alphas_used.append(0.0)
            continue

        tdd_eff_new = alpha * impl + (1 - alpha) * tdd_eff_prev
        vals.append(tdd_eff_new)
        alphas_used.append(alpha)
        tdd_eff_prev = tdd_eff_new

    nightly[col_eff]  = vals
    nightly[col_alpha] = alphas_used
    nightly[col_next] = nightly[col_eff].shift(1)  # next-day (no look-ahead)

# ── 7. Backtest All Methods ───────────────────────────────────────────────────
print("\n" + "═" * 72)
print("BACKTEST COMPARISON (+2h)")
print("═" * 72)

# Merge all next-day TDD_effective values to samples
merge_cols = ['date'] + [f'tdd_eff_next_{m}' for m, _, _ in methods]
v_bt = v.merge(nightly[merge_cols], on='date', how='left')

# Need at least one method's next-day value
v_bt = v_bt.dropna(subset=['tdd_eff_next_simple'])
print(f"  Samples with prior-night TDD_effective: {len(v_bt):,}")

# Compute predictions for each method
for mname, _, _ in methods:
    col_next = f'tdd_eff_next_{mname}'
    isf_col  = f'isf_{mname}'
    pred_col = f'pred_{mname}'
    err_col  = f'err_{mname}'

    v_bt[isf_col]  = (1700.0 / v_bt[col_next]) * LN_TARGET / v_bt['ln_bg']
    v_bt[pred_col] = v_bt['bg'] - v_bt['bg_drop_pred'] * (v_bt[isf_col] / v_bt['isf_v1'])
    v_bt[err_col]  = v_bt['actual_bg_2h'] - v_bt[pred_col]

# Also compute static 7D-TDD and v1 errors
v_bt['pred_7dtdd_s'] = v_bt['bg'] - v_bt['bg_drop_pred'] * (v_bt['isf_7dtdd'] / v_bt['isf_v1'])
v_bt['err_7dtdd_s']  = v_bt['actual_bg_2h'] - v_bt['pred_7dtdd_s']
v_bt['err_v1']       = v_bt['actual_bg_2h'] - v_bt['pred_iob_24']

# Summary table
print(f"\n  {'Method':<38s}  {'Bias':>6s}  {'MAE':>6s}  {'RMSE':>6s}  {'±1mmol':>7s}  {'±2mmol':>7s}")
print("  " + "─" * 78)

all_results = []

for label, err_col in [('v1 (Actual loop)', 'err_v1'),
                        ('7D-TDD (static)', 'err_7dtdd_s')]:
    e = v_bt[err_col].dropna()
    row = dict(method=label, bias=e.mean(), mae=e.abs().mean(),
               rmse=np.sqrt((e**2).mean()),
               w18=(e.abs()<=18).mean()*100, w36=(e.abs()<=36).mean()*100)
    all_results.append(row)
    print(f"  {label:<38s}  {row['bias']:+6.1f}  {row['mae']:6.1f}  {row['rmse']:6.1f}  "
          f"{row['w18']:6.1f}%  {row['w36']:6.1f}%")

for mname, _, mlabel in methods:
    err_col = f'err_{mname}'
    e = v_bt[err_col].dropna()
    row = dict(method=mlabel, bias=e.mean(), mae=e.abs().mean(),
               rmse=np.sqrt((e**2).mean()),
               w18=(e.abs()<=18).mean()*100, w36=(e.abs()<=36).mean()*100)
    all_results.append(row)
    print(f"  {mlabel:<38s}  {row['bias']:+6.1f}  {row['mae']:6.1f}  {row['rmse']:6.1f}  "
          f"{row['w18']:6.1f}%  {row['w36']:6.1f}%")

# ── 8. Per-Night Breakdown: Best Method vs Static ────────────────────────────
# Determine best adaptive method by MAE
best_method = min([(m, v_bt[f'err_{m}'].abs().mean()) for m, _, _ in methods], key=lambda x: x[1])
best_name = best_method[0]
print(f"\n  Best adaptive method by MAE: {best_name} ({best_method[1]:.1f} mg/dL)")

print(f"\n── Per-Night: Static 7D-TDD vs Best Adaptive ({best_name}) " + "─" * 20)
print(f"  {'Date':>12s}  {'n':>3s}  │ {'v1 MAE':>7s}  {'Stat MAE':>8s}  {'Adpt MAE':>8s}  │ "
      f"{'v1 bias':>7s}  {'Stat bias':>9s}  {'Adpt bias':>9s}  │ "
      f"{'TDD_7d':>6s}  {'TDD_eff':>7s}")
print("  " + "─" * 110)

nights_improved = 0
nights_total = 0
for date, grp in v_bt.groupby('date'):
    n = len(grp)
    tdd7 = grp['tdd_7day'].median()
    tdd_e = grp[f'tdd_eff_next_{best_name}'].median()
    mae_v1 = grp['err_v1'].abs().mean()
    mae_s  = grp['err_7dtdd_s'].abs().mean()
    mae_a  = grp[f'err_{best_name}'].abs().mean()
    bias_v1 = grp['err_v1'].mean()
    bias_s  = grp['err_7dtdd_s'].mean()
    bias_a  = grp[f'err_{best_name}'].mean()
    marker = " ✓" if mae_a < mae_s else ""
    if mae_a < mae_s:
        nights_improved += 1
    nights_total += 1
    print(f"  {str(date):>12s}  {n:3d}  │ "
          f"{mae_v1:7.1f}  {mae_s:8.1f}  {mae_a:8.1f}{marker:2s} │ "
          f"{bias_v1:+7.1f}  {bias_s:+9.1f}  {bias_a:+9.1f}  │ "
          f"{tdd7:6.1f}  {tdd_e:7.1f}")

print(f"\n  Adaptive wins on {nights_improved}/{nights_total} nights "
      f"({nights_improved/nights_total*100:.0f}%)")

# ── 9. Convergence and Stability ─────────────────────────────────────────────
print(f"\n── TDD_effective Convergence " + "─" * 45)
print(f"  {'Night':>5s}  {'Date':>12s}  {'TDD_7d':>6s}  │ "
      f"{'Simple':>7s}  {'Weight':>7s}  {'DeepN':>7s}  {'Adapt':>7s}  {'Direct':>7s}  │ "
      f"{'α_adapt':>7s}")
print("  " + "─" * 90)

for i, (_, row) in enumerate(nightly.iterrows()):
    if i >= 30:
        print(f"  ... ({len(nightly) - 30} more nights)")
        break
    print(f"  {i+1:5d}  {str(row['date']):>12s}  {row['tdd_7day_med']:6.1f}  │ "
          f"{row['tdd_eff_simple']:7.1f}  {row['tdd_eff_weighted']:7.1f}  "
          f"{row['tdd_eff_deepnight']:7.1f}  {row['tdd_eff_adaptive']:7.1f}  "
          f"{row['tdd_eff_direct']:7.1f}  │ "
          f"{row['alpha_adaptive']:7.3f}")

# ── 10. Sensitivity Offset Summary ───────────────────────────────────────────
print(f"\n── Sensitivity Offset (TDD_effective − TDD_7day) " + "─" * 25)
for mname, _, mlabel in methods:
    col = f'tdd_eff_{mname}'
    deltas = nightly[col] - nightly['tdd_7day_med']
    d = deltas.dropna()
    print(f"  {mlabel:<38s}  median={d.median():+.2f}  mean={d.mean():+.2f}  "
          f"std={d.std():.2f} U/day")

# ── 11. Figure ────────────────────────────────────────────────────────────────
print("\nGenerating figure...")

BG_C = '#0f0f0f'; PANEL = '#1a1a2e'; GRID = '#2a2a4a'; TXT = '#e0e0ff'
C_V1 = '#4fc3f7'; C_7D = '#ce93d8'; C_SIMP = '#66bb6a'; C_WGHT = '#ffb74d'
C_DEEP = '#f48fb1'; C_ADPT = '#80deea'; C_DIRE = '#fff176'

def style(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=TXT, labelsize=7)
    ax.set_title(title, color=TXT, fontsize=9, fontweight='bold')
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.grid(True, color=GRID, lw=0.5, ls='--', alpha=0.7)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.xaxis.label.set_fontsize(7); ax.yaxis.label.set_fontsize(7)

fig = plt.figure(figsize=(22, 24))
fig.patch.set_facecolor(BG_C)
gs = gridspec.GridSpec(5, 3, figure=fig, hspace=0.55, wspace=0.38)

dates_num = matplotlib.dates.date2num([pd.Timestamp(d) for d in nightly['date']])

# P1: TDD timeline — all methods (full width)
ax1 = fig.add_subplot(gs[0, :])
style(ax1, 'TDD_effective: All Methods vs TDD_7day')
ax1.plot(dates_num, nightly['tdd_7day_med'], 'o-', color=C_7D, lw=1.5, ms=3,
         label='TDD 7-day', alpha=0.7)
ax1.plot(dates_num, nightly['tdd_eff_simple'], '-', color=C_SIMP, lw=1.5,
         label='Simple median', alpha=0.8)
ax1.plot(dates_num, nightly['tdd_eff_weighted'], '-', color=C_WGHT, lw=2,
         label='Confidence-weighted', alpha=0.9)
ax1.plot(dates_num, nightly['tdd_eff_deepnight'], '-', color=C_DEEP, lw=1.5,
         label='Deep-night only', alpha=0.8)
ax1.plot(dates_num, nightly['tdd_eff_adaptive'], '-', color=C_ADPT, lw=1.5,
         label='Adaptive α', alpha=0.8)
ax1.plot(dates_num, nightly['tdd_eff_direct'], '--', color=C_DIRE, lw=1.5,
         label='Direct ratio', alpha=0.7)
ax1.xaxis_date(); ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))
ax1.set_ylabel('TDD (U/day)'); ax1.set_xlabel('')
ax1.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL, loc='best', ncol=3)
for tick in ax1.get_xticklabels(): tick.set_rotation(45); tick.set_fontsize(6)

# P2: BG-band ratio boxplot
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'Sensitivity Ratio by BG Band')
band_data = [v[v['bg_band'] == b]['ratio'].values for b in ['<90', '90-120', '120-150', '>150']]
bp = ax2.boxplot(band_data, labels=['<90', '90-120', '120-150', '>150'],
                  patch_artist=True, widths=0.6,
                  medianprops=dict(color='white', lw=2),
                  flierprops=dict(marker='.', ms=2, alpha=0.3))
colors_box = ['#66bb6a', '#4fc3f7', '#ffb74d', '#ef5350']
for patch, c in zip(bp['boxes'], colors_box):
    patch.set_facecolor(c); patch.set_alpha(0.6)
ax2.axhline(1.0, color='white', lw=0.8, ls='--', alpha=0.5)
ax2.set_ylabel('Ratio (actual/predicted drop)')
ax2.set_xlabel('BG Band (mg/dL)')

# P3: Dawn vs Deep Night ratio distributions
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'Sensitivity Ratio: Deep Night vs Pre-Dawn')
dn_vals = v[v['period'] == 'deep_night']['ratio'].clip(0, 3)
pd_vals = v[v['period'] == 'pre_dawn']['ratio'].clip(0, 3)
ax3.hist(dn_vals, bins=30, alpha=0.6, color=C_DEEP, density=True, label='Deep Night (00–03:30)')
ax3.hist(pd_vals, bins=30, alpha=0.6, color=C_ADPT, density=True, label='Pre-Dawn (03:30–07)')
ax3.axvline(dn_vals.median(), color=C_DEEP, lw=1.5, ls='-')
ax3.axvline(pd_vals.median(), color=C_ADPT, lw=1.5, ls='-')
ax3.axvline(1.0, color='white', lw=0.8, ls='--', alpha=0.5)
ax3.set_xlabel('Ratio'); ax3.set_ylabel('Density')
ax3.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P4: Adaptive α over time
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, 'Adaptive Learning Rate (α) Over Time')
ax4.plot(range(len(nightly)), nightly['alpha_adaptive'], 'o-', color=C_ADPT, ms=3, lw=1)
ax4.axhline(ALPHA_BASE, color='white', lw=0.8, ls='--', alpha=0.5,
            label=f'Base α={ALPHA_BASE}')
ax4.set_xlabel('Night #'); ax4.set_ylabel('α used')
ax4.set_ylim(0, 0.35)
ax4.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P5: MAE comparison bar chart
ax5 = fig.add_subplot(gs[2, 0])
style(ax5, 'MAE Comparison (+2h)')
labels_bar = ['v1\nActual', '7D-TDD\nStatic', 'Simple', 'Weighted', 'Deep\nNight', 'Adapt α', 'Direct\nRatio']
err_cols = ['err_v1', 'err_7dtdd_s', 'err_simple', 'err_weighted', 'err_deepnight', 'err_adaptive', 'err_direct']
maes = [v_bt[c].abs().mean() for c in err_cols]
cols_bar = [C_V1, C_7D, C_SIMP, C_WGHT, C_DEEP, C_ADPT, C_DIRE]
bars = ax5.bar(range(len(labels_bar)), maes, color=cols_bar, alpha=0.85)
for bar, val in zip(bars, maes):
    ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
             f'{val:.1f}', ha='center', fontsize=7, color=TXT)
ax5.set_xticks(range(len(labels_bar)))
ax5.set_xticklabels(labels_bar, fontsize=6)
ax5.set_ylabel('MAE (mg/dL)')

# P6: Bias comparison
ax6 = fig.add_subplot(gs[2, 1])
style(ax6, 'Prediction Bias (+2h)')
biases = [v_bt[c].mean() for c in err_cols]
bars = ax6.bar(range(len(labels_bar)), biases, color=cols_bar, alpha=0.85)
for bar, val in zip(bars, biases):
    ax6.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + (0.15 if val >= 0 else -1.0),
             f'{val:+.1f}', ha='center', fontsize=7, color=TXT)
ax6.axhline(0, color='white', lw=0.8, ls='--')
ax6.set_xticks(range(len(labels_bar)))
ax6.set_xticklabels(labels_bar, fontsize=6)
ax6.set_ylabel('Mean Error (mg/dL)')

# P7: ±1mmol comparison
ax7 = fig.add_subplot(gs[2, 2])
style(ax7, '% Within ±1 mmol/L (+2h)')
w18s = [(v_bt[c].abs() <= 18).mean() * 100 for c in err_cols]
bars = ax7.bar(range(len(labels_bar)), w18s, color=cols_bar, alpha=0.85)
for bar, val in zip(bars, w18s):
    ax7.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
             f'{val:.0f}%', ha='center', fontsize=7, color=TXT)
ax7.set_xticks(range(len(labels_bar)))
ax7.set_xticklabels(labels_bar, fontsize=6)
ax7.set_ylabel('%')

# P8: Error distributions — static vs best adaptive
ax8 = fig.add_subplot(gs[3, 0])
style(ax8, f'Error Distribution: Static vs {best_name}')
ax8.hist(v_bt['err_7dtdd_s'].clip(-80, 80), bins=35, alpha=0.5, color=C_7D,
         density=True, label='7D-TDD static')
ax8.hist(v_bt[f'err_{best_name}'].clip(-80, 80), bins=35, alpha=0.5, color=C_SIMP,
         density=True, label=f'7D-TDD {best_name}')
ax8.axvline(0, color='white', lw=0.8, ls='--')
ax8.set_xlabel('Pred Error (mg/dL)'); ax8.set_ylabel('Density')
ax8.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P9: Per-night MAE timeline — static vs best adaptive
ax9 = fig.add_subplot(gs[3, 1:3])
style(ax9, f'Per-Night MAE: Static vs {best_name} Adaptive')
night_dates_bt = []
night_mae_stat = []
night_mae_best = []
night_mae_v1 = []
for date, grp in v_bt.groupby('date'):
    night_dates_bt.append(date)
    night_mae_v1.append(grp['err_v1'].abs().mean())
    night_mae_stat.append(grp['err_7dtdd_s'].abs().mean())
    night_mae_best.append(grp[f'err_{best_name}'].abs().mean())
nd_num = matplotlib.dates.date2num([pd.Timestamp(d) for d in night_dates_bt])
ax9.plot(nd_num, night_mae_v1, 'o-', color=C_V1, lw=1, ms=3, label='v1 Actual', alpha=0.5)
ax9.plot(nd_num, night_mae_stat, 's-', color=C_7D, lw=1.5, ms=3, label='7D-TDD static')
ax9.plot(nd_num, night_mae_best, 'D-', color=C_SIMP, lw=2, ms=4, label=f'7D-TDD {best_name}')
ax9.xaxis_date(); ax9.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))
ax9.set_ylabel('MAE (mg/dL)'); ax9.set_xlabel('')
ax9.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)
for tick in ax9.get_xticklabels(): tick.set_rotation(45); tick.set_fontsize(6)

# P10: Sensitivity offset by method
ax10 = fig.add_subplot(gs[4, 0])
style(ax10, 'Sensitivity Offset by Method (TDD_eff − TDD_7day)')
method_labels = ['Simple', 'Weighted', 'Deep\nNight', 'Adapt α', 'Direct']
method_names = [m for m, _, _ in methods]
offsets_med = []
offsets_iqr = []
for mname in method_names:
    d = (nightly[f'tdd_eff_{mname}'] - nightly['tdd_7day_med']).dropna()
    offsets_med.append(d.median())
    offsets_iqr.append((d.quantile(0.25), d.quantile(0.75)))
cols_m = [C_SIMP, C_WGHT, C_DEEP, C_ADPT, C_DIRE]
bars = ax10.bar(range(5), offsets_med, color=cols_m, alpha=0.85)
for i, (bar, val) in enumerate(zip(bars, offsets_med)):
    ax10.errorbar(i, val, yerr=[[val - offsets_iqr[i][0]], [offsets_iqr[i][1] - val]],
                  color='white', lw=1, capsize=4)
    ax10.text(bar.get_x() + bar.get_width()/2, val - 0.3,
              f'{val:+.1f}', ha='center', fontsize=7, color=TXT)
ax10.axhline(0, color='white', lw=0.8, ls='--')
ax10.set_xticks(range(5)); ax10.set_xticklabels(method_labels, fontsize=7)
ax10.set_ylabel('Δ TDD (U/day)')

# P11: Ratio by hour of night
ax11 = fig.add_subplot(gs[4, 1])
style(ax11, 'Sensitivity Ratio by Hour of Night')
hourly = v.groupby('hour')['ratio'].agg(['median', lambda x: x.quantile(.25), lambda x: x.quantile(.75)])
hourly.columns = ['median', 'q25', 'q75']
ax11.fill_between(hourly.index, hourly['q25'], hourly['q75'], alpha=0.2, color=C_WGHT)
ax11.plot(hourly.index, hourly['median'], 'o-', color=C_WGHT, lw=2, ms=5)
ax11.axhline(1.0, color='white', lw=0.8, ls='--', alpha=0.5)
ax11.axvspan(DAWN_HOUR, 7, alpha=0.1, color=C_DEEP, label='Pre-dawn zone')
ax11.set_xlabel('Hour'); ax11.set_ylabel('Ratio')
ax11.set_xticks(range(8))
ax11.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P12: Nightly CV (confidence measure)
ax12 = fig.add_subplot(gs[4, 2])
style(ax12, 'Nightly Consistency (IQR/Median of TDD_implied)')
ax12.bar(range(len(nightly)), nightly['cv_night'].clip(0, 3), color=C_ADPT, alpha=0.7)
ax12.axhline(nightly['cv_night'].median(), color='white', lw=1, ls='--',
             label=f"median CV={nightly['cv_night'].median():.2f}")
ax12.set_xlabel('Night #'); ax12.set_ylabel('CV (IQR/median)')
ax12.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

fig.suptitle('TDD_effective v2: Refined Overnight Sensitivity Regression\n'
             'Confidence weighting, dawn stratification, adaptive α, BG-band analysis',
             color=TXT, fontsize=13, fontweight='bold', y=0.995)

plt.savefig('ns_tdd_effective_v2_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_tdd_effective_v2_results.png")

# ── 12. Save ──────────────────────────────────────────────────────────────────
nightly.to_csv('ns_tdd_effective_v2_nightly.csv', index=False)

# Summary text
lines = [
    "TDD_effective v2 — REFINED OVERNIGHT SENSITIVITY REGRESSION",
    "=" * 62,
    f"Input:           ns_backtest_overnight.csv",
    f"Valid samples:   {len(v):,}",
    f"Nights used:     {len(nightly)}",
    f"Base α:          {ALPHA_BASE}",
    f"Dawn split:      {DAWN_HOUR:.1f}h",
    "",
    "Methods compared:",
    "  1. Simple median (v1 approach, fixed α=0.15)",
    "  2. Confidence-weighted median (weight by |predicted BG drop|, fixed α=0.15)",
    "  3. Deep-night only (00:00–03:30, avoids Dawn Phenomenon, fixed α=0.15)",
    "  4. Adaptive α (weighted median, α scales by nightly consistency)",
    "  5. Direct ratio (TDD_eff = TDD_7day / ratio, no v1 ISF routing)",
    "",
]

lines.append("BG-Band ratio analysis:")
for band, grp in v.groupby('bg_band', observed=True):
    r = grp['ratio']
    lines.append(f"  {band:>10s}  n={len(grp):4d}  ratio_med={r.median():.3f}  "
                 f"IQR=[{r.quantile(.25):.3f}, {r.quantile(.75):.3f}]")

lines.append("")
lines.append("Dawn Phenomenon:")
for period in ['deep_night', 'pre_dawn']:
    grp = v[v['period'] == period]
    r = grp['ratio']
    label = "Deep Night (00:00–03:30)" if period == 'deep_night' else "Pre-Dawn   (03:30–07:00)"
    lines.append(f"  {label}  n={len(grp):,}  ratio_med={r.median():.3f}")

lines.append("")
lines.append("Backtest results (+2h):")
lines.append(f"  {'Method':<38s}  {'Bias':>6s}  {'MAE':>6s}  {'RMSE':>6s}  {'±1mmol':>7s}  {'±2mmol':>7s}")
lines.append("  " + "─" * 72)
for r in all_results:
    lines.append(f"  {r['method']:<38s}  {r['bias']:+6.1f}  {r['mae']:6.1f}  {r['rmse']:6.1f}  "
                 f"{r['w18']:6.1f}%  {r['w36']:6.1f}%")

lines.append("")
lines.append("Sensitivity offset (TDD_eff − TDD_7day):")
for mname, _, mlabel in methods:
    col = f'tdd_eff_{mname}'
    d = (nightly[col] - nightly['tdd_7day_med']).dropna()
    lines.append(f"  {mlabel:<38s}  median={d.median():+.2f}  mean={d.mean():+.2f}")

with open('ns_tdd_effective_v2_summary.txt', 'w') as f:
    f.write('\n'.join(lines))

print('\n' + '\n'.join(lines))
print("\nSaved: ns_tdd_effective_v2_nightly.csv, ns_tdd_effective_v2_summary.txt")
