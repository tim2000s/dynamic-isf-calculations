"""
TDD_effective — Overnight Sensitivity Regression
=================================================
Reads ns_backtest_overnight.csv (output of ns_overnight_backtest.py)
and implements a nightly rolling TDD correction:

  At 07:00 each morning, compare what the 7D-TDD formula predicted
  overnight (00:00–07:00) against what actually happened at +2h.
  Back-calculate the TDD that WOULD have made the prediction correct.
  Exponentially smooth that into TDD_effective for the next day.

Concept
-------
  ISF_7dtdd   = (1700 / TDD_7day) × ln(target/D+1) / ln(BG/D+1)

  The prediction residual at +2h tells us whether TDD_7day was
  too high (formula too aggressive) or too low (too conservative).

  ISF_eff(t)  = ISF_v1 × (BG_t − BG_actual_2h) / (BG_t − BG_pred_2h)
              = the ISF that would have given zero prediction error

  TDD_eff(t)  = 1700 × ln(target/D+1) / (ISF_eff(t) × ln(BG/D+1))
              = the TDD input that produces ISF_eff via the 7D-TDD formula

  Nightly:    TDD_implied_night = median( TDD_eff(t) ) over valid samples
  Update:     TDD_effective_new = α × TDD_implied_night + (1−α) × TDD_effective_prev

  Next day:   ISF_adaptive = (1700 / TDD_effective) × ln(target/D+1) / ln(BG/D+1)

  The formula structure is unchanged — only the TDD input is personalised.

Outputs
-------
  ns_tdd_effective_results.png
  ns_tdd_effective_summary.txt
  ns_tdd_effective_nightly.csv
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

D      = 82.0
TARGET = 99.0
LN_TARGET = np.log(TARGET / D + 1)

# Minimum IOB delta consumed over 2h window to trust the sample
MIN_IOB_DELTA  = 0.3   # units
# Minimum BG drop predicted — avoid near-zero denominators
MIN_BG_DROP    = 3.0    # mg/dL
# BG bounds — exclude extreme readings
BG_LO, BG_HI  = 72, 200

# Learning rates to compare
ALPHAS = [0.10, 0.15, 0.20, 0.30]
ALPHA_DEFAULT = 0.15

# ── 1. Load ───────────────────────────────────────────────────────────────────
print("Loading ns_backtest_overnight.csv ...")
df = pd.read_csv('ns_backtest_overnight.csv')
df['ts'] = pd.to_datetime(df['ts'], format='ISO8601')
df['date'] = pd.to_datetime(df['date']).dt.date
df['hour'] = df['ts'].dt.hour
print(f"  Total overnight rows: {len(df):,}")

# ── 2. Filter to valid 2h prediction samples ─────────────────────────────────
v = df.dropna(subset=['pred_iob_24', 'actual_bg_2h', 'isf_v1', 'tdd_7day']).copy()
v = v[(v['bg'] >= BG_LO) & (v['bg'] <= BG_HI)]

# Predicted and actual BG drop
v['bg_drop_pred'] = v['bg'] - v['pred_iob_24']       # loop's predicted drop
v['bg_drop_actual'] = v['bg'] - v['actual_bg_2h']     # what actually happened

# Filter: both drops must be in the same direction and meaningful
v = v[v['bg_drop_pred'].abs() >= MIN_BG_DROP]
v['ratio'] = v['bg_drop_actual'] / v['bg_drop_pred']
v = v[(v['ratio'] > 0) & (v['ratio'] < 5)]            # same direction, not pathological

# Filter: exclude samples where BG was rising sharply (possible missed carbs)
# Rising >0.5 mmol/L (9 mg/dL) over 2h when loop predicted a drop → suspicious
v = v[~((v['bg_drop_pred'] > 0) & (v['bg_drop_actual'] < -9))]

print(f"  Valid 2h samples after filtering: {len(v):,}")

# ── 3. Back-calculate TDD_eff per sample ──────────────────────────────────────
# ISF that would have zeroed the prediction error
v['isf_eff'] = v['isf_v1'] * v['ratio']

# TDD that produces ISF_eff via the 7D-TDD formula structure
# ISF_7dtdd = (1700/TDD) × ln_target / ln_bg
# → TDD = 1700 × ln_target / (ISF × ln_bg)
v['tdd_eff_sample'] = (1700.0 * LN_TARGET) / (v['isf_eff'] * v['ln_bg'])

# Sanity bound: TDD must be physiologically plausible
v = v[(v['tdd_eff_sample'] > 3) & (v['tdd_eff_sample'] < 120)]

print(f"  Valid samples after TDD bounds: {len(v):,}")

# ── 4. Nightly aggregation ────────────────────────────────────────────────────
nightly = v.groupby('date').agg(
    n_samples       = ('tdd_eff_sample', 'count'),
    tdd_implied_night = ('tdd_eff_sample', 'median'),
    tdd_eff_iqr25   = ('tdd_eff_sample', lambda x: x.quantile(0.25)),
    tdd_eff_iqr75   = ('tdd_eff_sample', lambda x: x.quantile(0.75)),
    tdd_7day_med    = ('tdd_7day', 'median'),
    ratio_med       = ('ratio', 'median'),
    ratio_iqr25     = ('ratio', lambda x: x.quantile(0.25)),
    ratio_iqr75     = ('ratio', lambda x: x.quantile(0.75)),
    isf_eff_med     = ('isf_eff', 'median'),
    isf_7dtdd_med   = ('isf_7dtdd', 'median'),
    bg_med          = ('bg', 'median'),
).reset_index().sort_values('date')

# Require minimum samples per night
MIN_SAMPLES = 5
nightly = nightly[nightly['n_samples'] >= MIN_SAMPLES]
print(f"  Nights with ≥{MIN_SAMPLES} valid samples: {len(nightly)}")

# ── 5. Rolling TDD_effective with exponential smoothing ───────────────────────
for alpha in ALPHAS:
    col = f'tdd_effective_a{int(alpha*100):02d}'
    vals = []
    tdd_eff_prev = None
    for _, row in nightly.iterrows():
        if tdd_eff_prev is None:
            # Bootstrap: use first night's tdd_7day as the prior
            tdd_eff_prev = row['tdd_7day_med']
        tdd_eff_new = alpha * row['tdd_implied_night'] + (1 - alpha) * tdd_eff_prev
        vals.append(tdd_eff_new)
        tdd_eff_prev = tdd_eff_new
    nightly[col] = vals

# Shift forward by 1 day — tonight's regression informs TOMORROW's TDD_effective
# This prevents look-ahead bias: the TDD_effective available on day N
# was computed from nights 1..N-1
for alpha in ALPHAS:
    col = f'tdd_effective_a{int(alpha*100):02d}'
    col_next = f'tdd_eff_nextday_a{int(alpha*100):02d}'
    nightly[col_next] = nightly[col].shift(1)

print("\n" + "═" * 72)
print("NIGHTLY TDD_effective REGRESSION")
print("═" * 72)

# ── 6. Per-night summary ─────────────────────────────────────────────────────
alpha_col = f'tdd_effective_a{int(ALPHA_DEFAULT*100):02d}'
alpha_next = f'tdd_eff_nextday_a{int(ALPHA_DEFAULT*100):02d}'

print(f"\n{'Date':>12s}  {'n':>3s}  {'TDD_7day':>8s}  {'TDD_impl':>8s}  {'TDD_eff':>8s}  "
      f"{'Δ eff-7d':>8s}  {'ratio':>7s}  {'BG_med':>6s}")
print("─" * 80)

for _, row in nightly.iterrows():
    tdd_eff = row[alpha_col]
    delta = tdd_eff - row['tdd_7day_med'] if pd.notna(row['tdd_7day_med']) else float('nan')
    print(f"{str(row['date']):>12s}  {row['n_samples']:3.0f}  "
          f"{row['tdd_7day_med']:8.1f}  {row['tdd_implied_night']:8.1f}  "
          f"{tdd_eff:8.1f}  {delta:+8.1f}  "
          f"{row['ratio_med']:7.2f}  {row['bg_med']:6.0f}")

# ── 7. Backtest: what if we had used TDD_effective? ───────────────────────────
# Merge TDD_effective (next-day shifted) back to samples
v_bt = v.merge(
    nightly[['date', alpha_next, 'tdd_implied_night']],
    on='date', how='left'
).dropna(subset=[alpha_next])

# Adaptive ISF using TDD_effective
v_bt['isf_adaptive'] = (1700.0 / v_bt[alpha_next]) * LN_TARGET / v_bt['ln_bg']

# Counterfactual predictions
# pred_adaptive = BG - bg_drop_pred × (ISF_adaptive / ISF_v1)
v_bt['pred_adaptive'] = v_bt['bg'] - v_bt['bg_drop_pred'] * (v_bt['isf_adaptive'] / v_bt['isf_v1'])
v_bt['err_adaptive']  = v_bt['actual_bg_2h'] - v_bt['pred_adaptive']

# Same for static 7D-TDD
v_bt['pred_7dtdd'] = v_bt['bg'] - v_bt['bg_drop_pred'] * (v_bt['isf_7dtdd'] / v_bt['isf_v1'])
v_bt['err_7dtdd']  = v_bt['actual_bg_2h'] - v_bt['pred_7dtdd']

# And for v1 actual
v_bt['pred_v1']  = v_bt['pred_iob_24']
v_bt['err_v1']   = v_bt['actual_bg_2h'] - v_bt['pred_v1']

print(f"\n── Backtest Comparison (+2h, α={ALPHA_DEFAULT}) ──────────────────────────")
print(f"  Samples with prior-night TDD_effective available: {len(v_bt):,}")
print(f"\n  {'Formula':<28s}  {'Bias':>6s}  {'MAE':>6s}  {'RMSE':>6s}  {'±1mmol':>7s}  {'±2mmol':>7s}  {'ISF_med':>7s}")
print("  " + "─" * 75)

for label, err_col, isf_col in [
    ('v1 (Actual)',       'err_v1',       'isf_v1'),
    ('7D-TDD (static)',   'err_7dtdd',    'isf_7dtdd'),
    ('7D-TDD (adaptive)', 'err_adaptive', 'isf_adaptive'),
]:
    e = v_bt[err_col].dropna()
    isf_med = v_bt[isf_col].median() if isf_col in v_bt.columns else float('nan')
    print(f"  {label:<28s}  {e.mean():+6.1f}  {e.abs().mean():6.1f}  "
          f"{np.sqrt((e**2).mean()):6.1f}  "
          f"{(e.abs() <= 18).mean()*100:6.1f}%  "
          f"{(e.abs() <= 36).mean()*100:6.1f}%  "
          f"{isf_med:7.1f}")

# ── 8. Per-night backtest comparison ──────────────────────────────────────────
print(f"\n── Per-Night Breakdown (+2h, α={ALPHA_DEFAULT}) " + "─" * 30)
print(f"  {'Date':>12s}  {'n':>3s}  │ {'v1 MAE':>7s}  {'7D MAE':>7s}  {'Adp MAE':>7s}  │ "
      f"{'v1 bias':>7s}  {'7D bias':>7s}  {'Adp bias':>8s}  │ "
      f"{'TDD_7d':>6s}  {'TDD_eff':>7s}")
print("  " + "─" * 100)

for date, grp in v_bt.groupby('date'):
    n = len(grp)
    tdd7 = grp['tdd_7day'].median()
    tdd_e = grp[alpha_next].median()
    print(f"  {str(date):>12s}  {n:3d}  │ "
          f"{grp['err_v1'].abs().mean():7.1f}  "
          f"{grp['err_7dtdd'].abs().mean():7.1f}  "
          f"{grp['err_adaptive'].abs().mean():7.1f}  │ "
          f"{grp['err_v1'].mean():+7.1f}  "
          f"{grp['err_7dtdd'].mean():+7.1f}  "
          f"{grp['err_adaptive'].mean():+8.1f}  │ "
          f"{tdd7:6.1f}  {tdd_e:7.1f}")

# ── 9. Learning rate comparison ───────────────────────────────────────────────
print(f"\n── Learning Rate Comparison (+2h) " + "─" * 40)
print(f"  {'α':>5s}  {'~Window':>7s}  {'Bias':>6s}  {'MAE':>6s}  {'RMSE':>6s}  {'±1mmol':>7s}")
print("  " + "─" * 50)

for alpha in ALPHAS:
    col_next = f'tdd_eff_nextday_a{int(alpha*100):02d}'
    tmp = v.merge(nightly[['date', col_next]], on='date', how='left').dropna(subset=[col_next])
    isf_a = (1700.0 / tmp[col_next]) * LN_TARGET / tmp['ln_bg']
    pred_a = tmp['bg'] - tmp['bg_drop_pred'] * (isf_a / tmp['isf_v1'])
    err_a = tmp['actual_bg_2h'] - pred_a
    e = err_a.dropna()
    window = f"~{1/alpha:.0f}n"
    marker = " ←" if alpha == ALPHA_DEFAULT else ""
    print(f"  {alpha:5.2f}  {window:>7s}  {e.mean():+6.1f}  {e.abs().mean():6.1f}  "
          f"{np.sqrt((e**2).mean()):6.1f}  {(e.abs()<=18).mean()*100:6.1f}%{marker}")

# ── 10. Convergence: how quickly does TDD_effective stabilise? ────────────────
print(f"\n── TDD_effective Convergence (α={ALPHA_DEFAULT}) " + "─" * 30)
print(f"  {'Night':>5s}  {'Date':>12s}  {'TDD_7day':>8s}  {'TDD_impl':>8s}  {'TDD_eff':>8s}  "
      f"{'Δ(night)':>8s}  {'Δ(cumul)':>8s}")
print("  " + "─" * 68)

first_tdd7 = nightly['tdd_7day_med'].iloc[0]
for i, (_, row) in enumerate(nightly.iterrows()):
    tdd_eff = row[alpha_col]
    delta_night = row['tdd_implied_night'] - tdd_eff  # how far tonight's implied was from current eff
    delta_cum = tdd_eff - first_tdd7                   # total drift from initial TDD_7day
    print(f"  {i+1:5d}  {str(row['date']):>12s}  {row['tdd_7day_med']:8.1f}  "
          f"{row['tdd_implied_night']:8.1f}  {tdd_eff:8.1f}  "
          f"{delta_night:+8.1f}  {delta_cum:+8.1f}")

# ── 11. The delta — the behavioural signature ────────────────────────────────
print(f"\n── Sensitivity Offset (TDD_effective − TDD_7day) " + "─" * 25)
deltas = nightly[alpha_col] - nightly['tdd_7day_med']
deltas_valid = deltas.dropna()
print(f"  Mean offset:   {deltas_valid.mean():+.2f} U/day")
print(f"  Median offset: {deltas_valid.median():+.2f} U/day")
print(f"  Std:           {deltas_valid.std():.2f} U/day")
print(f"  Range:         [{deltas_valid.min():+.1f}, {deltas_valid.max():+.1f}] U/day")
if deltas_valid.median() > 0.5:
    print("  → Persistently positive: patient less sensitive than TDD history implies")
elif deltas_valid.median() < -0.5:
    print("  → Persistently negative: patient more sensitive than TDD history implies")
else:
    print("  → Near zero: TDD_7day is a good predictor for this patient")

# ── 12. Figure ────────────────────────────────────────────────────────────────
print("\nGenerating figure...")

BG_C = '#0f0f0f'; PANEL = '#1a1a2e'; GRID = '#2a2a4a'; TXT = '#e0e0ff'
C_V1 = '#4fc3f7'; C_7D = '#ce93d8'; C_ADAPT = '#66bb6a'; C_IMPL = '#ffb74d'

def style(ax, title):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=TXT, labelsize=8)
    ax.set_title(title, color=TXT, fontsize=10, fontweight='bold')
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.grid(True, color=GRID, lw=0.5, ls='--', alpha=0.7)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.xaxis.label.set_fontsize(8); ax.yaxis.label.set_fontsize(8)

fig = plt.figure(figsize=(20, 20))
fig.patch.set_facecolor(BG_C)
gs = gridspec.GridSpec(4, 3, figure=fig, hspace=0.5, wspace=0.38)

dates_num = matplotlib.dates.date2num([pd.Timestamp(d) for d in nightly['date']])
date_labels = [str(d) for d in nightly['date']]

# P1: TDD timeline — 7day vs implied vs effective (full width)
ax1 = fig.add_subplot(gs[0, :])
style(ax1, f'Nightly TDD: 7-Day Rolling vs Implied vs Effective (α={ALPHA_DEFAULT})')
ax1.plot(dates_num, nightly['tdd_7day_med'], 'o-', color=C_7D, lw=1.5, ms=4, label='TDD 7-day')
ax1.plot(dates_num, nightly['tdd_implied_night'], 's-', color=C_IMPL, lw=1.5, ms=4, label='TDD implied (overnight)')
ax1.plot(dates_num, nightly[alpha_col], 'D-', color=C_ADAPT, lw=2, ms=5, label='TDD effective')
# IQR band for implied
ax1.fill_between(dates_num, nightly['tdd_eff_iqr25'], nightly['tdd_eff_iqr75'],
                 alpha=0.15, color=C_IMPL, label='Implied IQR')
ax1.xaxis_date(); ax1.xaxis.set_major_locator(matplotlib.dates.AutoDateLocator())
ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))
ax1.set_ylabel('TDD (U/day)'); ax1.set_xlabel('')
ax1.legend(fontsize=9, labelcolor=TXT, facecolor=PANEL, loc='best')
for tick in ax1.get_xticklabels(): tick.set_rotation(45); tick.set_fontsize(7)

# P2: Sensitivity offset (delta) over time
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'Sensitivity Offset (TDD_eff − TDD_7day)')
deltas_plot = nightly[alpha_col] - nightly['tdd_7day_med']
colors = ['#ef5350' if d > 0 else '#66bb6a' for d in deltas_plot]
ax2.bar(range(len(nightly)), deltas_plot, color=colors, alpha=0.8)
ax2.axhline(0, color='white', lw=0.8, ls='--')
ax2.set_xlabel('Night #'); ax2.set_ylabel('Δ TDD (U/day)')
ax2.annotate('+ve = less sensitive\n−ve = more sensitive', xy=(0.02, 0.98),
             xycoords='axes fraction', va='top', fontsize=7, color=TXT, alpha=0.7)

# P3: Overnight ratio distribution
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'Overnight Sensitivity Ratio Distribution')
ax3.hist(v['ratio'].clip(0, 3), bins=40, color=C_IMPL, alpha=0.7, density=True)
ax3.axvline(1.0, color='white', lw=1.2, ls='--', label='ratio=1 (perfect)')
ax3.axvline(v['ratio'].median(), color=C_ADAPT, lw=1.5, ls='-', label=f"median={v['ratio'].median():.2f}")
ax3.set_xlabel('Actual/Predicted BG drop'); ax3.set_ylabel('Density')
ax3.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P4: Learning rate comparison
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, 'TDD_effective by Learning Rate')
for alpha in ALPHAS:
    col = f'tdd_effective_a{int(alpha*100):02d}'
    lw = 2.5 if alpha == ALPHA_DEFAULT else 1.2
    al = 1.0 if alpha == ALPHA_DEFAULT else 0.6
    ax4.plot(range(len(nightly)), nightly[col], '-', lw=lw, alpha=al, label=f'α={alpha}')
ax4.plot(range(len(nightly)), nightly['tdd_7day_med'], '--', color='white', lw=0.8, alpha=0.5, label='TDD_7day')
ax4.set_xlabel('Night #'); ax4.set_ylabel('TDD_effective (U/day)')
ax4.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P5: Error distributions — v1 vs 7D-TDD static vs adaptive
ax5 = fig.add_subplot(gs[2, 0])
style(ax5, 'Prediction Error Distributions (+2h)')
for label, col, c in [('v1 Actual', 'err_v1', C_V1),
                        ('7D-TDD static', 'err_7dtdd', C_7D),
                        ('7D-TDD adaptive', 'err_adaptive', C_ADAPT)]:
    ax5.hist(v_bt[col].clip(-80, 80), bins=35, alpha=0.5, color=c, density=True, label=label)
ax5.axvline(0, color='white', lw=0.8, ls='--')
ax5.set_xlabel('Pred Error (mg/dL)'); ax5.set_ylabel('Density')
ax5.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P6: MAE comparison bar chart
ax6 = fig.add_subplot(gs[2, 1])
style(ax6, 'MAE Comparison (+2h)')
labels_bar = ['v1\nActual', '7D-TDD\nStatic', '7D-TDD\nAdaptive']
maes = [v_bt['err_v1'].abs().mean(), v_bt['err_7dtdd'].abs().mean(), v_bt['err_adaptive'].abs().mean()]
cols_bar = [C_V1, C_7D, C_ADAPT]
bars = ax6.bar(range(3), maes, color=cols_bar, alpha=0.85)
for bar, val in zip(bars, maes):
    ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{val:.1f}', ha='center', fontsize=9, color=TXT)
ax6.set_xticks(range(3)); ax6.set_xticklabels(labels_bar, fontsize=8)
ax6.set_ylabel('MAE (mg/dL)')

# P7: Bias comparison
ax7 = fig.add_subplot(gs[2, 2])
style(ax7, 'Prediction Bias (+2h)')
biases = [v_bt['err_v1'].mean(), v_bt['err_7dtdd'].mean(), v_bt['err_adaptive'].mean()]
bars = ax7.bar(range(3), biases, color=cols_bar, alpha=0.85)
for bar, val in zip(bars, biases):
    ax7.text(bar.get_x() + bar.get_width()/2,
             bar.get_height() + (0.3 if val >= 0 else -1.5),
             f'{val:+.1f}', ha='center', fontsize=9, color=TXT)
ax7.axhline(0, color='white', lw=0.8, ls='--')
ax7.set_xticks(range(3)); ax7.set_xticklabels(labels_bar, fontsize=8)
ax7.set_ylabel('Mean Error (mg/dL)')

# P8: Per-night MAE comparison — line chart
ax8 = fig.add_subplot(gs[3, 0:2])
style(ax8, 'Per-Night MAE: Static 7D-TDD vs Adaptive (+2h)')
night_dates = []
night_mae_static = []
night_mae_adapt = []
night_mae_v1 = []
for date, grp in v_bt.groupby('date'):
    night_dates.append(date)
    night_mae_v1.append(grp['err_v1'].abs().mean())
    night_mae_static.append(grp['err_7dtdd'].abs().mean())
    night_mae_adapt.append(grp['err_adaptive'].abs().mean())
night_dates_num = matplotlib.dates.date2num([pd.Timestamp(d) for d in night_dates])
ax8.plot(night_dates_num, night_mae_v1, 'o-', color=C_V1, lw=1.5, ms=4, label='v1 Actual')
ax8.plot(night_dates_num, night_mae_static, 's-', color=C_7D, lw=1.5, ms=4, label='7D-TDD static')
ax8.plot(night_dates_num, night_mae_adapt, 'D-', color=C_ADAPT, lw=2, ms=5, label='7D-TDD adaptive')
ax8.xaxis_date(); ax8.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))
ax8.set_ylabel('MAE (mg/dL)'); ax8.set_xlabel('')
ax8.legend(fontsize=9, labelcolor=TXT, facecolor=PANEL)
for tick in ax8.get_xticklabels(): tick.set_rotation(45); tick.set_fontsize(7)

# P9: ISF comparison scatter
ax9 = fig.add_subplot(gs[3, 2])
style(ax9, 'ISF: Static vs Adaptive 7D-TDD')
ax9.scatter(v_bt['isf_7dtdd'], v_bt['isf_adaptive'], alpha=0.3, s=8, color=C_ADAPT, edgecolors='none')
lims = [v_bt[['isf_7dtdd', 'isf_adaptive']].min().min(),
        v_bt[['isf_7dtdd', 'isf_adaptive']].max().max()]
ax9.plot(lims, lims, '--', color='white', lw=0.8, alpha=0.5)
ax9.set_xlabel('ISF 7D-TDD static'); ax9.set_ylabel('ISF 7D-TDD adaptive')

fig.suptitle(f'TDD_effective: Overnight Sensitivity Regression\n'
             f'7D-TDD formula with nightly adaptive correction  |  α={ALPHA_DEFAULT}  |  +2h horizon',
             color=TXT, fontsize=13, fontweight='bold', y=0.995)

plt.savefig('ns_tdd_effective_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_tdd_effective_results.png")

# ── 13. Save outputs ──────────────────────────────────────────────────────────
nightly.to_csv('ns_tdd_effective_nightly.csv', index=False)

# Summary text
lines = [
    "TDD_effective — OVERNIGHT SENSITIVITY REGRESSION",
    "=" * 62,
    f"Input:           ns_backtest_overnight.csv",
    f"Overnight rows:  {len(df):,}",
    f"Valid 2h samples:{len(v):,}  (after COB=0, ratio, BG bounds, IOB delta filters)",
    f"Nights used:     {len(nightly)} (≥{MIN_SAMPLES} samples each)",
    f"Learning rate:   α={ALPHA_DEFAULT}  (~{1/ALPHA_DEFAULT:.0f}-night effective window)",
    f"Formula:         ISF = (1700/TDD_eff) × ln({TARGET:.0f}/D+1) / ln(BG/D+1)",
    "",
    "Concept:",
    "  TDD_7day is a description of dose history.",
    "  TDD_effective is the TDD that WOULD HAVE produced correct overnight predictions.",
    "  Updated at 07:00 each morning from last night's 00:00–07:00 data.",
    "  The delta (TDD_eff − TDD_7day) is the patient's behavioural sensitivity signature.",
    "",
    f"Overall sensitivity ratio: median={v['ratio'].median():.3f}  "
    f"IQR=[{v['ratio'].quantile(.25):.3f}, {v['ratio'].quantile(.75):.3f}]",
    f"Sensitivity offset:        median={deltas_valid.median():+.2f} U/day  "
    f"mean={deltas_valid.mean():+.2f} U/day",
    "",
]

lines.append("Nightly detail:")
lines.append(f"  {'Date':>12s}  {'n':>3s}  {'TDD_7day':>8s}  {'TDD_impl':>8s}  {'TDD_eff':>8s}  {'Δ':>7s}  {'ratio':>7s}")
lines.append("  " + "─" * 60)
for _, row in nightly.iterrows():
    tdd_eff = row[alpha_col]
    delta = tdd_eff - row['tdd_7day_med'] if pd.notna(row['tdd_7day_med']) else float('nan')
    lines.append(f"  {str(row['date']):>12s}  {row['n_samples']:3.0f}  "
                 f"{row['tdd_7day_med']:8.1f}  {row['tdd_implied_night']:8.1f}  "
                 f"{tdd_eff:8.1f}  {delta:+7.1f}  {row['ratio_med']:7.2f}")

lines.append("")
lines.append("Backtest comparison (+2h):")
lines.append(f"  {'Formula':<28s}  {'Bias':>6s}  {'MAE':>6s}  {'RMSE':>6s}  {'±1mmol':>7s}")
lines.append("  " + "─" * 55)
for label, err_col in [('v1 (Actual)', 'err_v1'),
                        ('7D-TDD (static)', 'err_7dtdd'),
                        ('7D-TDD (adaptive)', 'err_adaptive')]:
    e = v_bt[err_col].dropna()
    lines.append(f"  {label:<28s}  {e.mean():+6.1f}  {e.abs().mean():6.1f}  "
                 f"{np.sqrt((e**2).mean()):6.1f}  {(e.abs()<=18).mean()*100:6.1f}%")

with open('ns_tdd_effective_summary.txt', 'w') as f:
    f.write('\n'.join(lines))

print('\n' + '\n'.join(lines))
print("\nSaved: ns_tdd_effective_nightly.csv, ns_tdd_effective_summary.txt")
