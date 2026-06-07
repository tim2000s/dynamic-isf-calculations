"""
k-Sensitivity Analysis for Power-Law DynamicISF BG Scaling
===========================================================
Sweeps k from 0.5 to 3.0 for ISF = (1700/TDD_7day) * (99/BG)^k
and reports per-BG-band accuracy metrics to find the best compromise.

Outputs
-------
  ns_bg_scaling_k_sensitivity.png
  ns_bg_scaling_k_sensitivity_summary.txt
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ── Constants ────────────────────────────────────────────────────────────────
D      = 82.0
TARGET = 99.0
LN_TARGET = np.log(TARGET / D + 1)
MIN_BG_DROP = 3.0
BG_LO, BG_HI = 72, 200

# Dark theme
BG_C   = '#0f0f0f'
PANEL  = '#1a1a2e'
GRID   = '#2a2a4a'
TXT    = '#e0e0ff'

# BG bands
BANDS = [
    ('<90',     72,  90),
    ('90-105',  90, 105),
    ('105-120', 105, 120),
    ('120-150', 120, 150),
    ('>150',    150, 200),
]
BAND_COLOURS = ['#00ccff', '#22dd88', '#ffcc00', '#ff8844', '#ff4466']

K_VALS = np.arange(0.5, 3.05, 0.1)

# ── 1. Load and filter (same as ns_bg_scaling.py) ───────────────────────────
print("Loading ns_backtest_overnight.csv ...")
df = pd.read_csv('ns_backtest_overnight.csv')
df['ts'] = pd.to_datetime(df['ts'], format='ISO8601')
df['date'] = pd.to_datetime(df['date']).dt.date
df['hour'] = df['ts'].dt.hour
print(f"  Total overnight rows: {len(df):,}")

v = df.dropna(subset=['pred_iob_24', 'actual_bg_2h', 'isf_v1', 'tdd_7day']).copy()
v = v[(v['bg'] >= BG_LO) & (v['bg'] <= BG_HI)]

v['bg_drop_pred']   = v['bg'] - v['pred_iob_24']
v['bg_drop_actual']  = v['bg'] - v['actual_bg_2h']

v = v[v['bg_drop_pred'].abs() >= MIN_BG_DROP]
v['ratio'] = v['bg_drop_actual'] / v['bg_drop_pred']
v = v[(v['ratio'] > 0) & (v['ratio'] < 5)]
v = v[~((v['bg_drop_pred'] > 0) & (v['bg_drop_actual'] < -9))]

v['isf_eff'] = v['isf_v1'] * v['ratio']
v = v[(v['isf_eff'] > 10) & (v['isf_eff'] < 800)]

print(f"  Valid samples: {len(v):,}")

# Pre-extract arrays for speed
bg_vals   = v['bg'].values
tdd_vals  = v['tdd_7day'].values
isf_v1    = v['isf_v1'].values
bg_drop   = v['bg_drop_pred'].values
actual_2h = v['actual_bg_2h'].values

# Band masks
band_masks = {}
for name, lo, hi in BANDS:
    band_masks[name] = (bg_vals >= lo) & (bg_vals < hi)

# ── 2. Sweep k values ───────────────────────────────────────────────────────
print(f"\nSweeping k = {K_VALS[0]:.1f} to {K_VALS[-1]:.1f} ...")

overall = {m: [] for m in ['mae', 'bias', 'rmse', 'w18']}
band_mae  = {b: [] for b, _, _ in BANDS}
band_bias = {b: [] for b, _, _ in BANDS}

for k in K_VALS:
    isf_cand = (1700.0 / tdd_vals) * (TARGET / bg_vals) ** k
    pred = bg_vals - bg_drop * (isf_cand / isf_v1)
    err  = actual_2h - pred

    overall['mae'].append(np.abs(err).mean())
    overall['bias'].append(err.mean())
    overall['rmse'].append(np.sqrt((err ** 2).mean()))
    overall['w18'].append((np.abs(err) <= 18).mean() * 100)

    for name, _, _ in BANDS:
        m = band_masks[name]
        if m.sum() > 0:
            band_mae[name].append(np.abs(err[m]).mean())
            band_bias[name].append(err[m].mean())
        else:
            band_mae[name].append(np.nan)
            band_bias[name].append(np.nan)

# Convert to arrays
for key in overall:
    overall[key] = np.array(overall[key])
for name in band_mae:
    band_mae[name]  = np.array(band_mae[name])
    band_bias[name] = np.array(band_bias[name])

# BG-band bias range per k
bias_range = np.array([
    max(band_bias[b][i] for b, _, _ in BANDS if not np.isnan(band_bias[b][i]))
    - min(band_bias[b][i] for b, _, _ in BANDS if not np.isnan(band_bias[b][i]))
    for i in range(len(K_VALS))
])

# ── 3. Find optimal k values ────────────────────────────────────────────────
idx_mae   = np.argmin(overall['mae'])
idx_bias  = np.argmin(np.abs(overall['bias']))
idx_range = np.argmin(bias_range)

# Best compromise: lowest MAE where band bias range < 10 mg/dL
mask_flat = bias_range < 10
if mask_flat.any():
    sub_mae = overall['mae'].copy()
    sub_mae[~mask_flat] = 999
    idx_comp = np.argmin(sub_mae)
else:
    # Relax: lowest MAE where range < 15
    mask_flat2 = bias_range < 15
    if mask_flat2.any():
        sub_mae = overall['mae'].copy()
        sub_mae[~mask_flat2] = 999
        idx_comp = np.argmin(sub_mae)
    else:
        idx_comp = idx_mae  # fallback

k_mae  = K_VALS[idx_mae]
k_bias = K_VALS[idx_bias]
k_range = K_VALS[idx_range]
k_comp  = K_VALS[idx_comp]

# ── 4. Current formula metrics (for comparison) ─────────────────────────────
isf_curr = (1700.0 / tdd_vals) * LN_TARGET / np.log(bg_vals / D + 1)
pred_curr = bg_vals - bg_drop * (isf_curr / isf_v1)
err_curr  = actual_2h - pred_curr
curr_mae  = np.abs(err_curr).mean()
curr_bias = err_curr.mean()
curr_rmse = np.sqrt((err_curr ** 2).mean())
curr_w18  = (np.abs(err_curr) <= 18).mean() * 100

# ── 5. Build summary ────────────────────────────────────────────────────────
lines = []
lines.append("=" * 72)
lines.append("k-SENSITIVITY ANALYSIS: ISF = (1700/TDD) * (99/BG)^k")
lines.append("=" * 72)
lines.append(f"  Valid 2h samples: {len(v):,}")
lines.append(f"  BG range: {BG_LO}-{BG_HI} mg/dL, min predicted drop: {MIN_BG_DROP} mg/dL")
lines.append(f"  k sweep: {K_VALS[0]:.1f} to {K_VALS[-1]:.1f}, step 0.1")
lines.append("")

lines.append("CURRENT FORMULA  ln(target/D+1) / ln(BG/D+1), D=82")
lines.append(f"  MAE={curr_mae:.2f}  bias={curr_bias:+.2f}  RMSE={curr_rmse:.2f}  ±1mmol={curr_w18:.1f}%")
lines.append("")

lines.append("OPTIMAL k VALUES")
lines.append(f"  Best MAE:          k={k_mae:.1f}   MAE={overall['mae'][idx_mae]:.2f}  "
             f"bias={overall['bias'][idx_mae]:+.2f}  band-range={bias_range[idx_mae]:.1f}")
lines.append(f"  Best |bias|:       k={k_bias:.1f}   MAE={overall['mae'][idx_bias]:.2f}  "
             f"bias={overall['bias'][idx_bias]:+.2f}  band-range={bias_range[idx_bias]:.1f}")
lines.append(f"  Best band range:   k={k_range:.1f}   MAE={overall['mae'][idx_range]:.2f}  "
             f"bias={overall['bias'][idx_range]:+.2f}  band-range={bias_range[idx_range]:.1f}")
lines.append(f"  Compromise:        k={k_comp:.1f}   MAE={overall['mae'][idx_comp]:.2f}  "
             f"bias={overall['bias'][idx_comp]:+.2f}  band-range={bias_range[idx_comp]:.1f}")
lines.append("")

lines.append("PER-BG-BAND METRICS AT EACH OPTIMAL k")
lines.append(f"  {'Band':<10s}  {'n':>5s}  " +
             "  ".join(f"bias@k={K_VALS[i]:.1f}" for i in [idx_mae, idx_bias, idx_range, idx_comp]))
for bi, (name, lo, hi) in enumerate(BANDS):
    n = band_masks[name].sum()
    vals = []
    for i in [idx_mae, idx_bias, idx_range, idx_comp]:
        vals.append(f"{band_bias[name][i]:+10.2f}")
    lines.append(f"  {name:<10s}  {n:5d}  " + "  ".join(vals))
lines.append("")

lines.append("FULL k SWEEP TABLE")
lines.append(f"  {'k':>4s}  {'MAE':>6s}  {'bias':>7s}  {'RMSE':>6s}  {'±1mmol':>6s}  {'BandRng':>7s}  " +
             "  ".join(f"{'b_'+b:<10s}" for b, _, _ in BANDS))
for i, k in enumerate(K_VALS):
    band_str = "  ".join(f"{band_bias[b][i]:+10.2f}" for b, _, _ in BANDS)
    lines.append(f"  {k:4.1f}  {overall['mae'][i]:6.2f}  {overall['bias'][i]:+7.2f}  "
                 f"{overall['rmse'][i]:6.2f}  {overall['w18'][i]:6.1f}  {bias_range[i]:7.1f}  {band_str}")

summary = "\n".join(lines)
print(summary)
with open('ns_bg_scaling_k_sensitivity_summary.txt', 'w') as f:
    f.write(summary)
print("\nSaved ns_bg_scaling_k_sensitivity_summary.txt")

# ── 6. Plot ──────────────────────────────────────────────────────────────────
print("Generating figure ...")

fig = plt.figure(figsize=(20, 14), facecolor=BG_C)
gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.30,
                       left=0.06, right=0.97, top=0.93, bottom=0.06)

def style_ax(ax, title):
    ax.set_facecolor(PANEL)
    ax.set_title(title, color=TXT, fontsize=12, fontweight='bold', pad=8)
    ax.tick_params(colors=TXT, labelsize=9)
    for sp in ax.spines.values():
        sp.set_color(GRID)
    ax.grid(True, color=GRID, alpha=0.4, linewidth=0.5)

# Mark helpers
def add_k_markers(ax):
    for ki, lab, col, ms in [(k_mae, 'MAE', '#00ccff', 'v'),
                              (k_bias, '|bias|', '#22dd88', 's'),
                              (k_range, 'range', '#ffcc00', 'D'),
                              (k_comp, 'comp', '#ff4466', 'o')]:
        ax.axvline(ki, color=col, alpha=0.3, linewidth=1, linestyle='--')

# Panel 1: MAE, bias, RMSE vs k
ax1 = fig.add_subplot(gs[0, 0])
style_ax(ax1, 'Overall Metrics vs k')
ax1.plot(K_VALS, overall['mae'],  color='#00ccff', linewidth=2, label='MAE')
ax1.plot(K_VALS, overall['bias'], color='#22dd88', linewidth=2, label='Bias')
ax1.plot(K_VALS, overall['rmse'], color='#ff8844', linewidth=2, label='RMSE')
ax1.axhline(0, color=TXT, alpha=0.3, linewidth=0.8)
add_k_markers(ax1)
ax1.set_xlabel('k', color=TXT, fontsize=10)
ax1.set_ylabel('mg/dL', color=TXT, fontsize=10)
ax1.legend(facecolor=PANEL, edgecolor=GRID, labelcolor=TXT, fontsize=9)

# Panel 2: Per-BG-band bias vs k
ax2 = fig.add_subplot(gs[0, 1])
style_ax(ax2, 'Per-BG-Band Bias vs k')
for bi, (name, _, _) in enumerate(BANDS):
    ax2.plot(K_VALS, band_bias[name], color=BAND_COLOURS[bi], linewidth=1.8, label=name)
ax2.axhline(0, color=TXT, alpha=0.3, linewidth=0.8)
add_k_markers(ax2)
ax2.set_xlabel('k', color=TXT, fontsize=10)
ax2.set_ylabel('Bias (mg/dL)', color=TXT, fontsize=10)
ax2.legend(facecolor=PANEL, edgecolor=GRID, labelcolor=TXT, fontsize=9, title='BG band',
           title_fontproperties={'weight': 'bold'})

# Panel 3: Per-BG-band MAE vs k
ax3 = fig.add_subplot(gs[0, 2])
style_ax(ax3, 'Per-BG-Band MAE vs k')
for bi, (name, _, _) in enumerate(BANDS):
    ax3.plot(K_VALS, band_mae[name], color=BAND_COLOURS[bi], linewidth=1.8, label=name)
add_k_markers(ax3)
ax3.set_xlabel('k', color=TXT, fontsize=10)
ax3.set_ylabel('MAE (mg/dL)', color=TXT, fontsize=10)
ax3.legend(facecolor=PANEL, edgecolor=GRID, labelcolor=TXT, fontsize=9, title='BG band',
           title_fontproperties={'weight': 'bold'})

# Panel 4: Band bias range vs k
ax4 = fig.add_subplot(gs[1, 0])
style_ax(ax4, 'BG-Band Bias Range vs k')
ax4.plot(K_VALS, bias_range, color='#ff4466', linewidth=2)
ax4.axhline(10, color='#ffcc00', alpha=0.5, linewidth=1, linestyle=':', label='10 mg/dL threshold')
# Mark optima
ax4.plot(k_mae,   bias_range[idx_mae],   'v', color='#00ccff', markersize=10, label=f'Best MAE k={k_mae:.1f}')
ax4.plot(k_bias,  bias_range[idx_bias],  's', color='#22dd88', markersize=10, label=f'Best |bias| k={k_bias:.1f}')
ax4.plot(k_range, bias_range[idx_range], 'D', color='#ffcc00', markersize=10, label=f'Best range k={k_range:.1f}')
ax4.plot(k_comp,  bias_range[idx_comp],  'o', color='#ff4466', markersize=10, label=f'Compromise k={k_comp:.1f}')
ax4.set_xlabel('k', color=TXT, fontsize=10)
ax4.set_ylabel('Band bias range (mg/dL)', color=TXT, fontsize=10)
ax4.legend(facecolor=PANEL, edgecolor=GRID, labelcolor=TXT, fontsize=8)

# Panel 5: ISF curves at selected k values + current formula
ax5 = fig.add_subplot(gs[1, 1])
style_ax(ax5, 'ISF Curves at Selected k Values')
bg_x = np.linspace(75, 200, 200)
tdd_ref = np.median(tdd_vals)  # Use median TDD for display

# Current formula
isf_c = (1700.0 / tdd_ref) * LN_TARGET / np.log(bg_x / D + 1)
ax5.plot(bg_x, isf_c, color='#888888', linewidth=2, linestyle='--', label=f'Current ln (D={D:.0f})')

for ki, lab, col in [(k_mae, f'Best MAE k={k_mae:.1f}', '#00ccff'),
                      (k_bias, f'Best |bias| k={k_bias:.1f}', '#22dd88'),
                      (k_range, f'Best range k={k_range:.1f}', '#ffcc00'),
                      (k_comp, f'Compromise k={k_comp:.1f}', '#ff4466')]:
    isf_k = (1700.0 / tdd_ref) * (TARGET / bg_x) ** ki
    ax5.plot(bg_x, isf_k, color=col, linewidth=1.8, label=lab)

ax5.axvline(TARGET, color=TXT, alpha=0.2, linewidth=0.8, linestyle=':')
ax5.set_xlabel('BG (mg/dL)', color=TXT, fontsize=10)
ax5.set_ylabel(f'ISF (mg/dL per U) @ TDD={tdd_ref:.0f}', color=TXT, fontsize=10)
ax5.legend(facecolor=PANEL, edgecolor=GRID, labelcolor=TXT, fontsize=8)

# Panel 6: ±1mmol% vs k
ax6 = fig.add_subplot(gs[1, 2])
style_ax(ax6, '±1 mmol/L Accuracy vs k')
ax6.plot(K_VALS, overall['w18'], color='#aa66ff', linewidth=2)
add_k_markers(ax6)
# Mark current formula w18
ax6.axhline(curr_w18, color='#888888', alpha=0.5, linewidth=1, linestyle=':', label=f'Current: {curr_w18:.1f}%')
ax6.plot(k_mae,  overall['w18'][idx_mae],  'v', color='#00ccff', markersize=10, label=f'Best MAE: {overall["w18"][idx_mae]:.1f}%')
ax6.plot(k_comp, overall['w18'][idx_comp], 'o', color='#ff4466', markersize=10, label=f'Compromise: {overall["w18"][idx_comp]:.1f}%')
ax6.set_xlabel('k', color=TXT, fontsize=10)
ax6.set_ylabel('±1 mmol/L accuracy (%)', color=TXT, fontsize=10)
ax6.legend(facecolor=PANEL, edgecolor=GRID, labelcolor=TXT, fontsize=8)

fig.suptitle('Power-Law k Sensitivity:  ISF = (1700 / TDD) × (99 / BG)^k',
             color=TXT, fontsize=15, fontweight='bold')
fig.savefig('ns_bg_scaling_k_sensitivity.png', dpi=180, facecolor=BG_C)
plt.close(fig)
print("Saved ns_bg_scaling_k_sensitivity.png")
print("Done.")
