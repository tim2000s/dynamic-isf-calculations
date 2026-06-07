"""
TDD Inflation Analysis
======================
Tests the argument: "Overestimated basal profiles inflate TDD,
causing TDD-based ISF formulas to produce values that are too low."

Three tests:
  1. Is 1700/TDD correct at target BG? (isolates TDD from BG scaling)
  2. Basal IOB analysis — is basal contributing excess IOB overnight?
  3. TDD sensitivity simulation — how does formula accuracy degrade
     with artificially inflated/deflated TDD?

Plus:
  4. What constant (instead of 1700) would be correct for this patient?
  5. Combined: power-law with corrected constant vs TDD inflation

Outputs
-------
  ns_tdd_inflation_results.png
  ns_tdd_inflation_summary.txt
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import minimize_scalar

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

v['isf_eff'] = v['isf_v1'] * v['ratio']
v = v[(v['isf_eff'] > 10) & (v['isf_eff'] < 800)]

print(f"  Valid samples: {len(v):,}")

# ── 2. TEST 1: Is 1700/TDD correct at target BG? ─────────────────────────────
print("\n" + "═" * 72)
print("TEST 1: IS 1700/TDD CORRECT AT TARGET BG?")
print("═" * 72)
print(f"  At target BG ({TARGET} mg/dL), all BG scaling cancels out.")
print(f"  ISF should = 1700/TDD if the 1700 constant and TDD are correct.")

# Samples near target BG (±10 mg/dL)
near_target = v[(v['bg'] >= TARGET - 10) & (v['bg'] <= TARGET + 10)].copy()
print(f"\n  Samples with BG {TARGET-10:.0f}–{TARGET+10:.0f}: {len(near_target):,}")

near_target['isf_1700_tdd'] = 1700.0 / near_target['tdd_7day']
near_target['isf_ratio_at_target'] = near_target['isf_eff'] / near_target['isf_1700_tdd']

print(f"\n  At BG ≈ target:")
print(f"    ISF_implied (what actually worked):    median = {near_target['isf_eff'].median():.1f} mg/dL/U")
print(f"    ISF from 1700/TDD:                     median = {near_target['isf_1700_tdd'].median():.1f} mg/dL/U")
print(f"    Ratio (implied / 1700÷TDD):            median = {near_target['isf_ratio_at_target'].median():.3f}")
print(f"                                           IQR = [{near_target['isf_ratio_at_target'].quantile(.25):.3f}, "
      f"{near_target['isf_ratio_at_target'].quantile(.75):.3f}]")

if near_target['isf_ratio_at_target'].median() > 1.05:
    print(f"\n  → ISF_implied > 1700/TDD by {(near_target['isf_ratio_at_target'].median()-1)*100:.0f}%")
    print(f"    This means 1700/TDD produces ISF values that are TOO LOW at target BG.")
    print(f"    Possible causes: TDD is inflated, OR 1700 is too small for this patient.")
elif near_target['isf_ratio_at_target'].median() < 0.95:
    print(f"\n  → ISF_implied < 1700/TDD by {(1-near_target['isf_ratio_at_target'].median())*100:.0f}%")
    print(f"    1700/TDD produces ISF values that are too HIGH at target BG.")
else:
    print(f"\n  → 1700/TDD is well-calibrated at target BG (within ±5%).")

# Wider BG bands to see if the problem is at target or away from it
print(f"\n  Ratio by BG proximity to target:")
for lo, hi, label in [(72, 85, 'BG 72-85 (well below)'),
                       (85, 95, 'BG 85-95 (below)'),
                       (95, 105, 'BG 95-105 (at target)'),
                       (105, 120, 'BG 105-120 (above)'),
                       (120, 150, 'BG 120-150 (high)'),
                       (150, 200, 'BG 150-200 (very high)')]:
    mask = (v['bg'] >= lo) & (v['bg'] < hi)
    if mask.sum() < 10:
        continue
    grp = v[mask]
    isf_tdd = 1700.0 / grp['tdd_7day']
    ratio = grp['isf_eff'] / isf_tdd
    print(f"    {label:30s}  n={mask.sum():4d}  ISF_eff/1700÷TDD = {ratio.median():.3f}  "
          f"[{ratio.quantile(.25):.3f}, {ratio.quantile(.75):.3f}]")

# ── 3. What constant SHOULD it be? ───────────────────────────────────────────
print("\n" + "═" * 72)
print("TEST 1b: WHAT CONSTANT FITS THIS PATIENT?")
print("═" * 72)

# At target BG, ISF = C/TDD. What C makes ISF = ISF_implied?
# C = ISF_implied × TDD
near_target['implied_constant'] = near_target['isf_eff'] * near_target['tdd_7day']
c_median = near_target['implied_constant'].median()
c_iqr25 = near_target['implied_constant'].quantile(0.25)
c_iqr75 = near_target['implied_constant'].quantile(0.75)

print(f"\n  At BG ≈ target, implied constant C where ISF = C/TDD:")
print(f"    median = {c_median:.0f}  IQR = [{c_iqr25:.0f}, {c_iqr75:.0f}]")
print(f"    (Current formula uses 1700; historical alternatives: 1500, 1700, 1800)")

# Also fit across ALL BG levels using the power-law formula
# ISF_eff = (C/TDD) × (target/BG)^k → C = ISF_eff × TDD / (target/BG)^k
for k in [1.5, 2.0, 2.3, 2.5]:
    v[f'implied_C_k{k}'] = v['isf_eff'] * v['tdd_7day'] / (TARGET / v['bg']) ** k
    c = v[f'implied_C_k{k}']
    c_clean = c[(c > 500) & (c < 5000)]
    print(f"    Power-law k={k}: implied C = {c_clean.median():.0f}  "
          f"IQR = [{c_clean.quantile(.25):.0f}, {c_clean.quantile(.75):.0f}]  "
          f"(n={len(c_clean):,})")

# ── 4. TEST 2: Basal IOB analysis ────────────────────────────────────────────
print("\n" + "═" * 72)
print("TEST 2: BASAL IOB ANALYSIS — IS BASAL OVERESTIMATED?")
print("═" * 72)

# iob_basal: the IOB attributable to basal insulin
# In a perfectly profiled patient, overnight basal IOB should be near zero
# (basal matches liver glucose output, so net effect is neutral)
# Positive iob_basal overnight means the loop is delivering MORE than profile
# basal (suggesting profile basal is already too high, or the loop is
# topping up), Negative means less.

has_basal = df.dropna(subset=['iob_basal', 'iob']).copy()
has_basal_on = has_basal[has_basal['hour'] < 8]

if len(has_basal_on) > 0:
    print(f"\n  Overnight samples with basal IOB data: {len(has_basal_on):,}")
    print(f"\n  IOB breakdown (overnight, all samples):")
    print(f"    Total IOB:     median = {has_basal_on['iob'].median():.3f} U  "
          f"IQR = [{has_basal_on['iob'].quantile(.25):.3f}, {has_basal_on['iob'].quantile(.75):.3f}]")
    print(f"    Basal IOB:     median = {has_basal_on['iob_basal'].median():.3f} U  "
          f"IQR = [{has_basal_on['iob_basal'].quantile(.25):.3f}, {has_basal_on['iob_basal'].quantile(.75):.3f}]")

    bolus_iob = has_basal_on['iob'] - has_basal_on['iob_basal']
    print(f"    Bolus IOB:     median = {bolus_iob.median():.3f} U  "
          f"IQR = [{bolus_iob.quantile(.25):.3f}, {bolus_iob.quantile(.75):.3f}]")

    # Basal as fraction of total IOB
    frac = has_basal_on['iob_basal'] / has_basal_on['iob'].replace(0, np.nan)
    frac_clean = frac.dropna()
    frac_clean = frac_clean[(frac_clean > -10) & (frac_clean < 10)]  # clip outliers
    print(f"    Basal/Total:   median = {frac_clean.median():.1%}  "
          f"IQR = [{frac_clean.quantile(.25):.1%}, {frac_clean.quantile(.75):.1%}]")

    # By hour
    print(f"\n  Basal IOB by hour of night:")
    print(f"    {'Hour':>4s}  {'n':>5s}  {'IOB_basal':>9s}  {'IOB_total':>9s}  {'Basal%':>7s}")
    print(f"    " + "─" * 40)
    for hour in range(8):
        h = has_basal_on[has_basal_on['hour'] == hour]
        if len(h) < 10:
            continue
        bf = h['iob_basal'] / h['iob'].replace(0, np.nan)
        bf_clean = bf.dropna()
        print(f"    {hour:4d}  {len(h):5d}  {h['iob_basal'].median():9.3f}  "
              f"{h['iob'].median():9.3f}  {bf_clean.median():7.1%}")

    # Interpretation
    if has_basal_on['iob_basal'].median() > 0.1:
        print(f"\n  → Basal IOB is significantly positive ({has_basal_on['iob_basal'].median():.3f} U).")
        print(f"    The loop is delivering excess basal → profile basal may be too high.")
        print(f"    This inflates TDD and could make ISF too low.")
    elif has_basal_on['iob_basal'].median() < -0.1:
        print(f"\n  → Basal IOB is negative ({has_basal_on['iob_basal'].median():.3f} U).")
        print(f"    The loop is withholding basal → profile basal may be too low.")
    else:
        print(f"\n  → Basal IOB is near zero ({has_basal_on['iob_basal'].median():.3f} U).")
        print(f"    Profile basal appears well-calibrated overnight.")
else:
    print("  No basal IOB data available.")

# ── 5. TEST 3: TDD sensitivity simulation ────────────────────────────────────
print("\n" + "═" * 72)
print("TEST 3: TDD INFLATION/DEFLATION SIMULATION")
print("═" * 72)
print("  If TDD is inflated by X%, how does formula accuracy change?")

bg_vals  = v['bg'].values
tdd_vals = v['tdd_7day'].values
isf_v1   = v['isf_v1'].values
bg_drop  = v['bg_drop_pred'].values
actual_2h = v['actual_bg_2h'].values

# Test TDD scaling from -50% to +100%
tdd_scales = np.arange(0.5, 2.05, 0.05)
sim_results = []

for scale in tdd_scales:
    tdd_scaled = tdd_vals * scale

    # Current formula with scaled TDD
    isf_ln = (1700.0 / tdd_scaled) * LN_TARGET / v['ln_bg'].values
    pred_ln = bg_vals - bg_drop * (isf_ln / isf_v1)
    err_ln = actual_2h - pred_ln

    # Power-law k=2.0 with scaled TDD
    isf_pw = (1700.0 / tdd_scaled) * (TARGET / bg_vals) ** 2.0
    pred_pw = bg_vals - bg_drop * (isf_pw / isf_v1)
    err_pw = actual_2h - pred_pw

    # Power-law k=2.3 with scaled TDD
    isf_pw23 = (1700.0 / tdd_scaled) * (TARGET / bg_vals) ** 2.3
    pred_pw23 = bg_vals - bg_drop * (isf_pw23 / isf_v1)
    err_pw23 = actual_2h - pred_pw23

    sim_results.append({
        'scale': scale,
        'tdd_pct': (scale - 1) * 100,
        'ln_mae': np.abs(err_ln).mean(),
        'ln_bias': err_ln.mean(),
        'pw20_mae': np.abs(err_pw).mean(),
        'pw20_bias': err_pw.mean(),
        'pw23_mae': np.abs(err_pw23).mean(),
        'pw23_bias': err_pw23.mean(),
    })

sim = pd.DataFrame(sim_results)

print(f"\n  {'TDD Δ':>7s}  │ {'ln MAE':>7s}  {'ln bias':>8s}  │ {'PL2.0 MAE':>9s}  {'PL2.0 bias':>10s}  │ "
      f"{'PL2.3 MAE':>9s}  {'PL2.3 bias':>10s}")
print("  " + "─" * 85)
for _, row in sim.iterrows():
    if row['tdd_pct'] % 10 == 0 or abs(row['tdd_pct']) <= 5:
        marker = " ←" if abs(row['tdd_pct']) < 1 else ""
        print(f"  {row['tdd_pct']:+6.0f}%  │ {row['ln_mae']:7.2f}  {row['ln_bias']:+8.2f}  │ "
              f"{row['pw20_mae']:9.2f}  {row['pw20_bias']:+10.2f}  │ "
              f"{row['pw23_mae']:9.2f}  {row['pw23_bias']:+10.2f}{marker}")

# Find optimal TDD scaling for each formula (where bias = 0)
for label, bias_col in [('Current ln', 'ln_bias'), ('Power-law k=2.0', 'pw20_bias'), ('Power-law k=2.3', 'pw23_bias')]:
    # Linear interpolation to find zero-bias crossing
    biases = sim[bias_col].values
    scales = sim['scale'].values
    for i in range(len(biases) - 1):
        if biases[i] * biases[i+1] <= 0 and biases[i] != biases[i+1]:
            zero_scale = scales[i] + (0 - biases[i]) / (biases[i+1] - biases[i]) * (scales[i+1] - scales[i])
            print(f"\n  {label}: zero bias at TDD × {zero_scale:.3f} ({(zero_scale-1)*100:+.1f}%)")
            break

# ── 6. TEST 4: What if we fit the constant AND k jointly? ────────────────────
print("\n" + "═" * 72)
print("TEST 4: JOINT OPTIMISATION OF CONSTANT AND k")
print("═" * 72)
print("  Fitting: ISF = (C/TDD) × (target/BG)^k")
print("  If C ≠ 1700, the TDD-to-ISF mapping is miscalibrated.")

from scipy.optimize import minimize

def obj_ck(params):
    C, k = params
    isf = (C / tdd_vals) * (TARGET / bg_vals) ** k
    pred = bg_vals - bg_drop * (isf / isf_v1)
    return np.abs(actual_2h - pred).mean()

def obj_ck_bias(params):
    C, k = params
    isf = (C / tdd_vals) * (TARGET / bg_vals) ** k
    pred = bg_vals - bg_drop * (isf / isf_v1)
    return abs((actual_2h - pred).mean())

# MAE optimisation
res_mae = minimize(obj_ck, x0=[1700, 2.0], bounds=[(800, 3000), (0.5, 4.0)], method='L-BFGS-B')
C_mae, k_mae = res_mae.x

# Bias optimisation
res_bias = minimize(obj_ck_bias, x0=[1700, 2.0], bounds=[(800, 3000), (0.5, 4.0)], method='L-BFGS-B')
C_bias, k_bias = res_bias.x

# Compute metrics for optimised
for label, C, k in [('Current (1700, ln)', 1700, None),
                     ('Power-law (1700, k=2.0)', 1700, 2.0),
                     ('Power-law (1700, k=2.3)', 1700, 2.3),
                     (f'Optimal MAE (C={C_mae:.0f}, k={k_mae:.2f})', C_mae, k_mae),
                     (f'Optimal bias (C={C_bias:.0f}, k={k_bias:.2f})', C_bias, k_bias)]:
    if k is None:
        isf = (C / tdd_vals) * LN_TARGET / v['ln_bg'].values
    else:
        isf = (C / tdd_vals) * (TARGET / bg_vals) ** k
    pred = bg_vals - bg_drop * (isf / isf_v1)
    err = actual_2h - pred
    print(f"  {label:45s}  MAE={np.abs(err).mean():.2f}  bias={err.mean():+.2f}  "
          f"RMSE={np.sqrt((err**2).mean()):.2f}")

print(f"\n  Interpretation:")
print(f"    If C_optimal >> 1700 → TDD is inflated (1700/TDD gives ISF too low)")
print(f"    If C_optimal ≈ 1700 → TDD is fine; the BG scaling was the issue")
print(f"    If C_optimal << 1700 → TDD is deflated (uncommon)")

c_ratio = C_mae / 1700
if c_ratio > 1.15:
    tdd_inflation = (c_ratio - 1) * 100
    print(f"\n  → C_optimal = {C_mae:.0f}, which is {tdd_inflation:.0f}% above 1700.")
    print(f"    This is equivalent to TDD being inflated by ~{tdd_inflation:.0f}%.")
    print(f"    The concern about basal overestimation has QUANTITATIVE SUPPORT.")
elif c_ratio < 0.85:
    print(f"\n  → C_optimal = {C_mae:.0f}, which is {(1-c_ratio)*100:.0f}% below 1700.")
    print(f"    TDD may be underestimated for this patient.")
else:
    print(f"\n  → C_optimal = {C_mae:.0f}, which is within 15% of 1700.")
    print(f"    TDD calibration is reasonable; the BG scaling was the dominant issue.")

# ── 7. TEST 5: Decompose — how much error is TDD vs BG scaling? ──────────────
print("\n" + "═" * 72)
print("TEST 5: ERROR DECOMPOSITION — TDD vs BG SCALING")
print("═" * 72)

# Error from TDD alone (at target BG, only TDD matters)
# ISF_implied = ISF_v1 × ratio
# At target BG: ISF should be C/TDD
# Error from TDD = ISF_implied(at target) vs 1700/TDD

# Use samples near target to estimate TDD error
near = v[(v['bg'] >= 90) & (v['bg'] <= 110)].copy()
tdd_error_pct = (near['isf_eff'] / (1700.0 / near['tdd_7day']) - 1).median() * 100

# Use ratio variation across BG bands to estimate BG scaling error
# At target BG, ratio should be 1.0 if TDD is correct
# Deviation from 1.0 across BG levels is the BG scaling error
ratio_at_target = near['ratio'].median()
ratio_at_low = v[v['bg'] < 90]['ratio'].median()
ratio_at_high = v[(v['bg'] >= 120) & (v['bg'] < 150)]['ratio'].median()

print(f"\n  Sensitivity ratio by BG zone:")
print(f"    Low BG (<90):        {ratio_at_low:.3f}")
print(f"    At target (90-110):  {ratio_at_target:.3f}")
print(f"    High BG (120-150):   {ratio_at_high:.3f}")
print(f"\n  ISF_implied / (1700/TDD) at target BG: {1 + tdd_error_pct/100:.3f} ({tdd_error_pct:+.1f}%)")

print(f"\n  Decomposition:")
print(f"    TDD component: ISF should be {tdd_error_pct:+.1f}% {'higher' if tdd_error_pct > 0 else 'lower'} at target BG")
print(f"    BG scaling component: ratio varies from {ratio_at_low:.3f} (low BG) to {ratio_at_high:.3f} (high BG)")
print(f"    BG scaling range: {abs(ratio_at_low - ratio_at_high):.3f}")
print(f"\n    If TDD error >> BG scaling range → TDD inflation is the primary problem")
print(f"    If BG scaling range >> TDD error → BG scaling is the primary problem")

tdd_mag = abs(tdd_error_pct)
bg_mag = abs(ratio_at_low - ratio_at_high) * 100  # convert to percentage terms
print(f"\n    TDD magnitude:       {tdd_mag:.1f}%")
print(f"    BG scaling magnitude: {bg_mag:.1f}%")
if tdd_mag > bg_mag * 1.5:
    print(f"    → TDD miscalibration dominates. The basal inflation concern is VALID.")
elif bg_mag > tdd_mag * 1.5:
    print(f"    → BG scaling dominates. TDD is a secondary factor.")
else:
    print(f"    → Both contribute comparably. Both need addressing.")

# ── 8. Figure ─────────────────────────────────────────────────────────────────
print("\nGenerating figure...")

BG_C = '#0f0f0f'; PANEL = '#1a1a2e'; GRID = '#2a2a4a'; TXT = '#e0e0ff'
C_LN = '#4fc3f7'; C_PW20 = '#66bb6a'; C_PW23 = '#ce93d8'; C_OPT = '#f44336'

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
gs = gridspec.GridSpec(4, 3, figure=fig, hspace=0.52, wspace=0.35)

# P1: TDD inflation — MAE vs TDD scale (full width)
ax1 = fig.add_subplot(gs[0, :])
style(ax1, 'Formula Accuracy vs TDD Inflation/Deflation')
ax1.plot(sim['tdd_pct'], sim['ln_mae'], '-', color=C_LN, lw=2, label='Current ln')
ax1.plot(sim['tdd_pct'], sim['pw20_mae'], '-', color=C_PW20, lw=2, label='Power-law k=2.0')
ax1.plot(sim['tdd_pct'], sim['pw23_mae'], '-', color=C_PW23, lw=2, label='Power-law k=2.3')
ax1.axvline(0, color='white', lw=0.8, ls='--', alpha=0.5, label='Actual TDD')
# Mark minima
for col, c in [('ln_mae', C_LN), ('pw20_mae', C_PW20), ('pw23_mae', C_PW23)]:
    idx = sim[col].idxmin()
    ax1.plot(sim.loc[idx, 'tdd_pct'], sim.loc[idx, col], 'o', color=c, ms=10, zorder=5)
    ax1.annotate(f"min at {sim.loc[idx, 'tdd_pct']:+.0f}%\nMAE={sim.loc[idx, col]:.1f}",
                 xy=(sim.loc[idx, 'tdd_pct'], sim.loc[idx, col]),
                 xytext=(sim.loc[idx, 'tdd_pct'] + 5, sim.loc[idx, col] + 0.5),
                 color=c, fontsize=7)
ax1.set_xlabel('TDD Change (%)'); ax1.set_ylabel('MAE (mg/dL)')
ax1.legend(fontsize=9, labelcolor=TXT, facecolor=PANEL, loc='upper left')

# P2: TDD inflation — bias vs TDD scale
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'Prediction Bias vs TDD Change')
ax2.plot(sim['tdd_pct'], sim['ln_bias'], '-', color=C_LN, lw=2, label='Current ln')
ax2.plot(sim['tdd_pct'], sim['pw20_bias'], '-', color=C_PW20, lw=2, label='PL k=2.0')
ax2.plot(sim['tdd_pct'], sim['pw23_bias'], '-', color=C_PW23, lw=2, label='PL k=2.3')
ax2.axhline(0, color='white', lw=1, ls='--', alpha=0.5)
ax2.axvline(0, color='white', lw=0.8, ls='--', alpha=0.3)
ax2.set_xlabel('TDD Change (%)'); ax2.set_ylabel('Bias (mg/dL)')
ax2.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P3: ISF_implied / (1700/TDD) by BG
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'ISF_implied / (1700/TDD) by BG')
bg_bins = np.arange(72, 160, 8)
v['bg_bin_c'] = pd.cut(v['bg'], bins=bg_bins)
ratio_profile = v.groupby('bg_bin_c', observed=True).apply(
    lambda g: pd.Series({
        'bg_mid': g['bg'].median(),
        'isf_ratio': (g['isf_eff'] / (1700.0 / g['tdd_7day'])).median(),
        'n': len(g),
    })
).reset_index(drop=True)
ax3.plot(ratio_profile['bg_mid'], ratio_profile['isf_ratio'], 'o-', color='#ffb74d', lw=2, ms=5)
ax3.axhline(1.0, color='white', lw=1, ls='--', alpha=0.5, label='1700/TDD is correct')
ax3.axvline(TARGET, color='white', lw=0.8, ls=':', alpha=0.5, label=f'Target {TARGET}')
ax3.set_xlabel('BG (mg/dL)'); ax3.set_ylabel('ISF_implied / (1700/TDD)')
ax3.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)
ax3.set_ylim(0, 3)

# P4: Basal IOB by hour
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, 'Basal IOB by Hour of Night')
if len(has_basal_on) > 0:
    hourly_basal = has_basal_on.groupby('hour').agg(
        iob_basal_med=('iob_basal', 'median'),
        iob_total_med=('iob', 'median'),
        iob_basal_q25=('iob_basal', lambda x: x.quantile(0.25)),
        iob_basal_q75=('iob_basal', lambda x: x.quantile(0.75)),
    )
    ax4.fill_between(hourly_basal.index, hourly_basal['iob_basal_q25'],
                     hourly_basal['iob_basal_q75'], alpha=0.2, color='#f48fb1')
    ax4.plot(hourly_basal.index, hourly_basal['iob_basal_med'], 'o-', color='#f48fb1',
             lw=2, ms=5, label='Basal IOB')
    ax4.plot(hourly_basal.index, hourly_basal['iob_total_med'], 's-', color=C_LN,
             lw=1.5, ms=4, label='Total IOB')
    ax4.axhline(0, color='white', lw=0.8, ls='--', alpha=0.5)
    ax4.set_xlabel('Hour'); ax4.set_ylabel('IOB (U)')
    ax4.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P5: ISF curves — 1700 vs optimal C
ax5 = fig.add_subplot(gs[2, 0])
style(ax5, f'ISF Curves: C=1700 vs C={C_mae:.0f} (optimised)')
bg_range = np.linspace(70, 200, 300)
tdd_med = np.median(tdd_vals)

isf_1700_k2 = (1700 / tdd_med) * (TARGET / bg_range) ** 2.0
isf_opt = (C_mae / tdd_med) * (TARGET / bg_range) ** k_mae
isf_current = (1700 / tdd_med) * LN_TARGET / np.log(bg_range / D + 1)

ax5.plot(bg_range, isf_current, '-', color=C_LN, lw=2, label=f'Current ln (1700)')
ax5.plot(bg_range, isf_1700_k2, '-', color=C_PW20, lw=2, label=f'PL k=2.0 (1700)')
ax5.plot(bg_range, isf_opt, '-', color=C_OPT, lw=2.5, label=f'PL k={k_mae:.1f} (C={C_mae:.0f})')
ax5.axvline(TARGET, color='white', lw=0.8, ls=':', alpha=0.5)
ax5.set_xlabel('BG (mg/dL)'); ax5.set_ylabel('ISF (mg/dL/U)')
ax5.set_ylim(0, 350)
ax5.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P6: Sensitivity — MAE vs constant C for power-law k=2.0
ax6 = fig.add_subplot(gs[2, 1])
style(ax6, 'MAE vs Constant C (power-law k=2.0)')
c_range = np.arange(1000, 2600, 50)
c_maes = []
c_biases = []
for C in c_range:
    isf = (C / tdd_vals) * (TARGET / bg_vals) ** 2.0
    pred = bg_vals - bg_drop * (isf / isf_v1)
    err = actual_2h - pred
    c_maes.append(np.abs(err).mean())
    c_biases.append(err.mean())

ax6.plot(c_range, c_maes, '-', color=C_PW20, lw=2)
ax6.axvline(1700, color='white', lw=1, ls='--', alpha=0.5, label='C=1700')
idx_min = np.argmin(c_maes)
ax6.plot(c_range[idx_min], c_maes[idx_min], 'o', color=C_OPT, ms=10, zorder=5)
ax6.annotate(f'Optimal C={c_range[idx_min]}\nMAE={c_maes[idx_min]:.2f}',
             xy=(c_range[idx_min], c_maes[idx_min]),
             xytext=(c_range[idx_min]+100, c_maes[idx_min]+0.3),
             color=C_OPT, fontsize=8,
             arrowprops=dict(arrowstyle='->', color=C_OPT, lw=1))
ax6.set_xlabel('Constant C'); ax6.set_ylabel('MAE (mg/dL)')
ax6.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P7: Bias vs C
ax7 = fig.add_subplot(gs[2, 2])
style(ax7, 'Prediction Bias vs Constant C (power-law k=2.0)')
ax7.plot(c_range, c_biases, '-', color=C_PW20, lw=2)
ax7.axhline(0, color='white', lw=1, ls='--', alpha=0.5)
ax7.axvline(1700, color='white', lw=1, ls='--', alpha=0.3, label='C=1700')
# Find zero crossing
for i in range(len(c_biases)-1):
    if c_biases[i] * c_biases[i+1] <= 0:
        c_zero = c_range[i] + (0 - c_biases[i]) / (c_biases[i+1] - c_biases[i]) * 50
        ax7.axvline(c_zero, color=C_OPT, lw=1, ls=':', label=f'Zero bias: C={c_zero:.0f}')
        break
ax7.set_xlabel('Constant C'); ax7.set_ylabel('Bias (mg/dL)')
ax7.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P8: Near-target ISF distribution — implied vs 1700/TDD
ax8 = fig.add_subplot(gs[3, 0])
style(ax8, f'ISF at Target BG: Implied vs 1700/TDD (BG {TARGET-10:.0f}-{TARGET+10:.0f})')
ax8.hist(near_target['isf_eff'].clip(0, 400), bins=40, alpha=0.6, color='#ffb74d',
         density=True, label='ISF_implied')
ax8.hist(near_target['isf_1700_tdd'].clip(0, 400), bins=40, alpha=0.6, color=C_LN,
         density=True, label='1700/TDD')
ax8.axvline(near_target['isf_eff'].median(), color='#ffb74d', lw=2, ls='-')
ax8.axvline(near_target['isf_1700_tdd'].median(), color=C_LN, lw=2, ls='-')
ax8.set_xlabel('ISF (mg/dL/U)'); ax8.set_ylabel('Density')
ax8.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P9: Implied constant C distribution
ax9 = fig.add_subplot(gs[3, 1])
style(ax9, 'Implied Constant C at Target BG')
c_vals = near_target['implied_constant'].clip(0, 5000)
ax9.hist(c_vals, bins=40, alpha=0.7, color='#ce93d8', density=True)
ax9.axvline(1700, color='white', lw=1.5, ls='--', label='C=1700')
ax9.axvline(c_median, color=C_OPT, lw=1.5, ls='-', label=f'Median={c_median:.0f}')
ax9.set_xlabel('Constant C'); ax9.set_ylabel('Density')
ax9.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P10: TDD inflation — which formula is more robust?
ax10 = fig.add_subplot(gs[3, 2])
style(ax10, 'TDD Robustness: MAE Degradation at ±20%')
# Show how much MAE worsens at ±20% TDD for each formula
baseline_row = sim[sim['tdd_pct'].abs() < 1].iloc[0]
degradation = []
for pct in [-20, -10, 0, 10, 20]:
    row = sim[(sim['tdd_pct'] - pct).abs() < 3].iloc[0]
    degradation.append({
        'pct': pct,
        'ln_degrad': row['ln_mae'] - baseline_row['ln_mae'],
        'pw20_degrad': row['pw20_mae'] - baseline_row['pw20_mae'],
        'pw23_degrad': row['pw23_mae'] - baseline_row['pw23_mae'],
    })
deg = pd.DataFrame(degradation)
x = np.arange(len(deg))
w = 0.25
ax10.bar(x - w, deg['ln_degrad'], w, color=C_LN, alpha=0.85, label='Current ln')
ax10.bar(x, deg['pw20_degrad'], w, color=C_PW20, alpha=0.85, label='PL k=2.0')
ax10.bar(x + w, deg['pw23_degrad'], w, color=C_PW23, alpha=0.85, label='PL k=2.3')
ax10.set_xticks(x); ax10.set_xticklabels([f'{p:+d}%' for p in deg['pct']], fontsize=8)
ax10.axhline(0, color='white', lw=0.8, ls='--')
ax10.set_xlabel('TDD Error'); ax10.set_ylabel('MAE Degradation (mg/dL)')
ax10.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

fig.suptitle('TDD Inflation Analysis: Is Basal Overestimation a Problem?\n'
             'Testing whether TDD miscalibration or BG scaling drives the formula error',
             color=TXT, fontsize=12, fontweight='bold', y=0.995)

plt.savefig('ns_tdd_inflation_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_tdd_inflation_results.png")

# ── 9. Summary ────────────────────────────────────────────────────────────────
lines = [
    "TDD INFLATION ANALYSIS",
    "=" * 62,
    f"Samples: {len(v):,}  |  Near-target samples: {len(near_target):,}",
    "",
    "QUESTION: Does basal overestimation inflate TDD and make",
    "TDD-based ISF formulas produce values that are too low?",
    "",
    "TEST 1 — Is 1700/TDD correct at target BG?",
    f"  ISF_implied at target: {near_target['isf_eff'].median():.1f} mg/dL/U",
    f"  1700/TDD at target:    {near_target['isf_1700_tdd'].median():.1f} mg/dL/U",
    f"  Ratio: {near_target['isf_ratio_at_target'].median():.3f}",
    f"  → {'YES, 1700/TDD underestimates ISF at target' if near_target['isf_ratio_at_target'].median() > 1.05 else 'NO, 1700/TDD is approximately correct at target' if near_target['isf_ratio_at_target'].median() > 0.95 else '1700/TDD overestimates ISF at target'}",
    "",
    "TEST 1b — Implied constant:",
    f"  C_implied at target: {c_median:.0f} (vs 1700)",
    f"  C_optimal (MAE): {C_mae:.0f}, k={k_mae:.2f}",
    f"  C_optimal (bias): {C_bias:.0f}, k={k_bias:.2f}",
    "",
    "TEST 2 — Basal IOB overnight:",
    f"  Median basal IOB: {has_basal_on['iob_basal'].median():.3f} U" if len(has_basal_on) > 0 else "  No data",
    "",
    "TEST 3 — TDD inflation sensitivity:",
    f"  Current ln:  zero bias at TDD scale ~{sim.loc[sim['ln_bias'].abs().idxmin(), 'scale']:.2f}",
    f"  PL k=2.0:    zero bias at TDD scale ~{sim.loc[sim['pw20_bias'].abs().idxmin(), 'scale']:.2f}",
    "",
    "TEST 5 — Error decomposition:",
    f"  TDD error at target: {tdd_error_pct:+.1f}%",
    f"  BG scaling range: {bg_mag:.1f}%",
    f"  Dominant factor: {'TDD' if tdd_mag > bg_mag * 1.5 else 'BG scaling' if bg_mag > tdd_mag * 1.5 else 'Both comparable'}",
    "",
    "CONCLUSION:",
]

if c_ratio > 1.15:
    lines.append(f"  The implied constant ({c_median:.0f}) is {(c_ratio-1)*100:.0f}% above 1700,")
    lines.append(f"  indicating TDD IS inflated for this patient. The basal overestimation")
    lines.append(f"  concern has quantitative support. However, the BG scaling error")
    lines.append(f"  (ratio varies {abs(ratio_at_low - ratio_at_high):.2f} across BG range)")
    lines.append(f"  is also significant. Both need addressing.")
elif c_ratio > 1.05:
    lines.append(f"  The implied constant ({c_median:.0f}) is modestly above 1700 ({(c_ratio-1)*100:.0f}%),")
    lines.append(f"  suggesting mild TDD inflation. The power-law BG scaling is the")
    lines.append(f"  larger correction needed.")
else:
    lines.append(f"  The implied constant ({c_median:.0f}) is close to 1700,")
    lines.append(f"  suggesting TDD calibration is adequate. The BG scaling term")
    lines.append(f"  is the primary source of error, not TDD inflation.")

with open('ns_tdd_inflation_summary.txt', 'w') as f:
    f.write('\n'.join(lines))

print('\n' + '\n'.join(lines))
print("\nSaved: ns_tdd_inflation_summary.txt")
