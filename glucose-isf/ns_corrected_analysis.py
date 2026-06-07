#!/usr/bin/env python3
"""
Corrected DynamicISF Analysis — Full 10-Month Dataset
=====================================================

Fixes from prior analysis:
  1. TDD computed from actual insulin deliveries (temp basals + boluses),
     NOT back-calculated from variable_sens.
  2. Formula uses the actual Boost TDD blending logic:
       tddWeighted = ((1.4 × tdd_last4h) + (0.6 × tdd_8to4h)) × 3
       TDD = 0.33 × tddWeighted + 0.34 × tdd7D + 0.33 × tdd1D
  3. variable_sens = sensNormalTarget × (1 - (1 - scaler) × velocity)
     where sensNormalTarget = 1800 / (TDD × ln(target/D+1))
     and scaler = ln(target/D+1) / ln(BG/D+1)

Data sources:
  - devicestatus records: BG, variable_sens, IOB, predictions
  - treatment records: temp basals + boluses for actual TDD
  - CGM entries: actual future BG for outcome matching
"""

import json
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import minimize_scalar, minimize
from pathlib import Path

warnings.filterwarnings('ignore')

D      = 82.0     # insulinDivisor: from peak=38 min
TARGET = 99.0     # normalTarget mg/dL (default in Boost)

HOME    = Path.home()
NS_WORK = HOME / 'Nightscout_Work'
OUT_DIR = HOME / 'Downloads'

# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════

def find_json(prefix):
    patterns = [str(NS_WORK / f'{prefix}_*.json'), str(NS_WORK / '*' / f'{prefix}_*.json'), str(HOME / f'{prefix}_*.json')]
    files = set()
    for pat in patterns:
        files.update(glob.glob(pat))
    return sorted(files)


def load_dedup(paths):
    combined = {}
    for p in paths:
        try:
            with open(p) as f:
                first = f.read(1); f.seek(0)
                if first == '[':
                    records = json.load(f)
                elif first == '{':
                    text = f.read()
                    if '}\n{' in text or '}{' in text:
                        text = '[' + text.replace('}\n{', '},{').replace('}{', '},{') + ']'
                    records = json.loads(text) if text.startswith('[') else [json.loads(text)]
                else:
                    continue
        except Exception as e:
            print(f"  SKIP {Path(p).name}: {e}")
            continue
        for r in records:
            if not isinstance(r, dict): continue
            raw_id = r.get('_id')
            if isinstance(raw_id, dict):
                key = raw_id.get('$oid', str(raw_id))
            elif raw_id:
                key = str(raw_id)
            else:
                key = r.get('created_at', '') + '_' + str(r.get('date', ''))
            combined[key] = r
        print(f"  {Path(p).name}: {len(records):,} records  (total unique: {len(combined):,})")
    return list(combined.values())


print("=" * 70)
print("CORRECTED DYNAMICISF ANALYSIS — ACTUAL TDD FROM TREATMENTS")
print("=" * 70)

print("\nLoading devicestatus...")
raw_ds = load_dedup(find_json('devicestatus'))
print(f"  Total unique: {len(raw_ds):,}")

print("\nLoading CGM entries...")
raw_entries = load_dedup(find_json('entries'))
print(f"  Total unique: {len(raw_entries):,}")

print("\nLoading treatments...")
raw_tx = load_dedup(find_json('treatments'))
print(f"  Total unique: {len(raw_tx):,}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. COMPUTE ACTUAL TDD FROM TREATMENTS
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Computing actual TDD from treatments ──")

# Parse temp basals and boluses with timestamps
insulin_events = []

for t in raw_tx:
    try:
        ts = pd.to_datetime(t['created_at'], utc=True)
    except:
        continue
    et = t.get('eventType', '')
    ts_epoch = int(ts.timestamp())

    if 'Temp Basal' in et:
        rate = t.get('rate')
        duration = t.get('duration')
        if rate is not None and duration is not None:
            insulin = float(rate) * float(duration) / 60.0
            insulin_events.append({'ts': ts, 'epoch': ts_epoch, 'insulin': insulin, 'type': 'basal'})

    if t.get('insulin') and float(t.get('insulin', 0)) > 0:
        insulin_events.append({'ts': ts, 'epoch': ts_epoch, 'insulin': float(t['insulin']), 'type': 'bolus'})

df_insulin = pd.DataFrame(insulin_events).sort_values('ts').reset_index(drop=True)
df_insulin['date'] = df_insulin['ts'].dt.date
print(f"  Insulin events: {len(df_insulin):,}")

# Daily actual TDD
daily_tdd = df_insulin.groupby('date')['insulin'].sum().reset_index().rename(columns={'insulin': 'tdd_actual'})
daily_tdd = daily_tdd.sort_values('date')
daily_tdd['tdd_7day'] = daily_tdd['tdd_actual'].rolling(7, min_periods=3).mean()
daily_tdd['tdd_1day'] = daily_tdd['tdd_actual']  # same-day total

print(f"  Daily TDD: median={daily_tdd['tdd_actual'].median():.1f}  mean={daily_tdd['tdd_actual'].mean():.1f}")
print(f"  7D rolling: median={daily_tdd['tdd_7day'].dropna().median():.1f}")

# Build epoch-indexed insulin array for computing tdd_last4h, tdd_8to4h per cycle
insulin_epochs = df_insulin['epoch'].values
insulin_amounts = df_insulin['insulin'].values


def compute_insulin_window(target_epoch, hours_back_start, hours_back_end):
    """Sum insulin delivered between (target - hours_back_start*3600) and (target - hours_back_end*3600)."""
    t_start = target_epoch - int(hours_back_start * 3600)
    t_end = target_epoch - int(hours_back_end * 3600)
    mask = (insulin_epochs >= t_start) & (insulin_epochs < t_end)
    return insulin_amounts[mask].sum()


def compute_boost_tdd(target_epoch, tdd_7d, tdd_1d):
    """Replicate the Boost TDD blending logic."""
    tdd_last4h = compute_insulin_window(target_epoch, 4, 0)
    tdd_8to4h = compute_insulin_window(target_epoch, 8, 4)
    tdd_weighted = ((1.4 * tdd_last4h) + (0.6 * tdd_8to4h)) * 3  # extrapolate to 24h

    if tdd_7d is None or tdd_7d == 0:
        return tdd_weighted, tdd_last4h, tdd_8to4h, tdd_weighted

    if (tdd_weighted < (0.75 * tdd_7d)) and (tdd_1d is not None):
        tdd = ((tdd_weighted + ((tdd_weighted / tdd_7d) * (tdd_7d - tdd_weighted))) * 0.34) + (tdd_1d * 0.33) + (tdd_weighted * 0.33)
    elif tdd_1d is not None:
        tdd = (tdd_weighted * 0.33) + (tdd_7d * 0.34) + (tdd_1d * 0.33)
    else:
        tdd = tdd_weighted

    return tdd, tdd_last4h, tdd_8to4h, tdd_weighted


# ══════════════════════════════════════════════════════════════════════════════
# 3. PARSE DEVICESTATUS → LOOP CYCLES
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Parsing loop cycles ──")
cycles = []
for r in raw_ds:
    try:
        ts = pd.to_datetime(r['created_at'], utc=True)
        sg = r.get('openaps', {}).get('suggested', {})
        if not sg or 'bg' not in sg:
            continue
        ib = r.get('openaps', {}).get('iob', {}) or {}
        pred_iob = sg.get('predBGs', {}).get('IOB') or []

        cycles.append({
            'ts':                ts,
            'bg':                sg.get('bg'),
            'variable_sens':     sg.get('variable_sens'),
            'sensitivity_ratio': sg.get('sensitivityRatio'),
            'cob':               sg.get('COB'),
            'iob':               sg.get('IOB'),
            'target_bg':         sg.get('targetBG'),
            'iob_basal':         ib.get('basaliob'),
            'pred_iob_12':       pred_iob[12] if len(pred_iob) > 12 else None,
            'pred_iob_24':       pred_iob[24] if len(pred_iob) > 24 else None,
        })
    except:
        continue

df = pd.DataFrame(cycles).dropna(subset=['ts', 'bg', 'variable_sens']).sort_values('ts').reset_index(drop=True)
print(f"  Parsed: {len(df):,} cycles  ({df['ts'].min().date()} → {df['ts'].max().date()})")


# ══════════════════════════════════════════════════════════════════════════════
# 4. PARSE CGM + FORWARD LOOKUP
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Parsing CGM entries ──")
cgm_rows = []
for e in raw_entries:
    try:
        sgv = e.get('sgv')
        if not sgv or not (40 <= sgv <= 400): continue
        ca = e.get('created_at')
        ts = pd.to_datetime(ca, utc=True) if ca else pd.to_datetime(int(e['date']), unit='ms', utc=True)
        cgm_rows.append({'ts': ts, 'sgv': float(sgv)})
    except:
        continue

df_cgm = pd.DataFrame(cgm_rows).sort_values('ts').reset_index(drop=True)
cgm_epochs = np.array([int(t.timestamp()) for t in df_cgm['ts']])
cgm_sgv = df_cgm['sgv'].values
print(f"  CGM readings: {len(df_cgm):,}")


def get_cgm_at(target_ts_s, tolerance_s=150):
    idx = np.searchsorted(cgm_epochs, target_ts_s)
    best_val, best_diff = np.nan, np.inf
    for k in (idx - 1, idx):
        if 0 <= k < len(cgm_epochs):
            diff = abs(cgm_epochs[k] - target_ts_s)
            if diff < best_diff:
                best_diff, best_val = diff, cgm_sgv[k]
    return best_val if best_diff <= tolerance_s else np.nan


# Parse bolus history for bolus_age filter
bolus_rows = []
for t in raw_tx:
    if 'Bolus' in t.get('eventType', ''):
        try:
            ins = t.get('insulin')
            ts = pd.to_datetime(t['created_at'], utc=True)
            if ins and float(ins) > 0:
                bolus_rows.append(int(ts.timestamp()))
        except:
            continue
bolus_epochs = np.array(sorted(bolus_rows))


def mins_since_bolus(target_ts_s):
    if len(bolus_epochs) == 0: return 9999.0
    idx = np.searchsorted(bolus_epochs, target_ts_s, side='right') - 1
    return (target_ts_s - bolus_epochs[idx]) / 60.0 if idx >= 0 else 9999.0


print("\n── Forward BG lookups + bolus age ──")
cycle_epochs = np.array([int(t.timestamp()) for t in df['ts']])

actual_2h = np.full(len(df), np.nan)
bolus_age = np.full(len(df), np.nan)

for i, t_s in enumerate(cycle_epochs):
    actual_2h[i] = get_cgm_at(t_s + 7200)
    bolus_age[i] = mins_since_bolus(t_s)
    if (i + 1) % 10000 == 0:
        print(f"  ... {i+1:,}/{len(df):,}")

df['actual_bg_2h'] = actual_2h
df['bolus_age_min'] = bolus_age
df['hour'] = df['ts'].dt.hour
df['date'] = df['ts'].dt.date

print(f"  2h coverage: {(~np.isnan(actual_2h)).sum():,}/{len(df):,}")


# ══════════════════════════════════════════════════════════════════════════════
# 5. ATTACH ACTUAL TDD TO EACH CYCLE
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Computing actual TDD per cycle ──")

# Merge daily TDD
daily_tdd['date'] = pd.to_datetime(daily_tdd['date']).dt.date
df = df.merge(daily_tdd[['date', 'tdd_actual', 'tdd_7day', 'tdd_1day']], on='date', how='left')

# Also compute Boost-style blended TDD for each cycle
# (expensive but necessary for accuracy)
tdd_boost = np.full(len(df), np.nan)
tdd_last4h_arr = np.full(len(df), np.nan)

print("  Computing Boost TDD blend per cycle (slow)...")
for i, (epoch, tdd7, tdd1) in enumerate(zip(cycle_epochs, df['tdd_7day'].values, df['tdd_1day'].values)):
    tdd7_v = tdd7 if not np.isnan(tdd7) else None
    tdd1_v = tdd1 if not np.isnan(tdd1) else None
    tdd, last4h, _, _ = compute_boost_tdd(epoch, tdd7_v, tdd1_v)
    tdd_boost[i] = tdd
    tdd_last4h_arr[i] = last4h
    if (i + 1) % 10000 == 0:
        print(f"  ... {i+1:,}/{len(df):,}")

df['tdd_boost'] = tdd_boost
df['tdd_last4h'] = tdd_last4h_arr

print(f"  TDD boost: median={np.nanmedian(tdd_boost):.1f}")
print(f"  TDD 7-day: median={df['tdd_7day'].median():.1f}")
print(f"  TDD actual daily: median={df['tdd_actual'].median():.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 6. OVERNIGHT FASTING FILTER
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Filtering overnight fasting ──")

mask = (
    (df['hour'] < 8) &
    (df['cob'].fillna(99) == 0) &
    (df['bg'] >= 72) & (df['bg'] <= 200) &
    (df['bolus_age_min'] >= 180) &
    (df['tdd_boost'] > 0) &
    (~np.isnan(df['tdd_7day']))
)
on = df[mask].copy()
print(f"  Overnight fasting: {len(on):,} cycles")

# Further filter for valid 2h outcomes with sufficient predicted drop
on = on.dropna(subset=['pred_iob_24', 'actual_bg_2h']).copy()
on['pred_drop'] = on['bg'] - on['pred_iob_24']
on['actual_drop'] = on['bg'] - on['actual_bg_2h']

strict = on[
    (np.abs(on['pred_drop']) > 3) &
    (on['actual_drop'] * on['pred_drop'] > 0) &  # same direction
    ((on['actual_drop'] / on['pred_drop']) < 5) &
    ((on['actual_drop'] / on['pred_drop']) > 0) &
    (on['bg'] - on['actual_bg_2h'] < 9 + on['pred_drop'])  # no suspected missed carbs
].copy()
print(f"  After strict filtering: {len(strict):,} valid 2h samples")


# ══════════════════════════════════════════════════════════════════════════════
# 7. FORMULA COMPARISON
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Computing formula ISFs ──")

ln_target = np.log(TARGET / D + 1)
bg = strict['bg'].values
ln_bg = np.log(bg / D + 1)
isf_actual = strict['variable_sens'].values

# A: Current formula (what the loop actually used) — ISF = variable_sens
# Already in the data

# B: Reproduced current formula using ACTUAL 7-day TDD
tdd7 = strict['tdd_7day'].values
isf_current_actual_tdd = 1800.0 / (tdd7 * ln_bg)

# C: Reproduced current formula using Boost blended TDD
tdd_b = strict['tdd_boost'].values
isf_current_boost_tdd = 1800.0 / (tdd_b * ln_bg)

# D: Power-law with actual 7-day TDD, k=2.0
k_default = 2.0
isf_pl_7d = (1800.0 / (tdd7 * ln_target)) * (TARGET / bg) ** k_default

# E: Power-law with actual 7-day TDD, fit k
def mae_powerlaw(k, tdd_arr, bg_arr, isf_act, pred_drop, actual_bg_2h):
    isf_f = (1800.0 / (tdd_arr * ln_target)) * (TARGET / bg_arr) ** k
    pred_f = bg_arr - pred_drop * (isf_f / isf_act)
    err = actual_bg_2h - pred_f
    return np.nanmean(np.abs(err))

pred_drop = strict['pred_drop'].values
actual_bg_2h = strict['actual_bg_2h'].values

res_k = minimize_scalar(mae_powerlaw, bounds=(0.5, 3.0), method='bounded',
                         args=(tdd7, bg, isf_actual, pred_drop, actual_bg_2h))
k_opt = res_k.x
print(f"  Optimal k (7-day TDD): {k_opt:.3f}")

isf_pl_kopt = (1800.0 / (tdd7 * ln_target)) * (TARGET / bg) ** k_opt


# Counterfactual predictions for each formula
def compute_errors(isf_formula, label):
    pred_f = bg - pred_drop * (isf_formula / isf_actual)
    err = actual_bg_2h - pred_f
    valid = ~np.isnan(err)
    e = err[valid]
    mae = np.abs(e).mean()
    bias = e.mean()
    rmse = np.sqrt((e**2).mean())
    w18 = (np.abs(e) <= 18).mean() * 100
    return {'label': label, 'mae': mae, 'bias': bias, 'rmse': rmse, 'w18': w18, 'n': len(e), 'errors': e}

formulas = [
    (isf_actual,            'A: Loop actual (variable_sens)'),
    (isf_current_actual_tdd,'B: Current ln, actual 7D-TDD'),
    (isf_current_boost_tdd, 'C: Current ln, Boost blended TDD'),
    (isf_pl_7d,             f'D: Power-law k=2.0, actual 7D-TDD'),
    (isf_pl_kopt,           f'E: Power-law k={k_opt:.2f}, actual 7D-TDD'),
]

print(f"\n{'═'*75}")
print(f"CORRECTED RESULTS — Overnight 00:00–07:00  (n={len(strict):,})")
print(f"  Actual TDD (7-day): median={np.median(tdd7):.1f}  mean={np.mean(tdd7):.1f}")
print(f"  Boost blended TDD:  median={np.median(tdd_b):.1f}  mean={np.mean(tdd_b):.1f}")
print(f"  TDD_implied (old):  {1800.0 / np.median(isf_actual * ln_bg):.1f}  (for comparison)")
print(f"{'═'*75}")
print(f"\n  {'Formula':<45s}  {'MAE':>5s}  {'Bias':>6s}  {'RMSE':>5s}  {'±1mmol':>6s}")
print(f"  {'─'*45}  {'─'*5}  {'─'*6}  {'─'*5}  {'─'*6}")

results = []
for isf_f, label in formulas:
    r = compute_errors(isf_f, label)
    results.append(r)
    print(f"  {label:<45s}  {r['mae']:5.1f}  {r['bias']:+6.1f}  {r['rmse']:5.1f}  {r['w18']:5.1f}%")


# ══════════════════════════════════════════════════════════════════════════════
# 8. BG-BAND ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n── BG-Band Analysis ──")
bands = [('<90', bg < 90), ('90-105', (bg >= 90) & (bg < 105)),
         ('105-120', (bg >= 105) & (bg < 120)), ('120-150', (bg >= 120) & (bg <= 150))]

for label_f, isf_f in [(formulas[1][1], formulas[1][0]), (formulas[3][1], formulas[3][0]), (formulas[4][1], formulas[4][0])]:
    print(f"\n  {label_f}:")
    pred_f = bg - pred_drop * (isf_f / isf_actual)
    err = actual_bg_2h - pred_f
    for bname, bmask in bands:
        e = err[bmask & ~np.isnan(err)]
        if len(e) > 0:
            print(f"    {bname:>8s}  n={len(e):4d}  MAE={np.abs(e).mean():5.1f}  bias={e.mean():+5.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 9. k SENSITIVITY SWEEP
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n── k Sensitivity Sweep ──")
print(f"  {'k':>4s}  {'MAE':>5s}  {'bias':>6s}  {'b_<90':>7s}  {'b_90-105':>8s}  {'b_105-120':>9s}  {'b_120-150':>9s}")

for k in np.arange(0.5, 3.05, 0.1):
    isf_k = (1800.0 / (tdd7 * ln_target)) * (TARGET / bg) ** k
    pred_k = bg - pred_drop * (isf_k / isf_actual)
    err_k = actual_bg_2h - pred_k
    valid = ~np.isnan(err_k)
    mae_k = np.abs(err_k[valid]).mean()
    bias_k = err_k[valid].mean()
    band_biases = []
    for _, bmask in bands:
        e = err_k[bmask & valid]
        band_biases.append(e.mean() if len(e) > 0 else np.nan)
    print(f"  {k:4.1f}  {mae_k:5.2f}  {bias_k:+6.2f}  {band_biases[0]:+7.2f}  {band_biases[1]:+8.2f}  {band_biases[2]:+9.2f}  {band_biases[3]:+9.2f}")


# ══════════════════════════════════════════════════════════════════════════════
# 10. TDD COMPARISON: ACTUAL vs IMPLIED
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n── TDD Comparison ──")
tdd_implied = 1800.0 / (isf_actual * ln_bg)
print(f"  TDD implied (back-calc from variable_sens): median={np.median(tdd_implied):.1f}")
print(f"  TDD actual 7-day:                           median={np.median(tdd7):.1f}")
print(f"  TDD Boost blended:                          median={np.median(tdd_b):.1f}")
print(f"  Ratio actual_7d / implied:                  {np.median(tdd7) / np.median(tdd_implied):.3f}")
print(f"  Ratio boost / implied:                      {np.median(tdd_b) / np.median(tdd_implied):.3f}")


# ══════════════════════════════════════════════════════════════════════════════
# 11. JOINT C + k OPTIMISATION
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n── Joint C + k Optimisation ──")

def mae_ck(params):
    C, k = params
    isf_f = (C / tdd7) * (TARGET / bg) ** k
    pred_f = bg - pred_drop * (isf_f / isf_actual)
    err = actual_bg_2h - pred_f
    return np.nanmean(np.abs(err))

def bias_ck(params):
    C, k = params
    isf_f = (C / tdd7) * (TARGET / bg) ** k
    pred_f = bg - pred_drop * (isf_f / isf_actual)
    err = actual_bg_2h - pred_f
    return abs(np.nanmean(err))

res_mae = minimize(mae_ck, [1800, 2.0], bounds=[(500, 3000), (0.5, 3.0)], method='L-BFGS-B')
res_bias = minimize(bias_ck, [1800, 2.0], bounds=[(500, 3000), (0.5, 3.0)], method='L-BFGS-B')

print(f"  MAE-optimal:  C={res_mae.x[0]:.0f}, k={res_mae.x[1]:.3f}  MAE={res_mae.fun:.2f}")
print(f"  Bias-optimal: C={res_bias.x[0]:.0f}, k={res_bias.x[1]:.3f}  |bias|={res_bias.fun:.3f}")

# Also check: does the old 1700/implied_tdd formula work differently?
print(f"\n  For comparison with 1800 constant:")
isf_1800_7d = 1800.0 / (tdd7 * ln_bg)
r1800 = compute_errors(isf_1800_7d, '1800/(TDD7d × ln(BG/D+1))')
print(f"    1800 current ln:  MAE={r1800['mae']:.2f}  bias={r1800['bias']:+.2f}")

isf_1700_7d = 1700.0 / (tdd7 * ln_bg)
r1700 = compute_errors(isf_1700_7d, '1700/(TDD7d × ln(BG/D+1))')
print(f"    1700 current ln:  MAE={r1700['mae']:.2f}  bias={r1700['bias']:+.2f}")


# ══════════════════════════════════════════════════════════════════════════════
# 12. FIGURE
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating figure...")

BG_C = '#0f0f0f'; PANEL = '#1a1a2e'; GRID = '#2a2a4a'; TXT = '#e0e0ff'

def style(ax, title):
    ax.set_facecolor(PANEL); ax.tick_params(colors=TXT, labelsize=8)
    ax.set_title(title, color=TXT, fontsize=9, fontweight='bold')
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.grid(True, color=GRID, lw=0.5, ls='--', alpha=0.7)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.xaxis.label.set_fontsize(8); ax.yaxis.label.set_fontsize(8)

fig = plt.figure(figsize=(20, 16))
fig.patch.set_facecolor(BG_C)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

colors = ['#4fc3f7', '#f48fb1', '#a5d6a7', '#ffb74d', '#ce93d8']

# P1: ISF curves at median TDD
ax1 = fig.add_subplot(gs[0, :])
tdd_med = np.median(tdd7)
bg_range = np.linspace(70, 200, 300)
ln_bg_r = np.log(bg_range / D + 1)

curve_current = 1800.0 / (tdd_med * ln_bg_r)
curve_pl20 = (1800.0 / (tdd_med * ln_target)) * (TARGET / bg_range) ** 2.0
curve_pl_kopt = (1800.0 / (tdd_med * ln_target)) * (TARGET / bg_range) ** k_opt

style(ax1, f'ISF vs BG at Actual TDD = {tdd_med:.1f} U/day (7-day median)')
ax1.plot(bg_range, curve_current, lw=2.5, color=colors[1], label=f'Current: 1800/(TDD×ln(BG/D+1))')
ax1.plot(bg_range, curve_pl20, lw=2.5, color=colors[3], label=f'Power-law k=2.0')
ax1.plot(bg_range, curve_pl_kopt, lw=2.5, color=colors[4], label=f'Power-law k={k_opt:.2f} (optimal)')
ax1.axvline(TARGET, color='white', lw=0.8, ls=':', alpha=0.6, label=f'Target {TARGET:.0f}')
ax1.set_xlabel('BG (mg/dL)'); ax1.set_ylabel('ISF (mg/dL/U)')
ax1.set_xlim(70, 200)
ax1.legend(fontsize=9, labelcolor=TXT, facecolor=PANEL, loc='upper right')

# P2: TDD comparison — actual vs implied
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'TDD: Actual (treatments) vs Implied (back-calc)')
ax2.hist(tdd_implied.clip(0, 50), bins=40, alpha=0.6, color=colors[0], label=f'Implied (med={np.median(tdd_implied):.1f})', density=True)
ax2.hist(tdd7[~np.isnan(tdd7)].clip(0, 50), bins=40, alpha=0.6, color=colors[1], label=f'Actual 7D (med={np.nanmedian(tdd7):.1f})', density=True)
ax2.set_xlabel('TDD (U/day)'); ax2.set_ylabel('Density')
ax2.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P3: MAE bar chart
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'MAE by Formula (+2h)')
labels_short = ['Loop\nactual', 'ln\n7D-TDD', 'ln\nBoost', f'PL\nk=2.0', f'PL\nk={k_opt:.1f}']
maes = [r['mae'] for r in results]
bars = ax3.bar(range(len(maes)), maes, color=colors, alpha=0.85)
for bar, v in zip(bars, maes):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2, f'{v:.1f}',
             ha='center', fontsize=7, color=TXT)
ax3.set_xticks(range(len(labels_short))); ax3.set_xticklabels(labels_short, fontsize=7)
ax3.set_ylabel('MAE (mg/dL)')

# P4: Bias bar chart
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, 'Prediction Bias (+2h)')
biases = [r['bias'] for r in results]
bars = ax4.bar(range(len(biases)), biases, color=colors, alpha=0.85)
for bar, v in zip(bars, biases):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + (0.3 if v >= 0 else -1.5),
             f'{v:+.1f}', ha='center', fontsize=7, color=TXT)
ax4.axhline(0, color='white', lw=0.8, ls='--')
ax4.set_xticks(range(len(labels_short))); ax4.set_xticklabels(labels_short, fontsize=7)
ax4.set_ylabel('Mean Error (mg/dL)')

# P5: Error distribution
ax5 = fig.add_subplot(gs[2, 0])
style(ax5, 'Prediction Error Distribution (+2h)')
for r, c in zip(results[1:], colors[1:]):
    ax5.hist(r['errors'].clip(-80, 80), bins=30, alpha=0.4, color=c,
             label=r['label'][:20], density=True)
ax5.axvline(0, color='white', lw=0.8, ls='--')
ax5.set_xlabel('Pred Error (mg/dL)'); ax5.set_ylabel('Density')
ax5.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

# P6: k sweep MAE
ax6 = fig.add_subplot(gs[2, 1])
style(ax6, 'k Sensitivity: MAE and Bias')
ks = np.arange(0.5, 3.05, 0.1)
maes_k, biases_k = [], []
for k in ks:
    isf_k = (1800.0 / (tdd7 * ln_target)) * (TARGET / bg) ** k
    pred_k = bg - pred_drop * (isf_k / isf_actual)
    err_k = actual_bg_2h - pred_k
    v = ~np.isnan(err_k)
    maes_k.append(np.abs(err_k[v]).mean())
    biases_k.append(err_k[v].mean())
ax6.plot(ks, maes_k, 'o-', color=colors[3], lw=1.5, ms=3, label='MAE')
ax6b = ax6.twinx()
ax6b.plot(ks, biases_k, 's-', color=colors[4], lw=1.5, ms=3, label='Bias')
ax6b.axhline(0, color='white', lw=0.5, ls='--')
ax6b.tick_params(colors=TXT, labelsize=7)
ax6b.yaxis.label.set_color(TXT)
ax6.set_xlabel('k'); ax6.set_ylabel('MAE (mg/dL)')
ax6b.set_ylabel('Bias (mg/dL)')
ax6.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL, loc='upper left')
ax6b.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL, loc='upper right')

# P7: Daily TDD over time
ax7 = fig.add_subplot(gs[2, 2])
style(ax7, 'Actual TDD Over Time')
dates = pd.to_datetime(daily_tdd['date'])
ax7.scatter(dates, daily_tdd['tdd_actual'], s=8, alpha=0.4, color=colors[0], label='Daily TDD')
ax7.plot(dates, daily_tdd['tdd_7day'], color=colors[3], lw=1.5, label='7-day rolling')
ax7.set_xlabel('Date'); ax7.set_ylabel('TDD (U/day)')
ax7.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)
ax7.tick_params(axis='x', rotation=30)

fig.suptitle(f'Corrected DynamicISF Analysis — Actual TDD from Treatment Records\n'
             f'Jun 2025 – Mar 2026  |  {len(strict):,} valid overnight samples  |  TDD median={tdd_med:.1f} U/day',
             color=TXT, fontsize=12, fontweight='bold', y=0.995)

plt.savefig(OUT_DIR / 'ns_corrected_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_corrected_results.png")


# ══════════════════════════════════════════════════════════════════════════════
# 13. SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
summary = f"""CORRECTED DYNAMICISF ANALYSIS — SUMMARY
{'='*60}
Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}
Dataset: Jun 2025 – Mar 2026 (10 months)
Valid overnight +2h samples: {len(strict):,}

KEY CORRECTION:
  Previous analysis used TDD back-calculated from variable_sens: ~14.9 U/day
  Actual TDD from treatment records (bolus + temp basal): {np.median(tdd7):.1f} U/day
  Ratio: actual/implied = {np.median(tdd7)/np.median(tdd_implied):.2f}x

  The previous analysis was using a synthetic TDD that was {100*(np.median(tdd7)/np.median(tdd_implied) - 1):+.0f}%
  different from the real insulin delivered.

BOOST FORMULA (actual code):
  sensNormalTarget = 1800 / (TDD * ln(target/D + 1))
  scaler = ln(target/D + 1) / ln(BG/D + 1)
  variable_sens = sensNormalTarget * (1 - (1 - scaler) * velocity)

  TDD = weighted blend of 7-day avg, 1-day total, and recent 4h/8h windows
  D = {D:.0f} (insulinDivisor, from peak time 38 min)
  target = {TARGET:.0f} mg/dL (normalTarget)

RESULTS (+2h overnight):
"""
for r in results:
    summary += f"  {r['label']:<45s}  MAE={r['mae']:.1f}  bias={r['bias']:+.1f}  ±1mmol={r['w18']:.1f}%\n"

summary += f"""
OPTIMAL k (power-law, actual 7D-TDD): {k_opt:.3f}

JOINT C + k OPTIMISATION:
  MAE-optimal:  C={res_mae.x[0]:.0f}, k={res_mae.x[1]:.3f}
  Bias-optimal: C={res_bias.x[0]:.0f}, k={res_bias.x[1]:.3f}
"""

with open(OUT_DIR / 'ns_corrected_summary.txt', 'w') as f:
    f.write(summary)
print(summary)
print("Saved: ns_corrected_summary.txt")
