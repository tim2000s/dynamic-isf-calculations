#!/usr/bin/env python3
"""
TDD Variant Analysis for Power-Law DynamicISF
==============================================

Tests multiple TDD computation methods with the power-law formula to determine
which TDD input produces the best overnight predictions and the most stable
(C, k) parameters.

TDD variants tested:
  1. Actual 7-day rolling average (simple, stable)
  2. Boost blended TDD (complex, overnight-deflated)
  3. Simple 50/50 blend: 0.5 × 7D + 0.5 × 1D
  4. Simple 70/30 blend: 0.7 × 7D + 0.3 × 1D
  5. 7-day TDD with recent-delivery discount: tdd7D × (1 - overnight_discount)
  6. Implied TDD (back-calculated from variable_sens) — as baseline/reference

For each variant: sweep k, find optimal k, do joint C+k optimisation,
and compute per-BG-band performance.
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

D      = 82.0
TARGET = 99.0
LN_TARGET = np.log(TARGET / D + 1)

HOME    = Path.home()
NS_WORK = HOME / 'Nightscout_Work'
OUT_DIR = HOME / 'Downloads'


# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA (same as ns_corrected_analysis.py)
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
print("TDD VARIANT ANALYSIS FOR POWER-LAW DYNAMICISF")
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

daily_tdd = df_insulin.groupby('date')['insulin'].sum().reset_index().rename(columns={'insulin': 'tdd_actual'})
daily_tdd = daily_tdd.sort_values('date')
daily_tdd['tdd_7day'] = daily_tdd['tdd_actual'].rolling(7, min_periods=3).mean()
daily_tdd['tdd_1day'] = daily_tdd['tdd_actual']

# Also compute per-day basal vs bolus
daily_basal = df_insulin[df_insulin['type'] == 'basal'].groupby('date')['insulin'].sum().reset_index().rename(columns={'insulin': 'basal_insulin'})
daily_bolus = df_insulin[df_insulin['type'] == 'bolus'].groupby('date')['insulin'].sum().reset_index().rename(columns={'insulin': 'bolus_insulin'})
daily_tdd['date_dt'] = pd.to_datetime(daily_tdd['date'])
daily_tdd = daily_tdd.merge(daily_basal.rename(columns={'date': 'date'}), on='date', how='left')
daily_tdd = daily_tdd.merge(daily_bolus.rename(columns={'date': 'date'}), on='date', how='left')
daily_tdd['basal_insulin'] = daily_tdd['basal_insulin'].fillna(0)
daily_tdd['bolus_insulin'] = daily_tdd['bolus_insulin'].fillna(0)

print(f"  Daily TDD: median={daily_tdd['tdd_actual'].median():.1f}  mean={daily_tdd['tdd_actual'].mean():.1f}")
print(f"  7D rolling: median={daily_tdd['tdd_7day'].dropna().median():.1f}")
print(f"  Basal: median={daily_tdd['basal_insulin'].median():.1f}  Bolus: median={daily_tdd['bolus_insulin'].median():.1f}")

insulin_epochs = df_insulin['epoch'].values
insulin_amounts = df_insulin['insulin'].values


def compute_insulin_window(target_epoch, hours_back_start, hours_back_end):
    t_start = target_epoch - int(hours_back_start * 3600)
    t_end = target_epoch - int(hours_back_end * 3600)
    mask = (insulin_epochs >= t_start) & (insulin_epochs < t_end)
    return insulin_amounts[mask].sum()


def compute_boost_tdd(target_epoch, tdd_7d, tdd_1d):
    tdd_last4h = compute_insulin_window(target_epoch, 4, 0)
    tdd_8to4h = compute_insulin_window(target_epoch, 8, 4)
    tdd_weighted = ((1.4 * tdd_last4h) + (0.6 * tdd_8to4h)) * 3

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
# 4. CGM + FORWARD LOOKUP
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
# 5. ATTACH TDD VARIANTS TO EACH CYCLE
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Computing TDD variants per cycle ──")

daily_tdd['date'] = pd.to_datetime(daily_tdd['date']).dt.date
df = df.merge(daily_tdd[['date', 'tdd_actual', 'tdd_7day', 'tdd_1day', 'basal_insulin', 'bolus_insulin']], on='date', how='left')

# Boost blended TDD per cycle
tdd_boost = np.full(len(df), np.nan)
tdd_last4h_arr = np.full(len(df), np.nan)
tdd_last24h_arr = np.full(len(df), np.nan)

print("  Computing Boost TDD blend + last24h per cycle...")
for i, (epoch, tdd7, tdd1) in enumerate(zip(cycle_epochs, df['tdd_7day'].values, df['tdd_1day'].values)):
    tdd7_v = tdd7 if not np.isnan(tdd7) else None
    tdd1_v = tdd1 if not np.isnan(tdd1) else None
    tdd, last4h, _, _ = compute_boost_tdd(epoch, tdd7_v, tdd1_v)
    tdd_boost[i] = tdd
    tdd_last4h_arr[i] = last4h
    tdd_last24h_arr[i] = compute_insulin_window(epoch, 24, 0)
    if (i + 1) % 10000 == 0:
        print(f"  ... {i+1:,}/{len(df):,}")

df['tdd_boost'] = tdd_boost
df['tdd_last4h'] = tdd_last4h_arr
df['tdd_last24h'] = tdd_last24h_arr

# Compute TDD variant arrays
df['tdd_50_50'] = 0.5 * df['tdd_7day'] + 0.5 * df['tdd_1day']
df['tdd_70_30'] = 0.7 * df['tdd_7day'] + 0.3 * df['tdd_1day']

# Simple overnight discount: reduce 7D-TDD by a fixed fraction overnight
# The rationale: overnight, basal-only delivery means actual recent insulin
# is lower than the daily average. A 15-25% discount approximates this.
df['tdd_7d_disc15'] = df['tdd_7day'] * 0.85
df['tdd_7d_disc25'] = df['tdd_7day'] * 0.75

# TDD from last 24h of actual delivery
df['tdd_24h_actual'] = tdd_last24h_arr

# Implied TDD (back-calculated — for reference only)
ln_bg_all = np.log(df['bg'].values / D + 1)
df['tdd_implied'] = 1800.0 / (df['variable_sens'].values * ln_bg_all)

print(f"\n  TDD variant medians (all cycles):")
for col in ['tdd_7day', 'tdd_1day', 'tdd_boost', 'tdd_50_50', 'tdd_70_30',
            'tdd_7d_disc15', 'tdd_7d_disc25', 'tdd_24h_actual', 'tdd_implied']:
    v = df[col].dropna()
    print(f"    {col:<18s}: median={v.median():.1f}  mean={v.mean():.1f}")


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

on = on.dropna(subset=['pred_iob_24', 'actual_bg_2h']).copy()
on['pred_drop'] = on['bg'] - on['pred_iob_24']
on['actual_drop'] = on['bg'] - on['actual_bg_2h']

strict = on[
    (np.abs(on['pred_drop']) > 3) &
    (on['actual_drop'] * on['pred_drop'] > 0) &
    ((on['actual_drop'] / on['pred_drop']) < 5) &
    ((on['actual_drop'] / on['pred_drop']) > 0) &
    (on['bg'] - on['actual_bg_2h'] < 9 + on['pred_drop'])
].copy()
print(f"  After strict filtering: {len(strict):,} valid 2h samples")

# Extract core arrays for speed
bg = strict['bg'].values
ln_bg = np.log(bg / D + 1)
isf_actual = strict['variable_sens'].values
pred_drop = strict['pred_drop'].values
actual_bg_2h = strict['actual_bg_2h'].values


# ══════════════════════════════════════════════════════════════════════════════
# 7. EVALUATE EACH TDD VARIANT
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("TDD VARIANT COMPARISON — POWER-LAW ISF")
print("=" * 80)

bands = [('<90', bg < 90), ('90-105', (bg >= 90) & (bg < 105)),
         ('105-120', (bg >= 105) & (bg < 120)), ('120-150', (bg >= 120) & (bg <= 150))]


def evaluate_tdd_variant(tdd_arr, name, C_base=1800.0):
    """For a given TDD array, find optimal k, joint (C,k), and per-band metrics."""
    valid_tdd = ~np.isnan(tdd_arr) & (tdd_arr > 0)
    tdd = np.where(valid_tdd, tdd_arr, np.nan)
    n_valid = valid_tdd.sum()

    if n_valid < 100:
        return None

    # ISF for current ln formula with this TDD
    isf_ln = 1800.0 / (tdd * ln_bg)
    pred_ln = bg - pred_drop * (isf_ln / isf_actual)
    err_ln = actual_bg_2h - pred_ln
    v_ln = ~np.isnan(err_ln)

    mae_ln = np.abs(err_ln[v_ln]).mean()
    bias_ln = err_ln[v_ln].mean()

    # Sweep k for power-law with this TDD
    def mae_pl(k):
        isf_f = (C_base / tdd) * (TARGET / bg) ** k
        pred_f = bg - pred_drop * (isf_f / isf_actual)
        err = actual_bg_2h - pred_f
        v = ~np.isnan(err)
        return np.abs(err[v]).mean()

    def bias_pl(k):
        isf_f = (C_base / tdd) * (TARGET / bg) ** k
        pred_f = bg - pred_drop * (isf_f / isf_actual)
        err = actual_bg_2h - pred_f
        v = ~np.isnan(err)
        return abs(err[v].mean())

    # Optimal k for MAE
    res_k_mae = minimize_scalar(mae_pl, bounds=(0.5, 4.0), method='bounded')
    k_opt_mae = res_k_mae.x
    mae_opt = res_k_mae.fun

    # Optimal k for bias
    res_k_bias = minimize_scalar(bias_pl, bounds=(0.5, 4.0), method='bounded')
    k_opt_bias = res_k_bias.x

    # Power-law at optimal k: full metrics
    isf_opt = (C_base / tdd) * (TARGET / bg) ** k_opt_mae
    pred_opt = bg - pred_drop * (isf_opt / isf_actual)
    err_opt = actual_bg_2h - pred_opt
    v_opt = ~np.isnan(err_opt)
    bias_at_kopt = err_opt[v_opt].mean()
    w18_opt = (np.abs(err_opt[v_opt]) <= 18).mean() * 100

    # Joint C + k optimisation
    def mae_ck(params):
        C, k = params
        isf_f = (C / tdd) * (TARGET / bg) ** k
        pred_f = bg - pred_drop * (isf_f / isf_actual)
        err = actual_bg_2h - pred_f
        return np.nanmean(np.abs(err))

    def bias_ck(params):
        C, k = params
        isf_f = (C / tdd) * (TARGET / bg) ** k
        pred_f = bg - pred_drop * (isf_f / isf_actual)
        err = actual_bg_2h - pred_f
        return abs(np.nanmean(err))

    res_ck_mae = minimize(mae_ck, [1800, 2.0], bounds=[(500, 5000), (0.5, 4.0)], method='L-BFGS-B')
    res_ck_bias = minimize(bias_ck, [1800, 2.0], bounds=[(500, 5000), (0.5, 4.0)], method='L-BFGS-B')

    # Per-BG-band at optimal k
    band_results = []
    for bname, bmask in bands:
        e = err_opt[bmask & v_opt]
        if len(e) > 0:
            band_results.append({'band': bname, 'n': len(e), 'mae': np.abs(e).mean(), 'bias': e.mean()})

    # k sweep for detailed view
    k_sweep = []
    for k in np.arange(0.5, 4.05, 0.25):
        isf_k = (C_base / tdd) * (TARGET / bg) ** k
        pred_k = bg - pred_drop * (isf_k / isf_actual)
        err_k = actual_bg_2h - pred_k
        v_k = ~np.isnan(err_k)
        k_sweep.append({
            'k': k,
            'mae': np.abs(err_k[v_k]).mean(),
            'bias': err_k[v_k].mean(),
            'w18': (np.abs(err_k[v_k]) <= 18).mean() * 100
        })

    return {
        'name': name,
        'tdd_median': np.nanmedian(tdd),
        'tdd_mean': np.nanmean(tdd),
        'n_valid': int(v_ln.sum()),
        'mae_ln': mae_ln,
        'bias_ln': bias_ln,
        'k_opt_mae': k_opt_mae,
        'k_opt_bias': k_opt_bias,
        'mae_at_kopt': mae_opt,
        'bias_at_kopt': bias_at_kopt,
        'w18_at_kopt': w18_opt,
        'joint_mae_C': res_ck_mae.x[0],
        'joint_mae_k': res_ck_mae.x[1],
        'joint_mae_val': res_ck_mae.fun,
        'joint_bias_C': res_ck_bias.x[0],
        'joint_bias_k': res_ck_bias.x[1],
        'joint_bias_val': res_ck_bias.fun,
        'band_results': band_results,
        'k_sweep': k_sweep,
    }


# Define TDD variants to test
tdd_variants = [
    ('tdd_7day',       'Actual 7-day rolling'),
    ('tdd_1day',       'Actual 1-day (same day)'),
    ('tdd_boost',      'Boost blended'),
    ('tdd_50_50',      '50/50 blend (7D+1D)'),
    ('tdd_70_30',      '70/30 blend (7D+1D)'),
    ('tdd_7d_disc15',  '7-day × 0.85 (15% discount)'),
    ('tdd_7d_disc25',  '7-day × 0.75 (25% discount)'),
    ('tdd_24h_actual', 'Last 24h actual delivery'),
    ('tdd_implied',    'Implied (back-calc, reference)'),
]

all_results = []
for col, name in tdd_variants:
    print(f"\n── {name} (median={strict[col].median():.1f}) ──")
    r = evaluate_tdd_variant(strict[col].values, name)
    if r:
        all_results.append(r)
        print(f"  ln formula:       MAE={r['mae_ln']:.2f}  bias={r['bias_ln']:+.2f}")
        print(f"  PL k_opt(MAE):    k={r['k_opt_mae']:.3f}  MAE={r['mae_at_kopt']:.2f}  bias={r['bias_at_kopt']:+.2f}  ±1mmol={r['w18_at_kopt']:.1f}%")
        print(f"  PL k_opt(bias):   k={r['k_opt_bias']:.3f}")
        print(f"  Joint MAE:        C={r['joint_mae_C']:.0f}  k={r['joint_mae_k']:.3f}  MAE={r['joint_mae_val']:.2f}")
        print(f"  Joint bias:       C={r['joint_bias_C']:.0f}  k={r['joint_bias_k']:.3f}  |bias|={r['joint_bias_val']:.3f}")


# ══════════════════════════════════════════════════════════════════════════════
# 8. ADDITIONAL: Power-law with Boost blended TDD — sweep k and joint C+k
# ══════════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 80)
print("DETAILED COMPARISON: POWER-LAW WITH DIFFERENT TDD INPUTS")
print("=" * 80)

# Also test: what if we use Boost blended but allow C to float?
# This tells us: does the blended TDD with a natural C produce good results?


# ══════════════════════════════════════════════════════════════════════════════
# 9. SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 80)
print("SUMMARY TABLE")
print("=" * 80)
print(f"\n  {'TDD Variant':<30s}  {'med':>5s}  {'ln MAE':>6s}  {'ln bias':>7s}  {'PL k*':>5s}  {'PL MAE':>6s}  {'PL bias':>7s}  {'±1mm':>5s}  {'C*':>5s}  {'k*(C)':>5s}  {'MAE*':>5s}")
print(f"  {'─'*30}  {'─'*5}  {'─'*6}  {'─'*7}  {'─'*5}  {'─'*6}  {'─'*7}  {'─'*5}  {'─'*5}  {'─'*5}  {'─'*5}")

for r in all_results:
    print(f"  {r['name']:<30s}  {r['tdd_median']:5.1f}  {r['mae_ln']:6.2f}  {r['bias_ln']:+7.2f}  {r['k_opt_mae']:5.2f}  {r['mae_at_kopt']:6.2f}  {r['bias_at_kopt']:+7.2f}  {r['w18_at_kopt']:5.1f}  {r['joint_mae_C']:5.0f}  {r['joint_mae_k']:5.2f}  {r['joint_mae_val']:5.2f}")


# ══════════════════════════════════════════════════════════════════════════════
# 10. PER-BG-BAND TABLE FOR TOP 3
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n── Per-BG-Band at optimal k (top variants) ──")
for r in all_results:
    print(f"\n  {r['name']} (k={r['k_opt_mae']:.2f}):")
    for b in r['band_results']:
        print(f"    {b['band']:>8s}  n={b['n']:4d}  MAE={b['mae']:5.1f}  bias={b['bias']:+5.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 11. STABILITY: how much does k vary across TDD inputs?
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n── Parameter stability ──")
# Exclude implied (reference) from stability analysis
real_results = [r for r in all_results if 'Implied' not in r['name']]
ks = [r['k_opt_mae'] for r in real_results]
Cs = [r['joint_mae_C'] for r in real_results]
maes = [r['mae_at_kopt'] for r in real_results]
joint_maes = [r['joint_mae_val'] for r in real_results]

print(f"  k_opt range: {min(ks):.2f} – {max(ks):.2f}  (spread: {max(ks)-min(ks):.2f})")
print(f"  C_opt range: {min(Cs):.0f} – {max(Cs):.0f}  (spread: {max(Cs)-min(Cs):.0f})")
print(f"  MAE range (at PL k_opt): {min(maes):.2f} – {max(maes):.2f}")
print(f"  MAE range (at joint C,k): {min(joint_maes):.2f} – {max(joint_maes):.2f}")

# Best joint MAE across all variants
best = min(all_results, key=lambda r: r['joint_mae_val'])
print(f"\n  Best overall (joint C,k optimised):")
print(f"    {best['name']}: C={best['joint_mae_C']:.0f}, k={best['joint_mae_k']:.2f}, MAE={best['joint_mae_val']:.2f}")

# Best MAE with C=1800 fixed
best_fixed = min(all_results, key=lambda r: r['mae_at_kopt'])
print(f"\n  Best with C=1800 fixed:")
print(f"    {best_fixed['name']}: k={best_fixed['k_opt_mae']:.2f}, MAE={best_fixed['mae_at_kopt']:.2f}")


# ══════════════════════════════════════════════════════════════════════════════
# 12. FIGURE
# ══════════════════════════════════════════════════════════════════════════════
print("\n\nGenerating figure...")

BG_C = '#0f0f0f'; PANEL = '#1a1a2e'; GRID = '#2a2a4a'; TXT = '#e0e0ff'

def style(ax, title):
    ax.set_facecolor(PANEL); ax.tick_params(colors=TXT, labelsize=8)
    ax.set_title(title, color=TXT, fontsize=9, fontweight='bold')
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.grid(True, color=GRID, lw=0.5, ls='--', alpha=0.7)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.xaxis.label.set_fontsize(8); ax.yaxis.label.set_fontsize(8)

fig = plt.figure(figsize=(22, 18))
fig.patch.set_facecolor(BG_C)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

colors = ['#4fc3f7', '#f48fb1', '#a5d6a7', '#ffb74d', '#ce93d8',
          '#80deea', '#ef9a9a', '#c5e1a5', '#ffe082']

# P1: Summary bar chart — MAE for each TDD variant (ln vs PL optimal)
ax1 = fig.add_subplot(gs[0, :])
style(ax1, 'MAE by TDD Variant: Current ln (C=1800) vs Power-Law at Optimal k')
x = np.arange(len(all_results))
w = 0.35
bars1 = ax1.bar(x - w/2, [r['mae_ln'] for r in all_results], w, color='#ef5350', alpha=0.8, label='Current ln')
bars2 = ax1.bar(x + w/2, [r['mae_at_kopt'] for r in all_results], w, color='#66bb6a', alpha=0.8, label='Power-law (k opt)')
for bar, v in zip(bars1, [r['mae_ln'] for r in all_results]):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, f'{v:.1f}',
             ha='center', fontsize=6, color=TXT)
for bar, v, r in zip(bars2, [r['mae_at_kopt'] for r in all_results], all_results):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, f'{v:.1f}\nk={r["k_opt_mae"]:.1f}',
             ha='center', fontsize=6, color=TXT)
ax1.set_xticks(x)
ax1.set_xticklabels([r['name'] for r in all_results], fontsize=7, rotation=25, ha='right')
ax1.set_ylabel('MAE (mg/dL)')
ax1.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P2: Bias comparison
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'Bias by TDD Variant (Power-Law at Optimal k)')
bars = ax2.bar(range(len(all_results)), [r['bias_at_kopt'] for r in all_results],
               color=colors[:len(all_results)], alpha=0.85)
for bar, v in zip(bars, [r['bias_at_kopt'] for r in all_results]):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + (0.2 if v >= 0 else -1),
             f'{v:+.1f}', ha='center', fontsize=6, color=TXT)
ax2.axhline(0, color='white', lw=0.8, ls='--')
ax2.set_xticks(range(len(all_results)))
ax2.set_xticklabels([r['name'][:15] for r in all_results], fontsize=6, rotation=45, ha='right')
ax2.set_ylabel('Bias (mg/dL)')

# P3: Optimal k by variant
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'Optimal k (MAE) by TDD Variant')
bars = ax3.bar(range(len(all_results)), [r['k_opt_mae'] for r in all_results],
               color=colors[:len(all_results)], alpha=0.85)
for bar, v in zip(bars, [r['k_opt_mae'] for r in all_results]):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
             f'{v:.2f}', ha='center', fontsize=6, color=TXT)
ax3.set_xticks(range(len(all_results)))
ax3.set_xticklabels([r['name'][:15] for r in all_results], fontsize=6, rotation=45, ha='right')
ax3.set_ylabel('Optimal k')

# P4: Joint optimal C by variant
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, 'Joint Optimal C (MAE) by TDD Variant')
bars = ax4.bar(range(len(all_results)), [r['joint_mae_C'] for r in all_results],
               color=colors[:len(all_results)], alpha=0.85)
for bar, v in zip(bars, [r['joint_mae_C'] for r in all_results]):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
             f'{v:.0f}', ha='center', fontsize=6, color=TXT)
ax4.set_xticks(range(len(all_results)))
ax4.set_xticklabels([r['name'][:15] for r in all_results], fontsize=6, rotation=45, ha='right')
ax4.set_ylabel('Optimal C')

# P5-P7: k sweep curves for the most interesting TDD variants
interesting = ['Actual 7-day rolling', 'Boost blended', '7-day × 0.75 (25% discount)',
               'Implied (back-calc, reference)']

ax5 = fig.add_subplot(gs[2, 0])
style(ax5, 'k Sweep: MAE by TDD Variant')
for i, r in enumerate(all_results):
    if r['name'] in interesting:
        ks = [s['k'] for s in r['k_sweep']]
        maes = [s['mae'] for s in r['k_sweep']]
        ax5.plot(ks, maes, 'o-', ms=3, lw=1.5, label=f"{r['name'][:20]} (med={r['tdd_median']:.0f})")
ax5.set_xlabel('k'); ax5.set_ylabel('MAE (mg/dL)')
ax5.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

ax6 = fig.add_subplot(gs[2, 1])
style(ax6, 'k Sweep: Bias by TDD Variant')
for i, r in enumerate(all_results):
    if r['name'] in interesting:
        ks = [s['k'] for s in r['k_sweep']]
        biases = [s['bias'] for s in r['k_sweep']]
        ax6.plot(ks, biases, 'o-', ms=3, lw=1.5, label=f"{r['name'][:20]} (med={r['tdd_median']:.0f})")
ax6.axhline(0, color='white', lw=0.8, ls='--')
ax6.set_xlabel('k'); ax6.set_ylabel('Bias (mg/dL)')
ax6.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

# P7: TDD distributions overnight
ax7 = fig.add_subplot(gs[2, 2])
style(ax7, 'Overnight TDD Distributions')
for col, name, c in [('tdd_7day', '7-day', colors[0]),
                      ('tdd_boost', 'Boost blended', colors[1]),
                      ('tdd_7d_disc25', '7D×0.75', colors[2]),
                      ('tdd_implied', 'Implied', colors[4])]:
    v = strict[col].dropna().values
    ax7.hist(v.clip(0, 50), bins=40, alpha=0.4, color=c,
             label=f'{name} (med={np.median(v):.1f})', density=True)
ax7.set_xlabel('TDD (U/day)'); ax7.set_ylabel('Density')
ax7.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

fig.suptitle(f'TDD Variant Analysis for Power-Law DynamicISF\n'
             f'Jun 2025 – Mar 2026  |  {len(strict):,} valid overnight samples',
             color=TXT, fontsize=12, fontweight='bold', y=0.995)

plt.savefig(OUT_DIR / 'ns_tdd_variants_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_tdd_variants_results.png")


# ══════════════════════════════════════════════════════════════════════════════
# 13. WRITE SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
summary_lines = [
    "TDD VARIANT ANALYSIS — SUMMARY",
    "=" * 60,
    f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}",
    f"Valid overnight samples: {len(strict):,}",
    "",
    "VARIANT COMPARISON (Power-law at optimal k, C=1800 fixed):",
    f"  {'Variant':<30s}  {'TDD med':>7s}  {'k_opt':>5s}  {'MAE':>5s}  {'bias':>7s}  {'±1mm':>5s}",
    f"  {'─'*30}  {'─'*7}  {'─'*5}  {'─'*5}  {'─'*7}  {'─'*5}",
]

for r in all_results:
    summary_lines.append(
        f"  {r['name']:<30s}  {r['tdd_median']:7.1f}  {r['k_opt_mae']:5.2f}  {r['mae_at_kopt']:5.2f}  {r['bias_at_kopt']:+7.2f}  {r['w18_at_kopt']:5.1f}"
    )

summary_lines += [
    "",
    "JOINT C + k OPTIMISATION (both free):",
    f"  {'Variant':<30s}  {'C':>5s}  {'k':>5s}  {'MAE':>5s}",
    f"  {'─'*30}  {'─'*5}  {'─'*5}  {'─'*5}",
]
for r in all_results:
    summary_lines.append(
        f"  {r['name']:<30s}  {r['joint_mae_C']:5.0f}  {r['joint_mae_k']:5.2f}  {r['joint_mae_val']:5.2f}"
    )

summary_lines += [
    "",
    "KEY FINDING:",
    f"  Best MAE (joint C,k): {best['name']} — C={best['joint_mae_C']:.0f}, k={best['joint_mae_k']:.2f}, MAE={best['joint_mae_val']:.2f}",
    f"  Best MAE (C=1800):    {best_fixed['name']} — k={best_fixed['k_opt_mae']:.2f}, MAE={best_fixed['mae_at_kopt']:.2f}",
    f"  k stability (real TDD variants): {min(ks):.2f} – {max(ks):.2f}",
    f"  C stability (joint optimised):   {min(Cs):.0f} – {max(Cs):.0f}",
]

summary_text = "\n".join(summary_lines)
with open(OUT_DIR / 'ns_tdd_variants_summary.txt', 'w') as f:
    f.write(summary_text)

print("\n" + summary_text)
print("\nSaved: ns_tdd_variants_summary.txt")
