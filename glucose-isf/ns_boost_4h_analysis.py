#!/usr/bin/env python3
"""
Boost Single-Patient + Synthetic Multi-Patient — End-of-IOB Analysis
=====================================================================

Extends the 2-hour Boost/AAPS analysis to the end-of-IOB prediction horizon,
matching the multi-site Trio analysis methodology.

Key changes from 2h version:
  - Uses pred_iob[-1] (last element) instead of pred_iob[24]
  - Looks up actual SGV at the dynamic prediction horizon
  - Includes falling/rising directional analysis
  - Runs synthetic 7-patient validation at end-of-IOB
"""

import json
import glob
import warnings
import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings('ignore')

D      = 82.0
TARGET = 99.0

HOME    = Path.home()
NS_WORK = HOME / 'Nightscout_Work'
OUT_DIR = HOME / 'Downloads' / '4 Hour analysis'
CACHE   = OUT_DIR / 'boost_4h_cache.pkl'


# ── ISF formulas ──
def isf_quartic(bg):
    G = np.asarray(bg, dtype=float)
    return 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4

def isf_full_diabeloop(bg):
    G = np.asarray(bg, dtype=float)
    quartic = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    quadratic = 98.03 - 1.077*G + 0.008868*G**2
    return np.where(G > 100, quartic, quadratic)

def isf_hybrid(bg):
    G = np.asarray(bg, dtype=float)
    poly = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    power = 75.8 * (105.0 / G) ** 3.5
    return np.where(G >= 105, poly, power)

QUARTIC_AT_99 = float(isf_quartic(99.0))    # ~81.6
FULL_DB_AT_99 = float(isf_full_diabeloop(99.0))  # ~81.6
HYBRID_AT_105 = 75.8


# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA (with caching)
# ══════════════════════════════════════════════════════════════════════════════

if CACHE.exists():
    print(f"Loading cached data from {CACHE.name}...")
    with open(CACHE, 'rb') as f:
        cache = pickle.load(f)
    strict = cache['strict']
    print(f"  {len(strict)} valid end-of-IOB samples loaded from cache")
else:
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
    print("BOOST SINGLE-PATIENT — END-OF-IOB ANALYSIS")
    print("=" * 70)

    print("\nLoading devicestatus...")
    raw_ds = load_dedup(find_json('devicestatus'))
    print("\nLoading CGM entries...")
    raw_entries = load_dedup(find_json('entries'))
    print("\nLoading treatments...")
    raw_tx = load_dedup(find_json('treatments'))

    # ── Compute actual TDD from treatments ──
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
    daily_tdd = df_insulin.groupby('date')['insulin'].sum().reset_index().rename(columns={'insulin': 'tdd_actual'})
    daily_tdd = daily_tdd.sort_values('date')
    daily_tdd['tdd_7day'] = daily_tdd['tdd_actual'].rolling(7, min_periods=3).mean()
    daily_tdd['tdd_1day'] = daily_tdd['tdd_actual']

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
            return tdd_weighted
        if (tdd_weighted < (0.75 * tdd_7d)) and (tdd_1d is not None):
            tdd = ((tdd_weighted + ((tdd_weighted / tdd_7d) * (tdd_7d - tdd_weighted))) * 0.34) + (tdd_1d * 0.33) + (tdd_weighted * 0.33)
        elif tdd_1d is not None:
            tdd = (tdd_weighted * 0.33) + (tdd_7d * 0.34) + (tdd_1d * 0.33)
        else:
            tdd = tdd_weighted
        return tdd

    # ── Parse loop cycles — store full pred_iob array ──
    print("\n── Parsing loop cycles (end-of-IOB) ──")
    cycles = []
    for r in raw_ds:
        try:
            ts = pd.to_datetime(r['created_at'], utc=True)
            sg = r.get('openaps', {}).get('suggested', {})
            if not sg or 'bg' not in sg: continue
            pred_iob = sg.get('predBGs', {}).get('IOB') or []
            if len(pred_iob) < 12:  # Need at least 1h of prediction
                continue
            cycles.append({
                'ts': ts, 'bg': sg.get('bg'),
                'variable_sens': sg.get('variable_sens'),
                'cob': sg.get('COB'), 'iob': sg.get('IOB'),
                'pred_iob_final': pred_iob[-1],
                'pred_horizon_s': (len(pred_iob) - 1) * 5 * 60,
            })
        except:
            continue

    df = pd.DataFrame(cycles).dropna(subset=['ts', 'bg', 'variable_sens']).sort_values('ts').reset_index(drop=True)
    print(f"  Parsed: {len(df):,} cycles  ({df['ts'].min().date()} → {df['ts'].max().date()})")
    pred_horizons = df['pred_horizon_s'].values / 60
    print(f"  Prediction horizon: median={np.median(pred_horizons):.0f} min, "
          f"range={pred_horizons.min():.0f}-{pred_horizons.max():.0f} min")

    # ── CGM + forward lookup at dynamic horizon ──
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

    print("\n── Forward SGV lookups at dynamic horizon ──")
    cycle_epochs = np.array([int(t.timestamp()) for t in df['ts']])
    actual_end = np.full(len(df), np.nan)
    bolus_age = np.full(len(df), np.nan)

    for i, t_s in enumerate(cycle_epochs):
        horizon = int(df.iloc[i]['pred_horizon_s'])
        actual_end[i] = get_cgm_at(t_s + horizon)
        bolus_age[i] = mins_since_bolus(t_s)
        if (i + 1) % 10000 == 0:
            print(f"  ... {i+1:,}/{len(df):,}")

    df['actual_bg_end'] = actual_end
    df['bolus_age_min'] = bolus_age
    df['hour'] = df['ts'].dt.hour
    df['date'] = df['ts'].dt.date

    print(f"  End-of-IOB coverage: {(~np.isnan(actual_end)).sum():,}/{len(df):,}")

    # ── Attach TDD ──
    print("\n── Attaching TDD ──")
    daily_tdd['date'] = pd.to_datetime(daily_tdd['date']).dt.date
    df = df.merge(daily_tdd[['date', 'tdd_actual', 'tdd_7day', 'tdd_1day']], on='date', how='left')

    tdd_boost = np.full(len(df), np.nan)
    print("  Computing Boost blend per cycle...")
    for i, (epoch, tdd7, tdd1) in enumerate(zip(cycle_epochs, df['tdd_7day'].values, df['tdd_1day'].values)):
        tdd7_v = tdd7 if not np.isnan(tdd7) else None
        tdd1_v = tdd1 if not np.isnan(tdd1) else None
        tdd_boost[i] = compute_boost_tdd(epoch, tdd7_v, tdd1_v)
        if (i + 1) % 10000 == 0:
            print(f"  ... {i+1:,}/{len(df):,}")

    df['tdd_boost'] = tdd_boost

    # ── Overnight fasting filter ──
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

    on = on.dropna(subset=['pred_iob_final', 'actual_bg_end']).copy()
    on['pred_drop'] = on['bg'] - on['pred_iob_final']
    on['actual_drop'] = on['bg'] - on['actual_bg_end']

    strict = on[
        (np.abs(on['pred_drop']) > 3) &
        ((on['actual_drop'] / on['pred_drop']).between(0, 5))
    ].copy().sort_values('ts').reset_index(drop=True)
    print(f"  After strict filtering: {len(strict):,} valid end-of-IOB samples")

    # Cache
    with open(CACHE, 'wb') as f:
        pickle.dump({'strict': strict}, f)
    print(f"  Cached to {CACHE.name}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. FORMULA COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("BOOST SINGLE-PATIENT — END-OF-IOB RESULTS")
print("=" * 80)

bg = strict['bg'].values
isf_actual = strict['variable_sens'].values
pred_drop = strict['pred_drop'].values
actual_bg_end = strict['actual_bg_end'].values
tdd_7d = strict['tdd_7day'].values
tdd_b = strict['tdd_boost'].values
timestamps = strict['ts'].values

ln_target = np.log(TARGET / D + 1)
ln_bg = np.log(bg / D + 1)

median_horizon_min = np.median(strict['pred_horizon_s'].values) / 60
real_tdd_median = np.median(tdd_7d)
real_isf_at_target = np.median(1800.0 / tdd_b)

print(f"\n  Samples: {len(strict):,}")
print(f"  Prediction horizon: median {median_horizon_min:.0f} min")
print(f"  TDD 7-day median: {real_tdd_median:.1f} U/day")
print(f"  Boost TDD median: {np.median(tdd_b):.1f} U/day")

# Compute ISF formulas
# A: Loop actual
# B: Current ln + Boost TDD
isf_current_boost = 1800.0 / (tdd_b * ln_bg)
# C: Power-law k=3.5 + Boost TDD
isf_pl_boost = (1800.0 / (tdd_b * ln_target)) * (TARGET / bg) ** 3.5
# D: Quartic (TDD-scaled)
S_q = (1800.0 / real_tdd_median) / QUARTIC_AT_99
isf_q = isf_quartic(bg) * S_q
# E: Full Diabeloop (TDD-scaled)
S_db = (1800.0 / real_tdd_median) / FULL_DB_AT_99
isf_db = isf_full_diabeloop(bg) * S_db
# F: Hybrid (TDD-scaled)
S_h = (1800.0 / real_tdd_median) / HYBRID_AT_105
isf_h = isf_hybrid(bg) * S_h

# Counterfactual predictions
pred_loop = bg - pred_drop * 1.0
pred_ln_boost = bg - pred_drop * (isf_current_boost / isf_actual)
pred_pl_boost = bg - pred_drop * (isf_pl_boost / isf_actual)
pred_q = bg - pred_drop * (isf_q / isf_actual)
pred_db = bg - pred_drop * (isf_db / isf_actual)
pred_h = bg - pred_drop * (isf_h / isf_actual)

formulas = [
    ('Loop actual',           pred_loop),
    ('Current ln + Boost',    pred_ln_boost),
    ('Power-law k=3.5+Boost', pred_pl_boost),
    ('Quartic (TDD-scaled)',  pred_q),
    ('Full Diabeloop',        pred_db),
    ('Hybrid (poly+PL)',      pred_h),
]

print(f"\n  {'Formula':<28s}  {'MAE':>5s}  {'Bias':>6s}")
print(f"  {'─'*28}  {'─'*5}  {'─'*6}")
for label, pred in formulas:
    err = actual_bg_end - pred
    valid = ~np.isnan(err)
    mae = np.abs(err[valid]).mean()
    bias = np.mean(pred[valid] - actual_bg_end[valid])
    print(f"  {label:<28s}  {mae:5.1f}  {bias:+6.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 3. FALLING vs RISING DIRECTIONAL ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("FALLING vs RISING — Pred−Actual (positive = over-predicts SGV)")
print("=" * 80)

combos = [
    ("<105 falling",  lambda b: b < 105,  lambda p: p > 0),
    ("<105 rising",   lambda b: b < 105,  lambda p: p < 0),
    ("≥105 falling",  lambda b: b >= 105, lambda p: p > 0),
    ("≥105 rising",   lambda b: b >= 105, lambda p: p < 0),
]

formula_preds = [
    ('Loop', pred_loop), ('Ln+Boost', pred_ln_boost), ('PL+Boost', pred_pl_boost),
    ('Quartic', pred_q), ('FullDB', pred_db), ('Hybrid', pred_h),
]

print(f"\n  {'Zone':<16s} {'N':>6s}", end='')
for fname, _ in formula_preds:
    print(f"  {fname:>10s}", end='')
print()
print(f"  {'─'*16} {'─'*6}", end='')
for _ in formula_preds:
    print(f"  {'─'*10}", end='')
print()

for label, bg_fn, dr_fn in combos:
    m = bg_fn(bg) & dr_fn(pred_drop)
    n = m.sum()
    if n == 0:
        print(f"  {label:<16s}    —")
        continue
    act_m = actual_bg_end[m]
    print(f"  {label:<16s} {n:6d}", end='')
    for fname, pred in formula_preds:
        bias = np.mean(pred[m]) - np.mean(act_m)
        print(f"  {bias:+10.1f}", end='')
    print()


# ══════════════════════════════════════════════════════════════════════════════
# 4. BG-BAND ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'─'*80}")
print("SGV-Band MAE and Bias (Pred−Actual)")
print(f"{'─'*80}")

bands = [('<90', bg < 90), ('90-105', (bg >= 90) & (bg < 105)),
         ('105-120', (bg >= 105) & (bg < 120)), ('120-150', (bg >= 120) & (bg <= 150))]

print(f"\n  {'Band':>8s} {'N':>5s}", end='')
for fname, _ in formula_preds:
    print(f"  {fname+' MAE':>12s} {fname+' Bias':>12s}", end='')
print()

for bname, bmask in bands:
    n = bmask.sum()
    act_b = actual_bg_end[bmask]
    print(f"  {bname:>8s} {n:5d}", end='')
    for fname, pred in formula_preds:
        err = pred[bmask] - act_b
        valid = ~np.isnan(err)
        mae = np.abs(err[valid]).mean()
        bias = err[valid].mean()
        print(f"  {mae:12.1f} {bias:+12.1f}", end='')
    print()


# ══════════════════════════════════════════════════════════════════════════════
# 5. BELOW-90 SAFETY DETAIL
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'─'*80}")
print("BELOW 90 — Predicted−Actual end-of-IOB SGV")
print(f"{'─'*80}")
m90 = bg < 90
n90 = m90.sum()
act90 = actual_bg_end[m90]
print(f"\n  Starting SGV < 90: n={n90}, mean actual end-of-IOB = {np.nanmean(act90):.0f}")
for fname, pred in formula_preds:
    bias = np.mean(pred[m90]) - np.mean(act90)
    print(f"  {fname:<28s}  Pred−Actual: {bias:+.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 6. SYNTHETIC MULTI-PATIENT VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("SYNTHETIC MULTI-PATIENT VALIDATION — END-OF-IOB")
print("=" * 80)

patients = [
    {'id': 'P1', 'label': 'Very resistant',     'scale': 0.50},
    {'id': 'P2', 'label': 'Resistant',           'scale': 0.67},
    {'id': 'P3', 'label': 'Moderately resistant', 'scale': 0.80},
    {'id': 'P4', 'label': 'Original patient',    'scale': 1.00},
    {'id': 'P5', 'label': 'Moderately sensitive', 'scale': 1.25},
    {'id': 'P6', 'label': 'Sensitive',           'scale': 1.60},
    {'id': 'P7', 'label': 'Very sensitive',      'scale': 2.00},
]


def compute_errors_synth(isf_formula, isf_true):
    pred_f = bg - pred_drop * (isf_formula / isf_true)
    err = pred_f - actual_bg_end  # Pred - Actual convention
    valid = ~np.isnan(err)
    e = err[valid]
    return np.abs(e).mean(), e.mean(), e


def simulate_phase1(patient_scale, isf_func, anchor_value):
    """Phase 1 (TDD-based) calibration for synthetic patient at end-of-IOB."""
    isf_true = isf_actual * patient_scale
    patient_tdd = real_tdd_median / patient_scale
    S_phase1 = (1800.0 / patient_tdd) / anchor_value

    isf_p1 = isf_func(bg) * S_phase1
    mae_p1, bias_p1, errors_p1 = compute_errors_synth(isf_p1, isf_true)

    # Uncalibrated
    isf_raw = isf_func(bg)
    mae_raw, bias_raw, _ = compute_errors_synth(isf_raw, isf_true)

    # Band analysis
    band_results = {}
    for bname, bmask in bands:
        isf_p1_b = isf_func(bg[bmask]) * S_phase1
        isf_true_b = isf_true[bmask]
        pred_b = bg[bmask] - pred_drop[bmask] * (isf_p1_b / isf_true_b)
        err_b = pred_b - actual_bg_end[bmask]
        valid = ~np.isnan(err_b)
        if valid.sum() > 0:
            band_results[bname] = {'mae': np.abs(err_b[valid]).mean(), 'bias': err_b[valid].mean(), 'n': valid.sum()}

    # Falling/rising for <105
    m_rising = (bg < 105) & (pred_drop < 0)
    if m_rising.sum() > 0:
        isf_p1_r = isf_func(bg[m_rising]) * S_phase1
        isf_true_r = isf_true[m_rising]
        pred_r = bg[m_rising] - pred_drop[m_rising] * (isf_p1_r / isf_true_r)
        bias_rising = np.mean(pred_r) - np.mean(actual_bg_end[m_rising])
    else:
        bias_rising = np.nan

    return {
        'tdd': patient_tdd, 'S_phase1': S_phase1,
        'mae_raw': mae_raw, 'bias_raw': bias_raw,
        'mae_p1': mae_p1, 'bias_p1': bias_p1,
        'bands': band_results,
        'bias_rising_105': bias_rising,
    }


synth_formulas = [
    ('Quartic',  isf_quartic,        QUARTIC_AT_99),
    ('FullDB',   isf_full_diabeloop, FULL_DB_AT_99),
    ('Hybrid',   isf_hybrid,         HYBRID_AT_105),
]

for formula_name, isf_func, anchor in synth_formulas:
    print(f"\n── {formula_name} — anchor={anchor:.1f} ──")
    print(f"  {'Patient':<6s}  {'Type':<22s}  {'TDD':>5s}  {'S':>6s}  "
          f"{'MAE_raw':>7s}  {'MAE_p1':>7s}  {'Bias_p1':>8s}  {'<105 ris':>8s}")
    print(f"  {'─'*6}  {'─'*22}  {'─'*5}  {'─'*6}  {'─'*7}  {'─'*7}  {'─'*8}  {'─'*8}")

    for p in patients:
        r = simulate_phase1(p['scale'], isf_func, anchor)
        print(f"  {p['id']:<6s}  {p['label']:<22s}  {r['tdd']:5.1f}  {r['S_phase1']:6.3f}  "
              f"{r['mae_raw']:7.1f}  {r['mae_p1']:7.1f}  {r['bias_p1']:+8.1f}  {r['bias_rising_105']:+8.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 7. COMPARISON TABLE: 2h vs End-of-IOB
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("KEY COMPARISON: BOOST SINGLE-PATIENT")
print(f"  End-of-IOB horizon: {median_horizon_min:.0f} min")
print(f"  Samples: {len(strict):,}")
print("=" * 80)

print(f"\n  Directional bias (Pred−Actual):")
for label, bg_fn, dr_fn in combos:
    m = bg_fn(bg) & dr_fn(pred_drop)
    n = m.sum()
    if n == 0: continue
    act_m = actual_bg_end[m]
    loop_b = np.mean(pred_loop[m]) - np.mean(act_m)
    q_b = np.mean(pred_q[m]) - np.mean(act_m)
    db_b = np.mean(pred_db[m]) - np.mean(act_m)
    h_b = np.mean(pred_h[m]) - np.mean(act_m)
    print(f"  {label:<16s} n={n:5d}  Loop:{loop_b:+6.1f}  Quartic:{q_b:+6.1f}  FullDB:{db_b:+6.1f}  Hybrid:{h_b:+6.1f}")


print("\n" + "=" * 80)
print("DONE — Boost end-of-IOB analysis complete")
print("=" * 80)
