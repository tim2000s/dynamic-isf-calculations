#!/usr/bin/env python3
"""
Full Polynomial vs Hybrid — Calibration Backtest
==================================================

Compares two ISF curve shapes through the same calibration pipeline:
  A: Full Diabeloop polynomial (all glucose values) × S
  B: Hybrid (polynomial ≥105, power-law <105) × S

Both use identical Phase 1 calibration: S = (1800 / TDD_7day) / anchor
Both tested across 7 synthetic patients (TDD 11–45 U/day).

The key question: does the polynomial's smooth, reasonable curve below 100
perform as well as the hybrid's aggressive power-law tail?
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
from pathlib import Path

warnings.filterwarnings('ignore')

D      = 82.0
TARGET = 99.0

HOME    = Path.home()
NS_WORK = HOME / 'Nightscout_Work'
OUT_DIR = HOME / 'Downloads'


# ── ISF formulas ──
def isf_polynomial(bg):
    """Full Diabeloop polynomial across all glucose values"""
    G = np.asarray(bg, dtype=float)
    return 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4


def isf_hybrid(bg):
    """Hybrid: polynomial ≥105, power-law <105"""
    G = np.asarray(bg, dtype=float)
    poly = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    power = 75.8 * (105.0 / G) ** 3.5
    return np.where(G >= 105, poly, power)


# Anchor points for S calibration: polynomial value at target
POLY_AT_TARGET = float(isf_polynomial(99.0))   # ~81.6
HYBRID_AT_105 = 75.8


# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA (same pipeline as polynomial backtest)
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
print("HYBRID ISF CALIBRATION BACKTEST — SYNTHETIC MULTI-PATIENT")
print("=" * 70)

print("\nLoading devicestatus...")
raw_ds = load_dedup(find_json('devicestatus'))

print("\nLoading CGM entries...")
raw_entries = load_dedup(find_json('entries'))

print("\nLoading treatments...")
raw_tx = load_dedup(find_json('treatments'))


# ── Compute actual TDD ──
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


# ── Parse loop cycles ──
print("\n── Parsing loop cycles ──")
cycles = []
for r in raw_ds:
    try:
        ts = pd.to_datetime(r['created_at'], utc=True)
        sg = r.get('openaps', {}).get('suggested', {})
        if not sg or 'bg' not in sg: continue
        pred_iob = sg.get('predBGs', {}).get('IOB') or []

        cycles.append({
            'ts': ts, 'bg': sg.get('bg'),
            'variable_sens': sg.get('variable_sens'),
            'cob': sg.get('COB'), 'iob': sg.get('IOB'),
            'pred_iob_24': pred_iob[24] if len(pred_iob) > 24 else None,
        })
    except:
        continue

df = pd.DataFrame(cycles).dropna(subset=['ts', 'bg', 'variable_sens']).sort_values('ts').reset_index(drop=True)
print(f"  Parsed: {len(df):,} cycles  ({df['ts'].min().date()} → {df['ts'].max().date()})")


# ── CGM + forward lookup ──
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


print("\n── Forward BG lookups ──")
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
strict = strict.sort_values('ts').reset_index(drop=True)
print(f"  After strict filtering: {len(strict):,} valid 2h samples")


# ══════════════════════════════════════════════════════════════════════════════
# 2. SYNTHETIC PATIENTS
# ══════════════════════════════════════════════════════════════════════════════

# Real patient baseline
bg = strict['bg'].values
isf_actual = strict['variable_sens'].values
pred_drop = strict['pred_drop'].values
actual_bg_2h = strict['actual_bg_2h'].values
tdd_7d = strict['tdd_7day'].values
timestamps = strict['ts'].values

real_tdd_median = np.median(tdd_7d)
real_isf_at_target = np.median(1800.0 / strict['tdd_boost'].values)

print(f"\n── Real patient baseline ──")
print(f"  Median 7D TDD: {real_tdd_median:.1f} U/day")
print(f"  Median ISF at target: {real_isf_at_target:.1f} mg/dL/U")

# Synthetic patient definitions
# scale_factor multiplies the patient's TRUE ISF (higher = more sensitive)
# TDD scales inversely (higher sensitivity = lower TDD)
patients = [
    {'id': 'P1', 'label': 'Very resistant',     'scale': 0.50},
    {'id': 'P2', 'label': 'Resistant',           'scale': 0.67},
    {'id': 'P3', 'label': 'Moderately resistant', 'scale': 0.80},
    {'id': 'P4', 'label': 'Original patient',    'scale': 1.00},
    {'id': 'P5', 'label': 'Moderately sensitive', 'scale': 1.25},
    {'id': 'P6', 'label': 'Sensitive',           'scale': 1.60},
    {'id': 'P7', 'label': 'Very sensitive',      'scale': 2.00},
]


def compute_errors_with_isf(isf_formula, isf_true):
    """Counterfactual: what if the loop had used isf_formula instead of isf_true?"""
    pred_f = bg - pred_drop * (isf_formula / isf_true)
    err = actual_bg_2h - pred_f
    valid = ~np.isnan(err)
    e = err[valid]
    mae = np.abs(e).mean()
    bias = e.mean()
    w18 = (np.abs(e) <= 18).mean() * 100
    return mae, bias, w18, e


# ══════════════════════════════════════════════════════════════════════════════
# 3. CALIBRATION SIMULATION
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'═'*80}")
print("CALIBRATION SIMULATION")
print(f"{'═'*80}")

# Time-based sample grouping for progressive calibration
# Group samples by date to simulate day-by-day auto-calibration
strict_dates = pd.to_datetime(strict['ts']).dt.date.values
unique_dates = np.unique(strict_dates)
n_dates = len(unique_dates)
print(f"\n  Data spans {n_dates} dates with fasting samples")

# For each date, find sample indices
date_to_idx = {}
for i, d in enumerate(strict_dates):
    date_to_idx.setdefault(d, []).append(i)


def simulate_calibration(patient_scale, isf_func, anchor_value):
    """
    Simulate the two-phase calibration for a synthetic patient.

    isf_func: the ISF curve function (isf_polynomial or isf_hybrid)
    anchor_value: the unscaled ISF value used for S calculation
    """
    # True ISF for this synthetic patient
    isf_true = isf_actual * patient_scale

    # True scaling factor (what auto-cal should converge to)
    true_isf_at_target = real_isf_at_target * patient_scale
    S_true = true_isf_at_target / anchor_value

    # Phase 1: from synthetic TDD
    patient_tdd = real_tdd_median / patient_scale
    S_phase1 = (1800.0 / patient_tdd) / anchor_value

    # Phase 1 errors
    isf_p1 = isf_func(bg) * S_phase1
    mae_p1, bias_p1, w18_p1, _ = compute_errors_with_isf(isf_p1, isf_true)

    # Phase 2: auto-calibration simulation
    S_current = S_phase1
    S_history = [S_current]
    mae_history = [mae_p1]
    bias_history = [bias_p1]
    calibration_ratios = []
    ROLLING_DAYS = 7

    for day_idx, date in enumerate(unique_dates):
        day_samples = date_to_idx[date]
        day_bg = bg[day_samples]
        day_isf_true = isf_true[day_samples]
        day_pred_drop = pred_drop[day_samples]
        day_actual_2h = actual_bg_2h[day_samples]
        day_isf_formula = isf_func(day_bg) * S_current

        day_pred_bg = day_bg - day_pred_drop * (day_isf_formula / day_isf_true)
        day_actual_drop = day_bg - day_actual_2h
        day_pred_drop_formula = day_bg - day_pred_bg

        valid = (np.abs(day_pred_drop_formula) > 3) & (~np.isnan(day_actual_drop))
        if valid.sum() > 0:
            ratios = day_actual_drop[valid] / day_pred_drop_formula[valid]
            sane = (ratios > 0.1) & (ratios < 10)
            if sane.sum() > 0:
                calibration_ratios.append((date, ratios[sane]))

        if len(calibration_ratios) > ROLLING_DAYS:
            calibration_ratios = calibration_ratios[-ROLLING_DAYS:]

        DAMPING = 0.30
        if len(calibration_ratios) >= 2:
            all_ratios = np.concatenate([r for _, r in calibration_ratios])
            median_ratio = np.median(all_ratios)
            S_indicated = S_current * median_ratio
            S_new = S_current + DAMPING * (S_indicated - S_current)
            S_current = np.clip(S_new, 0.5, 5.0)

        S_history.append(S_current)

        isf_current = isf_func(bg) * S_current
        mae_cur, bias_cur, _, _ = compute_errors_with_isf(isf_current, isf_true)
        mae_history.append(mae_cur)
        bias_history.append(bias_cur)

    # Final steady-state metrics
    isf_final = isf_func(bg) * S_current
    mae_final, bias_final, w18_final, _ = compute_errors_with_isf(isf_final, isf_true)

    # Uncalibrated (S=1.0)
    isf_raw = isf_func(bg)
    mae_raw, bias_raw, w18_raw, _ = compute_errors_with_isf(isf_raw, isf_true)

    # Perfect calibration
    isf_perfect = isf_func(bg) * S_true
    mae_perfect, bias_perfect, w18_perfect, _ = compute_errors_with_isf(isf_perfect, isf_true)

    return {
        'S_true': S_true, 'S_phase1': S_phase1, 'S_final': S_current,
        'S_history': S_history, 'mae_history': mae_history, 'bias_history': bias_history,
        'tdd_synthetic': real_tdd_median / patient_scale,
        'mae_raw': mae_raw, 'bias_raw': bias_raw, 'w18_raw': w18_raw,
        'mae_p1': mae_p1, 'bias_p1': bias_p1, 'w18_p1': w18_p1,
        'mae_final': mae_final, 'bias_final': bias_final, 'w18_final': w18_final,
        'mae_perfect': mae_perfect, 'bias_perfect': bias_perfect, 'w18_perfect': w18_perfect,
    }


# ── Run calibration for each patient × each formula ──
formulas_to_test = [
    ('poly', 'Full Polynomial', isf_polynomial, POLY_AT_TARGET),
    ('hybrid', 'Hybrid (poly≥105 + PL<105)', isf_hybrid, HYBRID_AT_105),
]

all_results = {}
for formula_key, formula_name, isf_func, anchor in formulas_to_test:
    print(f"\n── {formula_name} — anchor={anchor:.1f} ──\n")
    all_results[formula_key] = {}
    for p in patients:
        r = simulate_calibration(p['scale'], isf_func, anchor)
        all_results[formula_key][p['id']] = r
        print(f"  {p['id']} ({p['label']}, scale={p['scale']:.2f}):")
        print(f"    TDD={r['tdd_synthetic']:.1f}  S_true={r['S_true']:.3f}  "
              f"S_ph1={r['S_phase1']:.3f}  S_fin={r['S_final']:.3f}")
        print(f"    MAE: raw={r['mae_raw']:.1f}  ph1={r['mae_p1']:.1f}  "
              f"cal={r['mae_final']:.1f}  prf={r['mae_perfect']:.1f}")
        print(f"    Bias: raw={r['bias_raw']:+.1f}  ph1={r['bias_p1']:+.1f}  "
              f"cal={r['bias_final']:+.1f}  prf={r['bias_perfect']:+.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 4. SUMMARY TABLES
# ══════════════════════════════════════════════════════════════════════════════

for formula_key, formula_name, _, _ in formulas_to_test:
    results = all_results[formula_key]
    print(f"\n{'═'*100}")
    print(f"SUMMARY — {formula_name.upper()}")
    print(f"{'═'*100}")

    print(f"\n  {'Patient':<6s}  {'Type':<22s}  {'Scale':>5s}  {'TDD':>5s}  "
          f"{'S_true':>6s}  {'S_ph1':>6s}  {'S_fin':>6s}  "
          f"{'MAE_raw':>7s}  {'MAE_p1':>7s}  {'MAE_cal':>7s}  {'MAE_prf':>7s}")
    print(f"  {'─'*6}  {'─'*22}  {'─'*5}  {'─'*5}  "
          f"{'─'*6}  {'─'*6}  {'─'*6}  "
          f"{'─'*7}  {'─'*7}  {'─'*7}  {'─'*7}")

    for p in patients:
        r = results[p['id']]
        print(f"  {p['id']:<6s}  {p['label']:<22s}  {p['scale']:5.2f}  {r['tdd_synthetic']:5.1f}  "
              f"{r['S_true']:6.3f}  {r['S_phase1']:6.3f}  {r['S_final']:6.3f}  "
              f"{r['mae_raw']:7.1f}  {r['mae_p1']:7.1f}  {r['mae_final']:7.1f}  {r['mae_perfect']:7.1f}")

# ── Head-to-head comparison ──
print(f"\n{'═'*100}")
print("HEAD-TO-HEAD: FULL POLYNOMIAL vs HYBRID — Phase 1 (TDD-based)")
print(f"{'═'*100}")
print(f"\n  {'Patient':<6s}  {'TDD':>5s}  "
      f"{'Poly MAE':>8s}  {'Poly Bias':>9s}  "
      f"{'Hybr MAE':>8s}  {'Hybr Bias':>9s}  {'Δ MAE':>6s}")
print(f"  {'─'*6}  {'─'*5}  {'─'*8}  {'─'*9}  {'─'*8}  {'─'*9}  {'─'*6}")

for p in patients:
    rp = all_results['poly'][p['id']]
    rh = all_results['hybrid'][p['id']]
    delta = rp['mae_p1'] - rh['mae_p1']
    print(f"  {p['id']:<6s}  {rp['tdd_synthetic']:5.1f}  "
          f"{rp['mae_p1']:8.1f}  {rp['bias_p1']:+9.1f}  "
          f"{rh['mae_p1']:8.1f}  {rh['bias_p1']:+9.1f}  {delta:+6.1f}")

# ── BG-band head-to-head ──
print(f"\n── Glucose-Band MAE: Full Polynomial vs Hybrid (Phase 1, P4=original patient) ──\n")
bands = [('<90', bg < 90), ('90-105', (bg >= 90) & (bg < 105)),
         ('105-120', (bg >= 105) & (bg < 120)), ('120-150', (bg >= 120) & (bg <= 150))]

for formula_key, formula_name, isf_func, anchor in formulas_to_test:
    r = all_results[formula_key]['P4']
    isf_true = isf_actual * 1.0  # P4 = original
    isf_cal = isf_func(bg) * r['S_phase1']
    pred_f = bg - pred_drop * (isf_cal / isf_true)
    err = actual_bg_2h - pred_f

    print(f"  {formula_name}  (S={r['S_phase1']:.3f}):")
    for bname, bmask in bands:
        e = err[bmask & ~np.isnan(err)]
        if len(e) > 10:
            print(f"    {bname:>8s}  n={len(e):4d}  MAE={np.abs(e).mean():5.1f}  bias={e.mean():+6.1f}")
    print()


# ══════════════════════════════════════════════════════════════════════════════
# 5. FIGURE — Head-to-head comparison
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


fig = plt.figure(figsize=(24, 20))
fig.patch.set_facecolor(BG_C)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.42, wspace=0.35)

x = np.arange(len(patients))
poly_color = '#80cbc4'
hybrid_color = '#f48fb1'

# P1: ISF curves — full poly vs hybrid (both TDD-scaled for P4)
ax1 = fig.add_subplot(gs[0, 0])
style(ax1, 'ISF Curves: Full Polynomial vs Hybrid (P4 Phase 1 S)')
bg_range = np.linspace(60, 200, 300)
S_poly_p4 = all_results['poly']['P4']['S_phase1']
S_hybrid_p4 = all_results['hybrid']['P4']['S_phase1']
ax1.plot(bg_range, isf_polynomial(bg_range) * S_poly_p4, lw=2.5, color=poly_color,
         label=f'Full Poly ×{S_poly_p4:.2f}')
ax1.plot(bg_range, isf_hybrid(bg_range) * S_hybrid_p4, lw=2.5, color=hybrid_color,
         label=f'Hybrid ×{S_hybrid_p4:.2f}')
ax1.axvline(105, color='white', lw=0.8, ls=':', alpha=0.4)
ax1.axvline(TARGET, color='yellow', lw=0.8, ls=':', alpha=0.4, label=f'Target {TARGET:.0f}')
ax1.set_xlabel('Glucose (mg/dL)'); ax1.set_ylabel('ISF (mg/dL/U)')
ax1.set_xlim(60, 200); ax1.set_ylim(0, 250)
ax1.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P2: ISF curves zoomed into sub-105
ax2 = fig.add_subplot(gs[0, 1])
style(ax2, 'ISF Below 105 — Full Poly vs Hybrid (P4)')
bg_low = np.linspace(60, 110, 200)
ax2.plot(bg_low, isf_polynomial(bg_low) * S_poly_p4, lw=2.5, color=poly_color,
         label='Full Poly')
ax2.plot(bg_low, isf_hybrid(bg_low) * S_hybrid_p4, lw=2.5, color=hybrid_color,
         label='Hybrid')
ax2.axvline(90, color='#ff6e6e', lw=0.8, ls='--', alpha=0.6, label='90 mg/dL')
ax2.axvline(105, color='white', lw=0.8, ls=':', alpha=0.4, label='Junction 105')
ax2.set_xlabel('Glucose (mg/dL)'); ax2.set_ylabel('ISF (mg/dL/U)')
ax2.set_xlim(60, 110); ax2.set_ylim(50, 350)
ax2.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P3: Phase 1 MAE head-to-head
ax3 = fig.add_subplot(gs[0, 2])
style(ax3, 'Phase 1 MAE: Full Polynomial vs Hybrid')
w = 0.35
mae_poly = [all_results['poly'][p['id']]['mae_p1'] for p in patients]
mae_hybr = [all_results['hybrid'][p['id']]['mae_p1'] for p in patients]
ax3.bar(x - w/2, mae_poly, w, color=poly_color, alpha=0.85, label='Full Polynomial')
ax3.bar(x + w/2, mae_hybr, w, color=hybrid_color, alpha=0.85, label='Hybrid')
for xi, (v1, v2) in enumerate(zip(mae_poly, mae_hybr)):
    ax3.text(xi - w/2, v1 + 0.2, f'{v1:.1f}', ha='center', fontsize=7, color=TXT)
    ax3.text(xi + w/2, v2 + 0.2, f'{v2:.1f}', ha='center', fontsize=7, color=TXT)
ax3.set_xticks(x); ax3.set_xticklabels([p['id'] for p in patients], fontsize=8)
ax3.set_xlabel('Patient'); ax3.set_ylabel('MAE (mg/dL)')
ax3.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P4: Phase 1 Bias head-to-head
ax4 = fig.add_subplot(gs[1, 0])
style(ax4, 'Phase 1 Bias: Full Polynomial vs Hybrid')
bias_poly = [all_results['poly'][p['id']]['bias_p1'] for p in patients]
bias_hybr = [all_results['hybrid'][p['id']]['bias_p1'] for p in patients]
ax4.bar(x - w/2, bias_poly, w, color=poly_color, alpha=0.85, label='Full Polynomial')
ax4.bar(x + w/2, bias_hybr, w, color=hybrid_color, alpha=0.85, label='Hybrid')
ax4.axhline(0, color='white', lw=0.8, ls='--')
for xi, (v1, v2) in enumerate(zip(bias_poly, bias_hybr)):
    ax4.text(xi - w/2, v1 + 0.3, f'{v1:+.1f}', ha='center', fontsize=7, color=TXT)
    ax4.text(xi + w/2, v2 + 0.3, f'{v2:+.1f}', ha='center', fontsize=7, color=TXT)
ax4.set_xticks(x); ax4.set_xticklabels([p['id'] for p in patients], fontsize=8)
ax4.set_xlabel('Patient'); ax4.set_ylabel('Bias (mg/dL)')
ax4.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P5: Glucose-band MAE comparison (P4)
ax5 = fig.add_subplot(gs[1, 1])
style(ax5, 'Glucose-Band MAE (P4 Original Patient, Phase 1)')
band_names = ['<90', '90-105', '105-120', '120-150']
band_masks_list = [bg < 90, (bg >= 90) & (bg < 105), (bg >= 105) & (bg < 120), (bg >= 120) & (bg <= 150)]
bx = np.arange(len(band_names))

for fi, (fkey, fname, isf_func, anchor) in enumerate(formulas_to_test):
    r = all_results[fkey]['P4']
    isf_true = isf_actual
    isf_cal = isf_func(bg) * r['S_phase1']
    pred_f = bg - pred_drop * (isf_cal / isf_true)
    err = actual_bg_2h - pred_f
    band_maes = []
    for bmask in band_masks_list:
        e = err[bmask & ~np.isnan(err)]
        band_maes.append(np.abs(e).mean() if len(e) > 10 else np.nan)
    clr = poly_color if fkey == 'poly' else hybrid_color
    ax5.bar(bx + (fi - 0.5)*0.35, band_maes, 0.35, color=clr, alpha=0.85,
            label=fname[:20])
ax5.set_xticks(bx); ax5.set_xticklabels(band_names, fontsize=8)
ax5.set_xlabel('Glucose Band (mg/dL)'); ax5.set_ylabel('MAE (mg/dL)')
ax5.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P6: Glucose-band Bias comparison (P4)
ax6 = fig.add_subplot(gs[1, 2])
style(ax6, 'Glucose-Band Bias (P4 Original Patient, Phase 1)')
for fi, (fkey, fname, isf_func, anchor) in enumerate(formulas_to_test):
    r = all_results[fkey]['P4']
    isf_cal = isf_func(bg) * r['S_phase1']
    pred_f = bg - pred_drop * (isf_cal / isf_actual)
    err = actual_bg_2h - pred_f
    band_biases = []
    for bmask in band_masks_list:
        e = err[bmask & ~np.isnan(err)]
        band_biases.append(e.mean() if len(e) > 10 else np.nan)
    clr = poly_color if fkey == 'poly' else hybrid_color
    ax6.bar(bx + (fi - 0.5)*0.35, band_biases, 0.35, color=clr, alpha=0.85,
            label=fname[:20])
ax6.axhline(0, color='white', lw=0.8, ls='--')
ax6.set_xticks(bx); ax6.set_xticklabels(band_names, fontsize=8)
ax6.set_xlabel('Glucose Band (mg/dL)'); ax6.set_ylabel('Bias (mg/dL)')
ax6.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P7: MAE heatmap — poly vs hybrid across all patients (Phase 1)
ax7 = fig.add_subplot(gs[2, 0])
style(ax7, 'Phase 1 MAE Heatmap: Full Poly')
hm_poly = []
for p in patients:
    r = all_results['poly'][p['id']]
    isf_true = isf_actual * p['scale']
    isf_cal = isf_polynomial(bg) * r['S_phase1']
    pred_f = bg - pred_drop * (isf_cal / isf_true)
    err = actual_bg_2h - pred_f
    row = []
    for bmask in band_masks_list:
        e = err[bmask & ~np.isnan(err)]
        row.append(np.abs(e).mean() if len(e) > 10 else np.nan)
    hm_poly.append(row)
hm = np.array(hm_poly)
im = ax7.imshow(hm, aspect='auto', cmap='RdYlGn_r', vmin=8, vmax=30)
ax7.set_xticks(range(len(band_names))); ax7.set_xticklabels(band_names, fontsize=8)
ax7.set_yticks(range(len(patients))); ax7.set_yticklabels([p['id'] for p in patients], fontsize=8)
for i in range(len(patients)):
    for j in range(len(band_names)):
        if not np.isnan(hm[i, j]):
            ax7.text(j, i, f'{hm[i,j]:.1f}', ha='center', va='center',
                     fontsize=8, color='black' if hm[i,j] < 20 else 'white', fontweight='bold')
cbar = plt.colorbar(im, ax=ax7, fraction=0.046, pad=0.04)
cbar.ax.tick_params(colors=TXT, labelsize=7)

ax8 = fig.add_subplot(gs[2, 1])
style(ax8, 'Phase 1 MAE Heatmap: Hybrid')
hm_hybr = []
for p in patients:
    r = all_results['hybrid'][p['id']]
    isf_true = isf_actual * p['scale']
    isf_cal = isf_hybrid(bg) * r['S_phase1']
    pred_f = bg - pred_drop * (isf_cal / isf_true)
    err = actual_bg_2h - pred_f
    row = []
    for bmask in band_masks_list:
        e = err[bmask & ~np.isnan(err)]
        row.append(np.abs(e).mean() if len(e) > 10 else np.nan)
    hm_hybr.append(row)
hm2 = np.array(hm_hybr)
im2 = ax8.imshow(hm2, aspect='auto', cmap='RdYlGn_r', vmin=8, vmax=30)
ax8.set_xticks(range(len(band_names))); ax8.set_xticklabels(band_names, fontsize=8)
ax8.set_yticks(range(len(patients))); ax8.set_yticklabels([p['id'] for p in patients], fontsize=8)
for i in range(len(patients)):
    for j in range(len(band_names)):
        if not np.isnan(hm2[i, j]):
            ax8.text(j, i, f'{hm2[i,j]:.1f}', ha='center', va='center',
                     fontsize=8, color='black' if hm2[i,j] < 20 else 'white', fontweight='bold')
cbar2 = plt.colorbar(im2, ax=ax8, fraction=0.046, pad=0.04)
cbar2.ax.tick_params(colors=TXT, labelsize=7)

# P9: Delta MAE (poly - hybrid) heatmap
ax9 = fig.add_subplot(gs[2, 2])
style(ax9, 'Δ MAE (Poly − Hybrid) by Glucose Band')
delta_hm = np.array(hm_poly) - np.array(hm_hybr)
im3 = ax9.imshow(delta_hm, aspect='auto', cmap='RdBu_r', vmin=-10, vmax=10)
ax9.set_xticks(range(len(band_names))); ax9.set_xticklabels(band_names, fontsize=8)
ax9.set_yticks(range(len(patients))); ax9.set_yticklabels([p['id'] for p in patients], fontsize=8)
for i in range(len(patients)):
    for j in range(len(band_names)):
        if not np.isnan(delta_hm[i, j]):
            ax9.text(j, i, f'{delta_hm[i,j]:+.1f}', ha='center', va='center',
                     fontsize=8, color='black', fontweight='bold')
cbar3 = plt.colorbar(im3, ax=ax9, fraction=0.046, pad=0.04)
cbar3.ax.tick_params(colors=TXT, labelsize=7)
cbar3.set_label('Δ MAE (mg/dL)', color=TXT, fontsize=8)


fig.suptitle(f'Full Polynomial vs Hybrid — Phase 1 Calibration Comparison\n'
             f'Jun 2025 – Mar 2026  |  {len(strict):,} overnight samples  |  '
             f'7 synthetic patients (TDD 11–45 U/day)',
             color=TXT, fontsize=13, fontweight='bold', y=0.995)

plt.savefig(OUT_DIR / 'ns_poly_vs_hybrid_calibration.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_poly_vs_hybrid_calibration.png")


# ══════════════════════════════════════════════════════════════════════════════
# 6. SUMMARY FILE
# ══════════════════════════════════════════════════════════════════════════════

summary = f"""FULL POLYNOMIAL vs HYBRID — CALIBRATION COMPARISON
{'='*65}

Dataset: Jun 2025 – Mar 2026 (10 months)
Valid overnight +2h samples: {len(strict):,}
Days with fasting data: {n_dates}
Real patient TDD (7D median): {real_tdd_median:.1f} U/day
Polynomial anchor (at target {TARGET}): {POLY_AT_TARGET:.1f}
Hybrid anchor (at 105): {HYBRID_AT_105}

HEAD-TO-HEAD: Phase 1 (TDD-based calibration, no in-cycle TDD)
"""

summary += f"\n  {'Patient':<6s}  {'TDD':>5s}  "
summary += f"{'Poly MAE':>8s}  {'Poly Bias':>9s}  "
summary += f"{'Hybr MAE':>8s}  {'Hybr Bias':>9s}  {'Δ MAE':>6s}\n"
summary += f"  {'─'*6}  {'─'*5}  {'─'*8}  {'─'*9}  {'─'*8}  {'─'*9}  {'─'*6}\n"

for p in patients:
    rp = all_results['poly'][p['id']]
    rh = all_results['hybrid'][p['id']]
    delta = rp['mae_p1'] - rh['mae_p1']
    summary += (f"  {p['id']:<6s}  {rp['tdd_synthetic']:5.1f}  "
                f"{rp['mae_p1']:8.1f}  {rp['bias_p1']:+9.1f}  "
                f"{rh['mae_p1']:8.1f}  {rh['bias_p1']:+9.1f}  {delta:+6.1f}\n")

with open(OUT_DIR / 'ns_poly_vs_hybrid_summary.txt', 'w') as f:
    f.write(summary)

print("\n" + summary)
print("Saved: ns_poly_vs_hybrid_summary.txt")
