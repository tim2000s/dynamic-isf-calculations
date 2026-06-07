#!/usr/bin/env python3
"""
Polynomial ISF Backtest — ADA Poster Equation + Hybrid
=======================================================

Tests the polynomial ISF equation from the ADA scientific poster:
  ISF(G) = 272 - 3.121×G + 0.01511×G² - 3.305e-05×G³ + 2.69e-08×G⁴

Plus a hybrid formula:
  BG >= 105: polynomial (above)
  BG <  105: 75.8 × (105/BG)^3.5  (power-law anchored at poly value at 105)

Compared against:
  A: Loop actual (variable_sens)
  B: Current ln with Boost blended TDD
  C: Power-law k=3.5 with Boost blended TDD (our proposed formula)
  D: Polynomial (ADA) — no TDD
  E: Polynomial (ADA) — scaled to patient's ISF at target
  F: Hybrid (poly ≥105, power-law <105) — no TDD
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
LN_TARGET = np.log(TARGET / D + 1)

HOME    = Path.home()
NS_WORK = HOME / 'Nightscout_Work'
OUT_DIR = HOME / 'Downloads'


# ── The polynomial ISF equation ──
def isf_polynomial(bg):
    """ADA poster equation: ISF = 272 - 3.121G + 0.01511G² - 3.305e-05G³ + 2.69e-08G⁴"""
    G = np.asarray(bg, dtype=float)
    return 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4


def isf_hybrid(bg):
    """Hybrid: polynomial ≥105, power-law <105 anchored at poly(105)=75.8"""
    G = np.asarray(bg, dtype=float)
    poly = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    power = 75.8 * (105.0 / G) ** 3.5
    return np.where(G >= 105, poly, power)


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
print("POLYNOMIAL ISF BACKTEST — ADA POSTER EQUATION")
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
        ib = r.get('openaps', {}).get('iob', {}) or {}
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
print(f"  After strict filtering: {len(strict):,} valid 2h samples")


# ══════════════════════════════════════════════════════════════════════════════
# FORMULA COMPARISON
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Computing ISF for each formula ──")

bg = strict['bg'].values
ln_bg = np.log(bg / D + 1)
isf_actual = strict['variable_sens'].values
pred_drop = strict['pred_drop'].values
actual_bg_2h = strict['actual_bg_2h'].values
tdd_b = strict['tdd_boost'].values

# A: Loop actual (variable_sens)
isf_A = isf_actual

# B: Current ln with Boost blended TDD
isf_B = 1800.0 / (tdd_b * ln_bg)

# C: Power-law k=3.5 with Boost blended TDD
isf_C = (1800.0 / tdd_b) * (TARGET / bg) ** 3.5

# D: Polynomial (ADA) — no TDD, raw equation
isf_D = isf_polynomial(bg)

# E: Polynomial scaled to patient ISF at target
# Scale factor: patient's ISF at target / equation's ISF at target
isf_at_target_eq = isf_polynomial(TARGET)  # ~81.6
patient_isf_at_target = np.median(1800.0 / tdd_b)  # ~107
scale_factor = patient_isf_at_target / isf_at_target_eq
isf_E = isf_polynomial(bg) * scale_factor

print(f"  Polynomial ISF at target ({TARGET}): {isf_at_target_eq:.1f}")
print(f"  Patient ISF at target (1800/TDD_boost median): {patient_isf_at_target:.1f}")
print(f"  Scale factor: {scale_factor:.3f}")


def compute_errors(isf_formula, label):
    pred_f = bg - pred_drop * (isf_formula / isf_actual)
    err = actual_bg_2h - pred_f
    valid = ~np.isnan(err)
    e = err[valid]
    mae = np.abs(e).mean()
    bias = e.mean()
    w18 = (np.abs(e) <= 18).mean() * 100
    return {'label': label, 'mae': mae, 'bias': bias, 'w18': w18, 'n': len(e), 'errors': e}


# F: Hybrid (poly ≥105, power-law <105) — no TDD
isf_F = isf_hybrid(bg)

# G: Full polynomial scaled by TDD-based S (no hybrid cutoff)
# S = (1800 / TDD) / 75.8  per-cycle, using Boost blended TDD
S_tdd = (1800.0 / tdd_b) / 75.8
isf_G = isf_polynomial(bg) * S_tdd
S_tdd_median = np.median(S_tdd)
print(f"  TDD-based S (median): {S_tdd_median:.3f}")

# H: Full polynomial scaled by fixed S from median TDD (simpler)
isf_H = isf_polynomial(bg) * S_tdd_median

formulas = [
    (isf_A, 'A: Loop actual (variable_sens)'),
    (isf_B, 'B: Current ln, Boost blended TDD'),
    (isf_C, 'C: Power-law k=3.5, Boost blended TDD'),
    (isf_D, 'D: Polynomial — raw (no scaling)'),
    (isf_E, f'E: Polynomial — scaled at target ×{scale_factor:.2f}'),
    (isf_F, 'F: Hybrid (poly≥105 + PL<105) — no TDD'),
    (isf_G, 'G: Full polynomial — TDD-scaled per cycle'),
    (isf_H, f'H: Full polynomial — TDD-scaled fixed ×{S_tdd_median:.2f}'),
]

print(f"\n{'═'*80}")
print(f"RESULTS — Overnight 00:00–07:00  (n={len(strict):,})")
print(f"{'═'*80}")
print(f"\n  {'Formula':<50s}  {'MAE':>5s}  {'Bias':>7s}  {'±1mmol':>6s}")
print(f"  {'─'*50}  {'─'*5}  {'─'*7}  {'─'*6}")

results = []
for isf_f, label in formulas:
    r = compute_errors(isf_f, label)
    results.append(r)
    print(f"  {label:<50s}  {r['mae']:5.1f}  {r['bias']:+7.1f}  {r['w18']:5.1f}%")


# ── BG-Band Analysis ──
print(f"\n── BG-Band Analysis ──")
bands = [('<90', bg < 90), ('90-105', (bg >= 90) & (bg < 105)),
         ('105-120', (bg >= 105) & (bg < 120)), ('120-150', (bg >= 120) & (bg <= 150)),
         ('150-200', (bg >= 150) & (bg <= 200))]

for isf_f, label in formulas:
    pred_f = bg - pred_drop * (isf_f / isf_actual)
    err = actual_bg_2h - pred_f
    print(f"\n  {label}:")
    for bname, bmask in bands:
        e = err[bmask & ~np.isnan(err)]
        if len(e) > 10:
            print(f"    {bname:>8s}  n={len(e):4d}  MAE={np.abs(e).mean():5.1f}  bias={e.mean():+6.1f}")


# ── ISF curve comparison ──
print(f"\n── ISF Curve Values ──")
print(f"  {'BG':>5s}  {'ln+Boost':>9s}  {'PL k=3.5':>9s}  {'Poly raw':>9s}  {'Poly TDD':>9s}  {'Hybrid':>9s}")
tdd_med = np.median(tdd_b)
for g in [60, 70, 80, 90, 99, 105, 110, 120, 140, 160, 180, 200]:
    isf_ln = 1800.0 / (tdd_med * np.log(g / D + 1))
    isf_pl = (1800.0 / tdd_med) * (TARGET / g) ** 3.5
    isf_poly = isf_polynomial(g)
    isf_poly_tdd = isf_poly * S_tdd_median
    isf_hyb = isf_hybrid(g)
    print(f"  {g:5d}  {isf_ln:9.1f}  {isf_pl:9.1f}  {isf_poly:9.1f}  {isf_poly_tdd:9.1f}  {isf_hyb:9.1f}")


# ── Hourly analysis for polynomial ──
print(f"\n── Hourly Bias: Polynomial (no TDD) ──")
for h in range(8):
    hm = strict['hour'] == h
    if hm.sum() > 20:
        bg_h = strict.loc[hm, 'bg'].values
        isf_act_h = strict.loc[hm, 'variable_sens'].values
        pd_h = strict.loc[hm, 'pred_drop'].values
        abg_h = strict.loc[hm, 'actual_bg_2h'].values

        isf_poly_h = isf_polynomial(bg_h)
        pred_poly_h = bg_h - pd_h * (isf_poly_h / isf_act_h)
        err_h = abg_h - pred_poly_h
        v = ~np.isnan(err_h)
        print(f"  {h:02d}:00  n={v.sum():4d}  MAE={np.abs(err_h[v]).mean():5.1f}  bias={err_h[v].mean():+6.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE
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

fig = plt.figure(figsize=(22, 16))
fig.patch.set_facecolor(BG_C)
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

colors = ['#4fc3f7', '#f48fb1', '#a5d6a7', '#ffb74d', '#ce93d8', '#ff6e6e', '#80cbc4', '#fff176']

# P1: ISF curves comparison
ax1 = fig.add_subplot(gs[0, :])
style(ax1, f'ISF vs BG — Formula Comparison (Boost blended TDD median = {tdd_med:.1f})')
bg_range = np.linspace(60, 200, 300)
ln_bg_r = np.log(bg_range / D + 1)

curve_ln = 1800.0 / (tdd_med * ln_bg_r)
curve_pl = (1800.0 / tdd_med) * (TARGET / bg_range) ** 3.5
curve_poly = isf_polynomial(bg_range)
curve_poly_s = curve_poly * scale_factor

curve_hybrid = isf_hybrid(bg_range)
curve_poly_tdd = isf_polynomial(bg_range) * S_tdd_median

ax1.plot(bg_range, curve_ln, lw=2.5, color=colors[1], label='Current ln + Boost TDD')
ax1.plot(bg_range, curve_pl, lw=2.5, color=colors[2], label='Power-law k=3.5 + Boost TDD')
ax1.plot(bg_range, curve_poly, lw=2.5, color=colors[3], label='Polynomial — raw')
ax1.plot(bg_range, curve_poly_s, lw=2.5, color=colors[4], ls='--', label=f'Polynomial — scaled at target ×{scale_factor:.2f}')
ax1.plot(bg_range, curve_hybrid, lw=2.5, color=colors[5], ls='-', label='Hybrid (poly≥105 + PL<105)')
ax1.plot(bg_range, curve_poly_tdd, lw=2.5, color=colors[6], ls='-', label=f'Full polynomial — TDD-scaled ×{S_tdd_median:.2f}')
ax1.axvline(TARGET, color='white', lw=0.8, ls=':', alpha=0.6, label=f'Target {TARGET:.0f}')
ax1.set_xlabel('BG (mg/dL)'); ax1.set_ylabel('ISF (mg/dL/U)')
ax1.set_xlim(60, 200); ax1.set_ylim(0, 200)
ax1.legend(fontsize=9, labelcolor=TXT, facecolor=PANEL, loc='upper right')

# P2: MAE bar chart
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'MAE by Formula (+2h)')
labels_short = ['Loop\nactual', 'ln\nBoost', 'PL\nk=3.5', 'Poly\nraw', 'Poly\ntarget', 'Hybrid', 'Poly\nTDD/cyc', 'Poly\nTDD/fix']
maes = [r['mae'] for r in results]
bars = ax2.bar(range(len(maes)), maes, color=colors, alpha=0.85)
for bar, v in zip(bars, maes):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2, f'{v:.1f}',
             ha='center', fontsize=7, color=TXT)
ax2.set_xticks(range(len(labels_short))); ax2.set_xticklabels(labels_short, fontsize=7)
ax2.set_ylabel('MAE (mg/dL)')

# P3: Bias bar chart
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'Prediction Bias (+2h)')
biases = [r['bias'] for r in results]
bars = ax3.bar(range(len(biases)), biases, color=colors, alpha=0.85)
for bar, v in zip(bars, biases):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + (0.3 if v >= 0 else -1.5),
             f'{v:+.1f}', ha='center', fontsize=7, color=TXT)
ax3.axhline(0, color='white', lw=0.8, ls='--')
ax3.set_xticks(range(len(labels_short))); ax3.set_xticklabels(labels_short, fontsize=7)
ax3.set_ylabel('Mean Error (mg/dL)')

# P4: ±1mmol bar chart
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, '±1 mmol/L Accuracy (+2h)')
w18s = [r['w18'] for r in results]
bars = ax4.bar(range(len(w18s)), w18s, color=colors, alpha=0.85)
for bar, v in zip(bars, w18s):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3, f'{v:.1f}%',
             ha='center', fontsize=7, color=TXT)
ax4.set_xticks(range(len(labels_short))); ax4.set_xticklabels(labels_short, fontsize=7)
ax4.set_ylabel('% within ±18 mg/dL')

# P5: Error distributions
ax5 = fig.add_subplot(gs[2, 0])
style(ax5, 'Prediction Error Distribution (+2h)')
for r, c in zip(results[1:], colors[1:]):
    ax5.hist(r['errors'].clip(-80, 80), bins=40, alpha=0.35, color=c,
             label=r['label'][:25], density=True)
ax5.axvline(0, color='white', lw=0.8, ls='--')
ax5.set_xlabel('Pred Error (mg/dL)'); ax5.set_ylabel('Density')
ax5.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

# P6: BG-band MAE for key formulas
ax6 = fig.add_subplot(gs[2, 1])
style(ax6, 'MAE by BG Band')
band_names = ['<90', '90-105', '105-120', '120-150']
band_masks = [bg < 90, (bg >= 90) & (bg < 105), (bg >= 105) & (bg < 120), (bg >= 120) & (bg <= 150)]

x = np.arange(len(band_names))
w = 0.18
for i, (isf_f, label) in enumerate(formulas[1:]):  # skip loop actual
    pred_f = bg - pred_drop * (isf_f / isf_actual)
    err = actual_bg_2h - pred_f
    band_maes = []
    for bmask in band_masks:
        e = err[bmask & ~np.isnan(err)]
        band_maes.append(np.abs(e).mean() if len(e) > 10 else np.nan)
    ax6.bar(x + i*w - 1.5*w, band_maes, w, color=colors[i+1], alpha=0.85,
            label=label[:20])
ax6.set_xticks(x); ax6.set_xticklabels(band_names, fontsize=8)
ax6.set_xlabel('BG Band (mg/dL)'); ax6.set_ylabel('MAE (mg/dL)')
ax6.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

# P7: BG-band bias for key formulas
ax7 = fig.add_subplot(gs[2, 2])
style(ax7, 'Bias by BG Band')
for i, (isf_f, label) in enumerate(formulas[1:]):
    pred_f = bg - pred_drop * (isf_f / isf_actual)
    err = actual_bg_2h - pred_f
    band_biases = []
    for bmask in band_masks:
        e = err[bmask & ~np.isnan(err)]
        band_biases.append(e.mean() if len(e) > 10 else np.nan)
    ax7.bar(x + i*w - 1.5*w, band_biases, w, color=colors[i+1], alpha=0.85,
            label=label[:20])
ax7.axhline(0, color='white', lw=0.8, ls='--')
ax7.set_xticks(x); ax7.set_xticklabels(band_names, fontsize=8)
ax7.set_xlabel('BG Band (mg/dL)'); ax7.set_ylabel('Bias (mg/dL)')
ax7.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

fig.suptitle(f'Polynomial ISF Backtest — ADA Poster Equation + Hybrid\n'
             f'Jun 2025 – Mar 2026  |  {len(strict):,} valid overnight samples',
             color=TXT, fontsize=12, fontweight='bold', y=0.995)

plt.savefig(OUT_DIR / 'ns_polynomial_backtest.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_polynomial_backtest.png")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
summary = f"""POLYNOMIAL ISF BACKTEST — SUMMARY
{'='*60}
Source: ADA Scientific Poster
  https://ada.scientificposters.com/epsAbstractADA.cfm?id=1
Equation (extended to all BG):
  ISF = 272 - 3.121×G + 0.01511×G² - 3.305e-05×G³ + 2.69e-08×G⁴

Dataset: Jun 2025 – Mar 2026 (10 months)
Valid overnight +2h samples: {len(strict):,}

Polynomial ISF at target ({TARGET:.0f}): {isf_at_target_eq:.1f} mg/dL/U
Patient ISF at target (1800/TDD_boost): {patient_isf_at_target:.1f} mg/dL/U
Scale factor: {scale_factor:.3f}

RESULTS (+2h overnight):
"""

for r in results:
    summary += f"  {r['label']:<50s}  MAE={r['mae']:.1f}  bias={r['bias']:+.1f}  ±1mmol={r['w18']:.1f}%\n"

with open(OUT_DIR / 'ns_polynomial_backtest_summary.txt', 'w') as f:
    f.write(summary)

print("\n" + summary)
print("Saved: ns_polynomial_backtest_summary.txt")
