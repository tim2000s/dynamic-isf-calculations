#!/usr/bin/env python3
"""
TDD Variant Analysis — Daytime vs Overnight Comparison
=======================================================

Same framework as ns_tdd_variants_analysis.py but compares:
  - Overnight (00:00–07:00, COB=0, 3h since bolus) — fasting
  - Daytime fasting (08:00–23:59, COB=0, 3h since bolus) — fasting but awake
  - All hours fasting (COB=0, 3h since bolus) — combined

This tells us whether the TDD recommendation changes by time of day.
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
print("TDD VARIANT ANALYSIS — DAYTIME vs OVERNIGHT")
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
        return tdd_weighted, tdd_last4h, tdd_8to4h

    if (tdd_weighted < (0.75 * tdd_7d)) and (tdd_1d is not None):
        tdd = ((tdd_weighted + ((tdd_weighted / tdd_7d) * (tdd_7d - tdd_weighted))) * 0.34) + (tdd_1d * 0.33) + (tdd_weighted * 0.33)
    elif tdd_1d is not None:
        tdd = (tdd_weighted * 0.33) + (tdd_7d * 0.34) + (tdd_1d * 0.33)
    else:
        tdd = tdd_weighted

    return tdd, tdd_last4h, tdd_8to4h


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
            'sensitivity_ratio': sg.get('sensitivityRatio'),
            'cob': sg.get('COB'), 'iob': sg.get('IOB'),
            'pred_iob_24': pred_iob[24] if len(pred_iob) > 24 else None,
        })
    except:
        continue

df = pd.DataFrame(cycles).dropna(subset=['ts', 'bg', 'variable_sens']).sort_values('ts').reset_index(drop=True)
print(f"  Parsed: {len(df):,} cycles")


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


# ── Attach TDD variants ──
print("\n── Attaching TDD variants ──")
daily_tdd['date'] = pd.to_datetime(daily_tdd['date']).dt.date
df = df.merge(daily_tdd[['date', 'tdd_actual', 'tdd_7day', 'tdd_1day']], on='date', how='left')

tdd_boost = np.full(len(df), np.nan)
tdd_last24h_arr = np.full(len(df), np.nan)

print("  Computing Boost blend + last24h per cycle...")
for i, (epoch, tdd7, tdd1) in enumerate(zip(cycle_epochs, df['tdd_7day'].values, df['tdd_1day'].values)):
    tdd7_v = tdd7 if not np.isnan(tdd7) else None
    tdd1_v = tdd1 if not np.isnan(tdd1) else None
    tdd, _, _ = compute_boost_tdd(epoch, tdd7_v, tdd1_v)
    tdd_boost[i] = tdd
    tdd_last24h_arr[i] = compute_insulin_window(epoch, 24, 0)
    if (i + 1) % 10000 == 0:
        print(f"  ... {i+1:,}/{len(df):,}")

df['tdd_boost'] = tdd_boost
df['tdd_24h_actual'] = tdd_last24h_arr
df['tdd_7d_disc25'] = df['tdd_7day'] * 0.75

ln_bg_all = np.log(df['bg'].values / D + 1)
df['tdd_implied'] = 1800.0 / (df['variable_sens'].values * ln_bg_all)


# ── Fasting filter (all hours) ──
print("\n── Fasting filter ──")
fasting_mask = (
    (df['cob'].fillna(99) == 0) &
    (df['bg'] >= 72) & (df['bg'] <= 200) &
    (df['bolus_age_min'] >= 180) &
    (df['tdd_boost'] > 0) &
    (~np.isnan(df['tdd_7day']))
)
fasting = df[fasting_mask].copy()
print(f"  All fasting cycles: {len(fasting):,}")

fasting = fasting.dropna(subset=['pred_iob_24', 'actual_bg_2h']).copy()
fasting['pred_drop'] = fasting['bg'] - fasting['pred_iob_24']
fasting['actual_drop'] = fasting['bg'] - fasting['actual_bg_2h']

strict = fasting[
    (np.abs(fasting['pred_drop']) > 3) &
    (fasting['actual_drop'] * fasting['pred_drop'] > 0) &
    ((fasting['actual_drop'] / fasting['pred_drop']) < 5) &
    ((fasting['actual_drop'] / fasting['pred_drop']) > 0) &
    (fasting['bg'] - fasting['actual_bg_2h'] < 9 + fasting['pred_drop'])
].copy()
print(f"  After strict filtering: {len(strict):,} valid 2h samples")

# Time-of-day splits
overnight = strict[strict['hour'] < 8].copy()
daytime = strict[(strict['hour'] >= 8) & (strict['hour'] < 22)].copy()
dawn = strict[(strict['hour'] >= 4) & (strict['hour'] < 8)].copy()
deep_night = strict[strict['hour'] < 4].copy()

print(f"\n  Overnight (00–08): {len(overnight):,}")
print(f"    Deep night (00–04): {len(deep_night):,}")
print(f"    Dawn (04–08): {len(dawn):,}")
print(f"  Daytime (08–22): {len(daytime):,}")

# Hourly breakdown
print(f"\n  Hourly sample counts:")
for h in range(24):
    n = (strict['hour'] == h).sum()
    if n > 0:
        print(f"    {h:02d}:00  n={n:5d}")


# ══════════════════════════════════════════════════════════════════════════════
# EVALUATE FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

def evaluate_tdd(subset, tdd_col, name, C_base=1800.0):
    """Evaluate a TDD variant on a given subset."""
    bg = subset['bg'].values
    ln_bg = np.log(bg / D + 1)
    isf_actual = subset['variable_sens'].values
    pred_drop = subset['pred_drop'].values
    actual_bg_2h = subset['actual_bg_2h'].values
    tdd = subset[tdd_col].values

    valid_tdd = ~np.isnan(tdd) & (tdd > 0)
    if valid_tdd.sum() < 50:
        return None

    # ln formula
    isf_ln = 1800.0 / (tdd * ln_bg)
    pred_ln = bg - pred_drop * (isf_ln / isf_actual)
    err_ln = actual_bg_2h - pred_ln
    v_ln = ~np.isnan(err_ln)
    mae_ln = np.abs(err_ln[v_ln]).mean()
    bias_ln = err_ln[v_ln].mean()

    # Power-law: optimal k (MAE)
    def mae_pl(k):
        isf_f = (C_base / tdd) * (TARGET / bg) ** k
        pred_f = bg - pred_drop * (isf_f / isf_actual)
        err = actual_bg_2h - pred_f
        v = ~np.isnan(err)
        return np.abs(err[v]).mean()

    res_k = minimize_scalar(mae_pl, bounds=(0.5, 4.0), method='bounded')
    k_opt = res_k.x
    mae_opt = res_k.fun

    # Metrics at optimal k
    isf_opt = (C_base / tdd) * (TARGET / bg) ** k_opt
    pred_opt = bg - pred_drop * (isf_opt / isf_actual)
    err_opt = actual_bg_2h - pred_opt
    v_opt = ~np.isnan(err_opt)
    bias_opt = err_opt[v_opt].mean()
    w18 = (np.abs(err_opt[v_opt]) <= 18).mean() * 100

    # Joint C + k
    def mae_ck(params):
        C, k = params
        isf_f = (C / tdd) * (TARGET / bg) ** k
        pred_f = bg - pred_drop * (isf_f / isf_actual)
        err = actual_bg_2h - pred_f
        return np.nanmean(np.abs(err))

    res_ck = minimize(mae_ck, [1800, 2.0], bounds=[(500, 5000), (0.5, 4.0)], method='L-BFGS-B')

    return {
        'name': name,
        'tdd_median': np.nanmedian(tdd),
        'n': int(v_ln.sum()),
        'mae_ln': mae_ln, 'bias_ln': bias_ln,
        'k_opt': k_opt, 'mae_opt': mae_opt,
        'bias_opt': bias_opt, 'w18': w18,
        'joint_C': res_ck.x[0], 'joint_k': res_ck.x[1],
        'joint_mae': res_ck.fun,
    }


# ══════════════════════════════════════════════════════════════════════════════
# RUN COMPARISON ACROSS TIME WINDOWS
# ══════════════════════════════════════════════════════════════════════════════

tdd_cols = [
    ('tdd_7day',      'Actual 7-day'),
    ('tdd_boost',     'Boost blended'),
    ('tdd_7d_disc25', '7-day × 0.75'),
    ('tdd_24h_actual','Last 24h actual'),
    ('tdd_implied',   'Implied (ref)'),
]

time_windows = [
    ('Overnight 00–08', overnight),
    ('Deep night 00–04', deep_night),
    ('Dawn 04–08', dawn),
    ('Daytime 08–22', daytime),
    ('All fasting', strict),
]

print("\n\n" + "=" * 100)
print("RESULTS BY TIME WINDOW AND TDD VARIANT")
print("=" * 100)

all_window_results = {}

for window_name, subset in time_windows:
    print(f"\n{'─'*80}")
    print(f"  {window_name}  (n={len(subset):,})")
    print(f"{'─'*80}")
    print(f"  {'TDD Variant':<22s}  {'med':>5s}  {'ln MAE':>6s}  {'ln bias':>7s}  {'PL k':>5s}  {'PL MAE':>6s}  {'PL bias':>7s}  {'±1mm':>5s}  {'C*':>5s}  {'k*':>5s}  {'MAE*':>5s}")

    window_results = []
    for col, name in tdd_cols:
        r = evaluate_tdd(subset, col, name)
        if r:
            window_results.append(r)
            print(f"  {name:<22s}  {r['tdd_median']:5.1f}  {r['mae_ln']:6.2f}  {r['bias_ln']:+7.2f}  {r['k_opt']:5.2f}  {r['mae_opt']:6.2f}  {r['bias_opt']:+7.2f}  {r['w18']:5.1f}  {r['joint_C']:5.0f}  {r['joint_k']:5.2f}  {r['joint_mae']:5.2f}")

    all_window_results[window_name] = window_results

# ══════════════════════════════════════════════════════════════════════════════
# BOOST BLENDED TDD DISTRIBUTION BY HOUR
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'─'*80}")
print("BOOST BLENDED TDD BY HOUR OF DAY")
print(f"{'─'*80}")
print(f"  {'Hour':>4s}  {'n':>5s}  {'Boost med':>9s}  {'7-day med':>9s}  {'Ratio':>6s}  {'Implied med':>11s}")

for h in range(24):
    hm = strict['hour'] == h
    if hm.sum() > 10:
        boost_med = strict.loc[hm, 'tdd_boost'].median()
        tdd7_med = strict.loc[hm, 'tdd_7day'].median()
        impl_med = strict.loc[hm, 'tdd_implied'].median()
        ratio = boost_med / tdd7_med if tdd7_med > 0 else np.nan
        print(f"  {h:4d}  {hm.sum():5d}  {boost_med:9.1f}  {tdd7_med:9.1f}  {ratio:6.2f}  {impl_med:11.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# SENSITIVITY RATIO DISTRIBUTION BY HOUR
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'─'*80}")
print("SENSITIVITY RATIO BY HOUR (from loop data)")
print(f"{'─'*80}")
print(f"  {'Hour':>4s}  {'n':>5s}  {'median':>7s}  {'mean':>7s}  {'% at 1.0':>8s}")

for h in range(24):
    hm = strict['hour'] == h
    sr = strict.loc[hm, 'sensitivity_ratio'].dropna()
    if len(sr) > 10:
        at_one = (np.abs(sr - 1.0) < 0.01).mean() * 100
        print(f"  {h:4d}  {len(sr):5d}  {sr.median():7.2f}  {sr.mean():7.2f}  {at_one:7.1f}%")


# ══════════════════════════════════════════════════════════════════════════════
# KEY COMPARISON: same k across day and night?
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n\n{'='*80}")
print("KEY QUESTION: Does optimal k differ between day and night?")
print(f"{'='*80}")

for col, name in tdd_cols:
    r_night = evaluate_tdd(overnight, col, name)
    r_day = evaluate_tdd(daytime, col, name)
    r_all = evaluate_tdd(strict, col, name)
    if r_night and r_day and r_all:
        print(f"\n  {name}:")
        print(f"    Overnight: k={r_night['k_opt']:.2f}  C={r_night['joint_C']:.0f}  k*={r_night['joint_k']:.2f}  MAE={r_night['mae_opt']:.2f}")
        print(f"    Daytime:   k={r_day['k_opt']:.2f}  C={r_day['joint_C']:.0f}  k*={r_day['joint_k']:.2f}  MAE={r_day['mae_opt']:.2f}")
        print(f"    All:       k={r_all['k_opt']:.2f}  C={r_all['joint_C']:.0f}  k*={r_all['joint_k']:.2f}  MAE={r_all['mae_opt']:.2f}")
        k_diff = abs(r_night['k_opt'] - r_day['k_opt'])
        print(f"    k difference (night vs day): {k_diff:.2f}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE
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
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.5, wspace=0.38)

colors_tdd = ['#4fc3f7', '#f48fb1', '#a5d6a7', '#ffb74d', '#ce93d8']
colors_win = ['#4fc3f7', '#80deea', '#ffe082', '#ef9a9a', '#b39ddb']

# P1: Boost blended TDD by hour
ax1 = fig.add_subplot(gs[0, 0])
style(ax1, 'Boost Blended TDD by Hour')
hours_plot = []
boost_meds = []
tdd7_meds = []
for h in range(24):
    hm = strict['hour'] == h
    if hm.sum() > 10:
        hours_plot.append(h)
        boost_meds.append(strict.loc[hm, 'tdd_boost'].median())
        tdd7_meds.append(strict.loc[hm, 'tdd_7day'].median())
ax1.plot(hours_plot, boost_meds, 'o-', color=colors_tdd[1], lw=2, ms=5, label='Boost blended')
ax1.plot(hours_plot, tdd7_meds, 's-', color=colors_tdd[0], lw=2, ms=5, label='7-day actual')
ax1.set_xlabel('Hour'); ax1.set_ylabel('TDD (U/day)')
ax1.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P2: Sensitivity ratio by hour
ax2 = fig.add_subplot(gs[0, 1])
style(ax2, 'Sensitivity Ratio by Hour')
sr_meds = []
for h in hours_plot:
    sr = strict.loc[strict['hour'] == h, 'sensitivity_ratio'].dropna()
    sr_meds.append(sr.median() if len(sr) > 0 else np.nan)
ax2.plot(hours_plot, sr_meds, 'o-', color='#ffb74d', lw=2, ms=5)
ax2.axhline(1.0, color='white', lw=0.8, ls='--')
ax2.set_xlabel('Hour'); ax2.set_ylabel('sensitivityRatio (median)')

# P3: Implied TDD by hour
ax3 = fig.add_subplot(gs[0, 2])
style(ax3, 'Implied TDD by Hour')
impl_meds = []
for h in hours_plot:
    impl = strict.loc[strict['hour'] == h, 'tdd_implied'].dropna()
    impl_meds.append(impl.median() if len(impl) > 0 else np.nan)
ax3.plot(hours_plot, impl_meds, 'o-', color=colors_tdd[4], lw=2, ms=5, label='Implied')
ax3.plot(hours_plot, boost_meds, 's-', color=colors_tdd[1], lw=2, ms=5, label='Boost blended')
ax3.set_xlabel('Hour'); ax3.set_ylabel('TDD (U/day)')
ax3.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P4-P5: MAE comparison across time windows for each TDD variant
ax4 = fig.add_subplot(gs[1, 0])
style(ax4, 'Power-Law MAE at Optimal k by Window')
x = np.arange(len(tdd_cols))
w = 0.15
for i, (wname, _) in enumerate(time_windows[:4]):
    wresults = all_window_results.get(wname, [])
    if wresults:
        maes = [r['mae_opt'] for r in wresults]
        ax4.bar(x + i*w - 1.5*w, maes, w, color=colors_win[i], alpha=0.85, label=wname)
ax4.set_xticks(x)
ax4.set_xticklabels([n for _, n in tdd_cols], fontsize=7, rotation=25, ha='right')
ax4.set_ylabel('MAE (mg/dL)')
ax4.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

# P5: Optimal k by window
ax5 = fig.add_subplot(gs[1, 1])
style(ax5, 'Optimal k by Window')
for i, (wname, _) in enumerate(time_windows[:4]):
    wresults = all_window_results.get(wname, [])
    if wresults:
        ks = [r['k_opt'] for r in wresults]
        ax5.bar(x + i*w - 1.5*w, ks, w, color=colors_win[i], alpha=0.85, label=wname)
ax5.set_xticks(x)
ax5.set_xticklabels([n for _, n in tdd_cols], fontsize=7, rotation=25, ha='right')
ax5.set_ylabel('Optimal k')
ax5.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

# P6: Joint C by window
ax6 = fig.add_subplot(gs[1, 2])
style(ax6, 'Joint Optimal C by Window')
for i, (wname, _) in enumerate(time_windows[:4]):
    wresults = all_window_results.get(wname, [])
    if wresults:
        Cs = [r['joint_C'] for r in wresults]
        ax6.bar(x + i*w - 1.5*w, Cs, w, color=colors_win[i], alpha=0.85, label=wname)
ax6.set_xticks(x)
ax6.set_xticklabels([n for _, n in tdd_cols], fontsize=7, rotation=25, ha='right')
ax6.set_ylabel('Optimal C')
ax6.legend(fontsize=6, labelcolor=TXT, facecolor=PANEL)

# P7: Bias by window for Boost blended
ax7 = fig.add_subplot(gs[2, 0])
style(ax7, 'Bias by Window (Boost Blended TDD, PL at k_opt)')
win_names = []
win_biases_ln = []
win_biases_pl = []
for wname, _ in time_windows:
    wresults = all_window_results.get(wname, [])
    boost_r = next((r for r in wresults if 'Boost' in r['name']), None)
    if boost_r:
        win_names.append(wname)
        win_biases_ln.append(boost_r['bias_ln'])
        win_biases_pl.append(boost_r['bias_opt'])
xw = np.arange(len(win_names))
ax7.bar(xw - 0.15, win_biases_ln, 0.3, color='#ef5350', alpha=0.8, label='Current ln')
ax7.bar(xw + 0.15, win_biases_pl, 0.3, color='#66bb6a', alpha=0.8, label='Power-law')
ax7.axhline(0, color='white', lw=0.8, ls='--')
ax7.set_xticks(xw)
ax7.set_xticklabels(win_names, fontsize=7, rotation=25, ha='right')
ax7.set_ylabel('Bias (mg/dL)')
ax7.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P8: Sample count by hour
ax8 = fig.add_subplot(gs[2, 1])
style(ax8, 'Fasting Sample Count by Hour')
hour_counts = [strict[strict['hour'] == h].shape[0] for h in range(24)]
ax8.bar(range(24), hour_counts, color='#4fc3f7', alpha=0.8)
ax8.set_xlabel('Hour'); ax8.set_ylabel('n (fasting samples)')

# P9: Prediction error by hour (Boost blended, PL at k=3.5)
ax9 = fig.add_subplot(gs[2, 2])
style(ax9, 'Hourly Bias: Boost Blended PL k=3.5')
hourly_biases = []
for h in range(24):
    hm = strict['hour'] == h
    if hm.sum() > 20:
        bg_h = strict.loc[hm, 'bg'].values
        tdd_h = strict.loc[hm, 'tdd_boost'].values
        isf_act_h = strict.loc[hm, 'variable_sens'].values
        pd_h = strict.loc[hm, 'pred_drop'].values
        abg_h = strict.loc[hm, 'actual_bg_2h'].values

        isf_f = (1800.0 / tdd_h) * (TARGET / bg_h) ** 3.5
        pred_f = bg_h - pd_h * (isf_f / isf_act_h)
        err = abg_h - pred_f
        v = ~np.isnan(err)
        hourly_biases.append((h, err[v].mean(), len(err[v])))
    else:
        hourly_biases.append((h, np.nan, 0))

hours_b = [x[0] for x in hourly_biases if not np.isnan(x[1])]
biases_b = [x[1] for x in hourly_biases if not np.isnan(x[1])]
ax9.bar(hours_b, biases_b, color='#a5d6a7', alpha=0.8)
ax9.axhline(0, color='white', lw=0.8, ls='--')
ax9.set_xlabel('Hour'); ax9.set_ylabel('Bias (mg/dL)')

fig.suptitle(f'TDD Variant Analysis — Daytime vs Overnight\n'
             f'Jun 2025 – Mar 2026  |  {len(strict):,} fasting samples  |  Overnight: {len(overnight):,}  Daytime: {len(daytime):,}',
             color=TXT, fontsize=12, fontweight='bold', y=0.995)

plt.savefig(OUT_DIR / 'ns_tdd_daytime_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_tdd_daytime_results.png")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

for wname, _ in time_windows:
    wresults = all_window_results.get(wname, [])
    if not wresults: continue
    print(f"\n  {wname} (n={wresults[0]['n']:,}):")
    best_real = min([r for r in wresults if 'Implied' not in r['name']], key=lambda r: r['mae_opt'])
    best_joint = min([r for r in wresults if 'Implied' not in r['name']], key=lambda r: r['joint_mae'])
    print(f"    Best MAE (C=1800): {best_real['name']} — k={best_real['k_opt']:.2f}, MAE={best_real['mae_opt']:.2f}, bias={best_real['bias_opt']:+.2f}")
    print(f"    Best MAE (C,k free): {best_joint['name']} — C={best_joint['joint_C']:.0f}, k={best_joint['joint_k']:.2f}, MAE={best_joint['joint_mae']:.2f}")
