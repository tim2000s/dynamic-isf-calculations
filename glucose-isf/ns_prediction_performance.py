#!/usr/bin/env python3
"""
Prediction performance vs outcomes charts.

Uses cached per-sample data from ns_sigmoid_log_plots.py run.
Saves site data to JSON cache to avoid re-fetching.

Charts:
  1. Calibration plots (predicted vs actual, binned) — sigmoid & log
  2. Scatter: predicted vs actual with identity line
  3. Rolling error vs BG — how prediction bias shifts across glucose range
  4. Box plots of error distribution by BG band and formula
"""

import json, time, re, warnings, pickle
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

warnings.filterwarnings('ignore')

OUT_DIR = Path.home() / 'Downloads'
CACHE_FILE = OUT_DIR / 'multisite_sample_cache.pkl'

# ── ISF formulas ──
def isf_quartic(bg):
    G = np.asarray(bg, dtype=float)
    return 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4

def isf_full_diabeloop(bg):
    G = np.asarray(bg, dtype=float)
    quad = 98.03 - 1.077*G + 0.008868*G**2
    quart = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    return np.where(G <= 100, quad, quart)

def isf_hybrid(bg):
    G = np.asarray(bg, dtype=float)
    quart = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    power = 75.8 * (105.0 / G) ** 3.5
    return np.where(G >= 105, quart, power)

QUARTIC_AT_99 = float(isf_quartic(99.0))
FULL_DB_AT_99 = float(isf_full_diabeloop(99.0))
HYBRID_AT_105 = 75.8

SITES = [
    {'name': 'henny425',    'model': 'sigmoid', 'url': '***REDACTED-URL***',               'token': '***REDACTED-TOKEN***'},
    {'name': 'aadiabetes',  'model': 'sigmoid', 'url': '***REDACTED-URL***',     'token': '***REDACTED-TOKEN***'},
    {'name': 'diajesse',    'model': 'sigmoid', 'url': '***REDACTED-URL***',       'token': None},
    {'name': 'svns',        'model': 'sigmoid', 'url': '***REDACTED-URL***',  'token': None},
    {'name': 'andycgm',     'model': 'log',     'url': '***REDACTED-URL***',     'token': '***REDACTED-TOKEN***'},
    {'name': 'noahr',       'model': 'log',     'url': '***REDACTED-URL***',             'token': None},
    {'name': 'fuxchr',      'model': 'sigmoid', 'url': '***REDACTED-URL***',               'token': '***REDACTED-TOKEN***'},
    {'name': 'nightscout1', 'model': 'log',     'url': '***REDACTED-URL***',     'token': '***REDACTED-TOKEN***'},
    {'name': 'eli',         'model': 'log',     'url': '***REDACTED-URL***',            'token': '***REDACTED-TOKEN***'},
    {'name': 'mikens',      'model': 'sigmoid', 'url': '***REDACTED-URL***',          'token': '***REDACTED-TOKEN***'},
    {'name': 'ns_rot6',     'model': 'log',     'url': '***REDACTED-URL***',                      'token': '***REDACTED-TOKEN***'},
    {'name': 'kelseyhuss',  'model': 'log',     'url': '***REDACTED-URL***',  'token': '***REDACTED-TOKEN***'},
]

# ── API fetching (same as before) ──
def ns_fetch(base_url, endpoint, token=None, params=None, max_retries=3):
    if params is None: params = {}
    if token: params['token'] = token
    url = f"{base_url}/api/v1/{endpoint}?{urlencode(params, safe='[]$')}"
    for attempt in range(max_retries):
        try:
            req = Request(url, headers={'Accept': 'application/json',
                                         'User-Agent': 'Mozilla/5.0 (Nightscout-Backtest)'})
            with urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except HTTPError as e:
            if e.code in (401, 403): return None
            if attempt < max_retries - 1: time.sleep(2 ** attempt)
            else: return None
        except (URLError, TimeoutError):
            if attempt < max_retries - 1: time.sleep(2 ** attempt)
            else: return None

def fetch_all_paginated(base_url, endpoint, token=None, months_back=6, date_field='created_at'):
    all_records, seen_ids, page_size = [], set(), 10000
    now = datetime.utcnow()
    windows = []
    for m in range(months_back):
        end = now - timedelta(days=30 * m)
        start = now - timedelta(days=30 * (m + 1))
        windows.append((start.strftime('%Y-%m-%dT%H:%M:%S.000Z'),
                         end.strftime('%Y-%m-%dT%H:%M:%S.000Z')))
    for win_start, win_end in windows:
        oldest = win_end
        while True:
            params = {'count': page_size,
                      f'find[{date_field}][$gte]': win_start,
                      f'find[{date_field}][$lt]': oldest}
            batch = ns_fetch(base_url, endpoint, token, params)
            if not batch or len(batch) == 0: break
            new = 0
            for r in batch:
                rid = r.get('_id', r.get('identifier', id(r)))
                if rid not in seen_ids:
                    seen_ids.add(rid)
                    all_records.append(r)
                    new += 1
            if new == 0: break
            dates = [r.get(date_field) for r in batch if r.get(date_field)]
            if dates: oldest = min(dates)
            else: break
            if len(batch) < page_size: break
    return all_records

def parse_isf_from_reason(reason, bg_mgdl):
    m = re.search(r'ISF:\s*([\d.]+)\s*(?:→|->)\s*([\d.]+)', reason)
    if not m: return None
    dynamic_isf = float(m.group(2))
    target_m = re.search(r'Target:\s*([\d.]+)', reason)
    if target_m and float(target_m.group(1)) < 30:
        dynamic_isf *= 18.0
    return dynamic_isf

def parse_tdd_from_reason(reason):
    m = re.search(r'TDD:\s*([\d.]+)', reason)
    return float(m.group(1)) if m else None

def process_site(site):
    name, base_url, token, model = site['name'], site['url'], site['token'], site['model']
    print(f"  Fetching {name}...", end='', flush=True)
    test = ns_fetch(base_url, 'devicestatus.json', token, {'count': 1})
    if test is None:
        print(" SKIP")
        return None
    raw_ds = fetch_all_paginated(base_url, 'devicestatus.json', token, date_field='created_at')
    if not raw_ds or len(raw_ds) < 100:
        print(f" SKIP")
        return None
    raw_entries = fetch_all_paginated(base_url, 'entries.json', token, date_field='dateString')
    if not raw_entries or len(raw_entries) < 100:
        print(f" SKIP")
        return None
    raw_tx = fetch_all_paginated(base_url, 'treatments.json', token, date_field='created_at')

    cycles = []
    for r in raw_ds:
        try:
            ts = pd.to_datetime(r['created_at'], utc=True)
            sg = r.get('openaps', {}).get('suggested', {})
            if not sg or 'bg' not in sg: continue
            bg_val = sg['bg']
            isf_val = sg.get('ISF')
            reason = sg.get('reason', '')
            if not isf_val or isf_val <= 0:
                isf_val = parse_isf_from_reason(reason, bg_val)
            if not isf_val or isf_val <= 0: continue
            pred_iob = sg.get('predBGs', {}).get('IOB') or []
            tdd_val = sg.get('TDD')
            if tdd_val is None or tdd_val <= 0:
                tdd_val = parse_tdd_from_reason(reason)
            cycles.append({
                'ts': ts, 'bg': bg_val, 'isf_actual': float(isf_val),
                'cob': sg.get('COB'), 'iob': sg.get('IOB'),
                'tdd': float(tdd_val) if tdd_val else None,
                'pred_iob_24': pred_iob[24] if len(pred_iob) > 24 else None,
            })
        except Exception:
            continue
    if len(cycles) < 100:
        print(f" SKIP")
        return None
    df = pd.DataFrame(cycles).sort_values('ts').reset_index(drop=True)

    cgm_rows = []
    for e in raw_entries:
        try:
            sgv = e.get('sgv')
            if not sgv or not (40 <= sgv <= 400): continue
            ca = e.get('created_at') or e.get('dateString')
            ts = pd.to_datetime(ca, utc=True) if ca else pd.to_datetime(int(e['date']), unit='ms', utc=True)
            cgm_rows.append({'ts': ts, 'sgv': float(sgv)})
        except Exception:
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

    bolus_epochs = []
    for t in (raw_tx or []):
        try:
            ins = t.get('insulin')
            if ins and float(ins) > 0:
                ts = pd.to_datetime(t['created_at'], utc=True)
                bolus_epochs.append(int(ts.timestamp()))
        except Exception:
            continue
    bolus_epochs = np.array(sorted(bolus_epochs)) if bolus_epochs else np.array([])

    def mins_since_bolus(target_ts_s):
        if len(bolus_epochs) == 0: return 9999.0
        idx = np.searchsorted(bolus_epochs, target_ts_s, side='right') - 1
        return (target_ts_s - bolus_epochs[idx]) / 60.0 if idx >= 0 else 9999.0

    cycle_epochs = np.array([int(t.timestamp()) for t in df['ts']])
    actual_2h = np.full(len(df), np.nan)
    bolus_age = np.full(len(df), np.nan)
    for i, t_s in enumerate(cycle_epochs):
        actual_2h[i] = get_cgm_at(t_s + 7200)
        bolus_age[i] = mins_since_bolus(t_s)
    df['actual_bg_2h'] = actual_2h
    df['bolus_age_min'] = bolus_age
    df['hour'] = df['ts'].dt.hour
    df['date'] = df['ts'].dt.date
    tdd_valid = df[df['tdd'].notna() & (df['tdd'] > 0)]
    if len(tdd_valid) > 0:
        daily_tdd = tdd_valid.groupby('date')['tdd'].median().reset_index()
        daily_tdd.columns = ['date', 'tdd_daily']
        daily_tdd['tdd_7day'] = daily_tdd['tdd_daily'].rolling(7, min_periods=1).mean()
        df = df.merge(daily_tdd[['date', 'tdd_7day']], on='date', how='left')
    else:
        df['tdd_7day'] = df['tdd']

    mask = (
        (df['hour'] < 8) & (df['cob'].fillna(99) == 0) &
        (df['bg'] >= 72) & (df['bg'] <= 200) &
        (df['bolus_age_min'] >= 120) & (df['tdd'].notna()) & (df['tdd'] > 0)
    )
    on = df[mask].copy()
    on = on.dropna(subset=['pred_iob_24', 'actual_bg_2h']).copy()
    if len(on) < 10:
        print(f" SKIP")
        return None
    on['pred_drop'] = on['bg'] - on['pred_iob_24']
    on['actual_drop'] = on['bg'] - on['actual_bg_2h']
    strict = on[
        (np.abs(on['pred_drop']) > 3) &
        ((on['actual_drop'] / on['pred_drop']).between(0, 5))
    ].copy().sort_values('ts').reset_index(drop=True)
    if len(strict) < 10:
        print(f" SKIP")
        return None

    bg = strict['bg'].values
    isf_actual = strict['isf_actual'].values
    pred_drop = strict['pred_drop'].values
    actual_bg_2h = strict['actual_bg_2h'].values
    tdd_median = np.median(strict['tdd_7day'].values)

    S_q = (1800.0 / tdd_median) / QUARTIC_AT_99
    S_db = (1800.0 / tdd_median) / FULL_DB_AT_99
    S_h = (1800.0 / tdd_median) / HYBRID_AT_105

    isf_q = isf_quartic(bg) * S_q
    isf_db = isf_full_diabeloop(bg) * S_db
    isf_h = isf_hybrid(bg) * S_h

    pred_loop = bg - pred_drop * 1.0  # isf_actual / isf_actual = 1
    pred_q = bg - pred_drop * (isf_q / isf_actual)
    pred_db = bg - pred_drop * (isf_db / isf_actual)
    pred_h = bg - pred_drop * (isf_h / isf_actual)

    print(f" OK ({len(strict)} samples)")
    return {
        'name': name, 'model': model, 'n': len(strict),
        'bg': bg, 'actual_bg_2h': actual_bg_2h, 'pred_drop': pred_drop,
        'pred_loop': pred_loop, 'pred_q': pred_q, 'pred_db': pred_db, 'pred_h': pred_h,
    }


# ══════════════════════════════════════════════════════════════════════════════
# LOAD OR FETCH DATA
# ══════════════════════════════════════════════════════════════════════════════

if CACHE_FILE.exists():
    print(f"Loading cached data from {CACHE_FILE}")
    with open(CACHE_FILE, 'rb') as f:
        all_sites = pickle.load(f)
    print(f"Loaded {len(all_sites)} sites, {sum(s['n'] for s in all_sites)} samples")
else:
    print("Fetching data from all sites...")
    all_sites = []
    for site in SITES:
        r = process_site(site)
        if r is not None:
            all_sites.append(r)
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(all_sites, f)
    print(f"Cached to {CACHE_FILE}")

sigmoid_sites = [s for s in all_sites if s['model'] == 'sigmoid']
log_sites = [s for s in all_sites if s['model'] == 'log']

# Pool data per group
def pool(sites):
    bg = np.concatenate([s['bg'] for s in sites])
    act = np.concatenate([s['actual_bg_2h'] for s in sites])
    pd_ = np.concatenate([s['pred_drop'] for s in sites])
    pl = np.concatenate([s['pred_loop'] for s in sites])
    pq = np.concatenate([s['pred_q'] for s in sites])
    pdb = np.concatenate([s['pred_db'] for s in sites])
    ph = np.concatenate([s['pred_h'] for s in sites])
    return bg, act, pd_, pl, pq, pdb, ph

FORMULAS = [
    ('Loop (DynISF)', '#2E86AB'),
    ('Quartic',       '#E74C3C'),
    ('Full Diabeloop','#27AE60'),
    ('Hybrid',        '#F39C12'),
]


# ══════════════════════════════════════════════════════════════════════════════
# CHART 1: CALIBRATION PLOTS
# Each formula: bin predicted 2h BG, plot mean predicted vs mean actual
# Perfect calibration = points on the diagonal
# ══════════════════════════════════════════════════════════════════════════════

def plot_calibration(ax, actual, preds_dict, title):
    """Reliability diagram: binned mean predicted vs mean actual."""
    bins = np.arange(50, 220, 10)
    ax.plot([40, 250], [40, 250], 'k--', linewidth=1, alpha=0.5, label='Perfect')

    for (fname, color), pred in zip(FORMULAS, preds_dict.values()):
        bin_pred_means, bin_act_means, bin_ns = [], [], []
        for i in range(len(bins) - 1):
            mask = (pred >= bins[i]) & (pred < bins[i+1])
            if mask.sum() >= 5:
                bin_pred_means.append(pred[mask].mean())
                bin_act_means.append(actual[mask].mean())
                bin_ns.append(mask.sum())
        if bin_pred_means:
            ax.plot(bin_pred_means, bin_act_means, 'o-', color=color,
                    markersize=4, linewidth=1.5, label=fname, alpha=0.8)

    ax.set_xlabel('Mean Predicted 2h BG [mg/dL]', fontsize=10)
    ax.set_ylabel('Mean Actual 2h BG [mg/dL]', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_xlim(40, 220)
    ax.set_ylim(40, 220)
    ax.set_aspect('equal')
    ax.grid(alpha=0.3)


print("\nGenerating calibration plots...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))

for ax, (group_name, sites) in zip(axes, [('Sigmoid', sigmoid_sites), ('Log', log_sites)]):
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    n = len(bg)
    plot_calibration(ax, act, preds, f'{group_name} — Calibration\n({len(sites)} sites, {n:,} samples)')

fig.tight_layout()
fig.savefig(OUT_DIR / 'multisite_calibration.png', dpi=150, bbox_inches='tight')
print(f"Saved: multisite_calibration.png")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 2: SCATTER — predicted vs actual, coloured by starting BG zone
# ══════════════════════════════════════════════════════════════════════════════

def plot_scatter(axes_row, actual, preds_dict, bg, title_prefix):
    """Scatter of predicted vs actual for each formula, coloured by BG zone."""
    zone_colors = {
        '<90': '#E74C3C',
        '90-105': '#F39C12',
        '105-150': '#27AE60',
        '150+': '#2E86AB',
    }
    zone_masks = {
        '<90':     bg < 90,
        '90-105':  (bg >= 90) & (bg < 105),
        '105-150': (bg >= 105) & (bg < 150),
        '150+':    bg >= 150,
    }

    for ax, ((fname, fcolor), pred) in zip(axes_row, zip(FORMULAS, preds_dict.values())):
        ax.plot([40, 260], [40, 260], 'k--', linewidth=0.8, alpha=0.4)
        for zname, zmask in zone_masks.items():
            if zmask.sum() > 0:
                ax.scatter(pred[zmask], actual[zmask], s=3, alpha=0.15,
                          color=zone_colors[zname], label=zname, rasterized=True)
        ax.set_xlabel('Predicted 2h BG', fontsize=9)
        ax.set_ylabel('Actual 2h BG', fontsize=9)
        ax.set_title(f'{title_prefix} — {fname}', fontsize=10, fontweight='bold')
        ax.set_xlim(40, 260)
        ax.set_ylim(40, 260)
        ax.set_aspect('equal')
        ax.grid(alpha=0.2)
        if ax == axes_row[0]:
            ax.legend(fontsize=7, markerscale=3, title='Starting BG', title_fontsize=7)


print("Generating scatter plots...")
fig, axes = plt.subplots(2, 4, figsize=(20, 10))

for row, (group_name, sites) in enumerate([('Sigmoid', sigmoid_sites), ('Log', log_sites)]):
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    plot_scatter(axes[row], act, preds, bg, group_name)

fig.tight_layout()
fig.savefig(OUT_DIR / 'multisite_scatter.png', dpi=150, bbox_inches='tight')
print(f"Saved: multisite_scatter.png")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 3: ROLLING BIAS vs STARTING BG
# Shows how prediction error (pred - actual) varies across the BG range
# ══════════════════════════════════════════════════════════════════════════════

def plot_rolling_bias(ax, bg, actual, preds_dict, title):
    """Rolling mean of (predicted - actual) as a function of starting BG."""
    # Sort by BG and compute rolling mean with window
    order = np.argsort(bg)
    bg_sorted = bg[order]

    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.axvline(x=105, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)

    for (fname, color), pred in zip(FORMULAS, preds_dict.values()):
        err_sorted = (pred - actual)[order]
        # Use BG bins instead of rolling window for cleaner display
        bin_centers, bin_means = [], []
        for lo in range(72, 196, 4):
            hi = lo + 8
            mask = (bg_sorted >= lo) & (bg_sorted < hi)
            if mask.sum() >= 10:
                bin_centers.append((lo + hi) / 2)
                bin_means.append(err_sorted[mask].mean())
        if bin_centers:
            ax.plot(bin_centers, bin_means, color=color, linewidth=2, label=fname, alpha=0.85)

    ax.set_xlabel('Starting BG [mg/dL]', fontsize=10)
    ax.set_ylabel('Predicted − Actual 2h BG [mg/dL]', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Annotate
    ax.text(0.02, 0.98, '↑ Over-predicts glucose (conservative)',
            transform=ax.transAxes, fontsize=7, va='top', color='#666')
    ax.text(0.02, 0.02, '↓ Under-predicts glucose (aggressive)',
            transform=ax.transAxes, fontsize=7, va='bottom', color='#666')


print("Generating rolling bias plots...")
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))

for ax, (group_name, sites) in zip(axes, [('Sigmoid', sigmoid_sites), ('Log', log_sites)]):
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    n = len(bg)
    plot_rolling_bias(ax, bg, act, preds, f'{group_name} — Bias vs Starting BG\n({len(sites)} sites, {n:,} samples)')

fig.tight_layout()
fig.savefig(OUT_DIR / 'multisite_bias_vs_bg.png', dpi=150, bbox_inches='tight')
print(f"Saved: multisite_bias_vs_bg.png")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 4: BOX PLOTS — error distribution by BG band for each formula
# ══════════════════════════════════════════════════════════════════════════════

def plot_error_boxes(axes_row, bg, actual, preds_dict, title_prefix):
    """Box plots of (predicted - actual) by BG band for each formula."""
    bands = [('<90', 0, 90), ('90-105', 90, 105), ('105-120', 105, 120),
             ('120-150', 120, 150), ('150+', 150, 999)]

    for ax, ((fname, fcolor), pred) in zip(axes_row, zip(FORMULAS, preds_dict.values())):
        err = pred - actual
        data, labels = [], []
        for bname, lo, hi in bands:
            mask = (bg >= lo) & (bg < hi)
            if mask.sum() >= 5:
                data.append(err[mask])
                labels.append(f'{bname}\n(n={mask.sum()})')

        bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False,
                       medianprops=dict(color='black', linewidth=1.5),
                       whiskerprops=dict(linewidth=1),
                       boxprops=dict(facecolor=fcolor, alpha=0.4))

        ax.axhline(y=0, color='black', linewidth=0.8)
        ax.set_ylabel('Pred − Actual [mg/dL]', fontsize=9)
        ax.set_title(f'{title_prefix} — {fname}', fontsize=10, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        ax.tick_params(axis='x', labelsize=8)


print("Generating box plots...")
fig, axes = plt.subplots(2, 4, figsize=(20, 10))

for row, (group_name, sites) in enumerate([('Sigmoid', sigmoid_sites), ('Log', log_sites)]):
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    plot_error_boxes(axes[row], bg, act, preds, group_name)

fig.suptitle('Prediction Error Distribution by BG Band', fontsize=13, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(OUT_DIR / 'multisite_error_boxes.png', dpi=150, bbox_inches='tight')
print(f"Saved: multisite_error_boxes.png")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 5: COMBINED OVERVIEW — 2x2 grid for each model type
# Top-left: calibration, Top-right: bias vs BG
# Bottom: scatter for loop vs best alternative
# ══════════════════════════════════════════════════════════════════════════════

for group_name, sites in [('Sigmoid', sigmoid_sites), ('Log', log_sites)]:
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    n = len(bg)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Top-left: calibration
    plot_calibration(axes[0, 0], act, preds,
                     f'{group_name} — Calibration ({n:,} samples)')

    # Top-right: rolling bias
    plot_rolling_bias(axes[0, 1], bg, act, preds,
                      f'{group_name} — Bias vs Starting BG')

    # Bottom-left: scatter Loop
    ax = axes[1, 0]
    ax.plot([40, 260], [40, 260], 'k--', linewidth=0.8, alpha=0.4)
    falling = pd_ > 0
    ax.scatter(pl[falling], act[falling], s=4, alpha=0.15, color='#E74C3C',
               label='Falling', rasterized=True)
    ax.scatter(pl[~falling], act[~falling], s=4, alpha=0.15, color='#2E86AB',
               label='Rising', rasterized=True)
    ax.set_xlabel('Loop Predicted 2h BG', fontsize=10)
    ax.set_ylabel('Actual 2h BG', fontsize=10)
    ax.set_title(f'{group_name} — Loop (DynISF)', fontsize=11, fontweight='bold')
    ax.set_xlim(40, 260); ax.set_ylim(40, 260)
    ax.set_aspect('equal'); ax.grid(alpha=0.2)
    ax.legend(fontsize=8, markerscale=3)

    # Bottom-right: scatter Quartic (closest competitor)
    ax = axes[1, 1]
    ax.plot([40, 260], [40, 260], 'k--', linewidth=0.8, alpha=0.4)
    ax.scatter(pq[falling], act[falling], s=4, alpha=0.15, color='#E74C3C',
               label='Falling', rasterized=True)
    ax.scatter(pq[~falling], act[~falling], s=4, alpha=0.15, color='#2E86AB',
               label='Rising', rasterized=True)
    ax.set_xlabel('Quartic Predicted 2h BG', fontsize=10)
    ax.set_ylabel('Actual 2h BG', fontsize=10)
    ax.set_title(f'{group_name} — Quartic', fontsize=11, fontweight='bold')
    ax.set_xlim(40, 260); ax.set_ylim(40, 260)
    ax.set_aspect('equal'); ax.grid(alpha=0.2)
    ax.legend(fontsize=8, markerscale=3)

    fig.tight_layout()
    fname = f'multisite_{group_name.lower()}_overview.png'
    fig.savefig(OUT_DIR / fname, dpi=150, bbox_inches='tight')
    print(f"Saved: {fname}")

print("\nDONE — all charts saved")
