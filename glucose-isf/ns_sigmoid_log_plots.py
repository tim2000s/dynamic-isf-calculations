#!/usr/bin/env python3
"""
Sigmoid vs Log analysis with Diabeloop-poster-style ISF charts.

For each dynamic ISF type (sigmoid, log):
  - Bar chart: median actual ISF per BG bin, with IQR whiskers
  - Overlaid curves: Quartic, Full Diabeloop, Hybrid (all TDD-scaled per site)
  - Falling vs rising analysis table

Also produces per-site ISF charts for each site.
"""

import json, time, re, warnings
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

# ── Sites ──
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

# ── API fetching ──
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
        print(f" SKIP ({len(raw_ds) if raw_ds else 0} ds)")
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
        print(f" SKIP ({len(cycles)} cycles)")
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

    # Overnight fasting filter
    mask = (
        (df['hour'] < 8) &
        (df['cob'].fillna(99) == 0) &
        (df['bg'] >= 72) & (df['bg'] <= 200) &
        (df['bolus_age_min'] >= 120) &
        (df['tdd'].notna()) & (df['tdd'] > 0)
    )
    on = df[mask].copy()
    on = on.dropna(subset=['pred_iob_24', 'actual_bg_2h']).copy()
    if len(on) < 10:
        print(f" SKIP ({len(on)} samples)")
        return None

    on['pred_drop'] = on['bg'] - on['pred_iob_24']
    on['actual_drop'] = on['bg'] - on['actual_bg_2h']

    strict = on[
        (np.abs(on['pred_drop']) > 3) &
        ((on['actual_drop'] / on['pred_drop']).between(0, 5))
    ].copy()
    strict = strict.sort_values('ts').reset_index(drop=True)
    if len(strict) < 10:
        print(f" SKIP ({len(strict)} quality)")
        return None

    bg = strict['bg'].values
    isf_actual = strict['isf_actual'].values
    pred_drop = strict['pred_drop'].values
    actual_bg_2h = strict['actual_bg_2h'].values
    tdd_7d = strict['tdd_7day'].values
    tdd_median = np.median(tdd_7d)

    S_q = (1800.0 / tdd_median) / QUARTIC_AT_99
    S_db = (1800.0 / tdd_median) / FULL_DB_AT_99
    S_h = (1800.0 / tdd_median) / HYBRID_AT_105

    print(f" OK ({len(strict)} samples, TDD={tdd_median:.1f})")

    return {
        'name': name, 'model': model, 'n': len(strict),
        'tdd_median': tdd_median,
        'bg': bg, 'isf_actual': isf_actual,
        'pred_drop': pred_drop, 'actual_bg_2h': actual_bg_2h,
        'S_q': S_q, 'S_db': S_db, 'S_h': S_h,
    }


# ══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

BG_BINS = list(range(72, 201, 8))  # 72, 80, 88, ..., 200
BG_BIN_CENTERS = [(BG_BINS[i] + BG_BINS[i+1]) / 2 for i in range(len(BG_BINS)-1)]

COLORS = {
    'bars': '#2E86AB',
    'quartic': '#E74C3C',
    'full_db': '#27AE60',
    'hybrid': '#F39C12',
}


def plot_isf_chart(ax, bg_all, isf_all, S_q_weighted, S_db_weighted, S_h_weighted,
                   title, show_legend=True):
    """Plot Diabeloop-poster-style ISF chart on given axes."""

    # Bin the actual ISF data
    medians, q25s, q75s, counts = [], [], [], []
    for i in range(len(BG_BINS) - 1):
        lo, hi = BG_BINS[i], BG_BINS[i+1]
        mask = (bg_all >= lo) & (bg_all < hi)
        vals = isf_all[mask]
        if len(vals) >= 3:
            medians.append(np.median(vals))
            q25s.append(np.percentile(vals, 25))
            q75s.append(np.percentile(vals, 75))
            counts.append(len(vals))
        else:
            medians.append(np.nan)
            q25s.append(np.nan)
            q75s.append(np.nan)
            counts.append(0)

    medians = np.array(medians)
    q25s = np.array(q25s)
    q75s = np.array(q75s)

    valid = ~np.isnan(medians)
    bar_width = 6.5

    # Bars with IQR whiskers
    ax.bar(np.array(BG_BIN_CENTERS)[valid], medians[valid], width=bar_width,
           color=COLORS['bars'], alpha=0.6, label='Actual DynISF (median)', zorder=2)
    ax.errorbar(np.array(BG_BIN_CENTERS)[valid], medians[valid],
                yerr=[medians[valid] - q25s[valid], q75s[valid] - medians[valid]],
                fmt='none', ecolor='#1A5276', elinewidth=1.2, capsize=3, zorder=3)

    # Formula curves
    bg_curve = np.linspace(72, 200, 200)
    isf_q_curve = isf_quartic(bg_curve) * S_q_weighted
    isf_db_curve = isf_full_diabeloop(bg_curve) * S_db_weighted
    isf_h_curve = isf_hybrid(bg_curve) * S_h_weighted

    ax.plot(bg_curve, isf_q_curve, color=COLORS['quartic'], linewidth=2,
            linestyle='--', label='Quartic (TDD-scaled)', zorder=4)
    ax.plot(bg_curve, isf_db_curve, color=COLORS['full_db'], linewidth=2,
            linestyle='-.', label='Full Diabeloop (TDD-scaled)', zorder=4)
    ax.plot(bg_curve, isf_h_curve, color=COLORS['hybrid'], linewidth=2,
            linestyle=':', label='Hybrid (TDD-scaled)', zorder=4)

    # Vertical line at 105
    ax.axvline(x=105, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.text(106, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else 80,
            'BG=105', fontsize=7, color='gray', va='top')

    ax.set_xlabel('Blood Glucose [mg/dL]', fontsize=10)
    ax.set_ylabel('ISF [mg/dL/U]', fontsize=10)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlim(68, 204)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.grid(axis='y', alpha=0.3)

    if show_legend:
        ax.legend(fontsize=8, loc='upper right')


def plot_error_chart(ax, sites, title):
    """Plot falling/rising prediction error summary as grouped bars."""
    combos = [
        ("<105\nfalling",  lambda bg: bg < 105,  lambda pd: pd > 0),
        ("<105\nrising",   lambda bg: bg < 105,  lambda pd: pd < 0),
        ("≥105\nfalling",  lambda bg: bg >= 105, lambda pd: pd > 0),
        ("≥105\nrising",   lambda bg: bg >= 105, lambda pd: pd < 0),
    ]

    labels = []
    loop_vals, q_vals, db_vals, h_vals = [], [], [], []

    for label, bg_fn, dr_fn in combos:
        all_actual, all_pl, all_pq, all_pd, all_ph = [], [], [], [], []
        for s in sites:
            bg = s['bg']
            pd_ = s['pred_drop']
            m = bg_fn(bg) & dr_fn(pd_)
            if m.sum() == 0: continue

            isf_a = s['isf_actual']
            isf_q = isf_quartic(bg) * s['S_q']
            isf_db = isf_full_diabeloop(bg) * s['S_db']
            isf_h = isf_hybrid(bg) * s['S_h']

            pred_l = bg - pd_ * (isf_a / isf_a)
            pred_q = bg - pd_ * (isf_q / isf_a)
            pred_d = bg - pd_ * (isf_db / isf_a)
            pred_h = bg - pd_ * (isf_h / isf_a)

            act = s['actual_bg_2h']
            all_actual.extend(act[m])
            all_pl.extend(pred_l[m]); all_pq.extend(pred_q[m])
            all_pd.extend(pred_d[m]); all_ph.extend(pred_h[m])

        if len(all_actual) == 0:
            labels.append(label)
            loop_vals.append(0); q_vals.append(0); db_vals.append(0); h_vals.append(0)
            continue

        act = np.array(all_actual)
        labels.append(f"{label}\n(n={len(act)})")
        loop_vals.append(np.mean(all_pl) - act.mean())
        q_vals.append(np.mean(all_pq) - act.mean())
        db_vals.append(np.mean(all_pd) - act.mean())
        h_vals.append(np.mean(all_ph) - act.mean())

    x = np.arange(len(labels))
    w = 0.2
    ax.bar(x - 1.5*w, loop_vals, w, label='Loop (DynISF)', color=COLORS['bars'], alpha=0.8)
    ax.bar(x - 0.5*w, q_vals, w, label='Quartic', color=COLORS['quartic'], alpha=0.8)
    ax.bar(x + 0.5*w, db_vals, w, label='Full Diabeloop', color=COLORS['full_db'], alpha=0.8)
    ax.bar(x + 1.5*w, h_vals, w, label='Hybrid', color=COLORS['hybrid'], alpha=0.8)

    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Predicted − Actual 2h BG [mg/dL]', fontsize=9)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=7, loc='best')
    ax.grid(axis='y', alpha=0.3)

    # Annotate: above zero = over-predicts, below = under-predicts
    ylim = ax.get_ylim()
    ax.text(0.02, 0.98, '↑ Predicts HIGHER than actual', transform=ax.transAxes,
            fontsize=7, va='top', color='#666')
    ax.text(0.02, 0.02, '↓ Predicts LOWER than actual', transform=ax.transAxes,
            fontsize=7, va='bottom', color='#666')


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

print("Fetching data from all sites...")
all_sites = []
for site in SITES:
    r = process_site(site)
    if r is not None:
        all_sites.append(r)

sigmoid_sites = [s for s in all_sites if s['model'] == 'sigmoid']
log_sites = [s for s in all_sites if s['model'] == 'log']

print(f"\nSigmoid: {len(sigmoid_sites)} sites, {sum(s['n'] for s in sigmoid_sites)} samples")
print(f"Log: {len(log_sites)} sites, {sum(s['n'] for s in log_sites)} samples")


# ── Compute weighted-average scaling factors for group curves ──
def weighted_S(sites, key):
    total_n = sum(s['n'] for s in sites)
    return sum(s[key] * s['n'] / total_n for s in sites)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Sigmoid — ISF chart + error chart (2 panels)
# ══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Pool all sigmoid data
bg_sig = np.concatenate([s['bg'] for s in sigmoid_sites])
isf_sig = np.concatenate([s['isf_actual'] for s in sigmoid_sites])
S_q_sig = weighted_S(sigmoid_sites, 'S_q')
S_db_sig = weighted_S(sigmoid_sites, 'S_db')
S_h_sig = weighted_S(sigmoid_sites, 'S_h')

plot_isf_chart(axes[0], bg_sig, isf_sig, S_q_sig, S_db_sig, S_h_sig,
               f'Sigmoid DynISF — Actual vs Formulas\n({len(sigmoid_sites)} sites, {len(bg_sig):,} samples)')

plot_error_chart(axes[1], sigmoid_sites,
                 'Sigmoid — Prediction Error by Zone & Direction')

fig.tight_layout()
fig.savefig(OUT_DIR / 'multisite_sigmoid_analysis.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR / 'multisite_sigmoid_analysis.png'}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Log — ISF chart + error chart (2 panels)
# ══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

bg_log = np.concatenate([s['bg'] for s in log_sites])
isf_log = np.concatenate([s['isf_actual'] for s in log_sites])
S_q_log = weighted_S(log_sites, 'S_q')
S_db_log = weighted_S(log_sites, 'S_db')
S_h_log = weighted_S(log_sites, 'S_h')

plot_isf_chart(axes[0], bg_log, isf_log, S_q_log, S_db_log, S_h_log,
               f'Log DynISF — Actual vs Formulas\n({len(log_sites)} sites, {len(bg_log):,} samples)')

plot_error_chart(axes[1], log_sites,
                 'Log — Prediction Error by Zone & Direction')

fig.tight_layout()
fig.savefig(OUT_DIR / 'multisite_log_analysis.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'multisite_log_analysis.png'}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Per-site ISF charts (individual panels)
# ══════════════════════════════════════════════════════════════════════════════

n_sites = len(all_sites)
ncols = 3
nrows = (n_sites + ncols - 1) // ncols

fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
axes = axes.flatten()

for i, s in enumerate(sorted(all_sites, key=lambda x: (x['model'], x['name']))):
    plot_isf_chart(axes[i], s['bg'], s['isf_actual'],
                   s['S_q'], s['S_db'], s['S_h'],
                   f"{s['name']} ({s['model']})\n{s['n']} samples, TDD={s['tdd_median']:.1f}",
                   show_legend=(i == 0))
    # Set y-axis upper limit sensibly
    isf_p95 = np.percentile(s['isf_actual'], 95)
    hybrid_max = isf_hybrid(72) * s['S_h']
    y_max = min(max(isf_p95 * 1.5, 150), 500)
    axes[i].set_ylim(0, y_max)

# Hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle('Per-Site ISF Curves: Actual DynISF vs Diabeloop Formulas',
             fontsize=14, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(OUT_DIR / 'multisite_persite_isf.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'multisite_persite_isf.png'}")


# ══════════════════════════════════════════════════════════════════════════════
# TEXT SUMMARY — Sigmoid vs Log falling/rising
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("SIGMOID vs LOG — Falling/Rising Summary (Pred − Actual)")
print("=" * 80)

combos = [
    ("<105 falling",  lambda bg: bg < 105,  lambda pd: pd > 0),
    ("<105 rising",   lambda bg: bg < 105,  lambda pd: pd < 0),
    ("≥105 falling",  lambda bg: bg >= 105, lambda pd: pd > 0),
    ("≥105 rising",   lambda bg: bg >= 105, lambda pd: pd < 0),
]

for group_name, group_sites in [("SIGMOID", sigmoid_sites), ("LOG", log_sites)]:
    print(f"\n  {group_name} ({len(group_sites)} sites, {sum(s['n'] for s in group_sites)} samples)")
    print(f"  {'Zone':<16s} {'N':>6s}  {'Loop':>8s} {'Quartic':>8s} {'FullDB':>8s} {'Hybrid':>8s}")
    print(f"  {'─'*16} {'─'*6}  {'─'*8} {'─'*8} {'─'*8} {'─'*8}")

    for label, bg_fn, dr_fn in combos:
        all_act, all_pl, all_pq, all_pd, all_ph = [], [], [], [], []
        for s in group_sites:
            bg = s['bg']; pd_ = s['pred_drop']
            m = bg_fn(bg) & dr_fn(pd_)
            if m.sum() == 0: continue
            isf_a = s['isf_actual']
            isf_q = isf_quartic(bg) * s['S_q']
            isf_db = isf_full_diabeloop(bg) * s['S_db']
            isf_h = isf_hybrid(bg) * s['S_h']
            act = s['actual_bg_2h']
            pred_l = bg - pd_ * (isf_a / isf_a)
            pred_q = bg - pd_ * (isf_q / isf_a)
            pred_d = bg - pd_ * (isf_db / isf_a)
            pred_h = bg - pd_ * (isf_h / isf_a)
            all_act.extend(act[m]); all_pl.extend(pred_l[m])
            all_pq.extend(pred_q[m]); all_pd.extend(pred_d[m]); all_ph.extend(pred_h[m])

        if len(all_act) == 0:
            print(f"  {label:<16s}    —")
            continue
        act = np.mean(all_act)
        n = len(all_act)
        lv = np.mean(all_pl) - act
        qv = np.mean(all_pq) - act
        dv = np.mean(all_pd) - act
        hv = np.mean(all_ph) - act
        print(f"  {label:<16s} {n:6d}  {lv:+8.1f} {qv:+8.1f} {dv:+8.1f} {hv:+8.1f}")

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
