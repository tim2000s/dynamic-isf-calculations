#!/usr/bin/env python3
"""
4-Hour Multi-Site ISF Analysis — Combined Script
=================================================

Uses the LAST element of predBGs.IOB (= end of insulin action) as the
prediction horizon, and matches actual CGM at the same time offset.

Produces:
  1. Per-site backtest results + aggregate summary
  2. Falling vs rising analysis (<105 / ≥105)
  3. Sigmoid vs log ISF charts (Diabeloop-poster style)
  4. Prediction performance charts (calibration, scatter, bias, boxes)

All output goes to ~/Downloads/4 Hour analysis/
"""

import json, time, re, warnings, pickle, sys
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

OUT_DIR = Path.home() / 'Downloads' / '4 Hour analysis'
CACHE_FILE = OUT_DIR / 'multisite_4h_sample_cache.pkl'

# ══════════════════════════════════════════════════════════════════════════════
# ISF FORMULAS
# ══════════════════════════════════════════════════════════════════════════════

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
# Nightscout sites + tokens are loaded from ns_sites.json (gitignored — never commit credentials).
import json as _json, os as _os
SITES = _json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "ns_sites.json")))

# ══════════════════════════════════════════════════════════════════════════════
# API FETCHING
# ══════════════════════════════════════════════════════════════════════════════

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
        sys.stdout.write(f"\r    {len(all_records):,} records (back to {win_start[:10]})...")
        sys.stdout.flush()
    if all_records:
        sys.stdout.write(f"\r    {len(all_records):,} records total                    \n")
        sys.stdout.flush()
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


# ══════════════════════════════════════════════════════════════════════════════
# PROCESS SITE — uses last element of predBGs.IOB
# ══════════════════════════════════════════════════════════════════════════════

def process_site(site):
    name, base_url, token, model = site['name'], site['url'], site['token'], site['model']
    print(f"\n  Fetching {name}...", end='', flush=True)

    test = ns_fetch(base_url, 'devicestatus.json', token, {'count': 1})
    if test is None:
        print(" SKIP (connect fail)")
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

    # Parse loop cycles
    cycles = []
    pred_lens = []
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
            if len(pred_iob) < 12:  # Need at least 1h of prediction
                continue
            tdd_val = sg.get('TDD')
            if tdd_val is None or tdd_val <= 0:
                tdd_val = parse_tdd_from_reason(reason)

            pred_iob_final = pred_iob[-1]
            pred_horizon_s = (len(pred_iob) - 1) * 5 * 60
            pred_lens.append(len(pred_iob))

            cycles.append({
                'ts': ts, 'bg': bg_val, 'isf_actual': float(isf_val),
                'cob': sg.get('COB'), 'iob': sg.get('IOB'),
                'tdd': float(tdd_val) if tdd_val else None,
                'pred_iob_final': float(pred_iob_final),
                'pred_horizon_s': pred_horizon_s,
            })
        except Exception:
            continue

    if len(cycles) < 100:
        print(f" SKIP ({len(cycles)} cycles)")
        return None

    df = pd.DataFrame(cycles).sort_values('ts').reset_index(drop=True)
    median_horizon_min = np.median([c['pred_horizon_s'] for c in cycles]) / 60

    # Parse CGM
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

    # Forward BG at each cycle's own prediction horizon
    cycle_epochs = np.array([int(t.timestamp()) for t in df['ts']])
    actual_end = np.full(len(df), np.nan)
    bolus_age = np.full(len(df), np.nan)
    for i, t_s in enumerate(cycle_epochs):
        horizon = int(df.iloc[i]['pred_horizon_s'])
        actual_end[i] = get_cgm_at(t_s + horizon)
        bolus_age[i] = mins_since_bolus(t_s)
    df['actual_bg_end'] = actual_end
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
    on = on.dropna(subset=['pred_iob_final', 'actual_bg_end']).copy()
    if len(on) < 10:
        print(f" SKIP ({len(on)} samples)")
        return None

    on['pred_drop'] = on['bg'] - on['pred_iob_final']
    on['actual_drop'] = on['bg'] - on['actual_bg_end']

    strict = on[
        (np.abs(on['pred_drop']) > 3) &
        ((on['actual_drop'] / on['pred_drop']).between(0, 5))
    ].copy().sort_values('ts').reset_index(drop=True)
    if len(strict) < 10:
        print(f" SKIP ({len(strict)} quality)")
        return None

    bg = strict['bg'].values
    isf_actual = strict['isf_actual'].values
    pred_drop = strict['pred_drop'].values
    actual_bg_end = strict['actual_bg_end'].values
    tdd_7d = strict['tdd_7day'].values
    tdd_median = np.median(tdd_7d)

    S_q = (1800.0 / tdd_median) / QUARTIC_AT_99
    S_db = (1800.0 / tdd_median) / FULL_DB_AT_99
    S_h = (1800.0 / tdd_median) / HYBRID_AT_105

    isf_q = isf_quartic(bg) * S_q
    isf_db = isf_full_diabeloop(bg) * S_db
    isf_h = isf_hybrid(bg) * S_h

    pred_loop = bg - pred_drop * 1.0
    pred_q = bg - pred_drop * (isf_q / isf_actual)
    pred_db = bg - pred_drop * (isf_db / isf_actual)
    pred_h = bg - pred_drop * (isf_h / isf_actual)

    print(f" OK ({len(strict)} samples, TDD={tdd_median:.1f}, horizon={median_horizon_min:.0f}min)")

    return {
        'name': name, 'model': model, 'n': len(strict),
        'tdd_median': tdd_median, 'median_horizon_min': median_horizon_min,
        'bg': bg, 'isf_actual': isf_actual,
        'pred_drop': pred_drop, 'actual_bg_end': actual_bg_end,
        'pred_loop': pred_loop, 'pred_q': pred_q, 'pred_db': pred_db, 'pred_h': pred_h,
        'S_q': S_q, 'S_db': S_db, 'S_h': S_h,
    }


# ══════════════════════════════════════════════════════════════════════════════
# LOAD OR FETCH
# ══════════════════════════════════════════════════════════════════════════════

if CACHE_FILE.exists():
    print(f"Loading cached data from {CACHE_FILE}")
    with open(CACHE_FILE, 'rb') as f:
        all_sites = pickle.load(f)
    print(f"Loaded {len(all_sites)} sites, {sum(s['n'] for s in all_sites)} samples")
else:
    print("Fetching data from all sites (end-of-IOB prediction horizon)...")
    all_sites = []
    for site in SITES:
        r = process_site(site)
        if r is not None:
            all_sites.append(r)
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(all_sites, f)
    print(f"\nCached to {CACHE_FILE}")

sigmoid_sites = [s for s in all_sites if s['model'] == 'sigmoid']
log_sites = [s for s in all_sites if s['model'] == 'log']

total_n = sum(s['n'] for s in all_sites)
print(f"\nTotal: {len(all_sites)} sites, {total_n} samples")
print(f"Sigmoid: {len(sigmoid_sites)} sites, {sum(s['n'] for s in sigmoid_sites)} samples")
print(f"Log: {len(log_sites)} sites, {sum(s['n'] for s in log_sites)} samples")

for s in all_sites:
    print(f"  {s['name']:<14s} {s['model']:<8s} n={s['n']:5d}  TDD={s['tdd_median']:5.1f}  horizon={s['median_horizon_min']:.0f}min")


# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

FORMULAS = [
    ('Loop (DynISF)', '#2E86AB'),
    ('Quartic',       '#E74C3C'),
    ('Full Diabeloop','#27AE60'),
    ('Hybrid',        '#F39C12'),
]

def pool(sites):
    bg = np.concatenate([s['bg'] for s in sites])
    act = np.concatenate([s['actual_bg_end'] for s in sites])
    pd_ = np.concatenate([s['pred_drop'] for s in sites])
    pl = np.concatenate([s['pred_loop'] for s in sites])
    pq = np.concatenate([s['pred_q'] for s in sites])
    pdb = np.concatenate([s['pred_db'] for s in sites])
    ph = np.concatenate([s['pred_h'] for s in sites])
    return bg, act, pd_, pl, pq, pdb, ph

def weighted_S(sites, key):
    total_n = sum(s['n'] for s in sites)
    return sum(s[key] * s['n'] / total_n for s in sites)


# ══════════════════════════════════════════════════════════════════════════════
# 1. BACKTEST SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("BACKTEST SUMMARY — End-of-IOB Prediction Horizon")
print("=" * 80)

print(f"\n  {'Site':<14s} {'Model':<8s} {'N':>5s} {'TDD':>5s} {'Hrz':>4s}  {'Loop':>7s} {'Quart':>7s} {'FullDB':>7s} {'Hybrid':>7s}")
print(f"  {'─'*14} {'─'*8} {'─'*5} {'─'*5} {'─'*4}  {'─'*7} {'─'*7} {'─'*7} {'─'*7}")

for s in all_sites:
    bg, act = s['bg'], s['actual_bg_end']
    mae_l = np.abs(act - s['pred_loop']).mean()
    mae_q = np.abs(act - s['pred_q']).mean()
    mae_d = np.abs(act - s['pred_db']).mean()
    mae_h = np.abs(act - s['pred_h']).mean()
    maes = {'Loop': mae_l, 'Quart': mae_q, 'FullDB': mae_d, 'Hybrid': mae_h}
    best = min(maes, key=maes.get)
    marker = lambda k: '*' if k == best else ' '
    print(f"  {s['name']:<14s} {s['model']:<8s} {s['n']:5d} {s['tdd_median']:5.1f} {s['median_horizon_min']:4.0f}"
          f"  {mae_l:6.1f}{marker('Loop')} {mae_q:6.1f}{marker('Quart')} {mae_d:6.1f}{marker('FullDB')} {mae_h:6.1f}{marker('Hybrid')}")


# ══════════════════════════════════════════════════════════════════════════════
# 2. FALLING vs RISING ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("FALLING vs RISING — Pred−Actual (positive = over-predicts glucose)")
print("=" * 80)

combos = [
    ("<105 falling",  lambda bg: bg < 105,  lambda pd: pd > 0),
    ("<105 rising",   lambda bg: bg < 105,  lambda pd: pd < 0),
    ("≥105 falling",  lambda bg: bg >= 105, lambda pd: pd > 0),
    ("≥105 rising",   lambda bg: bg >= 105, lambda pd: pd < 0),
]

for group_name, group_sites in [("SIGMOID", sigmoid_sites), ("LOG", log_sites), ("ALL", all_sites)]:
    print(f"\n  {group_name} ({len(group_sites)} sites, {sum(s['n'] for s in group_sites)} samples)")
    print(f"  {'Zone':<16s} {'N':>6s}  {'Loop':>8s} {'Quartic':>8s} {'FullDB':>8s} {'Hybrid':>8s}")
    print(f"  {'─'*16} {'─'*6}  {'─'*8} {'─'*8} {'─'*8} {'─'*8}")

    for label, bg_fn, dr_fn in combos:
        all_act, all_pl, all_pq, all_pd, all_ph = [], [], [], [], []
        for s in group_sites:
            m = bg_fn(s['bg']) & dr_fn(s['pred_drop'])
            if m.sum() == 0: continue
            all_act.extend(s['actual_bg_end'][m])
            all_pl.extend(s['pred_loop'][m]); all_pq.extend(s['pred_q'][m])
            all_pd.extend(s['pred_db'][m]); all_ph.extend(s['pred_h'][m])

        if len(all_act) == 0:
            print(f"  {label:<16s}    —")
            continue
        act = np.mean(all_act); n = len(all_act)
        lv = np.mean(all_pl) - act; qv = np.mean(all_pq) - act
        dv = np.mean(all_pd) - act; hv = np.mean(all_ph) - act
        print(f"  {label:<16s} {n:6d}  {lv:+8.1f} {qv:+8.1f} {dv:+8.1f} {hv:+8.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# 3. BELOW-90 DETAIL
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("BELOW 90 — Mean Predicted vs Actual end-of-IOB BG")
print("  Pred−Actual: positive = over-predicts (conservative), negative = under-predicts")
print("=" * 80)

print(f"\n  {'Site':<14s} {'Model':<8s} {'N':>5s} {'Actual':>7s}  {'Loop':>7s} {'Quart':>7s} {'FullDB':>7s} {'Hybrid':>7s}")
print(f"  {'─'*14} {'─'*8} {'─'*5} {'─'*7}  {'─'*7} {'─'*7} {'─'*7} {'─'*7}")

for s in all_sites:
    m = s['bg'] < 90
    if m.sum() < 5: continue
    act = s['actual_bg_end'][m].mean()
    lv = s['pred_loop'][m].mean() - act
    qv = s['pred_q'][m].mean() - act
    dv = s['pred_db'][m].mean() - act
    hv = s['pred_h'][m].mean() - act
    print(f"  {s['name']:<14s} {s['model']:<8s} {m.sum():5d} {act:7.0f}  {lv:+7.1f} {qv:+7.1f} {dv:+7.1f} {hv:+7.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# CHARTS
# ══════════════════════════════════════════════════════════════════════════════

BG_BINS = list(range(72, 201, 8))
BG_BIN_CENTERS = [(BG_BINS[i] + BG_BINS[i+1]) / 2 for i in range(len(BG_BINS)-1)]
COLORS = {'bars': '#2E86AB', 'quartic': '#E74C3C', 'full_db': '#27AE60', 'hybrid': '#F39C12'}


def plot_isf_chart(ax, bg_all, isf_all, S_q, S_db, S_h, title, show_legend=True):
    medians, q25s, q75s = [], [], []
    for i in range(len(BG_BINS) - 1):
        lo, hi = BG_BINS[i], BG_BINS[i+1]
        vals = isf_all[(bg_all >= lo) & (bg_all < hi)]
        if len(vals) >= 3:
            medians.append(np.median(vals)); q25s.append(np.percentile(vals, 25)); q75s.append(np.percentile(vals, 75))
        else:
            medians.append(np.nan); q25s.append(np.nan); q75s.append(np.nan)
    medians, q25s, q75s = np.array(medians), np.array(q25s), np.array(q75s)
    valid = ~np.isnan(medians)
    ax.bar(np.array(BG_BIN_CENTERS)[valid], medians[valid], width=6.5,
           color=COLORS['bars'], alpha=0.6, label='Actual DynISF (median)', zorder=2)
    ax.errorbar(np.array(BG_BIN_CENTERS)[valid], medians[valid],
                yerr=[medians[valid]-q25s[valid], q75s[valid]-medians[valid]],
                fmt='none', ecolor='#1A5276', elinewidth=1.2, capsize=3, zorder=3)
    bg_curve = np.linspace(72, 200, 200)
    ax.plot(bg_curve, isf_quartic(bg_curve)*S_q, color=COLORS['quartic'], linewidth=2,
            linestyle='--', label='Quartic (TDD-scaled)', zorder=4)
    ax.plot(bg_curve, isf_full_diabeloop(bg_curve)*S_db, color=COLORS['full_db'], linewidth=2,
            linestyle='-.', label='Full Diabeloop (TDD-scaled)', zorder=4)
    ax.plot(bg_curve, isf_hybrid(bg_curve)*S_h, color=COLORS['hybrid'], linewidth=2,
            linestyle=':', label='Hybrid (TDD-scaled)', zorder=4)
    ax.axvline(x=105, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.set_xlabel('Sensor Glucose [mg/dL]'); ax.set_ylabel('ISF [mg/dL/U]')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlim(68, 204); ax.xaxis.set_major_locator(MultipleLocator(20)); ax.grid(axis='y', alpha=0.3)
    if show_legend: ax.legend(fontsize=8, loc='upper right')


def plot_error_bars(ax, sites, title):
    labels, loop_vals, q_vals, db_vals, h_vals = [], [], [], [], []
    for label, bg_fn, dr_fn in combos:
        all_act, all_pl, all_pq, all_pd, all_ph = [], [], [], [], []
        for s in sites:
            m = bg_fn(s['bg']) & dr_fn(s['pred_drop'])
            if m.sum() == 0: continue
            all_act.extend(s['actual_bg_end'][m])
            all_pl.extend(s['pred_loop'][m]); all_pq.extend(s['pred_q'][m])
            all_pd.extend(s['pred_db'][m]); all_ph.extend(s['pred_h'][m])
        if not all_act: labels.append(label); loop_vals.append(0); q_vals.append(0); db_vals.append(0); h_vals.append(0); continue
        act = np.mean(all_act)
        labels.append(f"{label}\n(n={len(all_act)})")
        loop_vals.append(np.mean(all_pl)-act); q_vals.append(np.mean(all_pq)-act)
        db_vals.append(np.mean(all_pd)-act); h_vals.append(np.mean(all_ph)-act)
    x = np.arange(len(labels)); w = 0.2
    ax.bar(x-1.5*w, loop_vals, w, label='Loop (DynISF)', color=COLORS['bars'], alpha=0.8)
    ax.bar(x-0.5*w, q_vals, w, label='Quartic', color=COLORS['quartic'], alpha=0.8)
    ax.bar(x+0.5*w, db_vals, w, label='Full Diabeloop', color=COLORS['full_db'], alpha=0.8)
    ax.bar(x+1.5*w, h_vals, w, label='Hybrid', color=COLORS['hybrid'], alpha=0.8)
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Predicted − Actual end-of-IOB SGV [mg/dL]', fontsize=9)
    ax.set_title(title, fontsize=12, fontweight='bold'); ax.legend(fontsize=7)
    ax.grid(axis='y', alpha=0.3)
    ax.text(0.02, 0.98, '↑ Predicts HIGHER than actual', transform=ax.transAxes, fontsize=7, va='top', color='#666')
    ax.text(0.02, 0.02, '↓ Predicts LOWER than actual', transform=ax.transAxes, fontsize=7, va='bottom', color='#666')


def plot_calibration(ax, actual, preds_dict, title):
    bins = np.arange(50, 220, 10)
    ax.plot([40, 250], [40, 250], 'k--', linewidth=1, alpha=0.5, label='Perfect')
    for (fname, color), pred in zip(FORMULAS, preds_dict.values()):
        bp, ba = [], []
        for i in range(len(bins)-1):
            mask = (pred >= bins[i]) & (pred < bins[i+1])
            if mask.sum() >= 5: bp.append(pred[mask].mean()); ba.append(actual[mask].mean())
        if bp: ax.plot(bp, ba, 'o-', color=color, markersize=4, linewidth=1.5, label=fname, alpha=0.8)
    ax.set_xlabel('Mean Predicted end-of-IOB SGV'); ax.set_ylabel('Mean Actual end-of-IOB SGV')
    ax.set_title(title, fontsize=11, fontweight='bold'); ax.legend(fontsize=8)
    ax.set_xlim(40, 220); ax.set_ylim(40, 220); ax.set_aspect('equal'); ax.grid(alpha=0.3)


def plot_rolling_bias(ax, bg, actual, preds_dict, title):
    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.axvline(x=105, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
    order = np.argsort(bg); bg_s = bg[order]
    for (fname, color), pred in zip(FORMULAS, preds_dict.values()):
        err_s = (pred - actual)[order]
        bc, bm = [], []
        for lo in range(72, 196, 4):
            mask = (bg_s >= lo) & (bg_s < lo+8)
            if mask.sum() >= 10: bc.append(lo+4); bm.append(err_s[mask].mean())
        if bc: ax.plot(bc, bm, color=color, linewidth=2, label=fname, alpha=0.85)
    ax.set_xlabel('Starting SGV [mg/dL]'); ax.set_ylabel('Predicted − Actual [mg/dL]')
    ax.set_title(title, fontsize=11, fontweight='bold'); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.text(0.02, 0.98, '↑ Over-predicts glucose (conservative)', transform=ax.transAxes, fontsize=7, va='top', color='#666')
    ax.text(0.02, 0.02, '↓ Under-predicts glucose', transform=ax.transAxes, fontsize=7, va='bottom', color='#666')


# ── Generate all charts ──
print("\nGenerating charts...")

# Chart 1+2: Sigmoid and Log ISF + error charts
for group_name, sites in [('Sigmoid', sigmoid_sites), ('Log', log_sites)]:
    bg_all = np.concatenate([s['bg'] for s in sites])
    isf_all = np.concatenate([s['isf_actual'] for s in sites])
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    plot_isf_chart(axes[0], bg_all, isf_all,
                   weighted_S(sites, 'S_q'), weighted_S(sites, 'S_db'), weighted_S(sites, 'S_h'),
                   f'{group_name} DynISF — Actual vs Formulas\n({len(sites)} sites, {len(bg_all):,} samples)')
    plot_error_bars(axes[1], sites, f'{group_name} — Prediction Error (end-of-IOB)')
    fig.tight_layout()
    fig.savefig(OUT_DIR / f'4h_{group_name.lower()}_analysis.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: 4h_{group_name.lower()}_analysis.png")

# Chart 3: Per-site ISF charts
n_sites = len(all_sites); ncols = 3; nrows = (n_sites + ncols - 1) // ncols
fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
axes_flat = axes.flatten()
for i, s in enumerate(sorted(all_sites, key=lambda x: (x['model'], x['name']))):
    plot_isf_chart(axes_flat[i], s['bg'], s['isf_actual'], s['S_q'], s['S_db'], s['S_h'],
                   f"{s['name']} ({s['model']})\n{s['n']} samples, TDD={s['tdd_median']:.1f}, {s['median_horizon_min']:.0f}min",
                   show_legend=(i==0))
    isf_p95 = np.percentile(s['isf_actual'], 95)
    axes_flat[i].set_ylim(0, min(max(isf_p95*1.5, 150), 500))
for j in range(i+1, len(axes_flat)): axes_flat[j].set_visible(False)
fig.suptitle('Per-Site ISF Curves — End-of-IOB Analysis', fontsize=14, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(OUT_DIR / '4h_persite_isf.png', dpi=150, bbox_inches='tight')
print(f"  Saved: 4h_persite_isf.png")

# Chart 4: Calibration
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
for ax, (gname, sites) in zip(axes, [('Sigmoid', sigmoid_sites), ('Log', log_sites)]):
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    plot_calibration(ax, act, preds, f'{gname} — Calibration ({len(bg):,} samples)')
fig.tight_layout()
fig.savefig(OUT_DIR / '4h_calibration.png', dpi=150, bbox_inches='tight')
print(f"  Saved: 4h_calibration.png")

# Chart 5: Rolling bias
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
for ax, (gname, sites) in zip(axes, [('Sigmoid', sigmoid_sites), ('Log', log_sites)]):
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    plot_rolling_bias(ax, bg, act, preds, f'{gname} — Bias vs Starting SGV (end-of-IOB)')
fig.tight_layout()
fig.savefig(OUT_DIR / '4h_bias_vs_bg.png', dpi=150, bbox_inches='tight')
print(f"  Saved: 4h_bias_vs_bg.png")

# Chart 6: Box plots
bands = [('<90', 0, 90), ('90-105', 90, 105), ('105-120', 105, 120), ('120-150', 120, 150), ('150+', 150, 999)]
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
for row, (gname, sites) in enumerate([('Sigmoid', sigmoid_sites), ('Log', log_sites)]):
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    for col, ((fname, fcolor), pred) in enumerate(zip(FORMULAS, [pl, pq, pdb, ph])):
        ax = axes[row, col]; err = pred - act
        data, labels = [], []
        for bname, lo, hi in bands:
            mask = (bg >= lo) & (bg < hi)
            if mask.sum() >= 5: data.append(err[mask]); labels.append(f'{bname}\n(n={mask.sum()})')
        ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False,
                   medianprops=dict(color='black', linewidth=1.5),
                   boxprops=dict(facecolor=fcolor, alpha=0.4))
        ax.axhline(y=0, color='black', linewidth=0.8)
        ax.set_ylabel('Pred − Actual [mg/dL]'); ax.set_title(f'{gname} — {fname}', fontsize=10, fontweight='bold')
        ax.grid(axis='y', alpha=0.3); ax.tick_params(axis='x', labelsize=8)
fig.suptitle('End-of-IOB Prediction Error Distribution by SGV Band', fontsize=13, fontweight='bold', y=1.01)
fig.tight_layout()
fig.savefig(OUT_DIR / '4h_error_boxes.png', dpi=150, bbox_inches='tight')
print(f"  Saved: 4h_error_boxes.png")

# Chart 7: Combined overview per model type
for gname, sites in [('Sigmoid', sigmoid_sites), ('Log', log_sites)]:
    bg, act, pd_, pl, pq, pdb, ph = pool(sites)
    preds = {'loop': pl, 'quartic': pq, 'full_db': pdb, 'hybrid': ph}
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    plot_calibration(axes[0,0], act, preds, f'{gname} — Calibration ({len(bg):,} samples)')
    plot_rolling_bias(axes[0,1], bg, act, preds, f'{gname} — Bias vs Starting SGV')
    for idx, (fname, pred, color) in enumerate([(f'Loop (DynISF)', pl, '#2E86AB'), ('Quartic', pq, '#E74C3C')]):
        ax = axes[1, idx]
        ax.plot([40, 260], [40, 260], 'k--', linewidth=0.8, alpha=0.4)
        falling = pd_ > 0
        ax.scatter(pred[falling], act[falling], s=4, alpha=0.15, color='#E74C3C', label='Falling', rasterized=True)
        ax.scatter(pred[~falling], act[~falling], s=4, alpha=0.15, color='#2E86AB', label='Rising', rasterized=True)
        ax.set_xlabel(f'{fname} Predicted SGV'); ax.set_ylabel('Actual SGV')
        ax.set_title(f'{gname} — {fname}', fontsize=11, fontweight='bold')
        ax.set_xlim(40, 260); ax.set_ylim(40, 260); ax.set_aspect('equal'); ax.grid(alpha=0.2)
        ax.legend(fontsize=8, markerscale=3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f'4h_{gname.lower()}_overview.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: 4h_{gname.lower()}_overview.png")

print("\n" + "=" * 80)
print("DONE — all 4h analysis complete")
print("=" * 80)
