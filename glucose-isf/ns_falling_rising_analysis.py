#!/usr/bin/env python3
"""
Analyse prediction errors split by falling vs rising glucose, <105 vs ≥105.

pred_drop = bg - pred_iob_24
  pred_drop > 0 → loop predicts glucose FALLS (active insulin)
  pred_drop < 0 → loop predicts glucose RISES (suspension / low insulin)

err = actual_bg_2h - pred_f   (positive = actual higher than predicted)

This script re-fetches the per-sample data from each site (using the same
pipeline as ns_multisite_backtest.py) and splits by direction + BG zone.
"""

import json, time, re, sys, warnings
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

OUT_DIR = Path.home() / 'Downloads'

# ── ISF formulas (same as backtest) ──
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

# ── API fetching (same as backtest) ──
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
            if dates:
                oldest = min(dates)
            else:
                break
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
        print(f" SKIP ({len(raw_entries) if raw_entries else 0} entries)")
        return None

    raw_tx = fetch_all_paginated(base_url, 'treatments.json', token, date_field='created_at')

    # Parse loop cycles
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

    # Bolus times
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

    # Forward BG
    cycle_epochs = np.array([int(t.timestamp()) for t in df['ts']])
    actual_2h = np.full(len(df), np.nan)
    bolus_age = np.full(len(df), np.nan)
    for i, t_s in enumerate(cycle_epochs):
        actual_2h[i] = get_cgm_at(t_s + 7200)
        bolus_age[i] = mins_since_bolus(t_s)
    df['actual_bg_2h'] = actual_2h
    df['bolus_age_min'] = bolus_age
    df['hour'] = df['ts'].dt.hour

    # TDD 7-day
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

    # Compute formula ISFs
    bg = strict['bg'].values
    isf_actual = strict['isf_actual'].values
    pred_drop = strict['pred_drop'].values
    actual_bg_2h = strict['actual_bg_2h'].values
    tdd_7d = strict['tdd_7day'].values
    tdd_median = np.median(tdd_7d)

    S_q = (1800.0 / tdd_median) / QUARTIC_AT_99
    S_db = (1800.0 / tdd_median) / FULL_DB_AT_99
    S_h = (1800.0 / tdd_median) / HYBRID_AT_105

    isf_q = isf_quartic(bg) * S_q
    isf_db = isf_full_diabeloop(bg) * S_db
    isf_h = isf_hybrid(bg) * S_h

    # Predicted 2h BG for each formula
    pred_loop = bg - pred_drop * (isf_actual / isf_actual)  # = pred_iob_24
    pred_q    = bg - pred_drop * (isf_q / isf_actual)
    pred_db   = bg - pred_drop * (isf_db / isf_actual)
    pred_h    = bg - pred_drop * (isf_h / isf_actual)

    print(f" OK ({len(strict)} samples)")

    return {
        'name': name, 'model': model, 'n': len(strict),
        'tdd_median': tdd_median,
        'bg': bg, 'pred_drop': pred_drop, 'actual_bg_2h': actual_bg_2h,
        'pred_loop': pred_loop, 'pred_q': pred_q, 'pred_db': pred_db, 'pred_h': pred_h,
    }


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 80)
print("FALLING vs RISING ANALYSIS — <105 and ≥105")
print("=" * 80)

all_sites = []
for site in SITES:
    r = process_site(site)
    if r is not None:
        all_sites.append(r)

print(f"\n\nLoaded {len(all_sites)} sites, {sum(s['n'] for s in all_sites)} total samples\n")

# ── Analysis function ──
def analyse_segment(label, sites, bg_mask_fn, drop_mask_fn):
    """Aggregate stats across sites for a BG zone + direction combo."""
    total_n = 0
    err_loop_all, err_q_all, err_db_all, err_h_all = [], [], [], []
    pred_loop_all, pred_q_all, pred_db_all, pred_h_all = [], [], [], []
    actual_all = []

    for s in sites:
        bg_m = bg_mask_fn(s['bg'])
        dr_m = drop_mask_fn(s['pred_drop'])
        m = bg_m & dr_m
        n = m.sum()
        if n == 0: continue
        total_n += n
        actual_all.extend(s['actual_bg_2h'][m])
        pred_loop_all.extend(s['pred_loop'][m])
        pred_q_all.extend(s['pred_q'][m])
        pred_db_all.extend(s['pred_db'][m])
        pred_h_all.extend(s['pred_h'][m])

    if total_n == 0:
        return None

    actual = np.array(actual_all)
    preds = {
        'Loop (DynISF)': np.array(pred_loop_all),
        'Quartic':       np.array(pred_q_all),
        'Full Diabeloop': np.array(pred_db_all),
        'Hybrid':        np.array(pred_h_all),
    }

    return total_n, actual, preds


def print_segment(label, total_n, actual, preds):
    """Print results for one segment."""
    mean_actual = actual.mean()
    print(f"\n  {label}  (n={total_n}, mean actual 2h BG = {mean_actual:.0f})")
    print(f"  {'Formula':<18s} {'Pred 2h':>8s} {'Pred−Act':>9s} {'|Err|':>6s}  Interpretation")
    print(f"  {'─'*18} {'─'*8} {'─'*9} {'─'*6}  {'─'*40}")
    for fname, pred in preds.items():
        mean_pred = pred.mean()
        diff = mean_pred - mean_actual  # positive = pred higher than actual
        mae = np.abs(pred - actual).mean()
        if diff > 2:
            interp = "predicts HIGHER than actual (over-estimate)"
        elif diff < -2:
            interp = "predicts LOWER than actual (under-estimate)"
        else:
            interp = "≈ matches actual"
        print(f"  {fname:<18s} {mean_pred:8.1f} {diff:+9.1f} {mae:6.1f}  {interp}")


def print_per_site(sites, bg_mask_fn, drop_mask_fn, direction_label):
    """Per-site breakdown for a segment."""
    print(f"\n  Per-site ({direction_label}):")
    print(f"  {'Site':<14s} {'Model':<8s} {'N':>5s} {'Actual':>7s}  {'Loop':>7s} {'Quart':>7s} {'FullDB':>7s} {'Hybrid':>7s}")
    print(f"  {'─'*14} {'─'*8} {'─'*5} {'─'*7}  {'─'*7} {'─'*7} {'─'*7} {'─'*7}")
    for s in sites:
        bg_m = bg_mask_fn(s['bg'])
        dr_m = drop_mask_fn(s['pred_drop'])
        m = bg_m & dr_m
        n = m.sum()
        if n < 3: continue
        act = s['actual_bg_2h'][m].mean()
        pl = s['pred_loop'][m].mean() - act
        pq = s['pred_q'][m].mean() - act
        pd_ = s['pred_db'][m].mean() - act
        ph = s['pred_h'][m].mean() - act
        print(f"  {s['name']:<14s} {s['model']:<8s} {n:5d} {act:7.0f}  {pl:+7.1f} {pq:+7.1f} {pd_:+7.1f} {ph:+7.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# BELOW 105
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 80)
print("BELOW 105 — GLUCOSE FALLING (pred_drop > 0, insulin active)")
print("  Loop predicts glucose will DROP. Actual 2h BG shows what happened.")
print("  Pred > Actual = over-estimated the drop (predicted too low)")
print("  Pred < Actual = under-estimated the drop (predicted too high)")
print("=" * 80)

r = analyse_segment("< 105 FALLING", all_sites,
                     lambda bg: bg < 105, lambda pd: pd > 0)
if r: print_segment("< 105, FALLING (insulin active)", *r)

print_per_site(all_sites, lambda bg: bg < 105, lambda pd: pd > 0, "<105 falling")

print("\n" + "=" * 80)
print("BELOW 105 — GLUCOSE RISING (pred_drop < 0, suspended/low insulin)")
print("  Loop predicts glucose will RISE. Actual 2h BG shows what happened.")
print("  Pred > Actual = over-estimated the rise (predicted too high)")
print("  Pred < Actual = under-estimated the rise (predicted too low)")
print("=" * 80)

r = analyse_segment("< 105 RISING", all_sites,
                     lambda bg: bg < 105, lambda pd: pd < 0)
if r: print_segment("< 105, RISING (suspended/low insulin)", *r)

print_per_site(all_sites, lambda bg: bg < 105, lambda pd: pd < 0, "<105 rising")


# ══════════════════════════════════════════════════════════════════════════════
# AT OR ABOVE 105
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("AT OR ABOVE 105 — GLUCOSE FALLING (pred_drop > 0, insulin active)")
print("=" * 80)

r = analyse_segment("≥ 105 FALLING", all_sites,
                     lambda bg: bg >= 105, lambda pd: pd > 0)
if r: print_segment("≥ 105, FALLING (insulin active)", *r)

print_per_site(all_sites, lambda bg: bg >= 105, lambda pd: pd > 0, "≥105 falling")

print("\n" + "=" * 80)
print("AT OR ABOVE 105 — GLUCOSE RISING (pred_drop < 0, suspended/low insulin)")
print("=" * 80)

r = analyse_segment("≥ 105 RISING", all_sites,
                     lambda bg: bg >= 105, lambda pd: pd < 0)
if r: print_segment("≥ 105, RISING (suspended/low insulin)", *r)

print_per_site(all_sites, lambda bg: bg >= 105, lambda pd: pd < 0, "≥105 rising")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 80)
print("SUMMARY: Mean (Predicted − Actual) 2h BG by zone and direction")
print("  Positive = formula predicts glucose will be HIGHER than it actually is")
print("  Negative = formula predicts glucose will be LOWER than it actually is")
print("=" * 80)

combos = [
    ("<105 falling",  lambda bg: bg < 105,  lambda pd: pd > 0),
    ("<105 rising",   lambda bg: bg < 105,  lambda pd: pd < 0),
    ("≥105 falling",  lambda bg: bg >= 105, lambda pd: pd > 0),
    ("≥105 rising",   lambda bg: bg >= 105, lambda pd: pd < 0),
]

print(f"\n  {'Zone + Direction':<18s} {'N':>6s}  {'Loop':>8s} {'Quartic':>8s} {'FullDB':>8s} {'Hybrid':>8s}")
print(f"  {'─'*18} {'─'*6}  {'─'*8} {'─'*8} {'─'*8} {'─'*8}")

for label, bg_fn, dr_fn in combos:
    r = analyse_segment(label, all_sites, bg_fn, dr_fn)
    if r is None:
        print(f"  {label:<18s} {'—':>6s}")
        continue
    n, actual, preds = r
    mean_act = actual.mean()
    vals = []
    for fname in ['Loop (DynISF)', 'Quartic', 'Full Diabeloop', 'Hybrid']:
        diff = preds[fname].mean() - mean_act
        vals.append(f"{diff:+8.1f}")
    print(f"  {label:<18s} {n:6d}  {vals[0]} {vals[1]} {vals[2]} {vals[3]}")

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
