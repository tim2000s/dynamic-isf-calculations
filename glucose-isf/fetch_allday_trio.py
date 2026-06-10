#!/usr/bin/env python3
"""
Rebuild Trio multisite cache with ALL HOURS (not just overnight).
Filter: COB=0, bolus_age >= 120min, BG 72-200, any hour.
Stores per-site numpy arrays + hour data.
"""

import json, time, re, warnings, pickle, sys
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

OUT_DIR = Path.home() / 'Downloads' / '4 Hour analysis'
CACHE_FILE = OUT_DIR / 'speculation' / 'multisite_allday_cache.pkl'

# ISF formulas
def isf_quartic(bg):
    G = np.asarray(bg, dtype=float)
    return 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4

QUARTIC_AT_99 = float(isf_quartic(99.0))

# Nightscout sites + tokens are loaded from ns_sites.json (gitignored — never commit credentials).
import json as _json, os as _os
SITES = _json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "ns_sites.json")))

ANON = {
    'henny425': 'User-A', 'aadiabetes': 'User-B', 'diajesse': 'User-C',
    'svns': 'User-D', 'fuxchr': 'User-E', 'mikens': 'User-F',
    'andycgm': 'User-G', 'noahr': 'User-H', 'nightscout1': 'User-I',
    'eli': 'User-J', 'ns_rot6': 'User-K', 'kelseyhuss': 'User-L',
}


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


def process_site(site):
    name, base_url, token, model = site['name'], site['url'], site['token'], site['model']
    anon_name = ANON.get(name, name)
    print(f"\n  Fetching {anon_name} ({name})...", end='', flush=True)

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
            if len(pred_iob) < 12: continue
            tdd_val = sg.get('TDD')
            if tdd_val is None or tdd_val <= 0:
                tdd_val = parse_tdd_from_reason(reason)

            cycles.append({
                'ts': ts, 'bg': bg_val, 'isf_actual': float(isf_val),
                'cob': sg.get('COB'), 'iob': sg.get('IOB'),
                'tdd': float(tdd_val) if tdd_val else None,
                'pred_iob_final': float(pred_iob[-1]),
                'pred_horizon_s': (len(pred_iob) - 1) * 5 * 60,
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

    bolus_epochs_list = []
    for t in (raw_tx or []):
        try:
            ins = t.get('insulin')
            if ins and float(ins) > 0:
                ts = pd.to_datetime(t['created_at'], utc=True)
                bolus_epochs_list.append(int(ts.timestamp()))
        except Exception:
            continue
    bolus_arr = np.array(sorted(bolus_epochs_list)) if bolus_epochs_list else np.array([])

    def mins_since_bolus(target_ts_s):
        if len(bolus_arr) == 0: return 9999.0
        idx = np.searchsorted(bolus_arr, target_ts_s, side='right') - 1
        return (target_ts_s - bolus_arr[idx]) / 60.0 if idx >= 0 else 9999.0

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

    # ── ALL-DAY FASTING FILTER (COB=0, bolus_age >= 120, any hour) ──
    mask_fasting = (
        (df['cob'].fillna(99) == 0) &
        (df['bg'] >= 72) & (df['bg'] <= 200) &
        (df['bolus_age_min'] >= 120) &
        (df['tdd'].notna()) & (df['tdd'] > 0)
    )
    fasting = df[mask_fasting].copy()
    fasting = fasting.dropna(subset=['pred_iob_final', 'actual_bg_end']).copy()
    if len(fasting) < 10:
        print(f" SKIP ({len(fasting)} fasting samples)")
        return None

    fasting['pred_drop'] = fasting['bg'] - fasting['pred_iob_final']
    fasting['actual_drop'] = fasting['bg'] - fasting['actual_bg_end']

    strict = fasting[
        (np.abs(fasting['pred_drop']) > 3) &
        ((fasting['actual_drop'] / fasting['pred_drop']).between(0, 5))
    ].copy().sort_values('ts').reset_index(drop=True)

    if len(strict) < 10:
        print(f" SKIP ({len(strict)} quality)")
        return None

    # Split into subsets
    overnight = strict[strict['hour'] < 8].copy().reset_index(drop=True)
    daytime = strict[strict['hour'] >= 8].copy().reset_index(drop=True)

    tdd_7d = strict['tdd_7day'].values
    tdd_median = np.median(tdd_7d[~np.isnan(tdd_7d)]) if np.any(~np.isnan(tdd_7d)) else np.median(strict['tdd'].values)

    def extract_arrays(subset):
        ts_epochs = np.array([int(t.timestamp()) for t in subset['ts']], dtype=np.int64)
        tdd_vals = subset['tdd'].values.astype(float) if 'tdd' in subset.columns else np.full(len(subset), np.nan)
        return {
            'bg': subset['bg'].values.astype(float),
            'isf_actual': subset['isf_actual'].values.astype(float),
            'pred_drop': subset['pred_drop'].values.astype(float),
            'actual_bg_end': subset['actual_bg_end'].values.astype(float),
            'hour': subset['hour'].values.astype(int),
            'ts_epoch': ts_epochs,
            'tdd': tdd_vals,
        }

    n_overnight = len(overnight)
    n_daytime = len(daytime)
    n_allday = len(strict)

    print(f" OK (allday={n_allday}, overnight={n_overnight}, daytime={n_daytime}, "
          f"TDD={tdd_median:.1f})")

    # Compute pred_loop for each subset
    def make_pred_loop(arrays):
        return arrays['bg'] - arrays['pred_drop']

    allday_arr = extract_arrays(strict)
    overnight_arr = extract_arrays(overnight) if n_overnight >= 10 else None
    daytime_arr = extract_arrays(daytime) if n_daytime >= 10 else None

    result = {
        'name': anon_name, 'original_name': name, 'model': model,
        'tdd_median': tdd_median,
        'allday': {**allday_arr, 'n': n_allday, 'pred_loop': make_pred_loop(allday_arr)},
    }
    if overnight_arr and n_overnight >= 10:
        result['overnight'] = {**overnight_arr, 'n': n_overnight, 'pred_loop': make_pred_loop(overnight_arr)}
    if daytime_arr and n_daytime >= 10:
        result['daytime'] = {**daytime_arr, 'n': n_daytime, 'pred_loop': make_pred_loop(daytime_arr)}

    return result


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

if CACHE_FILE.exists():
    print(f"Loading cached data from {CACHE_FILE}")
    with open(CACHE_FILE, 'rb') as f:
        all_sites = pickle.load(f)
    for s in all_sites:
        ad = s['allday']
        on = s.get('overnight', {})
        dt = s.get('daytime', {})
        print(f"  {s['name']:8s}  allday={ad['n']:5d}  "
              f"overnight={on.get('n', 0):5d}  daytime={dt.get('n', 0):5d}")
else:
    print("Fetching all-day fasting data from all sites...")
    print("(COB=0, bolus_age >= 120min, BG 72-200, ALL hours)")
    print("=" * 70)

    all_sites = []
    for site in SITES:
        r = process_site(site)
        if r is not None:
            all_sites.append(r)

    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(all_sites, f)
    print(f"\nCached to {CACHE_FILE}")

print(f"\n{len(all_sites)} sites loaded")
total_allday = sum(s['allday']['n'] for s in all_sites)
total_overnight = sum(s.get('overnight', {}).get('n', 0) for s in all_sites)
total_daytime = sum(s.get('daytime', {}).get('n', 0) for s in all_sites)
print(f"Total samples: allday={total_allday}, overnight={total_overnight}, daytime={total_daytime}")
print("DONE")
