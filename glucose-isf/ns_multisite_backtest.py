#!/usr/bin/env python3
"""
Multi-Site Nightscout Backtest — Diabeloop & Hybrid ISF Evaluation
===================================================================

Fetches data from multiple Nightscout (Trio) sites and evaluates:
  A: Loop actual ISF (the site's own sigmoid/log dynamic ISF)
  C: Quartic polynomial (Diabeloop >100, extended) — TDD-scaled
  D: Full Diabeloop (quadratic ≤100, quartic >100) — TDD-scaled
  E: Hybrid (quartic ≥105, power-law k=3.5 <105) — TDD-scaled

Each site's ISF predictions are compared against actual 2-hour glucose outcomes.
"""

import json
import warnings
import time
import sys
from urllib.parse import urlparse, parse_qs, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

warnings.filterwarnings('ignore')

OUT_DIR = Path.home() / 'Downloads'

# ── Site definitions ──
SITES = [
    {'name': 'henny425',    'model': 'sigmoid', 'url': '***REDACTED-URL***',               'token': '***REDACTED-TOKEN***'},
    {'name': 'aadiabetes',  'model': 'sigmoid', 'url': '***REDACTED-URL***',     'token': '***REDACTED-TOKEN***'},
    {'name': 'diajesse',    'model': 'sigmoid', 'url': '***REDACTED-URL***',       'token': None},
    {'name': 'svns',        'model': 'sigmoid', 'url': '***REDACTED-URL***',  'token': None},
    {'name': 'andycgm',     'model': 'log',     'url': '***REDACTED-URL***',     'token': '***REDACTED-TOKEN***'},
    {'name': 'noahr',       'model': 'log',     'url': '***REDACTED-URL***',             'token': None},
    {'name': 'fuxchr',      'model': 'sigmoid', 'url': '***REDACTED-URL***',               'token': '***REDACTED-TOKEN***'},
    {'name': 'taylor',      'model': 'sigmoid', 'url': '***REDACTED-URL***',      'token': None},
    {'name': 'nightscout1', 'model': 'log',     'url': '***REDACTED-URL***',     'token': '***REDACTED-TOKEN***'},
    {'name': 'eli',         'model': 'log',     'url': '***REDACTED-URL***',            'token': '***REDACTED-TOKEN***'},
    {'name': 'mikens',      'model': 'sigmoid', 'url': '***REDACTED-URL***',          'token': '***REDACTED-TOKEN***'},
    {'name': 'ns_rot6',     'model': 'log',     'url': '***REDACTED-URL***',                      'token': '***REDACTED-TOKEN***'},
    {'name': 'kelseyhuss',  'model': 'log',     'url': '***REDACTED-URL***',  'token': '***REDACTED-TOKEN***'},
]


# ══════════════════════════════════════════════════════════════════════════════
# ISF FORMULAS
# ══════════════════════════════════════════════════════════════════════════════

def isf_quartic(bg):
    """Diabeloop quartic (designed for >100, extended to all BG)"""
    G = np.asarray(bg, dtype=float)
    return 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4


def isf_full_diabeloop(bg):
    """Full Diabeloop: quadratic ≤100, quartic >100"""
    G = np.asarray(bg, dtype=float)
    quad = 98.03 - 1.077*G + 0.008868*G**2
    quart = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    return np.where(G <= 100, quad, quart)


def isf_hybrid(bg):
    """Hybrid: quartic ≥105, power-law <105"""
    G = np.asarray(bg, dtype=float)
    quart = 272.0 - 3.121*G + 0.01511*G**2 - 3.305e-05*G**3 + 2.69e-08*G**4
    power = 75.8 * (105.0 / G) ** 3.5
    return np.where(G >= 105, quart, power)


# Anchor points for TDD-based scaling
QUARTIC_AT_TARGET = float(isf_quartic(99.0))    # ~81.6
FULL_DB_AT_TARGET = float(isf_full_diabeloop(99.0))  # ~81.6 (quartic at 99)
HYBRID_AT_105 = 75.8


# ══════════════════════════════════════════════════════════════════════════════
# NIGHTSCOUT API FETCHING
# ══════════════════════════════════════════════════════════════════════════════

def ns_fetch(base_url, endpoint, token=None, params=None, max_retries=3):
    """Fetch from Nightscout API with pagination support."""
    if params is None:
        params = {}
    if token:
        params['token'] = token

    url = f"{base_url}/api/v1/{endpoint}?{urlencode(params, safe='[]$')}"
    for attempt in range(max_retries):
        try:
            req = Request(url, headers={
                'Accept': 'application/json',
                'User-Agent': 'Mozilla/5.0 (Nightscout-Backtest)',
            })
            with urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except HTTPError as e:
            if e.code in (401, 403):
                print(f"    AUTH FAIL ({e.code})")
                return None
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"    HTTP {e.code} after {max_retries} tries")
                return None
        except (URLError, TimeoutError) as e:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"    TIMEOUT/ERROR: {e}")
                return None


def fetch_all_paginated(base_url, endpoint, token=None, max_records=200000,
                        months_back=6, date_field='created_at'):
    """Fetch records using month-by-month date windows, then paginate within each.
    date_field: 'created_at' for devicestatus/treatments, 'dateString' for entries."""
    from datetime import datetime, timedelta

    all_records = []
    seen_ids = set()
    page_size = 10000

    now = datetime.utcnow()
    windows = []
    for m in range(months_back):
        end = now - timedelta(days=30 * m)
        start = now - timedelta(days=30 * (m + 1))
        windows.append((start.strftime('%Y-%m-%dT%H:%M:%S.000Z'),
                         end.strftime('%Y-%m-%dT%H:%M:%S.000Z')))

    for win_start, win_end in windows:
        oldest_in_window = win_end

        while True:
            params = {
                'count': page_size,
                f'find[{date_field}][$gte]': win_start,
                f'find[{date_field}][$lt]': oldest_in_window,
            }
            batch = ns_fetch(base_url, endpoint, token, params)
            if batch is None or len(batch) == 0:
                break

            new_records = []
            for r in batch:
                rid = r.get('_id')
                if rid and rid not in seen_ids:
                    seen_ids.add(rid)
                    new_records.append(r)

            if not new_records:
                break

            all_records.extend(new_records)

            dates = []
            for r in new_records:
                ca = r.get('created_at') or r.get('dateString')
                if ca:
                    dates.append(ca)
            if not dates:
                break

            new_oldest = min(dates)
            if new_oldest == oldest_in_window:
                break
            oldest_in_window = new_oldest

            if len(batch) < page_size * 0.9:
                break

            if len(all_records) >= max_records:
                break

        sys.stdout.write(f"\r    Fetched {len(all_records):,} records (back to {win_start[:10]})...")
        sys.stdout.flush()

        if len(all_records) >= max_records:
            break

    if all_records:
        sys.stdout.write(f"\r    Fetched {len(all_records):,} records total                              \n")
        sys.stdout.flush()

    return all_records


# ══════════════════════════════════════════════════════════════════════════════
# SITE DATA PROCESSING
# ══════════════════════════════════════════════════════════════════════════════

import re

def parse_isf_from_reason(reason, bg_mgdl):
    """Parse dynamic ISF from reason string. Returns ISF in mg/dL or None.
    Format: 'ISF: 110→108' (mg/dL) or 'ISF: 5→4.7' (mmol/L).
    The second number (after →) is the dynamic ISF."""
    m = re.search(r'ISF:\s*([\d.]+)\s*(?:→|->)\s*([\d.]+)', reason)
    if not m:
        return None
    dynamic_isf = float(m.group(2))
    # Detect mmol/L: check if Target in reason is in mmol range
    target_m = re.search(r'Target:\s*([\d.]+)', reason)
    if target_m and float(target_m.group(1)) < 30:
        # Target is in mmol/L, so ISF is too
        dynamic_isf *= 18.0
    return dynamic_isf


def parse_tdd_from_reason(reason):
    """Parse TDD from reason string. Format: 'TDD: 45.2 U'"""
    m = re.search(r'TDD:\s*([\d.]+)', reason)
    return float(m.group(1)) if m else None


def process_site(site):
    """Fetch data from a site and run the backtest. Returns results dict or None."""
    name = site['name']
    base_url = site['url']
    token = site['token']
    model = site['model']

    print(f"\n{'═'*70}")
    print(f"  {name} ({model})")
    print(f"{'═'*70}")

    # ── Quick connectivity check ──
    test = ns_fetch(base_url, 'devicestatus.json', token, {'count': 1})
    if test is None:
        print(f"  SKIPPING — cannot connect")
        return None

    # ── Fetch devicestatus ──
    print(f"  Fetching devicestatus...")
    raw_ds = fetch_all_paginated(base_url, 'devicestatus.json', token,
                                  date_field='created_at')
    if not raw_ds or len(raw_ds) < 100:
        print(f"  SKIPPING — only {len(raw_ds) if raw_ds else 0} devicestatus records")
        return None

    # ── Fetch CGM entries ──
    print(f"  Fetching entries...")
    raw_entries = fetch_all_paginated(base_url, 'entries.json', token,
                                       date_field='dateString')
    if not raw_entries or len(raw_entries) < 100:
        print(f"  SKIPPING — only {len(raw_entries) if raw_entries else 0} entries")
        return None

    # ── Fetch treatments ──
    print(f"  Fetching treatments...")
    raw_tx = fetch_all_paginated(base_url, 'treatments.json', token,
                                  date_field='created_at')

    # ── Parse loop cycles (Trio format) ──
    print(f"  Parsing loop cycles...")
    cycles = []
    isf_from_reason_count = 0
    tdd_from_reason_count = 0
    for r in raw_ds:
        try:
            ts = pd.to_datetime(r['created_at'], utc=True)
            sg = r.get('openaps', {}).get('suggested', {})
            if not sg or 'bg' not in sg:
                continue

            bg_val = sg['bg']
            isf_val = sg.get('ISF')
            reason = sg.get('reason', '')

            # Fallback: parse ISF from reason string if ISF=0 or missing
            if not isf_val or isf_val <= 0:
                isf_val = parse_isf_from_reason(reason, bg_val)
                if isf_val:
                    isf_from_reason_count += 1
            if not isf_val or isf_val <= 0:
                continue

            pred_iob = sg.get('predBGs', {}).get('IOB') or []
            tdd_val = sg.get('TDD')

            # Fallback: parse TDD from reason string
            if tdd_val is None or tdd_val <= 0:
                tdd_val = parse_tdd_from_reason(reason)
                if tdd_val:
                    tdd_from_reason_count += 1

            cycles.append({
                'ts': ts,
                'bg': bg_val,
                'isf_actual': float(isf_val),  # Trio's dynamic ISF in mg/dL
                'sensitivity_ratio': sg.get('sensitivityRatio'),
                'cob': sg.get('COB'),
                'iob': sg.get('IOB'),
                'tdd': float(tdd_val) if tdd_val else None,
                'pred_iob_24': pred_iob[24] if len(pred_iob) > 24 else None,
            })
        except Exception:
            continue

    if isf_from_reason_count:
        print(f"  (parsed ISF from reason string for {isf_from_reason_count:,} cycles)")
    if tdd_from_reason_count:
        print(f"  (parsed TDD from reason string for {tdd_from_reason_count:,} cycles)")

    if len(cycles) < 100:
        print(f"  SKIPPING — only {len(cycles)} valid loop cycles")
        return None

    df = pd.DataFrame(cycles).sort_values('ts').reset_index(drop=True)
    print(f"  Parsed: {len(df):,} cycles  ({df['ts'].min().date()} → {df['ts'].max().date()})")

    # ── Parse CGM ──
    cgm_rows = []
    for e in raw_entries:
        try:
            sgv = e.get('sgv')
            if not sgv or not (40 <= sgv <= 400):
                continue
            ca = e.get('created_at') or e.get('dateString')
            ts = pd.to_datetime(ca, utc=True) if ca else pd.to_datetime(int(e['date']), unit='ms', utc=True)
            cgm_rows.append({'ts': ts, 'sgv': float(sgv)})
        except Exception:
            continue

    df_cgm = pd.DataFrame(cgm_rows).sort_values('ts').reset_index(drop=True)
    cgm_epochs = np.array([int(t.timestamp()) for t in df_cgm['ts']])
    cgm_sgv = df_cgm['sgv'].values
    print(f"  CGM: {len(df_cgm):,} readings")

    def get_cgm_at(target_ts_s, tolerance_s=150):
        idx = np.searchsorted(cgm_epochs, target_ts_s)
        best_val, best_diff = np.nan, np.inf
        for k in (idx - 1, idx):
            if 0 <= k < len(cgm_epochs):
                diff = abs(cgm_epochs[k] - target_ts_s)
                if diff < best_diff:
                    best_diff, best_val = diff, cgm_sgv[k]
        return best_val if best_diff <= tolerance_s else np.nan

    # ── Parse bolus times ──
    bolus_epochs = []
    for t in (raw_tx or []):
        try:
            et = t.get('eventType', '')
            ins = t.get('insulin')
            if ins and float(ins) > 0:
                ts = pd.to_datetime(t['created_at'], utc=True)
                bolus_epochs.append(int(ts.timestamp()))
        except Exception:
            continue
    bolus_epochs = np.array(sorted(bolus_epochs)) if bolus_epochs else np.array([])

    def mins_since_bolus(target_ts_s):
        if len(bolus_epochs) == 0:
            return 9999.0
        idx = np.searchsorted(bolus_epochs, target_ts_s, side='right') - 1
        return (target_ts_s - bolus_epochs[idx]) / 60.0 if idx >= 0 else 9999.0

    # ── Forward BG lookups ──
    print(f"  Computing 2h forward BG...")
    cycle_epochs = np.array([int(t.timestamp()) for t in df['ts']])
    actual_2h = np.full(len(df), np.nan)
    bolus_age = np.full(len(df), np.nan)

    for i, t_s in enumerate(cycle_epochs):
        actual_2h[i] = get_cgm_at(t_s + 7200)
        bolus_age[i] = mins_since_bolus(t_s)

    df['actual_bg_2h'] = actual_2h
    df['bolus_age_min'] = bolus_age
    df['hour'] = df['ts'].dt.hour

    # ── Compute TDD from Trio's reported TDD (already in devicestatus) ──
    # Use 7-day rolling median of Trio's own TDD field
    df['date'] = df['ts'].dt.date
    tdd_valid = df[df['tdd'].notna() & (df['tdd'] > 0)]
    if len(tdd_valid) > 0:
        daily_tdd = tdd_valid.groupby('date')['tdd'].median().reset_index()
        daily_tdd.columns = ['date', 'tdd_daily']
        daily_tdd['tdd_7day'] = daily_tdd['tdd_daily'].rolling(7, min_periods=1).mean()
        df = df.merge(daily_tdd[['date', 'tdd_7day']], on='date', how='left')
    else:
        # Use per-cycle TDD directly if no daily aggregation possible
        df['tdd_7day'] = df['tdd']

    # ── Overnight fasting filter ──
    print(f"  Filtering overnight fasting...")
    mask = (
        (df['hour'] < 8) &
        (df['cob'].fillna(99) == 0) &
        (df['bg'] >= 72) & (df['bg'] <= 200) &
        (df['bolus_age_min'] >= 120) &
        (df['tdd'].notna()) & (df['tdd'] > 0)
    )
    on = df[mask].copy()
    print(f"  Overnight fasting: {len(on):,} cycles")

    if len(on) < 10:
        print(f"  SKIPPING — only {len(on)} overnight fasting cycles")
        return None

    on = on.dropna(subset=['pred_iob_24', 'actual_bg_2h']).copy()
    if len(on) < 10:
        print(f"  SKIPPING — only {len(on)} cycles with 2h forward BG")
        return None

    on['pred_drop'] = on['bg'] - on['pred_iob_24']
    on['actual_drop'] = on['bg'] - on['actual_bg_2h']

    # Quality filters (less strict than single-patient study)
    strict = on[
        (np.abs(on['pred_drop']) > 3) &
        ((on['actual_drop'] / on['pred_drop']).between(0, 5))
    ].copy()
    strict = strict.sort_values('ts').reset_index(drop=True)
    print(f"  After quality filtering: {len(strict):,} valid 2h samples")

    if len(strict) < 10:
        print(f"  SKIPPING — only {len(strict)} quality samples")
        return None

    # ── Extract arrays ──
    bg = strict['bg'].values
    isf_actual = strict['isf_actual'].values
    pred_drop = strict['pred_drop'].values
    actual_bg_2h_arr = strict['actual_bg_2h'].values
    tdd_7d = strict['tdd_7day'].values

    tdd_median = np.median(tdd_7d)

    def compute_errors(isf_formula, isf_true):
        pred_f = bg - pred_drop * (isf_formula / isf_true)
        err = actual_bg_2h_arr - pred_f
        valid = ~np.isnan(err)
        e = err[valid]
        if len(e) == 0:
            return np.nan, np.nan, np.nan
        mae = np.abs(e).mean()
        bias = e.mean()
        w18 = (np.abs(e) <= 18).mean() * 100
        return mae, bias, w18

    def compute_band_errors(isf_formula, isf_true):
        """MAE by glucose band."""
        pred_f = bg - pred_drop * (isf_formula / isf_true)
        err = actual_bg_2h_arr - pred_f
        bands = {}
        for bname, lo, hi in [('<90', 0, 90), ('90-105', 90, 105),
                               ('105-120', 105, 120), ('120-150', 120, 150),
                               ('150+', 150, 999)]:
            mask = (bg >= lo) & (bg < hi) & ~np.isnan(err)
            if mask.sum() >= 5:
                bands[bname] = {
                    'n': int(mask.sum()),
                    'mae': float(np.abs(err[mask]).mean()),
                    'bias': float(err[mask].mean()),
                }
        return bands

    # ── Formula A: Loop actual (baseline) ──
    mae_loop, bias_loop, w18_loop = compute_errors(isf_actual, isf_actual)

    # ── TDD-scaled formulas ──
    # Scale factor: S = (1800 / TDD_median) / anchor
    S_quartic = (1800.0 / tdd_median) / QUARTIC_AT_TARGET
    S_full_db = (1800.0 / tdd_median) / FULL_DB_AT_TARGET
    S_hybrid = (1800.0 / tdd_median) / HYBRID_AT_105

    isf_C = isf_quartic(bg) * S_quartic
    isf_D = isf_full_diabeloop(bg) * S_full_db
    isf_E = isf_hybrid(bg) * S_hybrid

    mae_C, bias_C, w18_C = compute_errors(isf_C, isf_actual)
    mae_D, bias_D, w18_D = compute_errors(isf_D, isf_actual)
    mae_E, bias_E, w18_E = compute_errors(isf_E, isf_actual)

    # Band errors
    bands_loop = compute_band_errors(isf_actual, isf_actual)
    bands_C = compute_band_errors(isf_C, isf_actual)
    bands_D = compute_band_errors(isf_D, isf_actual)
    bands_E = compute_band_errors(isf_E, isf_actual)

    # ── Mean predicted 2h BG below 90 ──
    mask_sub90 = bg < 90
    sub90_actual_mean = float(actual_bg_2h_arr[mask_sub90].mean()) if mask_sub90.sum() >= 5 else np.nan

    def mean_pred_2h_sub90(isf_formula, isf_true):
        if mask_sub90.sum() < 5:
            return np.nan
        pred_f = bg[mask_sub90] - pred_drop[mask_sub90] * (isf_formula[mask_sub90] / isf_true[mask_sub90])
        return float(np.nanmean(pred_f))

    pred_2h_loop_sub90 = mean_pred_2h_sub90(isf_actual, isf_actual)
    pred_2h_C_sub90 = mean_pred_2h_sub90(isf_C, isf_actual)
    pred_2h_D_sub90 = mean_pred_2h_sub90(isf_D, isf_actual)
    pred_2h_E_sub90 = mean_pred_2h_sub90(isf_E, isf_actual)

    result = {
        'name': name,
        'model': model,
        'n_cycles': len(df),
        'n_samples': len(strict),
        'date_range': f"{df['ts'].min().date()} → {df['ts'].max().date()}",
        'tdd_median': tdd_median,
        'S_quartic': S_quartic,
        'S_hybrid': S_hybrid,
        'overall': {
            'loop':    {'mae': mae_loop, 'bias': bias_loop, 'w18': w18_loop},
            'quartic': {'mae': mae_C, 'bias': bias_C, 'w18': w18_C},
            'full_db': {'mae': mae_D, 'bias': bias_D, 'w18': w18_D},
            'hybrid':  {'mae': mae_E, 'bias': bias_E, 'w18': w18_E},
        },
        'bands': {
            'loop': bands_loop, 'quartic': bands_C,
            'full_db': bands_D, 'hybrid': bands_E,
        },
        'sub90': {
            'n': int(mask_sub90.sum()),
            'actual_mean': sub90_actual_mean,
            'pred_loop': pred_2h_loop_sub90,
            'pred_quartic': pred_2h_C_sub90,
            'pred_full_db': pred_2h_D_sub90,
            'pred_hybrid': pred_2h_E_sub90,
        },
        'bg_mean': float(bg.mean()),
        'isf_actual_median': float(np.median(isf_actual)),
    }

    # Print summary
    print(f"\n  ── Results for {name} ({model}) ──")
    print(f"  Samples: {len(strict):,}  |  TDD median: {tdd_median:.1f} U/day  |  BG mean: {bg.mean():.0f}")
    print(f"  {'Formula':<20s}  {'MAE':>6s}  {'Bias':>7s}  {'W18%':>5s}")
    print(f"  {'─'*20}  {'─'*6}  {'─'*7}  {'─'*5}")
    print(f"  {'A: Loop actual':<20s}  {mae_loop:6.1f}  {bias_loop:+7.1f}  {w18_loop:5.1f}")
    print(f"  {'C: Quartic':<20s}  {mae_C:6.1f}  {bias_C:+7.1f}  {w18_C:5.1f}")
    print(f"  {'D: Full Diabeloop':<20s}  {mae_D:6.1f}  {bias_D:+7.1f}  {w18_D:5.1f}")
    print(f"  {'E: Hybrid':<20s}  {mae_E:6.1f}  {bias_E:+7.1f}  {w18_E:5.1f}")

    if bands_loop:
        print(f"\n  {'Band':<10s}  {'n':>4s}  {'Loop':>6s}  {'Quart':>6s}  {'FullDB':>6s}  {'Hybrid':>6s}")
        print(f"  {'─'*10}  {'─'*4}  {'─'*6}  {'─'*6}  {'─'*6}  {'─'*6}")
        for bname in ['<90', '90-105', '105-120', '120-150', '150+']:
            bl = bands_loop.get(bname, {})
            bc = bands_C.get(bname, {})
            bd = bands_D.get(bname, {})
            be = bands_E.get(bname, {})
            if bl:
                print(f"  {bname:<10s}  {bl['n']:4d}  {bl['mae']:6.1f}  "
                      f"{bc.get('mae', np.nan):6.1f}  {bd.get('mae', np.nan):6.1f}  "
                      f"{be.get('mae', np.nan):6.1f}")

    return result


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 70)
    print("MULTI-SITE NIGHTSCOUT BACKTEST")
    print("Diabeloop & Hybrid ISF Evaluation")
    print("=" * 70)

    results = []
    for site in SITES:
        try:
            r = process_site(site)
            if r:
                results.append(r)
        except Exception as e:
            print(f"\n  ERROR processing {site['name']}: {e}")
            import traceback
            traceback.print_exc()

    if not results:
        print("\nNo sites returned valid results.")
        sys.exit(1)

    # ══════════════════════════════════════════════════════════════════════════
    # AGGREGATE SUMMARY
    # ══════════════════════════════════════════════════════════════════════════
    print(f"\n\n{'═'*90}")
    print("AGGREGATE SUMMARY ACROSS ALL SITES")
    print(f"{'═'*90}")
    print(f"\nSites with valid data: {len(results)}/{len(SITES)}")
    print(f"Total valid samples: {sum(r['n_samples'] for r in results):,}")

    # Per-site summary table
    print(f"\n{'Site':<14s} {'Model':<8s} {'N':>5s} {'TDD':>5s}  "
          f"{'Loop':>6s} {'Quart':>6s} {'FullDB':>6s} {'Hybrid':>6s}  "
          f"{'Best':>8s}")
    print(f"{'─'*14} {'─'*8} {'─'*5} {'─'*5}  "
          f"{'─'*6} {'─'*6} {'─'*6} {'─'*6}  "
          f"{'─'*8}")

    for r in results:
        o = r['overall']
        maes = {'Loop': o['loop']['mae'], 'Quart': o['quartic']['mae'],
                'FullDB': o['full_db']['mae'], 'Hybrid': o['hybrid']['mae']}
        best = min(maes, key=maes.get)
        print(f"{r['name']:<14s} {r['model']:<8s} {r['n_samples']:5d} {r['tdd_median']:5.1f}  "
              f"{o['loop']['mae']:6.1f} {o['quartic']['mae']:6.1f} "
              f"{o['full_db']['mae']:6.1f} {o['hybrid']['mae']:6.1f}  "
              f"{best:<8s}")

    # Weighted average (by sample count)
    total_n = sum(r['n_samples'] for r in results)
    def weighted_avg(formula_key, metric):
        return sum(r['overall'][formula_key][metric] * r['n_samples'] for r in results) / total_n

    print(f"\n{'Weighted avg':<14s} {'':8s} {total_n:5d} {'':5s}  "
          f"{weighted_avg('loop', 'mae'):6.1f} {weighted_avg('quartic', 'mae'):6.1f} "
          f"{weighted_avg('full_db', 'mae'):6.1f} {weighted_avg('hybrid', 'mae'):6.1f}")

    # Below-90 summary
    sites_with_sub90 = [r for r in results if r['sub90']['n'] >= 5]
    if sites_with_sub90:
        print(f"\n── Below 90 mg/dL Summary ──")
        print(f"{'Site':<14s} {'N':>4s}  {'Actual':>6s}  {'Loop':>6s} {'Quart':>6s} {'FullDB':>6s} {'Hybrid':>6s}")
        print(f"{'─'*14} {'─'*4}  {'─'*6}  {'─'*6} {'─'*6} {'─'*6} {'─'*6}")
        for r in sites_with_sub90:
            s = r['sub90']
            print(f"{r['name']:<14s} {s['n']:4d}  {s['actual_mean']:6.0f}  "
                  f"{s['pred_loop']:6.0f} {s['pred_quartic']:6.0f} "
                  f"{s['pred_full_db']:6.0f} {s['pred_hybrid']:6.0f}")

    # ── Band-level aggregation ──
    print(f"\n── Glucose Band MAE (weighted average across sites) ──")
    print(f"{'Band':<10s}  {'N':>5s}  {'Loop':>6s}  {'Quart':>6s}  {'FullDB':>6s}  {'Hybrid':>6s}")
    print(f"{'─'*10}  {'─'*5}  {'─'*6}  {'─'*6}  {'─'*6}  {'─'*6}")
    for bname in ['<90', '90-105', '105-120', '120-150', '150+']:
        total_band_n = 0
        sums = {'loop': 0, 'quartic': 0, 'full_db': 0, 'hybrid': 0}
        for r in results:
            bl = r['bands']['loop'].get(bname, {})
            if not bl:
                continue
            n = bl['n']
            total_band_n += n
            for fkey in sums:
                b = r['bands'][fkey].get(bname, {})
                if b:
                    sums[fkey] += b['mae'] * n
        if total_band_n >= 5:
            print(f"{bname:<10s}  {total_band_n:5d}  "
                  f"{sums['loop']/total_band_n:6.1f}  "
                  f"{sums['quartic']/total_band_n:6.1f}  "
                  f"{sums['full_db']/total_band_n:6.1f}  "
                  f"{sums['hybrid']/total_band_n:6.1f}")

    # ── Save results to JSON ──
    out_json = OUT_DIR / 'multisite_backtest_results.json'
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to {out_json}")

    # ══════════════════════════════════════════════════════════════════════════
    # CHART
    # ══════════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Multi-Site Nightscout Backtest — Diabeloop & Hybrid ISF', fontsize=14, fontweight='bold')

    # Panel 1: Overall MAE by site
    ax = axes[0, 0]
    site_names = [r['name'] for r in results]
    x = np.arange(len(site_names))
    w = 0.2
    ax.bar(x - 1.5*w, [r['overall']['loop']['mae'] for r in results], w, label='Loop actual', color='#888888')
    ax.bar(x - 0.5*w, [r['overall']['quartic']['mae'] for r in results], w, label='Quartic', color='#e74c3c')
    ax.bar(x + 0.5*w, [r['overall']['full_db']['mae'] for r in results], w, label='Full Diabeloop', color='#e67e22')
    ax.bar(x + 1.5*w, [r['overall']['hybrid']['mae'] for r in results], w, label='Hybrid', color='#2ecc71')
    ax.set_xticks(x)
    ax.set_xticklabels(site_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('MAE (mg/dL)')
    ax.set_title('Overall MAE by Site')
    ax.legend(fontsize=8)

    # Panel 2: Overall bias by site
    ax = axes[0, 1]
    ax.bar(x - 1.5*w, [r['overall']['loop']['bias'] for r in results], w, label='Loop actual', color='#888888')
    ax.bar(x - 0.5*w, [r['overall']['quartic']['bias'] for r in results], w, label='Quartic', color='#e74c3c')
    ax.bar(x + 0.5*w, [r['overall']['full_db']['bias'] for r in results], w, label='Full Diabeloop', color='#e67e22')
    ax.bar(x + 1.5*w, [r['overall']['hybrid']['bias'] for r in results], w, label='Hybrid', color='#2ecc71')
    ax.set_xticks(x)
    ax.set_xticklabels(site_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Bias (mg/dL)')
    ax.set_title('Overall Bias by Site')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.legend(fontsize=8)

    # Panel 3: Weighted average MAE by glucose band
    ax = axes[1, 0]
    band_names = ['<90', '90-105', '105-120', '120-150', '150+']
    band_data = {fkey: [] for fkey in ['loop', 'quartic', 'full_db', 'hybrid']}
    band_ns = []
    for bname in band_names:
        total_n = 0
        sums = {fkey: 0 for fkey in band_data}
        for r in results:
            bl = r['bands']['loop'].get(bname, {})
            if not bl:
                continue
            n = bl['n']
            total_n += n
            for fkey in sums:
                b = r['bands'][fkey].get(bname, {})
                if b:
                    sums[fkey] += b['mae'] * n
        band_ns.append(total_n)
        for fkey in band_data:
            band_data[fkey].append(sums[fkey] / total_n if total_n > 0 else 0)

    bx = np.arange(len(band_names))
    ax.bar(bx - 1.5*w, band_data['loop'], w, label='Loop actual', color='#888888')
    ax.bar(bx - 0.5*w, band_data['quartic'], w, label='Quartic', color='#e74c3c')
    ax.bar(bx + 0.5*w, band_data['full_db'], w, label='Full Diabeloop', color='#e67e22')
    ax.bar(bx + 1.5*w, band_data['hybrid'], w, label='Hybrid', color='#2ecc71')
    ax.set_xticks(bx)
    ax.set_xticklabels([f"{bn}\n(n={n})" for bn, n in zip(band_names, band_ns)], fontsize=8)
    ax.set_ylabel('MAE (mg/dL)')
    ax.set_title('Weighted Average MAE by Glucose Band')
    ax.legend(fontsize=8)

    # Panel 4: Hybrid vs Loop MAE scatter
    ax = axes[1, 1]
    loop_maes = [r['overall']['loop']['mae'] for r in results]
    hybrid_maes = [r['overall']['hybrid']['mae'] for r in results]
    colors = ['#3498db' if r['model'] == 'sigmoid' else '#e74c3c' for r in results]
    ax.scatter(loop_maes, hybrid_maes, c=colors, s=80, zorder=5)
    for i, r in enumerate(results):
        ax.annotate(r['name'], (loop_maes[i], hybrid_maes[i]),
                    fontsize=7, ha='left', va='bottom')
    lims = [min(min(loop_maes), min(hybrid_maes)) - 2,
            max(max(loop_maes), max(hybrid_maes)) + 2]
    ax.plot(lims, lims, 'k--', alpha=0.3, linewidth=1)
    ax.set_xlabel('Loop Actual MAE (mg/dL)')
    ax.set_ylabel('Hybrid MAE (mg/dL)')
    ax.set_title('Hybrid vs Loop Actual MAE (blue=sigmoid, red=log)')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    plt.tight_layout()
    chart_path = OUT_DIR / 'multisite_backtest_results.png'
    plt.savefig(chart_path, dpi=150, bbox_inches='tight')
    print(f"Chart saved to {chart_path}")
    plt.close()

    print("\nDone.")
