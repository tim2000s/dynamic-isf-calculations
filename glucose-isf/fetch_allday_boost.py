#!/usr/bin/env python3
"""
Rebuild Boost cache with ALL HOURS (not just overnight).
Filter: COB=0, bolus_age >= 180min, BG 72-200.
Saves both 'allday' and 'overnight' subsets.
"""

import json, glob, warnings, pickle
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings('ignore')

D = 82.0
TARGET = 99.0
HOME = Path.home()
NS_WORK = HOME / 'Nightscout_Work'
OUT_DIR = HOME / 'Downloads' / '4 Hour analysis'
CACHE = OUT_DIR / 'speculation' / 'boost_allday_cache.pkl'


def find_json(prefix):
    patterns = [
        str(NS_WORK / f'{prefix}_*.json'),
        str(NS_WORK / '*' / f'{prefix}_*.json'),
        str(HOME / f'{prefix}_*.json'),
    ]
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
print("BOOST ALL-DAY FASTING CACHE BUILD")
print("=" * 70)

print("\nLoading devicestatus...")
raw_ds = load_dedup(find_json('devicestatus'))
print("\nLoading CGM entries...")
raw_entries = load_dedup(find_json('entries'))
print("\nLoading treatments...")
raw_tx = load_dedup(find_json('treatments'))

# ── Compute actual TDD from treatments ──
print("\n── Computing TDD from treatments ──")
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
        if len(pred_iob) < 12: continue
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

# ── CGM ──
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


print("\n── Forward SGV lookups ──")
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

# ── Attach TDD ──
print("\n── Attaching TDD ──")
daily_tdd['date'] = pd.to_datetime(daily_tdd['date']).dt.date
df = df.merge(daily_tdd[['date', 'tdd_actual', 'tdd_7day', 'tdd_1day']], on='date', how='left')

tdd_boost = np.full(len(df), np.nan)
for i, (epoch, tdd7, tdd1) in enumerate(zip(cycle_epochs, df['tdd_7day'].values, df['tdd_1day'].values)):
    tdd7_v = tdd7 if not np.isnan(tdd7) else None
    tdd1_v = tdd1 if not np.isnan(tdd1) else None
    tdd_boost[i] = compute_boost_tdd(epoch, tdd7_v, tdd1_v)
df['tdd_boost'] = tdd_boost

# ── ALL-DAY FASTING FILTER (COB=0, bolus_age >= 180, any hour) ──
print("\n── Filtering: all-day fasting (COB=0, bolus_age >= 180min, any hour) ──")
mask_allday = (
    (df['cob'].fillna(99) == 0) &
    (df['bg'] >= 72) & (df['bg'] <= 200) &
    (df['bolus_age_min'] >= 180) &
    (df['tdd_boost'] > 0) &
    (~np.isnan(df['tdd_7day']))
)
allday = df[mask_allday].copy()
allday = allday.dropna(subset=['pred_iob_final', 'actual_bg_end']).copy()
allday['pred_drop'] = allday['bg'] - allday['pred_iob_final']
allday['actual_drop'] = allday['bg'] - allday['actual_bg_end']

strict_allday = allday[
    (np.abs(allday['pred_drop']) > 3) &
    ((allday['actual_drop'] / allday['pred_drop']).between(0, 5))
].copy().sort_values('ts').reset_index(drop=True)

# Also create overnight subset for comparison
mask_overnight = mask_allday & (df['hour'] < 8)
overnight = df[mask_overnight].copy()
overnight = overnight.dropna(subset=['pred_iob_final', 'actual_bg_end']).copy()
overnight['pred_drop'] = overnight['bg'] - overnight['pred_iob_final']
overnight['actual_drop'] = overnight['bg'] - overnight['actual_bg_end']
strict_overnight = overnight[
    (np.abs(overnight['pred_drop']) > 3) &
    ((overnight['actual_drop'] / overnight['pred_drop']).between(0, 5))
].copy().sort_values('ts').reset_index(drop=True)

# Daytime-only subset
mask_daytime = mask_allday & (df['hour'] >= 8)
daytime = df[mask_daytime].copy()
daytime = daytime.dropna(subset=['pred_iob_final', 'actual_bg_end']).copy()
daytime['pred_drop'] = daytime['bg'] - daytime['pred_iob_final']
daytime['actual_drop'] = daytime['bg'] - daytime['actual_bg_end']
strict_daytime = daytime[
    (np.abs(daytime['pred_drop']) > 3) &
    ((daytime['actual_drop'] / daytime['pred_drop']).between(0, 5))
].copy().sort_values('ts').reset_index(drop=True)

print(f"  All-day fasting:  {len(strict_allday):,} samples")
print(f"  Overnight only:   {len(strict_overnight):,} samples")
print(f"  Daytime fasting:  {len(strict_daytime):,} samples")

# Hour distribution
print(f"\n  Hour distribution (all-day):")
hour_counts = strict_allday.groupby('hour').size()
for h, c in hour_counts.items():
    print(f"    {h:02d}:00  {c:4d} {'█' * (c // 20)}")

# Cache
with open(CACHE, 'wb') as f:
    pickle.dump({
        'allday': strict_allday,
        'overnight': strict_overnight,
        'daytime': strict_daytime,
    }, f)
print(f"\nSaved: {CACHE}")
print("DONE")
