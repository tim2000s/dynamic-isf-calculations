#!/usr/bin/env python3
"""
Simulate the Profile-Anchored Dynamic ISF with Zone Auto-Calibration.

Model: ISF(G) = ISF_profile × R(G)
  R(G) interpolates between control points at 76, 100, 130, 170 mg/dL.
  R(100) = 1.0 by definition.

Population defaults from 10-patient mean:
  R(76) = 1.15, R(100) = 1.00, R(130) = 0.80, R(170) = 0.70

Auto-calibration:
  Track prediction bias in 3 zones (<90, 90-130, >130).
  Every 7 days (or N observations), adjust zone ratios.
  Also track correction bolus ISF observations.

Test on all 13 subjects retrospectively.
"""

import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

OUT_DIR = Path(__file__).parent
DATA_DIR = OUT_DIR.parent
TRIO_CACHE = DATA_DIR / 'multisite_4h_sample_cache.pkl'
BOOST_CACHE = DATA_DIR / 'boost_4h_cache.pkl'

# ── Quartic (for comparison) ────────────────────────────────────────────────

def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

# ── The Anchored Dynamic ISF Model ──────────────────────────────────────────

# Control points: (glucose, default ratio)
DEFAULT_RATIOS = {76: 1.15, 100: 1.00, 130: 0.80, 170: 0.70}

# Bounds for auto-calibration
RATIO_BOUNDS = {76: (0.80, 2.00), 130: (0.40, 1.00), 170: (0.20, 1.00)}

CONTROL_GLUCOSES = sorted(DEFAULT_RATIOS.keys())  # [76, 100, 130, 170]


def ratio_function(g, ratios):
    """
    Compute R(G) by linear interpolation between control points.
    Extrapolate flat beyond the endpoints.
    """
    points = sorted(ratios.items())  # [(76, r1), (100, 1.0), (130, r2), (170, r3)]
    gs = [p[0] for p in points]
    rs = [p[1] for p in points]

    if g <= gs[0]:
        return rs[0]
    if g >= gs[-1]:
        return rs[-1]

    # Find segment
    for i in range(len(gs) - 1):
        if gs[i] <= g <= gs[i + 1]:
            t = (g - gs[i]) / (gs[i + 1] - gs[i])
            return rs[i] + t * (rs[i + 1] - rs[i])
    return 1.0


def anchored_isf(g, isf_profile, ratios):
    """Compute ISF at glucose g using the anchored model."""
    return isf_profile * ratio_function(g, ratios)


# ── Load data ───────────────────────────────────────────────────────────────

with open(TRIO_CACHE, 'rb') as f:
    trio_sites = pickle.load(f)

with open(BOOST_CACHE, 'rb') as f:
    boost_raw = pickle.load(f)

boost_df = boost_raw['strict']

# Anonymise site names
ANON_MAP = {
    'henny425': 'User-A', 'aadiabetes': 'User-B', 'diajesse': 'User-C',
    'svns': 'User-D', 'fuxchr': 'User-E', 'mikens': 'User-F',
    'andycgm': 'User-G', 'noahr': 'User-H', 'nightscout1': 'User-I',
    'eli': 'User-J', 'ns_rot6': 'User-K', 'kelseyhuss': 'User-L',
}
for s in trio_sites:
    s['name'] = ANON_MAP.get(s['name'], s['name'])

sites = []
for s in trio_sites:
    sites.append({
        'name': s['name'],
        'model': s['model'],
        'n': s['n'],
        'tdd': s['tdd_median'],
        'bg': s['bg'],
        'isf_actual': s['isf_actual'],
        'pred_drop': s['pred_drop'],
        'actual_bg_end': s['actual_bg_end'],
        'pred_loop': s['pred_loop'],
    })

boost_bg = boost_df['bg'].values.astype(float)
boost_isf = boost_df['variable_sens'].values.astype(float)
boost_pred_drop = boost_df['pred_drop'].values.astype(float)
boost_actual_end = boost_df['actual_bg_end'].values.astype(float)
boost_tdd = boost_df['tdd_7day'].median()

sites.append({
    'name': 'User-M',
    'model': 'AAPS',
    'n': len(boost_bg),
    'tdd': boost_tdd,
    'bg': boost_bg,
    'isf_actual': boost_isf,
    'pred_drop': boost_pred_drop,
    'actual_bg_end': boost_actual_end,
    'pred_loop': boost_bg - boost_pred_drop,
})

print(f"Loaded {len(sites)} sites, {sum(s['n'] for s in sites)} total samples\n")


# ── Compute each site's ISF_profile (median ISF near 100 mg/dL) ────────────

for s in sites:
    mask_ref = (s['bg'] >= 96) & (s['bg'] < 104)
    if mask_ref.sum() >= 5:
        s['isf_profile'] = np.median(s['isf_actual'][mask_ref])
    else:
        # Fallback: 1800/TDD
        s['isf_profile'] = 1800 / s['tdd']
    print(f"  {s['name']:15s}  ISF_profile={s['isf_profile']:6.1f}  TDD={s['tdd']:5.1f}  n={s['n']}")

print()


# ── Simulation 1: Population defaults (no calibration) ──────────────────────

print("=" * 80)
print("SIMULATION 1: Population defaults (R@76=1.15, R@130=0.80, R@170=0.70)")
print("=" * 80)


def evaluate_model(sites, get_ratios_fn, label="Model"):
    """Evaluate a model across all sites. get_ratios_fn(site) returns the ratios dict."""
    results = []
    for s in sites:
        bg = s['bg']
        isf_actual = s['isf_actual']
        pred_drop = s['pred_drop']
        actual_end = s['actual_bg_end']
        isf_profile = s['isf_profile']
        ratios = get_ratios_fn(s)

        # Compute anchored ISF for each sample
        isf_model = np.array([anchored_isf(g, isf_profile, ratios) for g in bg])

        # Counterfactual prediction
        pred_model = bg - pred_drop * (isf_model / isf_actual)

        # Zones
        zones = {
            '<105 falling': (bg < 105) & (pred_drop > 0),
            '<105 rising': (bg < 105) & (pred_drop <= 0),
            '>=105 falling': (bg >= 105) & (pred_drop > 0),
            '>=105 rising': (bg >= 105) & (pred_drop <= 0),
        }

        row = {'name': s['name'], 'n': s['n']}
        for zone, mask in zones.items():
            n = mask.sum()
            if n > 0:
                row[f'{zone}_bias'] = np.mean(pred_model[mask] - actual_end[mask])
                row[f'{zone}_mae'] = np.mean(np.abs(pred_model[mask] - actual_end[mask]))
                row[f'{zone}_n'] = n
            else:
                row[f'{zone}_bias'] = np.nan
                row[f'{zone}_mae'] = np.nan
                row[f'{zone}_n'] = 0

        row['overall_bias'] = np.mean(pred_model - actual_end)
        row['overall_mae'] = np.mean(np.abs(pred_model - actual_end))

        # Loop comparison
        pred_loop = s['pred_loop']
        row['loop_mae'] = np.mean(np.abs(pred_loop - actual_end))
        row['loop_bias'] = np.mean(pred_loop - actual_end)
        for zone, mask in zones.items():
            if mask.sum() > 0:
                row[f'{zone}_loop_bias'] = np.mean(pred_loop[mask] - actual_end[mask])

        results.append(row)
    return results


def print_results(results, label):
    print(f"\n  {label}:")
    print(f"  {'Site':15s} {'<105r':>8s} {'>=105f':>8s} {'MAE':>6s} {'Loop MAE':>9s}")
    print(f"  {'─' * 50}")
    for r in results:
        print(f"  {r['name']:15s} {r.get('<105 rising_bias', np.nan):+8.1f} "
              f"{r.get('>=105 falling_bias', np.nan):+8.1f} "
              f"{r['overall_mae']:6.1f} {r['loop_mae']:9.1f}")

    # Weighted aggregate
    total_n = sum(r['n'] for r in results)
    agg_mae = sum(r['overall_mae'] * r['n'] for r in results) / total_n
    loop_mae = sum(r['loop_mae'] * r['n'] for r in results) / total_n
    print(f"\n  Aggregate: MAE={agg_mae:.1f}, Loop MAE={loop_mae:.1f}")

    # Zone aggregates
    for zone in ['<105 falling', '<105 rising', '>=105 falling', '>=105 rising']:
        total_b = 0
        total_n_z = 0
        for r in results:
            n = r.get(f'{zone}_n', 0)
            if n > 0:
                total_b += r[f'{zone}_bias'] * n
                total_n_z += n
        if total_n_z > 0:
            print(f"  {zone:20s}  n={total_n_z:5d}  bias={total_b/total_n_z:+.1f}")


# Population defaults
pop_results = evaluate_model(sites, lambda s: dict(DEFAULT_RATIOS))
print_results(pop_results, "Population defaults")


# ── Simulation 2: Per-site optimal ratios ───────────────────────────────────

print(f"\n{'=' * 80}")
print("SIMULATION 2: Per-site optimal ratios (from binned ISF data)")
print("=" * 80)

BG_BINS = list(range(72, 201, 8))


def compute_optimal_ratios(s):
    """Compute the actual ISF ratios at control points for this site."""
    bg = s['bg']
    isf = s['isf_actual']
    isf_profile = s['isf_profile']

    ratios = {100: 1.0}
    for g_ctrl in [76, 130, 170]:
        # Find samples near this glucose
        window = 12  # ±12 mg/dL
        mask = (bg >= g_ctrl - window) & (bg < g_ctrl + window)
        if mask.sum() >= 5:
            median_isf = np.median(isf[mask])
            ratios[g_ctrl] = median_isf / isf_profile
        else:
            ratios[g_ctrl] = DEFAULT_RATIOS[g_ctrl]  # fallback to population
    return ratios


for s in sites:
    r = compute_optimal_ratios(s)
    print(f"  {s['name']:15s}  R@76={r.get(76, 'N/A'):5.2f}  R@130={r.get(130, 'N/A'):5.2f}  "
          f"R@170={r.get(170, 'N/A'):5.2f}")

opt_results = evaluate_model(sites, compute_optimal_ratios)
print_results(opt_results, "Per-site optimal ratios")


# ── Simulation 3: Auto-calibration from prediction data ─────────────────────

print(f"\n{'=' * 80}")
print("SIMULATION 3: Auto-calibration (weekly zone-based learning)")
print("=" * 80)

LEARNING_RATE = 0.05  # 5% per update
MIN_OBS = 20  # minimum observations per zone before adjusting
UPDATE_INTERVAL = 50  # update every N samples (approximating weekly)


def simulate_autocalibration(s):
    """
    Process samples sequentially. Track prediction bias in 3 zones.
    Periodically adjust ratios.
    """
    bg = s['bg']
    isf_actual = s['isf_actual']
    pred_drop = s['pred_drop']
    actual_end = s['actual_bg_end']
    isf_profile = s['isf_profile']
    n = len(bg)

    # Start with population defaults
    ratios = dict(DEFAULT_RATIOS)
    ratio_history = {76: [], 130: [], 170: []}

    # Zone accumulators
    zone_errors = {'low': [], 'mid': [], 'high': []}
    predictions = np.zeros(n)

    for i in range(n):
        g = bg[i]
        isf_m = anchored_isf(g, isf_profile, ratios)
        pred_sgv = g - pred_drop[i] * (isf_m / isf_actual[i])
        predictions[i] = pred_sgv

        error = pred_sgv - actual_end[i]

        # Assign to zone
        if g < 90:
            zone_errors['low'].append(error)
        elif g <= 130:
            zone_errors['mid'].append(error)
        else:
            zone_errors['high'].append(error)

        # Record ratio history
        ratio_history[76].append(ratios[76])
        ratio_history[130].append(ratios[130])
        ratio_history[170].append(ratios[170])

        # Periodic update
        if (i + 1) % UPDATE_INTERVAL == 0:
            # Zone low → adjust R(76)
            if len(zone_errors['low']) >= MIN_OBS:
                mean_bias = np.mean(zone_errors['low'][-MIN_OBS:])
                # Negative bias (under-predict) → ISF too high at low G → reduce R(76)
                # Positive bias (over-predict) → ISF too low → increase R(76)
                # But direction depends on falling vs rising...
                # Actually: if ISF_model > ISF_actual, pred_drop is amplified, pred_sgv is lower
                # So negative bias → ISF_model too high → reduce ratio
                # Positive bias → ISF_model too low → increase ratio
                adjustment = -mean_bias / isf_profile * LEARNING_RATE
                ratios[76] = np.clip(ratios[76] + adjustment,
                                     RATIO_BOUNDS[76][0], RATIO_BOUNDS[76][1])

            # Zone high → adjust R(130) and R(170) together
            if len(zone_errors['high']) >= MIN_OBS:
                mean_bias = np.mean(zone_errors['high'][-MIN_OBS:])
                adjustment = -mean_bias / isf_profile * LEARNING_RATE
                ratios[130] = np.clip(ratios[130] + adjustment,
                                      RATIO_BOUNDS[130][0], RATIO_BOUNDS[130][1])
                ratios[170] = np.clip(ratios[170] + adjustment * 0.8,  # slightly less at extreme
                                      RATIO_BOUNDS[170][0], RATIO_BOUNDS[170][1])

            # Zone mid → would adjust ISF_profile (skip for now, it's the anchor)

            # Reset accumulators
            zone_errors = {'low': [], 'mid': [], 'high': []}

    return {
        'predictions': predictions,
        'ratios_final': dict(ratios),
        'ratio_history': {k: np.array(v) for k, v in ratio_history.items()},
    }


auto_results_raw = {}
for s in sites:
    ar = simulate_autocalibration(s)
    auto_results_raw[s['name']] = ar
    rf = ar['ratios_final']
    print(f"  {s['name']:15s}  R@76: 1.15→{rf[76]:.2f}  "
          f"R@130: 0.80→{rf[130]:.2f}  R@170: 0.70→{rf[170]:.2f}")

# Evaluate auto-calibrated model
def get_auto_predictions(sites):
    results = []
    for s in sites:
        ar = auto_results_raw[s['name']]
        pred = ar['predictions']
        actual_end = s['actual_bg_end']
        bg = s['bg']
        pred_drop = s['pred_drop']

        zones = {
            '<105 falling': (bg < 105) & (pred_drop > 0),
            '<105 rising': (bg < 105) & (pred_drop <= 0),
            '>=105 falling': (bg >= 105) & (pred_drop > 0),
            '>=105 rising': (bg >= 105) & (pred_drop <= 0),
        }

        row = {'name': s['name'], 'n': s['n']}
        for zone, mask in zones.items():
            n = mask.sum()
            if n > 0:
                row[f'{zone}_bias'] = np.mean(pred[mask] - actual_end[mask])
                row[f'{zone}_n'] = n
        row['overall_mae'] = np.mean(np.abs(pred - actual_end))
        row['loop_mae'] = np.mean(np.abs(s['pred_loop'] - actual_end))
        results.append(row)
    return results

auto_results = get_auto_predictions(sites)
print_results(auto_results, "Auto-calibrated")


# ── Simulation 4: Quartic + TDD (for comparison) ───────────────────────────

print(f"\n{'=' * 80}")
print("SIMULATION 4: Quartic + TDD scaling (comparison baseline)")
print("=" * 80)

Q_REF = quartic(99)

def quartic_ratios(s):
    """Return ratios that reproduce the quartic shape."""
    q100 = quartic(100)
    return {76: quartic(76)/q100, 100: 1.0, 130: quartic(130)/q100, 170: quartic(170)/q100}

quartic_results = evaluate_model(sites, quartic_ratios)
print_results(quartic_results, "Quartic shape via anchored model")


# ══════════════════════════════════════════════════════════════════════════════
# AGGREGATE COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'=' * 80}")
print("AGGREGATE COMPARISON — All 13 subjects")
print("=" * 80)

models = {
    'Loop (calibrated)': None,
    'Quartic + TDD': quartic_results,
    'Anchored (pop defaults)': pop_results,
    'Anchored (optimal)': opt_results,
    'Anchored (auto-cal)': auto_results,
}

zones = ['<105 falling', '<105 rising', '>=105 falling', '>=105 rising']

print(f"\n  {'Model':<25s}", end='')
for z in zones:
    print(f" {z:>15s}", end='')
print(f" {'MAE':>8s}")
print(f"  {'─' * 95}")

for model_name, results in models.items():
    print(f"  {model_name:<25s}", end='')
    if model_name == 'Loop (calibrated)':
        # Compute from sites directly
        for z in zones:
            total_b = 0
            total_n = 0
            for s in sites:
                bg = s['bg']
                pd_ = s['pred_drop']
                ae = s['actual_bg_end']
                pl = s['pred_loop']
                if z == '<105 falling':
                    mask = (bg < 105) & (pd_ > 0)
                elif z == '<105 rising':
                    mask = (bg < 105) & (pd_ <= 0)
                elif z == '>=105 falling':
                    mask = (bg >= 105) & (pd_ > 0)
                else:
                    mask = (bg >= 105) & (pd_ <= 0)
                n = mask.sum()
                if n > 0:
                    total_b += np.mean(pl[mask] - ae[mask]) * n
                    total_n += n
            print(f" {total_b/total_n:+15.1f}", end='')
        total_mae = sum(np.mean(np.abs(s['pred_loop'] - s['actual_bg_end'])) * s['n']
                       for s in sites) / sum(s['n'] for s in sites)
        print(f" {total_mae:8.1f}")
    else:
        for z in zones:
            total_b = 0
            total_n = 0
            for r in results:
                n = r.get(f'{z}_n', 0)
                if n > 0:
                    total_b += r[f'{z}_bias'] * n
                    total_n += n
            print(f" {total_b/total_n:+15.1f}" if total_n > 0 else f" {'N/A':>15s}", end='')
        total_n_all = sum(r['n'] for r in results)
        agg_mae = sum(r['overall_mae'] * r['n'] for r in results) / total_n_all
        print(f" {agg_mae:8.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# ── Figure 1: The model — ratio function with patient data ──────────────────

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

bg_range = np.linspace(72, 200, 200)

# Panel 1: Population default ratio curve vs quartic
ax = axes[0]
pop_curve = [ratio_function(g, DEFAULT_RATIOS) for g in bg_range]
q100 = quartic(100)
quartic_curve = [quartic(g) / q100 for g in bg_range]
ax.plot(bg_range, pop_curve, 'b-', linewidth=2.5, label='Anchored (population default)')
ax.plot(bg_range, quartic_curve, 'r--', linewidth=1.5, label='Quartic (normalised)')
# Control points
for g, r in DEFAULT_RATIOS.items():
    ax.plot(g, r, 'bo', markersize=8, zorder=5)
ax.axhline(1.0, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('R(G) — ISF ratio relative to ISF at target')
ax.set_title('ISF Ratio Function: Population Default vs Quartic')
ax.legend(fontsize=9)
ax.set_xlim(72, 200)
ax.set_ylim(0.2, 1.8)

# Panel 2: Patient actual ratios overlaid
ax = axes[1]
ax.plot(bg_range, pop_curve, 'b-', linewidth=2.5, alpha=0.5, label='Population default')
BG_BINS_C = list(range(72, 201, 8))
BG_CENTERS = [(BG_BINS_C[i] + BG_BINS_C[i+1])/2 for i in range(len(BG_BINS_C)-1)]
ref_idx = 3  # ~100 mg/dL bin

for s in sites:
    if s['name'] in ('User-I', 'User-H', 'User-J'):
        continue
    bg = s['bg']
    isf = s['isf_actual']
    medians = []
    for i in range(len(BG_BINS_C) - 1):
        mask = (bg >= BG_BINS_C[i]) & (bg < BG_BINS_C[i+1])
        medians.append(np.median(isf[mask]) if mask.sum() >= 5 else np.nan)
    medians = np.array(medians)
    ref_isf = medians[ref_idx]
    if np.isnan(ref_isf) or ref_isf <= 0:
        continue
    normalised = medians / ref_isf
    style = '-' if s['model'] == 'sigmoid' else '--'
    ax.plot(BG_CENTERS, normalised, style, linewidth=1, alpha=0.6, label=s['name'])

ax.axhline(1.0, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('Normalised ISF')
ax.set_title('Population Default vs Individual Patient Curves')
ax.legend(fontsize=6, loc='upper right', ncol=2)
ax.set_xlim(72, 200)
ax.set_ylim(0.2, 2.0)

# Panel 3: Auto-calibrated ratios per site
ax = axes[2]
site_names = [s['name'] for s in sites]
x = np.arange(len(sites))
width = 0.25

r76 = [auto_results_raw[s['name']]['ratios_final'][76] for s in sites]
r130 = [auto_results_raw[s['name']]['ratios_final'][130] for s in sites]
r170 = [auto_results_raw[s['name']]['ratios_final'][170] for s in sites]

ax.bar(x - width, r76, width, label='R(76)', color='#4A90D9', alpha=0.8)
ax.bar(x, r130, width, label='R(130)', color='#F5A623', alpha=0.8)
ax.bar(x + width, r170, width, label='R(170)', color='#D94A4A', alpha=0.8)
ax.axhline(1.15, color='#4A90D9', linestyle='--', alpha=0.4)
ax.axhline(0.80, color='#F5A623', linestyle='--', alpha=0.4)
ax.axhline(0.70, color='#D94A4A', linestyle='--', alpha=0.4)
ax.set_xticks(x)
ax.set_xticklabels(site_names, rotation=45, ha='right', fontsize=7)
ax.set_ylabel('Ratio value')
ax.set_title('Auto-Calibrated Ratios Per Site\n(dashed = population defaults)')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(OUT_DIR / 'anchored_model_overview.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR / 'anchored_model_overview.png'}")


# ── Figure 2: Auto-calibration convergence per site ─────────────────────────

fig, axes = plt.subplots(3, 5, figsize=(20, 10))
axes_flat = axes.flatten()

for i, s in enumerate(sites):
    if i >= 13:
        break
    ax = axes_flat[i]
    ar = auto_results_raw[s['name']]

    ax.plot(ar['ratio_history'][76], linewidth=0.8, color='#4A90D9', label='R(76)')
    ax.plot(ar['ratio_history'][130], linewidth=0.8, color='#F5A623', label='R(130)')
    ax.plot(ar['ratio_history'][170], linewidth=0.8, color='#D94A4A', label='R(170)')

    ax.axhline(1.15, color='#4A90D9', linestyle=':', alpha=0.3)
    ax.axhline(0.80, color='#F5A623', linestyle=':', alpha=0.3)
    ax.axhline(0.70, color='#D94A4A', linestyle=':', alpha=0.3)
    ax.axhline(1.0, color='grey', linestyle=':', alpha=0.2)

    ax.set_title(f"{s['name']} (n={s['n']})", fontsize=8)
    ax.set_ylim(0.1, 2.1)
    if i >= 8:
        ax.set_xlabel('Sample', fontsize=7)
    if i % 5 == 0:
        ax.set_ylabel('Ratio', fontsize=7)
    if i == 0:
        ax.legend(fontsize=6)

for j in range(len(sites), len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.suptitle('Auto-Calibration Convergence: Zone Ratios Over Time', fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'anchored_convergence.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'anchored_convergence.png'}")


# ── Figure 3: Per-site ISF curves — actual vs anchored models ───────────────

fig, axes = plt.subplots(3, 5, figsize=(20, 10))
axes_flat = axes.flatten()

for i, s in enumerate(sites):
    if i >= 13:
        break
    ax = axes_flat[i]

    # Actual ISF bars
    bg = s['bg']
    isf = s['isf_actual']
    medians = []
    for bi in range(len(BG_BINS_C) - 1):
        mask = (bg >= BG_BINS_C[bi]) & (bg < BG_BINS_C[bi+1])
        medians.append(np.median(isf[mask]) if mask.sum() >= 5 else np.nan)

    ax.bar(BG_CENTERS, medians, width=7, alpha=0.3, color='steelblue',
           edgecolor='navy', linewidth=0.5, label='Actual')

    isf_p = s['isf_profile']

    # Population default
    y_pop = [anchored_isf(g, isf_p, DEFAULT_RATIOS) for g in bg_range]
    ax.plot(bg_range, y_pop, 'g-', linewidth=1.5, alpha=0.8, label='Pop default')

    # Auto-calibrated
    rf = auto_results_raw[s['name']]['ratios_final']
    y_auto = [anchored_isf(g, isf_p, rf) for g in bg_range]
    ax.plot(bg_range, y_auto, color='orange', linewidth=1.5, alpha=0.8, label='Auto-cal')

    # Quartic + TDD
    S = (1800 / s['tdd']) / quartic(99)
    y_quartic = [quartic(g) * S for g in bg_range]
    ax.plot(bg_range, y_quartic, 'r--', linewidth=1, alpha=0.6, label='Quartic+TDD')

    ax.set_title(f"{s['name']} (TDD={s['tdd']:.0f})", fontsize=8)
    ax.set_xlim(72, 200)
    ax.set_ylim(0, min(max(np.nanmax(medians) * 1.3, 50), 400))
    if i >= 8:
        ax.set_xlabel('BG', fontsize=7)
    if i % 5 == 0:
        ax.set_ylabel('ISF', fontsize=7)
    ax.legend(fontsize=5, loc='upper right')

for j in range(len(sites), len(axes_flat)):
    axes_flat[j].set_visible(False)

plt.suptitle('Per-Site ISF: Actual vs Population Default vs Auto-Calibrated vs Quartic',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'anchored_isf_per_site.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'anchored_isf_per_site.png'}")


# ── Figure 4: Bias comparison bar chart ─────────────────────────────────────

fig, axes = plt.subplots(1, 4, figsize=(18, 5))
model_labels = ['Loop', 'Quartic+TDD', 'Anchored\n(pop)', 'Anchored\n(optimal)', 'Anchored\n(auto-cal)']
colors = ['#4A90D9', '#D94A4A', '#7ED321', '#F5A623', '#9B59B6']

all_model_results = [None, quartic_results, pop_results, opt_results, auto_results]

for zi, zone in enumerate(zones):
    ax = axes[zi]
    biases = []
    for mi, (label, results) in enumerate(zip(model_labels, all_model_results)):
        if results is None:
            # Loop
            total_b = 0
            total_n = 0
            for s in sites:
                bg = s['bg']
                pd_ = s['pred_drop']
                ae = s['actual_bg_end']
                pl = s['pred_loop']
                if zone == '<105 falling':
                    mask = (bg < 105) & (pd_ > 0)
                elif zone == '<105 rising':
                    mask = (bg < 105) & (pd_ <= 0)
                elif zone == '>=105 falling':
                    mask = (bg >= 105) & (pd_ > 0)
                else:
                    mask = (bg >= 105) & (pd_ <= 0)
                n = mask.sum()
                if n > 0:
                    total_b += np.mean(pl[mask] - ae[mask]) * n
                    total_n += n
            biases.append(total_b / total_n if total_n > 0 else 0)
        else:
            total_b = 0
            total_n = 0
            for r in results:
                n = r.get(f'{zone}_n', 0)
                if n > 0:
                    total_b += r[f'{zone}_bias'] * n
                    total_n += n
            biases.append(total_b / total_n if total_n > 0 else 0)

    bars = ax.bar(range(len(biases)), biases, color=colors, alpha=0.8,
                  edgecolor='black', linewidth=0.5)
    ax.axhline(0, color='grey', linewidth=0.5)
    ax.set_xticks(range(len(biases)))
    ax.set_xticklabels(model_labels, fontsize=7)
    ax.set_ylabel('Bias (mg/dL)')
    ax.set_title(zone, fontsize=10)
    for bar, val in zip(bars, biases):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{val:+.1f}', ha='center',
                va='bottom' if val >= 0 else 'top', fontsize=7)

plt.suptitle('Prediction Bias by Zone: All 13 Subjects\n'
             '(predicted SGV at end-of-IOB minus actual SGV)', fontsize=11)
plt.tight_layout()
plt.savefig(OUT_DIR / 'anchored_bias_comparison.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'anchored_bias_comparison.png'}")

print(f"\nDONE")
