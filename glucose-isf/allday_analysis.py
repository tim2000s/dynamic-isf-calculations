#!/usr/bin/env python3
"""
All-Day Fasting Analysis: Compare overnight vs daytime vs combined
for all ISF models across 13 subjects.

Uses the allday caches (COB=0, bolus_age >= 120/180min, all hours).
"""

import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path(__file__).parent
TRIO_CACHE = OUT_DIR / 'multisite_allday_cache.pkl'
BOOST_CACHE = OUT_DIR / 'boost_allday_cache.pkl'


def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

Q_REF = quartic(100)


# ── Load data ──────────────────────────────────────────────────────────────

with open(TRIO_CACHE, 'rb') as f:
    trio_sites = pickle.load(f)

with open(BOOST_CACHE, 'rb') as f:
    boost_cache = pickle.load(f)

print(f"Loaded {len(trio_sites)} Trio sites")
print(f"Boost: allday={len(boost_cache['allday'])}, overnight={len(boost_cache['overnight'])}, "
      f"daytime={len(boost_cache['daytime'])}")


# ── Build unified site list for each time period ──────────────────────────

def build_site_list(trio_sites, boost_cache, period='allday'):
    """Build a unified list of sites for a given period (allday/overnight/daytime)."""
    sites = []
    for s in trio_sites:
        data = s.get(period)
        if data is None or data['n'] < 10:
            continue
        bg = data['bg']
        isf = data['isf_actual']
        m100 = (bg >= 96) & (bg < 104)
        isf100 = np.median(isf[m100]) if m100.sum() >= 5 else np.nan
        sites.append({
            'name': s['name'], 'model': s['model'],
            'tdd': s['tdd_median'], 'n': data['n'],
            'bg': bg, 'isf_actual': isf,
            'pred_drop': data['pred_drop'],
            'actual_bg_end': data['actual_bg_end'],
            'pred_loop': data['pred_loop'],
            'hour': data.get('hour'),
            'isf_true': isf100,
            'isf_tdd': 1800 / s['tdd_median'],
        })

    # Boost
    boost_df = boost_cache.get(period)
    if boost_df is not None and len(boost_df) >= 10:
        bb = boost_df['bg'].values.astype(float)
        bi = boost_df['variable_sens'].values.astype(float)
        m100 = (bb >= 96) & (bb < 104)
        isf100 = np.median(bi[m100]) if m100.sum() >= 5 else np.nan
        tdd = boost_df['tdd_7day'].median()
        sites.append({
            'name': 'User-M', 'model': 'AAPS',
            'tdd': tdd, 'n': len(bb),
            'bg': bb, 'isf_actual': bi,
            'pred_drop': boost_df['pred_drop'].values.astype(float),
            'actual_bg_end': boost_df['actual_bg_end'].values.astype(float),
            'pred_loop': bb - boost_df['pred_drop'].values.astype(float),
            'hour': boost_df['hour'].values.astype(int),
            'isf_true': isf100,
            'isf_tdd': 1800 / tdd,
        })

    return sites


def evaluate_models(sites, label=""):
    """Evaluate all models on a set of sites. Returns results dict."""
    if not sites:
        return None

    valid = [s for s in sites if not np.isnan(s['isf_true'])]
    if not valid:
        return None

    weights = np.array([s['n'] for s in valid])

    results = {}

    # 1. Loop actual
    loop_maes = [np.mean(np.abs(s['pred_loop'] - s['actual_bg_end'])) for s in valid]
    results['Loop (tuned)'] = {'wmae': np.average(loop_maes, weights=weights),
                                'per_site': {s['name']: m for s, m in zip(valid, loop_maes)}}

    # 2. Quartic + 1800/TDD
    q_maes = []
    for s in valid:
        isf_model = np.array([s['isf_tdd'] * quartic(g) / Q_REF for g in s['bg']])
        pred = s['bg'] - s['pred_drop'] * (isf_model / s['isf_actual'])
        q_maes.append(np.mean(np.abs(pred - s['actual_bg_end'])))
    results['Quartic+TDD'] = {'wmae': np.average(q_maes, weights=weights),
                               'per_site': {s['name']: m for s, m in zip(valid, q_maes)}}

    # 3. Profile-anchored + population ratios
    RATIOS = {76: 1.15, 100: 1.00, 130: 0.80, 170: 0.70}

    def ratio_fn(g):
        points = sorted(RATIOS.items())
        gs = [p[0] for p in points]; rs = [p[1] for p in points]
        if g <= gs[0]: return rs[0]
        if g >= gs[-1]: return rs[-1]
        for i in range(len(gs) - 1):
            if gs[i] <= g <= gs[i + 1]:
                t = (g - gs[i]) / (gs[i + 1] - gs[i])
                return rs[i] + t * (rs[i + 1] - rs[i])
        return 1.0

    anch_maes = []
    for s in valid:
        isf_model = np.array([s['isf_true'] * ratio_fn(g) for g in s['bg']])
        pred = s['bg'] - s['pred_drop'] * (isf_model / s['isf_actual'])
        anch_maes.append(np.mean(np.abs(pred - s['actual_bg_end'])))
    results['Anchored+pop'] = {'wmae': np.average(anch_maes, weights=weights),
                                'per_site': {s['name']: m for s, m in zip(valid, anch_maes)}}

    # 4. Profile ISF + quartic
    profq_maes = []
    for s in valid:
        isf_model = np.array([s['isf_true'] * quartic(g) / Q_REF for g in s['bg']])
        pred = s['bg'] - s['pred_drop'] * (isf_model / s['isf_actual'])
        profq_maes.append(np.mean(np.abs(pred - s['actual_bg_end'])))
    results['Profile+Quartic'] = {'wmae': np.average(profq_maes, weights=weights),
                                   'per_site': {s['name']: m for s, m in zip(valid, profq_maes)}}

    # 5. Auto-learning anchor (rolling median, quartic shape)
    auto_maes_overall = []
    auto_maes_last200 = []
    for s in valid:
        bg = s['bg']; isf_actual = s['isf_actual']
        pred_drop = s['pred_drop']; actual_end = s['actual_bg_end']
        isf_est = s['isf_tdd']
        implied = []
        errors = []

        for i in range(len(bg)):
            g = bg[i]
            ratio = quartic(g) / Q_REF
            isf_model = isf_est * ratio
            pred = g - pred_drop[i] * (isf_model / isf_actual[i])
            errors.append(abs(pred - actual_end[i]))

            actual_drop = g - actual_end[i]
            if abs(actual_drop) > 10 and abs(pred_drop[i]) > 5:
                correct_isf = isf_actual[i] * actual_drop / pred_drop[i]
                implied_isf_100 = correct_isf / ratio if ratio > 0 else isf_est
                if 3 < implied_isf_100 < 600:
                    implied.append(implied_isf_100)

            if len(implied) >= 10:
                isf_est = np.median(implied[-100:])

        errors = np.array(errors)
        auto_maes_overall.append(np.mean(errors))
        auto_maes_last200.append(np.mean(errors[-200:]) if len(errors) >= 200 else np.mean(errors[-100:]))

    results['Auto-learn (overall)'] = {
        'wmae': np.average(auto_maes_overall, weights=weights),
        'per_site': {s['name']: m for s, m in zip(valid, auto_maes_overall)}}
    results['Auto-learn (last200)'] = {
        'wmae': np.average(auto_maes_last200, weights=weights),
        'per_site': {s['name']: m for s, m in zip(valid, auto_maes_last200)}}

    return results, valid, weights


# ── Run analysis for each time period ──────────────────────────────────────

periods = ['allday', 'overnight', 'daytime']
all_results = {}

for period in periods:
    print(f"\n{'=' * 80}")
    print(f"PERIOD: {period.upper()}")
    print(f"{'=' * 80}")

    sites = build_site_list(trio_sites, boost_cache, period)
    print(f"  Sites: {len(sites)}, Samples: {sum(s['n'] for s in sites)}")

    for s in sites:
        print(f"    {s['name']:8s}  n={s['n']:5d}  model={s['model']:7s}  "
              f"TDD={s['tdd']:5.1f}  ISF@100={'%.0f' % s['isf_true'] if not np.isnan(s['isf_true']) else 'N/A':>5s}")

    result = evaluate_models(sites)
    if result is None:
        print("  SKIP (insufficient data)")
        continue

    results, valid, weights = result
    all_results[period] = results

    print(f"\n  {'Model':25s} {'Wt MAE':>8s}")
    print(f"  {'-' * 35}")
    for model, data in results.items():
        print(f"  {model:25s} {data['wmae']:8.1f}")


# ── Comparison table ───────────────────────────────────────────────────────

print(f"\n\n{'=' * 80}")
print("COMPARISON: Overnight vs Daytime vs All-Day")
print(f"{'=' * 80}\n")

models_to_compare = ['Loop (tuned)', 'Quartic+TDD', 'Anchored+pop', 'Profile+Quartic',
                     'Auto-learn (overall)', 'Auto-learn (last200)']

header = f"{'Model':25s}"
for period in periods:
    header += f" {period:>12s}"
print(header)
print("-" * 65)

for model in models_to_compare:
    row = f"{model:25s}"
    for period in periods:
        if period in all_results and model in all_results[period]:
            row += f" {all_results[period][model]['wmae']:12.1f}"
        else:
            row += f" {'N/A':>12s}"
    print(row)


# ── Per-site overnight vs daytime comparison ───────────────────────────────

print(f"\n\n{'=' * 80}")
print("PER-SITE: Loop MAE overnight vs daytime")
print(f"{'=' * 80}\n")

if 'overnight' in all_results and 'daytime' in all_results:
    on_sites = all_results['overnight']['Loop (tuned)']['per_site']
    dt_sites = all_results['daytime']['Loop (tuned)']['per_site']

    print(f"{'Site':8s} {'Overnight':>10s} {'Daytime':>10s} {'Diff':>8s}")
    print("-" * 40)
    for name in sorted(set(on_sites.keys()) & set(dt_sites.keys())):
        on_mae = on_sites[name]
        dt_mae = dt_sites[name]
        diff = dt_mae - on_mae
        print(f"{name:8s} {on_mae:10.1f} {dt_mae:10.1f} {diff:+8.1f}")


# ── Charts ─────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("All-Day Fasting Analysis: Overnight vs Daytime ISF Model Performance",
             fontsize=13, fontweight='bold')

# Chart 1: Model comparison across periods
ax = axes[0, 0]
models_short = ['Loop', 'Q+TDD', 'Anch+pop', 'Prof+Q', 'Auto(all)', 'Auto(L200)']
x = np.arange(len(models_short))
width = 0.25
colors = {'overnight': 'tab:blue', 'daytime': 'tab:orange', 'allday': 'tab:green'}

for i, period in enumerate(periods):
    if period not in all_results: continue
    vals = [all_results[period].get(m, {}).get('wmae', np.nan) for m in models_to_compare]
    ax.bar(x + (i - 1) * width, vals, width, label=period.capitalize(),
           color=colors[period], alpha=0.8)

ax.set_xticks(x)
ax.set_xticklabels(models_short, rotation=30, ha='right', fontsize=8)
ax.set_ylabel('Weighted Mean MAE (mg/dL)')
ax.set_title('A. Model Performance by Time Period')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

# Chart 2: Per-site Loop MAE overnight vs daytime
ax = axes[0, 1]
if 'overnight' in all_results and 'daytime' in all_results:
    on_sites = all_results['overnight']['Loop (tuned)']['per_site']
    dt_sites = all_results['daytime']['Loop (tuned)']['per_site']
    common = sorted(set(on_sites.keys()) & set(dt_sites.keys()))
    x_pos = np.arange(len(common))
    on_vals = [on_sites[n] for n in common]
    dt_vals = [dt_sites[n] for n in common]
    ax.bar(x_pos - 0.15, on_vals, 0.3, label='Overnight', color='tab:blue', alpha=0.8)
    ax.bar(x_pos + 0.15, dt_vals, 0.3, label='Daytime', color='tab:orange', alpha=0.8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(common, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Loop MAE (mg/dL)')
    ax.set_title('B. Loop MAE: Overnight vs Daytime per Site')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')

# Chart 3: Sample distribution by hour
ax = axes[1, 0]
allday_sites = build_site_list(trio_sites, boost_cache, 'allday')
all_hours = np.concatenate([s['hour'] for s in allday_sites if s.get('hour') is not None])
hour_counts = np.bincount(all_hours.astype(int), minlength=24)
ax.bar(range(24), hour_counts, color='tab:green', alpha=0.7)
ax.axvline(7.5, color='red', linestyle='--', alpha=0.5, label='Overnight/daytime split')
ax.set_xlabel('Hour of day')
ax.set_ylabel('Number of fasting samples')
ax.set_title('C. Fasting Sample Distribution by Hour')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

# Chart 4: Daytime improvement — which model benefits most?
ax = axes[1, 1]
if 'overnight' in all_results and 'allday' in all_results:
    models_for_delta = ['Loop (tuned)', 'Quartic+TDD', 'Anchored+pop', 'Auto-learn (overall)']
    labels = ['Loop', 'Q+TDD', 'Anch+pop', 'Auto-learn']
    on_wmaes = [all_results['overnight'].get(m, {}).get('wmae', np.nan) for m in models_for_delta]
    all_wmaes = [all_results['allday'].get(m, {}).get('wmae', np.nan) for m in models_for_delta]
    deltas = [a - o for a, o in zip(all_wmaes, on_wmaes)]
    bars = ax.bar(range(len(labels)), deltas, color=['green' if d < 0 else 'red' for d in deltas], alpha=0.7)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('MAE change (allday − overnight)')
    ax.set_title('D. Adding Daytime Data: MAE Change')
    ax.axhline(0, color='black', linewidth=0.5)
    for bar, val in zip(bars, deltas):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{val:+.1f}', ha='center', fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(OUT_DIR / 'allday_analysis.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR / 'allday_analysis.png'}")
print("\nDONE")
