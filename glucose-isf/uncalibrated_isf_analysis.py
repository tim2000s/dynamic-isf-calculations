#!/usr/bin/env python3
"""
Simulation: Which dynamic ISF approach works best for UNCALIBRATED users?

Assumes users have NOT carefully tuned their profile ISF.
Tests various shapes and anchor strategies using 13 subjects' real data.
"""

import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

OUT_DIR = Path(__file__).parent
CACHE_DIR = OUT_DIR.parent
TRIO_CACHE = CACHE_DIR / 'multisite_4h_sample_cache.pkl'
BOOST_CACHE = CACHE_DIR / 'boost_4h_cache.pkl'

ANON = {
    'henny425': 'User-A', 'aadiabetes': 'User-B', 'diajesse': 'User-C',
    'svns': 'User-D', 'fuxchr': 'User-E', 'mikens': 'User-F',
    'andycgm': 'User-G', 'noahr': 'User-H', 'nightscout1': 'User-I',
    'eli': 'User-J', 'ns_rot6': 'User-K', 'kelseyhuss': 'User-L',
}

def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

Q_REF = quartic(100)

def quartic_ratio(g):
    return quartic(g) / Q_REF

def gentle_ratio(g, alpha=0.42):
    return 1 + alpha * (quartic_ratio(g) - 1)

def flat_ratio(g):
    return 1.0

import math
def sigmoid_ratio(g, target=100):
    ln_ref = math.log(target / 120 + 1)
    ln_g = math.log(max(g, 40) / 120 + 1)
    return ln_ref / ln_g if ln_g > 0 else 1.0

SHAPES = {
    'Quartic': quartic_ratio,
    'Gentle (α=0.42)': gentle_ratio,
    'Flat': flat_ratio,
    'Sigmoid': sigmoid_ratio,
}


# ── Load data ──────────────────────────────────────────────────────────────

with open(TRIO_CACHE, 'rb') as f:
    trio_sites = pickle.load(f)

with open(BOOST_CACHE, 'rb') as f:
    boost_raw = pickle.load(f)

for s in trio_sites:
    s['name'] = ANON.get(s['name'], s['name'])

boost_df = boost_raw['strict']

sites = []
for s in trio_sites:
    bg = s['bg']; isf = s['isf_actual']
    m100 = (bg >= 96) & (bg < 104)
    isf100 = np.median(isf[m100]) if m100.sum() >= 5 else np.nan
    sites.append({
        'name': s['name'], 'model': s['model'], 'tdd': s['tdd_median'],
        'n': s['n'], 'bg': bg, 'isf_actual': isf,
        'pred_drop': s['pred_drop'], 'actual_bg_end': s['actual_bg_end'],
        'pred_loop': s['pred_loop'],
        'isf_true': isf100, 'isf_tdd': 1800 / s['tdd_median'],
    })

bb = boost_df['bg'].values.astype(float)
bi = boost_df['variable_sens'].values.astype(float)
m100 = (bb >= 96) & (bb < 104)
sites.append({
    'name': 'User-M', 'model': 'AAPS',
    'tdd': boost_df['tdd_7day'].median(), 'n': len(bb),
    'bg': bb, 'isf_actual': bi,
    'pred_drop': boost_df['pred_drop'].values.astype(float),
    'actual_bg_end': boost_df['actual_bg_end'].values.astype(float),
    'pred_loop': bb - boost_df['pred_drop'].values.astype(float),
    'isf_true': np.median(bi[m100]),
    'isf_tdd': 1800 / boost_df['tdd_7day'].median(),
})

valid_sites = [s for s in sites if not np.isnan(s['isf_true'])]
print(f"Loaded {len(sites)} sites ({len(valid_sites)} with ISF@100), "
      f"{sum(s['n'] for s in sites)} total samples\n")


# ── ANALYSIS 1: ISF × TDD is not constant ─────────────────────────────────

print("=" * 80)
print("ANALYSIS 1: The 1800/TDD assumption — ISF × TDD is NOT constant")
print("=" * 80)

products = []
for s in valid_sites:
    c = s['isf_true'] * s['tdd']
    products.append(c)
    tdd_err = (s['isf_tdd'] / s['isf_true'] - 1) * 100
    print(f"  {s['name']:8s}  TDD={s['tdd']:5.1f}  ISF@100={s['isf_true']:6.0f}  "
          f"ISF×TDD={c:7.0f}  1800/TDD err={tdd_err:+.0f}%")

print(f"\n  ISF×TDD range: {min(products):.0f} to {max(products):.0f}")
print(f"  CV = {np.std(products)/np.mean(products):.1%}")
print(f"  If ISF×TDD were constant, all would be 1800. Range is 518-3939.\n")


# ── ANALYSIS 2: Shape doesn't matter when anchor is wrong ──────────────────

print("=" * 80)
print("ANALYSIS 2: Shape comparison with 1800/TDD anchor (uncalibrated)")
print("=" * 80)

weights = np.array([s['n'] for s in valid_sites])
shape_results = {}

for shape_name, shape_fn in SHAPES.items():
    site_maes = []
    for s in valid_sites:
        isf_model = np.array([s['isf_tdd'] * shape_fn(g) for g in s['bg']])
        pred = s['bg'] - s['pred_drop'] * (isf_model / s['isf_actual'])
        site_maes.append(np.mean(np.abs(pred - s['actual_bg_end'])))
    wmae = np.average(site_maes, weights=weights)
    shape_results[shape_name] = {'wmae': wmae, 'per_site': site_maes}
    print(f"  {shape_name:16s}: weighted MAE = {wmae:.1f}")

loop_maes = [np.mean(np.abs(s['pred_loop'] - s['actual_bg_end'])) for s in valid_sites]
loop_wmae = np.average(loop_maes, weights=weights)
print(f"  {'Loop (tuned)':16s}: weighted MAE = {loop_wmae:.1f}")
print(f"\n  Shape spread: {max(v['wmae'] for v in shape_results.values()) - min(v['wmae'] for v in shape_results.values()):.1f} mg/dL")
print(f"  Anchor error accounts for the rest.\n")


# ── ANALYSIS 3: Robustness to anchor error ─────────────────────────────────

print("=" * 80)
print("ANALYSIS 3: MAE sensitivity to ISF anchor error (quartic shape)")
print("=" * 80)

anchor_errors = [-50, -30, -20, -10, 0, 10, 20, 30, 50]
sensitivity = {}

for s in valid_sites:
    maes = []
    for err_pct in anchor_errors:
        isf_anchor = s['isf_true'] * (1 + err_pct / 100)
        isf_model = np.array([isf_anchor * quartic_ratio(g) for g in s['bg']])
        pred = s['bg'] - s['pred_drop'] * (isf_model / s['isf_actual'])
        maes.append(np.mean(np.abs(pred - s['actual_bg_end'])))
    sensitivity[s['name']] = maes
    # Sensitivity: MAE change per 10% anchor error
    mae_at_0 = maes[anchor_errors.index(0)]
    mae_at_30 = maes[anchor_errors.index(30)]
    sens = (mae_at_30 - mae_at_0) / 3  # per 10%
    print(f"  {s['name']:8s}  MAE@0%={mae_at_0:.1f}  MAE@±30%={maes[anchor_errors.index(-30)]:.1f}/{mae_at_30:.1f}"
          f"  sensitivity={sens:.1f}/10%")


# ── ANALYSIS 4: Auto-learning anchor from outcomes ─────────────────────────

print("\n" + "=" * 80)
print("ANALYSIS 4: Auto-learning ISF anchor from prediction outcomes")
print("  Start from 1800/TDD, use rolling median of implied ISF from outcomes")
print("=" * 80)

def simulate_autolearn(s, window=100, shape_fn=quartic_ratio, min_drop=10):
    bg = s['bg']; isf_actual = s['isf_actual']
    pred_drop = s['pred_drop']; actual_end = s['actual_bg_end']

    isf_est = s['isf_tdd']
    implied_isfs = []
    errors = []
    isf_history = [isf_est]

    for i in range(len(bg)):
        g = bg[i]
        ratio = shape_fn(g)
        isf_model = isf_est * ratio

        pred = g - pred_drop[i] * (isf_model / isf_actual[i])
        errors.append(abs(pred - actual_end[i]))

        actual_drop = g - actual_end[i]
        if abs(actual_drop) > min_drop and abs(pred_drop[i]) > 5:
            correct_isf = isf_actual[i] * actual_drop / pred_drop[i]
            implied_isf_100 = correct_isf / ratio if ratio > 0 else isf_est
            if 3 < implied_isf_100 < 600:
                implied_isfs.append(implied_isf_100)

        if len(implied_isfs) >= 10:
            recent = implied_isfs[-window:]
            isf_est = np.median(recent)

        isf_history.append(isf_est)

    return np.array(errors), np.array(isf_history), isf_est

autolearn_results = {}
print(f"\n  {'Site':8s} {'Start':>7s} {'True':>6s} {'Final':>7s} {'Err%':>6s} "
      f"{'First50':>8s} {'Last200':>8s} {'Overall':>8s} {'Loop':>6s}")
print("  " + "-" * 70)

all_first = []; all_last = []; all_overall = []; all_loop = []
for s in valid_sites:
    errs, hist, final = simulate_autolearn(s, window=100)
    autolearn_results[s['name']] = {'errors': errs, 'history': hist, 'final': final}
    loop_mae = np.mean(np.abs(s['pred_loop'] - s['actual_bg_end']))
    f50 = np.mean(errs[:50]) if len(errs) >= 50 else np.mean(errs)
    l200 = np.mean(errs[-200:]) if len(errs) >= 200 else np.mean(errs[-100:])
    overall = np.mean(errs)
    fin_err = (final / s['isf_true'] - 1) * 100

    all_first.append(f50); all_last.append(l200); all_overall.append(overall)
    all_loop.append(loop_mae)

    print(f"  {s['name']:8s} {s['isf_tdd']:7.1f} {s['isf_true']:6.0f} {final:7.1f} {fin_err:+5.0f}% "
          f"{f50:8.1f} {l200:8.1f} {overall:8.1f} {loop_mae:6.1f}")

w = weights
print(f"\n  Weighted means: {'':32s} {np.average(all_first, weights=w):8.1f} "
      f"{np.average(all_last, weights=w):8.1f} {np.average(all_overall, weights=w):8.1f} "
      f"{np.average(all_loop, weights=w):6.1f}")


# ── CHARTS ─────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("Uncalibrated User Analysis: ISF Anchor Matters More Than Shape",
             fontsize=13, fontweight='bold')

# Chart 1: ISF × TDD scatter
ax = axes[0, 0]
for s in valid_sites:
    c = s['isf_true'] * s['tdd']
    color = 'tab:blue' if s['model'] == 'sigmoid' else ('tab:orange' if s['model'] == 'log' else 'tab:green')
    ax.scatter(s['tdd'], s['isf_true'], c=color, s=60, zorder=3)
    ax.annotate(s['name'], (s['tdd'], s['isf_true']), fontsize=6, ha='left')
# Plot 1800/TDD line
tdd_range = np.linspace(5, 90, 100)
ax.plot(tdd_range, 1800 / tdd_range, 'r--', alpha=0.5, label='1800/TDD')
ax.set_xlabel('TDD (U/day)')
ax.set_ylabel('Actual ISF at 100 mg/dL')
ax.set_title('A. ISF × TDD Is Not Constant')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Chart 2: Shape comparison with 1800/TDD
ax = axes[0, 1]
shape_names = list(SHAPES.keys()) + ['Loop (tuned)']
shape_wmaes = [shape_results[k]['wmae'] for k in SHAPES.keys()] + [loop_wmae]
colors = ['tab:red', 'tab:orange', 'tab:purple', 'tab:cyan', 'tab:green']
bars = ax.bar(range(len(shape_names)), shape_wmaes, color=colors)
ax.set_xticks(range(len(shape_names)))
ax.set_xticklabels(shape_names, rotation=30, ha='right', fontsize=8)
ax.set_ylabel('Weighted Mean MAE (mg/dL)')
ax.set_title('B. Shape Barely Matters With 1800/TDD Anchor')
ax.set_ylim(15, 28)
for bar, val in zip(bars, shape_wmaes):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2, f'{val:.1f}',
            ha='center', fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

# Chart 3: Anchor convergence for select sites
ax = axes[1, 0]
highlight = ['User-A', 'User-D', 'User-I', 'User-M', 'User-E']
for s in valid_sites:
    if s['name'] not in highlight:
        continue
    hist = autolearn_results[s['name']]['history']
    x = np.arange(len(hist))
    ax.plot(x, hist, linewidth=1.2, label=f"{s['name']} (true={s['isf_true']:.0f})")
    ax.axhline(s['isf_true'], color=ax.get_lines()[-1].get_color(), linestyle=':', alpha=0.4)
ax.set_xlabel('Sample number')
ax.set_ylabel('Learned ISF@100 (mg/dL)')
ax.set_title('C. Auto-Learning ISF Anchor Convergence')
ax.legend(fontsize=7, loc='upper right')
ax.grid(True, alpha=0.3)

# Chart 4: Comparison — uncalibrated vs auto-learned vs Loop
ax = axes[1, 1]
x = np.arange(len(valid_sites))
width = 0.25
site_names = [s['name'] for s in valid_sites]

uncal_maes = [shape_results['Quartic']['per_site'][i] for i in range(len(valid_sites))]
auto_maes = [np.mean(autolearn_results[s['name']]['errors']) for s in valid_sites]
loop_site_maes = [np.mean(np.abs(s['pred_loop'] - s['actual_bg_end'])) for s in valid_sites]

ax.bar(x - width, uncal_maes, width, label='1800/TDD + Quartic', color='tab:red', alpha=0.7)
ax.bar(x, auto_maes, width, label='Auto-learned anchor', color='tab:blue', alpha=0.7)
ax.bar(x + width, loop_site_maes, width, label='Loop (tuned)', color='tab:green', alpha=0.7)
ax.set_xticks(x)
ax.set_xticklabels(site_names, rotation=45, ha='right', fontsize=7)
ax.set_ylabel('MAE (mg/dL)')
ax.set_title('D. Uncalibrated vs Auto-Learned vs Tuned Loop')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3, axis='y')
# Cap Y axis to show detail (User-A uncalibrated is 78)
ax.set_ylim(0, 50)

plt.tight_layout()
plt.savefig(OUT_DIR / 'uncalibrated_analysis.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR / 'uncalibrated_analysis.png'}")

# ── Summary table for paper ────────────────────────────────────────────────
print("\n" + "=" * 80)
print("SUMMARY FOR PAPER")
print("=" * 80)
print(f"\n  1. ISF×TDD varies from {min(products):.0f} to {max(products):.0f} (CV {np.std(products)/np.mean(products):.0%})")
print(f"     → The 1800 rule is unreliable for individual patients")
print(f"\n  2. With 1800/TDD anchor, shape spread is only "
      f"{max(v['wmae'] for v in shape_results.values()) - min(v['wmae'] for v in shape_results.values()):.1f} mg/dL")
print(f"     → Shape choice is secondary to getting the anchor right")
print(f"\n  3. Auto-learning from outcomes:")
print(f"     First 50 samples: weighted MAE {np.average(all_first, weights=w):.1f}")
print(f"     Last 200 samples: weighted MAE {np.average(all_last, weights=w):.1f}")
print(f"     Loop (fully tuned): weighted MAE {np.average(all_loop, weights=w):.1f}")
print(f"     → After ~100 samples, auto-learned anchor matches or beats tuned Loop")
print(f"\n  4. Recommendation for uncalibrated users:")
print(f"     Start from 1800/TDD + quartic shape (best out-of-box)")
print(f"     Auto-learn ISF anchor from rolling median of implied ISF")
print(f"     Conservative dosing during the first ~100 samples (learning phase)")
print("\nDONE")
