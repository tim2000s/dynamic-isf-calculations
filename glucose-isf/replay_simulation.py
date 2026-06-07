#!/usr/bin/env python3
"""
Closed-Loop Replay Simulation
==============================
Replays actual fasting windows and simulates what would have happened to
glucose if the AID system had used a different ISF model.  Unlike the
counterfactual method (which only changes the prediction), this propagates
altered insulin delivery through a simple glucose model.

Core logic:
    1. Infer IOB from the system's own prediction: IOB = pred_drop / isf_actual
    2. Compute what the system would have predicted with the test ISF
    3. Compute the difference in correction insulin delivery
    4. Apply that delta-insulin to the observed glucose outcome

Data: allday fasting caches (COB=0, no recent bolus, BG 72-200).
"""

import pickle, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
OUT_DIR = Path('/Users/tims/Downloads/4 Hour analysis')
CACHE_DIR = OUT_DIR / 'daytime analysis'
TRIO_CACHE = CACHE_DIR / 'multisite_allday_cache.pkl'
BOOST_CACHE = CACHE_DIR / 'boost_allday_cache.pkl'

# ── Anonymisation ─────────────────────────────────────────────────────────────
ANON_MAP = {
    'henny425': 'User-A', 'aadiabetes': 'User-B', 'diajesse': 'User-C',
    'svns': 'User-D', 'fuxchr': 'User-E', 'mikens': 'User-F',
    'andycgm': 'User-G', 'noahr': 'User-H', 'nightscout1': 'User-I',
    'eli': 'User-J', 'ns_rot6': 'User-K', 'kelseyhuss': 'User-L',
}

# ── ISF model functions ──────────────────────────────────────────────────────

def quartic(g):
    return 272 - 3.121*g + 0.01511*g**2 - 3.305e-5*g**3 + 2.69e-8*g**4

Q_REF = quartic(100)

RATIOS_POP = {76: 1.15, 100: 1.00, 130: 0.80, 170: 0.70}

def ratio_fn(g, ratios=None):
    if ratios is None:
        ratios = RATIOS_POP
    points = sorted(ratios.items())
    gs = [p[0] for p in points]; rs = [p[1] for p in points]
    if g <= gs[0]: return rs[0]
    if g >= gs[-1]: return rs[-1]
    for i in range(len(gs) - 1):
        if gs[i] <= g <= gs[i+1]:
            t = (g - gs[i]) / (gs[i+1] - gs[i])
            return rs[i] + t * (rs[i+1] - rs[i])
    return 1.0

def sigmoid_ratio(g, target=100):
    ln_ref = math.log(target / 120 + 1)
    ln_g = math.log(max(g, 40) / 120 + 1)
    return ln_ref / ln_g if ln_g > 0 else 1.0

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading caches ...")
with open(TRIO_CACHE, 'rb') as f:
    trio_sites = pickle.load(f)
with open(BOOST_CACHE, 'rb') as f:
    boost_cache = pickle.load(f)

# ── Build unified site list ──────────────────────────────────────────────────

def build_sites(trio_sites, boost_cache, period='allday'):
    sites = []
    for s in trio_sites:
        data = s.get(period)
        if data is None or data['n'] < 10:
            continue
        bg = data['bg']; isf = data['isf_actual']
        m100 = (bg >= 96) & (bg < 104)
        isf100 = np.median(isf[m100]) if m100.sum() >= 5 else np.nan
        name = ANON_MAP.get(s['name'], s['name'])
        sites.append({
            'name': name, 'model': s['model'],
            'tdd': s['tdd_median'], 'n': data['n'],
            'bg': bg, 'isf_actual': isf,
            'pred_drop': data['pred_drop'],
            'actual_bg_end': data['actual_bg_end'],
            'hour': data.get('hour'),
            'isf_true': isf100,
            'isf_tdd': 1800 / s['tdd_median'],
        })

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
            'hour': boost_df['hour'].values.astype(int),
            'isf_true': isf100, 'isf_tdd': 1800 / tdd,
        })
    return sites


# ══════════════════════════════════════════════════════════════════════════════
# REPLAY SIMULATION ENGINE
# ══════════════════════════════════════════════════════════════════════════════

TARGET = 100  # mg/dL — typical AID target

def replay_site(site, isf_test_arr):
    """
    For each fasting sample, simulate what glucose would have been if the
    system had used isf_test instead of isf_actual.

    Returns array of simulated end-glucose values.
    """
    bg = site['bg']
    isf_actual = site['isf_actual']
    pred_drop = site['pred_drop']
    actual_bg_end = site['actual_bg_end']

    # Infer effective IOB from system's prediction
    iob = pred_drop / isf_actual

    # What the system predicted with its own ISF
    eventual_orig = bg - iob * isf_actual  # = bg - pred_drop

    # What the system would predict with the test ISF
    eventual_test = bg - iob * isf_test_arr

    # Correction insulin: (eventual - target) / ISF
    # Positive = deliver more insulin; negative = suspend / reduce
    corr_orig = (eventual_orig - TARGET) / isf_actual
    corr_test = (eventual_test - TARGET) / isf_test_arr

    # Delta insulin delivered
    delta_insulin = corr_test - corr_orig

    # Clip to +-2U (system rate limits over ~4h horizon)
    delta_insulin = np.clip(delta_insulin, -2.0, 2.0)

    # Estimate true ISF from observed outcomes
    actual_drop = bg - actual_bg_end
    with np.errstate(divide='ignore', invalid='ignore'):
        true_isf = np.where(np.abs(iob) > 0.01,
                            actual_drop / iob,
                            isf_actual)
    true_isf = np.clip(true_isf, 5, 500)

    # Simulated end glucose
    simulated_end = actual_bg_end - delta_insulin * true_isf

    # Physiological clamp
    simulated_end = np.clip(simulated_end, 40, 400)

    return simulated_end


# ── ISF model definitions ────────────────────────────────────────────────────
# Each returns an array of ISF values for the site's BG array

MODEL_DEFS = {
    'Original (validation)': lambda s: s['isf_actual'],
    'Quartic + 1800/TDD': lambda s: np.array([s['isf_tdd'] * quartic(g) / Q_REF for g in s['bg']]),
    'Profile + quartic': lambda s: np.array([s['isf_true'] * quartic(g) / Q_REF for g in s['bg']]),
    'Profile + pop-ratio': lambda s: np.array([s['isf_true'] * ratio_fn(g) for g in s['bg']]),
    'Profile + sigmoid': lambda s: np.array([s['isf_true'] * sigmoid_ratio(g) for g in s['bg']]),
    'Profile + flat': lambda s: np.full(len(s['bg']), s['isf_true']),
    '1800/TDD + flat': lambda s: np.full(len(s['bg']), s['isf_tdd']),
}


# ── Run simulation ───────────────────────────────────────────────────────────

sites = build_sites(trio_sites, boost_cache, 'allday')
# Filter to sites with valid isf_true
sites = [s for s in sites if not np.isnan(s['isf_true'])]

print(f"\nLoaded {len(sites)} sites with valid ISF@100")
total_n = sum(s['n'] for s in sites)
print(f"Total fasting samples: {total_n:,}\n")

# Store results: {model_name: {site_name: simulated_end_array}}
results = {}

for model_name, model_fn in MODEL_DEFS.items():
    results[model_name] = {}
    for s in sites:
        isf_test = model_fn(s)
        sim_end = replay_site(s, isf_test)
        results[model_name][s['name']] = sim_end


# ══════════════════════════════════════════════════════════════════════════════
# COMPUTE METRICS
# ══════════════════════════════════════════════════════════════════════════════

def compute_metrics(bg_arr):
    """Compute TIR, TITR (normoglycaemia), time<70, time<54, time>140, mean, SD."""
    n = len(bg_arr)
    if n == 0:
        return {'tir': np.nan, 'titr': np.nan, 'below70': np.nan, 'below54': np.nan,
                'above140': np.nan, 'mean': np.nan, 'sd': np.nan, 'n': 0}
    tir = np.mean((bg_arr >= 70) & (bg_arr <= 180)) * 100
    titr = np.mean((bg_arr >= 70) & (bg_arr <= 140)) * 100   # normoglycaemia
    below70 = np.mean(bg_arr < 70) * 100
    below54 = np.mean(bg_arr < 54) * 100
    above140 = np.mean(bg_arr > 140) * 100
    return {
        'tir': tir, 'titr': titr, 'below70': below70, 'below54': below54,
        'above140': above140, 'mean': np.mean(bg_arr), 'sd': np.std(bg_arr), 'n': n,
    }


# Actual outcomes (baseline)
actual_metrics = {}
for s in sites:
    actual_metrics[s['name']] = compute_metrics(s['actual_bg_end'])

# Simulated outcomes per model
sim_metrics = {}  # {model: {site: metrics_dict}}
for model_name in MODEL_DEFS:
    sim_metrics[model_name] = {}
    for s in sites:
        sim_metrics[model_name][s['name']] = compute_metrics(results[model_name][s['name']])


# Weighted aggregation
def weighted_aggregate(metric_dict, sites):
    """Compute sample-weighted averages across sites."""
    weights = np.array([s['n'] for s in sites])
    names = [s['name'] for s in sites]
    agg = {}
    for key in ['tir', 'titr', 'below70', 'below54', 'above140', 'mean', 'sd']:
        vals = np.array([metric_dict[n][key] for n in names])
        agg[key] = np.average(vals, weights=weights)
    agg['n'] = sum(weights)
    return agg

actual_agg = weighted_aggregate(actual_metrics, sites)


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATION: Original model should closely match actual
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 90)
print("VALIDATION: Original ISF replay vs actual outcomes")
print("=" * 90)

orig_agg = weighted_aggregate(sim_metrics['Original (validation)'], sites)
print(f"\n  {'Metric':20s} {'Actual':>10s} {'Sim (orig)':>10s} {'Delta':>10s}")
print("  " + "-" * 55)
for key in ['tir', 'titr', 'below70', 'below54', 'above140', 'mean', 'sd']:
    fmt = '.1f' if key != 'n' else 'd'
    suffix = '%' if key in ('tir', 'titr', 'below70', 'below54', 'above140') else ' mg/dL'
    a = actual_agg[key]; o = orig_agg[key]
    print(f"  {key:20s} {a:10.1f}{suffix:>6s} {o:10.1f}{suffix:>6s} {o-a:+10.1f}")

print("\n  Per-site validation (TIR):")
print(f"  {'Site':8s} {'Actual TIR':>10s} {'Sim TIR':>10s} {'Delta':>8s}")
print("  " + "-" * 40)
for s in sites:
    a_tir = actual_metrics[s['name']]['tir']
    s_tir = sim_metrics['Original (validation)'][s['name']]['tir']
    print(f"  {s['name']:8s} {a_tir:10.1f}% {s_tir:10.1f}% {s_tir-a_tir:+8.1f}")


# ══════════════════════════════════════════════════════════════════════════════
# MODEL COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("MODEL COMPARISON: Simulated glucose outcomes")
print("=" * 90)

model_aggs = {}
for model_name in MODEL_DEFS:
    model_aggs[model_name] = weighted_aggregate(sim_metrics[model_name], sites)

print(f"\n  {'Model':25s} {'TIR%':>8s} {'TITR%':>8s} {'>140%':>8s} {'<70%':>8s} {'<54%':>8s} {'Mean':>8s} {'SD':>8s}")
print("  " + "-" * 86)
# Print actual baseline first
print(f"  {'** Actual outcomes **':25s} {actual_agg['tir']:8.1f} {actual_agg['titr']:8.1f} "
      f"{actual_agg['above140']:8.1f} {actual_agg['below70']:8.1f} "
      f"{actual_agg['below54']:8.1f} {actual_agg['mean']:8.1f} {actual_agg['sd']:8.1f}")
print("  " + "-" * 86)
for model_name in MODEL_DEFS:
    a = model_aggs[model_name]
    print(f"  {model_name:25s} {a['tir']:8.1f} {a['titr']:8.1f} "
          f"{a['above140']:8.1f} {a['below70']:8.1f} "
          f"{a['below54']:8.1f} {a['mean']:8.1f} {a['sd']:8.1f}")


# ── Per-site breakdown for top models ─────────────────────────────────────────

TOP_MODELS = ['Original (validation)', 'Quartic + 1800/TDD',
              'Profile + quartic', 'Profile + pop-ratio']

print("\n\n" + "=" * 90)
print("PER-SITE TIR: Top models")
print("=" * 90)

header = f"  {'Site':8s} {'Actual':>8s}"
for m in TOP_MODELS:
    short = m[:12]
    header += f" {short:>12s}"
print(header)
print("  " + "-" * (8 + 8 + 12 * len(TOP_MODELS) + len(TOP_MODELS)))

for s in sites:
    line = f"  {s['name']:8s} {actual_metrics[s['name']]['tir']:8.1f}"
    for m in TOP_MODELS:
        tir = sim_metrics[m][s['name']]['tir']
        delta = tir - actual_metrics[s['name']]['tir']
        line += f" {tir:7.1f}({delta:+.0f})"
    print(line)


# ── Per-site safety ──────────────────────────────────────────────────────────

print("\n\n" + "=" * 90)
print("PER-SITE SAFETY: Time below 70 mg/dL")
print("=" * 90)

header = f"  {'Site':8s} {'Actual':>8s}"
for m in TOP_MODELS:
    short = m[:12]
    header += f" {short:>12s}"
print(header)
print("  " + "-" * (8 + 8 + 12 * len(TOP_MODELS) + len(TOP_MODELS)))

for s in sites:
    line = f"  {s['name']:8s} {actual_metrics[s['name']]['below70']:8.1f}"
    for m in TOP_MODELS:
        b70 = sim_metrics[m][s['name']]['below70']
        delta = b70 - actual_metrics[s['name']]['below70']
        line += f" {b70:7.1f}({delta:+.0f})"
    print(line)


# ══════════════════════════════════════════════════════════════════════════════
# CHARTS
# ══════════════════════════════════════════════════════════════════════════════

MODEL_COLORS = {
    'Original (validation)': '#555555',
    'Quartic + 1800/TDD': '#2196F3',
    'Profile + quartic': '#4CAF50',
    'Profile + pop-ratio': '#FF9800',
    'Profile + sigmoid': '#9C27B0',
    'Profile + flat': '#F44336',
    '1800/TDD + flat': '#795548',
}
MODEL_SHORT = {
    'Original (validation)': 'Original',
    'Quartic + 1800/TDD': 'Q+TDD',
    'Profile + quartic': 'Prof+Q',
    'Profile + pop-ratio': 'Prof+Pop',
    'Profile + sigmoid': 'Prof+Sig',
    'Profile + flat': 'Prof+Flat',
    '1800/TDD + flat': 'TDD+Flat',
}


# ── Chart 1: TIR and TITR comparison ─────────────────────────────────────────

fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(16, 5.5))
fig1.suptitle("Closed-Loop Replay: Simulated Glycaemic Outcomes\n"
              "Fasting samples only, weighted across patients",
              fontsize=12, fontweight='bold')

model_names = list(MODEL_DEFS.keys())
x = np.arange(len(model_names) + 1)  # +1 for actual baseline
labels = ['Actual'] + [MODEL_SHORT[m] for m in model_names]
colors = ['#333333'] + [MODEL_COLORS[m] for m in model_names]

# Panel A: TIR (70-180)
tir_vals = [actual_agg['tir']] + [model_aggs[m]['tir'] for m in model_names]
bars_a = ax1a.bar(x, tir_vals, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5)
for bar, val in zip(bars_a, tir_vals):
    ax1a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{val:.1f}%', ha='center', fontsize=8, fontweight='bold')
ax1a.set_xticks(x)
ax1a.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
ax1a.set_ylabel('Time in Range (%)')
ax1a.set_title('A. Time in Range (70–180 mg/dL)')
ax1a.set_ylim(min(tir_vals) - 5, max(tir_vals) + 3)
ax1a.axhline(actual_agg['tir'], color='#333333', linestyle='--', alpha=0.4, linewidth=1)
ax1a.grid(True, alpha=0.2, axis='y')

# Panel B: TITR (70-140) — normoglycaemia
titr_vals = [actual_agg['titr']] + [model_aggs[m]['titr'] for m in model_names]
bars_b = ax1b.bar(x, titr_vals, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5)
for bar, val in zip(bars_b, titr_vals):
    ax1b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
             f'{val:.1f}%', ha='center', fontsize=8, fontweight='bold')
ax1b.set_xticks(x)
ax1b.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
ax1b.set_ylabel('Time in Tight Range (%)')
ax1b.set_title('B. Time in Normoglycaemia (70–140 mg/dL)')
ax1b.set_ylim(min(titr_vals) - 5, max(titr_vals) + 3)
ax1b.axhline(actual_agg['titr'], color='#333333', linestyle='--', alpha=0.4, linewidth=1)
ax1b.grid(True, alpha=0.2, axis='y')

plt.tight_layout()
plt.savefig(OUT_DIR / 'replay_tir_comparison.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR / 'replay_tir_comparison.png'}")


# ── Chart 2: Safety — Time below 70 and below 54 ─────────────────────────────

fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(12, 5))
fig2.suptitle("Closed-Loop Replay: Hypoglycaemia Risk by ISF Model\n"
              "Fasting samples, weighted across patients",
              fontsize=11, fontweight='bold')

# Time below 70
labels_m = [MODEL_SHORT[m] for m in model_names]
b70_vals = [model_aggs[m]['below70'] for m in model_names]
b70_colors = [MODEL_COLORS[m] for m in model_names]
bars2a = ax2a.bar(range(len(model_names)), b70_vals, color=b70_colors, alpha=0.85)
ax2a.axhline(actual_agg['below70'], color='red', linestyle='--', alpha=0.5, label='Actual')
for bar, val in zip(bars2a, b70_vals):
    ax2a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
             f'{val:.2f}%', ha='center', fontsize=7)
ax2a.set_xticks(range(len(model_names)))
ax2a.set_xticklabels(labels_m, rotation=35, ha='right', fontsize=8)
ax2a.set_ylabel('Time < 70 mg/dL (%)')
ax2a.set_title('A. Time Below Range')
ax2a.legend(fontsize=8)
ax2a.grid(True, alpha=0.2, axis='y')

# Time below 54
b54_vals = [model_aggs[m]['below54'] for m in model_names]
bars2b = ax2b.bar(range(len(model_names)), b54_vals, color=b70_colors, alpha=0.85)
ax2b.axhline(actual_agg['below54'], color='red', linestyle='--', alpha=0.5, label='Actual')
for bar, val in zip(bars2b, b54_vals):
    ax2b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
             f'{val:.2f}%', ha='center', fontsize=7)
ax2b.set_xticks(range(len(model_names)))
ax2b.set_xticklabels(labels_m, rotation=35, ha='right', fontsize=8)
ax2b.set_ylabel('Time < 54 mg/dL (%)')
ax2b.set_title('B. Clinically Significant Hypo')
ax2b.legend(fontsize=8)
ax2b.grid(True, alpha=0.2, axis='y')

plt.tight_layout()
plt.savefig(OUT_DIR / 'replay_safety.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'replay_safety.png'}")


# ── Chart 3: Per-site TIR heatmap ────────────────────────────────────────────

display_models = TOP_MODELS
site_names = [s['name'] for s in sites]

fig3, ax3 = plt.subplots(figsize=(10, max(5, len(sites) * 0.45 + 1)))
fig3.suptitle("Closed-Loop Replay: Per-Site TIR Change vs Actual\n"
              "(Positive = improvement, Negative = worse)",
              fontsize=11, fontweight='bold')

# Build delta-TIR matrix
delta_matrix = np.zeros((len(sites), len(display_models)))
for j, m in enumerate(display_models):
    for i, s in enumerate(sites):
        delta_matrix[i, j] = (sim_metrics[m][s['name']]['tir'] -
                              actual_metrics[s['name']]['tir'])

im = ax3.imshow(delta_matrix, cmap='RdYlGn', aspect='auto',
                vmin=-5, vmax=5)
ax3.set_xticks(range(len(display_models)))
ax3.set_xticklabels([MODEL_SHORT[m] for m in display_models], rotation=35, ha='right', fontsize=9)
ax3.set_yticks(range(len(sites)))
ax3.set_yticklabels(site_names, fontsize=8)

# Annotate cells
for i in range(len(sites)):
    for j in range(len(display_models)):
        val = delta_matrix[i, j]
        color = 'white' if abs(val) > 3 else 'black'
        ax3.text(j, i, f'{val:+.1f}', ha='center', va='center', fontsize=7, color=color)

plt.colorbar(im, ax=ax3, label='TIR change (pp)')
plt.tight_layout()
plt.savefig(OUT_DIR / 'replay_per_site.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'replay_per_site.png'}")


# ── Chart 4: Glucose distribution violin/box plots ───────────────────────────

fig4, ax4 = plt.subplots(figsize=(12, 5))
fig4.suptitle("Closed-Loop Replay: Simulated End-Glucose Distribution\n"
              "Pooled across all patients (fasting samples)",
              fontsize=11, fontweight='bold')

# Collect pooled distributions
all_models_for_violin = ['Original (validation)', 'Quartic + 1800/TDD',
                         'Profile + quartic', 'Profile + pop-ratio',
                         'Profile + flat', '1800/TDD + flat']

# Start with actual
violin_data = [np.concatenate([s['actual_bg_end'] for s in sites])]
violin_labels = ['Actual']
violin_colors = ['#333333']

for m in all_models_for_violin:
    pooled = np.concatenate([results[m][s['name']] for s in sites])
    violin_data.append(pooled)
    violin_labels.append(MODEL_SHORT[m])
    violin_colors.append(MODEL_COLORS[m])

positions = range(len(violin_data))
parts = ax4.violinplot(violin_data, positions=positions, showmedians=True,
                       showextrema=False)

for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(violin_colors[i])
    pc.set_alpha(0.6)
parts['cmedians'].set_color('black')

# Add IQR boxes
for i, data in enumerate(violin_data):
    q1, med, q3 = np.percentile(data, [25, 50, 75])
    ax4.vlines(i, q1, q3, color='black', linewidth=3, alpha=0.5)

# Range lines
ax4.axhline(70, color='red', linestyle='--', alpha=0.4, linewidth=1, label='Hypo (70)')
ax4.axhline(180, color='orange', linestyle='--', alpha=0.4, linewidth=1, label='Hyper (180)')
ax4.axhline(TARGET, color='green', linestyle=':', alpha=0.4, linewidth=1, label='Target (100)')

ax4.set_xticks(positions)
ax4.set_xticklabels(violin_labels, rotation=30, ha='right', fontsize=9)
ax4.set_ylabel('End Glucose (mg/dL)')
ax4.set_ylim(40, 300)
ax4.legend(fontsize=8, loc='upper right')
ax4.grid(True, alpha=0.2, axis='y')

plt.tight_layout()
plt.savefig(OUT_DIR / 'replay_glucose_distribution.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'replay_glucose_distribution.png'}")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("SUMMARY: Replay Simulation Results")
print("=" * 90)

# Rank models by TIR improvement
ranked = sorted(MODEL_DEFS.keys(),
                key=lambda m: model_aggs[m]['tir'], reverse=True)

print(f"\n  Models ranked by simulated TIR (actual baseline: {actual_agg['tir']:.1f}%):\n")
for i, m in enumerate(ranked):
    a = model_aggs[m]
    delta_tir = a['tir'] - actual_agg['tir']
    delta_titr = a['titr'] - actual_agg['titr']
    delta_b70 = a['below70'] - actual_agg['below70']
    safe_flag = ' [!HYPO]' if delta_b70 > 0.5 else ''
    print(f"    {i+1}. {m:25s}  TIR={a['tir']:.1f}% ({delta_tir:+.1f})  "
          f"TITR={a['titr']:.1f}% ({delta_titr:+.1f})  "
          f"<70={a['below70']:.2f}% ({delta_b70:+.2f}){safe_flag}  "
          f"Mean={a['mean']:.0f}  SD={a['sd']:.0f}")

print("\n\nDONE")
