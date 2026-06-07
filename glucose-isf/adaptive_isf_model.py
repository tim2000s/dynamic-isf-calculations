#!/usr/bin/env python3
"""
Adaptive ISF Model: Two-parameter dynamic ISF using population polynomial + per-patient shape.

Model: ISF(G) = S × [1 + α × (Q_norm(G) - 1)]
  S = 1800 / TDD           (magnitude, from TDD)
  Q_norm(G) = quartic(G) / quartic(G_ref)  (population shape, normalised to 1.0 at G_ref)
  α = per-patient shape factor (0 = flat, 1 = full quartic)

Simulations:
  1. Optimal α per site (best achievable)
  2. Population mean α (one-size-fits-all shape dampening)
  3. Adaptive α (start at 0.5, learn from prediction errors)
  4. Comparison against Loop, Quartic, Full Diabeloop

All 13 subjects: 12 Trio + 1 AAPS/Boost.
"""

import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

OUT_DIR = Path.home() / 'Downloads' / '4 Hour analysis'
TRIO_CACHE = OUT_DIR / 'multisite_4h_sample_cache.pkl'
BOOST_CACHE = OUT_DIR / 'boost_4h_cache.pkl'

BG_BINS = list(range(72, 201, 8))
BG_BIN_CENTERS = [(BG_BINS[i] + BG_BINS[i + 1]) / 2 for i in range(len(BG_BINS) - 1)]
G_REF = 100  # Reference glucose for normalisation

# ── Formulas ────────────────────────────────────────────────────────────────

def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

def full_diabeloop(g):
    return (98.03 - 1.077 * g + 0.008868 * g**2) if g <= 100 else quartic(g)

Q_REF = quartic(G_REF)

def quartic_norm(g):
    """Normalised quartic: 1.0 at G_REF."""
    return quartic(g) / Q_REF

def adaptive_isf(g, S, alpha):
    """Two-parameter adaptive ISF."""
    return S * Q_REF * (1 + alpha * (quartic_norm(g) - 1))


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

# Build unified site list
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

# Boost user
boost_bg = boost_df['bg'].values.astype(float)
boost_isf = boost_df['variable_sens'].values.astype(float)
boost_pred_drop = boost_df['pred_drop'].values.astype(float)
boost_actual_end = boost_df['actual_bg_end'].values.astype(float)
boost_tdd = boost_df['tdd_7day'].median()
# Loop prediction for Boost = current_bg - pred_drop (its own model)
boost_pred_loop = boost_bg - boost_pred_drop

sites.append({
    'name': 'User-M',
    'model': 'AAPS',
    'n': len(boost_bg),
    'tdd': boost_tdd,
    'bg': boost_bg,
    'isf_actual': boost_isf,
    'pred_drop': boost_pred_drop,
    'actual_bg_end': boost_actual_end,
    'pred_loop': boost_pred_loop,
})

print(f"Loaded {len(sites)} sites, {sum(s['n'] for s in sites)} total samples\n")


# ── Compute normalised ISF slopes and optimal α ────────────────────────────

def bin_isf(bg, isf):
    medians = []
    for i in range(len(BG_BINS) - 1):
        mask = (bg >= BG_BINS[i]) & (bg < BG_BINS[i + 1])
        n = mask.sum()
        medians.append(np.median(isf[mask]) if n >= 5 else np.nan)
    return np.array(medians)

# Quartic normalised slope
q_x = np.array(BG_BIN_CENTERS)
q_y = np.array([quartic_norm(g) for g in q_x])
Q_SLOPE = stats.linregress(q_x, q_y).slope

print(f"Quartic normalised slope: {Q_SLOPE*1000:.2f} ×10⁻³/mg/dL\n")

ref_idx = None
for i in range(len(BG_BINS) - 1):
    if BG_BINS[i] <= G_REF < BG_BINS[i + 1]:
        ref_idx = i
        break

for s in sites:
    medians = bin_isf(s['bg'], s['isf_actual'])
    ref_isf = medians[ref_idx]

    if np.isnan(ref_isf) or ref_isf <= 0 or s['name'] == 'User-I':
        s['slope'] = np.nan
        s['alpha_optimal'] = np.nan
        s['ref_isf'] = ref_isf if not np.isnan(ref_isf) else np.nan
        s['normalised'] = np.full_like(medians, np.nan)
        reason = 'constant ISF' if s['name'] == 'User-I' else 'no ref ISF'
        print(f"  {s['name']:15s}  EXCLUDED from shape analysis ({reason})")
        continue

    normalised = medians / ref_isf
    valid = ~np.isnan(normalised)
    if valid.sum() < 4:
        s['slope'] = np.nan
        s['alpha_optimal'] = np.nan
        s['ref_isf'] = ref_isf
        s['normalised'] = normalised
        print(f"  {s['name']:15s}  EXCLUDED (insufficient bins: {valid.sum()})")
        continue

    x = np.array(BG_BIN_CENTERS)[valid]
    y = normalised[valid]
    slope = stats.linregress(x, y).slope
    alpha = slope / Q_SLOPE

    s['slope'] = slope
    s['alpha_optimal'] = alpha
    s['ref_isf'] = ref_isf
    s['normalised'] = normalised

    print(f"  {s['name']:15s}  model={s['model']:7s}  TDD={s['tdd']:5.1f}  "
          f"slope={slope*1000:+6.2f}  α_opt={alpha:+5.2f}  ISF@100={ref_isf:.0f}")

evaluable = [s for s in sites if not np.isnan(s.get('alpha_optimal', np.nan))]
all_alphas = [s['alpha_optimal'] for s in evaluable]
ALPHA_MEAN = np.mean(all_alphas)
ALPHA_MEDIAN = np.median(all_alphas)

print(f"\nEvaluable sites: {len(evaluable)}")
print(f"α range: {min(all_alphas):.2f} to {max(all_alphas):.2f}")
print(f"α mean:  {ALPHA_MEAN:.2f}, median: {ALPHA_MEDIAN:.2f}, std: {np.std(all_alphas):.2f}")


# ── Counterfactual predictions with adaptive model ──────────────────────────

def compute_bias(sites, alpha_mode='optimal'):
    """Compute prediction bias for the adaptive model at different α settings."""
    results = []
    for s in sites:
        bg = s['bg']
        isf_actual = s['isf_actual']
        pred_drop = s['pred_drop']
        actual_end = s['actual_bg_end']

        # TDD scaling
        anchor = quartic(99)
        S = (1800 / s['tdd']) / anchor

        # Choose α
        if alpha_mode == 'optimal':
            alpha = s.get('alpha_optimal', ALPHA_MEAN)
            if np.isnan(alpha):
                alpha = ALPHA_MEAN
        elif alpha_mode == 'mean':
            alpha = ALPHA_MEAN
        elif alpha_mode == 'quartic':  # α=1.0, original quartic
            alpha = 1.0
        elif alpha_mode == 'flat':  # α=0.0, flat ISF
            alpha = 0.0
        else:
            alpha = float(alpha_mode)

        # Compute adaptive ISF for each sample
        isf_adaptive = np.array([adaptive_isf(g, S, alpha) for g in bg])

        # Counterfactual prediction
        pred_adaptive = bg - pred_drop * (isf_adaptive / isf_actual)

        # Bias by zone
        zones = {
            '<105 falling': (bg < 105) & (pred_drop > 0),
            '<105 rising': (bg < 105) & (pred_drop <= 0),
            '>=105 falling': (bg >= 105) & (pred_drop > 0),
            '>=105 rising': (bg >= 105) & (pred_drop <= 0),
        }

        row = {'name': s['name'], 'model': s['model'], 'n': s['n'],
               'tdd': s['tdd'], 'alpha': alpha}
        for zone, mask in zones.items():
            n = mask.sum()
            if n > 0:
                bias = np.mean(pred_adaptive[mask] - actual_end[mask])
                mae = np.mean(np.abs(pred_adaptive[mask] - actual_end[mask]))
            else:
                bias = np.nan
                mae = np.nan
            row[f'{zone}_n'] = n
            row[f'{zone}_bias'] = bias
            row[f'{zone}_mae'] = mae

        # Overall
        bias_all = np.mean(pred_adaptive - actual_end)
        mae_all = np.mean(np.abs(pred_adaptive - actual_end))
        row['overall_bias'] = bias_all
        row['overall_mae'] = mae_all

        # Loop for comparison
        pred_loop = s['pred_loop']
        row['loop_bias'] = np.mean(pred_loop - actual_end)
        row['loop_mae'] = np.mean(np.abs(pred_loop - actual_end))

        for zone, mask in zones.items():
            n = mask.sum()
            if n > 0:
                row[f'{zone}_loop_bias'] = np.mean(pred_loop[mask] - actual_end[mask])
            else:
                row[f'{zone}_loop_bias'] = np.nan

        results.append(row)
    return results


# ── Adaptive learning simulation ────────────────────────────────────────────

def simulate_adaptive(s, alpha_init=0.5, learning_rate=0.01, window=50):
    """
    Simulate adaptive α learning on chronologically ordered data.
    After each batch of `window` samples, adjust α based on observed bias.

    Strategy: If predictions are systematically high at low BG (α too high → ISF too high
    at low glucose → predicts too much drop → predicts too low... wait, let me think about
    the direction carefully.

    ISF higher → predicts larger drop → predicted SGV lower.
    If predicted < actual (negative bias) at low BG, ISF is too high at low BG → α too high.
    If predicted > actual (positive bias) at low BG, ISF is too low at low BG → α too low.

    At high BG, ISF is lower when α is higher.
    If predicted < actual at high BG, ISF too low at high BG → need higher ISF → lower α.
    If predicted > actual at high BG, ISF too high at high BG → need lower ISF → higher α.

    Simpler approach: adjust α based on the correlation between prediction error and glucose.
    If error trends negative with glucose (over-predicts at low, under-predicts at high),
    α is too high. If error trends positive with glucose, α is too low.
    """
    bg = s['bg']
    isf_actual = s['isf_actual']
    pred_drop = s['pred_drop']
    actual_end = s['actual_bg_end']
    n = len(bg)

    anchor = quartic(99)
    S = (1800 / s['tdd']) / anchor

    alpha = alpha_init
    alpha_history = []
    bias_history = []
    pred_all = np.zeros(n)

    for i in range(n):
        g = bg[i]
        isf_a = adaptive_isf(g, S, alpha)
        pred_sgv = g - pred_drop[i] * (isf_a / isf_actual[i])
        pred_all[i] = pred_sgv
        alpha_history.append(alpha)

        # Update α every `window` samples
        if (i + 1) % window == 0 and i >= window:
            recent_bg = bg[i - window + 1:i + 1]
            recent_pred = pred_all[i - window + 1:i + 1]
            recent_actual = actual_end[i - window + 1:i + 1]
            recent_error = recent_pred - recent_actual

            # Compute slope of error vs glucose
            if len(recent_bg) > 10:
                slope, _, r, p, _ = stats.linregress(recent_bg, recent_error)
                # If error trends negative with glucose → α too high → decrease
                # If error trends positive with glucose → α too low → increase
                # The quartic_norm decreases with glucose, so:
                # Higher α → more negative ISF slope → more negative error slope
                # We want error slope = 0
                alpha_adjustment = slope * learning_rate * 100  # scale factor
                alpha = np.clip(alpha + alpha_adjustment, -1.0, 2.0)

            # Also adjust based on overall bias
            mean_bias = np.mean(recent_error)
            # If overall bias positive (predicting too high) and mostly low BG → α too low
            # Actually, overall bias is more about S than α. Skip this.

        bias_history.append(pred_all[i] - actual_end[i])

    return {
        'alpha_history': np.array(alpha_history),
        'alpha_final': alpha,
        'predictions': pred_all,
        'bias_history': np.array(bias_history),
    }


# ══════════════════════════════════════════════════════════════════════════════
# RUN ALL ANALYSES
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'='*80}")
print("ANALYSIS 1: Prediction bias by α setting")
print(f"{'='*80}\n")

for mode_name, mode_val in [('Original Quartic (α=1.0)', 'quartic'),
                             ('Flat ISF (α=0.0)', 'flat'),
                             (f'Population mean (α={ALPHA_MEAN:.2f})', 'mean'),
                             ('Per-site optimal α', 'optimal')]:
    results = compute_bias(sites, alpha_mode=mode_val)
    print(f"\n  {mode_name}:")
    print(f"  {'Site':15s} {'α':>5s} {'<105r bias':>10s} {'>=105f bias':>11s} "
          f"{'MAE':>6s} {'Loop MAE':>9s}")
    print(f"  {'─'*60}")
    for r in results:
        print(f"  {r['name']:15s} {r['alpha']:5.2f} "
              f"{r.get('<105 rising_bias', np.nan):+10.1f} "
              f"{r.get('>=105 falling_bias', np.nan):+11.1f} "
              f"{r['overall_mae']:6.1f} {r['loop_mae']:9.1f}")

    # Aggregate
    total_bias_105r = []
    total_n_105r = 0
    total_bias_105f = []
    total_n_105f = 0
    total_mae = []
    loop_mae = []
    for r in results:
        n_r = r.get('<105 rising_n', 0)
        if n_r > 0:
            total_bias_105r.extend([r['<105 rising_bias']] * n_r)
            total_n_105r += n_r
        n_f = r.get('>=105 falling_n', 0)
        if n_f > 0:
            total_bias_105f.extend([r['>=105 falling_bias']] * n_f)
            total_n_105f += n_f
        total_mae.append(r['overall_mae'] * r['n'])
        loop_mae.append(r['loop_mae'] * r['n'])
    total_n = sum(r['n'] for r in results)


print(f"\n{'='*80}")
print("ANALYSIS 2: Adaptive α learning simulation")
print(f"{'='*80}\n")

adaptive_results = {}
for s in sites:
    result = simulate_adaptive(s, alpha_init=0.5, learning_rate=0.005, window=30)
    adaptive_results[s['name']] = result

    # Compute final bias
    final_pred = result['predictions']
    actual_end = s['actual_bg_end']
    bg = s['bg']
    pred_drop = s['pred_drop']

    mask_105r = (bg < 105) & (pred_drop <= 0)
    mask_105f = (bg >= 105) & (pred_drop > 0)

    bias_105r = np.mean(final_pred[mask_105r] - actual_end[mask_105r]) if mask_105r.sum() > 0 else np.nan
    bias_105f = np.mean(final_pred[mask_105f] - actual_end[mask_105f]) if mask_105f.sum() > 0 else np.nan
    overall_mae = np.mean(np.abs(final_pred - actual_end))

    opt_alpha = s.get('alpha_optimal', np.nan)
    print(f"  {s['name']:15s}  α: 0.50 → {result['alpha_final']:+.2f}  "
          f"(optimal={opt_alpha:+.2f} if known)  "
          f"<105r={bias_105r:+.1f}  >=105f={bias_105f:+.1f}  MAE={overall_mae:.1f}")


# ── Aggregate comparison table ──────────────────────────────────────────────

print(f"\n{'='*80}")
print("ANALYSIS 3: Aggregate comparison — all 13 subjects combined")
print(f"{'='*80}\n")

zones = ['<105 falling', '<105 rising', '>=105 falling', '>=105 rising']
models = {
    'Loop': {},
    'Quartic (α=1.0)': {},
    f'Damped (α={ALPHA_MEAN:.2f})': {},
    'Per-site optimal α': {},
    'Adaptive (learned α)': {},
}

for zone in zones:
    # Compute weighted mean bias for each model
    for model_name in models:
        total_bias_sum = 0
        total_n = 0
        for s in sites:
            bg = s['bg']
            pred_drop = s['pred_drop']
            actual_end = s['actual_bg_end']
            isf_actual = s['isf_actual']

            if zone == '<105 falling':
                mask = (bg < 105) & (pred_drop > 0)
            elif zone == '<105 rising':
                mask = (bg < 105) & (pred_drop <= 0)
            elif zone == '>=105 falling':
                mask = (bg >= 105) & (pred_drop > 0)
            else:
                mask = (bg >= 105) & (pred_drop <= 0)

            n = mask.sum()
            if n == 0:
                continue

            if model_name == 'Loop':
                pred = s['pred_loop'][mask]
            elif model_name == 'Adaptive (learned α)':
                pred = adaptive_results[s['name']]['predictions'][mask]
            else:
                anchor = quartic(99)
                S = (1800 / s['tdd']) / anchor
                if model_name == 'Quartic (α=1.0)':
                    alpha = 1.0
                elif 'Damped' in model_name:
                    alpha = ALPHA_MEAN
                elif model_name == 'Per-site optimal α':
                    alpha = s.get('alpha_optimal', ALPHA_MEAN)
                    if np.isnan(alpha):
                        alpha = ALPHA_MEAN

                isf_a = np.array([adaptive_isf(g, S, alpha) for g in bg[mask]])
                pred = bg[mask] - pred_drop[mask] * (isf_a / isf_actual[mask])

            bias = np.mean(pred - actual_end[mask])
            total_bias_sum += bias * n
            total_n += n

        models[model_name][zone] = total_bias_sum / total_n if total_n > 0 else np.nan

print(f"  {'Model':<25s}", end='')
for zone in zones:
    print(f"  {zone:>15s}", end='')
print()
print(f"  {'─'*90}")
for model_name, zone_biases in models.items():
    print(f"  {model_name:<25s}", end='')
    for zone in zones:
        print(f"  {zone_biases.get(zone, np.nan):+15.1f}", end='')
    print()

# Compute overall MAE for each model
print(f"\n  Overall MAE:")
for model_name in models:
    total_ae = 0
    total_n = 0
    for s in sites:
        actual_end = s['actual_bg_end']
        bg = s['bg']
        pred_drop = s['pred_drop']
        isf_actual = s['isf_actual']
        n = s['n']

        if model_name == 'Loop':
            pred = s['pred_loop']
        elif model_name == 'Adaptive (learned α)':
            pred = adaptive_results[s['name']]['predictions']
        else:
            anchor = quartic(99)
            S = (1800 / s['tdd']) / anchor
            if model_name == 'Quartic (α=1.0)':
                alpha = 1.0
            elif 'Damped' in model_name:
                alpha = ALPHA_MEAN
            elif model_name == 'Per-site optimal α':
                alpha = s.get('alpha_optimal', ALPHA_MEAN)
                if np.isnan(alpha):
                    alpha = ALPHA_MEAN

            isf_a = np.array([adaptive_isf(g, S, alpha) for g in bg])
            pred = bg - pred_drop * (isf_a / isf_actual)

        total_ae += np.sum(np.abs(pred - actual_end))
        total_n += n

    print(f"  {model_name:<25s}  MAE = {total_ae/total_n:.1f} mg/dL")


# ══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════════════════════

# ── Figure 1: Adaptive α convergence per site ───────────────────────────────

fig, axes = plt.subplots(3, 5, figsize=(20, 10))
axes = axes.flatten()
for i, s in enumerate(sites):
    if i >= 13:
        break
    ax = axes[i]
    ar = adaptive_results[s['name']]
    ax.plot(ar['alpha_history'], linewidth=0.8, color='steelblue')
    opt = s.get('alpha_optimal', np.nan)
    if not np.isnan(opt):
        ax.axhline(opt, color='red', linestyle='--', alpha=0.6, linewidth=1,
                   label=f'Optimal α={opt:.2f}')
    ax.axhline(0.5, color='grey', linestyle=':', alpha=0.3)
    ax.set_title(f"{s['name']}\n(n={s['n']}, TDD={s['tdd']:.0f})", fontsize=8)
    ax.set_ylim(-1, 2)
    ax.set_ylabel('α', fontsize=8)
    if i >= 8:
        ax.set_xlabel('Sample', fontsize=8)
    if not np.isnan(opt):
        ax.legend(fontsize=6)

for j in range(len(sites), len(axes)):
    axes[j].set_visible(False)

plt.suptitle('Adaptive α Convergence Per Site\n(blue = learned α over time, red dashed = optimal α)',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'adaptive_alpha_convergence.png', dpi=150, bbox_inches='tight')
print(f"\nSaved: {OUT_DIR / 'adaptive_alpha_convergence.png'}")


# ── Figure 2: Comparison bar chart — bias by zone and model ─────────────────

fig, axes = plt.subplots(1, 4, figsize=(18, 5))
model_names = list(models.keys())
colors = ['#4A90D9', '#D94A4A', '#F5A623', '#7ED321', '#9B59B6']

for i, zone in enumerate(zones):
    ax = axes[i]
    biases = [models[m].get(zone, 0) for m in model_names]
    bars = ax.bar(range(len(model_names)), biases, color=colors, alpha=0.8, edgecolor='black',
                  linewidth=0.5)
    ax.axhline(0, color='grey', linestyle='-', linewidth=0.5)
    ax.set_xticks(range(len(model_names)))
    ax.set_xticklabels([m.split('(')[0].strip() for m in model_names],
                       rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('Bias (mg/dL)')
    ax.set_title(zone, fontsize=10)
    # Add value labels
    for bar, val in zip(bars, biases):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f'{val:+.1f}', ha='center', va='bottom' if val >= 0 else 'top',
                fontsize=7)

plt.suptitle('Prediction Bias by Zone: All 13 Subjects Combined\n'
             '(predicted SGV at end-of-IOB minus actual SGV)', fontsize=11)
plt.tight_layout()
plt.savefig(OUT_DIR / 'adaptive_model_bias_comparison.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'adaptive_model_bias_comparison.png'}")


# ── Figure 3: Per-site α values with optimal and adaptive final ─────────────

fig, ax = plt.subplots(figsize=(14, 5))
x = np.arange(len(sites))
width = 0.35

opt_alphas = [s.get('alpha_optimal', np.nan) for s in sites]
adapt_alphas = [adaptive_results[s['name']]['alpha_final'] for s in sites]
site_names = [s['name'] for s in sites]

# Plot
bars1 = ax.bar(x - width/2, opt_alphas, width, label='Optimal α (from data)',
               color='#4A90D9', alpha=0.8, edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x + width/2, adapt_alphas, width, label='Adaptive α (learned)',
               color='#F5A623', alpha=0.8, edgecolor='black', linewidth=0.5)
ax.axhline(1.0, color='green', linestyle='--', alpha=0.5, label='Original quartic (α=1.0)')
ax.axhline(ALPHA_MEAN, color='red', linestyle='--', alpha=0.5, label=f'Population mean (α={ALPHA_MEAN:.2f})')
ax.axhline(0, color='grey', linestyle=':', alpha=0.3)

ax.set_xticks(x)
ax.set_xticklabels(site_names, rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Shape factor α')
ax.set_title('Per-Site Shape Factor: Optimal vs Adaptive Learning')
ax.legend(fontsize=8, loc='upper right')

plt.tight_layout()
plt.savefig(OUT_DIR / 'adaptive_alpha_per_site.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'adaptive_alpha_per_site.png'}")


# ── Figure 4: Model-corrected ISF curves overlaid on patient data ───────────

fig, axes = plt.subplots(3, 5, figsize=(20, 10))
axes = axes.flatten()
bg_smooth = np.linspace(74, 198, 200)

for i, s in enumerate(sites):
    if i >= 13:
        break
    ax = axes[i]

    # Patient actual ISF bars
    medians = bin_isf(s['bg'], s['isf_actual'])
    ax.bar(BG_BIN_CENTERS, medians, width=7, alpha=0.3, color='steelblue',
           edgecolor='navy', linewidth=0.5, label='Actual ISF')

    # TDD scaling
    anchor = quartic(99)
    S = (1800 / s['tdd']) / anchor

    # Original quartic
    y_quartic = [quartic(g) * S for g in bg_smooth]
    ax.plot(bg_smooth, y_quartic, 'r--', linewidth=1, alpha=0.7, label='Quartic (α=1.0)')

    # Damped quartic (population mean α)
    y_damped = [adaptive_isf(g, S, ALPHA_MEAN) for g in bg_smooth]
    ax.plot(bg_smooth, y_damped, 'g-', linewidth=1.5, alpha=0.8, label=f'Damped (α={ALPHA_MEAN:.2f})')

    # Optimal α
    opt = s.get('alpha_optimal', np.nan)
    if not np.isnan(opt):
        y_opt = [adaptive_isf(g, S, opt) for g in bg_smooth]
        ax.plot(bg_smooth, y_opt, color='orange', linewidth=1.5, alpha=0.8,
                label=f'Optimal (α={opt:.2f})')

    ax.set_title(f"{s['name']} (TDD={s['tdd']:.0f})", fontsize=8)
    ax.set_xlim(72, 200)
    ax.set_ylim(0, min(max(np.nanmax(medians) * 1.3, 50), 400))
    if i >= 8:
        ax.set_xlabel('BG (mg/dL)', fontsize=7)
    if i % 5 == 0:
        ax.set_ylabel('ISF (mg/dL/U)', fontsize=7)
    ax.legend(fontsize=5, loc='upper right')

for j in range(len(sites), len(axes)):
    axes[j].set_visible(False)

plt.suptitle('Per-Site ISF Curves: Actual vs Quartic vs Damped vs Optimal',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'adaptive_isf_curves_per_site.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'adaptive_isf_curves_per_site.png'}")


# ── Figure 5: Rolling prediction bias comparison ────────────────────────────

fig, axes = plt.subplots(3, 5, figsize=(20, 10))
axes = axes.flatten()

for i, s in enumerate(sites):
    if i >= 13:
        break
    ax = axes[i]
    bg = s['bg']
    actual_end = s['actual_bg_end']

    # Sort by BG for rolling window
    idx = np.argsort(bg)
    bg_sorted = bg[idx]

    # Loop bias
    loop_err = (s['pred_loop'] - actual_end)[idx]
    # Adaptive bias
    adapt_err = (adaptive_results[s['name']]['predictions'] - actual_end)[idx]

    # Rolling mean (window=50)
    w = min(50, len(bg_sorted) // 4)
    if w > 5:
        loop_roll = pd.Series(loop_err).rolling(w, center=True).mean().values
        adapt_roll = pd.Series(adapt_err).rolling(w, center=True).mean().values

        ax.plot(bg_sorted, loop_roll, 'b-', linewidth=1, alpha=0.7, label='Loop')
        ax.plot(bg_sorted, adapt_roll, color='orange', linewidth=1, alpha=0.7, label='Adaptive')
        ax.axhline(0, color='grey', linestyle=':', alpha=0.3)

    ax.set_title(f"{s['name']}", fontsize=8)
    ax.set_xlim(72, 200)
    ax.set_ylim(-60, 60)
    if i >= 8:
        ax.set_xlabel('Starting BG', fontsize=7)
    if i % 5 == 0:
        ax.set_ylabel('Bias (mg/dL)', fontsize=7)
    ax.legend(fontsize=6)

for j in range(len(sites), len(axes)):
    axes[j].set_visible(False)

plt.suptitle('Rolling Prediction Bias vs Starting Glucose\n(Loop vs Adaptive model)',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(OUT_DIR / 'adaptive_bias_vs_bg.png', dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_DIR / 'adaptive_bias_vs_bg.png'}")


# ── Summary statistics for the paper ────────────────────────────────────────

print(f"\n{'='*80}")
print("SUMMARY FOR PAPER")
print(f"{'='*80}")

print(f"\n  Total subjects: {len(sites)}")
print(f"  Total samples: {sum(s['n'] for s in sites)}")
print(f"  Evaluable for shape analysis: {len(evaluable)} (excluded: User-I, User-H, User-J)")
print(f"  α range: {min(all_alphas):.2f} to {max(all_alphas):.2f}")
print(f"  α mean: {ALPHA_MEAN:.2f}, median: {ALPHA_MEDIAN:.2f}")
print(f"  Quartic slope: {Q_SLOPE*1000:.1f} ×10⁻³")
print(f"  r(TDD, α) = {stats.pearsonr([s['tdd'] for s in evaluable], all_alphas)[0]:.3f}")

# Adaptive convergence summary
print(f"\n  Adaptive α convergence:")
for s in sites:
    ar = adaptive_results[s['name']]
    opt = s.get('alpha_optimal', np.nan)
    final = ar['alpha_final']
    diff = abs(final - opt) if not np.isnan(opt) else np.nan
    print(f"    {s['name']:15s}  final={final:+.2f}  optimal={opt:+.2f}  |diff|={diff:.2f}"
          if not np.isnan(opt) else
          f"    {s['name']:15s}  final={final:+.2f}  optimal=N/A")

print(f"\n  DONE")
