#!/usr/bin/env python3
"""
Investigate WHY ISF curve shapes differ between patients
and whether a shape correction factor could fix the generic model.

Questions:
1. Does the normalised ISF slope correlate with TDD, model type, or other measurables?
2. If we add a per-patient shape factor α to the quartic, what α does each site need?
3. Can α be derived from anything measurable, or does it require per-patient calibration?
4. What would prediction bias look like with a shape-corrected quartic?
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

OUT_DIR = Path.home() / 'Downloads' / '4 Hour analysis'
CACHE_FILE = OUT_DIR / 'multisite_4h_sample_cache.pkl'

BG_BINS = list(range(72, 201, 8))
BG_BIN_CENTERS = [(BG_BINS[i] + BG_BINS[i + 1]) / 2 for i in range(len(BG_BINS) - 1)]
REF_BG = 100

def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

def full_diabeloop(g):
    return (98.03 - 1.077 * g + 0.008868 * g**2) if g <= 100 else quartic(g)

# ── Load cache ──────────────────────────────────────────────────────────────

with open(CACHE_FILE, 'rb') as f:
    all_sites = pickle.load(f)

# Anonymise site names
ANON_MAP = {
    'henny425': 'User-A', 'aadiabetes': 'User-B', 'diajesse': 'User-C',
    'svns': 'User-D', 'fuxchr': 'User-E', 'mikens': 'User-F',
    'andycgm': 'User-G', 'noahr': 'User-H', 'nightscout1': 'User-I',
    'eli': 'User-J', 'ns_rot6': 'User-K', 'kelseyhuss': 'User-L',
}
for s in all_sites:
    s['name'] = ANON_MAP.get(s['name'], s['name'])

print(f"Loaded {len(all_sites)} sites\n")

# ── Compute per-site normalised ISF and slope ───────────────────────────────

ref_idx = None
for i in range(len(BG_BINS) - 1):
    if BG_BINS[i] <= REF_BG < BG_BINS[i + 1]:
        ref_idx = i
        break

site_data = []
for s in all_sites:
    bg = s['bg']
    isf = s['isf_actual']
    # Bin ISF
    medians = []
    for i in range(len(BG_BINS) - 1):
        mask = (bg >= BG_BINS[i]) & (bg < BG_BINS[i + 1])
        n = mask.sum()
        medians.append(np.median(isf[mask]) if n >= 5 else np.nan)
    medians = np.array(medians)

    ref_isf = medians[ref_idx]
    if np.isnan(ref_isf) or ref_isf <= 0:
        continue
    if s['name'] == 'User-I':
        continue  # constant ISF

    normalised = medians / ref_isf
    valid = ~np.isnan(normalised)
    if valid.sum() < 4:
        continue

    x = np.array(BG_BIN_CENTERS)[valid]
    y = normalised[valid]
    slope, intercept, r, p, se = stats.linregress(x, y)

    site_data.append({
        'name': s['name'],
        'model': s['model'],
        'tdd': s['tdd_median'],
        'n': s['n'],
        'slope': slope,
        'r2': r**2,
        'ref_isf': ref_isf,
        'normalised': normalised,
    })

# Quartic normalised slope
q_ref = quartic(REF_BG)
q_x = np.array(BG_BIN_CENTERS)
q_y = np.array([quartic(g) / q_ref for g in q_x])
q_slope = stats.linregress(q_x, q_y).slope

print(f"Quartic normalised slope: {q_slope*1000:.1f} ×10⁻³/mg/dL\n")

# ── Analysis 1: Does slope correlate with TDD? ─────────────────────────────

print("=" * 70)
print("ANALYSIS 1: Does ISF curve shape correlate with TDD?")
print("=" * 70)

names = [d['name'] for d in site_data]
slopes = np.array([d['slope'] for d in site_data])
tdds = np.array([d['tdd'] for d in site_data])
models = [d['model'] for d in site_data]

r_tdd, p_tdd = stats.pearsonr(tdds, slopes)
print(f"\n  Pearson r(TDD, slope) = {r_tdd:.3f}, p = {p_tdd:.3f}")
print(f"  {'Significant' if p_tdd < 0.05 else 'Not significant'} at p < 0.05")

for d in site_data:
    print(f"    {d['name']:15s}  TDD={d['tdd']:5.1f}  slope={d['slope']*1000:+6.1f}×10⁻³  "
          f"model={d['model']}")

# ── Analysis 2: Does slope differ by model type? ───────────────────────────

print(f"\n{'=' * 70}")
print("ANALYSIS 2: Slope by model type (sigmoid vs logarithmic)")
print("=" * 70)

sig_slopes = [d['slope'] for d in site_data if d['model'] == 'sigmoid']
log_slopes = [d['slope'] for d in site_data if d['model'] == 'log']

print(f"\n  Sigmoid (n={len(sig_slopes)}): mean={np.mean(sig_slopes)*1000:.1f}, "
      f"range=[{min(sig_slopes)*1000:.1f}, {max(sig_slopes)*1000:.1f}] ×10⁻³")
print(f"  Log     (n={len(log_slopes)}): mean={np.mean(log_slopes)*1000:.1f}, "
      f"range=[{min(log_slopes)*1000:.1f}, {max(log_slopes)*1000:.1f}] ×10⁻³")

if len(sig_slopes) >= 2 and len(log_slopes) >= 2:
    t, p_model = stats.ttest_ind(sig_slopes, log_slopes)
    print(f"  t-test: t={t:.2f}, p={p_model:.3f}")
    print(f"  {'Significant' if p_model < 0.05 else 'Not significant'} difference between model types")

# ── Analysis 3: Shape correction factor α ───────────────────────────────────

print(f"\n{'=' * 70}")
print("ANALYSIS 3: Per-site shape correction factor")
print("=" * 70)
print(f"\n  If quartic_corrected(G) = S × [1 + α × (quartic_norm(G) - 1)]")
print(f"  where α = patient_slope / quartic_slope scales the curve's deviation from flat")
print(f"  α = 1.0 → original quartic shape")
print(f"  α = 0.0 → flat ISF (no glucose dependence)")
print(f"  α < 0   → ISF increases with glucose (inverted)")

alphas = []
for d in site_data:
    alpha = d['slope'] / q_slope
    alphas.append(alpha)
    print(f"    {d['name']:15s}  slope={d['slope']*1000:+6.1f}  α={alpha:+5.2f}")

print(f"\n  α range: {min(alphas):.2f} to {max(alphas):.2f}")
print(f"  α mean:  {np.mean(alphas):.2f}")
print(f"  α std:   {np.std(alphas):.2f}")

# Does α correlate with TDD?
r_alpha_tdd, p_alpha_tdd = stats.pearsonr(tdds, alphas)
print(f"\n  Pearson r(TDD, α) = {r_alpha_tdd:.3f}, p = {p_alpha_tdd:.3f}")
print(f"  {'Could derive α from TDD' if p_alpha_tdd < 0.05 else 'Cannot derive α from TDD — requires per-patient calibration'}")

# ── Analysis 4: Simulate shape-corrected quartic prediction bias ────────────

print(f"\n{'=' * 70}")
print("ANALYSIS 4: Simulated prediction bias with shape-corrected quartic")
print("=" * 70)

# For each site, compute what the prediction bias would be if we used
# a quartic with that site's own α (best case: we know the right α)
print(f"\n  Using per-site optimal α (best possible with this approach):")

for i, s in enumerate(all_sites):
    # Find matching site_data entry
    sd = next((d for d in site_data if d['name'] == s['name']), None)
    if sd is None:
        continue

    alpha = sd['slope'] / q_slope
    bg = s['bg']
    isf_actual = s['isf_actual']
    pred_drop = s['pred_drop']
    actual_end = s['actual_bg_end']

    # TDD scaling factor
    anchor = quartic(99)
    S = (1800 / sd['tdd']) / anchor

    # Original quartic ISF
    isf_quartic = np.array([quartic(g) * S for g in bg])

    # Shape-corrected quartic ISF
    q_norm_vals = np.array([quartic(g) / q_ref for g in bg])
    isf_corrected = np.array([S * q_ref * (1 + alpha * (qn - 1)) for qn in q_norm_vals])

    # Counterfactual predictions
    current_sgv = bg
    pred_q = current_sgv - pred_drop * (isf_quartic / isf_actual)
    pred_c = current_sgv - pred_drop * (isf_corrected / isf_actual)

    # Bias in <105 rising
    mask_105r = (bg < 105) & (pred_drop < 0)
    if mask_105r.sum() > 10:
        bias_q = np.mean(pred_q[mask_105r] - actual_end[mask_105r])
        bias_c = np.mean(pred_c[mask_105r] - actual_end[mask_105r])
        bias_loop = np.mean(current_sgv[mask_105r] - pred_drop[mask_105r] - actual_end[mask_105r])
        print(f"    {sd['name']:15s}  α={alpha:+5.2f}  <105r: Loop={bias_loop:+6.1f}  "
              f"Quartic={bias_q:+6.1f}  Corrected={bias_c:+6.1f}")

# ── Analysis 5: What if we use the MEAN α? ─────────────────────────────────

print(f"\n  Using population mean α = {np.mean(alphas):.2f} (all sites get same α):")

alpha_mean = np.mean(alphas)
for i, s in enumerate(all_sites):
    sd = next((d for d in site_data if d['name'] == s['name']), None)
    if sd is None:
        continue

    bg = s['bg']
    isf_actual = s['isf_actual']
    pred_drop = s['pred_drop']
    actual_end = s['actual_bg_end']

    anchor = quartic(99)
    S = (1800 / sd['tdd']) / anchor

    q_norm_vals = np.array([quartic(g) / q_ref for g in bg])
    isf_corrected = np.array([S * q_ref * (1 + alpha_mean * (qn - 1)) for qn in q_norm_vals])

    current_sgv = bg
    pred_c = current_sgv - pred_drop * (isf_corrected / isf_actual)

    mask_105r = (bg < 105) & (pred_drop < 0)
    if mask_105r.sum() > 10:
        bias_c = np.mean(pred_c[mask_105r] - actual_end[mask_105r])
        print(f"    {sd['name']:15s}  <105r Corrected={bias_c:+6.1f}")

# ── PLOTS ───────────────────────────────────────────────────────────────────

# ── Analysis 6: Real-world differences between sigmoid and log users ────────

print(f"\n{'=' * 70}")
print("ANALYSIS 6: Real-world differences between sigmoid and log model users")
print("=" * 70)
print(f"\n  Do users on different model types have different actual glucose outcomes?")
print(f"  (What the model choice means in practice)\n")

print(f"  {'Site':15s} {'Model':7s} {'TDD':>5s} {'N':>5s} {'Mean BG':>8s} {'Med BG':>7s} "
      f"{'Mean End':>9s} {'<80 %':>6s} {'>140 %':>7s} {'ISF@100':>8s} {'Slope':>8s}")
print("  " + "─" * 100)

for s in all_sites:
    sd = next((d for d in site_data if d['name'] == s['name']), None)
    bg = s['bg']
    actual_end = s['actual_bg_end']
    pct_below80 = 100 * np.mean(bg < 80)
    pct_above140 = 100 * np.mean(bg > 140)

    slope_str = f"{sd['slope']*1000:+6.1f}" if sd else "  N/A"
    isf_ref = f"{sd['ref_isf']:6.1f}" if sd else "  N/A"

    print(f"  {s['name']:15s} {s['model']:7s} {s['tdd_median']:5.1f} {s['n']:5d} "
          f"{np.mean(bg):8.1f} {np.median(bg):7.1f} {np.mean(actual_end):9.1f} "
          f"{pct_below80:6.1f} {pct_above140:7.1f} {isf_ref:>8s} {slope_str:>8s}")

# Compare actual outcomes
sig_sites = [s for s in all_sites if s['model'] == 'sigmoid']
log_sites = [s for s in all_sites if s['model'] == 'log']

sig_all_bg = np.concatenate([s['bg'] for s in sig_sites])
log_all_bg = np.concatenate([s['bg'] for s in log_sites])
sig_all_end = np.concatenate([s['actual_bg_end'] for s in sig_sites])
log_all_end = np.concatenate([s['actual_bg_end'] for s in log_sites])

print(f"\n  Sigmoid aggregate: n={len(sig_all_bg)}, mean start BG={np.mean(sig_all_bg):.1f}, "
      f"mean end BG={np.mean(sig_all_end):.1f}")
print(f"  Log aggregate:     n={len(log_all_bg)}, mean start BG={np.mean(log_all_bg):.1f}, "
      f"mean end BG={np.mean(log_all_end):.1f}")

# ISF at 100 by model type
sig_isf = [d['ref_isf'] for d in site_data if d['model'] == 'sigmoid']
log_isf = [d['ref_isf'] for d in site_data if d['model'] == 'log']
print(f"\n  ISF at 100 mg/dL — Sigmoid: {np.mean(sig_isf):.1f} ± {np.std(sig_isf):.1f} "
      f"(range {min(sig_isf):.0f}-{max(sig_isf):.0f})")
print(f"  ISF at 100 mg/dL — Log:     {np.mean(log_isf):.1f} ± {np.std(log_isf):.1f} "
      f"(range {min(log_isf):.0f}-{max(log_isf):.0f})")

# TDD by model type
sig_tdd = [d['tdd'] for d in site_data if d['model'] == 'sigmoid']
log_tdd = [d['tdd'] for d in site_data if d['model'] == 'log']
print(f"\n  TDD — Sigmoid: {np.mean(sig_tdd):.1f} ± {np.std(sig_tdd):.1f} "
      f"(range {min(sig_tdd):.0f}-{max(sig_tdd):.0f})")
print(f"  TDD — Log:     {np.mean(log_tdd):.1f} ± {np.std(log_tdd):.1f} "
      f"(range {min(log_tdd):.0f}-{max(log_tdd):.0f})")

# Prediction accuracy by model type (Loop's own predictions)
sig_pred = np.concatenate([s['pred_loop'] for s in sig_sites])
log_pred = np.concatenate([s['pred_loop'] for s in log_sites])
sig_bias = np.mean(sig_pred - sig_all_end)
log_bias = np.mean(log_pred - log_all_end)
sig_mae = np.mean(np.abs(sig_pred - sig_all_end))
log_mae = np.mean(np.abs(log_pred - log_all_end))

print(f"\n  Loop prediction accuracy (its own model):")
print(f"  Sigmoid: bias={sig_bias:+.1f} mg/dL, MAE={sig_mae:.1f} mg/dL")
print(f"  Log:     bias={log_bias:+.1f} mg/dL, MAE={log_mae:.1f} mg/dL")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Slope vs TDD
ax = axes[0, 0]
for d in site_data:
    color = '#4A90D9' if d['model'] == 'sigmoid' else '#D94A4A'
    marker = 'o' if d['model'] == 'sigmoid' else 's'
    ax.scatter(d['tdd'], d['slope'] * 1000, color=color, marker=marker,
              s=80, zorder=5, edgecolors='black', linewidth=0.5)
    ax.annotate(d['name'], (d['tdd'], d['slope'] * 1000),
               fontsize=7, xytext=(5, 5), textcoords='offset points')
ax.axhline(q_slope * 1000, linestyle='--', color='green', alpha=0.6, label=f'Quartic ({q_slope*1000:.1f})')
ax.axhline(0, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('TDD (U/day)')
ax.set_ylabel('Normalised ISF slope (×10⁻³/mg/dL)')
ax.set_title(f'Curve Slope vs TDD\nr={r_tdd:.2f}, p={p_tdd:.3f}')
ax.legend(fontsize=8)

# Plot 2: α values per site
ax = axes[0, 1]
colors = ['#4A90D9' if d['model'] == 'sigmoid' else '#D94A4A' for d in site_data]
bars = ax.bar(range(len(site_data)), alphas, color=colors, alpha=0.8, edgecolor='navy')
ax.axhline(1.0, linestyle='--', color='green', alpha=0.6, label='α=1.0 (original quartic)')
ax.axhline(0.0, linestyle=':', color='grey', alpha=0.3, label='α=0.0 (flat)')
ax.axhline(np.mean(alphas), linestyle='--', color='orange', alpha=0.6,
           label=f'Mean α={np.mean(alphas):.2f}')
ax.set_xticks(range(len(site_data)))
ax.set_xticklabels([d['name'] for d in site_data], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Shape factor α')
ax.set_title('Per-Site Shape Correction Factor\n(α = patient slope / quartic slope)')
ax.legend(fontsize=8)

# Plot 3: Corrected vs original normalised curves
ax = axes[1, 0]
bg_smooth = np.linspace(74, 198, 200)
for d in site_data:
    alpha = d['slope'] / q_slope
    q_norm_smooth = np.array([quartic(g) / q_ref for g in bg_smooth])
    corrected = 1 + alpha * (q_norm_smooth - 1)
    style = '-' if d['model'] == 'sigmoid' else '--'
    ax.plot(bg_smooth, corrected, style, linewidth=1, alpha=0.5, label=d['name'])
# Original quartic for reference
q_norm_smooth = np.array([quartic(g) / q_ref for g in bg_smooth])
ax.plot(bg_smooth, q_norm_smooth, 'k-', linewidth=2.5, label='Original quartic', zorder=10)
ax.axhline(1.0, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('Normalised ISF')
ax.set_title('Shape-Corrected Quartic per Site\n(each site gets its optimal α)')
ax.legend(fontsize=6, loc='upper right', ncol=2)
ax.set_ylim(0, 2.5)
ax.set_xlim(72, 200)

# Plot 4: α vs TDD
ax = axes[1, 1]
for d, alpha in zip(site_data, alphas):
    color = '#4A90D9' if d['model'] == 'sigmoid' else '#D94A4A'
    marker = 'o' if d['model'] == 'sigmoid' else 's'
    ax.scatter(d['tdd'], alpha, color=color, marker=marker,
              s=80, zorder=5, edgecolors='black', linewidth=0.5)
    ax.annotate(d['name'], (d['tdd'], alpha),
               fontsize=7, xytext=(5, 5), textcoords='offset points')
ax.axhline(1.0, linestyle='--', color='green', alpha=0.6, label='Quartic shape')
ax.axhline(0.0, linestyle=':', color='grey', alpha=0.3)
ax.set_xlabel('TDD (U/day)')
ax.set_ylabel('Shape factor α')
ax.set_title(f'Shape Factor vs TDD\nr={r_alpha_tdd:.2f}, p={p_alpha_tdd:.3f}')
ax.legend(fontsize=8)

plt.tight_layout()
out = OUT_DIR / 'shape_correction_analysis.png'
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out}")

print(f"\n{'=' * 70}")
print("SUMMARY")
print("=" * 70)
print(f"\n  Patient normalised slopes: {min(slopes)*1000:+.1f} to {max(slopes)*1000:+.1f} ×10⁻³/mg/dL")
print(f"  Quartic normalised slope:  {q_slope*1000:.1f} ×10⁻³/mg/dL")
print(f"  Shape factor α range:      {min(alphas):.2f} to {max(alphas):.2f}")
print(f"  Correlation r(TDD, α):     {r_alpha_tdd:.3f} (p={p_alpha_tdd:.3f})")
print(f"  Correlation r(TDD, slope): {r_tdd:.3f} (p={p_tdd:.3f})")
if p_alpha_tdd >= 0.05:
    print(f"\n  CONCLUSION: Shape factor α does NOT correlate with TDD.")
    print(f"  A generic model needs per-patient shape calibration, not just TDD scaling.")
    print(f"  This is functionally what Trio's adjustment factor and model selection provide.")
else:
    print(f"\n  CONCLUSION: Shape factor α correlates with TDD — a two-parameter")
    print(f"  model using TDD for both magnitude and shape might work.")
print()
