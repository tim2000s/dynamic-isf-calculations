#!/usr/bin/env python3
"""
Analyse how ISF-glucose curve shapes differ between patients.

Loads the cached per-site data, normalises each site's ISF curve by its
value at a reference glucose (~100 mg/dL) to remove magnitude differences,
then quantifies and visualises the shape variation across patients.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path.home() / 'Downloads' / '4 Hour analysis'
CACHE_FILE = OUT_DIR / 'multisite_4h_sample_cache.pkl'

# Glucose bins (must match the main analysis)
BG_BINS = list(range(72, 201, 8))
BG_BIN_CENTERS = [(BG_BINS[i] + BG_BINS[i + 1]) / 2 for i in range(len(BG_BINS) - 1)]
REF_BG = 100  # Reference glucose for normalisation

# Diabeloop formulas (unscaled, for shape comparison)
def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

def full_diabeloop(g):
    if g <= 100:
        return 98.03 - 1.077 * g + 0.008868 * g**2
    return quartic(g)

def hybrid(g):
    if g < 105:
        return 75.8 * (105 / g) ** 3.5
    return quartic(g)


# ── Load cache ───────────────────────────────────────────────────────────────

with open(CACHE_FILE, 'rb') as f:
    all_sites = pickle.load(f)

print(f"Loaded {len(all_sites)} sites from cache\n")


# ── Compute per-site median ISF by glucose bin ───────────────────────────────

def bin_isf(bg_arr, isf_arr):
    """Return median ISF per glucose bin and the count per bin."""
    medians = []
    counts = []
    for i in range(len(BG_BINS) - 1):
        lo, hi = BG_BINS[i], BG_BINS[i + 1]
        mask = (bg_arr >= lo) & (bg_arr < hi)
        n = mask.sum()
        counts.append(n)
        if n >= 5:
            medians.append(np.median(isf_arr[mask]))
        else:
            medians.append(np.nan)
    return np.array(medians), np.array(counts)


site_curves = []
for s in all_sites:
    medians, counts = bin_isf(s['bg'], s['isf_actual'])
    site_curves.append({
        'name': s['name'],
        'model': s['model'],
        'tdd': s['tdd_median'],
        'n': s['n'],
        'medians': medians,
        'counts': counts,
    })

# ── Normalise by ISF at reference glucose ────────────────────────────────────

# Find the bin containing the reference glucose
ref_idx = None
for i in range(len(BG_BINS) - 1):
    if BG_BINS[i] <= REF_BG < BG_BINS[i + 1]:
        ref_idx = i
        break

print(f"Reference bin: {BG_BINS[ref_idx]}-{BG_BINS[ref_idx+1]} mg/dL "
      f"(centre {BG_BIN_CENTERS[ref_idx]})\n")

for sc in site_curves:
    ref_isf = sc['medians'][ref_idx]
    if np.isnan(ref_isf) or ref_isf <= 0:
        sc['normalised'] = np.full_like(sc['medians'], np.nan)
        print(f"  {sc['name']}: no valid ISF at reference — skipping normalisation")
    else:
        sc['normalised'] = sc['medians'] / ref_isf
        print(f"  {sc['name']:15s}  model={sc['model']:7s}  TDD={sc['tdd']:5.1f}  "
              f"ISF@ref={ref_isf:.1f}  range={np.nanmin(sc['normalised']):.2f}-"
              f"{np.nanmax(sc['normalised']):.2f}")

# Normalise the Diabeloop formulas the same way
bg_smooth = np.linspace(74, 198, 200)
q_vals = np.array([quartic(g) for g in bg_smooth])
db_vals = np.array([full_diabeloop(g) for g in bg_smooth])
h_vals = np.array([hybrid(g) for g in bg_smooth])

q_ref = quartic(REF_BG)
db_ref = full_diabeloop(REF_BG)
h_ref = hybrid(REF_BG)

q_norm = q_vals / q_ref
db_norm = db_vals / db_ref
h_norm = h_vals / h_ref


# ── Quantify shape variation ─────────────────────────────────────────────────

# Stack all normalised curves (exclude nightscout1 — constant ISF, not dynamic)
all_norm = np.array([sc['normalised'] for sc in site_curves
                     if sc['name'] != 'nightscout1'])

# Per-bin statistics across sites
bin_mean = np.nanmean(all_norm, axis=0)
bin_std = np.nanstd(all_norm, axis=0)
bin_cv = bin_std / bin_mean  # Coefficient of variation
bin_min = np.nanmin(all_norm, axis=0)
bin_max = np.nanmax(all_norm, axis=0)
bin_range = bin_max - bin_min

print(f"\n{'Glucose':>8s} {'Mean':>6s} {'SD':>6s} {'CV':>6s} {'Min':>6s} {'Max':>6s} {'Range':>6s}")
print("─" * 50)
for i, gc in enumerate(BG_BIN_CENTERS):
    valid = np.sum(~np.isnan(all_norm[:, i]))
    if valid >= 3:
        print(f"{gc:8.0f} {bin_mean[i]:6.2f} {bin_std[i]:6.2f} {bin_cv[i]:6.2f} "
              f"{bin_min[i]:6.2f} {bin_max[i]:6.2f} {bin_range[i]:6.2f}")

# Split by model type
# Exclude nightscout1 from log group — its ISF is constant (13 mg/dL/U in every bin),
# meaning dynamic ISF is not active or is clamped at autosens_min.
sig_norm = np.array([sc['normalised'] for sc in site_curves if sc['model'] == 'sigmoid'])
log_norm = np.array([sc['normalised'] for sc in site_curves
                     if sc['model'] == 'log' and sc['name'] != 'nightscout1'])

# Flag it
for sc in site_curves:
    if sc['name'] == 'nightscout1':
        print(f"\n  NOTE: nightscout1 excluded from log mean — ISF is constant "
              f"({sc['medians'][0]:.0f} mg/dL/U in every bin, dynamic ISF not active or clamped)")


# ── PLOT 1: All normalised curves overlaid ───────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel 1: All sites
ax = axes[0]
for sc in site_curves:
    style = '-' if sc['model'] == 'sigmoid' else '--'
    ax.plot(BG_BIN_CENTERS, sc['normalised'], style, linewidth=1.5, alpha=0.7,
            label=f"{sc['name']} ({sc['model'][:3]})")
ax.axhline(1.0, color='grey', linestyle=':', alpha=0.5)
ax.axvline(105, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('Normalised ISF (ISF / ISF@100)')
ax.set_title('All Sites — Normalised ISF Curves')
ax.legend(fontsize=7, loc='upper right')
ax.set_ylim(0, 3.5)
ax.set_xlim(72, 200)

# Panel 2: Sigmoid vs Log mean ± SD
ax = axes[1]
sig_mean = np.nanmean(sig_norm, axis=0)
sig_std = np.nanstd(sig_norm, axis=0)
log_mean = np.nanmean(log_norm, axis=0)
log_std = np.nanstd(log_norm, axis=0)

ax.plot(BG_BIN_CENTERS, sig_mean, 'b-', linewidth=2, label='Sigmoid mean')
ax.fill_between(BG_BIN_CENTERS, sig_mean - sig_std, sig_mean + sig_std,
                alpha=0.2, color='blue')
ax.plot(BG_BIN_CENTERS, log_mean, 'r-', linewidth=2, label='Log mean')
ax.fill_between(BG_BIN_CENTERS, log_mean - log_std, log_mean + log_std,
                alpha=0.2, color='red')
ax.axhline(1.0, color='grey', linestyle=':', alpha=0.5)
ax.axvline(105, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('Normalised ISF (ISF / ISF@100)')
ax.set_title('Sigmoid vs Logarithmic — Mean ± SD')
ax.legend(fontsize=9)
ax.set_ylim(0, 3.5)
ax.set_xlim(72, 200)

# Panel 3: Diabeloop formula shapes vs patient envelope
ax = axes[2]
ax.fill_between(BG_BIN_CENTERS, bin_min, bin_max, alpha=0.15, color='grey',
                label='Patient range')
ax.plot(BG_BIN_CENTERS, bin_mean, 'k-', linewidth=2, label='Patient mean')
ax.plot(bg_smooth, q_norm, 'g--', linewidth=1.5, label='Quartic')
ax.plot(bg_smooth, db_norm, 'r--', linewidth=1.5, label='Full Diabeloop')
ax.plot(bg_smooth, h_norm, color='orange', linestyle='--', linewidth=1.5, label='Hybrid')
ax.axhline(1.0, color='grey', linestyle=':', alpha=0.5)
ax.axvline(105, color='grey', linestyle=':', alpha=0.3)
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('Normalised ISF (ISF / ISF@100)')
ax.set_title('Diabeloop Shapes vs Patient Envelope')
ax.legend(fontsize=9)
ax.set_ylim(0, 3.5)
ax.set_xlim(72, 200)

plt.tight_layout()
out1 = OUT_DIR / 'curve_shape_comparison.png'
plt.savefig(out1, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out1}")


# ── PLOT 2: Shape variation metrics ──────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: CV across glucose range
ax = axes[0]
ax.bar(BG_BIN_CENTERS, bin_cv * 100, width=7, color='steelblue', alpha=0.7, edgecolor='navy')
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('Coefficient of Variation (%)')
ax.set_title('Shape Variation Across Glucose Range\n(CV of normalised ISF, dynamic sites only)')
ax.axvline(105, color='grey', linestyle=':', alpha=0.3)
ax.set_xlim(72, 200)

# Panel 2: Min-max spread
ax = axes[1]
ax.fill_between(BG_BIN_CENTERS, bin_min, bin_max, alpha=0.3, color='steelblue',
                label='Min-Max range')
ax.plot(BG_BIN_CENTERS, bin_mean, 'k-', linewidth=2, label='Mean')
ax.plot(BG_BIN_CENTERS, bin_min, 'b--', linewidth=1, alpha=0.5)
ax.plot(BG_BIN_CENTERS, bin_max, 'b--', linewidth=1, alpha=0.5)
ax.set_xlabel('Sensor Glucose (mg/dL)')
ax.set_ylabel('Normalised ISF (ISF / ISF@100)')
ax.set_title('Normalised ISF Spread Across Patients')
ax.axhline(1.0, color='grey', linestyle=':', alpha=0.5)
ax.axvline(105, color='grey', linestyle=':', alpha=0.3)
ax.legend(fontsize=9)
ax.set_ylim(0, 3.5)
ax.set_xlim(72, 200)

plt.tight_layout()
out2 = OUT_DIR / 'curve_shape_variation.png'
plt.savefig(out2, dpi=150, bbox_inches='tight')
print(f"Saved: {out2}")


# ── PLOT 3: Individual site normalised curves with slope metric ──────────────

# Compute shape metric using the slope of normalised ISF across available bins.
# Use linear regression of normalised ISF vs glucose across all bins with data.
# Slope = how much normalised ISF changes per mg/dL of glucose.
# Negative slope = ISF decreases as glucose rises (steeper dynamic response).
# Near-zero slope = flat ISF curve (less glucose-dependent).

print(f"\nShape metric: Linear slope of normalised ISF vs glucose")
print(f"{'Site':>15s} {'Model':>7s} {'Slope':>10s} {'R²':>6s} {'Bins':>5s} {'Low':>6s} {'High':>6s}")
print("─" * 62)

shape_data = []
for sc in site_curves:
    if sc['name'] == 'nightscout1':
        print(f"  {'nightscout1':>15s} {'log':>7s}    EXCLUDED — constant ISF, dynamic ISF not active")
        continue
    valid = ~np.isnan(sc['normalised'])
    n_valid = valid.sum()
    if n_valid >= 4:
        x = np.array(BG_BIN_CENTERS)[valid]
        y = sc['normalised'][valid]
        # Linear fit
        coeffs = np.polyfit(x, y, 1)
        slope = coeffs[0]
        # R²
        y_pred = np.polyval(coeffs, x)
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        shape_data.append({
            'name': sc['name'], 'model': sc['model'],
            'slope': slope, 'r2': r2, 'n_bins': n_valid,
            'low': np.min(y), 'high': np.max(y)
        })
        print(f"{sc['name']:>15s} {sc['model']:>7s} {slope:10.5f} {r2:6.3f} "
              f"{n_valid:5d} {np.min(y):6.2f} {np.max(y):6.2f}")

# Same for Diabeloop formulas (evaluated at bin centres)
for name, func, ref in [('Quartic', quartic, q_ref), ('Full DB', full_diabeloop, db_ref),
                         ('Hybrid', hybrid, h_ref)]:
    x = np.array(BG_BIN_CENTERS)
    y = np.array([func(g) / ref for g in x])
    coeffs = np.polyfit(x, y, 1)
    slope = coeffs[0]
    y_pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    print(f"{name:>15s} {'formula':>7s} {slope:10.5f} {r2:6.3f} "
          f"{len(x):5d} {np.min(y):6.2f} {np.max(y):6.2f}")

# Bar chart of the slope metric
fig, ax = plt.subplots(figsize=(12, 5))
names = [d['name'] for d in shape_data]
slopes = [d['slope'] * 1000 for d in shape_data]  # ×1000 for readability
colors = ['#4A90D9' if d['model'] == 'sigmoid' else '#D94A4A' for d in shape_data]

bars = ax.bar(range(len(names)), slopes, color=colors, alpha=0.8, edgecolor='navy')

# Add Diabeloop formula slopes
for name, func, ref in [('Quartic', quartic, q_ref),
                         ('Full DB', full_diabeloop, db_ref),
                         ('Hybrid', hybrid, h_ref)]:
    x = np.array(BG_BIN_CENTERS)
    y = np.array([func(g) / ref for g in x])
    s = np.polyfit(x, y, 1)[0] * 1000
    ax.axhline(s, linestyle='--', alpha=0.6,
               color='green' if 'Quartic' in name else ('red' if 'Full' in name else 'orange'),
               label=f'{name} ({s:.1f})')

ax.set_xticks(range(len(names)))
ax.set_xticklabels(names, rotation=45, ha='right')
ax.set_ylabel('Normalised ISF slope (×1000 per mg/dL)')
ax.set_title('Curve Steepness: Linear Slope of Normalised ISF vs Glucose\n'
             'Blue = sigmoid sites, Red = logarithmic sites, Dashed = Diabeloop formulas\n'
             'More negative = ISF decreases faster as glucose rises')
ax.axhline(0, color='grey', linestyle=':', alpha=0.3)
ax.legend(fontsize=9, loc='lower left')

plt.tight_layout()
out3 = OUT_DIR / 'curve_steepness_comparison.png'
plt.savefig(out3, dpi=150, bbox_inches='tight')
print(f"\nSaved: {out3}")

print("\nDONE")
