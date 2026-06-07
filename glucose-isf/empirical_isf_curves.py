#!/usr/bin/env python3
"""
Empirical ISF-Glucose Curves from Observed Outcomes
====================================================
Computes the ISF value that WOULD have produced the correct prediction
for each fasting sample, then bins by starting glucose to reveal each
patient's true ISF-glucose relationship.

ISF_observed = isf_actual * (actual_drop / pred_drop)

Charts:
  1. empirical_isf_raw.png       — 4x4 small-multiples, scatter + binned medians
  2. empirical_isf_normalised.png — all 13 normalised ratio curves vs models
  3. empirical_isf_model_fit.png  — RMSE bar chart + population mean curve
  4. empirical_isf_overnight_vs_daytime.png — period comparison
"""

import pickle
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path(__file__).parent
TRIO_CACHE = OUT_DIR / 'daytime analysis' / 'multisite_allday_cache.pkl'
BOOST_CACHE = OUT_DIR / 'daytime analysis' / 'boost_allday_cache.pkl'

ANON_MAP = {
    'henny425': 'User-A', 'aadiabetes': 'User-B', 'diajesse': 'User-C',
    'svns': 'User-D', 'fuxchr': 'User-E', 'mikens': 'User-F',
    'andycgm': 'User-G', 'noahr': 'User-H', 'nightscout1': 'User-I',
    'eli': 'User-J', 'ns_rot6': 'User-K', 'kelseyhuss': 'User-L',
}

# ── Binning parameters ──────────────────────────────────────────────────────

BIN_EDGES = list(range(72, 201, 10))  # 72, 82, 92, ..., 192, 202
BIN_CENTRES = [e + 5 for e in BIN_EDGES[:-1]]  # 77, 87, ..., 195
MIN_BIN_N = 5

# ── Model curves ────────────────────────────────────────────────────────────

def quartic(g):
    return 272 - 3.121 * g + 0.01511 * g**2 - 3.305e-5 * g**3 + 2.69e-8 * g**4

Q_REF = quartic(100)

def quartic_ratio(g):
    return quartic(g) / Q_REF

def sigmoid_ratio(g, target=100):
    return math.log(target / 120 + 1) / math.log(max(g, 40) / 120 + 1)

def pop_ratio(g):
    points = [(76, 1.15), (100, 1.00), (130, 0.80), (170, 0.70)]
    gs = [p[0] for p in points]; rs = [p[1] for p in points]
    if g <= gs[0]: return rs[0]
    if g >= gs[-1]: return rs[-1]
    for i in range(len(gs) - 1):
        if gs[i] <= g <= gs[i + 1]:
            t = (g - gs[i]) / (gs[i + 1] - gs[i])
            return rs[i] + t * (rs[i + 1] - rs[i])
    return 1.0


# ── Load data ───────────────────────────────────────────────────────────────

with open(TRIO_CACHE, 'rb') as f:
    trio_sites = pickle.load(f)

with open(BOOST_CACHE, 'rb') as f:
    boost_cache = pickle.load(f)


def build_patient(name, model, bg, isf_actual, pred_drop, actual_bg_end, hour=None):
    """Compute ISF_observed and apply quality filters. Returns dict or None."""
    actual_drop = bg - actual_bg_end
    n_total = len(bg)

    # Filter 1: |actual_drop| > 10
    m1 = np.abs(actual_drop) > 10
    # Filter 2: |pred_drop| > 5
    m2 = np.abs(pred_drop) > 5
    # Filter 3: ratio in (0.1, 5) — same direction, reasonable
    ratio = np.where(m2, actual_drop / pred_drop, 0.0)
    m3 = (ratio > 0.1) & (ratio < 5.0)
    # Combined mask so far
    mask = m1 & m2 & m3

    isf_obs = np.where(mask, isf_actual * ratio, np.nan)

    # Filter 4: ISF_observed between 3 and 600
    m4 = (isf_obs >= 3) & (isf_obs <= 600)
    mask = mask & m4

    n_after_f1 = int(m1.sum())
    n_after_f2 = int((m1 & m2).sum())
    n_after_f3 = int((m1 & m2 & m3).sum())
    n_final = int(mask.sum())

    if n_final < MIN_BIN_N:
        return None

    return {
        'name': name, 'model': model,
        'bg': bg[mask], 'isf_obs': isf_obs[mask],
        'hour': hour[mask] if hour is not None else None,
        'n_total': n_total, 'n_after_f1': n_after_f1,
        'n_after_f2': n_after_f2, 'n_after_f3': n_after_f3,
        'n_final': n_final,
    }


def build_all_patients(period='allday'):
    """Build patient list from both caches for a given period."""
    patients = []

    for s in trio_sites:
        data = s.get(period)
        if data is None or data['n'] < 10:
            continue
        raw_name = s['name']
        anon_name = ANON_MAP.get(raw_name, raw_name)
        p = build_patient(
            name=anon_name, model=s['model'],
            bg=data['bg'], isf_actual=data['isf_actual'],
            pred_drop=data['pred_drop'], actual_bg_end=data['actual_bg_end'],
            hour=data.get('hour'),
        )
        if p is not None:
            patients.append(p)

    boost_df = boost_cache.get(period)
    if boost_df is not None and len(boost_df) >= 10:
        p = build_patient(
            name='User-M', model='AAPS',
            bg=boost_df['bg'].values.astype(float),
            isf_actual=boost_df['variable_sens'].values.astype(float),
            pred_drop=boost_df['pred_drop'].values.astype(float),
            actual_bg_end=boost_df['actual_bg_end'].values.astype(float),
            hour=boost_df['hour'].values.astype(int),
        )
        if p is not None:
            patients.append(p)

    return patients


def bin_patient(bg, values):
    """Bin values by bg into 10 mg/dL bins. Returns medians, q25, q75, counts."""
    medians = np.full(len(BIN_CENTRES), np.nan)
    q25 = np.full(len(BIN_CENTRES), np.nan)
    q75 = np.full(len(BIN_CENTRES), np.nan)
    counts = np.zeros(len(BIN_CENTRES), dtype=int)

    for i, (lo, hi) in enumerate(zip(BIN_EDGES[:-1], BIN_EDGES[1:])):
        mask = (bg >= lo) & (bg < hi)
        n = mask.sum()
        counts[i] = n
        if n >= MIN_BIN_N:
            vals = values[mask]
            medians[i] = np.median(vals)
            q25[i] = np.percentile(vals, 25)
            q75[i] = np.percentile(vals, 75)

    return medians, q25, q75, counts


def normalise_curve(medians):
    """Normalise binned medians by the value at ~100 mg/dL (bin centre 97)."""
    # Find bin containing 100 mg/dL
    ref_idx = None
    for i, c in enumerate(BIN_CENTRES):
        if 95 <= c <= 105:
            ref_idx = i
            break
    if ref_idx is None or np.isnan(medians[ref_idx]):
        # Fall back: interpolate from neighbours
        for i, c in enumerate(BIN_CENTRES):
            if 90 <= c <= 110 and not np.isnan(medians[i]):
                ref_idx = i
                break
    if ref_idx is None or np.isnan(medians[ref_idx]):
        return None, None
    ref_val = medians[ref_idx]
    ratio = medians / ref_val
    return ratio, ref_val


def rmse_vs_model(ratio_curve, model_fn):
    """RMSE between empirical ratio curve and model ratio, ignoring NaN bins."""
    model_vals = np.array([model_fn(c) for c in BIN_CENTRES])
    valid = ~np.isnan(ratio_curve)
    if valid.sum() < 3:
        return np.nan
    return np.sqrt(np.mean((ratio_curve[valid] - model_vals[valid])**2))


# ══════════════════════════════════════════════════════════════════════════════
# BUILD DATA
# ══════════════════════════════════════════════════════════════════════════════

patients = build_all_patients('allday')
patients.sort(key=lambda p: p['name'])

print("=" * 90)
print("FILTER CASCADE")
print("=" * 90)
for p in patients:
    print(f"  {p['name']:8s}  total={p['n_total']:5d}  "
          f"|drop|>10={p['n_after_f1']:5d}  "
          f"|pred|>5={p['n_after_f2']:5d}  "
          f"ratio_ok={p['n_after_f3']:5d}  "
          f"final={p['n_final']:5d}")
print(f"\n  TOTAL: {sum(p['n_total'] for p in patients):,} -> "
      f"{sum(p['n_final'] for p in patients):,} samples")

# Bin each patient
for p in patients:
    meds, q25, q75, cts = bin_patient(p['bg'], p['isf_obs'])
    p['bin_medians'] = meds
    p['bin_q25'] = q25
    p['bin_q75'] = q75
    p['bin_counts'] = cts
    ratio, ref_val = normalise_curve(meds)
    p['ratio'] = ratio
    p['isf_at_100'] = ref_val

# Compute per-patient RMSE vs each model
flat_fn = lambda g: 1.0
models = {
    'quartic': quartic_ratio,
    'sigmoid': sigmoid_ratio,
    'pop-ratio': pop_ratio,
    'flat': flat_fn,
}

print("\n" + "=" * 90)
print("PER-PATIENT SUMMARY")
print("=" * 90)
header = (f"{'Name':8s}  {'Model':7s}  {'N':>5s}  {'ISF@100':>7s}  "
          f"{'RMSE_q':>7s}  {'RMSE_s':>7s}  {'RMSE_p':>7s}  {'RMSE_f':>7s}")
print(header)
print("-" * len(header))

for p in patients:
    rmses = {}
    for mname, mfn in models.items():
        rmses[mname] = rmse_vs_model(p['ratio'], mfn) if p['ratio'] is not None else np.nan
    p['rmses'] = rmses
    isf_str = f"{p['isf_at_100']:.0f}" if p['isf_at_100'] is not None else 'N/A'
    print(f"  {p['name']:8s}  {p['model']:7s}  {p['n_final']:5d}  {isf_str:>7s}  "
          f"{rmses['quartic']:7.3f}  {rmses['sigmoid']:7.3f}  "
          f"{rmses['pop-ratio']:7.3f}  {rmses['flat']:7.3f}")

# Model ranking
print("\n" + "=" * 90)
print("MODEL RANKING (weighted mean RMSE, weighted by n_final)")
print("=" * 90)
for mname in ['quartic', 'sigmoid', 'pop-ratio', 'flat']:
    weights = []; vals = []
    for p in patients:
        r = p['rmses'].get(mname, np.nan)
        if not np.isnan(r):
            weights.append(p['n_final'])
            vals.append(r)
    if weights:
        wmean = np.average(vals, weights=weights)
        print(f"  {mname:12s}  weighted RMSE = {wmean:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 1: Raw ISF scatter + binned medians (4x4 grid)
# ══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(4, 4, figsize=(16, 14), constrained_layout=True)
axes_flat = axes.flatten()

for idx, p in enumerate(patients):
    ax = axes_flat[idx]
    ax.scatter(p['bg'], p['isf_obs'], s=3, alpha=0.05, color='grey', rasterized=True)

    valid = ~np.isnan(p['bin_medians'])
    centres = np.array(BIN_CENTRES)
    ax.plot(centres[valid], p['bin_medians'][valid], 'b-', lw=1.5, zorder=5)
    ax.fill_between(centres[valid], p['bin_q25'][valid], p['bin_q75'][valid],
                     alpha=0.25, color='cornflowerblue', zorder=4)

    ax.set_title(f"{p['name']}  n={p['n_final']}  ({p['model']})", fontsize=9)
    ax.set_xlim(65, 210)
    y_upper = min(np.nanpercentile(p['isf_obs'], 99), 400)
    ax.set_ylim(0, y_upper)
    ax.tick_params(labelsize=7)
    if idx >= 12:
        ax.set_xlabel('BG (mg/dL)', fontsize=8)
    if idx % 4 == 0:
        ax.set_ylabel('ISF_obs (mg/dL/U)', fontsize=8)

# Hide unused panels
for idx in range(len(patients), 16):
    axes_flat[idx].set_visible(False)

fig.suptitle('Empirical ISF vs Glucose — Raw Observations + Binned Medians (IQR)', fontsize=13)
fig.savefig(OUT_DIR / 'empirical_isf_raw.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"\nSaved: {OUT_DIR / 'empirical_isf_raw.png'}")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 2: Normalised ratio curves (all patients + model curves)
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 6))

cmap = plt.cm.tab20
centres = np.array(BIN_CENTRES)
g_smooth = np.linspace(72, 200, 200)

for idx, p in enumerate(patients):
    if p['ratio'] is None:
        continue
    valid = ~np.isnan(p['ratio'])
    colour = cmap(idx / max(len(patients) - 1, 1))
    ax.plot(centres[valid], p['ratio'][valid], '-', lw=1.0, color=colour,
            alpha=0.7, label=p['name'])

# Model curves
ax.plot(g_smooth, [quartic_ratio(g) for g in g_smooth], 'g--', lw=2.5,
        label='Quartic (Diabeloop)', zorder=10)
ax.plot(g_smooth, [sigmoid_ratio(g) for g in g_smooth], 'b--', lw=2.5,
        label='Sigmoid (Trio)', zorder=10)
ax.plot(g_smooth, [pop_ratio(g) for g in g_smooth], color='orange', ls='--', lw=2.5,
        label='Pop-ratio (piecewise)', zorder=10)
ax.axhline(1.0, color='grey', ls=':', lw=0.8, alpha=0.5)

ax.set_xlabel('Starting Glucose (mg/dL)')
ax.set_ylabel('ISF Ratio (normalised to ISF@100)')
ax.set_title('Empirical ISF-Glucose Ratio Curves — All Patients vs Models')
ax.set_xlim(72, 200)
ax.set_ylim(0.3, 2.0)
ax.legend(fontsize=7, ncol=3, loc='upper right')
ax.grid(True, alpha=0.3)

fig.savefig(OUT_DIR / 'empirical_isf_normalised.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {OUT_DIR / 'empirical_isf_normalised.png'}")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 3: Model fit — RMSE bars + population mean curve
# ══════════════════════════════════════════════════════════════════════════════

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1.2, 1]})

# Panel A: Grouped bar chart of per-patient RMSE
model_names = ['quartic', 'sigmoid', 'pop-ratio', 'flat']
model_colours = ['green', 'blue', 'orange', 'grey']
n_models = len(model_names)
n_patients = len(patients)
x = np.arange(n_patients)
bar_w = 0.8 / n_models

for mi, (mname, mcol) in enumerate(zip(model_names, model_colours)):
    rmse_vals = [p['rmses'].get(mname, 0) for p in patients]
    ax_a.bar(x + mi * bar_w, rmse_vals, bar_w, label=mname, color=mcol, alpha=0.75)

ax_a.set_xticks(x + bar_w * (n_models - 1) / 2)
ax_a.set_xticklabels([p['name'] for p in patients], rotation=45, ha='right', fontsize=8)
ax_a.set_ylabel('RMSE (ratio units)')
ax_a.set_title('A. Per-Patient RMSE: Empirical vs Model Curves')
ax_a.legend(fontsize=8)
ax_a.grid(axis='y', alpha=0.3)

# Panel B: Population mean empirical ratio ± SD + model curves
all_ratios = []
for p in patients:
    if p['ratio'] is not None:
        all_ratios.append(p['ratio'])

ratio_stack = np.array(all_ratios)  # (n_patients, n_bins)
pop_mean = np.nanmean(ratio_stack, axis=0)
pop_sd = np.nanstd(ratio_stack, axis=0)

valid = ~np.isnan(pop_mean)
ax_b.plot(centres[valid], pop_mean[valid], 'k-', lw=2, label='Population mean')
ax_b.fill_between(centres[valid],
                   (pop_mean - pop_sd)[valid], (pop_mean + pop_sd)[valid],
                   alpha=0.2, color='black', label='±1 SD')

g_smooth = np.linspace(72, 200, 200)
# Compute weighted RMSE for annotation
for mname, mfn, mcol, mls in [
    ('Quartic', quartic_ratio, 'green', '--'),
    ('Sigmoid', sigmoid_ratio, 'blue', '--'),
    ('Pop-ratio', pop_ratio, 'orange', '--'),
]:
    model_at_centres = np.array([mfn(c) for c in BIN_CENTRES])
    wrmse = np.sqrt(np.nanmean((pop_mean[valid] - model_at_centres[valid])**2))
    ax_b.plot(g_smooth, [mfn(g) for g in g_smooth], ls=mls, color=mcol, lw=2,
              label=f'{mname} (RMSE={wrmse:.3f})')

ax_b.axhline(1.0, color='grey', ls=':', lw=0.8, alpha=0.5)
ax_b.set_xlabel('Starting Glucose (mg/dL)')
ax_b.set_ylabel('ISF Ratio')
ax_b.set_title('B. Population Mean Empirical Curve vs Models')
ax_b.set_xlim(72, 200)
ax_b.set_ylim(0.3, 1.8)
ax_b.legend(fontsize=8)
ax_b.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(OUT_DIR / 'empirical_isf_model_fit.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {OUT_DIR / 'empirical_isf_model_fit.png'}")


# ══════════════════════════════════════════════════════════════════════════════
# CHART 4: Overnight vs Daytime
# ══════════════════════════════════════════════════════════════════════════════

patients_on = build_all_patients('overnight')
patients_day = build_all_patients('daytime')

# Bin and normalise each period
def get_ratio_stack(plist):
    """Bin, normalise, return (name->ratio) dict and ratio matrix."""
    name_ratio = {}
    for p in plist:
        meds, _, _, _ = bin_patient(p['bg'], p['isf_obs'])
        ratio, _ = normalise_curve(meds)
        if ratio is not None:
            name_ratio[p['name']] = ratio
    return name_ratio

on_ratios = get_ratio_stack(patients_on)
day_ratios = get_ratio_stack(patients_day)

# Find patients present in both
common_names = sorted(set(on_ratios.keys()) & set(day_ratios.keys()))
on_stack = np.array([on_ratios[n] for n in common_names])
day_stack = np.array([day_ratios[n] for n in common_names])

on_mean = np.nanmean(on_stack, axis=0)
on_sd = np.nanstd(on_stack, axis=0)
day_mean = np.nanmean(day_stack, axis=0)
day_sd = np.nanstd(day_stack, axis=0)

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Population mean curves
valid_on = ~np.isnan(on_mean)
valid_day = ~np.isnan(day_mean)

ax_a.plot(centres[valid_on], on_mean[valid_on], 'b-', lw=2, label='Overnight')
ax_a.fill_between(centres[valid_on],
                   (on_mean - on_sd)[valid_on], (on_mean + on_sd)[valid_on],
                   alpha=0.15, color='blue')
ax_a.plot(centres[valid_day], day_mean[valid_day], 'r-', lw=2, label='Daytime')
ax_a.fill_between(centres[valid_day],
                   (day_mean - day_sd)[valid_day], (day_mean + day_sd)[valid_day],
                   alpha=0.15, color='red')
ax_a.axhline(1.0, color='grey', ls=':', lw=0.8)
ax_a.set_xlabel('Starting Glucose (mg/dL)')
ax_a.set_ylabel('ISF Ratio')
ax_a.set_title('A. Population Mean: Overnight vs Daytime')
ax_a.set_xlim(72, 200)
ax_a.set_ylim(0.3, 1.8)
ax_a.legend()
ax_a.grid(True, alpha=0.3)

# Panel B: Per-patient delta (daytime - overnight)
for i, name in enumerate(common_names):
    delta = day_ratios[name] - on_ratios[name]
    valid = ~np.isnan(delta)
    colour = cmap(i / max(len(common_names) - 1, 1))
    ax_b.plot(centres[valid], delta[valid], '-', lw=1, color=colour, alpha=0.7, label=name)

ax_b.axhline(0, color='grey', ls=':', lw=0.8)
ax_b.set_xlabel('Starting Glucose (mg/dL)')
ax_b.set_ylabel('Delta (Daytime - Overnight)')
ax_b.set_title('B. Per-Patient: Daytime minus Overnight Ratio')
ax_b.set_xlim(72, 200)
ax_b.legend(fontsize=7, ncol=2, loc='best')
ax_b.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig(OUT_DIR / 'empirical_isf_overnight_vs_daytime.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f"Saved: {OUT_DIR / 'empirical_isf_overnight_vs_daytime.png'}")

# Overnight vs daytime slope difference
print("\n" + "=" * 90)
print("OVERNIGHT VS DAYTIME: SLOPE DIFFERENCE")
print("=" * 90)
slopes_on = []; slopes_day = []
for name in common_names:
    for ratios, slope_list in [(on_ratios[name], slopes_on), (day_ratios[name], slopes_day)]:
        valid_idx = [i for i in range(len(BIN_CENTRES)) if not np.isnan(ratios[i])]
        if len(valid_idx) >= 3:
            xs = np.array([BIN_CENTRES[i] for i in valid_idx])
            ys = np.array([ratios[i] for i in valid_idx])
            slope = np.polyfit(xs, ys, 1)[0]
            slope_list.append(slope)
        else:
            slope_list.append(np.nan)

slopes_on = np.array(slopes_on); slopes_day = np.array(slopes_day)
valid = ~(np.isnan(slopes_on) | np.isnan(slopes_day))
if valid.sum() > 0:
    mean_diff = np.mean(slopes_day[valid] - slopes_on[valid])
    print(f"  Mean slope (overnight): {np.mean(slopes_on[valid]):.5f} per mg/dL")
    print(f"  Mean slope (daytime):   {np.mean(slopes_day[valid]):.5f} per mg/dL")
    print(f"  Mean slope difference (daytime - overnight): {mean_diff:.5f} per mg/dL")

print("\nDone.")
