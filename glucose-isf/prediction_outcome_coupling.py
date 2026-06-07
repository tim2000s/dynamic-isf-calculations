#!/usr/bin/env python3
"""
Paper 9 – "How Tightly Are Predictions Linked to Outcomes?"
============================================================
Multi-step feedback replay to quantify the coupling between ISF prediction
accuracy and closed-loop glycaemic outcomes.  Extends the single-step replay
from Paper 8 with sub-step feedback correction that models how AID systems
re-evaluate and adjust dosing as glucose moves.

Four levels of analysis:
  1. Model-Site aggregate (MAE vs TIR/TITR scatter + correlation)
  2. Per-sample prediction error vs outcome deviation (density)
  3. Error direction (over- vs under-prediction)
  4. Magnitude threshold (bins of prediction error)

Plus: dampening quantification and Boost sequential validation (appendix).
"""

import pickle, math, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from scipy import stats

warnings.filterwarnings('ignore', category=RuntimeWarning)

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


MODEL_DEFS = {
    'Original': lambda s: s['isf_actual'],
    'Quartic + 1800/TDD': lambda s: np.array([s['isf_tdd'] * quartic(g) / Q_REF for g in s['bg']]),
    'Profile + quartic': lambda s: np.array([s['isf_true'] * quartic(g) / Q_REF for g in s['bg']]),
    'Profile + pop-ratio': lambda s: np.array([s['isf_true'] * ratio_fn(g) for g in s['bg']]),
    'Profile + sigmoid': lambda s: np.array([s['isf_true'] * sigmoid_ratio(g) for g in s['bg']]),
    'Profile + flat': lambda s: np.full(len(s['bg']), s['isf_true']),
    '1800/TDD + flat': lambda s: np.full(len(s['bg']), s['isf_tdd']),
}

MODEL_SHORT = {
    'Original': 'Original',
    'Quartic + 1800/TDD': 'Q+TDD',
    'Profile + quartic': 'Prof+Q',
    'Profile + pop-ratio': 'Prof+Pop',
    'Profile + sigmoid': 'Prof+Sig',
    'Profile + flat': 'Prof+Flat',
    '1800/TDD + flat': 'TDD+Flat',
}

MODEL_COLORS = {
    'Original': '#555555',
    'Quartic + 1800/TDD': '#2196F3',
    'Profile + quartic': '#4CAF50',
    'Profile + pop-ratio': '#FF9800',
    'Profile + sigmoid': '#9C27B0',
    'Profile + flat': '#F44336',
    '1800/TDD + flat': '#795548',
}


# ══════════════════════════════════════════════════════════════════════════════
# INSULIN ACTIVITY CURVE  (bilinear, DIA=300 min, peak at DIA/3)
# ══════════════════════════════════════════════════════════════════════════════

def iob_fraction_remaining(t_min, DIA=300):
    """Fraction of IOB remaining at time t_min (bilinear activity curve).
    Activity: linear ramp 0 -> peak at DIA/3, linear ramp peak -> 0 at DIA.
    IOB_remaining = 1 - integral(activity, 0, t) / integral(activity, 0, DIA)."""
    peak = DIA / 3.0
    # Total area under activity curve = 0.5 * peak * h + 0.5 * (DIA-peak) * h
    # = 0.5 * h * DIA  (where h = peak height). We can set h=1, total = DIA/2.
    # Integral from 0 to t:
    if t_min <= 0:
        return 1.0
    if t_min >= DIA:
        return 0.0
    # Activity height at peak: h = 2/DIA (normalised so integral = 1)
    h = 2.0 / DIA
    if t_min <= peak:
        # Ramp up: activity = h * t / peak
        # Integral = h * t^2 / (2*peak)
        acted = h * t_min**2 / (2.0 * peak)
    else:
        # Full ramp-up integral + ramp-down portion
        up_integral = h * peak / 2.0
        # Ramp down: activity at peak+dt = h * (1 - dt/(DIA-peak))
        dt = t_min - peak
        down_len = DIA - peak
        # Integral of ramp down from 0 to dt = h * (dt - dt^2/(2*down_len))
        down_integral = h * (dt - dt**2 / (2.0 * down_len))
        acted = up_integral + down_integral
    return max(0.0, 1.0 - acted)


def insulin_acted_fraction(t_start, t_end, DIA=300):
    """Fraction of a bolus that acts between t_start and t_end minutes."""
    return iob_fraction_remaining(t_start, DIA) - iob_fraction_remaining(t_end, DIA)


# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA & BUILD SITES
# ══════════════════════════════════════════════════════════════════════════════

print("Loading caches ...")
with open(TRIO_CACHE, 'rb') as f:
    trio_sites = pickle.load(f)
with open(BOOST_CACHE, 'rb') as f:
    boost_cache = pickle.load(f)


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
            'ts': boost_df['ts'].values,  # keep timestamps for sequential validation
            'isf_true': isf100, 'isf_tdd': 1800 / tdd,
        })
    return sites


# ══════════════════════════════════════════════════════════════════════════════
# REPLAY ENGINES: SINGLE-STEP & MULTI-STEP
# ══════════════════════════════════════════════════════════════════════════════

TARGET = 100  # mg/dL


def replay_single_step(site, isf_test_arr):
    """Single-step replay (same as Paper 8)."""
    bg = site['bg']
    isf_actual = site['isf_actual']
    pred_drop = site['pred_drop']
    actual_bg_end = site['actual_bg_end']

    iob = pred_drop / isf_actual
    eventual_orig = bg - iob * isf_actual
    eventual_test = bg - iob * isf_test_arr
    corr_orig = (eventual_orig - TARGET) / isf_actual
    corr_test = (eventual_test - TARGET) / isf_test_arr
    delta_insulin = np.clip(corr_test - corr_orig, -2.0, 2.0)

    actual_drop = bg - actual_bg_end
    true_isf = np.where(np.abs(iob) > 0.01, actual_drop / iob, isf_actual)
    true_isf = np.clip(true_isf, 5, 500)

    simulated_end = np.clip(actual_bg_end - delta_insulin * true_isf, 40, 400)
    return simulated_end, delta_insulin, true_isf, iob


def replay_multi_step(site, isf_test_arr, n_steps=3):
    """Multi-step feedback replay with sub-step BG correction.

    At each sub-step the system re-evaluates with updated BG, leading to
    dampened total delta insulin compared to single-step.
    """
    bg = site['bg']
    isf_actual = site['isf_actual']
    pred_drop = site['pred_drop']
    actual_bg_end = site['actual_bg_end']
    n_samples = len(bg)

    DIA = 300  # minutes
    window_minutes = 240  # 4-hour window
    step_duration = window_minutes / n_steps

    # Infer IOB and true ISF (same as single-step)
    iob_initial = pred_drop / isf_actual
    actual_drop = bg - actual_bg_end
    true_isf = np.where(np.abs(iob_initial) > 0.01,
                        actual_drop / iob_initial, isf_actual)
    true_isf = np.clip(true_isf, 5, 500)

    # Per-sample multi-step loop
    total_delta_insulin = np.zeros(n_samples)
    bg_current = bg.copy().astype(float)

    # Track delta insulin delivered at each prior sub-step so we can
    # compute how much of it is still on board (IOB) at later steps.
    delta_history = []  # list of (delivery_time, delta_amount) tuples

    for step in range(n_steps):
        t_start = step * step_duration

        # IOB from original delivery remaining at this sub-step
        iob_orig_remaining = iob_initial * iob_fraction_remaining(t_start, DIA)

        # IOB from our own prior delta corrections still on board
        iob_delta_remaining = np.zeros(n_samples)
        for (t_delivered, delta_amt) in delta_history:
            elapsed = t_start - t_delivered
            iob_delta_remaining += delta_amt * iob_fraction_remaining(elapsed, DIA)

        # Total IOB the system sees: original remaining + our prior deltas
        total_iob_test = iob_orig_remaining + iob_delta_remaining
        total_iob_orig = iob_orig_remaining  # original system has no deltas

        # What the test ISF predicts at current BG
        eventual_test = bg_current - total_iob_test * isf_test_arr
        corr_test = (eventual_test - TARGET) / isf_test_arr

        # What the original ISF predicts at current BG
        eventual_orig = bg_current - total_iob_orig * isf_actual
        corr_orig = (eventual_orig - TARGET) / isf_actual

        # Delta insulin this sub-step (rate-limited)
        delta_step = np.clip(corr_test - corr_orig, -0.67, 0.67)
        total_delta_insulin += delta_step
        delta_history.append((t_start, delta_step))

        # BG evolves: delta insulin acts with true ISF over this sub-step
        action_frac = insulin_acted_fraction(0, step_duration, DIA)
        bg_current = bg_current - delta_step * true_isf * action_frac
        bg_current = np.clip(bg_current, 40, 400)

    # Final simulated end: actual outcome adjusted by total delta insulin
    simulated_end = np.clip(actual_bg_end - total_delta_insulin * true_isf, 40, 400)
    return simulated_end, total_delta_insulin, true_isf


def compute_metrics(bg_arr):
    n = len(bg_arr)
    if n == 0:
        return {'tir': np.nan, 'titr': np.nan, 'below70': np.nan,
                'below54': np.nan, 'above140': np.nan, 'mean': np.nan,
                'sd': np.nan, 'n': 0}
    return {
        'tir': np.mean((bg_arr >= 70) & (bg_arr <= 180)) * 100,
        'titr': np.mean((bg_arr >= 70) & (bg_arr <= 140)) * 100,
        'below70': np.mean(bg_arr < 70) * 100,
        'below54': np.mean(bg_arr < 54) * 100,
        'above140': np.mean(bg_arr > 140) * 100,
        'mean': np.mean(bg_arr), 'sd': np.std(bg_arr), 'n': n,
    }


# ══════════════════════════════════════════════════════════════════════════════
# RUN ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

PERIODS = ['allday', 'overnight', 'daytime']

all_results = {}  # period -> {model -> {site_name -> dict of arrays}}

for period in PERIODS:
    sites = build_sites(trio_sites, boost_cache, period)
    sites = [s for s in sites if not np.isnan(s['isf_true'])]

    print(f"\n{'='*80}")
    print(f"Period: {period.upper()} — {len(sites)} sites, "
          f"{sum(s['n'] for s in sites):,} samples")
    print(f"{'='*80}")

    period_results = {}
    for model_name, model_fn in MODEL_DEFS.items():
        period_results[model_name] = {}
        for s in sites:
            isf_test = model_fn(s)
            sim_end_ss, delta_ins_ss, true_isf, iob = replay_single_step(s, isf_test)
            sim_end_ms, delta_ins_ms, _ = replay_multi_step(s, isf_test, n_steps=3)

            # Counterfactual prediction: what the model would predict as end BG
            counterfactual_end = s['bg'] - (s['pred_drop'] / s['isf_actual']) * isf_test
            counterfactual_mae = np.mean(np.abs(counterfactual_end - s['actual_bg_end']))

            period_results[model_name][s['name']] = {
                'sim_end_ss': sim_end_ss,
                'sim_end_ms': sim_end_ms,
                'delta_ins_ss': delta_ins_ss,
                'delta_ins_ms': delta_ins_ms,
                'true_isf': true_isf,
                'iob': iob,
                'counterfactual_end': counterfactual_end,
                'counterfactual_mae': counterfactual_mae,
                'metrics_actual': compute_metrics(s['actual_bg_end']),
                'metrics_ss': compute_metrics(sim_end_ss),
                'metrics_ms': compute_metrics(sim_end_ms),
            }

    all_results[period] = {'sites': sites, 'models': period_results}


# ══════════════════════════════════════════════════════════════════════════════
# LEVEL 1: MODEL-SITE SCATTER  (MAE vs TIR / TITR)
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("LEVEL 1: Prediction Accuracy vs Glycaemic Outcomes (Model-Site)")
print("=" * 90)

for period in PERIODS:
    pr = all_results[period]
    sites = pr['sites']
    models = pr['models']

    print(f"\n--- {period.upper()} ---")
    print(f"{'Metric Pair':35s} {'Pearson r':>10s} {'p-val':>10s} {'Spearman r':>10s} {'p-val':>10s}")
    print("-" * 80)

    # Collect (MAE, TIR, TITR, below70) for each model-site point
    mae_list, tir_ss, titr_ss, b70_ss = [], [], [], []
    tir_ms, titr_ms, b70_ms = [], [], []
    model_labels = []

    for model_name in MODEL_DEFS:
        for s in sites:
            r = models[model_name][s['name']]
            mae_list.append(r['counterfactual_mae'])
            tir_ss.append(r['metrics_ss']['tir'])
            titr_ss.append(r['metrics_ss']['titr'])
            b70_ss.append(r['metrics_ss']['below70'])
            tir_ms.append(r['metrics_ms']['tir'])
            titr_ms.append(r['metrics_ms']['titr'])
            b70_ms.append(r['metrics_ms']['below70'])
            model_labels.append(model_name)

    mae_arr = np.array(mae_list)
    pairs = [
        ('MAE vs TIR (single-step)', mae_arr, np.array(tir_ss)),
        ('MAE vs TITR (single-step)', mae_arr, np.array(titr_ss)),
        ('MAE vs <70% (single-step)', mae_arr, np.array(b70_ss)),
        ('MAE vs TIR (multi-step)', mae_arr, np.array(tir_ms)),
        ('MAE vs TITR (multi-step)', mae_arr, np.array(titr_ms)),
        ('MAE vs <70% (multi-step)', mae_arr, np.array(b70_ms)),
    ]
    for label, x, y in pairs:
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 5:
            continue
        pr_r, pr_p = stats.pearsonr(x[mask], y[mask])
        sp_r, sp_p = stats.spearmanr(x[mask], y[mask])
        print(f"  {label:35s} {pr_r:10.3f} {pr_p:10.4f} {sp_r:10.3f} {sp_p:10.4f}")

    # Store for chart
    all_results[period]['l1_mae'] = mae_arr
    all_results[period]['l1_tir_ss'] = np.array(tir_ss)
    all_results[period]['l1_titr_ss'] = np.array(titr_ss)
    all_results[period]['l1_tir_ms'] = np.array(tir_ms)
    all_results[period]['l1_titr_ms'] = np.array(titr_ms)
    all_results[period]['l1_labels'] = model_labels


# ══════════════════════════════════════════════════════════════════════════════
# LEVEL 2: PER-SAMPLE CORRELATION
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("LEVEL 2: Per-Sample Prediction Error vs Outcome Deviation")
print("=" * 90)

ALL_MODELS = list(MODEL_DEFS.keys())
FOCUS_MODELS = ['Profile + pop-ratio', 'Quartic + 1800/TDD']  # for detailed charts
period = 'allday'
pr = all_results[period]
sites = pr['sites']
models = pr['models']

l2_data = {}  # model -> (pred_error, dev_ss, dev_ms)

for model_name in ALL_MODELS:
    pred_errors, dev_ss, dev_ms = [], [], []
    for s in sites:
        r = models[model_name][s['name']]
        pe = np.abs(r['counterfactual_end'] - s['actual_bg_end'])
        dss = np.abs(r['sim_end_ss'] - s['actual_bg_end'])
        dms = np.abs(r['sim_end_ms'] - s['actual_bg_end'])
        pred_errors.append(pe)
        dev_ss.append(dss)
        dev_ms.append(dms)

    pred_errors = np.concatenate(pred_errors)
    dev_ss = np.concatenate(dev_ss)
    dev_ms = np.concatenate(dev_ms)
    l2_data[model_name] = (pred_errors, dev_ss, dev_ms)

    mask = np.isfinite(pred_errors) & np.isfinite(dev_ss) & np.isfinite(dev_ms)
    pr_ss, _ = stats.pearsonr(pred_errors[mask], dev_ss[mask])
    pr_ms, _ = stats.pearsonr(pred_errors[mask], dev_ms[mask])
    sp_ss, _ = stats.spearmanr(pred_errors[mask], dev_ss[mask])
    sp_ms, _ = stats.spearmanr(pred_errors[mask], dev_ms[mask])

    print(f"\n  {model_name} (n={mask.sum():,}):")
    print(f"    Single-step: Pearson r={pr_ss:.3f}, Spearman r={sp_ss:.3f}")
    print(f"    Multi-step:  Pearson r={pr_ms:.3f}, Spearman r={sp_ms:.3f}")


# ══════════════════════════════════════════════════════════════════════════════
# LEVEL 3: ERROR DIRECTION
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("LEVEL 3: Outcome Deviation by Prediction Error Direction")
print("=" * 90)

l3_data = {}

for model_name in ALL_MODELS:
    signed_errors, dev_ss_signed, dev_ms_signed = [], [], []
    for s in sites:
        r = models[model_name][s['name']]
        # Signed: counterfactual_end - actual_end
        # Negative = over-predicting drop (predicted lower than actual)
        se = r['counterfactual_end'] - s['actual_bg_end']
        dss = r['sim_end_ss'] - s['actual_bg_end']  # signed deviation
        dms = r['sim_end_ms'] - s['actual_bg_end']
        signed_errors.append(se)
        dev_ss_signed.append(dss)
        dev_ms_signed.append(dms)

    signed_errors = np.concatenate(signed_errors)
    dev_ss_signed = np.concatenate(dev_ss_signed)
    dev_ms_signed = np.concatenate(dev_ms_signed)
    l3_data[model_name] = (signed_errors, dev_ss_signed, dev_ms_signed)

    over = signed_errors < 0  # over-predicting drop
    under = signed_errors > 0  # under-predicting drop

    print(f"\n  {model_name}:")
    print(f"    {'Direction':25s} {'N':>8s} {'Mean dev SS':>12s} {'Mean dev MS':>12s} {'|Mean| SS':>12s} {'|Mean| MS':>12s}")
    print(f"    {'-'*75}")
    for label, mask in [('Over-predict drop', over), ('Under-predict drop', under)]:
        if mask.sum() == 0:
            continue
        print(f"    {label:25s} {mask.sum():8d} "
              f"{np.mean(dev_ss_signed[mask]):12.2f} {np.mean(dev_ms_signed[mask]):12.2f} "
              f"{np.mean(np.abs(dev_ss_signed[mask])):12.2f} {np.mean(np.abs(dev_ms_signed[mask])):12.2f}")


# ══════════════════════════════════════════════════════════════════════════════
# LEVEL 4: MAGNITUDE THRESHOLD
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("LEVEL 4: Outcome Deviation by Prediction Error Magnitude")
print("=" * 90)

BINS = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 999)]
BIN_LABELS = ['0-5', '5-10', '10-20', '20-40', '40+']

l4_data = {}

for model_name in ALL_MODELS:
    pred_err, dev_ss, dev_ms = l2_data[model_name]

    print(f"\n  {model_name}:")
    print(f"    {'Bin (mg/dL)':>12s} {'N':>8s} {'Mean|dev| SS':>14s} {'Med|dev| SS':>14s} "
          f"{'Mean|dev| MS':>14s} {'Med|dev| MS':>14s}")
    print(f"    {'-'*75}")

    bin_stats = []
    for (lo, hi), lab in zip(BINS, BIN_LABELS):
        mask = (pred_err >= lo) & (pred_err < hi)
        n = mask.sum()
        if n == 0:
            bin_stats.append((lab, n, np.nan, np.nan, np.nan, np.nan))
            continue
        mean_ss = np.mean(dev_ss[mask])
        med_ss = np.median(dev_ss[mask])
        mean_ms = np.mean(dev_ms[mask])
        med_ms = np.median(dev_ms[mask])
        bin_stats.append((lab, n, mean_ss, med_ss, mean_ms, med_ms))
        print(f"    {lab:>12s} {n:8d} {mean_ss:14.2f} {med_ss:14.2f} "
              f"{mean_ms:14.2f} {med_ms:14.2f}")

    l4_data[model_name] = bin_stats


# ══════════════════════════════════════════════════════════════════════════════
# DAMPENING ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("DAMPENING: Multi-step vs Single-step Delta Insulin")
print("=" * 90)

dampening_data = {}

for model_name in MODEL_DEFS:
    if model_name == 'Original':
        continue
    all_ss, all_ms = [], []
    for s in sites:
        r = models[model_name][s['name']]
        all_ss.append(r['delta_ins_ss'])
        all_ms.append(r['delta_ins_ms'])
    all_ss = np.concatenate(all_ss)
    all_ms = np.concatenate(all_ms)

    # Dampening ratio: where |single-step| > 0.01
    mask = np.abs(all_ss) > 0.01
    if mask.sum() > 0:
        ratios = all_ms[mask] / all_ss[mask]
        dampening_data[model_name] = ratios
        print(f"  {MODEL_SHORT[model_name]:12s}: median ratio = {np.median(ratios):.3f}, "
              f"mean = {np.mean(ratios):.3f}, "
              f"IQR = [{np.percentile(ratios, 25):.3f}, {np.percentile(ratios, 75):.3f}], "
              f"n = {mask.sum():,}")


# ══════════════════════════════════════════════════════════════════════════════
# SENSITIVITY: N_STEPS = 1, 2, 3
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("SENSITIVITY: Convergence across N_STEPS = 1, 2, 3")
print("=" * 90)

test_model = 'Profile + pop-ratio'
print(f"\n  Model: {test_model}")
print(f"  {'N_STEPS':>8s} {'TIR (wt)':>10s} {'TITR (wt)':>10s} {'<70% (wt)':>10s} {'Mean delta_ins':>16s}")
print(f"  {'-'*60}")

for ns in [1, 2, 3]:
    tir_w, titr_w, b70_w, weights = [], [], [], []
    all_di = []
    for s in sites:
        isf_test = MODEL_DEFS[test_model](s)
        sim_end, di, _ = replay_multi_step(s, isf_test, n_steps=ns)
        m = compute_metrics(sim_end)
        tir_w.append(m['tir']); titr_w.append(m['titr']); b70_w.append(m['below70'])
        weights.append(s['n'])
        all_di.append(di)
    weights = np.array(weights, dtype=float)
    tir_agg = np.average(tir_w, weights=weights)
    titr_agg = np.average(titr_w, weights=weights)
    b70_agg = np.average(b70_w, weights=weights)
    mean_di = np.mean(np.abs(np.concatenate(all_di)))
    print(f"  {ns:8d} {tir_agg:10.2f} {titr_agg:10.2f} {b70_agg:10.3f} {mean_di:16.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# BOOST SEQUENTIAL VALIDATION (APPENDIX)
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("BOOST SEQUENTIAL VALIDATION (User-M)")
print("=" * 90)

boost_site = None
for s in all_results['allday']['sites']:
    if s['name'] == 'User-M':
        boost_site = s
        break

boost_seq_results = {}

if boost_site is not None and 'ts' in boost_site:
    ts = boost_site['ts']
    bg = boost_site['bg']
    actual_end = boost_site['actual_bg_end']

    # Convert timestamps to minutes for gap calculation
    ts_epoch = ts.astype('datetime64[s]').astype(float) / 60.0  # minutes

    # Find sequential pairs: gap 4-6 hours (240-360 min)
    pairs = []
    for i in range(len(ts) - 1):
        for j in range(i + 1, min(i + 50, len(ts))):
            gap = ts_epoch[j] - ts_epoch[i]
            if 240 <= gap <= 360:
                pairs.append((i, j))
                break
            elif gap > 360:
                break

    print(f"  Found {len(pairs)} sequential pairs (4-6 hour gap)")

    if len(pairs) >= 10:
        for model_name in FOCUS_MODELS:
            r = all_results['allday']['models'][model_name][boost_site['name']]
            sim_end_ss = r['sim_end_ss']
            sim_end_ms = r['sim_end_ms']

            # For each pair: first sample's simulated end vs second sample's actual start
            sim_starts_ss = np.array([sim_end_ss[i] for i, j in pairs])
            sim_starts_ms = np.array([sim_end_ms[i] for i, j in pairs])
            actual_starts = np.array([bg[j] for i, j in pairs])

            mae_ss = np.mean(np.abs(sim_starts_ss - actual_starts))
            mae_ms = np.mean(np.abs(sim_starts_ms - actual_starts))
            r_ss, _ = stats.pearsonr(sim_starts_ss, actual_starts)
            r_ms, _ = stats.pearsonr(sim_starts_ms, actual_starts)

            boost_seq_results[model_name] = {
                'sim_ss': sim_starts_ss, 'sim_ms': sim_starts_ms,
                'actual': actual_starts,
                'mae_ss': mae_ss, 'mae_ms': mae_ms,
                'r_ss': r_ss, 'r_ms': r_ms,
            }

            print(f"\n  {model_name}:")
            print(f"    Single-step: MAE = {mae_ss:.1f} mg/dL, r = {r_ss:.3f}")
            print(f"    Multi-step:  MAE = {mae_ms:.1f} mg/dL, r = {r_ms:.3f}")
    else:
        print("  Insufficient sequential pairs for validation.")
else:
    print("  Boost site not found or no timestamps.")


# ══════════════════════════════════════════════════════════════════════════════
# CHARTS
# ══════════════════════════════════════════════════════════════════════════════

print("\n\nGenerating charts ...")

# ── Chart 1: Level 1 Scatter (allday) ─────────────────────────────────────────

fig1, axes1 = plt.subplots(2, 2, figsize=(14, 11))
fig1.suptitle("Level 1: Prediction Accuracy vs Glycaemic Outcomes (All-Day)\n"
              "Each point = one ISF model applied to one site",
              fontsize=12, fontweight='bold')

ad = all_results['allday']
mae = ad['l1_mae']
labels = ad['l1_labels']

plot_pairs = [
    (axes1[0, 0], ad['l1_tir_ss'], 'TIR % (single-step)', 'A'),
    (axes1[0, 1], ad['l1_titr_ss'], 'TITR % (single-step)', 'B'),
    (axes1[1, 0], ad['l1_tir_ms'], 'TIR % (multi-step)', 'C'),
    (axes1[1, 1], ad['l1_titr_ms'], 'TITR % (multi-step)', 'D'),
]

for ax, y_arr, ylabel, panel in plot_pairs:
    for model_name in MODEL_DEFS:
        mask = np.array([l == model_name for l in labels])
        ax.scatter(mae[mask], y_arr[mask], c=MODEL_COLORS[model_name],
                   s=30, alpha=0.7, label=MODEL_SHORT[model_name], edgecolors='white',
                   linewidth=0.3)

    # Regression line
    m_ok = np.isfinite(mae) & np.isfinite(y_arr)
    if m_ok.sum() > 5:
        slope, intercept, r_val, p_val, _ = stats.linregress(mae[m_ok], y_arr[m_ok])
        x_line = np.linspace(mae[m_ok].min(), mae[m_ok].max(), 100)
        ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5, linewidth=1)
        ax.text(0.05, 0.05, f'r = {r_val:.3f} (p = {p_val:.4f})',
                transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel('Counterfactual MAE (mg/dL)')
    ax.set_ylabel(ylabel)
    ax.set_title(f'{panel}. {ylabel}')
    ax.grid(True, alpha=0.2)

# Single legend
handles, leg_labels = axes1[0, 0].get_legend_handles_labels()
fig1.legend(handles, leg_labels, loc='lower center', ncol=4, fontsize=8,
            bbox_to_anchor=(0.5, -0.02))

plt.tight_layout(rect=[0, 0.04, 1, 0.95])
fig1.savefig(OUT_DIR / 'coupling_level1_scatter.png', dpi=150, bbox_inches='tight')
print(f"  Saved: coupling_level1_scatter.png")


# ── Chart 2: Level 1 Periods (overnight + daytime) ──────────────────────────

fig2, axes2 = plt.subplots(2, 2, figsize=(14, 11))
fig2.suptitle("Level 1: Prediction Accuracy vs Outcomes — Overnight vs Daytime",
              fontsize=12, fontweight='bold')

for col, period in enumerate(['overnight', 'daytime']):
    pd_data = all_results.get(period)
    if pd_data is None or 'l1_mae' not in pd_data:
        continue
    p_mae = pd_data['l1_mae']
    p_labels = pd_data['l1_labels']

    for row, (key, ylabel) in enumerate([('l1_tir_ms', 'TIR %'), ('l1_titr_ms', 'TITR %')]):
        ax = axes2[row, col]
        y_arr = pd_data[key]
        for model_name in MODEL_DEFS:
            mask = np.array([l == model_name for l in p_labels])
            ax.scatter(p_mae[mask], y_arr[mask], c=MODEL_COLORS[model_name],
                       s=30, alpha=0.7, label=MODEL_SHORT[model_name],
                       edgecolors='white', linewidth=0.3)

        m_ok = np.isfinite(p_mae) & np.isfinite(y_arr)
        if m_ok.sum() > 5:
            slope, intercept, r_val, p_val, _ = stats.linregress(p_mae[m_ok], y_arr[m_ok])
            x_line = np.linspace(p_mae[m_ok].min(), p_mae[m_ok].max(), 100)
            ax.plot(x_line, slope * x_line + intercept, 'k--', alpha=0.5, linewidth=1)
            ax.text(0.05, 0.05, f'r = {r_val:.3f}',
                    transform=ax.transAxes, fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        ax.set_xlabel('Counterfactual MAE (mg/dL)')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{period.capitalize()}: {ylabel} (multi-step)')
        ax.grid(True, alpha=0.2)

handles, leg_labels = axes2[0, 0].get_legend_handles_labels()
fig2.legend(handles, leg_labels, loc='lower center', ncol=4, fontsize=8,
            bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.04, 1, 0.95])
fig2.savefig(OUT_DIR / 'coupling_level1_periods.png', dpi=150, bbox_inches='tight')
print(f"  Saved: coupling_level1_periods.png")


# ── Chart 3: Level 2 Density (hexbin) ────────────────────────────────────────

fig3, axes3 = plt.subplots(len(FOCUS_MODELS), 2, figsize=(14, 5 * len(FOCUS_MODELS)))
fig3.suptitle("Level 2: Per-Sample Prediction Error vs Outcome Deviation\n"
              "(Pooled across sites, all-day)",
              fontsize=12, fontweight='bold')

if len(FOCUS_MODELS) == 1:
    axes3 = axes3.reshape(1, -1)

for row, model_name in enumerate(FOCUS_MODELS):
    pred_err, dev_ss, dev_ms = l2_data[model_name]

    for col, (dev, label) in enumerate([(dev_ss, 'Single-step'), (dev_ms, 'Multi-step')]):
        ax = axes3[row, col]
        mask = np.isfinite(pred_err) & np.isfinite(dev) & (pred_err < 100) & (dev < 100)
        hb = ax.hexbin(pred_err[mask], dev[mask], gridsize=40, cmap='YlOrRd',
                       mincnt=1, linewidths=0.1)
        plt.colorbar(hb, ax=ax, label='Count')

        # Diagonal reference
        lim = max(pred_err[mask].max(), dev[mask].max())
        ax.plot([0, lim], [0, lim], 'k--', alpha=0.3, linewidth=1, label='1:1')

        pr_r, _ = stats.pearsonr(pred_err[mask], dev[mask])
        ax.text(0.05, 0.95, f'r = {pr_r:.3f}\nn = {mask.sum():,}',
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_xlabel('|Prediction Error| (mg/dL)')
        ax.set_ylabel(f'|Outcome Deviation| (mg/dL) — {label}')
        ax.set_title(f'{MODEL_SHORT[model_name]}: {label}')
        ax.grid(True, alpha=0.2)

plt.tight_layout()
fig3.savefig(OUT_DIR / 'coupling_level2_density.png', dpi=150, bbox_inches='tight')
print(f"  Saved: coupling_level2_density.png")


# ── Chart 4: Level 3 Direction ────────────────────────────────────────────────

fig4, axes4 = plt.subplots(1, len(FOCUS_MODELS), figsize=(7 * len(FOCUS_MODELS), 6))
fig4.suptitle("Level 3: Outcome Deviation by Prediction Error Direction",
              fontsize=12, fontweight='bold')

if len(FOCUS_MODELS) == 1:
    axes4 = [axes4]

for idx, model_name in enumerate(FOCUS_MODELS):
    ax = axes4[idx]
    signed_err, dev_ss_s, dev_ms_s = l3_data[model_name]

    over = signed_err < 0
    under = signed_err > 0

    # Absolute deviations grouped by direction
    data_to_plot = []
    labels_bp = []
    colors_bp = []

    for label, mask, col in [('Over-predict\n(single)', over, '#2196F3'),
                              ('Over-predict\n(multi)', over, '#64B5F6'),
                              ('Under-predict\n(single)', under, '#F44336'),
                              ('Under-predict\n(multi)', under, '#EF9A9A')]:
        if 'single' in label:
            d = np.abs(dev_ss_s[mask])
        else:
            d = np.abs(dev_ms_s[mask])
        data_to_plot.append(d)
        labels_bp.append(label)
        colors_bp.append(col)

    bp = ax.boxplot(data_to_plot, tick_labels=labels_bp, patch_artist=True,
                    showfliers=False, widths=0.6)
    for patch, color in zip(bp['boxes'], colors_bp):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax.set_ylabel('|Outcome Deviation| (mg/dL)')
    ax.set_title(MODEL_SHORT[model_name])
    ax.grid(True, alpha=0.2, axis='y')

    # Add median values as text
    for i, d in enumerate(data_to_plot):
        med = np.median(d)
        ax.text(i + 1, med + 0.5, f'{med:.1f}', ha='center', fontsize=8, fontweight='bold')

plt.tight_layout()
fig4.savefig(OUT_DIR / 'coupling_level3_direction.png', dpi=150, bbox_inches='tight')
print(f"  Saved: coupling_level3_direction.png")


# ── Chart 5: Level 4 Threshold ───────────────────────────────────────────────

fig5, axes5 = plt.subplots(1, len(FOCUS_MODELS), figsize=(7 * len(FOCUS_MODELS), 5.5))
fig5.suptitle("Level 4: Outcome Deviation by Prediction Error Magnitude\n"
              "(How much does ISF accuracy matter at each error level?)",
              fontsize=12, fontweight='bold')

if len(FOCUS_MODELS) == 1:
    axes5 = [axes5]

for idx, model_name in enumerate(FOCUS_MODELS):
    ax = axes5[idx]
    bin_stats = l4_data[model_name]

    x_pos = np.arange(len(BIN_LABELS))
    mean_ss = [b[2] for b in bin_stats]
    med_ss = [b[3] for b in bin_stats]
    mean_ms = [b[4] for b in bin_stats]
    med_ms = [b[5] for b in bin_stats]

    ax.plot(x_pos, mean_ss, 'o-', color='#F44336', linewidth=2, markersize=8,
            label='Mean |dev| single-step')
    ax.plot(x_pos, mean_ms, 's-', color='#2196F3', linewidth=2, markersize=8,
            label='Mean |dev| multi-step')
    ax.plot(x_pos, med_ss, 'o--', color='#F44336', linewidth=1, markersize=6,
            alpha=0.6, label='Median |dev| single-step')
    ax.plot(x_pos, med_ms, 's--', color='#2196F3', linewidth=1, markersize=6,
            alpha=0.6, label='Median |dev| multi-step')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(BIN_LABELS)
    ax.set_xlabel('|Prediction Error| bin (mg/dL)')
    ax.set_ylabel('|Outcome Deviation| (mg/dL)')
    ax.set_title(MODEL_SHORT[model_name])
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    # Annotate sample counts
    for i, b in enumerate(bin_stats):
        ax.text(i, max(mean_ss[i] or 0, mean_ms[i] or 0) + 1,
                f'n={b[1]:,}', ha='center', fontsize=7, alpha=0.6)

plt.tight_layout()
fig5.savefig(OUT_DIR / 'coupling_level4_threshold.png', dpi=150, bbox_inches='tight')
print(f"  Saved: coupling_level4_threshold.png")


# ── Chart 6: Dampening Distribution ──────────────────────────────────────────

fig6, ax6 = plt.subplots(figsize=(10, 5.5))
fig6.suptitle("Feedback Dampening: Multi-step / Single-step Delta Insulin Ratio\n"
              "(Values < 1 = feedback attenuates the dosing change)",
              fontsize=12, fontweight='bold')

violin_data = []
violin_labels = []
violin_colors = []
for model_name in dampening_data:
    ratios = dampening_data[model_name]
    # Clip extreme ratios for visualisation
    clipped = np.clip(ratios, -2, 3)
    violin_data.append(clipped)
    violin_labels.append(MODEL_SHORT[model_name])
    violin_colors.append(MODEL_COLORS[model_name])

if violin_data:
    positions = range(len(violin_data))
    parts = ax6.violinplot(violin_data, positions=positions, showmedians=True,
                           showextrema=False)
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(violin_colors[i])
        pc.set_alpha(0.6)
    parts['cmedians'].set_color('black')

    # Add IQR
    for i, d in enumerate(violin_data):
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax6.vlines(i, q1, q3, color='black', linewidth=3, alpha=0.5)

    ax6.axhline(1.0, color='red', linestyle='--', alpha=0.4, linewidth=1,
                label='No dampening (ratio = 1)')
    ax6.axhline(0.0, color='grey', linestyle=':', alpha=0.3, linewidth=1)
    ax6.set_xticks(positions)
    ax6.set_xticklabels(violin_labels, rotation=30, ha='right', fontsize=9)
    ax6.set_ylabel('Delta Insulin Ratio (multi / single)')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.2, axis='y')

plt.tight_layout()
fig6.savefig(OUT_DIR / 'coupling_dampening.png', dpi=150, bbox_inches='tight')
print(f"  Saved: coupling_dampening.png")


# ── Chart 7: Boost Sequential Validation ─────────────────────────────────────

fig7, axes7 = plt.subplots(1, len(boost_seq_results) if boost_seq_results else 1,
                            figsize=(7 * max(len(boost_seq_results), 1), 6))
fig7.suptitle("Boost Sequential Validation (User-M)\n"
              "Simulated end of sample i vs actual start of sample i+1 (4-6h gap)",
              fontsize=12, fontweight='bold')

if not isinstance(axes7, np.ndarray):
    axes7 = [axes7]

if boost_seq_results:
    for idx, model_name in enumerate(boost_seq_results):
        ax = axes7[idx]
        bsr = boost_seq_results[model_name]
        actual = bsr['actual']

        # Single-step
        ax.scatter(bsr['sim_ss'], actual, s=15, alpha=0.4, c='#F44336',
                   label=f'Single-step (r={bsr["r_ss"]:.3f}, MAE={bsr["mae_ss"]:.1f})',
                   edgecolors='none')
        # Multi-step
        ax.scatter(bsr['sim_ms'], actual, s=15, alpha=0.4, c='#2196F3',
                   label=f'Multi-step (r={bsr["r_ms"]:.3f}, MAE={bsr["mae_ms"]:.1f})',
                   edgecolors='none')

        # Identity line
        lo = min(actual.min(), bsr['sim_ss'].min(), bsr['sim_ms'].min()) - 10
        hi = max(actual.max(), bsr['sim_ss'].max(), bsr['sim_ms'].max()) + 10
        ax.plot([lo, hi], [lo, hi], 'k--', alpha=0.3, linewidth=1)

        ax.set_xlabel('Simulated End BG (mg/dL)')
        ax.set_ylabel('Next Sample Start BG (mg/dL)')
        ax.set_title(MODEL_SHORT[model_name])
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)
else:
    axes7[0].text(0.5, 0.5, 'No sequential validation data available',
                  ha='center', va='center', transform=axes7[0].transAxes)

plt.tight_layout()
fig7.savefig(OUT_DIR / 'coupling_boost_validation.png', dpi=150, bbox_inches='tight')
print(f"  Saved: coupling_boost_validation.png")


# ══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

print("\n\n" + "=" * 90)
print("SUMMARY")
print("=" * 90)

print("\nKey findings:")
print("  1. Level 1 (Model-Site): Correlations between prediction MAE and outcome metrics")
print("  2. Level 2 (Per-Sample): How individual prediction errors propagate to outcomes")
print("  3. Level 3 (Direction): Asymmetry in over- vs under-prediction impact")
print("  4. Level 4 (Threshold): The error magnitude at which outcomes start diverging")
print("  5. Dampening: How much feedback attenuates ISF model differences")

# Quick summary stats
ad = all_results['allday']
mae = ad['l1_mae']
for key, label in [('l1_tir_ms', 'TIR'), ('l1_titr_ms', 'TITR')]:
    y = ad[key]
    m_ok = np.isfinite(mae) & np.isfinite(y)
    if m_ok.sum() > 5:
        r, p = stats.pearsonr(mae[m_ok], y[m_ok])
        print(f"\n  All-day MAE vs {label} (multi-step): r = {r:.3f}, p = {p:.4f}")

if dampening_data:
    all_ratios = np.concatenate(list(dampening_data.values()))
    print(f"\n  Overall dampening ratio: median = {np.median(all_ratios):.3f}, "
          f"mean = {np.mean(all_ratios):.3f}")
    print(f"  Fraction with ratio < 0.9 (meaningful dampening): "
          f"{np.mean(all_ratios < 0.9)*100:.1f}%")

# ── Comprehensive all-model summary table ─────────────────────────────────────

print("\n\n" + "=" * 120)
print("ALL-MODEL SUMMARY TABLE (for Paper 9)")
print("=" * 120)
print(f"\n  {'Model':25s} {'r(SS)':>7s} {'r(MS)':>7s} "
      f"{'Over med':>10s} {'Under med':>10s} {'Asym':>6s} "
      f"{'40+ mean':>10s} {'40+ med':>10s}")
print("  " + "-" * 110)

for model_name in ALL_MODELS:
    # Level 2: per-sample correlation
    pe, dss, dms = l2_data[model_name]
    mask = np.isfinite(pe) & np.isfinite(dss) & np.isfinite(dms)
    r_ss, _ = stats.pearsonr(pe[mask], dss[mask]) if mask.sum() > 5 else (np.nan, 1)
    r_ms, _ = stats.pearsonr(pe[mask], dms[mask]) if mask.sum() > 5 else (np.nan, 1)

    # Level 3: direction asymmetry
    se, _, dms_s = l3_data[model_name]
    over = se < 0
    under = se > 0
    over_med = np.median(np.abs(dms_s[over])) if over.sum() > 0 else np.nan
    under_med = np.median(np.abs(dms_s[under])) if under.sum() > 0 else np.nan
    asym = over_med / under_med if under_med > 0 else np.nan

    # Level 4: 40+ bin
    bin_40 = [b for b in l4_data[model_name] if b[0] == '40+']
    mean_40 = bin_40[0][4] if bin_40 and not np.isnan(bin_40[0][4]) else np.nan  # mean MS
    med_40 = bin_40[0][5] if bin_40 and not np.isnan(bin_40[0][5]) else np.nan   # med MS

    print(f"  {model_name:25s} {r_ss:7.3f} {r_ms:7.3f} "
          f"{over_med:10.2f} {under_med:10.2f} {asym:6.2f} "
          f"{mean_40:10.2f} {med_40:10.2f}")

print("\n\nDONE — all charts saved to:", OUT_DIR)
