#!/usr/bin/env python3
"""
TDD Drift Analysis — Tests whether 7-day TDD can serve as a
relative drift correction or change-detection trigger for dynamic ISF.

Hypotheses tested:
  Test 1: Relative Drift Correction (TDD_baseline / TDD_7day scaling)
  Test 4: Change Detection Trigger (re-anchor ISF when TDD shifts >15%)

Data: Boost/AAPS cache (User-M), allday DataFrame.
"""

import pickle
import warnings
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MaxNLocator

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ────────────────────────────────────────────────────────────
BASE = "/Users/tims/Downloads/4 Hour analysis"
CACHE = f"{BASE}/daytime analysis/boost_allday_cache.pkl"

# ── Load data ────────────────────────────────────────────────────────
with open(CACHE, "rb") as fh:
    raw = pickle.load(fh)
df = raw["allday"].copy()
df = df.sort_values("ts").reset_index(drop=True)
df["ts"] = pd.to_datetime(df["ts"])

print(f"Loaded {len(df):,} samples  |  {df['ts'].min():%Y-%m-%d} → {df['ts'].max():%Y-%m-%d}")
print(f"BG range: {df['bg'].min():.0f}–{df['bg'].max():.0f} mg/dL")
print(f"TDD_7day range: {df['tdd_7day'].min():.1f}–{df['tdd_7day'].max():.1f} U")
print()

# ── Helpers ──────────────────────────────────────────────────────────

def ratio_fn(bg_values):
    """Piecewise-linear population ratio curve.
    Anchors: {76: 1.15, 100: 1.00, 130: 0.80, 170: 0.70}
    Extrapolate flat outside the knots.
    """
    knots_bg  = np.array([76.0, 100.0, 130.0, 170.0])
    knots_rat = np.array([1.15,  1.00,  0.80,  0.70])
    return np.interp(bg_values, knots_bg, knots_rat)


def counterfactual_end(bg, pred_drop, isf_model, variable_sens):
    """predicted_end = bg - pred_drop * (ISF_model / variable_sens)"""
    return bg - pred_drop * (isf_model / variable_sens)


def mae(predicted, actual):
    mask = np.isfinite(predicted) & np.isfinite(actual)
    return np.mean(np.abs(predicted[mask] - actual[mask]))


def rolling_mae(predicted, actual, window=200):
    err = np.abs(predicted - actual)
    s = pd.Series(err)
    return s.rolling(window, min_periods=50).mean()


# ══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ══════════════════════════════════════════════════════════════════════

# ISF anchor: median variable_sens where BG ≈ 100
mask_anchor = (df["bg"] >= 96) & (df["bg"] <= 104)
ISF_ANCHOR = df.loc[mask_anchor, "variable_sens"].median()
print(f"ISF_anchor (median ISF @ BG 96-104): {ISF_ANCHOR:.1f} mg/dL per U")
print(f"  (computed from {mask_anchor.sum():,} samples)")

# TDD baseline: global median of tdd_7day
TDD_BASELINE = df["tdd_7day"].median()
print(f"TDD_baseline (global median tdd_7day): {TDD_BASELINE:.2f} U")
print()

# ══════════════════════════════════════════════════════════════════════
# TEST 1 — Relative Drift Correction
# ══════════════════════════════════════════════════════════════════════

print("=" * 72)
print("TEST 1: RELATIVE DRIFT CORRECTION")
print("=" * 72)

# --- Model A: Static anchor (no TDD adjustment) ---
df["isf_static"] = ISF_ANCHOR * ratio_fn(df["bg"].values)
df["pred_end_static"] = counterfactual_end(
    df["bg"].values, df["pred_drop"].values,
    df["isf_static"].values, df["variable_sens"].values
)

# --- Model B: TDD drift-corrected (global baseline) ---
df["isf_tdd_global"] = ISF_ANCHOR * (TDD_BASELINE / df["tdd_7day"].values) * ratio_fn(df["bg"].values)
df["pred_end_tdd_global"] = counterfactual_end(
    df["bg"].values, df["pred_drop"].values,
    df["isf_tdd_global"].values, df["variable_sens"].values
)

# --- Model C: TDD drift-corrected (30-day rolling baseline) ---
df["tdd_rolling_30d_med"] = (
    df.set_index("ts")["tdd_7day"]
    .rolling("30D", min_periods=10)
    .median()
    .values
)
df["isf_tdd_rolling"] = ISF_ANCHOR * (df["tdd_rolling_30d_med"] / df["tdd_7day"]) * ratio_fn(df["bg"].values)
df["pred_end_tdd_rolling"] = counterfactual_end(
    df["bg"].values, df["pred_drop"].values,
    df["isf_tdd_rolling"].values, df["variable_sens"].values
)

# --- Model D: System's actual variable_sens (reference) ---
df["pred_end_system"] = counterfactual_end(
    df["bg"].values, df["pred_drop"].values,
    df["variable_sens"].values, df["variable_sens"].values
)
# This simplifies to bg - pred_drop, which is the system's own prediction

mae_static     = mae(df["pred_end_static"].values, df["actual_bg_end"].values)
mae_tdd_global = mae(df["pred_end_tdd_global"].values, df["actual_bg_end"].values)
mae_tdd_roll   = mae(df["pred_end_tdd_rolling"].values, df["actual_bg_end"].values)
mae_system     = mae(df["pred_end_system"].values, df["actual_bg_end"].values)

print(f"\n{'Model':<42} {'MAE (mg/dL)':>12}")
print("-" * 56)
print(f"{'System (actual variable_sens)':<42} {mae_system:>12.2f}")
print(f"{'Static anchor (no TDD)':<42} {mae_static:>12.2f}")
print(f"{'TDD drift — global baseline':<42} {mae_tdd_global:>12.2f}")
print(f"{'TDD drift — 30-day rolling baseline':<42} {mae_tdd_roll:>12.2f}")
print()

# Breakdown by TDD quartile
df["tdd_quartile"] = pd.qcut(df["tdd_7day"], 4, labels=["Q1 (low)", "Q2", "Q3", "Q4 (high)"])
print("MAE by TDD_7day quartile:")
print(f"{'Quartile':<14} {'TDD range':>16} {'Static':>10} {'TDD-global':>12} {'TDD-roll':>10} {'System':>10}")
print("-" * 74)
for q in ["Q1 (low)", "Q2", "Q3", "Q4 (high)"]:
    m = df["tdd_quartile"] == q
    tdd_lo = df.loc[m, "tdd_7day"].min()
    tdd_hi = df.loc[m, "tdd_7day"].max()
    print(f"{q:<14} {tdd_lo:>6.1f}–{tdd_hi:<6.1f}   "
          f"{mae(df.loc[m,'pred_end_static'].values, df.loc[m,'actual_bg_end'].values):>10.2f}"
          f"{mae(df.loc[m,'pred_end_tdd_global'].values, df.loc[m,'actual_bg_end'].values):>12.2f}"
          f"{mae(df.loc[m,'pred_end_tdd_rolling'].values, df.loc[m,'actual_bg_end'].values):>10.2f}"
          f"{mae(df.loc[m,'pred_end_system'].values, df.loc[m,'actual_bg_end'].values):>10.2f}")
print()

# ══════════════════════════════════════════════════════════════════════
# TEST 4 — Change Detection Trigger
# ══════════════════════════════════════════════════════════════════════

print("=" * 72)
print("TEST 4: CHANGE DETECTION TRIGGER")
print("=" * 72)

# Detect shift events: tdd_7day deviates >15% from trailing 30-day median
df["tdd_30d_med"] = (
    df.set_index("ts")["tdd_7day"]
    .rolling("30D", min_periods=10)
    .median()
    .values
)
df["tdd_pct_shift"] = (df["tdd_7day"] - df["tdd_30d_med"]) / df["tdd_30d_med"] * 100
df["shift_event"] = df["tdd_pct_shift"].abs() > 15

n_shifts = df["shift_event"].sum()
print(f"\nShift events detected (|TDD shift| > 15%): {n_shifts:,} / {len(df):,} samples "
      f"({100*n_shifts/len(df):.1f}%)")

# Compute implied ISF for qualifying samples
df["actual_drop"] = df["bg"] - df["actual_bg_end"]
df["isf_implied"] = np.nan
qual = (df["actual_drop"].abs() > 10) & (df["pred_drop"].abs() > 5)
df.loc[qual, "isf_implied"] = df.loc[qual, "variable_sens"] * (
    df.loc[qual, "actual_drop"] / df.loc[qual, "pred_drop"]
)
# Normalise implied ISF to BG=100 by dividing out the ratio curve
df.loc[qual, "isf_implied_at100"] = df.loc[qual, "isf_implied"] / ratio_fn(df.loc[qual, "bg"].values)

print(f"Qualifying samples for implied ISF (|actual_drop|>10, |pred_drop|>5): {qual.sum():,}")

# --- Model E: Triggered re-anchor ---
# Walk through chronologically; re-anchor when shift detected
anchor_e = ISF_ANCHOR
isf_triggered = np.full(len(df), np.nan)
anchors_over_time = []
last_shift_idx = -999

for i in range(len(df)):
    if df["shift_event"].iloc[i] and (i - last_shift_idx > 20):
        # Re-anchor: median of last 20 qualifying implied ISF@100
        lookback = df["isf_implied_at100"].iloc[max(0, i-200):i].dropna()
        if len(lookback) >= 5:
            recent_20 = lookback.tail(20)
            anchor_e = recent_20.median()
            last_shift_idx = i
            anchors_over_time.append((df["ts"].iloc[i], anchor_e, "re-anchor"))

    isf_triggered[i] = anchor_e * ratio_fn(np.array([df["bg"].iloc[i]]))[0]

df["isf_triggered"] = isf_triggered
df["pred_end_triggered"] = counterfactual_end(
    df["bg"].values, df["pred_drop"].values,
    df["isf_triggered"].values, df["variable_sens"].values
)

# --- Model F: Continuous auto-learning (rolling 100 implied ISF) ---
df["isf_implied_at100_filled"] = df["isf_implied_at100"].copy()
# Forward-fill for rolling, but use raw for median
rolling_anchor = df["isf_implied_at100"].rolling(200, min_periods=20).median()
# Use ISF_ANCHOR for initial period where rolling isn't available
rolling_anchor = rolling_anchor.fillna(ISF_ANCHOR)
df["isf_continuous"] = rolling_anchor * ratio_fn(df["bg"].values)
df["pred_end_continuous"] = counterfactual_end(
    df["bg"].values, df["pred_drop"].values,
    df["isf_continuous"].values, df["variable_sens"].values
)

# MAE comparison for Test 4
mae_triggered  = mae(df["pred_end_triggered"].values, df["actual_bg_end"].values)
mae_continuous = mae(df["pred_end_continuous"].values, df["actual_bg_end"].values)

print(f"\n{'Model':<42} {'MAE (mg/dL)':>12}")
print("-" * 56)
print(f"{'System (actual variable_sens)':<42} {mae_system:>12.2f}")
print(f"{'Static anchor (never updates)':<42} {mae_static:>12.2f}")
print(f"{'TDD drift — global baseline (Test 1)':<42} {mae_tdd_global:>12.2f}")
print(f"{'TDD drift — 30-day rolling (Test 1)':<42} {mae_tdd_roll:>12.2f}")
print(f"{'Triggered re-anchor (Test 4)':<42} {mae_triggered:>12.2f}")
print(f"{'Continuous auto-learning (Test 4)':<42} {mae_continuous:>12.2f}")
print()

# ══════════════════════════════════════════════════════════════════════
# GRAND COMPARISON TABLE
# ══════════════════════════════════════════════════════════════════════

print("=" * 72)
print("GRAND COMPARISON — ALL MODELS")
print("=" * 72)

models = [
    ("System (actual variable_sens)", mae_system),
    ("Static anchor (no TDD)", mae_static),
    ("TDD drift — global baseline", mae_tdd_global),
    ("TDD drift — 30-day rolling baseline", mae_tdd_roll),
    ("Triggered re-anchor (shift >15%)", mae_triggered),
    ("Continuous auto-learning (roll 200)", mae_continuous),
]
models_sorted = sorted(models, key=lambda x: x[1])

print(f"\n{'Rank':<6} {'Model':<42} {'MAE':>10} {'vs System':>10}")
print("-" * 70)
for rank, (name, m) in enumerate(models_sorted, 1):
    delta = m - mae_system
    sign = "+" if delta >= 0 else ""
    print(f"{rank:<6} {name:<42} {m:>10.2f} {sign}{delta:>9.2f}")
print()

# Additional insight: correlation between TDD change and ISF error
df["tdd_change_pct"] = (df["tdd_7day"] - TDD_BASELINE) / TDD_BASELINE * 100
df["static_error"] = df["pred_end_static"] - df["actual_bg_end"]
corr = df[["tdd_change_pct", "static_error"]].dropna().corr().iloc[0, 1]
print(f"Correlation between TDD change (%) and static-model error: {corr:.3f}")
print("  (Positive = when TDD rises, static model over-predicts drops)")
print()

# ══════════════════════════════════════════════════════════════════════
# CHART 1 — Time series: TDD, ISF@100, shift events
# ══════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True,
                         gridspec_kw={"hspace": 0.12})

# Panel 1: TDD over time
ax = axes[0]
ax.plot(df["ts"], df["tdd_7day"], color="#2196F3", alpha=0.5, lw=0.6, label="TDD 7-day")
ax.plot(df["ts"], df["tdd_30d_med"], color="#0D47A1", lw=1.8, label="TDD 30-day median")
ax.set_ylabel("TDD (U/day)")
ax.legend(loc="upper right", fontsize=9)
ax.set_title("TDD Drift Analysis — Time Series", fontsize=13, fontweight="bold")
ax.grid(True, alpha=0.25)

# Panel 2: ISF at BG=100 (rolling implied vs anchor)
ax = axes[1]
isf_at100_system = df["variable_sens"] / ratio_fn(df["bg"].values)  # back out ISF@100
roll_isf100 = isf_at100_system.rolling(100, min_periods=20).median()
ax.plot(df["ts"], roll_isf100, color="#4CAF50", lw=1.4, label="System ISF@100 (roll 100)")
ax.axhline(ISF_ANCHOR, color="#888888", ls="--", lw=1.2, label=f"Static anchor ({ISF_ANCHOR:.0f})")
# Show rolling anchor from continuous model
ax.plot(df["ts"], rolling_anchor, color="#FF9800", lw=1.2, alpha=0.8,
        label="Auto-learned anchor (roll 200)")
ax.set_ylabel("ISF @ BG 100\n(mg/dL per U)")
ax.legend(loc="upper right", fontsize=9)
ax.grid(True, alpha=0.25)

# Panel 3: Shift events
ax = axes[2]
shift_mask = df["shift_event"]
ax.fill_between(df["ts"], 0, 1, where=shift_mask,
                color="#F44336", alpha=0.35, transform=ax.get_xaxis_transform(),
                label="TDD shift >15%")
ax.plot(df["ts"], df["tdd_pct_shift"], color="#9C27B0", lw=0.7, alpha=0.7)
ax.axhline(15, color="#F44336", ls=":", lw=1)
ax.axhline(-15, color="#F44336", ls=":", lw=1)
ax.set_ylabel("TDD shift (%)")
ax.set_xlabel("Date")
ax.legend(loc="upper right", fontsize=9)
ax.grid(True, alpha=0.25)

ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
plt.xticks(rotation=30, ha="right")

plt.tight_layout()
fig.savefig(f"{BASE}/tdd_drift_timeseries.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved: tdd_drift_timeseries.png")

# ══════════════════════════════════════════════════════════════════════
# CHART 2 — Bar chart MAE comparison
# ══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 6))

names = [m[0] for m in models_sorted]
maes  = [m[1] for m in models_sorted]
colors = []
for n in names:
    if "System" in n:
        colors.append("#2196F3")
    elif "Static" in n:
        colors.append("#9E9E9E")
    elif "global" in n:
        colors.append("#FF9800")
    elif "rolling" in n.lower() and "TDD" in n:
        colors.append("#FF5722")
    elif "Triggered" in n:
        colors.append("#4CAF50")
    elif "Continuous" in n:
        colors.append("#8BC34A")
    else:
        colors.append("#607D8B")

bars = ax.barh(range(len(names)), maes, color=colors, edgecolor="white", height=0.6)
ax.set_yticks(range(len(names)))
ax.set_yticklabels(names, fontsize=10)
ax.set_xlabel("MAE (mg/dL)", fontsize=11)
ax.set_title("Prediction MAE — All Models", fontsize=13, fontweight="bold")
ax.invert_yaxis()
ax.grid(True, axis="x", alpha=0.25)

for bar, val in zip(bars, maes):
    ax.text(val + 0.3, bar.get_y() + bar.get_height() / 2,
            f"{val:.2f}", va="center", fontsize=10, fontweight="bold")

plt.tight_layout()
fig.savefig(f"{BASE}/tdd_drift_mae_comparison.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved: tdd_drift_mae_comparison.png")

# ══════════════════════════════════════════════════════════════════════
# CHART 3 — Rolling MAE over time
# ══════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(14, 6))

window = 200
pairs = [
    ("System",              df["pred_end_system"].values,      "#2196F3", 2.0),
    ("Static anchor",       df["pred_end_static"].values,      "#9E9E9E", 1.4),
    ("TDD global baseline", df["pred_end_tdd_global"].values,  "#FF9800", 1.4),
    ("TDD rolling baseline",df["pred_end_tdd_rolling"].values, "#FF5722", 1.4),
    ("Triggered re-anchor", df["pred_end_triggered"].values,   "#4CAF50", 1.4),
    ("Continuous learning",  df["pred_end_continuous"].values,  "#8BC34A", 1.4),
]

for label, pred, col, lw in pairs:
    rmae = rolling_mae(pred, df["actual_bg_end"].values, window)
    ax.plot(df["ts"], rmae, color=col, lw=lw, alpha=0.85, label=label)

ax.set_ylabel(f"Rolling {window}-sample MAE (mg/dL)", fontsize=11)
ax.set_xlabel("Date", fontsize=11)
ax.set_title(f"Rolling MAE Over Time (window = {window})", fontsize=13, fontweight="bold")
ax.legend(loc="upper right", fontsize=9, ncol=2)
ax.grid(True, alpha=0.25)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
plt.xticks(rotation=30, ha="right")

plt.tight_layout()
fig.savefig(f"{BASE}/tdd_drift_rolling_error.png", dpi=180, bbox_inches="tight")
plt.close(fig)
print("Saved: tdd_drift_rolling_error.png")

print("\n" + "=" * 72)
print("ANALYSIS COMPLETE")
print("=" * 72)
