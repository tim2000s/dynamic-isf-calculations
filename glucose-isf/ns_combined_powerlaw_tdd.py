#!/usr/bin/env python3
"""
Combined Power-Law BG Scaling + TDD_effective Analysis
=======================================================
Tests whether power-law scaling and overnight TDD calibration are additive.

Variants:
  A  Current loop ISF  (isf_v1 = variable_sens)
  B  Power-law + static TDD_7day
  C  ln-scaling + TDD_effective (overnight regression)
  D  Power-law + TDD_effective
"""

import warnings, textwrap, pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

warnings.filterwarnings("ignore")

# ── Theme ──────────────────────────────────────────────────────────
BG_C   = "#0f0f0f"
PANEL  = "#1a1a2e"
GRID   = "#2a2a4a"
TXT    = "#e0e0ff"
COLS   = {"A": "#ff6b6b", "B": "#4ecdc4", "C": "#f7b731", "D": "#a29bfe"}

OUT_FIG = pathlib.Path("/Users/tims/Downloads/ns_combined_powerlaw_tdd.png")
OUT_TXT = pathlib.Path("/Users/tims/Downloads/ns_combined_powerlaw_tdd_summary.txt")
CSV     = pathlib.Path("/Users/tims/Downloads/ns_backtest_overnight.csv")

TARGET = 99.0
D_CONST = 82.0
K_FIXED = 2.0
REF_FACTOR = 1700.0
LN_REF = np.log(TARGET / D_CONST + 1)   # ln(99/82+1) ≈ ln(2.207)
ALPHA = 0.15   # EMA smoothing for TDD_effective

# ── Helpers ────────────────────────────────────────────────────────
def ln_scale(bg):
    return LN_REF / np.log(bg / D_CONST + 1)

def power_scale(bg, k):
    return (TARGET / bg) ** k

def isf_ln_tdd(tdd, bg):
    return (REF_FACTOR / tdd) * ln_scale(bg)

def isf_pow_tdd(tdd, bg, k):
    return (REF_FACTOR / tdd) * power_scale(bg, k)


# ── 1. Load & filter ──────────────────────────────────────────────
raw = pd.read_csv(CSV)
raw["ts"] = pd.to_datetime(raw["ts"], format="ISO8601", utc=True)
raw["date"] = pd.to_datetime(raw["date"])

df = raw.copy()
df = df[df["actual_bg_2h"].notna()]
df = df[(df["bg"] >= 72) & (df["bg"] <= 200)]
df = df[df["tdd_7day"].notna()]

# predicted BG drop using loop ISF
df["bg_drop_pred"] = df["bg"] - df["pred_iob_24"]
df = df[df["bg_drop_pred"].abs() >= 3]

# ratio filter
df["ratio"] = (df["bg"] - df["actual_bg_2h"]) / (df["bg"] - df["pred_iob_24"])
df = df[(df["ratio"] >= 0) & (df["ratio"] <= 5)]

print(f"Filtered samples: {len(df)}")

# ── 2. Compute TDD_effective ──────────────────────────────────────
# Per-sample implied effective TDD from overnight regression ratio
df["isf_eff"] = df["isf_v1"] * df["ratio"]
df["tdd_eff_sample"] = REF_FACTOR * LN_REF / (df["isf_eff"] * np.log(df["bg"] / D_CONST + 1))
df["tdd_eff_sample"] = df["tdd_eff_sample"].where(
    (df["tdd_eff_sample"] >= 3) & (df["tdd_eff_sample"] <= 120)
)

# Nightly median (overnight hours already filtered in data: 0-7)
nightly = (
    df.dropna(subset=["tdd_eff_sample"])
    .groupby("date")["tdd_eff_sample"]
    .agg(["median", "count"])
    .rename(columns={"median": "tdd_implied_night", "count": "n_samples"})
)
nightly = nightly[nightly["n_samples"] >= 5]
nightly = nightly.sort_index()

# EMA with 1-day forward shift (no look-ahead)
tdd_eff_series = {}
prev = nightly["tdd_implied_night"].iloc[0]  # seed with first night
for dt, row in nightly.iterrows():
    prev = ALPHA * row["tdd_implied_night"] + (1 - ALPHA) * prev
    tdd_eff_series[dt] = prev

tdd_eff_df = pd.DataFrame.from_dict(tdd_eff_series, orient="index", columns=["tdd_effective"])
tdd_eff_df.index = pd.to_datetime(tdd_eff_df.index)
tdd_eff_df = tdd_eff_df.sort_index()
# Shift forward 1 day — value computed from night N is available on day N+1
tdd_eff_df["tdd_effective_shifted"] = tdd_eff_df["tdd_effective"].shift(1)
tdd_eff_df = tdd_eff_df.dropna(subset=["tdd_effective_shifted"])

df = df.merge(tdd_eff_df[["tdd_effective_shifted"]], left_on="date", right_index=True, how="inner")
df.rename(columns={"tdd_effective_shifted": "tdd_eff"}, inplace=True)

print(f"Samples with TDD_effective: {len(df)}")
print(f"TDD_effective range: {df['tdd_eff'].min():.1f} – {df['tdd_eff'].max():.1f}, "
      f"median {df['tdd_eff'].median():.1f}")

# ── 3. Compute ISF variants & counterfactual predictions ──────────
df["isf_A"] = df["isf_v1"]
df["isf_B"] = isf_pow_tdd(df["tdd_7day"], df["bg"], K_FIXED)
df["isf_C"] = isf_ln_tdd(df["tdd_eff"], df["bg"])
df["isf_D"] = isf_pow_tdd(df["tdd_eff"], df["bg"], K_FIXED)

for v in "ABCD":
    col = f"isf_{v}"
    df[f"pred_{v}"] = df["bg"] - df["bg_drop_pred"] * (df[col] / df["isf_v1"])
    df[f"err_{v}"] = df["actual_bg_2h"] - df[f"pred_{v}"]

# ── 4. Aggregate metrics ──────────────────────────────────────────
def metrics(errs):
    return {"MAE": np.abs(errs).mean(), "Bias": errs.mean(),
            "RMSE": np.sqrt((errs**2).mean()), "N": len(errs)}

summary_rows = []
for v in "ABCD":
    m = metrics(df[f"err_{v}"])
    m["Variant"] = v
    summary_rows.append(m)
summary = pd.DataFrame(summary_rows).set_index("Variant")

# BG bands
bands = [(0, 90, "<90"), (90, 105, "90-105"), (105, 120, "105-120"), (120, 150, "120-150")]
band_stats = []
for lo, hi, label in bands:
    sub = df[(df["bg"] >= lo) & (df["bg"] < hi)]
    if len(sub) < 10:
        continue
    for v in "ABCD":
        m = metrics(sub[f"err_{v}"])
        m["Variant"] = v
        m["Band"] = label
        band_stats.append(m)
band_df = pd.DataFrame(band_stats)

# ── 5. k grid search (with and without TDD_effective) ─────────────
ks = np.arange(1.0, 3.01, 0.1)
grid_results = []
for k in ks:
    # B-like: power-law + TDD_7day
    isf_b = isf_pow_tdd(df["tdd_7day"], df["bg"], k)
    pred_b = df["bg"] - df["bg_drop_pred"] * (isf_b / df["isf_v1"])
    err_b = df["actual_bg_2h"] - pred_b

    # D-like: power-law + TDD_eff
    isf_d = isf_pow_tdd(df["tdd_eff"], df["bg"], k)
    pred_d = df["bg"] - df["bg_drop_pred"] * (isf_d / df["isf_v1"])
    err_d = df["actual_bg_2h"] - pred_d

    grid_results.append({
        "k": k,
        "MAE_static": np.abs(err_b).mean(),
        "Bias_static": err_b.mean(),
        "MAE_tddeff": np.abs(err_d).mean(),
        "Bias_tddeff": err_d.mean(),
    })
grid = pd.DataFrame(grid_results)

best_static = grid.loc[grid["MAE_static"].idxmin()]
best_tddeff = grid.loc[grid["MAE_tddeff"].idxmin()]

# ── 6. Per-night MAE timeline (B vs D) ───────────────────────────
night_mae = []
for dt, grp in df.groupby("date"):
    if len(grp) < 3:
        continue
    row = {"date": dt}
    for v in ["B", "D"]:
        row[f"MAE_{v}"] = np.abs(grp[f"err_{v}"]).mean()
    night_mae.append(row)
night_df = pd.DataFrame(night_mae)

# ── 7. Summary text ──────────────────────────────────────────────
lines = []
lines.append("=" * 70)
lines.append("  Combined Power-Law + TDD_effective Analysis")
lines.append("=" * 70)
lines.append("")
lines.append("Variants:")
lines.append("  A  Current loop ISF (isf_v1 = variable_sens)")
lines.append(f"  B  Power-law k={K_FIXED:.1f} + static TDD_7day")
lines.append(f"  C  ln-scaling + TDD_effective (alpha={ALPHA})")
lines.append(f"  D  Power-law k={K_FIXED:.1f} + TDD_effective")
lines.append("")
lines.append(f"Filtered samples: {len(df)}")
lines.append(f"TDD_7day  median: {df['tdd_7day'].median():.2f}")
lines.append(f"TDD_eff   median: {df['tdd_eff'].median():.2f}")
lines.append("")
lines.append("─── Overall Metrics ───")
lines.append(summary.to_string())
lines.append("")

# Improvement vs A
for v in "BCD":
    mae_imp = (summary.loc["A", "MAE"] - summary.loc[v, "MAE"]) / summary.loc["A", "MAE"] * 100
    bias_imp = abs(summary.loc["A", "Bias"]) - abs(summary.loc[v, "Bias"])
    lines.append(f"  {v} vs A:  MAE {mae_imp:+.1f}%  |  Bias improvement {bias_imp:+.1f} mg/dL")

# Additivity check
mae_A = summary.loc["A", "MAE"]
mae_B = summary.loc["B", "MAE"]
mae_C = summary.loc["C", "MAE"]
mae_D = summary.loc["D", "MAE"]
imp_B = mae_A - mae_B
imp_C = mae_A - mae_C
imp_D = mae_A - mae_D
lines.append("")
lines.append("─── Additivity Check ───")
lines.append(f"  Individual improvement B (power-law):     {imp_B:.2f} mg/dL")
lines.append(f"  Individual improvement C (TDD_eff):       {imp_C:.2f} mg/dL")
lines.append(f"  Sum of individual improvements (B+C):     {imp_B + imp_C:.2f} mg/dL")
lines.append(f"  Actual combined improvement D:            {imp_D:.2f} mg/dL")
additivity = imp_D / (imp_B + imp_C) * 100 if (imp_B + imp_C) != 0 else float("nan")
lines.append(f"  Additivity ratio: {additivity:.0f}%  (100% = perfectly additive)")

lines.append("")
lines.append("─── k Grid Search ───")
lines.append(f"  Best k (static TDD_7day): {best_static['k']:.1f}  MAE={best_static['MAE_static']:.2f}  Bias={best_static['Bias_static']:+.2f}")
lines.append(f"  Best k (TDD_effective):   {best_tddeff['k']:.1f}  MAE={best_tddeff['MAE_tddeff']:.2f}  Bias={best_tddeff['Bias_tddeff']:+.2f}")

lines.append("")
lines.append("─── Per-BG-Band Metrics ───")
for label in [b[2] for b in bands]:
    sub = band_df[band_df["Band"] == label]
    if sub.empty:
        continue
    lines.append(f"\n  Band {label} mg/dL:")
    for _, r in sub.iterrows():
        lines.append(f"    {r['Variant']}  MAE={r['MAE']:.1f}  Bias={r['Bias']:+.1f}  N={int(r['N'])}")

lines.append("")
lines.append("=" * 70)

summary_text = "\n".join(lines)
print(summary_text)
OUT_TXT.write_text(summary_text)
print(f"\nSummary saved to {OUT_TXT}")

# ── 8. Figure ─────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor": BG_C, "axes.facecolor": PANEL,
    "axes.edgecolor": GRID, "axes.labelcolor": TXT,
    "text.color": TXT, "xtick.color": TXT, "ytick.color": TXT,
    "grid.color": GRID, "grid.alpha": 0.5,
    "legend.facecolor": PANEL, "legend.edgecolor": GRID,
    "font.size": 10,
})

fig = plt.figure(figsize=(22, 26))
gs = GridSpec(4, 2, figure=fig, hspace=0.35, wspace=0.28,
              left=0.06, right=0.96, top=0.96, bottom=0.03)

# Panel 1: ISF curves (full width)
ax1 = fig.add_subplot(gs[0, :])
bg_range = np.linspace(72, 200, 200)
tdd_med = df["tdd_7day"].median()
tdd_eff_med = df["tdd_eff"].median()
ax1.plot(bg_range, isf_ln_tdd(tdd_med, bg_range), color=COLS["A"], lw=2.2, label=f"A  ln-scale, TDD_7d={tdd_med:.1f}")
ax1.plot(bg_range, isf_pow_tdd(tdd_med, bg_range, K_FIXED), color=COLS["B"], lw=2.2, label=f"B  power k={K_FIXED}, TDD_7d={tdd_med:.1f}")
ax1.plot(bg_range, isf_ln_tdd(tdd_eff_med, bg_range), color=COLS["C"], lw=2.2, ls="--", label=f"C  ln-scale, TDD_eff={tdd_eff_med:.1f}")
ax1.plot(bg_range, isf_pow_tdd(tdd_eff_med, bg_range, K_FIXED), color=COLS["D"], lw=2.2, ls="--", label=f"D  power k={K_FIXED}, TDD_eff={tdd_eff_med:.1f}")
ax1.set_xlabel("BG (mg/dL)")
ax1.set_ylabel("ISF (mg/dL per U)")
ax1.set_title("ISF Curves — All 4 Variants at Median TDD", fontsize=13, fontweight="bold")
ax1.legend(fontsize=9)
ax1.grid(True)

# Panel 2: MAE bars
ax2 = fig.add_subplot(gs[1, 0])
variants = list("ABCD")
mae_vals = [summary.loc[v, "MAE"] for v in variants]
bars = ax2.bar(variants, mae_vals, color=[COLS[v] for v in variants], edgecolor="white", linewidth=0.5)
for b, val in zip(bars, mae_vals):
    ax2.text(b.get_x() + b.get_width()/2, val + 0.15, f"{val:.2f}", ha="center", fontsize=10, color=TXT)
ax2.set_ylabel("MAE (mg/dL)")
ax2.set_title("Overall MAE", fontsize=12, fontweight="bold")
ax2.grid(axis="y")

# Panel 3: Bias bars
ax3 = fig.add_subplot(gs[1, 1])
bias_vals = [summary.loc[v, "Bias"] for v in variants]
bars = ax3.bar(variants, bias_vals, color=[COLS[v] for v in variants], edgecolor="white", linewidth=0.5)
for b, val in zip(bars, bias_vals):
    yoff = 0.15 if val >= 0 else -0.6
    ax3.text(b.get_x() + b.get_width()/2, val + yoff, f"{val:+.2f}", ha="center", fontsize=10, color=TXT)
ax3.axhline(0, color=TXT, lw=0.8, ls="--")
ax3.set_ylabel("Bias (mg/dL)")
ax3.set_title("Overall Bias", fontsize=12, fontweight="bold")
ax3.grid(axis="y")

# Panel 4: Per-band MAE grouped bars
ax4 = fig.add_subplot(gs[2, 0])
band_labels = [b[2] for b in bands if b[2] in band_df["Band"].values]
x = np.arange(len(band_labels))
w = 0.18
for i, v in enumerate(variants):
    vals = [band_df[(band_df["Band"] == bl) & (band_df["Variant"] == v)]["MAE"].values[0]
            for bl in band_labels]
    ax4.bar(x + i * w - 1.5 * w, vals, w, color=COLS[v], label=v, edgecolor="white", linewidth=0.4)
ax4.set_xticks(x)
ax4.set_xticklabels(band_labels)
ax4.set_xlabel("BG Band (mg/dL)")
ax4.set_ylabel("MAE (mg/dL)")
ax4.set_title("Per-BG-Band MAE", fontsize=12, fontweight="bold")
ax4.legend(fontsize=9)
ax4.grid(axis="y")

# Panel 5: Per-band Bias grouped bars
ax5 = fig.add_subplot(gs[2, 1])
for i, v in enumerate(variants):
    vals = [band_df[(band_df["Band"] == bl) & (band_df["Variant"] == v)]["Bias"].values[0]
            for bl in band_labels]
    ax5.bar(x + i * w - 1.5 * w, vals, w, color=COLS[v], label=v, edgecolor="white", linewidth=0.4)
ax5.axhline(0, color=TXT, lw=0.8, ls="--")
ax5.set_xticks(x)
ax5.set_xticklabels(band_labels)
ax5.set_xlabel("BG Band (mg/dL)")
ax5.set_ylabel("Bias (mg/dL)")
ax5.set_title("Per-BG-Band Bias", fontsize=12, fontweight="bold")
ax5.legend(fontsize=9)
ax5.grid(axis="y")

# Panel 6: k optimisation
ax6 = fig.add_subplot(gs[3, 0])
ax6.plot(grid["k"], grid["MAE_static"], "o-", color=COLS["B"], markersize=4, lw=1.8, label="Static TDD_7day")
ax6.plot(grid["k"], grid["MAE_tddeff"], "s-", color=COLS["D"], markersize=4, lw=1.8, label="TDD_effective")
ax6.axvline(best_static["k"], color=COLS["B"], ls=":", lw=1, alpha=0.7)
ax6.axvline(best_tddeff["k"], color=COLS["D"], ls=":", lw=1, alpha=0.7)
ax6.set_xlabel("Power-law exponent k")
ax6.set_ylabel("MAE (mg/dL)")
ax6.set_title("k Optimisation — MAE vs k", fontsize=12, fontweight="bold")
ax6.legend(fontsize=9)
ax6.grid(True)

# Panel 7: Error distributions
ax7 = fig.add_subplot(gs[3, 1])
bins = np.linspace(-60, 60, 61)
for v in variants:
    ax7.hist(df[f"err_{v}"], bins=bins, color=COLS[v], alpha=0.45, label=v, edgecolor="none")
ax7.axvline(0, color=TXT, lw=1, ls="--")
ax7.set_xlabel("Prediction Error (mg/dL)")
ax7.set_ylabel("Count")
ax7.set_title("Error Distributions", fontsize=12, fontweight="bold")
ax7.legend(fontsize=9)
ax7.grid(axis="y")

# ── Use remaining space for per-night timeline as annotation ──
# We'll add a small inset-style panel
ax8 = fig.add_axes([0.06, -0.08, 0.90, 0.10])  # below main grid
ax8.set_facecolor(PANEL)
if len(night_df) > 0:
    ax8.plot(night_df["date"], night_df["MAE_B"], color=COLS["B"], lw=1.3, alpha=0.8, label="B (power+TDD7d)")
    ax8.plot(night_df["date"], night_df["MAE_D"], color=COLS["D"], lw=1.3, alpha=0.8, label="D (power+TDDeff)")
    ax8.set_ylabel("MAE", fontsize=9)
    ax8.set_title("Per-Night MAE Timeline — B vs D", fontsize=10, fontweight="bold")
    ax8.legend(fontsize=8)
    ax8.grid(True, color=GRID, alpha=0.5)
    ax8.tick_params(labelsize=8)

# Better approach: put timeline in gs by expanding layout
# Actually let's just place it properly — redo as 5 rows
plt.close(fig)

fig = plt.figure(figsize=(22, 30))
gs = GridSpec(5, 2, figure=fig, hspace=0.35, wspace=0.28,
              left=0.06, right=0.96, top=0.97, bottom=0.03,
              height_ratios=[1, 1, 1, 1, 0.8])

# Re-draw all panels
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(bg_range, isf_ln_tdd(tdd_med, bg_range), color=COLS["A"], lw=2.2, label=f"A  ln-scale, TDD_7d={tdd_med:.1f}")
ax1.plot(bg_range, isf_pow_tdd(tdd_med, bg_range, K_FIXED), color=COLS["B"], lw=2.2, label=f"B  power k={K_FIXED}, TDD_7d={tdd_med:.1f}")
ax1.plot(bg_range, isf_ln_tdd(tdd_eff_med, bg_range), color=COLS["C"], lw=2.2, ls="--", label=f"C  ln-scale, TDD_eff={tdd_eff_med:.1f}")
ax1.plot(bg_range, isf_pow_tdd(tdd_eff_med, bg_range, K_FIXED), color=COLS["D"], lw=2.2, ls="--", label=f"D  power k={K_FIXED}, TDD_eff={tdd_eff_med:.1f}")
ax1.set_xlabel("BG (mg/dL)"); ax1.set_ylabel("ISF (mg/dL per U)")
ax1.set_title("ISF Curves — All 4 Variants at Median TDD", fontsize=13, fontweight="bold")
ax1.legend(fontsize=9); ax1.grid(True)

ax2 = fig.add_subplot(gs[1, 0])
bars = ax2.bar(variants, mae_vals, color=[COLS[v] for v in variants], edgecolor="white", linewidth=0.5)
for b, val in zip(bars, mae_vals):
    ax2.text(b.get_x()+b.get_width()/2, val+0.15, f"{val:.2f}", ha="center", fontsize=10, color=TXT)
ax2.set_ylabel("MAE (mg/dL)"); ax2.set_title("Overall MAE", fontsize=12, fontweight="bold"); ax2.grid(axis="y")

ax3 = fig.add_subplot(gs[1, 1])
bars = ax3.bar(variants, bias_vals, color=[COLS[v] for v in variants], edgecolor="white", linewidth=0.5)
for b, val in zip(bars, bias_vals):
    yoff = 0.15 if val >= 0 else -0.6
    ax3.text(b.get_x()+b.get_width()/2, val+yoff, f"{val:+.2f}", ha="center", fontsize=10, color=TXT)
ax3.axhline(0, color=TXT, lw=0.8, ls="--")
ax3.set_ylabel("Bias (mg/dL)"); ax3.set_title("Overall Bias", fontsize=12, fontweight="bold"); ax3.grid(axis="y")

ax4 = fig.add_subplot(gs[2, 0])
for i, v in enumerate(variants):
    vals = [band_df[(band_df["Band"]==bl)&(band_df["Variant"]==v)]["MAE"].values[0] for bl in band_labels]
    ax4.bar(x+i*w-1.5*w, vals, w, color=COLS[v], label=v, edgecolor="white", linewidth=0.4)
ax4.set_xticks(x); ax4.set_xticklabels(band_labels)
ax4.set_xlabel("BG Band (mg/dL)"); ax4.set_ylabel("MAE (mg/dL)")
ax4.set_title("Per-BG-Band MAE", fontsize=12, fontweight="bold"); ax4.legend(fontsize=9); ax4.grid(axis="y")

ax5 = fig.add_subplot(gs[2, 1])
for i, v in enumerate(variants):
    vals = [band_df[(band_df["Band"]==bl)&(band_df["Variant"]==v)]["Bias"].values[0] for bl in band_labels]
    ax5.bar(x+i*w-1.5*w, vals, w, color=COLS[v], label=v, edgecolor="white", linewidth=0.4)
ax5.axhline(0, color=TXT, lw=0.8, ls="--"); ax5.set_xticks(x); ax5.set_xticklabels(band_labels)
ax5.set_xlabel("BG Band (mg/dL)"); ax5.set_ylabel("Bias (mg/dL)")
ax5.set_title("Per-BG-Band Bias", fontsize=12, fontweight="bold"); ax5.legend(fontsize=9); ax5.grid(axis="y")

ax6 = fig.add_subplot(gs[3, 0])
ax6.plot(grid["k"], grid["MAE_static"], "o-", color=COLS["B"], markersize=4, lw=1.8, label="Static TDD_7day")
ax6.plot(grid["k"], grid["MAE_tddeff"], "s-", color=COLS["D"], markersize=4, lw=1.8, label="TDD_effective")
ax6.axvline(best_static["k"], color=COLS["B"], ls=":", lw=1, alpha=0.7)
ax6.axvline(best_tddeff["k"], color=COLS["D"], ls=":", lw=1, alpha=0.7)
ax6.annotate(f"k*={best_static['k']:.1f}", xy=(best_static["k"], best_static["MAE_static"]),
             xytext=(10, 10), textcoords="offset points", color=COLS["B"], fontsize=9)
ax6.annotate(f"k*={best_tddeff['k']:.1f}", xy=(best_tddeff["k"], best_tddeff["MAE_tddeff"]),
             xytext=(10, -15), textcoords="offset points", color=COLS["D"], fontsize=9)
ax6.set_xlabel("Power-law exponent k"); ax6.set_ylabel("MAE (mg/dL)")
ax6.set_title("k Optimisation — MAE vs k", fontsize=12, fontweight="bold"); ax6.legend(fontsize=9); ax6.grid(True)

ax7 = fig.add_subplot(gs[3, 1])
bins = np.linspace(-60, 60, 61)
for v in variants:
    ax7.hist(df[f"err_{v}"], bins=bins, color=COLS[v], alpha=0.45, label=v, edgecolor="none")
ax7.axvline(0, color=TXT, lw=1, ls="--")
ax7.set_xlabel("Prediction Error (mg/dL)"); ax7.set_ylabel("Count")
ax7.set_title("Error Distributions", fontsize=12, fontweight="bold"); ax7.legend(fontsize=9); ax7.grid(axis="y")

ax8 = fig.add_subplot(gs[4, :])
if len(night_df) > 0:
    ax8.plot(night_df["date"], night_df["MAE_B"], color=COLS["B"], lw=1.3, alpha=0.85, label="B (power + TDD_7day)")
    ax8.plot(night_df["date"], night_df["MAE_D"], color=COLS["D"], lw=1.3, alpha=0.85, label="D (power + TDD_eff)")
    # 7-day rolling
    if len(night_df) > 7:
        ax8.plot(night_df["date"], night_df["MAE_B"].rolling(7, min_periods=3).mean(),
                 color=COLS["B"], lw=2.5, alpha=0.6, ls="--")
        ax8.plot(night_df["date"], night_df["MAE_D"].rolling(7, min_periods=3).mean(),
                 color=COLS["D"], lw=2.5, alpha=0.6, ls="--")
ax8.set_ylabel("MAE (mg/dL)"); ax8.set_xlabel("Date")
ax8.set_title("Per-Night MAE Timeline — B vs D  (dashed = 7-day rolling)", fontsize=12, fontweight="bold")
ax8.legend(fontsize=9); ax8.grid(True)
ax8.tick_params(axis="x", rotation=30)

fig.savefig(OUT_FIG, dpi=160, bbox_inches="tight")
print(f"\nFigure saved to {OUT_FIG}")
