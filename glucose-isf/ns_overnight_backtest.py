"""
Overnight DynamicISF Backtest  (00:00 – 07:00)
================================================
Compares five ISF formulas against actual overnight BG outcomes.

Formulas tested
---------------
  v1 (Actual)  : variable_sens as logged — IS the v1 formula output
  v2           : 115000 / (TDD² × ln(BG/D+1))  [TDD back-calculated from v1]
  No-TDD       : profileISF(h) × ln(target/D+1) / ln(BG/D+1)   ← proposed
  7D-TDD       : (1700/TDD_7day) × ln(target/D+1) / ln(BG/D+1)  ← proposed
  Flat         : profileISF(h)  [no BG scaling — baseline]

TDD back-calculation
--------------------
v1 formula: variable_sens = 1800 / (TDD × ln(BG/D+1))
Rearranged: TDD_implied    = 1800 / (variable_sens × ln(BG/D+1))

This gives the effective TDD the loop was using each cycle (absorbs any
adjustmentFactor already baked into variable_sens). v2 ISF is then computed
using this same TDD, showing what the quadratic formula would have produced.

7D-TDD derivation
-----------------
  K = (1700/TDD_7day) × ln(normalTarget/D+1)
  ISF = K / ln(BG/D+1)
  At target BG: ISF = 1700/TDD_7day  (recovers 1700 rule)
  TDD_7day = 7-day rolling mean of daily median tdd_implied

Counterfactual prediction
-------------------------
  bg_drop_pred = BG_t − predBGs.IOB[horizon]  (loop's predicted drop)
  pred_f       = BG_t − bg_drop_pred × (ISF_f / ISF_actual)

Outputs
-------
  ns_backtest_overnight.csv
  ns_backtest_results.png
  ns_backtest_summary.txt
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

warnings.filterwarnings('ignore')

D      = 82.0    # insulinDivisor: peak=38 min, Boost plugin → (90-38)+30 = 82
TARGET = 99.0    # normalTarget mg/dL

# ── 1. Load & filter ───────────────────────────────────────────────────────────
print("Loading pipeline data...")
df = pd.read_csv('ns_pipeline_analysis.csv')
df['ts'] = pd.to_datetime(df['ts'], format='ISO8601')
df['hour'] = df['ts'].dt.hour

on = df[df['hour'] < 8].copy()
print(f"  00–07h fasting cycles: {len(on):,}")

# ── 2. Derive per-hour profile ISF ─────────────────────────────────────────────
eq = on.query("75 <= bg <= 140 and 30 <= profile_isf_implied <= 400").copy()
hourly_profile_isf = eq.groupby('hour')['profile_isf_implied'].median()
overall_profile_isf = eq['profile_isf_implied'].median()
print(f"  Overall median profile ISF: {overall_profile_isf:.1f} mg/dL/U")

def get_profile_isf(hour):
    return hourly_profile_isf[hour] if (hour in hourly_profile_isf.index and hourly_profile_isf[hour] > 0) else overall_profile_isf

on['profile_isf_h'] = on['hour'].map(get_profile_isf)

# ── 3. Back-calculate TDD from v1, then compute all formula ISFs ───────────────
ln_target = np.log(TARGET / D + 1)

# Store ln_bg as a column so it survives the upcoming merge intact
on['ln_bg'] = np.log(on['bg'] / D + 1)

# TDD implied by v1 (absorbs any adjustmentFactor already in variable_sens)
on['tdd_implied'] = 1800.0 / (on['variable_sens'] * on['ln_bg'])

tdd_median = on['tdd_implied'].median()
tdd_mean   = on['tdd_implied'].mean()
print(f"  Implied TDD: median={tdd_median:.1f}  mean={tdd_mean:.1f}  "
      f"IQR=[{on['tdd_implied'].quantile(.25):.1f}, {on['tdd_implied'].quantile(.75):.1f}] U/day")

# 7-day rolling TDD: daily median of tdd_implied, then 7-day rolling mean
on['date'] = on['ts'].dt.date
daily_tdd = (on.groupby('date')['tdd_implied'].median()
               .reset_index().rename(columns={'tdd_implied': 'tdd_daily_med'})
               .sort_values('date'))
daily_tdd['tdd_7day'] = daily_tdd['tdd_daily_med'].rolling(7, min_periods=3).mean()
on = on.merge(daily_tdd[['date', 'tdd_7day']], on='date', how='left')
tdd_7day_median = on['tdd_7day'].median()
print(f"  7-Day rolling TDD: median={tdd_7day_median:.1f} U/day  "
      f"(days with data: {daily_tdd['tdd_7day'].notna().sum()})")

# Formula ISFs (reference on['ln_bg'] — survives merge correctly)
on['isf_v1']    = on['variable_sens']                                          # actual
on['isf_v2']    = 115000.0 / (on['tdd_implied']**2 * on['ln_bg'])             # quadratic TDD²
on['isf_notdd'] = on['profile_isf_h'] * ln_target / on['ln_bg']               # proposed no-TDD
on['isf_7dtdd'] = (1700.0 / on['tdd_7day']) * ln_target / on['ln_bg']         # proposed 7D-TDD
on['isf_flat']  = on['profile_isf_h']                                          # baseline

# Cap all formulas to same safety range to avoid division blow-ups
isf_lo, isf_hi = on['isf_v1'].quantile(0.005), on['isf_v1'].quantile(0.995)
for col in ['isf_v2', 'isf_notdd', 'isf_7dtdd', 'isf_flat']:
    on[col] = on[col].clip(isf_lo, isf_hi)

# ── 4. Counterfactual predictions ──────────────────────────────────────────────
w1 = on.dropna(subset=['pred_iob_12', 'actual_bg_1h']).copy()
w2 = on.dropna(subset=['pred_iob_24', 'actual_bg_2h']).copy()
print(f"  Clean 1h rows: {len(w1):,}   Clean 2h rows: {len(w2):,}")

formulas = [
    ('v1',    'v1 — Actual (variable_sens)', '#4fc3f7'),
    ('v2',    'v2 — Quadratic TDD²',         '#ffb74d'),
    ('notdd', 'Proposed No-TDD',             '#f48fb1'),
    ('7dtdd', 'Proposed 7D-TDD',             '#ce93d8'),
    ('flat',  'Flat profileISF',             '#a5d6a7'),
]

for df_w in [w1, w2]:
    isf_act = df_w['isf_v1'].values
    bg      = df_w['bg'].values
    for fname, _, _ in formulas:
        for horizon, pred_col, actual_col in [('1h','pred_iob_12','actual_bg_1h'),
                                               ('2h','pred_iob_24','actual_bg_2h')]:
            if pred_col not in df_w.columns or actual_col not in df_w.columns:
                continue
            isf_f  = df_w[f'isf_{fname}'].values
            pred_l = df_w[pred_col].values
            bg_drop = bg - pred_l
            pred_f  = bg - bg_drop * (isf_f / isf_act)
            df_w[f'pred_{fname}_{horizon}'] = pred_f
            df_w[f'err_{fname}_{horizon}']  = df_w[actual_col].values - pred_f

# ── 5. Dosing direction test ───────────────────────────────────────────────────
for df_w in [w1, w2]:
    actual_1h = df_w['actual_bg_1h'].values
    isf_act   = df_w['isf_v1'].values
    bg_ended_high = actual_1h > TARGET
    for fname, _, _ in formulas[1:]:   # skip v1 (it IS actual)
        isf_f = df_w[f'isf_{fname}'].values
        more_agg = isf_f < isf_act
        df_w[f'dir_correct_{fname}'] = (bg_ended_high & more_agg) | (~bg_ended_high & ~more_agg)

# ── 6. Statistical summary ─────────────────────────────────────────────────────
print("\n" + "═"*72)
print("BACKTEST RESULTS — OVERNIGHT 00:00–07:00")
print("═"*72)

summary_rows = []
for df_w, horizon, actual_col in [(w1,'1h','actual_bg_1h'), (w2,'2h','actual_bg_2h')]:
    print(f"\n── Horizon +{horizon}  (n={len(df_w):,}) " + "─"*40)
    actual = df_w[actual_col].values
    for fname, flabel, _ in formulas:
        err_col = f'err_{fname}_{horizon}'
        if err_col not in df_w.columns:
            continue
        err  = df_w[err_col].values
        mask = ~np.isnan(err)          # some formulas have NaN rows (e.g. first days of rolling TDD)
        e    = err[mask]
        n_e  = len(e)
        mae  = np.abs(e).mean()
        rmse = np.sqrt((e**2).mean())
        bias = e.mean()
        w18  = (np.abs(e) <= 18).mean() * 100
        w36  = (np.abs(e) <= 36).mean() * 100
        ppos = (e > 0).mean() * 100
        ddir = df_w[f'dir_correct_{fname}'].mean()*100 if f'dir_correct_{fname}' in df_w.columns else np.nan
        n_str = f'' if n_e == len(df_w) else f'  n={n_e:,}'
        print(f"  {flabel}")
        print(f"    bias={bias:+.1f}  MAE={mae:.1f}  RMSE={rmse:.1f}  mg/dL{n_str}")
        print(f"    within ±1mmol={w18:.1f}%  within ±2mmol={w36:.1f}%")
        print(f"    % BG > pred: {ppos:.1f}%" + (f"  dir_correct: {ddir:.1f}%" if not np.isnan(ddir) else ""))
        summary_rows.append(dict(formula=flabel, horizon=horizon,
            bias=bias, mae=mae, rmse=rmse, w18=w18, w36=w36, ppos=ppos, ddir=ddir, n=n_e))

df_sum = pd.DataFrame(summary_rows)

print("\n── ISF Magnitude (all overnight cycles) " + "─"*32)
for fname, flabel, _ in formulas:
    col = f'isf_{fname}'
    print(f"  {flabel:35s}  median={on[col].median():.1f}  mean={on[col].mean():.1f}  "
          f"5–95pct=[{on[col].quantile(.05):.0f}, {on[col].quantile(.95):.0f}]")

print(f"\n── Actual Overnight BG Outcomes " + "─"*40)
for df_w, h, ac in [(w1,'1h','actual_bg_1h'),(w2,'2h','actual_bg_2h')]:
    a = df_w[ac].values
    print(f"  +{h}: TIR(70-180)={(( a>=70)&(a<=180)).mean()*100:.1f}%  "
          f"TIR(70-140)={((a>=70)&(a<=140)).mean()*100:.1f}%  "
          f"TAR={(a>180).mean()*100:.1f}%  TBR={(a<70).mean()*100:.1f}%  "
          f"mean={a.mean():.1f} mg/dL")

# ── 7. Figure ──────────────────────────────────────────────────────────────────
print("\nGenerating figure...")

BG_C = '#0f0f0f'; PANEL = '#1a1a2e'; GRID = '#2a2a4a'; TXT = '#e0e0ff'
FCOLS = {f: c for f, _, c in formulas}

def style(ax, title):
    ax.set_facecolor(PANEL); ax.tick_params(colors=TXT, labelsize=8)
    ax.set_title(title, color=TXT, fontsize=9, fontweight='bold')
    for sp in ax.spines.values(): sp.set_edgecolor(GRID)
    ax.grid(True, color=GRID, lw=0.5, ls='--', alpha=0.7)
    ax.xaxis.label.set_color(TXT); ax.yaxis.label.set_color(TXT)
    ax.xaxis.label.set_fontsize(8); ax.yaxis.label.set_fontsize(8)

fig = plt.figure(figsize=(20, 18))
fig.patch.set_facecolor(BG_C)
gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.38)

# P1: ISF curves vs BG at median TDD (like Figure 6)
ax1 = fig.add_subplot(gs[0, :])   # full-width top panel
style(ax1, f'ISF vs BG — All Four Formulas at Median TDD = {tdd_median:.1f} U/day  (D={D:.0f} mg/dL)')
bg_range = np.linspace(70, 200, 300)
ln_bg_r  = np.log(bg_range / D + 1)
ln_tgt   = np.log(TARGET / D + 1)

curves = {
    'v1':    1800.0 / (tdd_median * ln_bg_r),
    'v2':    115000.0 / (tdd_median**2 * ln_bg_r),
    'notdd': overall_profile_isf * ln_tgt / ln_bg_r,
    '7dtdd': (1700.0 / tdd_7day_median) * ln_tgt / ln_bg_r,
    'flat':  np.full_like(bg_range, overall_profile_isf),
}

for fname, flabel, fc in formulas:
    ax1.plot(bg_range, curves[fname], lw=2.5, color=fc, label=flabel)

# Add shaded IQR bands from observed data
for fname, _, fc in formulas:
    col = f'isf_{fname}'
    q25, q75 = on[col].quantile(.25), on[col].quantile(.75)
    ax1.axhspan(q25, q75, alpha=0.07, color=fc)

ax1.axvline(TARGET, color='white', lw=0.8, ls=':', alpha=0.6, label=f'Target {TARGET:.0f} mg/dL')
ax1.axhline(overall_profile_isf, color='white', lw=0.8, ls=':', alpha=0.4)
ax1.set_xlabel('BG (mg/dL)'); ax1.set_ylabel('ISF (mg/dL per U)')
ax1.set_xlim(70, 200); ax1.set_ylim(0, min(600, curves['v2'].max() * 1.1))
ax1.legend(fontsize=9, labelcolor=TXT, facecolor=PANEL, loc='upper right')
# Annotate at BG=99 (target)
idx99 = np.argmin(np.abs(bg_range - 99))
for fname, _, fc in formulas:
    y = curves[fname][idx99]
    ax1.annotate(f'{y:.0f}', xy=(99, y), xytext=(105, y), color=fc, fontsize=8,
                 arrowprops=dict(arrowstyle='->', color=fc, lw=0.8))

# P2: Error distributions +1h
ax2 = fig.add_subplot(gs[1, 0])
style(ax2, 'Pred Error Distribution (+1h)')
for fname, flabel, fc in formulas:
    if f'err_{fname}_1h' in w1.columns:
        ax2.hist(w1[f'err_{fname}_1h'].clip(-100,100), bins=30, alpha=0.5, color=fc,
                 label=flabel.split('—')[0].strip(), density=True)
ax2.axvline(0, color='white', lw=0.8, ls='--')
ax2.set_xlabel('Pred Error (mg/dL)'); ax2.set_ylabel('Density')
ax2.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P3: MAE bar chart
ax3 = fig.add_subplot(gs[1, 1])
style(ax3, 'MAE by Formula (+1h and +2h)')
x = np.arange(len(formulas)); w = 0.35
for i, horizon in enumerate(['1h','2h']):
    maes = [df_sum.query(f"formula=='{fl}' and horizon=='{horizon}'")['mae'].values[0]
            for _, fl, _ in formulas]
    bars = ax3.bar(x + (i-0.5)*w, maes, w,
                   color=[c for _,_,c in formulas], alpha=0.85 if i==0 else 0.5,
                   label=f'+{horizon}')
    for bar, v in zip(bars, maes):
        ax3.text(bar.get_x()+bar.get_width()/2, bar.get_height()+.3,
                 f'{v:.1f}', ha='center', fontsize=6, color=TXT)
ax3.set_xticks(x); ax3.set_xticklabels(['v1','v2','No-TDD','7D-TDD','Flat'], fontsize=8)
ax3.set_ylabel('MAE (mg/dL)'); ax3.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P4: Bias bar chart
ax4 = fig.add_subplot(gs[1, 2])
style(ax4, 'Prediction Bias (+1h and +2h)')
for i, horizon in enumerate(['1h','2h']):
    biases = [df_sum.query(f"formula=='{fl}' and horizon=='{horizon}'")['bias'].values[0]
              for _, fl, _ in formulas]
    bars = ax4.bar(x + (i-0.5)*w, biases, w,
                   color=[c for _,_,c in formulas], alpha=0.85 if i==0 else 0.5,
                   label=f'+{horizon}')
    for bar, v in zip(bars, biases):
        ax4.text(bar.get_x()+bar.get_width()/2,
                 bar.get_height()+(0.5 if v>=0 else -2.5),
                 f'{v:+.1f}', ha='center', fontsize=6, color=TXT)
ax4.axhline(0, color='white', lw=0.8, ls='--')
ax4.set_xticks(x); ax4.set_xticklabels(['v1','v2','No-TDD','7D-TDD','Flat'], fontsize=8)
ax4.set_ylabel('Mean Error (mg/dL)'); ax4.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

# P5: ISF distributions
ax5 = fig.add_subplot(gs[2, 0])
style(ax5, 'ISF Distribution — Overnight Cycles')
for fname, flabel, fc in formulas:
    ax5.hist(on[f'isf_{fname}'].clip(0,500), bins=50, alpha=0.55, color=fc,
             label=flabel.split('—')[0].strip(), density=True)
ax5.set_xlabel('ISF (mg/dL/U)'); ax5.set_ylabel('Density')
ax5.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P6: Bias by hour
ax6 = fig.add_subplot(gs[2, 1])
style(ax6, 'Bias by Hour of Night (+1h)')
for fname, flabel, fc in formulas:
    if f'err_{fname}_1h' in w1.columns:
        hb = w1.groupby('hour')[f'err_{fname}_1h'].mean()
        ax6.plot(hb.index, hb.values, 'o-', color=fc, lw=1.5, ms=4,
                 label=flabel.split('—')[0].strip())
ax6.axhline(0, color='white', lw=0.8, ls='--')
ax6.set_xlabel('Hour'); ax6.set_ylabel('Mean Error (mg/dL)')
ax6.legend(fontsize=7, labelcolor=TXT, facecolor=PANEL)

# P7: Within ±1mmol bar
ax7 = fig.add_subplot(gs[2, 2])
style(ax7, '% Within ±1 mmol/L Tolerance')
for i, horizon in enumerate(['1h','2h']):
    vals = [df_sum.query(f"formula=='{fl}' and horizon=='{horizon}'")['w18'].values[0]
            for _, fl, _ in formulas]
    bars = ax7.bar(x + (i-0.5)*w, vals, w,
                   color=[c for _,_,c in formulas], alpha=0.85 if i==0 else 0.5,
                   label=f'+{horizon}')
    for bar, v in zip(bars, vals):
        ax7.text(bar.get_x()+bar.get_width()/2, bar.get_height()+.3,
                 f'{v:.0f}%', ha='center', fontsize=6, color=TXT)
ax7.set_xticks(x); ax7.set_xticklabels(['v1','v2','No-TDD','7D-TDD','Flat'], fontsize=8)
ax7.set_ylabel('%'); ax7.legend(fontsize=8, labelcolor=TXT, facecolor=PANEL)

fig.suptitle(f'Overnight DynamicISF Backtest — 00:00–07:00  (Sep 2025–Mar 2026)\n'
             f'Median TDD = {tdd_median:.1f} U/day  |  Profile ISF = {overall_profile_isf:.1f} mg/dL/U  |  D = {D:.0f} mg/dL',
             color=TXT, fontsize=12, fontweight='bold', y=0.995)

plt.savefig('ns_backtest_results.png', dpi=150, bbox_inches='tight', facecolor=BG_C)
plt.close()
print("Saved: ns_backtest_results.png")

# ── 8. Save outputs ────────────────────────────────────────────────────────────
on.to_csv('ns_backtest_overnight.csv', index=False)

lines = [
    "OVERNIGHT DYNAMICISF BACKTEST — 00:00–07:00",
    "="*62,
    f"Data:            Sep 2025 – Mar 2026",
    f"Filter:          COB=0, bolus_age≥3h, BG 72–200 mg/dL, hour 0–7",
    f"Total cycles:    {len(on):,}",
    f"Clean 1h pairs:  {len(w1):,}   Clean 2h pairs: {len(w2):,}",
    f"Implied TDD:     median={tdd_median:.1f}  mean={tdd_mean:.1f}  IQR=[{on['tdd_implied'].quantile(.25):.1f}, {on['tdd_implied'].quantile(.75):.1f}] U/day",
    f"Profile ISF:     {overall_profile_isf:.1f} mg/dL/U (per-hour medians, BG 75–140 equilibrium filter)",
    "",
    "Formula definitions:",
    f"  v1 (Actual) = variable_sens  [= 1800/(TDD×ln(BG/D+1)), D={D:.0f}]",
    f"  v2          = 115000/(TDD²×ln(BG/D+1))  [TDD back-calculated from v1]",
    f"  No-TDD      = profileISF(h)×ln({TARGET:.0f}/D+1)/ln(BG/D+1)",
    f"  7D-TDD      = (1700/TDD_7day)×ln({TARGET:.0f}/D+1)/ln(BG/D+1)  [TDD_7day={tdd_7day_median:.1f} U/day median]",
    f"  Flat        = profileISF(h)",
    "",
]
for h in ['1h','2h']:
    lines.append(f"── +{h} ──────────────────────────────────────────────────────")
    for _, fl, _ in formulas:
        r = df_sum.query(f"formula=='{fl}' and horizon=='{h}'").iloc[0]
        lines.append(f"  {fl}")
        lines.append(f"    bias={r['bias']:+.1f}  MAE={r['mae']:.1f}  RMSE={r['rmse']:.1f}  within±1mmol={r['w18']:.1f}%  within±2mmol={r['w36']:.1f}%")
        if not np.isnan(r['ddir']):
            lines.append(f"    dosing direction correct: {r['ddir']:.1f}%")
    lines.append("")
lines.append("ISF magnitude (all overnight cycles):")
for fn, fl, _ in formulas:
    col = f'isf_{fn}'
    lines.append(f"  {fl:38s} median={on[col].median():.1f}  5-95pct=[{on[col].quantile(.05):.0f}, {on[col].quantile(.95):.0f}]")
with open('ns_backtest_summary.txt', 'w') as f:
    f.write('\n'.join(lines))

print('\n' + '\n'.join(lines))
print("\nSaved: ns_backtest_overnight.csv, ns_backtest_summary.txt")
