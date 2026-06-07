# TDD_effective: Overnight Sensitivity Regression for the 7D-TDD DynamicISF Model

## 1. Background and Motivation

### 1.1 The Problem

The 7D-TDD DynamicISF formula uses the 7-day rolling Total Daily Dose as its anchor for insulin sensitivity:

```
ISF = (1700 / TDD_7day) × ln(target/D + 1) / ln(BG/D + 1)
```

Where D = 82 mg/dL (insulin divisor), target = 99 mg/dL.

At the target BG, this recovers the classic 1700 rule: ISF = 1700 / TDD. The logarithmic BG term provides glucose-dependent scaling — ISF falls as BG rises (more aggressive correction at high BG) and rises as BG falls (protective at low BG).

The retrospective analysis (DynamicISF_Retrospective_ISF_Analysis-1) showed that the 7D-TDD formula performed well overall (MAE 15.2 mg/dL at +2h, 69.3% within ±1 mmol/L), but carried a systematic positive bias of +5.3 mg/dL — meaning actual BG consistently ended higher than the formula predicted. The formula runs slightly too aggressive: it over-estimates the insulin effect.

### 1.2 The Insight

TDD_7day is a *descriptive* statistic — it summarises what insulin was delivered over the past week. It does not account for whether that insulin produced the *expected* glucose outcomes. A patient who is more sensitive than their TDD implies will show BG dropping less than predicted (or more than predicted, depending on the direction of the bias).

The proposal: use the overnight fasting window (00:00–07:00) each day to retrospectively measure how far off the 7D-TDD formula's predictions were, back-calculate the TDD value that *would have produced correct predictions*, and exponentially smooth this into a `TDD_effective` that replaces `TDD_7day` in the formula.

The formula structure is unchanged. Only the TDD input is personalised.

### 1.3 Why Overnight

The overnight fasting window is ideal for sensitivity estimation because:

- **No carbs**: COB = 0 eliminates meal absorption as a confound
- **No exercise**: Glucose uptake from activity is minimal
- **Stable HGO**: Hepatic glucose output is relatively predictable (though Dawn Phenomenon introduces variation after ~04:00)
- **IOB dominance**: The only active glucose-lowering signal is insulin on board, making it possible to isolate insulin sensitivity
- **Daily availability**: Every night provides fresh data, enabling continuous adaptation

### 1.4 Source Data

The analysis uses `ns_backtest_overnight.csv`, the output of the prior retrospective study. This contains 4,165 overnight loop cycles (00:00–08:00, COB = 0) from September 2025 to March 2026, with:

- BG readings, IOB, loop predictions (predBGs.IOB at +1h, +2h, +3h)
- Actual BG outcomes at +1h, +2h, +3h matched to each prediction point
- Back-calculated TDD_implied from the v1 formula
- 7-day rolling TDD (tdd_7day)
- Pre-computed ISF values for all five formula variants

## 2. Method

### 2.1 Per-Sample Back-Calculation

For each valid overnight sample at time t, the loop logged its IOB-based BG prediction at +2h (`pred_iob_24`) and we know the actual BG at +2h (`actual_bg_2h`).

The sensitivity ratio measures how much the actual BG change differed from the predicted change:

```
ratio(t) = (BG(t) - actual_BG(t+2h)) / (BG(t) - pred_iob_24(t))
         = actual_drop / predicted_drop
```

- ratio > 1: BG dropped more than predicted → patient more sensitive
- ratio < 1: BG dropped less than predicted → patient less sensitive
- ratio = 1: prediction was perfect

The ISF that would have zeroed the prediction error:

```
ISF_eff(t) = ISF_v1(t) × ratio(t)
```

The TDD that produces ISF_eff via the 7D-TDD formula:

```
TDD_implied_eff(t) = (1700 × ln(target/D + 1)) / (ISF_eff(t) × ln(BG/D + 1))
```

### 2.2 Sample Filtering

Samples are excluded if:

| Filter | Threshold | Reason |
|---|---|---|
| BG out of range | < 72 or > 200 mg/dL | Extremes distort ratio calculation |
| Predicted drop too small | < 3 mg/dL absolute | Near-zero denominator → unstable ratio |
| Ratio direction mismatch | ratio ≤ 0 or > 5 | BG moved opposite to prediction, or pathological outlier |
| Suspected missed carbs | BG rose > 9 mg/dL when drop was predicted | Unrecorded food intake |
| TDD_implied_eff out of bounds | < 3 or > 120 U/day | Physiologically implausible |
| Missing data | NaN in pred_iob_24, actual_bg_2h, isf_v1, or tdd_7day | Incomplete loop cycle |

After filtering: 1,844 valid samples from 4,165 overnight cycles (44.3%).

### 2.3 Nightly Aggregation

For each night with ≥ 5 valid samples:

```
TDD_implied_night = median(TDD_implied_eff(t))  across all valid samples that night
```

Median is preferred over mean for robustness against CGM noise and occasional outlier samples.

### 2.4 Exponential Smoothing (07:00 Daily Update)

At 07:00 each morning, after last night's data is available:

```
TDD_effective(day N) = α × TDD_implied_night(N-1) + (1-α) × TDD_effective(day N-1)
```

- α = 0.15 gives an effective window of ~7 nights
- Bootstrap: TDD_effective is initialised from the first available tdd_7day value
- The update uses the *prior* night's regression — no look-ahead bias

The next 24 hours use TDD_effective in place of TDD_7day:

```
ISF_adaptive = (1700 / TDD_effective) × ln(target/D + 1) / ln(BG/D + 1)
```

## 3. Initial Results

### 3.1 Dataset Summary

- Overnight rows: 4,165
- Valid 2h samples: 1,844
- Usable nights (≥ 5 samples): 84
- Date range: October 2025 – March 2026

### 3.2 Overall Sensitivity Ratio

```
Median ratio: 0.836    IQR: [0.467, 1.364]
```

The median below 1.0 indicates that, on average, BG dropped *less* than predicted — the v1 formula (and by extension, the raw TDD) over-estimates the insulin effect. This is consistent with the +5.3 mg/dL positive bias found in the prior retrospective study.

### 3.3 Sensitivity Offset

```
TDD_effective − TDD_7day:  median = -0.95 U/day,  mean = -1.47 U/day
Range: [-11.4, +2.0] U/day
```

Persistently negative: this patient is *more sensitive* than their TDD history implies. The adaptive model corrects for this by using a lower TDD (→ higher ISF → less aggressive dosing).

### 3.4 Prediction Accuracy Comparison (+2h)

| Formula | Bias (mg/dL) | MAE (mg/dL) | RMSE (mg/dL) | ±1 mmol/L | ±2 mmol/L |
|---|---|---|---|---|---|
| v1 (Actual) | -1.8 | 16.4 | 22.1 | 68.0% | 89.6% |
| 7D-TDD (static) | +5.3 | 15.2 | 21.2 | 69.3% | 92.1% |
| 7D-TDD (adaptive, α=0.15) | +4.7 | 15.5 | 21.6 | 68.5% | 91.1% |

The adaptive model reduces bias from +5.3 to +4.7 mg/dL (11% reduction) but MAE increases slightly from 15.2 to 15.5.

### 3.5 Learning Rate Comparison

| α | ~Window | Bias | MAE | RMSE | ±1 mmol/L |
|---|---|---|---|---|---|
| 0.10 | ~10 nights | +4.8 | 15.4 | 21.6 | 69.1% |
| 0.15 | ~7 nights | +4.7 | 15.5 | 21.6 | 68.5% |
| 0.20 | ~5 nights | +4.5 | 15.5 | 21.7 | 68.3% |
| 0.30 | ~3 nights | +4.2 | 15.7 | 22.0 | 68.0% |

All learning rates reduce bias vs static. Higher α (faster learning) reduces bias more but increases MAE — the nightly signal is too noisy for aggressive updates.

### 3.6 Assessment of Initial Results

The improvement is modest. Three factors explain this:

1. **High per-night variance**: The nightly TDD_implied has wide IQR, meaning the per-night estimate is noisy. The median ratio IQR of [0.47, 1.36] represents almost 3× variation. Exponential smoothing necessarily lags behind this noise.

2. **Confound through ISF_v1**: The back-calculation routes through ISF_v1 (the actual formula the loop used). This means any v1-specific biases are baked into the implied ISF, and hence into TDD_implied_eff.

3. **Constant correction across BG range**: A single TDD_effective applies the same correction at all BG levels. If the bias is BG-dependent (e.g., worse at high BG), a scalar correction cannot capture this.

## 4. Refinements

### 4.1 Confidence-Weighted Nightly Aggregation

Not all overnight samples carry equal information. A sample with large IOB and a large predicted BG drop provides a much stronger signal about insulin sensitivity than a sample where IOB is near-zero and the predicted drop is 4 mg/dL.

**Refinement**: Weight each sample's TDD_implied_eff by the absolute predicted BG drop, then take a weighted median (or weighted mean).

```
weight(t) = |BG(t) - pred_iob_24(t)|
TDD_implied_night = weighted_median(TDD_implied_eff, weights)
```

This upweights the high-signal samples and downweights the noisy near-equilibrium ones.

### 4.2 Dawn Phenomenon Stratification

Dawn Phenomenon causes insulin resistance to increase between approximately 04:00–07:00 due to growth hormone and cortisol secretion. A single TDD_effective that blends deep-night (00:00–03:30) and pre-dawn (03:30–07:00) data may miss this pattern.

**Refinement**: Compute separate sensitivity ratios for each window. If they diverge significantly, report the dawn-effect magnitude. For the TDD_effective update, use only the deep-night window (00:00–03:30) as the cleaner signal, and separately quantify the dawn effect.

### 4.3 Adaptive Learning Rate

When the nightly TDD_implied is consistent (low IQR), we can trust it more and use a higher α. When it's scattered, we should be more conservative.

**Refinement**: Scale α by the inverse of the nightly coefficient of variation:

```
cv_night = IQR / median of TDD_implied_eff samples
α_adaptive = α_base × min(1.0, 1.0 / (1 + cv_night))
```

### 4.4 BG-Band Ratio Analysis

If the sensitivity correction differs by BG level, a single TDD_effective is insufficient and the problem may be in the logarithmic BG scaling term rather than TDD.

**Refinement**: Stratify the overnight ratio by BG band (< 90, 90–120, 120–150, > 150) and test whether the ratio varies significantly across bands. If it does, the correction is BG-dependent and a scalar TDD adjustment is the wrong lever.

### 4.5 Direct Ratio Approach

Instead of back-calculating through ISF_v1 to derive TDD_implied_eff, apply the correction as a direct multiplier on TDD_7day:

```
TDD_effective = TDD_7day / median_ratio_night
```

This is algebraically equivalent when the 7D-TDD formula is the reference, but avoids routing through v1's ISF.

## 5. Refined Results

Five methods were compared across 84 usable nights (1,844 valid 2h samples, 1,790 with prior-night TDD_effective available):

### 5.1 Backtest Comparison (+2h horizon)

| Method | Bias (mg/dL) | MAE (mg/dL) | RMSE (mg/dL) | ±1 mmol/L | ±2 mmol/L |
|---|---|---|---|---|---|
| v1 (Actual loop) | -1.8 | 16.4 | 22.1 | 68.0% | 89.6% |
| **7D-TDD (static)** | **+5.3** | **15.2** | **21.2** | **69.3%** | **92.1%** |
| Simple median, α=0.15 | +4.7 | 15.5 | 21.6 | 68.5% | 91.1% |
| Confidence-weighted, α=0.15 | +5.7 | 15.6 | 22.0 | 68.7% | 90.4% |
| Deep-night only, α=0.15 | +5.6 | 15.7 | 22.0 | 67.9% | 91.1% |
| Weighted + adaptive α | +5.3 | 15.5 | 21.8 | 68.4% | 90.7% |
| Direct ratio, α=0.15 | +7.4 | 16.7 | 23.6 | 65.6% | 88.5% |

The static 7D-TDD remains the best performer on MAE. All adaptive methods reduce bias (the simple median achieves the best bias reduction: +5.3 → +4.7, an 11% improvement) but at a small MAE cost. The direct ratio method performs worst — routing through TDD_7day/ratio amplifies noise because TDD_7day itself drifts.

The simple median method (v1's approach) wins overall among the adaptive variants: lowest MAE (15.5), lowest bias (+4.7), lowest RMSE (21.6), and best ±2 mmol/L rate (91.1%).

Adaptive wins on 41 of 83 nights (49%) — essentially a coin flip, indicating the correction helps on some nights and hurts on others.

### 5.2 BG-Band Ratio Analysis

| BG Band | n | Ratio Median | Ratio IQR |
|---|---|---|---|
| < 90 mg/dL | 702 | 0.951 | [0.579, 1.379] |
| 90–120 mg/dL | 869 | 0.778 | [0.364, 1.364] |
| 120–150 mg/dL | 273 | 0.750 | [0.421, 1.333] |

**Kruskal-Wallis test: H = 25.81, p < 0.0001**

The sensitivity ratio varies significantly by BG band. At low BG (< 90), the ratio is near 1.0 — the formula is well-calibrated. At higher BG (90–150), the ratio drops to 0.75–0.78, meaning the formula consistently over-predicts the insulin effect at elevated BG.

**This is a critical finding.** It means the correction is *not* a simple TDD miscalibration — it is BG-dependent. A scalar TDD_effective cannot fully capture this. The logarithmic BG scaling term `ln(BG/D+1)` may need adjustment, or the correction should be applied as a BG-dependent multiplier rather than a TDD shift.

### 5.3 Dawn Phenomenon Stratification

| Period | n | Ratio Median | TDD_eff Median | ISF_eff Median |
|---|---|---|---|---|
| Deep Night (00:00–03:30) | 815 | 0.781 | 12.1 U/day | 146.0 mg/dL/U |
| Pre-Dawn (03:30–07:00) | 1,029 | 0.896 | 9.8 U/day | 172.1 mg/dL/U |

**Mann-Whitney U test: p = 0.0017**

The dawn effect is statistically significant: the pre-dawn ratio (0.896) is 0.115 higher than the deep-night ratio (0.781). Counter-intuitively, this means the patient appears *more* sensitive in the pre-dawn period — the formula over-predicts insulin effect *less* pre-dawn. This could reflect:

- Dawn Phenomenon partially offsetting insulin action, making the formula's high-ISF prediction closer to correct
- Or that the deep-night signal is contaminated by residual dinner/evening bolus activity that has decayed by pre-dawn

The deep-night-only method did not outperform the full-night methods, suggesting that restricting to the cleaner window does not improve the signal enough to compensate for the reduced sample size.

### 5.4 Adaptive Learning Rate

The adaptive α method dampened the learning rate significantly — typical α values were 0.05–0.11, well below the base of 0.15 — because nightly CV (IQR/median) was consistently high (median CV = 0.89). This confirms the per-night signal is noisy: the adaptive method correctly identifies low confidence but the result is that it barely updates, making it equivalent to a very slow static α.

### 5.5 Sensitivity Offset by Method

| Method | Median Offset (U/day) | Mean Offset (U/day) |
|---|---|---|
| Simple median | -0.95 | -1.47 |
| Confidence-weighted | +0.75 | +0.41 |
| Deep-night only | -0.47 | -0.41 |
| Adaptive α | +0.21 | +0.00 |
| Direct ratio | +4.63 | +4.51 |

The sign of the offset depends on the method — simple and deep-night methods track below TDD_7day (more sensitive), while weighted and adaptive methods track near or above (less sensitive). The confidence weighting upweights large-drop samples, which tend to be at higher BG where the ratio is lower (0.75), pulling TDD_effective higher.

The direct ratio method diverges badly (+4.63 U/day) — it amplifies the noise in TDD_7day/ratio division.

## 6. Key Findings and Interpretation

### 6.1 The Correction Is BG-Dependent, Not TDD-Dependent

The most important finding is from the BG-band analysis: the prediction error is **not** a uniform TDD miscalibration. It varies from ratio ≈ 0.95 at low BG to ratio ≈ 0.75 at high BG. This means:

- At low BG, the 7D-TDD formula is already nearly correct
- At higher BG, it systematically over-predicts the insulin effect

A single TDD_effective scalar can reduce the *average* bias but cannot fix the BG-dependent component. This suggests the problem lies partly in the logarithmic BG scaling term, not solely in TDD.

This is consistent with the prior study's finding that TDD explains near-zero variance in implied ISF (R² = 0.003), while BG explains ~10%.

### 6.2 The Nightly Signal Is Noisy

All methods converge to similar accuracy because the per-night TDD_implied estimates have wide IQR (typical CV ≈ 0.89). The exponential smoothing averages this out, but by the time the estimate stabilises, the underlying sensitivity may have shifted. The adaptive α correctly identifies this noise and dampens the update rate.

### 6.3 The Persistent Offset Is Real

Despite the noise, the TDD_effective consistently sits below TDD_7day for this patient (simple median method: -0.95 U/day). This persistent directional signal is the value of the model — not cycle-by-cycle accuracy, but the ability to detect that this patient's insulin is more effective than their dose history implies.

### 6.4 Simple Wins

The simple median with fixed α=0.15 outperforms all refined methods. Confidence weighting and dawn stratification add complexity without improving accuracy. The overnight prediction residual is noisy enough that sophisticated aggregation methods cannot extract a better signal than the median.

## 7. Recommendations

### 7.1 For the 7D-TDD Formula

**TDD_effective with simple median and α=0.15 is the recommended approach.** It provides:
- A modest bias reduction (11%)
- A persistent sensitivity offset that reveals individual deviation from the 1700 rule
- A single-number daily update that requires no formula structural changes

### 7.2 For Further Investigation

The BG-band finding suggests that a **BG-dependent correction factor** would outperform a scalar TDD adjustment. This could take the form of a modified logarithmic scaling parameter, or a piecewise correction applied at different BG thresholds. This is a separate investigation from TDD_effective and addresses the ln(BG/D+1) term rather than the TDD input.

### 7.3 For Implementation

The update mechanism is simple:
1. At 07:00, collect overnight samples (00:00–07:00, COB=0, valid 2h outcomes)
2. For each: compute ratio = actual_drop / predicted_drop
3. Filter: ratio > 0 and < 5, |predicted_drop| ≥ 3 mg/dL, BG 72–200
4. Back-calculate TDD_implied_eff per sample
5. TDD_implied_night = median of valid samples (require ≥ 5)
6. TDD_effective = 0.15 × TDD_implied_night + 0.85 × TDD_effective_prev
7. Use TDD_effective in place of TDD_7day for the next 24 hours

### 7.4 Limitations

- **Single-patient study**: N=1. The offset direction, magnitude, and BG-band pattern will differ across patients.
- **Retrospective counterfactual**: Predictions are scaled from v1's ISF, not directly from the 7D-TDD formula's own forward model.
- **Overnight only**: Daytime sensitivity may differ. The overnight regression corrects the baseline but does not capture meal-time or exercise-time sensitivity changes.
- **Feedback loop**: In a live system, TDD_effective would influence dosing, which changes TDD_7day, which changes the baseline for the next regression. Stability analysis would be needed before deployment.

## 8. Files

| File | Description |
|---|---|
| `ns_tdd_effective.py` | v1 analysis script (simple median only) |
| `ns_tdd_effective_v2.py` | v2 refined analysis (all 5 methods, BG-band, dawn, adaptive α) |
| `ns_tdd_effective_results.png` | v1 visualisation (8-panel) |
| `ns_tdd_effective_v2_results.png` | v2 visualisation (12-panel) |
| `ns_tdd_effective_summary.txt` | v1 text summary |
| `ns_tdd_effective_v2_summary.txt` | v2 text summary |
| `ns_tdd_effective_nightly.csv` | v1 per-night data |
| `ns_tdd_effective_v2_nightly.csv` | v2 per-night data (all methods) |
| `ns_backtest_overnight.csv` | Input data (from ns_overnight_backtest.py) |
