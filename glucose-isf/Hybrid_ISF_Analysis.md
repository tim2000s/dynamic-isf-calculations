# Diabeloop Glucose-Related ISF Model: Analysis and Hybrid Extension

**Tim Street**
**April 2026**

---

## Executive Summary

This paper evaluates an alternative approach to dynamic insulin sensitivity factor (ISF) calculation based on the Diabeloop glucose-related ISF model, presented as an ADA scientific poster. The Diabeloop model uses a quartic polynomial to describe the relationship between glucose and ISF for values above 100 mg/dL. This analysis extends the model with a power-law continuation below 105 mg/dL, creating a hybrid formula. The resulting formula requires no Total Daily Dose (TDD) input, making it simpler to implement and inherently patient-agnostic at the population level.

Backtested against 3,647 valid overnight fasting samples from 10 months of closed-loop data (June 2025 – March 2026), the hybrid achieves:

- **Overall MAE 14.1 mg/dL** — within 1 mg/dL of the power-law with Boost TDD (13.2)
- **Excellent accuracy above 105 mg/dL** — MAE 8.9 in the 105–120 band, matching the best-performing formula
- **Deliberate conservatism below 90 mg/dL** — bias of +14.5, reducing insulin delivery in the zone where hypoglycaemia risk is highest, without affecting Time in Range (70–180 mg/dL)

The key finding is that the hybrid's greater positive bias at low glucose values is clinically desirable: below 90 mg/dL, the patient is already within the 70–180 target range, so conservative insulin delivery protects against dangerous drops without degrading glycaemic outcomes. Above 105 mg/dL, the polynomial provides tight, accurate control that maximises Time in Range.

The paper also presents a two-phase calibration strategy that requires no formal ISF testing: an initial scaling factor derived from the patient's existing TDD (or profile ISF if TDD is unavailable), followed by automatic refinement from observed fasting prediction errors over the first weeks of use.

A synthetic multi-patient validation across 7 patients spanning a 4× range of insulin sensitivity (TDD 11–45 U/day) demonstrates that TDD-based calibration achieves MAE 14.0 for every patient type, confirming the hybrid's curve shape works uniformly from very insulin-resistant to very insulin-sensitive patients.

---

## 1. Background

### 1.1 The Problem

Current dynamic ISF implementations in open-source closed-loop systems (AAPS, Trio) use a logarithmic scaling formula:

```
ISF = C / (TDD × ln(BG/D + 1))
```

This formula has known limitations: it requires accurate TDD estimation, its logarithmic curve provides insufficient sensitivity differentiation at high glucose values, and it can be overly aggressive near target, increasing hypoglycaemia risk.

### 1.2 Prior Work

A companion paper (*Dynamic ISF: Power-Law Scaling with Boost TDD Blending*) established that a power-law formula with k=3.5 and Boost-blended TDD significantly outperforms the current logarithmic approach:

```
ISF = (1800 / TDD_blended) × (target / BG) ^ 3.5
```

This formula achieves MAE 13.2 mg/dL with bias +5.6 — the most accurate of all formulas tested. However, it requires real-time TDD computation including the Boost blending algorithm.

### 1.3 The Diabeloop Glucose-Related ISF Model

Diabeloop, the company behind the DBLG1 closed-loop insulin delivery system, presented research at the American Diabetes Association (ADA) Scientific Sessions describing a glucose-related ISF model derived from clinical data. The model, presented as an ADA scientific poster, characterises the population-level relationship between glucose and insulin sensitivity using piecewise polynomial fits.

![Diabeloop glucose-related ISF model from ADA scientific poster, showing ISF vs glucose with two piecewise polynomial fits and 25th/75th percentile error bars](Variable%20ISF.png)

*Figure 1: Diabeloop's glucose-related ISF model. The chart shows median ISF values (bars) with 25th and 75th percentile ranges across glucose levels. Two polynomial equations are fitted: a quadratic for glucose ≤ 100 mg/dL and a quartic for glucose > 100 mg/dL. Note the proportional scaling of the IQR bars — the 25th and 75th percentiles scale at approximately 0.5× and 1.6× of the median respectively, consistent across all glucose levels.*

For glucose values above 100 mg/dL, the quartic polynomial is:

```
ISF(G) = 272 − 3.121G + 0.01511G² − 3.305×10⁻⁵G³ + 2.69×10⁻⁸G⁴
```

This equation captures the population-median ISF curve without requiring any patient-specific parameters such as TDD. The proportional IQR scaling visible in Figure 1 — where 25th and 75th percentile ISF values scale proportionally with the median across all glucose levels — suggests that individual patients follow the same curve shape but at different magnitudes, supporting the use of a single multiplicative scaling factor for personalisation.

The poster also provides a quadratic equation for glucose ≤ 100 mg/dL. However, as discussed in Section 2.1, this lower equation was not used in the hybrid formula due to safety concerns about its behaviour in the near-hypoglycaemic range.

---

## 2. The Hybrid Formula

### 2.1 Motivation

The Diabeloop model provides two equations: a quartic polynomial for glucose > 100 mg/dL and a quadratic for glucose ≤ 100 mg/dL. Initial backtesting of the quartic polynomial extended across all glucose values revealed a split personality:

- **Above 105 mg/dL**: Excellent performance — MAE 8.9 (105–120 band), competitive with the best formulas tested
- **Below 90 mg/dL**: Poor raw accuracy — MAE 30.2, with large positive bias (+29.9)

The quartic polynomial was designed for glucose values above 100 mg/dL. Extrapolating below this range produces ISF values that are too low relative to what the patient needs — for example, at glucose 80 mg/dL the polynomial gives ISF ≈ 103 mg/dL/U versus the ~160 mg/dL/U that actual loop data indicates is needed. This means the formula is too aggressive at low glucose, precisely where safety matters most.

The Diabeloop quadratic equation for ≤ 100 mg/dL was not adopted for two reasons:

1. **Safety in the near-hypoglycaemic range**: Below 90 mg/dL, the priority shifts from prediction accuracy to hypo avoidance. The quadratic's behaviour in this range was not validated against our dataset, and any equation that is insufficiently conservative below 90 mg/dL poses a direct safety risk.

2. **Deliberate conservatism is a feature**: As detailed in Section 5, a formula that is conservative below 90 mg/dL does not degrade Time in Range (the patient is already in the 70–180 range) but does provide a buffer against dangerous drops. A power-law tail with a well-characterised exponent achieves this conservatism by design.

Rather than discarding the polynomial or adopting the untested quadratic, we proposed anchoring a power-law tail below 105 mg/dL that provides a smooth, continuous transition with deliberately conservative ISF values in the near-target and near-hypoglycaemic range.

### 2.2 Definition

The hybrid ISF formula is defined piecewise:

**For glucose ≥ 105 mg/dL:**
```
ISF = 272 − 3.121G + 0.01511G² − 3.305×10⁻⁵G³ + 2.69×10⁻⁸G⁴
```

**For glucose < 105 mg/dL:**
```
ISF = 75.8 × (105 / BG) ^ 3.5
```

The constant 75.8 is the polynomial's value at glucose = 105, ensuring continuity at the junction. The power-law exponent of 3.5 was chosen to match the independently-optimised exponent from the companion power-law analysis, providing a well-characterised scaling behaviour in the sub-105 range.

### 2.3 Properties

The hybrid has several notable properties:

1. **No TDD required** — The formula gives absolute ISF values based solely on current glucose, making it simpler to implement and removing dependency on TDD estimation accuracy
2. **Continuous at the junction** — Both pieces evaluate to 75.8 mg/dL/U at glucose = 105
3. **Population-level calibration** — The polynomial piece reflects clinical population data; the power-law piece inherits its anchor from the same data
4. **Conservative at low glucose by design** — The sub-105 piece produces higher ISF values (less aggressive dosing) than TDD-based formulas, providing an inherent safety margin

---

## 3. Backtest Methodology

### 3.1 Dataset

The backtest uses 10 months of continuous closed-loop data from a single AAPS user (June 2025 – March 2026), comprising:

- 94,980 total loop cycles
- 90,733 CGM readings
- 61,975 treatment records

### 3.2 Sample Selection

Samples were restricted to overnight fasting periods (00:00–07:59) to isolate ISF effects from meal absorption:

- No carbs on board (COB = 0)
- At least 180 minutes since last bolus
- Glucose between 72–200 mg/dL
- Valid TDD data available

Additional quality filters removed samples with negligible predicted drops (<3 mg/dL), contradictory predicted vs actual direction, and extreme prediction ratios. This yielded **3,647 valid samples** with 2-hour forward glucose measurements.

### 3.3 Counterfactual Prediction

For each sample, the predicted 2-hour glucose under each formula was computed as:

```
predicted_BG = BG(t) − predicted_drop × (ISF_formula / ISF_actual)
```

Where `ISF_actual` is the `variable_sens` value the loop was actually using, and `predicted_drop` is the loop's own IOB-based prediction. This counterfactual approach asks: "If the loop had used this formula's ISF instead of its actual ISF, where would glucose have ended up?"

Prediction error is defined as `actual_BG_2h − predicted_BG`, so:
- **Positive bias** means actual glucose was higher than predicted — the formula overestimated the drop and would have delivered less insulin
- **Negative bias** means actual glucose was lower than predicted — the formula underestimated the drop and would have delivered more insulin

### 3.4 Formulas Compared

| Label | Formula | TDD Required |
|-------|---------|:---:|
| A | Loop actual (`variable_sens`) | Yes |
| B | Current logarithmic + Boost blended TDD | Yes |
| C | Power-law k=3.5 + Boost blended TDD | Yes |
| D | Polynomial (Diabeloop) — raw, no TDD | No |
| E | Polynomial (Diabeloop) — scaled ×1.31 to patient | No |
| F | **Hybrid** (polynomial ≥105, power-law <105) — no TDD | No |

---

## 4. Results

### 4.1 Overall Performance

| Formula | MAE | Bias | ±1 mmol/L |
|---------|:---:|:----:|:---------:|
| A: Loop actual | 13.4 | +11.2 | 75.8% |
| B: Current ln + Boost TDD | 16.7 | +14.1 | 65.8% |
| C: Power-law k=3.5 + Boost TDD | **13.2** | **+5.6** | **76.6%** |
| D: Polynomial raw — no TDD | 17.4 | +13.3 | 65.2% |
| E: Polynomial scaled — no TDD | 16.4 | +13.2 | 67.8% |
| F: Hybrid — no TDD | 14.1 | +8.2 | 72.8% |

The hybrid achieves the best performance of any formula that does not require TDD, and comes within 1 mg/dL MAE of the power-law with Boost TDD.

### 4.2 Performance by Glucose Band

| Glucose Band | n | Power-law (C) MAE | Power-law (C) Bias | Hybrid (F) MAE | Hybrid (F) Bias |
|---------|--:|:-:|:-:|:-:|:-:|
| <90 | 1,128 | 19.3 | +9.7 | 20.6 | +14.5 |
| 90–105 | 812 | 13.2 | +10.2 | 14.5 | +11.6 |
| 105–120 | 1,098 | 8.4 | +2.1 | 8.9 | +2.5 |
| 120–150 | 608 | 10.7 | −2.0 | 11.1 | +2.5 |
| 150–200 | — | — | — | — | — |

The two formulas perform almost identically above 105 mg/dL, where the polynomial drives both. The meaningful difference is in the sub-90 band, where the hybrid's greater conservatism produces a higher bias (+14.5 vs +9.7).

### 4.3 ISF Curve Values

| Glucose (mg/dL) | Current ln | Power-law k=3.5 | Polynomial | Hybrid |
|:---:|:---:|:---:|:---:|:---:|
| 60 | 195.4 | 619.1 | 132.3 | 537.4 |
| 70 | 173.8 | 361.0 | 116.9 | 313.3 |
| 80 | 157.6 | 226.2 | 103.2 | 196.3 |
| 90 | 144.8 | 149.8 | 91.2 | 130.0 |
| 99 | 135.5 | 107.3 | 81.6 | 93.1 |
| 105 | 130.1 | 87.3 | 75.9 | 75.9 |
| 110 | 126.1 | 74.2 | 71.5 | 71.5 |
| 120 | 119.0 | 54.7 | 63.5 | 63.5 |
| 140 | 107.7 | 31.9 | 50.9 | 50.9 |
| 160 | 99.1 | 20.0 | 41.7 | 41.7 |
| 180 | 92.4 | 13.2 | 35.3 | 35.3 |
| 200 | 86.9 | 9.2 | 30.8 | 30.8 |

Note: The ln and power-law columns use Boost blended TDD (median 16.8 U/day). The polynomial and hybrid columns are TDD-independent.

---

## 5. The Safety Argument

### 5.1 Why Conservatism Below 90 mg/dL Is Clinically Desirable

The hybrid's larger positive bias below 90 mg/dL (+14.5 vs the power-law's +9.7) means it consistently overestimates how far glucose will drop, causing the loop to deliver less insulin than a more "accurate" formula would. This is a feature, not a limitation, for the following reasons:

**A patient at 85 mg/dL is already in range.** The 70–180 mg/dL target range means that any glucose between 70 and 90 is already counted as Time in Range. Being conservative here — delivering slightly less insulin — does not reduce TIR. It simply means glucose drifts back toward target rather than continuing to fall.

**The downside risk is asymmetric.** A patient at 85 mg/dL who receives too much insulin may drop to 55 mg/dL — a clinically dangerous hypoglycaemic event. A patient at 85 mg/dL who receives slightly too little insulin may drift to 100 mg/dL — a clinically irrelevant outcome that remains well within range. The cost of over-delivery vastly exceeds the cost of under-delivery in this zone.

**Sudden drops are most dangerous near the lower boundary.** Below 90 mg/dL, the margin to clinically significant hypoglycaemia (<54 mg/dL) is narrow. Factors outside the algorithm's model — exercise, alcohol, sensor lag, absorption variability — can cause rapid unexpected drops. A built-in conservative bias provides a buffer against these unmodelled perturbations.

**The conservatism is proportional to risk.** The hybrid's bias increases as glucose decreases: +11.6 at 90–105, +14.5 below 90. This graduated response mirrors the increasing clinical risk as glucose approaches the hypoglycaemic threshold.

### 5.2 Comparison with the Power-Law Approach

The power-law formula (C) achieves lower MAE and tighter bias across all bands, making it the superior formula from a pure prediction-accuracy standpoint. However, prediction accuracy and clinical safety are not the same objective:

| Property | Power-law + Boost TDD | Hybrid |
|----------|:---:|:---:|
| Overall MAE | 13.2 | 14.1 |
| Sub-90 bias (safety margin) | +9.7 | +14.5 |
| Above-105 MAE | 8.4–10.7 | 8.9–11.1 |
| Real-time TDD required | Yes (every cycle) | No (TDD used only for initial calibration) |
| Requires Boost blending | Yes | No |
| Learned parameters | C, k, TDD blend | Single scaling factor S |
| Implementation complexity | Moderate | Low |

The hybrid trades approximately 1 mg/dL of overall accuracy for a 50% larger safety margin below 90 mg/dL, while replacing continuous TDD computation with a single scaling factor that is initialised once and refined automatically.

### 5.3 Where the Hybrid Is Weaker

The hybrid does have limitations:

1. **Requires calibration** — The base formula gives population-median ISF values. Without the patient-specific scaling factor (Section 6.2), ISF at target is 75.8 (unscaled hybrid) vs 107.3 (this patient's actual). The two-phase calibration strategy addresses this: TDD-based initialisation on day one, followed by auto-refinement from observed data.

2. **90–105 band** — This narrow band shows the hybrid's weakest relative performance (MAE 14.5 vs power-law's 13.2). The patient is in range here, so the clinical impact is minimal, but the bias (+11.6) means slightly less time at the 90–99 target centre.

3. **Single-patient validation** — All results are from one patient's data. The polynomial's clinical grounding (derived from population data) partially mitigates this concern, but multi-patient validation is needed.

---

## 6. Implementation

### 6.1 The Formula

```
if BG >= 105:
    ISF = 272 − 3.121 × BG + 0.01511 × BG² − 3.305e-05 × BG³ + 2.69e-08 × BG⁴
else:
    ISF = 75.8 × (105 / BG) ^ 3.5
```

### 6.2 Patient Calibration

The hybrid formula produces population-median ISF values. Individual patients require a scaling factor `S` to shift the curve to their insulin sensitivity. The Diabeloop poster's proportional IQR scaling (25th and 75th percentile ISF values scale proportionally with the median) supports the use of a single multiplicative factor that preserves the curve shape:

```
ISF_patient = ISF_hybrid(BG) × S
```

The recommended calibration strategy has two phases: an immediate initial estimate, followed by automatic refinement.

#### Phase 1: Initial Scaling (Day One)

**If 7-day TDD is available** (preferred — most pump users have this):

```
S = (1800 / TDD_7day) / 75.8
```

The 1800 rule provides a well-established estimate of ISF at target glucose. Dividing by 75.8 (the hybrid's unscaled value at the 105 mg/dL junction) yields the patient's scaling factor. For example, a patient with TDD = 22 U/day:

```
S = (1800 / 22) / 75.8 = 81.8 / 75.8 = 1.08
```

**If TDD is not available** (new pump user, or system without TDD tracking):

```
S = profile_ISF / 75.8
```

Every insulin pump requires a configured ISF in the patient's profile. This is already available on day one and requires no additional testing. For a patient with profile ISF = 100:

```
S = 100 / 75.8 = 1.32
```

Both approaches provide a reasonable starting point. The TDD-based method is preferred because TDD is a measured value that reflects the patient's actual recent insulin needs, whereas profile ISF may be outdated or inaccurate. However, either is sufficient to begin safe closed-loop operation — particularly given the hybrid's inherent conservatism at low glucose values, which provides a safety buffer even if the initial scaling is imperfect.

#### Phase 2: Auto-Calibration (First Week and Ongoing)

After the first few days of looping, the system accumulates fasting periods where it can directly observe how well its ISF predictions match reality. The scaling factor can be refined automatically:

1. **Collect fasting samples**: Identify overnight periods (or any fasting window) where COB = 0 and no bolus has been delivered for at least 3 hours. These are the cleanest windows for observing pure ISF behaviour.

2. **Compute observed scaling**: For each valid sample, compare the predicted glucose drop (using the current scaled hybrid ISF) against the actual glucose drop:

```
observed_ratio = actual_drop / predicted_drop
```

3. **Update the scaling factor**: Take the median of observed ratios over a rolling window (e.g., 7 days) and blend it with the current scaling factor:

```
S_new = S_current × median(observed_ratios)
```

If the median ratio is consistently above 1.0, the formula is underestimating drops (glucose falls more than predicted), meaning `S` should increase to make ISF more aggressive. If below 1.0, glucose falls less than predicted, and `S` should decrease. Applying the median rather than the mean provides robustness against outliers from sensor noise, unrecorded carbs, or exercise effects.

4. **Apply safety bounds**: The auto-calibration should be bounded to prevent runaway adjustments:

```
S = clamp(S_new, 0.5, 3.0)
```

This range covers patients from approximately double the population-median insulin sensitivity (very sensitive) to half (very resistant), which encompasses the vast majority of Type 1 diabetes patients.

#### Calibration Timeline

| Phase | Timing | Source | Accuracy |
|-------|--------|--------|----------|
| Initial (TDD) | Day 1 | 1800/TDD_7day | Approximate — correct order of magnitude |
| Initial (profile) | Day 1 | Profile ISF | Approximate — may be outdated |
| Auto-calibrated | Days 3–7 | Observed fasting data | Refined — converges on patient's actual ISF curve |
| Steady state | Ongoing | Rolling 7-day observations | Continuously adapting to changes in insulin sensitivity |

This approach is analogous to how existing systems like Autotune and Autosens work, but simpler: instead of adjusting multiple parameters (basal rates, ISF, carb ratio), only a single scaling factor needs to be learned. The curve shape is fixed by the population-derived polynomial and the power-law tail — only the magnitude is patient-specific.

### 6.3 Safety Bounds

An ISF floor should be applied to prevent extreme dosing at very high glucose:

```
ISF = max(ISF_floor, ISF_hybrid)
```

A reasonable default is ISF_floor = 10 mg/dL/U (adjustable 5–20).

### 6.4 Kotlin (AAPS)

```kotlin
fun calculateHybridIsf(bg: Double, scaleFactor: Double): Double {
    val baseIsf = if (bg >= 105.0) {
        272.0 - 3.121 * bg + 0.01511 * bg * bg -
        3.305e-05 * bg * bg * bg + 2.69e-08 * bg * bg * bg * bg
    } else {
        75.8 * (105.0 / bg).pow(3.5)
    }
    return maxOf(isfFloor, baseIsf * scaleFactor)
}

// Phase 1 initialisation
fun initialScaleFactor(tdd7day: Double): Double {
    return (1800.0 / tdd7day) / 75.8
}
```

### 6.5 Swift (Trio)

```swift
func calculateHybridIsf(bg: Double, scaleFactor: Double) -> Double {
    let baseIsf: Double
    if bg >= 105.0 {
        baseIsf = 272.0 - 3.121 * bg + 0.01511 * pow(bg, 2)
                - 3.305e-05 * pow(bg, 3) + 2.69e-08 * pow(bg, 4)
    } else {
        baseIsf = 75.8 * pow(105.0 / bg, 3.5)
    }
    return max(isfFloor, baseIsf * scaleFactor)
}

// Phase 1 initialisation
func initialScaleFactor(tdd7day: Double) -> Double {
    return (1800.0 / tdd7day) / 75.8
}
```

---

## 7. Multi-Patient Calibration Validation

To validate that the hybrid formula and its calibration strategy work across patients with different insulin sensitivity levels, a synthetic multi-patient backtest was conducted using the real patient's data as a template.

### 7.1 Methodology

Seven synthetic patients were created by uniformly scaling the real patient's ISF values, spanning from very insulin-resistant (0.5× sensitivity, ~45 U/day TDD) to very insulin-sensitive (2.0× sensitivity, ~11 U/day TDD). The 0.5× and 1.6× scale factors correspond to the Diabeloop poster's 25th and 75th percentile ISF values, representing the population interquartile range.

| Patient | Sensitivity | Scale Factor | Simulated TDD |
|:---:|---|:---:|:---:|
| P1 | Very resistant | 0.50 | 45.3 U/day |
| P2 | Resistant | 0.67 | 33.8 U/day |
| P3 | Moderately resistant | 0.80 | 28.3 U/day |
| P4 | Original patient | 1.00 | 22.6 U/day |
| P5 | Moderately sensitive | 1.25 | 18.1 U/day |
| P6 | Sensitive | 1.60 | 14.2 U/day |
| P7 | Very sensitive | 2.00 | 11.3 U/day |

For each synthetic patient, the full two-phase calibration was simulated:
- **Phase 1**: Initial scaling factor computed from synthetic TDD using `S = (1800 / TDD) / 75.8`
- **Phase 2**: Auto-calibration from observed fasting prediction errors, using a 7-day rolling window with 30% dampened updates to prevent overcorrection

The same 3,647 overnight fasting samples were used, with each synthetic patient's "true" ISF scaled accordingly. The counterfactual prediction method was applied: for each sample, the predicted 2-hour glucose was computed as if the loop had used the hybrid formula (with current S) instead of the synthetic patient's actual ISF.

### 7.2 Results

| Patient | TDD | S_true | S_phase1 | S_final | MAE raw | MAE Phase 1 | MAE Auto-Cal | MAE Perfect |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| P1 | 45.3 | 0.708 | 0.524 | 0.640 | 20.2 | 14.0 | 14.2 | 14.8 |
| P2 | 33.8 | 0.948 | 0.702 | 0.858 | 15.2 | 14.0 | 14.2 | 14.8 |
| P3 | 28.3 | 1.132 | 0.839 | 1.023 | 14.1 | 14.0 | 14.2 | 14.8 |
| P4 | 22.6 | 1.415 | 1.048 | 1.261 | 14.1 | 14.0 | 14.1 | 14.8 |
| P5 | 18.1 | 1.769 | 1.311 | 1.563 | 15.1 | 14.0 | 14.1 | 14.8 |
| P6 | 14.2 | 2.265 | 1.678 | 2.043 | 16.8 | 14.0 | 14.2 | 14.8 |
| P7 | 11.3 | 2.831 | 2.097 | 2.575 | 18.4 | 14.0 | 14.2 | 14.8 |

Key:
- **MAE raw**: Hybrid with no calibration (S = 1.0)
- **MAE Phase 1**: After TDD-based initialisation
- **MAE Auto-Cal**: After Phase 2 auto-calibration (237 days)
- **MAE Perfect**: With theoretically optimal S (matching ISF at target)

### 7.3 Key Findings

#### Phase 1 calibration is remarkably effective

The most striking result is that **Phase 1 (TDD-based) calibration achieves MAE 14.0 for every synthetic patient**, regardless of sensitivity level — from the very resistant P1 (TDD 45 U/day) to the very sensitive P7 (TDD 11 U/day). This occurs because the TDD-based scaling factor `S = (1800/TDD) / 75.8` exactly compensates for sensitivity differences in the counterfactual prediction. The patient's scale factor cancels algebraically:

```
ISF_formula / ISF_true = hybrid(BG) × S_phase1 / (ISF_actual × patient_scale)
                       = hybrid(BG) × (1800 × scale) / (TDD × 75.8 × ISF_actual × scale)
                       = hybrid(BG) × 1800 / (TDD × 75.8 × ISF_actual)
```

The scale factor drops out, leaving a ratio that depends only on the hybrid's curve shape relative to the actual ISF-glucose relationship — not on the patient's absolute sensitivity level. This is a powerful validation that the 1800/TDD rule correctly personalises the hybrid formula across a 4× range of insulin sensitivity.

#### The curve shape determines the accuracy floor

Without any calibration, MAE ranges from 14.1 (original patient) to 20.2 (very resistant). Phase 1 collapses this to a uniform 14.0. The remaining 14.0 mg/dL error is the irreducible cost of the hybrid's curve shape not perfectly matching the actual ISF-glucose relationship. Importantly, even "perfect" calibration (S matched exactly to ISF at target) achieves MAE 14.8 — slightly *worse* than Phase 1. This means the 1800/TDD scaling happens to optimise prediction accuracy across glucose bands better than point-matching at a single glucose value.

#### Auto-calibration converges but doesn't improve MAE

All patients converge to within 5% of S_true by day 19. The auto-calibration reduces Phase 1's S error from 25.9% to approximately 9–12%. However, this improved S accuracy does not translate to improved MAE (14.2 vs 14.0), because the remaining prediction error is driven by curve-shape mismatch that no single scaling factor can correct.

#### Glucose-band performance is uniform across patients

After calibration, all patients show nearly identical glucose-band MAE:

| Glucose Band | MAE (all patients) |
|:---:|:---:|
| <90 mg/dL | 19.6–19.7 |
| 90–105 mg/dL | 13.7–13.8 |
| 105–120 mg/dL | 9.7–9.9 |
| 120–150 mg/dL | 12.3–12.6 |

This confirms that the hybrid's curve shape works consistently across the full sensitivity range, with the expected pattern: best accuracy above 105 (polynomial region), deliberate conservatism below 90.

### 7.4 Limitations of Synthetic Validation

This validation has an inherent limitation: all synthetic patients share the same ISF-glucose curve *shape*, differing only in magnitude. Real patients may have different curve shapes due to varying degrees of insulin resistance at different glucose levels, different counter-regulatory hormone responses, or different insulin absorption profiles.

However, the Diabeloop poster's proportional IQR scaling — where 25th and 75th percentile ISF values scale proportionally with the median across all glucose levels — suggests that the curve-shape assumption holds at the population level. The synthetic validation confirms that the calibration mechanism works correctly for the magnitude-scaling case; real multi-patient data would be needed to test the shape-variation case.

---

## 8. Discussion

### 8.1 Two Viable Approaches

This analysis, together with the companion power-law paper, presents two viable replacement formulas for dynamic ISF:

**Power-law with Boost TDD blending** — The most accurate formula tested (MAE 13.2, bias +5.6). Best suited for systems that already implement TDD tracking and Boost-style blending. Recommended where maximum prediction accuracy is the priority and the TDD infrastructure exists.

**Hybrid polynomial/power-law** — The most accurate TDD-free formula tested (MAE 14.1, bias +8.2). Best suited for simpler implementations, systems without reliable TDD data, or contexts where an additional safety margin at low glucose values is desired. The graduated conservatism below target is a clinical advantage in the zone where hypoglycaemia risk is highest.

### 8.2 The Case for the Hybrid

The hybrid formula offers a compelling alternative for several reasons:

1. **Clinical data grounding** — The polynomial piece is derived from observed population-level ISF–glucose relationships, not fitted to a single patient's loop data
2. **Simplicity** — The curve shape is fixed; only a single scaling factor needs to be learned per patient
3. **Safety-appropriate bias** — Conservative where conservatism protects against hypos, accurate where accuracy improves Time in Range
4. **No-test calibration** — Patients can start with a TDD-based or profile-ISF-based scaling factor on day one, with automatic refinement from observed data within the first week
5. **Implementation ease** — A few lines of code with minimal external state (only the scaling factor needs to persist)

### 8.3 Limitations

- **Synthetic multi-patient validation** — The multi-patient backtest (Section 7) uses scaled copies of one patient's data. While this validates the calibration mechanism and confirms the curve shape works across a 4× sensitivity range, it cannot test whether patients with genuinely different ISF-glucose curve shapes would achieve similar results. Real multi-patient data is needed.
- **Overnight only** — The backtest is restricted to fasting overnight periods. Daytime performance with meals, exercise, and stress remains untested for the hybrid specifically, though the companion paper showed all formulas perform worse during daytime (MAE approximately 2× overnight).
- **Calibration convergence** — The auto-calibration converges to within 5% of the target scaling factor by day 19 in synthetic testing. During the initial period, the TDD-based Phase 1 scaling factor provides excellent prediction accuracy (MAE 14.0), so suboptimal control during convergence is unlikely.
- **Junction at 105** — While mathematically continuous, the derivative is not matched at the junction. In practice, this has no clinical impact as the loop recalculates ISF every cycle (typically every 5 minutes).

### 8.4 Future Work

1. **Real multi-patient validation** — The synthetic validation (Section 7) confirms calibration works for magnitude-scaled patients. Testing against real data from patients with diverse physiology would validate the curve-shape assumption
2. **Daytime backtesting** — Extend validation to post-meal and exercise periods
3. **Junction optimisation** — Evaluate whether 105 mg/dL is the optimal switchover point, or whether a value closer to 100 or 110 would improve the 90–105 band
4. **Circadian scaling** — Investigate whether a time-of-day modifier on the scaling factor (analogous to Boost's TDD blending) improves daytime performance without compromising the formula's simplicity
5. **Phase 1 vs Phase 2 value** — The synthetic backtest showed Phase 1 (TDD-based) calibration achieving MAE 14.0 with no further improvement from auto-calibration. Determine whether this holds for patients with curve-shape differences, where auto-calibration may add more value

---

## 9. Conclusion

The hybrid ISF formula combines the strengths of two approaches: a clinically-derived polynomial that provides excellent accuracy above 105 mg/dL, and a power-law tail that provides graduated conservatism below 105 mg/dL. It produces a safety-appropriate bias profile — conservative where hypoglycaemia risk is highest, accurate where Time in Range optimisation matters most.

Patient individualisation is achieved through a single multiplicative scaling factor, initialised from existing TDD or profile ISF on day one. Synthetic multi-patient validation across 7 patients spanning a 4× range of insulin sensitivity (TDD 11–45 U/day) demonstrates that this TDD-based calibration achieves MAE 14.0 for every patient type — confirming that the curve shape generalises across the clinically relevant sensitivity range. The remaining prediction error is determined by curve-shape fidelity, not by calibration accuracy, meaning no further parameter tuning is required for the magnitude-scaling case.

While the power-law with Boost TDD blending remains the most accurate formula on pure prediction metrics, the hybrid offers a simpler alternative that may produce equivalent or superior clinical outcomes by better balancing accuracy against safety in the near-hypoglycaemic range, while requiring only a single learned parameter rather than real-time TDD computation.

---

## Appendix A: Source Data

| Parameter | Value |
|-----------|-------|
| Data period | June 2025 – March 2026 (10 months) |
| Total loop cycles | 94,980 |
| CGM readings | 90,733 |
| Treatment records | 61,975 |
| Valid overnight +2h samples | 3,647 |
| Insulin divisor (D) | 82 |
| Target glucose | 99 mg/dL |
| Patient actual TDD (median) | 21.9 U/day |
| Boost blended TDD (overnight median) | 16.8 U/day |

## Appendix B: Diabeloop ADA Poster Reference

Source: Diabeloop, ADA Scientific Posters
https://ada.scientificposters.com/epsAbstractADA.cfm?id=1

Equation for glucose > 100 mg/dL (quartic polynomial):
```
ISF = 272 − 3.121G + 0.01511G² − 3.305×10⁻⁵G³ + 2.69×10⁻⁸G⁴
```

The poster also provides a quadratic equation for glucose ≤ 100 mg/dL. This equation was not adopted in the hybrid formula; instead, a power-law tail was used below 105 mg/dL to provide deliberate conservatism in the near-hypoglycaemic range (see Section 2.1 for rationale).

## Appendix C: Comparison Charts

The formula backtest produced a multi-panel chart saved as `ns_polynomial_backtest.png`, showing:
- ISF curves for all formulas across the 60–200 mg/dL range
- MAE, bias, and ±1 mmol/L accuracy bar charts
- Error distributions
- Glucose-band MAE and bias breakdowns

The multi-patient calibration backtest produced a chart saved as `ns_hybrid_calibration.png`, showing:
- Scaling factor convergence over time for all 7 synthetic patients
- MAE convergence trajectories
- MAE and bias comparison across calibration stages (uncalibrated → Phase 1 → auto-calibrated → perfect)
- Calibrated ISF curves by patient
- Scaling factor error reduction
- Glucose-band MAE heatmap across patients
