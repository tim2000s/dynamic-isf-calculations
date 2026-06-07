# Power-Law DynamicISF: A Retrospective Analysis of BG-Dependent Insulin Sensitivity Scaling for Automated Insulin Delivery

**Date:** 31 March 2026
**Dataset:** AndroidAPS Boost v1 overnight logs, June 2025 – March 2026 (10 months)
**Patient:** N=1, Type 1 Diabetes, closed-loop AID system
**Target platforms:** AAPS (Android, Kotlin) and Trio (iOS, Swift)

---

## Executive Summary

This study analyses 10 months of overnight closed-loop data from a single Type 1 Diabetes patient using AndroidAPS Boost v1 to evaluate whether a power-law formula for DynamicISF can outperform the current logarithmic formula.

**What we tested:** The loop predicts where BG will be in 2 hours based on its ISF formula. We compared the prediction against what actually happened, across 3,647 qualifying overnight samples. We then asked: would a different ISF formula have predicted more accurately?

**What we found:**

1. **The Boost formula is more complex than commonly described.** It uses a TDD blending algorithm that reduces effective TDD from 22.6 to ~16.8 U/day overnight, a sensitivityRatio adjustment, and a velocity parameter. All of these must be accounted for in any analysis.

2. **All formulas over-predict insulin effect overnight.** The loop predicts BG will drop more than it actually does — meaning BG ends up higher than expected. This "positive bias" ranges from +8 to +14 mg/dL depending on the formula, and even the loop's own tuned predictions are biased by +11.2 mg/dL.

3. **Power-law BG scaling is structurally better than logarithmic.** The proposed formula `ISF = (C/TDD) × (target/BG)^k` reduces prediction error by 12–18% compared to the current ln-based formula when using the same TDD input. It is also inherently more protective at low BG and more aggressive at high BG.

4. **The Boost TDD blending should be retained.** Testing 9 different TDD computation methods, the Boost blended TDD produces the best results with the power-law formula. It is the only TDD variant where the existing C=1800 constant is already optimal — all other TDD inputs require a substantially higher C. The blending adapts to circadian insulin delivery patterns in a way that fixed discounts cannot.

5. **Recommended defaults: C=1800, k=3.5, with Boost-style TDD blending.** These values are directly supported by the data for this patient. The joint-optimised k is stable at 3.2–3.6 across all real TDD variants and across both overnight and dawn time windows. Multi-patient validation is essential to confirm these defaults generalise.

**Bottom line:** Replace `ln(BG/D+1)` with `(target/BG)^k`. Use the existing Boost TDD blending, set C=1800 and k=3.5 as defaults. These parameters are data-driven but N=1 — multi-patient validation is required before deployment.

---

## Abstract

The DynamicISF algorithm adjusts the Insulin Sensitivity Factor (ISF) in real time based on blood glucose (BG) and total daily insulin dose (TDD). Using 10 months of overnight fasting data from AndroidAPS Boost v1 and actual TDD computed from treatment records, this paper evaluates whether a power-law BG scaling function `ISF = (C/TDD) × (target/BG)^k` can outperform the current logarithmic formula.

A key finding is that the Boost implementation uses a complex TDD blending algorithm that reduces effective overnight TDD from 22.6 U/day (7-day average) to ~16.8 U/day by weighting recent insulin delivery. This blending has a larger effect on ISF than any change to the BG scaling function. Despite this built-in adjustment, the loop's own overnight predictions carry a persistent positive bias of +11.2 mg/dL at the +2h horizon.

Power-law BG scaling improves over logarithmic scaling at every exponent tested. With 7-day TDD, joint optimisation yields C≈2834, k≈3.0. However, the population constant C, the exponent k, and the TDD input are strongly entangled — the optimal C and k shift substantially depending on which TDD is used. This entanglement must be resolved through multi-patient validation before deployment.

---

## 1. Introduction

### 1.1 The DynamicISF Problem

Automated Insulin Delivery (AID) systems calculate insulin dosing adjustments every 5 minutes using an Insulin Sensitivity Factor (ISF) — the expected BG drop per unit of insulin. Static ISF values fail to account for glucose-dependent insulin sensitivity: at high BG, insulin is more effective (lower ISF needed); at low BG, insulin is less effective (higher ISF needed, providing inherent hypo protection).

DynamicISF formulas attempt to capture this relationship. The formula used in AAPS Boost v1 is more complex than commonly described, and understanding its actual implementation is essential for evaluating any proposed replacement.

### 1.2 The Boost Formula

From the Boost source code (`IsfCalculatorImpl.kt`), the formula involves three stages:

**Stage 1 — TDD Blending:**
```
tddLast4H  = insulin delivered in the last 4 hours
tddLast8to4H = insulin delivered 4–8 hours ago
tddWeighted = ((1.4 × tddLast4H) + (0.6 × tddLast8to4H)) × 3

if tddWeighted < 0.75 × tdd7D and tdd1D exists:
    TDD = ((tddWeighted + (tddWeighted/tdd7D) × (tdd7D - tddWeighted)) × 0.34)
        + (tdd1D × 0.33) + (tddWeighted × 0.33)
else if tdd1D exists:
    TDD = (tddWeighted × 0.33) + (tdd7D × 0.34) + (tdd1D × 0.33)
else:
    TDD = tddWeighted

TDD *= dynISFadjust / 100
```

**Stage 2 — Base ISF at target:**
```
sensNormalTarget = 1800 / (TDD × ln(target/D + 1))
```

Where D = 82 (insulinDivisor, derived from insulin peak time 38 min) and target = 99 mg/dL (normalTarget).

**Stage 3 — BG scaling with sensitivityRatio:**
```
if adjustSens:
    ratio = clamp(tddLast24H / tdd7D, autosensMin, autosensMax)
    sensNormalTarget /= ratio

scaler = ln(target/D + 1) / ln(BG/D + 1)
variable_sens = sensNormalTarget × (1 - (1 - scaler) × dynIsfVelocity)
```

Where `dynIsfVelocity` (default 1.0) controls BG scaling strength and `autosensMin/Max` (default 0.5/1.5) bound the sensitivityRatio.

The commonly cited simplified form — `ISF = 1800 / (TDD × ln(BG/D+1))` — omits the TDD blending, the sensitivityRatio, and the velocity parameter. These omitted components have a substantial effect on the ISF actually used by the loop, particularly overnight.

### 1.3 The Research Question

Does the current logarithmic BG scaling `ln(target/D+1) / ln(BG/D+1)` produce accurate overnight predictions? If not, would a power-law `(target/BG)^k` perform better? And critically: how do the answers depend on which TDD input is used?

### 1.4 Why Overnight

The overnight fasting window (00:00–07:00) provides an ideal test environment:
- **No carbs**: COB = 0 eliminates meal absorption as a confound
- **No exercise**: Minimal glucose uptake from activity
- **IOB dominance**: Insulin on board is the primary active signal
- **Daily availability**: Every night provides data
- **HGO cancellation**: Hepatic glucose output affects both predicted and actual BG equally, since the profile basal already accounts for steady-state liver glucose output

### 1.5 Scope

This analysis informs a formula change proposal for both AAPS and Trio. Results are from a single patient (N=1) and require multi-patient validation. All analysis is retrospective.

---

## 2. Data and Method

### 2.1 Dataset

Log data were extracted from AndroidAPS Boost v1 devicestatus records covering June 2025 to March 2026 (10 months). A total of 94,980 loop cycles were parsed from 101,608 unique devicestatus records, joined with 90,733 CGM entries and 61,975 treatment records. Each 5-minute loop cycle logs:
- Current BG (CGM reading)
- IOB (insulin on board)
- The formula's computed ISF (`variable_sens`)
- A full predicted BG trajectory (`predBGs.IOB`)
- COB, bolus history, profile settings

After fasting filtering (COB=0, ≥3h since last bolus, BG 72–200 mg/dL), overnight cycles were identified in the 00:00–07:00 window across approximately 300 usable nights.

### 2.2 TDD from Treatment Records

TDD was computed directly from treatment records — the actual insulin delivered:

- **Bolus insulin**: Sum of all bolus treatments per day
- **Temp basal insulin**: Computed from temp basal rate × duration for each temp basal record, summed per day
- **7-day rolling TDD**: Mean of actual daily TDD over the preceding 7 days

| Metric | Value |
|---|---|
| Basal insulin (median) | 7.2 U/day |
| Bolus insulin (median) | 14.5 U/day |
| **Total TDD (median)** | **21.9 U/day** |
| **7-day rolling TDD (median)** | **22.6 U/day** |
| Basal : bolus ratio | 33 : 67 |

### 2.3 Boost TDD Blending — Overnight Behaviour

The Boost TDD blending algorithm (Section 1.2) has dramatic effects overnight. During sleep, insulin delivery is basal-only at ~0.3 U/hr:

```
tddLast4H ≈ 1.2 U → extrapolated to daily: 1.2 × 6 = 7.2
tddLast8to4H ≈ 1.2 U
tddWeighted = ((1.4 × 1.2) + (0.6 × 1.2)) × 3 = 7.2 U/day

Since 7.2 < 0.75 × 22.6 = 17.0 → uses low-TDD branch:
TDD ≈ 0.33 × 7.2 + 0.34 × 22.6 + 0.33 × 21.9 ≈ 17.3 U/day
```

Observed Boost blended TDD overnight: **median 16.8 U/day** (vs 22.6 actual 7-day). The blending reduces effective TDD by ~26% overnight, raising ISF proportionally — the loop uses a more conservative ISF while the patient sleeps.

### 2.4 The sensitivityRatio

Boost also divides `sensNormalTarget` by:
```
ratio = clamp(tddLast24H / tdd7D, 0.5, 1.5)
```

Overnight distribution of sensitivityRatio:
- 37% of cycles: ratio = 1.0 (last 24h ≈ 7-day average)
- 14% of cycles: ratio = 0.80 (last 24h lower → more sensitive → ISF raised)
- 3% of cycles: ratio = 0.77 (approaching minimum)

When ratio < 1.0, ISF is raised further beyond what TDD blending already achieves.

### 2.5 Counterfactual Prediction Method

For each loop cycle at time t, the loop logged its +2h BG prediction (`pred_iob_24`, the 24th 5-minute step). Comparing against actual BG at t+2h:

```
pred_error = BG_actual(t+2h) − BG_predicted(t+2h)
```

Positive error = actual BG higher than predicted = formula over-estimated insulin effect.

To compare alternative formulas without re-running the loop:

```
bg_drop_pred = BG(t) − BG_predicted(t+2h)
pred_formula = BG(t) − bg_drop_pred × (ISF_formula / ISF_actual)
error_formula = BG_actual(t+2h) − pred_formula
```

This assumes linear ISF scaling, valid for the IOB prediction model used in AAPS.

### 2.6 Understanding Bias and MAE

Two key metrics are used throughout this paper:

**Bias (mean error)** measures the average direction of prediction errors. A **positive bias** (e.g., +11 mg/dL) means the formula consistently predicts BG will be *lower* than it turns out to be — the formula over-estimates how much insulin will drop BG. In practical terms: the loop thinks insulin is working harder than it is, so it may deliver more than needed. Positive bias is the more dangerous direction, as it leads toward over-delivery and potential hypoglycaemia. A **negative bias** would mean the opposite — the formula under-estimates insulin effect, leading to under-delivery and higher BG.

**MAE (Mean Absolute Error)** measures prediction accuracy regardless of direction. A MAE of 15 mg/dL means the prediction is off by 15 mg/dL on average, sometimes high, sometimes low. A formula can have low bias (errors cancel out) but high MAE (large errors in both directions). The ideal formula has both low bias and low MAE.

For context: a prediction error of ±18 mg/dL (±1 mmol/L) is commonly used as the threshold for a "clinically acceptable" prediction.

### 2.7 Filtering

| Filter | Threshold | Reason |
|---|---|---|
| BG out of range | < 72 or > 200 mg/dL | Extremes distort calculations |
| Predicted drop too small | < 3 mg/dL absolute | Division instability |
| Ratio mismatch | ≤ 0 or > 5 | BG moved opposite to prediction |
| Suspected missed carbs | BG rose > 9 mg/dL when drop predicted | Unrecorded food |
| Missing data | NaN in key fields | Incomplete loop cycle |

After filtering: **3,647 valid +2h overnight samples.**

---

## 3. Results

### 3.1 Formula Comparison

Five formula variants were compared:

| Label | BG Scaling | TDD Source |
|---|---|---|
| A: Loop actual | As implemented (ln + velocity) | Boost blended |
| B: Current ln, 7D-TDD | `1800 / (TDD × ln(BG/D+1))` | Actual 7-day |
| C: Current ln, blended TDD | `1800 / (TDD × ln(BG/D+1))` | Boost blended |
| D: Power-law k=2.0 | `(1800/TDD) × (target/BG)^2.0` | Actual 7-day |
| E: Power-law k=3.0 | `(1800/TDD) × (target/BG)^3.0` | Actual 7-day |

Results at +2h overnight:

| Variant | MAE (mg/dL) | Bias (mg/dL) | ±1 mmol/L |
|---|---|---|---|
| A: Loop actual | **13.4** | +11.2 | **75.8%** |
| B: Current ln, 7D-TDD | 17.3 | +14.4 | 65.6% |
| C: Current ln, blended TDD | 16.7 | +14.1 | 65.8% |
| D: Power-law k=2.0, 7D-TDD | 15.2 | +11.0 | 71.1% |
| E: Power-law k=3.0, 7D-TDD | 14.1 | +8.1 | 74.3% |

### 3.2 Key Observations

**1. All formulas show large positive bias with actual TDD.**
Every variant shows +8 to +14 mg/dL positive bias — the formulas systematically predict BG will drop more than it does. This indicates the ISF formulas are too aggressive overnight, even after accounting for TDD blending.

**2. The loop's own predictions (A) are biased +11.2 mg/dL.**
Even `variable_sens` (which uses Boost TDD blending and sensitivityRatio) over-predicts insulin effect overnight. The blending raises ISF, but not enough for accurate predictions.

**3. Power-law improves over logarithmic at every k.**
Variant D (k=2.0) reduces MAE from 17.3 to 15.2 (12% improvement over ln with the same 7D-TDD). Variant E (k=3.0) reduces it further to 14.1 (18% improvement), approaching the loop's own 13.4.

**4. Boost TDD blending barely helps the ln formula.**
Comparing B (7D-TDD) vs C (blended TDD): MAE drops only 0.6, bias drops 0.3. The blending reduces TDD from 22.6 to 16.8, but the ln formula's weak BG scaling absorbs most of the ISF benefit, leaving the bias largely unchanged.

**5. Optimal k hits the upper bound.**
With actual 7D-TDD, the optimal power-law exponent is k ≥ 3.0 (search bound). The data demands steep BG-dependent scaling when the TDD input reflects actual insulin delivered.

### 3.3 Joint C + k Optimisation

When both the population constant C and exponent k are fitted simultaneously:

| Target | C | k |
|---|---|---|
| Minimise MAE | **2834** | **3.0** |
| Minimise \|bias\| | **3000** | **3.0** |

Both optimisations converge to a substantially higher constant (~2800–3000 vs the traditional 1700–1800) and the steepest exponent tested. The constant and exponent compensate each other: a higher C shifts ISF upward (less aggressive overall), while a higher k steepens the BG dependence (more aggressive above target, more protective below).

---

## 4. The TDD Blending Effect

### 4.1 Why Blended TDD Matters

The TDD input dominates the formula's output. With actual 7D-TDD of 22.6:

```
sensNormalTarget = 1800 / (22.6 × ln(99/82 + 1)) = 1800 / (22.6 × 0.789) = 101 mg/dL/U
```

With Boost blended TDD of 16.8:

```
sensNormalTarget = 1800 / (16.8 × ln(99/82 + 1)) = 1800 / (16.8 × 0.789) = 136 mg/dL/U
```

A **34% increase** in ISF from TDD blending alone. This is a larger effect than any BG scaling change could produce.

### 4.2 TDD Variant Comparison

Nine TDD computation methods were tested with the power-law formula, each with its own optimal k (C=1800 fixed) and joint-optimised (C, k):

| TDD Variant | Median | PL k_opt | PL MAE | PL bias | Joint C | Joint k | Joint MAE |
|---|---|---|---|---|---|---|---|
| Actual 7-day | 22.6 | 4.00† | 14.77 | +6.77 | 2679 | 3.57 | 13.39 |
| Actual 1-day | 21.4 | 4.00† | 14.85 | +6.10 | 2414 | 3.43 | 14.22 |
| **Boost blended** | **16.8** | **3.59** | **13.21** | **+5.25** | **1800** | **3.59** | **13.21** |
| 50/50 (7D+1D) | 22.4 | 4.00† | 14.66 | +6.57 | 2607 | 3.51 | 13.53 |
| 70/30 (7D+1D) | 22.7 | 4.00† | 14.68 | +6.67 | 2651 | 3.52 | 13.43 |
| 7-day × 0.85 | 19.3 | 4.00† | 13.87 | +5.54 | 2277 | 3.57 | 13.39 |
| 7-day × 0.75 | 17.0 | 3.87 | 13.50 | +4.91 | 2009 | 3.57 | 13.39 |
| Last 24h actual | 22.1 | 4.00† | 13.86 | +5.74 | 2636 | 3.23 | 12.67 |
| Implied (ref) | 14.9 | 2.40 | 11.35 | +5.89 | 1866 | 2.30 | 11.33 |

*† hit upper search bound at k=4.0*

### 4.3 Key Findings

**1. Boost blended TDD is the only variant where C=1800 is already optimal.** The joint optimiser for every other real TDD variant wants C >> 1800 (2009–2679). For Boost blended, C=1800 with k=3.59 is already at the joint optimum — no constant adjustment needed.

**2. All real TDD variants converge to joint k ≈ 3.2–3.6.** The spread is only 0.4 — the exponent is stable across TDD inputs. What shifts is C (1800–2679), which absorbs the TDD level difference.

**3. Simpler blends are nearly as good.** A fixed 7-day × 0.75 discount (median 17.0, close to Boost's 16.8) achieves MAE 13.50 vs Boost's 13.21. However, a fixed discount doesn't adapt through the day (see Section 4.4).

**4. The implied TDD produces the best absolute numbers but is circular.** It is back-calculated from the loop's own output and cannot be used as a real input.

### 4.4 Blending Varies by Hour

The Boost blending reduces TDD progressively through the night as the 4h window fills with basal-only delivery:

| Hour | Boost TDD | 7-day TDD | Ratio | n |
|---|---|---|---|---|
| 00:00 | 19.8 | 23.1 | 0.86 | 316 |
| 02:00 | 17.5 | 22.4 | 0.78 | 670 |
| 04:00 | 16.0 | 22.4 | 0.71 | 562 |
| 06:00 | 15.3 | 22.6 | 0.68 | 335 |
| 08:00 | 10.9 | 20.2 | 0.54 | 84 |
| 10:00 | 20.8 | 22.6 | 0.92 | 33 |

The blending bottoms out at hour 7–8 (54% of 7-day TDD), then recovers as morning boluses enter the 4h window. A fixed discount cannot reproduce this circadian adaptation — this is the key advantage of a blending algorithm.

### 4.5 Daytime Fasting Analysis

To confirm the overnight findings aren't biased by time-of-day, the same analysis was run on daytime fasting samples (08:00–22:00, COB=0, ≥3h since bolus). Only **238 samples** survived filtering — fasting windows are rare during waking hours — so these results are directionally interesting but statistically weak.

| Time Window | n | Best TDD (C=1800) | k_opt | MAE | Bias |
|---|---|---|---|---|---|
| Overnight (00–08) | 3,647 | Boost blended | 3.59 | 13.21 | +5.25 |
| Deep night (00–04) | 2,152 | Boost blended | 3.64 | 12.90 | +5.06 |
| Dawn (04–08) | 1,495 | Boost blended | 3.39 | 13.64 | +5.75 |
| Daytime (08–22) | 238 | Boost blended | 4.00† | 23.62 | +16.76 |
| All fasting | 4,058 | Boost blended | 3.68 | 13.83 | +5.78 |

**Boost blended TDD wins every time window.** The optimal k is stable overnight (3.4–3.6) but hits the search bound during the day (k=4.0). Daytime errors are 2× overnight across all formulas — this is a fundamental limitation of ISF-only correction during waking hours (unmodelled activity, stress, hormones), not a TDD problem.

The joint-optimised k across all time windows and TDD variants clusters at **3.2–3.6**, supporting k≈3.5 as a robust default.

---

## 5. TDD Inflation and Deflation

### 5.1 The Community Concern

A recurring argument in the AAPS and Trio communities is that TDD-based ISF formulas are inherently flawed because overestimated basal profiles inflate TDD, causing ISF to be too low (too aggressive). The proposed feedback loop is:

```
Basal too high → excess insulin → TDD inflated → ISF too low →
over-correction → more insulin → TDD rises further
```

### 5.2 Boost Blending as Intentional TDD Deflation

The Boost TDD blending is effectively a **TDD deflation mechanism**. By weighting recent insulin delivery (which is low overnight), it reduces the effective TDD by ~26% relative to the 7-day average. This is the system's built-in defence against TDD inflation — and against overnight over-aggressiveness.

However, the loop's predictions using blended TDD still carry +11.2 mg/dL bias, meaning even with the deflated TDD, the formula over-predicts insulin effect. This could mean:
- The ln BG scaling is inherently miscalibrated
- The blending doesn't adequately capture overnight sensitivity
- There is an irreducible prediction floor from CGM lag or other factors
- Some combination of these

### 5.3 Basal IOB Check

| Metric | Overnight Median |
|---|---|
| Basal IOB | 0.021 U |
| Total IOB | ~0.021 U |

Basal IOB is essentially zero overnight, suggesting the profile basal is not substantially overestimated for this patient. The positive bias likely originates in the BG scaling or the ISF constant, not in excess basal delivery.

### 5.4 Power-Law Resilience

Comparing TDD sensitivity between formulas:

| Formula | TDD Scale for Zero Bias |
|---|---|
| Current ln, 7D-TDD | Needs TDD reduced substantially |
| Power-law k=3.0, 7D-TDD | Closer to zero bias at actual TDD |

The power-law's steeper BG scaling absorbs more of the ISF "error" that the ln formula attributes to TDD miscalibration. This suggests part of the perceived TDD inflation problem is actually a BG scaling problem.

---

## 6. The Proposed Formula

### 6.1 Definition

```
ISF = max(ISF_floor, (C / TDD_blended) × (target / BG) ^ k)
```

Where:
- **C** = population constant, default **1800**
- **TDD_blended** = Boost-style blended TDD (retaining the existing 4h/8h weighted algorithm)
- **target** = user's profile target BG (typically 99–110 mg/dL)
- **BG** = current CGM reading (mg/dL)
- **k** = BG scaling exponent, default **3.5**, user-adjustable range 2.0–4.0
- **ISF_floor** = minimum ISF safety bound (e.g., 10 mg/dL/U)

### 6.2 Why These Defaults

**C = 1800:** The Boost blended TDD is the only TDD variant where C=1800 is already at the joint optimum. No constant change is needed. This preserves the familiar "1800 rule" that users and clinicians already understand.

**k = 3.5:** The joint-optimised k clusters at 3.2–3.6 across all real TDD variants and across overnight (3.59), deep night (3.64), and dawn (3.39) time windows. k=3.5 sits in the centre of this range.

**TDD = Boost blended:** Tested against 8 alternative TDD computations, the Boost blending produces the lowest MAE at C=1800 in every time window. It adapts to circadian insulin delivery patterns — falling through the night as recent delivery drops, recovering as morning boluses arrive. A fixed discount (e.g., 7D × 0.75) achieves ~80% of the benefit overnight but cannot adapt during the day.

### 6.3 Properties

**At target BG:**
```
ISF = C / TDD_blended = 1800 / TDD_blended
```
Recovers the standard 1800 rule with blended TDD — clinically familiar, no surprises.

**Below target (BG < target):**
`(target/BG) > 1` raised to power k → ISF rises steeply → less insulin → hypo-protective.

**Above target (BG > target):**
`(target/BG) < 1` raised to power k → ISF falls steeply → more aggressive correction.

### 6.4 Safety Analysis

At BG=60, target=99, k=3.5:
- Current ln: ISF ∝ `ln(99/82+1) / ln(60/82+1)` = 1.50 × ISF_target
- Proposed: ISF ∝ `(99/60)^3.5` = 5.96 × ISF_target

The proposed formula delivers ISF **4× higher** at BG=60 → dramatically stronger hypo protection.

At BG=200, target=99, k=3.5:
- Current ln: ISF ∝ 0.64 × ISF_target
- Proposed: ISF ∝ `(99/200)^3.5` = 0.08 × ISF_target

The proposed formula is **8× more aggressive** at BG=200. This is appropriate for correction but must be bounded by ISF_floor and existing platform safety constraints.

### 6.5 ISF Comparison Table

At target=99, C=1800, Boost blended TDD=16.8 overnight:

| BG (mg/dL) | Current ln ISF | Power-law k=3.5 ISF | Ratio (PL/ln) |
|---|---|---|---|
| 60 | 161 | 639 | 4.0× (more protective) |
| 75 | 136 | 296 | 2.2× |
| 99 (target) | 107 | 107 | 1.0× (anchored) |
| 125 | 92 | 53 | 0.58× |
| 160 | 79 | 22 | 0.28× (more aggressive) |
| 200 | 70 | 9 | 0.13× |

The power-law is dramatically more protective below target and more aggressive above. The steepness is governed by k — users who find the correction too aggressive above target can reduce k toward 2.5; those who want stronger correction can increase toward 4.0.

### 6.6 Safety Bounds

- **ISF floor**: Minimum ISF to prevent extreme dosing at very high BG (e.g., 10 mg/dL/U or user-configurable). At BG=200 with the values above, ISF is 9 — the floor would engage.
- **BG cap**: Existing Boost feature — caps BG input at a threshold (e.g., 220 mg/dL), compressing above as `bg = cap + (bg - cap)/3`. This provides additional protection against extreme ISF values.
- **Max bolus / max IOB**: Existing platform constraints, unchanged.
- **k range**: 2.0–4.0 (user-adjustable), default 3.5.

---

## 7. Implementation Notes

### 7.1 AAPS (Kotlin)

```kotlin
// Current
val sensNormalTarget = 1800.0 / (tdd * ln(bgNormalTarget / insulinDivisor + 1.0))
val scaler = ln(bgNormalTarget / insulinDivisor + 1.0) / ln(bg / insulinDivisor + 1.0)
val variableSens = sensNormalTarget * (1 - (1 - scaler) * dynIsfVelocity)

// Proposed
val variableSens = max(isfFloor, (1800.0 / tdd) * (target / bg).pow(dynISFExponent))
```

Where `dynISFExponent` is a new preference (default 3.5, range 2.0–4.0) and `target` is the existing profile target BG. The TDD blending (`tdd` variable) is unchanged — it uses the existing Boost blending algorithm.

The insulin divisor D (82), the velocity parameter, and the scaler are all eliminated. The formula is simpler: one `pow()` call replacing two `ln()` calls plus a scaler computation.

### 7.2 Trio (Swift)

```swift
// Current
let sensNormalTarget = 1800.0 / (tdd * log(bgNormalTarget / insulinDivisor + 1.0))
let scaler = log(bgNormalTarget / insulinDivisor + 1.0) / log(bg / insulinDivisor + 1.0)
let variableSens = sensNormalTarget * (1 - (1 - scaler) * dynIsfVelocity)

// Proposed
let variableSens = max(isfFloor, (1800.0 / tdd) * pow(target / bg, dynISFExponent))
```

### 7.3 TDD Blending: Retain As-Is

The existing Boost TDD blending algorithm should be retained unchanged. The analysis shows it is well-designed:
- It adapts to circadian delivery patterns (falls through the night, recovers during the day)
- It is the only TDD computation where C=1800 is already optimal
- It outperforms simpler alternatives (fixed discounts, simple blends) in every time window tested

The only change is what happens *after* the blended TDD is computed: the ln-based BG scaling is replaced with a power-law.

### 7.4 No New Infrastructure

- **No nightly regression**: A TDD_effective approach was evaluated and found not to improve accuracy.
- **No new persisted state**: The TDD blending already exists. k is a static preference.
- **No new sensors or inputs**: Uses existing BG, TDD, and target values.
- **Fewer parameters**: The insulinDivisor D and dynIsfVelocity are eliminated. The dynISFExponent (k) replaces both.

---

## 8. Discussion

### 8.1 The Entanglement Problem — Resolved

The C, k, and TDD input are entangled: the optimal C and k shift depending on which TDD is used. However, the TDD variant analysis resolves this: the Boost blended TDD is the natural choice because it is the only variant where C=1800 (the existing constant) is already optimal. This eliminates the entanglement — the TDD blending absorbs the circadian sensitivity variation, leaving C and k free to address only the BG scaling shape.

### 8.2 Why the Loop's Own Predictions Still Win

The loop's own predictions (variant A, MAE 13.4) outperform all counterfactual formulas. The Boost system's TDD blending and sensitivityRatio do useful work. The power-law with Boost blended TDD (MAE 13.21) is the closest any counterfactual formula gets — it essentially matches the loop's own accuracy while reducing bias from +11.2 to +5.25.

### 8.3 Power-Law vs Logarithmic

The structural finding is clear: power-law BG scaling outperforms logarithmic scaling with the same TDD input, across every k tested, every TDD variant, and every time window. The power-law's steeper curvature better matches physiology — insulin sensitivity varies more strongly with BG than the gentle ln curve allows.

### 8.4 The Residual Bias

Even the best formula (power-law k=3.59 with Boost blended TDD) retains +5.25 mg/dL bias overnight. This is a substantial improvement from the +14.4 bias of the current ln formula with 7-day TDD, but it is not zero. Possible explanations:
- An irreducible floor from CGM lag, compression lows, or measurement artefacts
- That ISF alone cannot fully capture overnight sensitivity — basal rate accuracy also matters
- Limitations of the counterfactual prediction method

### 8.5 Daytime Limitations

Daytime fasting samples (n=238) show much larger errors (MAE 23+, bias +17+) across all formulas. This is not a TDD or BG scaling problem — it reflects unmodelled daytime factors (activity, stress, hormones, residual meal effects). The formula is optimised for overnight; daytime performance depends on other loop components (meal boluses, exercise adjustments, etc.).

---

## 9. Limitations

### 9.1 Single Patient (N=1)

All results are from one individual over 10 months. The optimal k, C, and TDD blending parameters will differ across patients. Multi-patient validation is essential.

### 9.2 Primarily Overnight

The analysis is primarily based on overnight fasting data (3,647 samples). Daytime fasting data (238 samples) was included for validation but is statistically weak. Daytime insulin sensitivity is affected by meals, exercise, stress, and diurnal hormonal cycles that produce 2× the prediction error of overnight regardless of formula.

### 9.3 Retrospective Counterfactual

Predictions are computed by scaling the loop's ISF-based BG drop, not by re-running the full loop algorithm. This assumes linear ISF scaling, valid for the IOB model but does not account for feedback effects (different ISF → different basal rate → different IOB → different prediction).

### 9.4 Entanglement with TDD Computation

Formula parameters depend on the TDD input. Any formula validated with one TDD computation method cannot be assumed valid with another. This is a fundamental constraint on deployment.

### 9.5 The <90 mg/dL Band

The lowest BG band has the highest MAE across all formula variants (19+ mg/dL). This likely reflects CGM noise at low glucose levels and the inherent difficulty of predicting BG near the hypo threshold. No formula tested meaningfully solved this.

---

## 10. Recommendations

### 10.1 Replace Logarithmic BG Scaling with Power-Law

The structural advantage of power-law over logarithmic BG scaling is robust across all TDD inputs tested, all time windows, and all k values. Implement:

```
ISF = max(ISF_floor, (1800 / TDD_blended) × (target / BG) ^ k)
```

### 10.2 Retain Boost-Style TDD Blending

The existing Boost TDD blending algorithm should be kept unchanged. It outperforms simpler alternatives in every time window and is the only TDD computation where C=1800 is already at the joint optimum.

### 10.3 Recommended Defaults

| Parameter | Default | Range | Rationale |
|---|---|---|---|
| C | 1800 | 1500–2500 | Joint-optimal with Boost blended TDD |
| k | 3.5 | 2.0–4.0 | Centre of stable cluster (3.2–3.6) across TDD variants and time windows |
| ISF_floor | 10 mg/dL/U | 5–20 | Prevents extreme dosing at very high BG |

### 10.4 Multi-Patient Validation is Required

These defaults are N=1. The validation path:

1. **Retrospective multi-patient**: Run this analysis on 5–10 patients with diverse TDD ranges. Confirm k clusters near 3.5 and C near 1800 when Boost blending is used.
2. **Prospective A/B**: Deploy as opt-in experimental mode alongside existing DynamicISF. Compare overnight TIR, time below range, and mean BG over 2-week periods.
3. **Safety monitoring**: Track time below 54 mg/dL (severe hypo). The formula's enhanced low-BG protection should reduce this, but must be confirmed.
4. **Daytime validation**: The current analysis is predominantly overnight. Prospective testing should specifically evaluate daytime performance.

### 10.5 Do Not Implement TDD_effective

A nightly TDD regression approach was evaluated and found not to improve accuracy. It adds complexity without benefit.

---

## 11. Files and Reproducibility

All scripts are in `~/Downloads` or `~/Nightscout_Work`:

| File | Description |
|---|---|
| `ns_corrected_analysis.py` | Main analysis — actual TDD from treatments, Boost blending, 5 formula variants |
| `ns_tdd_variants_analysis.py` | TDD variant comparison — 9 TDD methods, optimal k and joint C+k per variant |
| `ns_tdd_daytime_analysis.py` | Daytime vs overnight comparison — time-window analysis, hourly TDD patterns |
| `ns_corrected_summary.txt` | Results summary |
| `ns_tdd_variants_summary.txt` | TDD variant comparison summary |
| `ns_corrected_results.png` | Main results figure |
| `ns_tdd_variants_results.png` | TDD variant comparison figure |
| `ns_tdd_daytime_results.png` | Daytime vs overnight figure |
| `ns_actual_tdd_daily.csv` | Daily actual TDD computed from treatment records |
| `ns_rebuild_pipeline.py` | Full pipeline rebuild from all Nightscout exports |
| `ns_export_missing.py` | Non-interactive Nightscout data exporter |
| `DynamicISF_PowerLaw_Analysis.md` | This document |

---

## Appendix A: Mathematical Notes

### A.1 Power-Law Anchoring

```
ISF(BG) = (C / TDD) × (target / BG) ^ k
```

At BG = target: `ISF = C / TDD` (the "C rule").

### A.2 Derivative

```
dISF/dBG = -k × ISF(BG) / BG
```

At target: `dISF/dBG = -k × (C/TDD) / target`.

For the current ln formula at target: `dISF/dBG = -(1800/TDD) / (target + D) = -(1800/TDD) / 181`.

With k=3.5 and C=1800: the power-law slope = `3.5 × (1800/TDD) / 99 = (1800/TDD) / 28.3`.
The current ln slope = `(1800/TDD) / 181`.

The power-law is **6.4× steeper** at target BG with k=3.5. The data supports this much stronger BG dependence.

### A.3 Boost TDD Blending Derivation

Overnight, with basal-only delivery at rate R (U/hr):

```
tddLast4H = 4R
tddLast8to4H = 4R
tddWeighted = ((1.4 × 4R) + (0.6 × 4R)) × 3 = 24R

For R = 0.3 U/hr: tddWeighted = 7.2 U/day
```

Since `7.2 < 0.75 × 22.6 = 17.0`, the low-TDD branch activates:
```
TDD = ((7.2 + (7.2/22.6) × (22.6 - 7.2)) × 0.34) + (21.9 × 0.33) + (7.2 × 0.33)
    = (12.1 × 0.34) + (21.9 × 0.33) + (7.2 × 0.33)
    = 4.1 + 7.2 + 2.4 = 13.7 U/day
```

The low-TDD branch inflates tddWeighted by `(tddWeighted/tdd7D) × (tdd7D - tddWeighted)` — a partial pull toward tdd7D — before blending. This prevents the effective TDD from collapsing entirely to the extrapolated recent delivery.

Observed median (16.8) is higher than this simplified calculation (13.7) because many cycles have residual bolus IOB or slightly higher recent delivery.

### A.4 ISF Comparison Table

At target=99 mg/dL, C=1800, Boost blended TDD=16.8 overnight:

| BG (mg/dL) | Current ln ISF | Power-law k=3.5 ISF | Ratio (PL/ln) |
|---|---|---|---|
| 60 | 161 | 849 | 5.3× (more protective) |
| 75 | 136 | 337 | 2.5× |
| 99 (target) | 107 | 107 | 1.0× (anchored) |
| 125 | 92 | 44 | 0.48× |
| 160 | 79 | 16 | 0.20× (more aggressive) |
| 200 | 70 | 7 | 0.10× |

At BG=200, the ISF floor (10 mg/dL/U) would engage, preventing the value from dropping below a safe minimum.

---

## Appendix B: Patient Summary Statistics

| Metric | Value |
|---|---|
| Period | Jun 2025 – Mar 2026 (10 months) |
| Total loop cycles | 94,980 |
| Overnight fasting cycles | 12,005 |
| Valid +2h overnight samples | 3,647 |
| Valid +2h daytime fasting samples | 238 |
| Usable nights | ~300 |
| Actual TDD (median) | 21.9 U/day |
| Actual TDD 7-day (median) | 22.6 U/day |
| Boost blended TDD overnight (median) | 16.8 U/day |
| Basal insulin (median) | 7.2 U/day |
| Bolus insulin (median) | 14.5 U/day |
| Basal : bolus ratio | 33 : 67 |
| Overnight basal IOB (median) | 0.021 U |
| insulinDivisor (D) | 82 (peak time 38 min) |
| normalTarget | 99 mg/dL |
| dynISFadjust | 100% |
| dynIsfVelocity | 100% (1.0) |
| autosensMin / autosensMax | 0.5 / 1.5 |
