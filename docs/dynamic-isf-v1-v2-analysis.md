# Dynamic ISF — analysis of the v1 and v2 equations

**2026-06-07** · Tim Street / Claude · Data: 171 people using open-source AID systems

---

## Summary

Dynamic ISF sets correction sensitivity from total daily dose (TDD) and current glucose.
The original equation (**v1**, Chris Wilson) makes the sensitivity anchor inversely
proportional to TDD; a later revision of the maths (**v2**) makes it inversely
proportional to TDD squared. We generated each person's dynamic ISF under both equations
from their real glucose and insulin-delivery history — roughly 9 million glucose readings
across 171 people — and compared them.

Three findings:

1. **v1 and v2 differ only through their TDD terms.** Compared at the same glucose and
   TDD, every other term cancels in the ratio between them, leaving
   ISF(v2)/ISF(v1) = **63.9/TDD** — independent of glucose. Below ~64 U/day v2 yields a
   higher ISF (weaker corrections); above it, a lower ISF (stronger corrections).
2. **77% of this cohort sits below the crossover**, so v2 weakens correction dosing for
   most people — by roughly 3× for someone using 20 U/day — and strengthens it only for
   the heaviest insulin users.
3. **Neither TDD power law matches observation, and v2 is the worse fit.** Sensitivity
   calculated independently from each person's own data follows ISF ∝ TDD^−0.5 or
   shallower. v1 assumes TDD^−1; v2 assumes TDD^−2. Against calculated sensitivity, v2 is
   the worst-fitting equation tested — worse than v1 and worse than the historical
   1700-rule.

The revision from v1 to v2 increased the steepness of the TDD dependence. The data shows
the true dependence is in fact **shallower than v1**, not steeper — so v2 moved in the
wrong direction.

---

## 1. The two equations

Both equations derive a blended TDD, then a sensitivity anchor at normal target, then
scale it by current glucose. The blend and the glucose scaler are identical between them.
They differ only in the anchor:

| | sensitivity anchor at normal target | implied law |
|---|---|---|
| **v1** | `1800 / (TDD · ln(target/divisor + 1))` | ISF ∝ 1/TDD |
| **v2** | `2300 / (ln(target/divisor + 1) · TDD² · 0.02)` | ISF ∝ 1/TDD² |

In full long form (defaults: target 99 mg/dL, divisor 75, glucose capped at 210):

```
v1:   ISF(BG) = 1800    / ( TDD   · ln(bg_capped/75 + 1) )
v2:   ISF(BG) = 115 000 / ( TDD²  · ln(bg_capped/75 + 1) )
```

Because the glucose term and the divisor are identical in both, they cancel **in the ratio
between the two equations** (within each equation they apply fully and the ISF falls as
glucose rises):

```
ISF_v2 / ISF_v1 = 2300 / (0.02 · 1800 · TDD) = 63.9 / TDD
```

The crossover where the equations agree is TDD = 63.9 U/day:

- 20 U/day → v2 ISF 3.2× higher → corrections ~3× smaller
- 40 U/day → 1.6× higher → corrections ~40% smaller
- 64 U/day → identical
- 100 U/day → 0.64× → corrections ~1.6× larger
- 150 U/day → 0.43× → corrections ~2.3× larger

---

## 2. Method in brief

For each person, every glucose reading was passed through both equations using the TDD
their device would have computed at that moment: device-logged TDD for Trio; for AAPS and
OpenAPS users, TDD reconstructed from raw delivery records (boluses plus temp-basal
segments over the profile basal schedule, on a 5-minute grid) and run through the same
five-window blend the equations use. Relative timestamps were re-anchored to absolute time
and validated against the recorded hour-of-day (median join coverage 99.4%; 5 of 148
reconstructed users flagged uncertain). The equation implementations carry 18 unit tests
against hand-computed fixtures. Full methodology is in the companion methodology paper.

For 114 people we also have sensitivity calculated independently from their own data (a
regression of glucose change on insulin absorbed over fasting windows), used here as the
ground truth for which TDD law reality follows.

---

## 3. Results

### 3.1 The crossover, observed

![Observed per-person median ISF ratio vs TDD, with the theoretical 63.9/TDD curve](charts/inv008/fig_crossover.png)

Per-person median ratios lie on the theoretical curve across the full TDD range
(13.6–272.9 U/day, median 47.5). **131 of 170 people (77%) fall below the 64 U/day
crossover** — for them, v2 computes weaker corrections than v1; for the 23% above it,
stronger.

### 3.2 Which TDD law does observed sensitivity follow?

![Log-log ISF vs TDD: observed points, v1 slope −1, v2 slope −2, fitted slope −0.56](charts/inv008/fig_tdd_loglog.png)

On a log-log plot, v1 is a line of slope −1 and v2 a line of slope −2. The calculated
sensitivities follow a fitted slope of about **−0.56** — *shallower* than v1, and far from
v2. Both equations over-steepen the TDD dependence: they over-estimate ISF for people
using smaller daily doses and under-estimate it for people using larger ones. v2, with
twice the log-space slope, does so roughly twice as hard.

### 3.3 Agreement with calculated sensitivity

![Equation ISF vs calculated sensitivity, v1 and v2 panels](charts/inv008/fig_empirical.png)

| | v1 (TDD⁻¹) | v2 (TDD⁻²) |
|---|---|---|
| Median absolute error (mg/dL per U) | **21.8** | 36.7 |
| Closer of the two | **112/138 (81%)** | 26/138 (19%) |
| Error, TDD < 64 U/day (n=110) | **24.4** | 53.1 |
| Error, TDD ≥ 64 U/day (n=28) | 10.9 | **9.5** |

v2's marginal advantage in the high-TDD band is a curve-crossing artefact, not evidence of
a better model: a slope-−2 line through this data must cross the observed cloud somewhere,
and it happens to do so around 60–100 U/day.

### 3.4 Implementation validation against device-calculated ISF

Independently of the comparison, we checked that the v1 implementation reproduces what
devices actually computed. Trio logs its own per-cycle ISF, giving ground truth for its
dynamic-ISF users. We expect, per person, a positive log-log correlation (curve shape
tracks the device) and a tight, roughly constant multiplicative offset — the device
additionally applies an adjustment factor, an insulin divisor, and an autosensitivity
ratio, which the replay does not model — rather than a ratio of exactly one.

All nine users on the logarithmic dynamic-ISF form track their devices with stable offsets
(median per-reading log-correlation 0.60; median ratio 0.52–0.90; median ratio
interquartile range 0.29). This corroborates that the replayed ISF is the equation's true
output, to within unmodelled per-person settings. (A data-quality note: four users had
switched the *units* of their logged ISF mid-history; a per-reading correction resolves
all four and affects only this validation read.)

### 3.5 Per-person view

![Example per-person page: ISF–glucose curves, two-week time series, ratio distribution](charts/inv008/users/U073.png)

A page like this exists for each of the 170 people: the ISF–glucose curves under each
equation at their TDD, a two-week sample of both dynamic-ISF traces over their real
glucose, and their per-reading ratio distribution.

---

## 4. Reading the result

The revision from v1 to v2 was a change of one quantity: the TDD exponent, from −1 to −2 —
a steeper dependence of sensitivity on total daily dose. The data points the other way.
Sensitivity *calculated from people's own glucose and insulin records* increases with TDD
far more gently than even v1 assumes — a slope near −0.5. Against that benchmark:

- **v1 is too steep**, but only moderately, and it remains the closer of the two for 81%
  of people;
- **v2 is much too steep**, and is the worst-fitting equation of every option tested
  against calculated sensitivity (see the companion proposal document for the full
  candidate comparison);
- in practice, adopting v2 *weakens* corrections for the 77% of people below ~64 U/day,
  and — set against measured sensitivity — over-estimates their ISF more severely than v1
  already does, exactly where the cohort is most concentrated.

The natural conclusion is that the correct TDD dependence is shallower than v1, not
steeper than it. That is the basis for the v-next proposal.

---

## 5. Caveats

1. **Counterfactual replay.** These people ran their own AID configurations; we compare
   the ISF each equation *would have computed*, not closed-loop outcomes.
2. **Basal approximation for one platform.** AAPS exports carry no temp-basal records, so
   basal TDD for those 39 users uses the profile schedule.
3. **Calculated-sensitivity benchmark.** The regression estimate is per-person and may be
   biased low by unrecorded carbohydrate or endogenous-glucose effects; it tests the
   between-person TDD law, not within-person glucose scaling (which is identical across
   equations).
4. **Single cohort.** Open-source AID users, mostly 2016–2023; n = 114/138 for the
   calculated-sensitivity analyses.

---

## Reproducibility

- Implementation + tests: `inv008/dynisf.py`, `inv008/tests/`
- Pipeline: `inv008/` (TDD reconstruction → ISF replay → figures)
- Per-person tables and figures: `charts/inv008/`; summary `charts/inv008/cohort_summary.json`
- Device validation: `inv008/validate_device_isf.py`
- Repository: `github.com/tim2000s/dynamic-isf-calculations`
