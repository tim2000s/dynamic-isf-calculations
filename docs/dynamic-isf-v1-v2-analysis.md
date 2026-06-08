# Dynamic ISF — analysis of the v1 and v2 equations

**2026-06-07** · Tim Street / Claude · Data: 171 people using open-source AID systems

---

## Summary

Dynamic ISF sets correction sensitivity from total daily dose (TDD) and current glucose. v1
(Chris Wilson) makes the sensitivity anchor inversely proportional to TDD. v2 makes it
inversely proportional to TDD squared, and uses the glucose term `ln(BG/divisor)` with glucose
floored at `divisor+1`. We ran each person's real glucose and insulin history through both
equations, about 9 million glucose readings across 171 people, and compared the sensitivity
each one would have produced.

The two equations differ on both axes. Because v1 keeps a `+1` in its glucose log and v2 does
not, the ratio between them changes with glucose: v2 produces a far higher ISF at low glucose
and only a modest one when high. Across the cohort v2 is the gentler equation almost
everywhere, giving a weaker correction than v1 on 92% of readings (a median of 3.0×, rising to
roughly 53× below 80 mg/dL and easing to about 1.5× above 200). That low-glucose behaviour is
sensible hypo protection.

The TDD exponent is the problem. Sensitivity calculated independently from each person's own
data scales as roughly TDD^−0.5, or shallower; v2 assumes TDD^−2. As a between-person predictor
of sensitivity, v2 is the worst of every form we tested (median error around 171 mg/dL per unit
against measured sensitivity), and its overall level sits well above where observed sensitivity
lies. So v2's glucose behaviour points the right way while its TDD scaling and level are off.

---

## 1. The two equations

Both derive a blended TDD, then a sensitivity anchor at normal target, then scale that by
glucose.

| | sensitivity anchor at normal target | implied TDD law |
|---|---|---|
| v1 | `1800 / (TDD · ln(target/divisor + 1))` | ISF ∝ 1/TDD |
| v2 | `2300 / (ln(target/divisor) · TDD² · 0.02)` | ISF ∝ 1/TDD² |

v2 uses `ln(BG/divisor)` (no `+1`) in both the anchor and the glucose scaler, with glucose
floored at `divisor+1` so the log stays positive. In long form, at target 99 mg/dL, divisor 75,
high cap 210:

```
v1:   ISF(BG) = 1800   / ( TDD  · ln(BG_capped/75 + 1) )
v2:   ISF(BG) = 115000 / ( TDD² · ln(BG_floored/75) )       BG_floored = max(BG, 76)
```

The glucose terms are not the same (`ln(BG/75+1)` for v1, `ln(BG/75)` for v2), so the ratio
between the two equations depends on glucose. It is largest at low glucose, where the v2 log
falls toward zero as BG approaches its floor, and shrinks to a small margin when glucose is
high. The cohort numbers for this are in §3.1.

---

## 2. Method in brief

For each person, every glucose reading was passed through both equations using the TDD their
device would have computed: device-logged TDD for Trio, and for AAPS and OpenAPS users a TDD
reconstructed from raw delivery records (boluses plus temp-basal segments over the profile
basal, on a 5-minute grid) through the same five-window blend. Relative timestamps were
re-anchored to absolute time and checked against recorded hour-of-day, with median join
coverage of 99.4%. The equation implementations carry 25 unit tests against hand-computed
fixtures, covering the v2 glucose floor, the glucose-dependent v1/v2 ratio, both branches of
the TDD blend, and the missing-data gates. The companion methodology paper has the full
pipeline.

For 114 people we also have sensitivity calculated directly from their own data, by regressing
glucose change on insulin absorbed over fasting windows. That serves as the ground truth for
the TDD law.

---

## 3. Results

### 3.1 The v2/v1 ratio depends on glucose

![v2 vs v1: ratio vs glucose (left) and per-user ratio vs TDD (right)](charts/inv008/fig_v1_v2.png)

Median ISF_v2/ISF_v1 by glucose band (170 users, 9.5M readings):

| glucose band | median ISF_v2 / ISF_v1 |
|---|---|
| 40–80 | 53× |
| 80–100 | 6.2× |
| 100–120 | 3.5× |
| 120–150 | 2.5× |
| 150–200 | 1.9× |
| 200–360 | 1.5× |

v2 gives a higher ISF, meaning a weaker correction, on 92% of readings, a median of 3.0×. The
very high ISF below about 100 mg/dL (close to no correction at all) comes from the
`ln(BG/divisor)` term and the glucose floor, and amounts to strong hypo protection.

### 3.2 Which TDD law does observed sensitivity follow?

![Log-log ISF vs TDD: observed points, v1 slope −1, v2 slope −2, fitted slope −0.56](charts/inv008/fig_tdd_loglog.png)

On a log-log plot v1 has slope −1 and v2 slope −2, while the calculated sensitivities follow a
fitted slope near −0.56: shallower than v1 and far from v2. Both equations over-steepen the
TDD dependence, and v2, with twice the log-space slope, does so twice as hard. The choice of
glucose term has no bearing on this; the TDD exponent is still −2.

### 3.3 Agreement with calculated sensitivity (and tuned profiles)

Scoring each form as a between-person ISF predictor under leave-one-user-out cross-validation,
against measured sensitivity (n=114) and tuned-profile ISF (n=138):

| candidate | median \|err\| vs measured | log-err | within ±30% |
|---|---|---|---|
| K/√TDD (v-next) | 6.2 | 0.30 | 45% |
| 1700-rule | 16.2 | 0.61 | 15% |
| v1 (TDD⁻¹) | 26.0 | 0.81 | 7% |
| v2 (TDD⁻²) | 171 | 2.32 | 2% |

v2 is the worst between-person predictor by a wide margin. On top of the over-steep TDD
exponent, its level sits far above observed sensitivity. The ranking holds against tuned
profiles too (v2 median error around 124 against 12.8 for √TDD).

### 3.4 Implementation validation against device-calculated ISF

The v1 implementation reproduces what devices actually computed. For the nine Trio users on the
logarithmic form, replayed v1 tracks the device's logged ISF with stable per-person offsets
(median per-reading log-correlation 0.60, ratio IQR 0.29), the offsets reflecting the
unmodelled adjustment factor, divisor and autosens. Four users had switched the units of their
logged ISF mid-history, and a per-reading correction resolves all four.

### 3.5 Per-person view

![Example per-person page](charts/inv008/users/U073.png)

A page like this exists for each of the 170 people: the ISF–glucose curves at their TDD, a
two-week sample of dynamic-ISF traces over real glucose, and the per-reading ratio.

---

## 4. Reading the result

v2's glucose behaviour is defensible. More ISF, and so less insulin, at low glucose is the hypo
protection a correction curve should provide, and it lines up with the glucose-dependent ISF
seen in the Diabeloop and power-law work. Two problems remain.

The TDD exponent is too steep, −2 against an observed −0.5, which makes v2 the worst-fitting
between-person predictor of any form tested. And the overall level is about 3× too high, so
corrections come out weak across the board: not only for the lighter-dosing majority but at
high glucose too, where v2 is still around 1.5× gentler than v1 above 200 mg/dL and more
aggression is usually wanted.

So v2 offers good hypo protection but its TDD scaling and level are wrong. The v-next proposal
keeps the useful idea, strong glucose-dependent protection at low BG, but delivers it through a
√TDD level and the Diabeloop glucose curve, anchored to each person.

---

## 5. Caveats

1. Counterfactual replay: we compare the ISF each equation would have computed, not closed-loop
   outcomes.
2. Basal approximation: AAPS exports lack temp-basal records, so basal TDD for those 39 users
   uses the profile schedule.
3. The calculated-sensitivity benchmark is per-person and may be biased low by unrecorded
   carbohydrate or endogenous-glucose effects; it tests the between-person TDD law.
4. Single cohort: open-source AID users, mostly 2016–2023, with n = 114/138 for the
   calculated-sensitivity analyses.

---

## Reproducibility

- Implementation and tests: `inv008/dynisf.py`, `inv008/tests/`
- v1 vs v2 comparison: `inv008/compare_v1_v2.py` → `results/v1_v2_comparison.*`, `charts/inv008/fig_v1_v2.png`
- Pipeline and candidate search: `inv008/`, `fit_best_isf.py`
- Repository: `github.com/tim2000s/dynamic-isf-calculations`
