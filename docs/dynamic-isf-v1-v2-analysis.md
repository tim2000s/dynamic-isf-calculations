# Dynamic ISF — analysis of the v1 and v2 equations

**2026-06-07** · Tim Street / Claude · Data: 171 people using open-source AID systems

---

## Summary

Dynamic ISF sets correction sensitivity from total daily dose (TDD) and current glucose.
**v1** (Chris Wilson) makes the sensitivity anchor inversely proportional to TDD; **v2**
makes it inversely proportional to TDD squared and uses a glucose term `ln(BG/divisor)`
with glucose floored at `divisor+1`. We generated each person's dynamic ISF under both
equations from their real glucose and insulin history — ~9 million glucose readings across
171 people — and compared them.

Three findings:

1. **v1 and v2 differ on *both* axes — TDD and glucose.** They use different glucose terms
   (v1 `ln(BG/divisor+1)`, v2 `ln(BG/divisor)`), so the v2/v1 ratio is **glucose-dependent**:
   v2 produces a far higher ISF at low glucose, tapering at high glucose.
2. **v2 is much gentler than v1, and strongly hypo-protective.** It gives a higher ISF
   (weaker corrections) than v1 on **92% of readings** — a median **3.0×** v1, rising to
   **~53× below 80 mg/dL** (near-zero correction when low) and falling to **~1.5× above 200**.
3. **The v2 TDD exponent is too steep.** Sensitivity calculated independently from each
   person's own data follows ISF ∝ TDD^−0.5 or shallower; v2 is TDD^−2. As a between-person
   ISF predictor, v2 is by far the worst of every form tested (median |error| ≈ 171 mg/dL/U
   vs measured sensitivity) — its level sits well above where observed sensitivity lies.

So v2's *glucose* behaviour is in a sensible (hypo-protective) direction, but its *TDD*
exponent is too steep and its overall level sits well above both v1 and observed
sensitivity.

---

## 1. The two equations

Both derive a blended TDD, a sensitivity anchor at normal target, then scale by glucose.

| | sensitivity anchor at normal target | implied TDD law |
|---|---|---|
| **v1** | `1800 / (TDD · ln(target/divisor + 1))` | ISF ∝ 1/TDD |
| **v2** | `2300 / (ln(target/divisor) · TDD² · 0.02)` | ISF ∝ 1/TDD² |

v2 uses `ln(BG/divisor)` (no `+1`) in both the anchor and the glucose scaler, with glucose
**floored at `divisor+1`** so `ln(BG/divisor)` stays positive. In full long form (target 99
mg/dL, divisor 75, high cap 210):

```
v1:   ISF(BG) = 1800    / ( TDD  · ln(BG_capped/75 + 1) )
v2:   ISF(BG) = 115 000 / ( TDD² · ln(BG_floored/75) )      BG_floored = max(BG, 76)
```

The two glucose terms differ (v1 `ln(BG/75+1)`, v2 `ln(BG/75)`), so the between-equation
ratio is glucose-dependent:

```
ISF_v2 / ISF_v1 = (63.9 / TDD) · ln(BG/75 + 1) / ln(BG_floored/75)
```

The bracket is large at low glucose (the v2 log → 0 as BG approaches the floor) and ~1.3 at
high glucose — so v2 is dramatically more protective when low and modestly gentler when high.

---

## 2. Method in brief

For each person, every glucose reading was passed through both equations using the TDD
their device would have computed: device-logged TDD for Trio; for AAPS and OpenAPS users,
TDD reconstructed from raw delivery records (boluses + temp-basal segments over the profile
basal, on a 5-minute grid) through the same five-window blend. Relative timestamps were
re-anchored to absolute time and validated against recorded hour-of-day (median join
coverage 99.4%). The equation implementations carry 20 unit tests against hand-computed
fixtures (including the v2 collapse, the glucose floor, and the BG-dependent ratio).
Full methodology in the companion methodology paper.

For 114 people we also have sensitivity calculated independently from their own data (a
regression of glucose change on insulin absorbed over fasting windows), used as the ground
truth for the TDD law.

---

## 3. Results

### 3.1 The v2/v1 ratio is now glucose-dependent

![v2 vs v1: ratio vs glucose (left) and per-user ratio vs TDD (right)](charts/inv008/fig_v2updated.png)

Median ISF_v2/ISF_v1 by glucose band (170 users, 9.5M readings):

| glucose band | median ISF_v2 / ISF_v1 |
|---|---|
| 40–80 | **53×** |
| 80–100 | 6.2× |
| 100–120 | 3.5× |
| 120–150 | 2.5× |
| 150–200 | 1.9× |
| 200–360 | 1.5× |

v2 is weaker than v1 (higher ISF) on **92%** of readings, median **3.0×**. The strong
low-glucose ISF (near-zero correction below ~100 mg/dL) comes from the `ln(BG/divisor)`
term and the glucose floor — strong hypo protection.

### 3.2 Which TDD law does observed sensitivity follow?

![Log-log ISF vs TDD: observed points, v1 slope −1, v2 slope −2, fitted slope −0.56](charts/inv008/fig_tdd_loglog.png)

On a log-log plot, v1 is slope −1 and v2 slope −2; the calculated sensitivities follow a
fitted slope of about **−0.56** — shallower than v1, far from v2. Both equations
over-steepen the TDD dependence; v2, with twice the log-space slope, does so twice as hard.
The glucose-term change does not affect this — the TDD exponent is still −2.

### 3.3 Agreement with calculated sensitivity (and tuned profiles)

Scoring each form as a between-person ISF predictor (leave-one-user-out), against measured
sensitivity (n=114) and tuned-profile ISF (n=138):

| candidate | median \|err\| vs measured | log-err | within ±30% |
|---|---|---|---|
| K/√TDD (v-next) | **6.2** | **0.30** | **45%** |
| 1700-rule | 16.2 | 0.61 | 15% |
| v1 (TDD⁻¹) | 26.0 | 0.81 | 7% |
| **v2 (TDD⁻²)** | **171** | 2.32 | 2% |

v2 is overwhelmingly the worst between-person predictor — on top of the too-steep TDD
exponent, its level sits far above observed sensitivity. (Against tuned profiles the
ranking is the same: v2 median |err| ≈ 124 vs √TDD 12.8.)

### 3.4 Implementation validation against device-calculated ISF

Independently, the v1 implementation reproduces what devices actually computed: for the nine
Trio users on the logarithmic form, replayed v1 tracks the device's logged ISF with stable
per-person offsets (median per-reading log-correlation 0.60; ratio IQR 0.29) — the offsets
being unmodelled adjustment factor / divisor / autosens. (Four users had switched the
*units* of logged ISF mid-history; a per-reading correction resolves all four.)

### 3.5 Per-person view

![Example per-person page](charts/inv008/users/U073.png)

A page like this exists for each of the 170 people: ISF–glucose curves at their TDD, a
two-week sample of dynamic-ISF traces over real glucose, and the per-reading ratio.

---

## 4. Reading the result

v2's glucose behaviour is in a defensible direction — much more ISF (less insulin) at low
glucose is exactly the hypo protection a correction curve should provide, and it echoes the
glucose-dependent ISF established elsewhere (Diabeloop / power-law work). But two problems
remain:

- **The TDD exponent is too steep** (−2 vs observed −0.5), and v2 is the worst-fitting
  between-person predictor of any form tested.
- **The overall level is ~3× too high**, so corrections are very weak across the board —
  not just for the lighter-dosing majority but at high glucose too (still ~1.5× gentler than
  v1 above 200), where more aggression is usually wanted.

So v2 provides good hypo-protection but its TDD scaling and level are off. The v-next
proposal keeps the good idea (strong, glucose-dependent low-BG protection) but via a √TDD
level and a validated power-law glucose curve, anchored per-patient.

---

## 5. Caveats

1. **Counterfactual replay** — we compare the ISF each equation *would have computed*, not
   closed-loop outcomes.
2. **Basal approximation** — AAPS exports lack temp-basal records, so basal TDD for those 39
   users uses the profile schedule.
3. **Calculated-sensitivity benchmark** is per-person and may be biased low by unrecorded
   carbohydrate / endogenous-glucose effects; it tests the between-person TDD law.
4. **Single cohort** — open-source AID users, mostly 2016–2023; n = 114/138 for the
   calculated-sensitivity analyses.

---

## Reproducibility

- Implementation + tests: `inv008/dynisf.py`, `inv008/tests/`
- v1-vs-v2 comparison: `inv008/compare_v1_v2updated.py` → `results/v1_v2updated_comparison.*`, `charts/inv008/fig_v2updated.png`
- Pipeline / candidate search: `inv008/`, `fit_best_isf.py`
- Repository: `github.com/tim2000s/dynamic-isf-calculations`
