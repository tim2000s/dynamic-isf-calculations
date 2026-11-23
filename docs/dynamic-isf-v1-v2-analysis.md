# Dynamic ISF — analysis of the v1 and v2 equations

**2026-06-07** (same-window outcome test §3.6 added 2026-06-09) · Tim Street / Claude · Data: 171 people using open-source AID systems

> Audit update, 16 September 2026: the V2 replay used the intended no-`+1` equation with
> divisor/floor 75/76. The delta-IOB estimator is not independent physiological ground truth,
> so the claimed physiological TDD exponent and square-root recommendation are withdrawn. The
> equation replay and within-model prediction comparison remain valid. See
> `DYNAMIC-ISF-AUDIT-2026-09.md`.
>
> Correction, 23 November 2026: the figures below were rerun after the daily-dose readers were
> corrected and a unit conversion in the selected-window analysis was fixed (REVIEW.md item 20).

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
everywhere, giving a weaker correction than v1 on 98% of readings (a median of 3.3×, rising to
roughly 60× below 80 mg/dL and easing to about 1.6× above 200). That low-glucose behaviour is
sensible hypo protection.

The TDD exponent is the problem within the tested proxy. The delta-IOB sensitivity proxy scales
as roughly TDD^−0.4 to TDD^−0.7, depending on the fit; v2 assumes TDD^−2. As a between-person
predictor of the delta-IOB proxy, v2 is the worst of every form we tested (median error around
207 mg/dL per unit). Its overall level sits well above the proxy. This comparison does not establish which
TDD exponent describes physiological insulin sensitivity.

A separate within-loop prediction test also places v2 last. The loop's own glucose prediction was
rescaled to each candidate ISF on 62,824 selected overnight windows. A person's tuned static ISF
predicted the realised drop nearly as well as the loop itself. On the 57,151 windows where both
dynamic equations could be scored, v1's median absolute error was 2.2 mg/dL above static ISF's
(95% interval −0.2 to 4.9), so v1 was not shown to be less accurate than static ISF, while v2's
was 39.0 mg/dL above it (23.7 to 59.4). Because this reuses the loop's prediction model and selects windows using
future glucose, it is predictive evidence rather than an independent physiological measurement.

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
coverage of 99.4%. The equation implementations carry 27 unit tests against hand-computed
fixtures, covering the v2 glucose floor, the glucose-dependent v1/v2 ratio, both branches of
the TDD blend, and the missing-data gates. The companion methodology paper has the full
pipeline.

For 116 people the repository also fits glucose change against delta IOB over fasting windows.
This provides an observational comparison target. It does not serve as physiological ground truth
because in-window delivery changes IOB and is selected by the closed-loop controller.

---

## 3. Results

### 3.1 The v2/v1 ratio depends on glucose

![v2 vs v1: ratio vs glucose (left) and per-user ratio vs TDD (right)](charts/inv008/fig_v1_v2.png)

Median ISF_v2/ISF_v1 by glucose band (170 users, 9.5M readings):

| glucose band | median ISF_v2 / ISF_v1 |
|---|---|
| 40–80 | 60× |
| 80–100 | 6.8× |
| 100–120 | 3.8× |
| 120–150 | 2.8× |
| 150–200 | 2.1× |
| 200–360 | 1.6× |

v2 gives a higher ISF, meaning a weaker correction, on 98% of readings, a median of 3.3×. The
very high ISF below about 100 mg/dL (close to no correction at all) comes from the
`ln(BG/divisor)` term and the glucose floor, and amounts to strong hypo protection.

### 3.2 Which TDD law does the delta-IOB proxy follow?

![Log-log delta-IOB proxy vs TDD: proxy points, v1 slope −1, v2 slope −2, fitted slope −0.72](charts/inv008/fig_tdd_loglog.png)

On a log-log plot v1 has slope −1 and v2 slope −2, while the delta-IOB outcome proxy follows a
fitted slope near −0.72. This is shallower than v1 and far from v2. It shows how the equations
compare with that proxy. It does not show which TDD law physiological insulin sensitivity follows,
because new delivery changes IOB during the interval and the controller selects delivery in
response to glucose.

### 3.3 Agreement with the delta-IOB proxy and tuned profiles

Scoring each form as a between-person ISF predictor under leave-one-user-out cross-validation,
against the delta-IOB proxy (n=116) and tuned-profile ISF (n=141):

| candidate | median \|err\| vs proxy | log-err | within ±30% |
|---|---|---|---|
| K/√TDD (v-next) | 6.0 | 0.30 | 41% |
| 1700-rule | 18.6 | 0.63 | 14% |
| v1 (TDD⁻¹) | 29.0 | 0.86 | 5% |
| v2 (TDD⁻²) | 207 | 2.43 | 0% |

v2 is the worst between-person predictor of the delta-IOB proxy by a wide margin. Its level sits
far above that proxy. The ranking holds against tuned profiles too (v2 median error around 159
against 12.6 for √TDD).

The next check concerns software implementation, not physiological validity.

### 3.4 Device-log implementation validation

The v1 implementation reproduces what devices actually computed. For the nine Trio users on the
logarithmic form, replayed v1 tracks the device's logged ISF with stable per-person offsets
(median per-reading log-correlation 0.60, ratio IQR 0.29), the offsets reflecting the
unmodelled adjustment factor, divisor and autosens. Four users had switched the units of their
logged ISF mid-history, and a per-reading correction resolves all four.

### 3.5 Per-person view

A page for each of the 170 people (the curves of ISF against glucose at their TDD, a two-week sample of
dynamic-ISF traces over real glucose, and the per-reading ratio) is drawn by
`inv008/stage3_plots.py` into the private data folder. Per-person results are not published.

### 3.6 Same-window outcome test: which ISF predicts the realised drop?

The results above compare the ISF each equation *would have computed*. A separate within-loop
prediction test asks which rescaled ISF best predicts the glucose drop that *actually happened*.
The loop's IOB-based glucose
prediction is linear in ISF — predicted drop = ISF × an activity integral that does not depend
on ISF — so on any window we can take the loop's own prediction (made with the ISF it ran) and
rescale it to any candidate ISF, then compare to the observed end glucose. Every ISF form is
tested on the *same* window, which removes the between-person confound that the cross-validation
in §3.3 cannot: each person's outcomes only ever score their own equation there, whereas here
all four forms are scored on one shared set of windows.

Over 62,824 overnight, future-rise-screened, four-hour windows from 89 people (v1's `+1` makes
its prediction well defined throughout; AAPS users excluded for the same TDD-reconciliation
reason as elsewhere), scoring error as observed end glucose minus predicted end glucose. V1 and V2
were scored on the 57,151 of these windows that had a usable blended TDD:

![Prediction error by ISF form (left) and bias vs glucose (right)](charts/inv008/fig_head_to_head.png)

| ISF form | median \|error\| (mg/dL) | bias |
|---|---|---|
| loop (what ran) | 18.6 | +7.4 |
| static (tuned profile ISF) | 20.5 | +8.1 |
| v1 (TDD⁻¹) | 23.3 | −1.3 |
| v2 (TDD⁻²) | 60.1 | +52.5 |

The person's tuned static ISF is within 1.9 mg/dL of the loop's own calibrated ISF. On the 57,151
common windows static ISF's median absolute error was 21.1 mg/dL. V1's was 2.2 mg/dL above it
(95% interval −0.2 to 4.9, resampling people), an interval that includes zero, so this test does
not show v1 to be less accurate than static ISF. V2's was 39.0 mg/dL above it (23.7 to 59.4), about
three times static ISF's error (`results/common_window.md`). Counting which form is the single best
predictor per person (80 people with enough windows), the loop wins for 35, static for 22, v1 for
19 and v2 for only 4.

The bias-by-glucose panel shows how the dynamic forms' errors depend on glucose within this model. V1's error is
strongly glucose-dependent. With error defined as observed end glucose minus predicted end
glucose, its positive error at low glucose means it predicts too large a drop, while its negative
error at high glucose means it predicts too small a drop. V2 has positive error almost
everywhere, so it predicts too large a drop. A
near-constant static ISF carries less glucose-linked error in these selected windows.

This is a predictive comparison inside the loop's own linear IOB model, not an independent
measurement of insulin action. Within that scope, v1 was not distinguishable from a well-tuned
static level and v2 was worse than both.

---

## 4. Reading the result

V2 provides strong low-glucose protection. More ISF produces less insulin near the configured
floor. The observational analysis does not establish that its glucose curve represents
physiological sensitivity.

The TDD exponent is steep relative to the fitted proxy, −2 against between about −0.4 and −0.7 depending on the fit, which makes v2
the worst-fitting between-person predictor of that proxy among the forms tested. Its overall
level is about 3× the proxy level, so
corrections come out weak across the board: not only for the lighter-dosing majority but at
high glucose too, where v2 is still around 1.6× gentler than v1 above 200 mg/dL and more
aggression is usually wanted.

V2 offers strong low-glucose protection, but its TDD scaling and level performed poorly against
the tested outcomes. The square-root form fitted the observational targets better. This does not
establish a physiological TDD law or support a dosing recommendation.

---

## 5. Caveats

1. Counterfactual replay: the main comparison (§3.1–3.5) is of the ISF each equation would have
   computed, not closed-loop outcomes. The same-window test in §3.6 partly addresses this by
   scoring each ISF against the realised glucose drop, but it still rescales the loop's own
   linear prediction rather than re-running the controller, so second-order effects (changed
   insulin delivery feeding back into later IOB and glucose) are not captured.
2. Basal approximation: until 23 November 2026 the reconstruction for the 39 AAPS users did not
   read the export's temporary-basal records and used the profile schedule for basal. The
   figures here were rerun with those records read.
3. The delta-IOB benchmark is a per-person outcome proxy affected by in-window delivery,
   unrecorded carbohydrate, endogenous-glucose effects and controller feedback. It does not
   identify a physiological TDD law.
4. Single cohort: open-source AID users, mostly 2016–2023, with n = 116/141 for the
   proxy/profile analyses.

---

## Reproducibility

- Implementation and tests: `inv008/dynisf.py`, `inv008/tests/`
- v1 vs v2 comparison: `inv008/compare_v1_v2.py` → `results/v1_v2_comparison.*`, `charts/inv008/fig_v1_v2.png`
- Same-window outcome test: `inv008/head_to_head.py` → `results/head_to_head.{json,md}`, `results/head_to_head_windows.parquet`, `charts/inv008/fig_head_to_head.png`
- Pipeline and candidate search: `inv008/`, `fit_best_isf.py`
- Repository: `github.com/tim2000s/dynamic-isf-calculations`
