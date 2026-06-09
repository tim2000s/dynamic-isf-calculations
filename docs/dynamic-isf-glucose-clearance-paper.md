# The net insulin sensitivity factor does not fall with glucose: a same-window analysis of dynamic ISF in open-source automated insulin delivery

**Tim Street, with Claude (Anthropic)** · 2026-06-09 (revised)
*Data: ~62,000 overnight and ~64,000 daytime fasting correction windows from 73–119 people using open-source automated insulin delivery (predominantly oref0; some Trio)*

---

## Abstract

Dynamic insulin sensitivity factor (ISF) algorithms lower the correction factor as glucose rises, on
the premise — well supported by molecular physiology and by a large clinical population model
(Diabeloop) — that hyperglycaemia induces insulin resistance. Yet on outcome data, dynamic equations
rarely beat a well-set static ISF. We examine this paradox in open-source automated insulin delivery
(AID) data using a *same-window* method: each window's **effective ISF** is computed from the realised
four-hour glucose drop and compared to candidate forms on the identical window, removing the
between-person confound exactly. We compute the effective ISF two independent ways — by rescaling the
loop's own insulin-on-board prediction, and, to remove any dependence on the loop's insulin-action
model, by conservation from the insulin that actually acted (ΔIOB plus delivered insulin). Both agree
on the central result: **the *net* effective ISF — the quantity an AID controller predicts and doses
against — does not fall with glucose.** It is suppressed near target (counterregulation) and flat to
mildly rising above it, whereas the v1, v2, and Diabeloop dynamic forms all fall steeply; transplanted
onto this cohort the dynamic shapes degrade prediction, the steeper the worse. The recoverable
individual signal lives in the per-user *level*, not a glucose curve: a best-fit model assigns a
+7 mg/dL gain to per-user scale and an optimal shared glucose steepness of zero. We considered whether
a hidden high-glucose *insulin* resistance, offset by glucose-rising insulin-independent clearance,
could reconcile the resistance premise with the flat net response; we show this decomposition is **not
testable** in fasting AID data — clean insulin-free high-glucose windows are essentially absent and the
result flips sign with the window selection — and we therefore present it only as a hypothesis. We
recommend a per-user-adapted *net* effective ISF (K/√TDD prior, online scale adaptation) with a
near-target easing clamp, and no glucose-dependent correction term.

---

## 1. Introduction

The insulin sensitivity factor (ISF) — expected glucose fall per unit of insulin — is the core
parameter of every correction-dosing decision in AID. Static profiles use one value; dynamic ISF
(dynISF) algorithms vary it with current glucose and total daily dose (TDD). The clinical rationale for
the glucose dependence is strong: chronic hyperglycaemia drives insulin resistance through
well-characterised mechanisms — hexokinase-linked glycolytic overload [1], disruption of the
IRS-1/PI3K/Akt/GLUT4 cascade [2], oxidative and endoplasmic-reticulum stress and mTOR/S6K1 feedback
[3,4], ectopic lipid deposition [5,6], and glucotoxic β-cell decline [7,9]. The Diabeloop DBLG1 system,
using a refactored oref controller, published a population glucose-ISF curve derived from thousands of
patients; an independent single-patient analysis recovered a similar steep dependence,
ISF = (C/TDD)·(target/BG)^k with k≈3.5.

Against this, outcome evidence is unfavourable to dynamic equations: across open-source AID cohorts a
well-tuned static ISF predicts realised glucose drops about as well as the loop itself and better than
the v1 (ISF∝1/TDD) and v2 (∝1/TDD²) forms, whose TDD exponents are far steeper than the empirically
observed ≈−0.5. We ask a precise version of the question: does the ISF a controller should use — the
factor relating its insulin to the *realised net glucose drop* — actually fall with glucose?

A note on scope and honesty. An earlier draft of this work claimed to *resolve* the paradox by showing
that high-glucose insulin resistance is real but cancelled by insulin-independent clearance. A
subsequent critical audit (§6) found that decomposition is not reliably estimable in this data. This
revision restricts its claims to what survives that audit, and treats the resistance question as open.

## 2. Data and methods

**Cohort.** Per-tick closed-loop logs were loaded into TimescaleDB: `oref_v5` (Trio) and `oref_v7`
(oref0); an AAPS-classic table was excluded for non-reconciling IOB/ISF accounting. Analyses use
overnight (23:00–02:00) and daytime (09:00–16:00) fasting windows, carbohydrate-screened by rejecting
any window containing a glucose rise above 2 mg/dL per 5-min step, over a four-hour horizon. Units were
cleaned (mmol-scale ISF ×18.018; verified no uncleaned values remained).

**Effective ISF on the same window.** Each window's effective ISF is the ISF that explains its realised
drop: `effective_ISF = realised_drop / insulin_activity`. The candidate-ISF error reduces exactly to
`err(candidate) = insulin_activity × (candidate_ISF − effective_ISF)` (verified to 1.1×10⁻¹³), so
scoring a candidate by its error is identical to comparing it to the effective ISF, weighted by the
insulin acting — the correct precision weighting, since the raw effective ISF is a noise-amplifying
ratio. The comparison is not circular: the ISF the loop ran cancels.

We compute the insulin activity two ways:
- **Loop-prediction:** `activity = predicted_drop / sug_isf`, from the loop's IOB prediction
  (`reason_IOBpredBG`). Exact and convenient, but conditional on the loop's DIA/peak insulin-action
  model.
- **Conservation (model-independent):** `insulin_acted = (IOB_start − IOB_end) + SMBs delivered +
  ∫(temp_basal − profile_basal) dt`. No activity curve: ΔIOB is the observed insulin absorbed, plus the
  insulin delivered during the window. Because by four hours most of a fast-insulin dose has acted
  (≈85–93%), the residual IOB_end is small and the result is robust to the curve. This requires basal
  profiles, available predominantly for oref0 users.

**Confound controls.** The realised drop is not pure insulin action; it reflects insulin-independent
clearance (renal above ~180 mg/dL, and glucose effectiveness), endogenous glucose production, and
counterregulation near target. We control these through per-user centring, entry-trajectory
stratification, explicit correction-magnitude terms, and (where testable) a data-derived clearance
estimate. Out-of-sample evaluation uses grouped (leave-one-user-out) cross-validation; per-user effects
are pooled by DerSimonian–Laird meta-analysis. All follow-on analyses run single-process.

## 3. Results

### 3.1 Static beats the dynamic equations

On identical windows, median absolute prediction error (mg/dL): loop 18.6, static profile ISF 20.3, v1
24.6, v2 49.7. A well-set static ISF essentially ties the loop and beats both dynamic forms; v2 — the
steepest TDD dependence — is far the worst. The empirical ISF–TDD relationship scales as ≈TDD^−0.5.

### 3.2 The net effective ISF does not fall with glucose

Computed both ways, the net effective ISF rises gently with glucose, while every dynamic form falls
(Table 1; each ÷ the user's own profile ISF to strip the between-person level).

**Table 1. Net effective ISF versus candidate forms, by glucose (median, ÷ profile ISF).**

| glucose (mg/dL) | effective ISF (loop-pred) | effective ISF (conservation) | v1 | v2 | Diabeloop quartic |
|---|---|---|---|---|---|
| 100–120 | 0.71 | 0.20 | 0.93 | 2.89 | 0.89 |
| 120–145 | 0.82 | 0.37 | 0.87 | 2.35 | 0.68 |
| 145–175 | 0.91 | 0.47 | 0.72 | 1.57 | 0.52 |
| 175–205 | 0.89 | 0.53 | 0.60 | 1.08 | 0.41 |
| 205–260 | 0.94 | 0.45 | 0.50 | 0.78 | 0.33 |

The two effective-ISF estimates differ in *level* — the conservation estimate is biased low, because
the temp-basal-deviation term over-counts insulin where profile basal is mis-set — but agree in
*direction*: both rise (implied slope k≈−0.2 to −0.4), opposite to the dynamic forms (quartic k≈+1.3).
The rise is robust to entry trajectory: restricting to flat-entry windows gives [0.24, 0.37, 0.43,
0.56, 0.56], essentially unchanged. The shape is a near-target suppression (counterregulation blunting
the fall as glucose approaches normal) that recovers above target — not a decline. Transplanting the
dynamic shapes onto the cohort confirms the mismatch: anchored to each user's profile and scored
out-of-user, every shape degraded prediction, the steeper the worse (power-law k=3.5 reached MAE 102 at
205–260 mg/dL versus 38 for static), because each predicts a drop the net data does not deliver.

![Net effective ISF versus glucose, computed model-independently (conservation, ΔIOB-based) and via the loop's prediction — both rise — against the Diabeloop quartic, which falls.](charts/inv008/fig_effective_isf_independent.png)

### 3.3 The actionable signal is the per-user level, not a glucose curve

A best-fit individualised model — `actual_drop ≈ a_u + s_u·[pred_drop·(100/BG)^k]`, per-user intercept
and scale, shared steepness k, within-user cross-validated — put the recoverable individual signal in
the per-user **level** (static 24.2 → 17.0 MAE, a 7 mg/dL gain) and set the optimal shared glucose
steepness to **k = 0**. The per-user scale is not predictable from TDD (ρ=0.03) or profile ISF, so it
must be learned online. A nested cross-validation found that learning a per-user glucose steepness beat
a flat model for only 24% of users out-of-sample (below the 50% chance line): the apparent minority
preferring a glucose curve was selection overfit.

### 3.4 What inflates the apparent glucose dependence

The realised drop is dominated by **correction magnitude**, not glucose: `actual_drop ≈ 0.4·pred_drop +
25`, affine, platform-invariant, with the loop over-predicting large drops by ~2×. Glucose and
correction magnitude are correlated (r≈0.54), and in a jointly-controlled gradient-boosted model
magnitude dominates attribution while a glucose-only model performs *worse* than predicting each user's
mean. The over-prediction of large corrections is properly a property of the insulin-action model and
basal, not the ISF.

### 3.5 Daytime

In daytime fasting windows the net effective ISF was flat and indistinguishable from overnight
(k=−0.09 versus −0.10), so the conclusion is not specific to the overnight regime. Active postprandial
dynamics are excluded by the carbohydrate screen and are a meal-bolus problem, not a correction-ISF
one.

### 3.6 Two individual case studies (external Nightscout data)

To test the method on data entirely outside the cohort, we applied it to two individual open-source
AID users via their live Nightscout instances, computing the effective ISF both ways (loop-prediction
and model-independent conservation) over the available history. The two gave opposite results — and
the difference is informative.

**User A — AAPS DynamicISF, fully closed-loop without announced carbohydrate (UAM).** Twelve months,
409 fasting windows. The effective ISF **falls steeply with glucose**, and critically this is
confirmed by the model-independent conservation estimate (slope k≈2.0 both ways), so it is **not** the
insulin-action-model artifact. A glucose-dependent ISF predicts this person's realised drops far
better than a static profile (median error: Diabeloop quartic 25, the dynISF they run 28, static
profile 49 mg/dL). For this individual, dynamic ISF demonstrably earns its keep.

**User B — oref with autosens, carbohydrate-aware.** Five months available, 610 fasting windows. The
effective ISF is **flat to mildly rising** (loop k≈0, conservation k≈−0.8), matching the cohort. A
static or autosens-adjusted ISF and the gentle v1 form all predict well (median error ≈15–19 mg/dL),
while the steep Diabeloop quartic and v2 forms are worse (21 and 69). For this individual a glucose
curve does not help.

These do not contradict the cohort; they illustrate the heterogeneity it contains, and the split
tracks **carbohydrate-announcement behaviour rather than algorithm**. User B announces carbohydrate,
so its fasting windows are genuinely carb-free and the effective ISF is flat. User A is UAM, so
*uncovered carbohydrate load is present in every window*; uncovered carbs at high glucose depress the
drop-per-unit-insulin and read as a falling ISF. The most likely interpretation is therefore that User
A's apparent glucose dependence is substantially uncovered-carb dynamics that its dynISF usefully
*compensates* for — not pure glucotoxic resistance. Each case is N=1 (User A 12 months, User B 5
months) and the conservation *level* is biased low (basal under-counting), so we read only the
slopes, which are robust.

## 4. The resistance question is open, not resolved

Molecular physiology and the Diabeloop population model say hyperglycaemia causes insulin resistance,
yet the net effective ISF does not fall with glucose. One reconciliation is that high-glucose *insulin*
resistance is real but offset by glucose-rising insulin-independent clearance, leaving the net flat.
This is physiologically plausible — the loop predicts and doses against the net, so a flat net is what a
controller needs regardless — but **we cannot confirm it in this data.** Estimating the clearance
requires fasting windows with high glucose and essentially no insulin acting; these are only ~2% of
windows, are increasingly contaminated by rising-entry (uncovered-carbohydrate-suspect) windows as
glucose climbs (11%→44%), and the truly steady, insulin-free, carbohydrate-free high-glucose window is
essentially absent (n≈16 at 205–260 mg/dL). The implied resistance flips sign with the window selection
(corrected slope +0.85 with the contaminated estimate, −4.22 with the cleaner but unstable one). We
therefore present the resistance-offset-by-clearance account as a hypothesis only, and base no
recommendation on it. The actionable point is unaffected: the net does not fall, and that is what an AID
controller acts on.

## 5. Identification limits

Observational closed-loop data only weakly identifies the causal ISF. A direct dose-response test —
regressing the realised drop on delivered correction insulin — returned physically impossible negative
sensitivities, because the controller's dose is *reactive* (it delivers more correction precisely when
glucose is not responding; in the 145–175 band, low-correction windows dropped 58 mg/dL, high-
correction only 38). The dose is endogenous, so a dose→outcome relationship cannot identify ISF; this
caution extends to any naive realised-ISF = drop/dose ratio, which reactive dosing would bias downward
at high glucose. We therefore relied on the loop's counterfactual prediction (rescaled) and the
model-independent conservation estimate, neither of which regresses on dose. The insulin/non-insulin
decomposition (§4) is the part this data cannot support; a definitive causal answer would need a
fixed-dose natural experiment or prospective shadow testing.

## 6. Audit

The findings above were re-derived from the raw data and stress-tested. Sign conventions, the
effective-ISF identity (1.1×10⁻¹³), unit cleaning, the static-versus-dynamic ranking, the magnitude
bias, and the net-effective-ISF rise (robust across both computation methods and across entry
trajectory) all held. The one result that did not survive was the clearance-versus-resistance
decomposition (§4), whose high-glucose estimate proved selection-dependent and unstable; it has been
demoted to a hypothesis accordingly. No residual code defects were found beyond two already fixed (an
array-broadcast error and a coefficient-scaling error).

## 7. Recommended algorithm

> **ISF = s_u · (K/√TDD) · (100/BG)^{k}  + near-target easing clamp**, where K/√TDD is the population
> cold-start level at BG 100; **s_u is a per-user effective-sensitivity scale (with baseline), adapted
> online from the user's own outcomes** (the dominant ~7 mg/dL lever, not predictable cold-start);
> **k is scheduled, ≈0.75 at cold-start fading to 0** as s_u is learned; and the near-target clamp is a
> one-sided hypo-safety guardrail (raise ISF approaching target), motivated by counterregulation and
> risk asymmetry, not fitted to the realised drop.

In steady state this is an individually-calibrated, glucose-flat *net* ISF with a safety floor. No
glucose-dependent correction term is included: the net effective ISF does not fall with glucose, so a
declining ISF would predict drops the data does not deliver and over-dose hyperglycaemia. The residual
magnitude over-prediction belongs in the insulin-action model and basal.

## 8. Limitations

Fasting windows only; active postprandial dynamics excluded. The model-independent estimate is
predominantly oref0 (Trio basal profiles were not assembled) and its *level* is biased low by
basal-deviation over-counting — directions are sound, absolute levels are not. The resistance/clearance
decomposition is not testable here (§4). Causal identification is weak (§5). The recommendation should
be confirmed prospectively in shadow mode before deployment.

## 9. Conclusion

The ISF that an AID controller should use — the factor relating its insulin to the realised *net*
glucose drop — does not fall with glucose; it is suppressed near target by counterregulation and flat
to mildly rising above. Dynamic equations that lower ISF as glucose rises therefore predict drops the
net data does not deliver and degrade prediction. Whether a hidden insulin-only resistance exists
beneath an offsetting clearance is a reasonable hypothesis this data cannot settle. The actionable
conclusion is independent of that question: individualise the *net* ISF level online, keep a near-target
safety clamp, and do not add a glucose curve by default. Two external individual case studies (§3.6)
sharpen this: a carbohydrate-aware user matched the cohort (flat, model-confirmed), while a UAM user
showed a genuine falling effective ISF that a dynISF helps — most plausibly by compensating for routine
uncovered-carbohydrate load rather than for glucotoxic resistance. A glucose term thus earns its place
mainly where uncovered carbohydrate is routine; per-user adaptation should learn it for those
individuals while defaulting to flat for the carbohydrate-aware majority.

---

## References

1. Rabbani N, Thornalley P (2024). *Front Endocrinol.* doi:10.3389/fendo.2023.1268308
2. Khalid M, Alkaabi J, Khan M, et al. (2021). *Int J Mol Sci.* doi:10.3390/ijms22168590
3. Zhao X, An X, Yang C, et al. (2023). *Front Endocrinol.* doi:10.3389/fendo.2023.1149239
4. Simon-Szabó L, Lizák B, Sturm G, et al. (2024). *Int J Mol Sci.* doi:10.3390/ijms25169113
5. Galicia-Garcia U, Benito-Vicente A, Jebari S, et al. (2020). *Int J Mol Sci.* doi:10.3390/ijms21176275
6. Wang Y (2025). *Int J Mol Sci.* doi:10.3390/ijms26209910
7. Beaupere C, Liboz A, Fève B, et al. (2021). *Int J Mol Sci.* doi:10.3390/ijms22020623
8. Allocca S, Monda A, Messina A, et al. (2025). *Healthcare.* doi:10.3390/healthcare13121437
9. Młynarska E, Czarnik W, Dzięca N, et al. (2025). *Int J Mol Sci.* doi:10.3390/ijms26031094
10. Yang B, Sherman A (2025). *bioRxiv.* doi:10.1101/2025.02.24.639807

## Reproducibility

Code: `github.com/tim2000s/dynamic-isf-calculations`, package `inv008/` (`python -m inv008.<name>`).
Stages: `head_to_head`, `effective_isf_independent` (conservation), `bridge_diabeloop`, `magnitude_bias`,
`gradient_isf_fit`, `clearance_independent` / `clearance_corrected_isf` (the untestable decomposition,
retained for transparency), `daytime_clearance`, `adaptive_k_nestedcv`, `dose_response_db`.
