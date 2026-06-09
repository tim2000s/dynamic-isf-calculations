# Glucose-dependent insulin sensitivity in automated insulin delivery: resistance is real, but clearance cancels it

**Tim Street, with Claude (Anthropic)** · 2026-06-09
*Data: ~62,000 overnight and ~64,000 daytime fasting correction windows from 88–119 people using open-source automated insulin delivery (oref0 and Trio)*

---

## Abstract

Dynamic insulin sensitivity factor (ISF) algorithms lower the correction factor as glucose rises, on
the premise — well supported by molecular physiology — that hyperglycaemia induces insulin
resistance. Population data from thousands of patients (the Diabeloop glucose-related ISF model) and a
single-patient power-law analysis both recover a steep glucose dependence (ISF ∝ (target/BG)^k,
k≈3.5). Yet on outcome data, dynamic equations rarely beat a well-set static ISF. We resolve this
paradox. Using a *same-window* method — computing each window's effective ISF from the realised four-hour drop
and scoring every candidate against it (equivalently, rescaling the loop's own insulin-on-board
prediction), which removes the between-person confound exactly — we test static, the v1/v2 dynamic
equations, and the Diabeloop/power-law shapes on identical fasting windows. The dynamic
shapes degrade prediction; a per-user *level* fit captures essentially all recoverable individual
signal while the best shared glucose steepness is k=0. Decomposing the realised glucose drop into its
insulin and non-insulin components — estimating insulin-independent clearance directly from windows
where the loop expects no insulin action — shows that **high-glucose insulin resistance is real**
(insulin-only effective ISF falls with glucose, k≈0.44) **but is offset by insulin-independent
clearance (renal excretion and mass action), which rises with glucose, leaving the *net* realised ISF
flat.** The loop predicts and doses against the net, so a flat net effective ISF is correct; the
power law is right about insulin physiology but wrong as a predictor because it omits the offsetting
clearance and double-counts. The cancellation holds overnight and by day. We further show that
observational closed-loop data only *weakly identifies* the causal ISF — a dose-response analysis is
broken by reactive-dosing endogeneity. We recommend a per-user-adapted *net* effective ISF
(K/√TDD prior, online scale adaptation) with a near-target easing clamp, and no glucose-dependent
correction term.

---

## 1. Introduction

The insulin sensitivity factor (ISF) — expected glucose fall per unit of insulin — is the core
parameter of every correction-dosing decision in automated insulin delivery (AID). Static profiles
use one value; dynamic ISF (dynISF) algorithms vary it with current glucose and total daily dose
(TDD). The clinical rationale for the glucose dependence is strong: chronic hyperglycaemia drives
insulin resistance through well-characterised mechanisms — hexokinase-linked glycolytic overload
[1], disruption of the IRS-1/PI3K/Akt/GLUT4 cascade [2], oxidative and endoplasmic-reticulum stress
and mTOR/S6K1 feedback [3,4], ectopic lipid deposition [5,6], and glucotoxic β-cell decline [7,9].
The Diabeloop DBLG1 system, using a refactored oref controller, published a population glucose-ISF
curve derived from thousands of patients; an independent single-patient power-law analysis recovered
a similar steep dependence, ISF = (C/TDD)·(target/BG)^k with k≈3.5, improving two-hour prediction over
the inherited logarithmic scaler.

Against this, outcome evidence is stubbornly unfavourable to dynamic equations: across open-source AID
cohorts, a well-tuned static ISF predicts realised glucose drops about as well as the loop itself and
better than the v1 (ISF∝1/TDD) and v2 (∝1/TDD²) dynamic forms, whose TDD exponents are far steeper
than the empirically observed ≈−0.5. The paradox — strong physiological and population support for a
glucose-dependent ISF, weak outcome support for deploying one — motivates this work.

We show the paradox dissolves once the realised glucose drop is decomposed into its insulin and
non-insulin components.

## 2. Data and methods

**Cohort.** Per-tick closed-loop logs from open-source AID users were loaded into a TimescaleDB
instance: `oref_v5` (Trio) and `oref_v7` (oref0); a third table (AAPS-classic) was excluded because
its IOB/ISF accounting did not reconcile. Analyses use overnight (23:00–02:00) and daytime (09:00–
16:00) fasting windows, carbohydrate-screened by rejecting any window containing a glucose rise above
2 mg/dL per 5-min step, over a four-hour horizon. Units were cleaned (mmol-scale ISF ×18.018).

**Same-window head-to-head.** The loop's IOB-based glucose prediction is *linear in ISF*: predicted
drop = ISF × an activity integral that is independent of ISF. We therefore take the loop's own
prediction on each window (made with the ISF it ran) and rescale it to any candidate ISF, comparing
each to the realised end glucose **on the identical window**. This removes the between-person confound
exactly and works for every user. For each window we record glucose, TDD, IOB, hour, the loop's
predicted and realised drops, delivered correction insulin (super-micro-bolus units), the entry
trajectory (a 30-minute backward slope), and the realised ISF — the ISF that would have made the
prediction exact.

This is mathematically identical to computing each window's **effective ISF directly from the
realised drop** and comparing it to the candidate forms — the method a reader might expect. The
effective ISF is `effective_ISF = realised_drop / activity_integral`, where `activity_integral =
predicted_drop / sug_isf` is the insulin action the loop's model attributes to the window. The
candidate-ISF error then reduces exactly to
`err(candidate) = activity_integral × (candidate_ISF − effective_ISF)` (an identity we verified to
machine precision): each candidate is scored by its distance from the effective ISF, weighted by the
insulin that was actually acting. Two points follow. First, the comparison is **not circular** —
`sug_isf` cancels, so the effective ISF is independent of which ISF the loop ran. Second, the error
(mg/dL) framing is the *correctly precision-weighted* form of "effective ISF versus candidate": the
raw effective ISF is a noise-amplifying ratio that diverges when little insulin was acting, and
weighting by the activity integral down-weights exactly those uninformative windows. The one
assumption both framings share is the loop's insulin-action model, which supplies the activity
integral; §3.4 shows that model over-states action, so every ISF estimate here is conditional on it,
and the only model-free alternative (attributing the drop to delivered insulin, §5) is defeated by
reactive-dosing endogeneity.

**Confound controls.** The realised drop is not pure insulin action; it also reflects insulin-
independent clearance (renal excretion above the ~180 mg/dL threshold, and glucose effectiveness /
mass action), endogenous glucose production, and counterregulation near target. We control these
through: per-user centring (between-person resistance), entry-trajectory stratification (mean
reversion), explicit correction-magnitude terms (the loop over-trusts large corrections), and — the
key step — a data-derived estimate of insulin-independent clearance (§3.5).

**Statistics.** Out-of-sample evaluation uses grouped (leave-one-user-out) cross-validation with
users weighted equally; per-user effects are pooled by DerSimonian–Laird random-effects meta-analysis;
gradient-boosted models use SHAP attribution. All follow-on analyses run single-process on the cached
window dataset or one user at a time from the database.

## 3. Results

### 3.1 Static beats the dynamic equations

On identical windows, median absolute prediction error (mg/dL) was: loop 18.6, static profile ISF
20.3, v1 24.6, v2 49.7. A well-set static ISF essentially ties the loop and beats both dynamic forms;
v2 — the steepest TDD dependence — is far the worst. The empirical ISF–TDD relationship scales as
≈TDD^−0.5, consistent with a √TDD law and inconsistent with v1's −1 or v2's −2.

Equivalently, the window-level effective ISF (computed from the realised drop; §2) diverges from
every dynamic form as glucose rises (Table 1; each candidate normalised to the user's own profile to
remove the between-person level). The effective (net) ISF is flat-to-rising with glucose, whereas v1,
v2, and the Diabeloop quartic all fall steeply — i.e. the dynamic forms predict a declining ISF that
the net realised data does not show.

**Table 1. Window-level effective ISF versus candidate forms, by glucose (median, ÷ each user's
profile ISF).**

| glucose (mg/dL) | effective ISF | v1 (dynISF) | v2 (dynISF) | Diabeloop quartic |
|---|---|---|---|---|
| 100–120 | 0.71 | 0.93 | 2.89 | 0.87 |
| 120–145 | 0.82 | 0.87 | 2.35 | 0.68 |
| 145–175 | 0.91 | 0.72 | 1.57 | 0.53 |
| 175–205 | 0.89 | 0.60 | 1.08 | 0.42 |
| 205–260 | 0.94 | 0.50 | 0.78 | 0.34 |

(The effective ISF here is the *net*; §3.5 shows its insulin-only component does fall with glucose,
as the dynamic forms expect — but is offset by clearance.)

### 3.2 Every observational glucose slope is a different confound

The realised-ISF-versus-glucose relationship changed sign depending on the conditioning choice. A
minimum-predicted-drop filter flips the near-target sign (it selects rare high-IOB-at-low-BG windows
where counterregulation blunts the fall). Correction magnitude — the loop's predicted drop — is
correlated with glucose (r≈0.59) and, in a jointly-controlled gradient-boosted model, dominates
attribution (mean |SHAP| 2.12 versus 0.81 for glucose); a glucose-only model performs *worse* than
predicting each user's mean. Entering trajectory (mean reversion) inflates the high-glucose ratio.
Each confound manufactures a different apparent slope.

### 3.3 The dynamic shapes do not transfer; the level does

Applying the Diabeloop quartic and power-law shapes to the cohort (anchored to each user's profile
ISF at 100 mg/dL, curvature only), scored out-of-user: every shape degraded prediction, the steeper
the worse (power-law k=3.5 reached MAE 102 at 205–260 mg/dL versus 38 for static). A glucose-blind
correction of the loop's magnitude bias beat all of them. A best-fit individualised model
(actual_drop ≈ a_u + s_u·[pred_drop·(100/BG)^k], per-user intercept and scale, shared steepness k,
with k=0 a special case) put the entire recoverable individual signal in the per-user **level**
(static→17.0 MAE, a 7 mg/dL gain) and set the optimal shared glucose steepness to **k=0**.

### 3.4 The magnitude bias is affine, platform-invariant, and not an ISF term

The loop systematically over-predicts large drops: actual_drop ≈ 0.51·predicted_drop + 21 mg/dL,
affine with no consistent saturation, identical across platforms (Trio↔oref0 cross-applied MAE within
0.6). IOB adds nothing over predicted-drop magnitude. This is a constant insulin-action over-scale
plus a baseline drift — properties of the insulin-action model and basal, not the ISF formula.

### 3.5 Decomposition: resistance is real, but clearance offsets it

To separate insulin from non-insulin disposal we estimated the insulin-independent flux *from the
data*: in windows where the loop expects essentially no insulin action (|cgm − predicted settled
glucose| < 5 mg/dL), the observed four-hour change is the non-insulin flux at that glucose. It rises
monotonically and accelerates above ~180 mg/dL — the renal-threshold signature — reaching ~81 mg/dL
at 205–260 (Table 1).

**Table 2. Realised-ISF ratio (realised ÷ profile; <1 = insulin did less than expected), overnight.**

| glucose (mg/dL) | non-insulin flux (mg/dL/4h) | raw (net) ratio | clearance-corrected (insulin-only) ratio |
|---|---|---|---|
| 100–120 | 10 | 0.59 | 0.30 |
| 120–145 | 23 | 0.84 | 0.31 |
| 145–175 | 44 | 0.86 | 0.21 |
| 175–205 | 67 | 0.87 | 0.17 |
| 205–260 | 81 | 0.89 | 0.25 |

The **raw (net) ratio is flat** with glucose (k=−0.1). But once the data-derived clearance is removed,
the **insulin-only ratio falls with glucose (k≈0.44)**: per unit of insulin, glucose drops less at
high BG. **Insulin resistance is real.** The clearance required to reproduce the power law's
high-glucose aggression (~51–75 mg/dL) is precisely what the data shows (67–81), so the power law is
physiologically achievable — for the *insulin* component.

![Data-derived insulin-independent flux (left) rises with glucose and accelerates past ~180 mg/dL; the raw net realised-ISF ratio is flat while the clearance-corrected insulin-only ratio falls with glucose (right).](charts/inv008/fig_clearance_corrected_isf.png)

### 3.6 The cancellation holds by day, and per-user adaptation adds nothing

In daytime fasting windows the net realised ISF was flat and identical to overnight (k=−0.09 versus
−0.10), with a comparable clearance curve; the insulin-only daytime resistance looked steeper
(k≈2.8, noisy, consistent with the power law's reported daytime steepening) but was again offset in
the net. A nested cross-validation of per-user glucose steepness found that learning an individual k
beat a flat model for only 24% of users out-of-sample (below the 50% chance line), with zero median
gain: the apparent minority preferring a glucose curve was selection overfit.

![Net realised ISF by glucose, day versus overnight (left): both flat and near-identical. The non-insulin clearance flux (right) behaves the same in both regimes.](charts/inv008/fig_daytime_clearance.png)

## 4. The reconciliation

The central result is a cancellation:

> At high glucose, **insulin resistance is real** (insulin does less per unit), **but insulin-
> independent clearance — renal excretion and mass action — rises with glucose and roughly cancels it
> in the net glucose drop the loop acts on.** The *net* effective ISF is therefore flat.

This makes every prior result consistent. The Diabeloop population model and the power law measured
the *insulin* component and correctly found resistance. The best-fit individualised model found flat
k=0 best for *prediction*, because the loop predicts and doses against the *net*. The transplant
experiment found the power law *degrades* prediction, because it encodes the resistance but omits the
offsetting clearance — it double-counts, predicting too small a drop at high glucose and over-dosing.
The debate is not static-versus-dynamic: the physiology is resistant, the net is flat, a flat net
effective ISF predicts the net, and the power law is right about insulin and wrong as a predictor.

A corollary is clinically important: the cancellation is a population *average* that depends on renal
handling. Where clearance exceeds resistance — for example SGLT2-inhibitor users with enhanced
glucosuria — the net effective ISF at high glucose is even higher (insulin need lower); where
clearance is impaired — chronic kidney disease, dehydration — resistance dominates and the net becomes
genuinely glucose-dependent. A fixed glucose curve cannot track this; per-user online adaptation of
the net level can.

## 5. Identification limits

Observational closed-loop data only weakly identifies the causal ISF. A direct dose-response test —
regressing the realised drop on delivered correction insulin — returned physically impossible negative
sensitivities, because the controller's dose is *reactive*: it delivers more correction precisely when
glucose is not responding (in the 145–175 band, low-correction windows dropped 58 mg/dL, high-
correction only 38). The dose is endogenous, so a dose→outcome relationship cannot identify ISF. This
caution extends to any method that derives ISF from a closed loop's own doses, including a naive
realised-ISF = drop/dose ratio, which reactive dosing would bias downward at high glucose and so
could manufacture a power-law-shaped curve. We therefore relied on the loop's *counterfactual
prediction* (rescaled) rather than its doses, and on a data-derived clearance estimate; even so, the
clearance and resistance terms are partially entangled, so the resistance magnitude is an upper bound
rather than a point estimate. A definitive causal answer would require a fixed-dose natural experiment
or prospective shadow testing.

## 6. Recommended algorithm

The deployable ISF for an oref-family controller is:

> **ISF = s_u · (K/√TDD) · (100/BG)^k  + near-target easing clamp**, where K/√TDD is the population
> cold-start level at 100 mg/dL; **s_u is a per-user effective-sensitivity scale (with baseline)
> adapted online from the user's own outcomes** — the dominant, ~7 mg/dL lever, and not predictable
> cold-start from TDD or profile; **k is scheduled, ≈0.75 at cold-start fading to 0** as s_u is learned
> (a Bayesian shrinkage from population prior to individual); and the near-target clamp is a one-sided
> hypo-safety guardrail (raise ISF approaching target), motivated by counterregulation and risk
> asymmetry, not fitted to the realised drop.

In steady state this is an individually-calibrated, glucose-flat net ISF with a safety floor. No
glucose-dependent correction term is included: resistance is offset by clearance in the net, so a
resistance curve would double-count and over-dose hyperglycaemia. The residual error the analysis did
identify (the affine magnitude over-scale) belongs in the insulin-action model and basal, not the ISF.

## 7. Limitations

The analysis is restricted to fasting windows; active postprandial dynamics are excluded by the
carbohydrate screen and are a meal-bolus problem rather than a correction-ISF one. The clearance
estimate over-subtracts to an unknown degree (the corrected insulin-only ratio is low across the
range), so the resistance magnitude is bounded above. High-glucose low-insulin windows are sparse.
The cohort is open-source AID users, predominantly 2016–2023. Causal identification is weak for the
reasons in §5; the recommendation should be confirmed prospectively in shadow mode before deployment.

## 8. Conclusion

Glucose-dependent insulin resistance is real and physiologically well-founded, and it is visible in
closed-loop data once insulin-independent clearance is removed. But clearance rises with glucose and
cancels resistance in the net glucose response that an AID controller predicts and doses against. A
static — better, a per-user-adapted — net effective ISF is therefore not a failure to model the
physiology; it is the correct net of two opposing, glucose-dependent processes. The actionable lever
is online individualisation of the net level, not a glucose curve.

---

## References

1. Rabbani N, Thornalley P (2024). Hexokinase-linked glycolytic overload and unscheduled glycolysis in hyperglycemia-induced pathogenesis of insulin resistance, beta-cell glucotoxicity, and diabetic vascular complications. *Front Endocrinol.* doi:10.3389/fendo.2023.1268308
2. Khalid M, Alkaabi J, Khan M, et al. (2021). Insulin Signal Transduction Perturbations in Insulin Resistance. *Int J Mol Sci.* doi:10.3390/ijms22168590
3. Zhao X, An X, Yang C, et al. (2023). The crucial role and mechanism of insulin resistance in metabolic disease. *Front Endocrinol.* doi:10.3389/fendo.2023.1149239
4. Simon-Szabó L, Lizák B, Sturm G, et al. (2024). Molecular Aspects in the Development of Type 2 Diabetes and Possible Preventive and Complementary Therapies. *Int J Mol Sci.* doi:10.3390/ijms25169113
5. Galicia-Garcia U, Benito-Vicente A, Jebari S, et al. (2020). Pathophysiology of Type 2 Diabetes Mellitus. *Int J Mol Sci.* doi:10.3390/ijms21176275
6. Wang Y (2025). Triglycerides, Glucose Metabolism, and Type 2 Diabetes. *Int J Mol Sci.* doi:10.3390/ijms26209910
7. Beaupere C, Liboz A, Fève B, et al. (2021). Molecular Mechanisms of Glucocorticoid-Induced Insulin Resistance. *Int J Mol Sci.* doi:10.3390/ijms22020623
8. Allocca S, Monda A, Messina A, et al. (2025). Endocrine and Metabolic Mechanisms Linking Obesity to Type 2 Diabetes: Implications for Targeted Therapy. *Healthcare.* doi:10.3390/healthcare13121437
9. Młynarska E, Czarnik W, Dzięca N, et al. (2025). Type 2 Diabetes Mellitus: New Pathogenetic Mechanisms, Treatment and the Most Important Complications. *Int J Mol Sci.* doi:10.3390/ijms26031094
10. Yang B, Sherman A (2025). Crafting Mathematical Models for Type 2 Diabetes Progression: Leveraging Longitudinal Data. *bioRxiv.* doi:10.1101/2025.02.24.639807

## Reproducibility

All analysis code is open: `github.com/tim2000s/dynamic-isf-calculations` (package `inv008/`, run
each stage with `python -m inv008.<name>`). Key stages: `head_to_head` (same-window dataset),
`bridge_diabeloop` (shape transfer), `magnitude_bias`, `gradient_isf_fit` (individualised best fit),
`clearance_corrected_isf` (decomposition), `daytime_clearance`, `adaptive_k_nestedcv`,
`dose_response_db` (the inconclusive dose-response). The full working record, including counterfactuals
and the confound table, is in `docs/OREF-INV-008-Glucose-ISF-Investigation.md`.
