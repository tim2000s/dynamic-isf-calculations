# OREF-INV-008 — Does ISF depend on glucose? A same-window investigation

**Date:** 2026-06-09 · Tim Street / Claude · Working record of process, methods, findings, and counterfactuals

---

## 0. TL;DR

We set out to settle, on real open-source AID data, whether the insulin sensitivity factor (ISF)
should vary with glucose — and if so, how to build a deployable update for oref instantiations
(AAPS / Trio / oref0). The arc:

- A **same-window method** (rescaling the loop's own IOB prediction to any candidate ISF) lets us
  test every ISF form on the *identical* window, removing the between-user confound exactly.
- **Static profile ISF predicts the realised drop about as well as the loop itself, and far better
  than the dynamic v1/v2 equations** (which over-steepen the TDD dependence).
- The apparent glucose dependence of ISF is **dominated by confounds**: between-user level,
  a minimum-predicted-drop selection effect near target, **correction-magnitude** (`pred_drop`),
  and **mean reversion** at high glucose. Each confound manufactures a *different* glucose slope.
- A properly **individualised best-fit** model (shared glucose shape × per-user scale, k=0 a special
  case) puts essentially all the recoverable signal in the **per-user level/scale (+7 mg/dL)** and
  sets the average glucose steepness to **k≈0**. 36% of users individually prefer k≥1, but that is
  selection-inflated and not predictable from TDD/profile.
- The **Diabeloop n≈1000s power law** (real, large, within-user, *same algorithm and method*) is the
  serious counter-evidence. We have **conceded** several of our earlier objections to it. The live
  disagreement is now narrow and specific: at **high glucose**, our realised-ISF curve **rises**
  (insulin stays effective) where the power law needs it to **fall** (resistance).
- **The high-BG question is RESOLVED into a clearance-vs-resistance cancellation (§3.13).** Insulin
  resistance at high glucose is **real** (corrected insulin-only ISF falls, k≈0.44), but insulin-
  independent clearance (renal + mass-action) rises with glucose and **roughly cancels it in the net
  drop the loop acts on** — so the *net* effective ISF is flat. The power law is right about insulin
  physiology but wrong as a *predictor* (it omits the offsetting clearance and double-counts). Test #2
  (dose-response) was inconclusive — broken by reactive-dosing endogeneity. Net: a flat per-user net
  ISF is correct for the loop; do not add the insulin-only resistance curve.

**Current recommendation (provisional, overnight-fasting scope):**
> ISF = a **per-user-adapted effective-sensitivity scale** (`K/√TDD` as the cold-start prior →
> adapt scale + baseline online from the user's own outcomes) **+ a near-target easing safety
> clamp.** Glucose steepness `k=0` by default, with a gentle cold-start `k≈0.75` that fades as the
> per-user scale is learned. A genuine high-BG resistance term remains *possible* and is under test.

---

## 1. The question and why it is hard

ISF is the expected mg/dL fall in glucose per unit of insulin. Static profiles use one value; dynamic
ISF (dynISF, Diabeloop) vary it with glucose and total daily dose (TDD). The clinical claim is that
**at high glucose insulin is less effective (glucotoxic resistance) so ISF should be lower (dose
harder), and near/below target insulin should be eased off for hypo safety.**

The difficulty is that the realised glucose drop in a fasting window is **not pure insulin action**:

```
ΔBG (fall) = insulin disposal  +  insulin-independent clearance (renal >~180, glucose effectiveness)
             − endogenous glucose production (EGP; hepatic, dawn)  ± residual carb absorption
           all gated by counterregulation as BG approaches normal
```

ISF is asked to absorb the whole sum. So a naive "realised ISF vs glucose" curve conflates true
insulin sensitivity with these glucose-dependent non-insulin fluxes — and with how the loop and the
person co-vary glucose with insulin dose. Every observational metric is confounded; they are
confounded *differently*, which is why they disagree.

---

## 2. Method foundation — the same-window head-to-head

**Key property:** the loop's IOB-based BG prediction is *linear in ISF* — predicted drop =
ISF × (an activity integral that does not depend on ISF). So on any window we can take the loop's own
prediction (made with the ISF it actually ran, `sug_isf`) and **rescale it to any candidate ISF**,
then compare to the realised end glucose **on the identical window**. This removes the between-user
confound completely and works for every user.

- `inv008/head_to_head.py` — builds the per-window dataset over **overnight (23:00–02:00),
  carbohydrate-screened, 4-hour-horizon** windows. Per window it stores: `bg`, `tdd`, `iob`, `hour`,
  `sug_isf`, `profile_isf`, `isf_v1`, `isf_v2`, `realised_isf`, the prediction errors of each form,
  `start_slope` (15-min forward), `pre_slope` (30-min backward — momentum *entering* the window),
  and `bg_end`. Units cleaned (mmol-scale ISF ×18.018). v6 excluded (iob/isf accounting did not
  reconcile). Output: `results/head_to_head_windows.parquet` (**~62,300 windows, 80–89 users,
  oref_v5 Trio + oref_v7 oref0**).
- `realised_isf` = the ISF that would have made the prediction exact = `sug_isf · actual_drop /
  predicted_drop`. It is a noise-amplifying ratio (tails to ±thousands) — used only via robust
  medians or bounded transforms.

**The crash (process note).** The first two runs *hard-froze the Mac* (2026-06-08, two reboots, no
panic file = swap-death). Cause: `ai = pdl / sugisf` used the full ~615k-element array instead of
`sugisf[i]`, so `e_static/e_v1/e_v2` became arrays appended into every row → ~40 GB in one worker ×
12 workers → 64 GB exhausted. Fixed to `sugisf[i]` → 0.57 GB/worker. This produced a standing
**safe-run protocol**: measure one worker before going parallel; prefer the cached parquet (most
follow-ons are single-process, tens of MB) over re-running the DB pipeline.

---

## 3. What we looked at, in order

### 3.1 Static vs dynamic (v1/v2), same windows
`head_to_head.py` scored each ISF form by the realised-drop MAE on the same windows.

| form | MAE (mg/dL) | bias | per-user best |
|---|---|---|---|
| loop (what ran) | 18.6 | +7.4 | 35 |
| **static profile ISF** | **20.3** | +6.6 | 24 |
| v1 (∝1/TDD) | 24.6 | −10.4 | 17 |
| v2 (∝1/TDD²) | 49.7 | +25.2 | 4 |

**Finding:** a well-set static ISF essentially ties the loop and beats both dynamic equations; v2 is
far the worst; v1 is systematically BG-dependent (the steep 1/TDD coupling). This corroborates the
independent cohort finding that empirical ISF∼TDD scales like **TDD^−0.5 (√TDD)**, not v1's −1 or
v2's −2.

### 3.2 Characterising the static-ISF error vs glucose
`inv008/err_curve.py`, `inv008/err_consistency.py`. The error `err_static = actual_end −
predicted_end(static)` (>0 = static over-predicted the drop), in two views:

- **Absolute (mg/dL):** U-shaped — large near target, minimum mid-range, an up-turn at high BG.
- **Fractional (effective-sensitivity gap = 1 − realised/profile):** monotone — 83% at 80–100 →
  ~13% at 175–260. The two disagree by construction: corrections are bigger at high BG, so a flat
  proportional mismatch shows as a larger mg/dL error there.

Random-effects pooling (`err_consistency.py`): the **near-target effect is a strong cross-user law**
(fractional slope p≈0, 96% of users share the sign after confound control); the **high-BG up-turn is
inconsistent** (curvature 55–64% same-sign, n.s. after control). Per-user curves overlaid show the
spread.

### 3.3 The direction debate and the confound reconciliation
Two observational metrics disagreed on the glucose *direction*:

- **Calibration audit** (all windows, loop ISF, pooled): error rises from −3.6 near target to +6 at
  high BG ⇒ realised ISF *falls* with glucose ⇒ the **dynISF/power-law direction**.
- **Head-to-head** (predicted-drop ≥ 8 filter, per-user): over-prediction worst near target ⇒ the
  *opposite* gradient.

We traced the disagreement to a **minimum-predicted-drop selection effect**: near target the loop
usually predicts ~0 drop, so requiring `pred_drop ≥ 8` keeps only the rare high-IOB-at-low-BG windows
where counterregulation then blunts the fall (median predicted 16.8, actual 4.0). The calibration
audit includes the dropped windows. So the two are not measuring the same sub-population near target;
neither cleanly settles direction. **Lesson:** every univariate glucose slope here is manufactured by
a *different* conditioning choice.

### 3.4 Joint-control ML and the physiological reframe
`inv008/isf_pattern_ml.py` — a leave-one-user-out gradient-boosted model, target centred within user
(so only shared, actionable structure can be learned), glucose competing against momentum
(`pre_slope`, the IOB-orthogonalised endogenous momentum), `iob`, `pred_drop`, `hour`.

- `pred_drop` (correction magnitude) **dominates** (mean|SHAP| 2.12 vs bg 0.81).
- bg survives joint control but adds only ~0.2 mg/dL; **bg-only is worse than the per-user-mean
  baseline** (a glucose-alone ISF — structurally what dynISF is — underperforms doing nothing).
- bg's SHAP curve is monotone (+2.69 near target → −1.34 at high BG): net of confounds, static ISF
  over-predicts near target, under-predicts when high.

**Physiology:** the near-target over-prediction is **counterregulation** (the body defends against
hypo, blunting the fall) — a *safety brake, not a sensitivity*. The high-BG under-prediction is
**mass-action + renal clearance** outweighing glucotoxic resistance in the acute overnight window.
**Critical safety point:** the low *realised* ISF near target must NOT be fitted (it would dose more
into a defended low); it is a signal to respect, hence the easing clamp.

### 3.5 The Diabeloop / inv004 g(BG) work (n≈1000s)
The prior body of work (other Drive, `Dynamic ISF data/`) established glucose-dependent ISF
externally:

- **Diabeloop ADA poster:** ISF-vs-glucose from **thousands of patients**, piecewise polynomial
  (quartic >100: `272 − 3.121G + 0.01511G² − 3.305e-5G³ + 2.69e-8G⁴`), individualised by a single
  multiplicative per-patient scale. **Algorithm-fit** (loop-correction accuracy), control algorithm
  = **refactored oref**.
- **N=1 power law** (Tim's 3,647 overnight samples): `ISF = (C/TDD)(target/BG)^k`, **k≈3.5**, 12–18%
  better 2h-prediction than the log scaler. Found the *same* overnight over-prediction we did
  (**+11.2 mg/dL** vs our +6–7).
- Hybrid (no-TDD): Diabeloop quartic ≥105 + power-law tail `75.8·(105/BG)^3.5` below; conservative
  <90 by design (safety, not fit).

### 3.6 Does the Diabeloop / power-law SHAPE transfer to oref?
`inv008/bridge_diabeloop.py` — apply each shape anchored to each user's profile ISF at 100 mg/dL
(curvature only), score out-of-user, with a magnitude-only control.

| candidate | MAE | vs static | vs magnitude-only |
|---|---|---|---|
| static | 24.2 | — | — |
| power-law k=2 | 34.8 | −10.6 | −15.7 |
| power-law k=3.5 | 43.9 | −19.7 | −24.8 |
| diabeloop hybrid | 29.8 | −5.7 | −10.8 |
| **magnitude-only (glucose-blind)** | **19.1** | — | — |
| bg-only | 19.4 | | |
| bg + pred_drop | 15.7 | | |

The transplanted shapes **degrade** oref badly (steeper = worse; k=3.5 MAE 102 at 205–260 — it
predicts a tiny drop where realised drops are large). A glucose-blind magnitude correction beats them
all; bg-only ≈ magnitude-only; bg+pred_drop adds a modest interaction. *(Caveat raised later: anchor
choice and bg↔magnitude collinearity make this not a clean refutation of glucose — see §4.)*

### 3.7 The magnitude bias — is it a deployable ISF term?
`inv008/magnitude_bias.py`. `actual_drop ≈ 21 + 0.51·pred_drop`. **Affine** (curvature
pooled-significant but inconsistent across users, 70% same-sign, I²=96% → no universal saturation).
`m(pred_drop)` reclaims 24.2→19.1; **`m(iob)` adds 0.0 over pred_drop** (IOB is not the carrier).
**Platform-invariant** (Trio v5 affine applied to oref0 v7 and vice versa, MAE within ~0.6).

**Read:** the biggest lever (~5 mg/dL) is **not** an IOB-conditioned ISF — it is a roughly constant
insulin-action over-scale (b≈0.5; partly regression attenuation, so stated as an empirically
validated transferable correction, not literally "insulin 2× too strong") plus a baseline/
mean-reversion drift (intercept ~21, BG-driven). Both belong in the insulin-action model / √TDD level
/ basal, **not** the ISF formula.

### 3.8 Best-fit individualised model (the gradient fit)
`inv008/gradient_isf_fit.py`. `actual_drop ≈ a_u + s_u·(pred_drop·(100/BG)^k)`; per-user `(a_u,s_u)`
fit by within-user 5-fold CV; shared steepness `k` searched; **k=0 a special case so the data
chooses.**

| model | out-of-user MAE |
|---|---|
| static (no fit) | 24.2 |
| per-user scale+baseline, **k=0** | **17.0** |
| best-fit glucose **k\*=0.0** | 17.0 |
| Diabeloop shape, per-user scaled | 22.5 |
| cold-start population shape (no adaptation, best k=0.75) | 22.4 |

**Findings:** per-user level adaptation is the dominant lever (**+7.2 mg/dL**); the optimal shared
glucose steepness is **k=0** (glucose adds 0.0 beyond the level); cold-start mildly prefers a gentle
`k≈0.75` (worth +1.7 mg/dL) that should *fade* as the level is learned; **36% of users individually
prefer k≥1** but the per-user scale is **not predictable from TDD (ρ=0.03)** or profile ISF
(ρ=−0.17), so individualisation must be learned online.

### 3.9 Glucotoxicity literature (provided by Tim)
Ten papers establishing glucose→insulin-resistance mechanisms (Rabbani 2024 HK2; Khalid 2021
IRS-1/PI3K/Akt/GLUT4; Zhao 2023; Simon-Szabó 2024 ER-stress/mTOR; Galicia-Garcia 2020; Wang 2025
ectopic lipid; Beaupère 2021; Allocca 2025 SGLT2; Młynarska 2025; Yang & Sherman 2025). **Accepted:
the relationship is real.** Almost all of it is *chronic* glucotoxicity (days–years) → a *between-
person / slow-drift* property → captured in the per-user level / √TDD. The ISF *term* acts on
minute–hour excursions, where chronic resistance competes with acute clearance.

### 3.10 Tim's challenges and our concessions
Tim corrected several points, and we conceded:

1. **Diabeloop's control algorithm is refactored oref** → "different algorithm, doesn't transfer" is
   weakened. *Conceded.*
2. **The Diabeloop curve is the average of thousands of *in-user* curves**, not a cross-sectional
   fit → our "between-vs-within / ecological fallacy" argument is **wrong**. *Conceded fully.*
3. **Same method** (fasting windows; correction-dose drop vs static-ISF expectation) → results should
   be comparable, so a discrepancy is a real puzzle, not a sample-size issue.
4. **bg and `pred_drop` are collinear (r≈0.59), and bg-only ≈ magnitude-only.** Our gradient fit
   parameterised magnitude as the base and glucose as the increment, so "k=0" partly reflects that
   *ordering*. "Glucose adds nothing" was **too strong**. *Conceded.*

This narrowed the disagreement to one specific place: **at high glucose, our raw realised-ISF curve
rises (insulin stays effective) where the power law needs it to fall (resistance).**

### 3.11 Removing mean reversion (physiological)
`inv008/insulin_dose_response.py`. In a fasting T1D system glucose falls because insulin drives it
down; "mean reversion" of a high spike is not a statistical force. So we isolated windows where the
drop had to be *initiated* by insulin — glucose **flat or rising at entry** (`pre_slope`), not
already falling — and recomputed realised ÷ profile ISF by glucose.

| BG | falling (mean-reverting) | flat (clean) | rising (clean) |
|---|---|---|---|
| 100–120 | 0.59 | 0.84 | 0.65 |
| 145–175 | 0.84 | 0.98 | 0.90 |
| 205–260 | 0.95 | 0.93 | 0.91 |

**Removing mean reversion does NOT recover the power law.** Even the flat-entry (clean) curve is
flat-to-rising (implied k ≈ −0.1); insulin is ~fully effective at high BG (0.93–0.98) and *least*
effective near target (0.84, counterregulation, identical across strata). The faint dip at the very
top (0.98→0.93) is the only hint of resistance.

---

### 3.12 Marginal correction-dose response (test #2) — INCONCLUSIVE (important)
`inv008/dose_response_db.py` (model-free: regress realised drop on *delivered* correction insulin
`sug_smb_units`, controlling starting IOB, per glucose band; the slope = effective ISF; renal/
mass-action clearance ∝ BG falls into the intercept, not the slope). 89 users, 81,106 windows,
single-process.

**Result: the effective-ISF slopes came out negative (−0.4 to −1.7 mg/dL/U) in every band — physically
impossible.** Diagnosis (145–175 band): lo-correction windows dropped 58 mg/dL, hi-correction 38 —
*more correction insulin, less drop.* This is **confounding-by-indication / reverse causation**: in a
closed loop the correction dose is *reactive* (the controller keeps dosing precisely when glucose
isn't responding), so dose ⊥ outcome is violated and **a dose→drop relationship cannot identify ISF.**
The negative coefficient is the signature.

**Consequences:**
- Test #2 does **not** answer the resistance question. The method is broken for closed-loop data.
- **This caution generalises to any method that reads ISF from the loop's own doses — including a
  naive realised-ISF = drop/dose ratio vs BG.** Reactive dosing pushes drop/dose *down* exactly where
  doses are biggest (high BG), which would manufacture a power-law-shaped "resistance" curve
  *artefactually.* Whether the Diabeloop derivation controls for this is unknown to us — if it does,
  their curve stands; if it is a naive ratio, it is suspect. We cannot verify either way, so we do
  **not** claim their curve is an artefact — only that the same confound that broke our test could
  inflate a naive one.
- The `head_to_head` `realised_isf` is **less** exposed: it uses the loop's counterfactual
  *prediction* rescaled to candidate ISFs, not a dose-response slope. Its remaining confound is the
  non-insulin clearance (renal / glucose effectiveness), which test #1 targets and which does *not*
  depend on dose exogeneity.

**Net epistemic state:** the causal ISF-vs-glucose relationship is only *weakly identified* from
observational closed-loop data. Neither our flat curve nor the power law is confound-free — ours by
clearance, a naive dose/ratio method by reactive dosing. A fully clean answer likely needs a
fixed-dose natural experiment or an instrument. The robust, identification-clean step we *can* take is
test #1 (clearance subtraction on the prediction-based realised ISF).

### 3.13 Clearance-corrected ISF (test #1) — resistance IS real, but offset by clearance
`inv008/clearance_corrected_isf.py`. We estimate the insulin-INDEPENDENT flux *from the data* (the
4h drop in windows where the loop expects ~no insulin action, `|cgm − reason_IOBpredBG| < 5`), then
subtract it from the insulin-active windows and recompute the realised-ISF ratio.

| BG band | non-insulin flux (mg/dL/4h) | raw ratio | clearance-corrected ratio |
|---|---|---|---|
| 100–120 | 10 | 0.59 | 0.30 |
| 145–175 | 44 | 0.86 | 0.21 |
| 175–205 | 67 | 0.87 | 0.17 |
| 205–260 | 81 | 0.89 | 0.25 |

- The data-derived flux **rises with glucose and accelerates past ~180** — the renal-threshold
  signature, from data not literature constants.
- **Raw ratio is flat (k=−0.1); the clearance-corrected (insulin-only) ratio FALLS with glucose
  (k≈0.44).** Insulin per unit genuinely does less at high BG — **resistance is real.**
- The power law is *achievable*: bending the high-BG ratio to its level needs ~51–75 mg/dL of
  clearance; the data shows 67–81 — right where it sits.

**THE RECONCILIATION (the central result of this investigation):**
> At high glucose, **insulin resistance is real** (insulin does less per unit), **but insulin-
> independent clearance (renal + mass-action) rises with glucose too and roughly *cancels* it in the
> net observed drop.**

This makes every earlier result consistent: Diabeloop/the power law measured the *insulin* component
(resistance — physiologically correct); our gradient fit found flat `k=0` best for *prediction*
(because the loop predicts/doses against the **net**, where clearance offsets resistance); the
bridging experiment found the power law *degrades* prediction (it encodes resistance but **omits the
offsetting clearance — double-counting** and over-correcting highs). So it is not static-vs-dynamic:
the physiology is resistant, the net is flat, a flat effective ISF predicts the net, and the power
law is right about insulin but wrong as a predictor.

**Bound (honesty):** the corrected ratio is low everywhere (0.17–0.31), implying the loop over-states
pure insulin ~3–5× — a strong claim; the clearance estimate may over-subtract (absorb basal/
under-modelled insulin), so the resistance *magnitude* (k=0.44) is an **upper bound**. High-BG
low-insulin windows are sparse (n≈100–220). The cancellation is an *average* — it depends on renal
function / hydration / individual variation (e.g. **SGLT2-inhibitor users** have enhanced clearance →
net even flatter; **impaired-renal** users have less → net resistant), which is precisely why
**online per-user net adaptation** beats a fixed glucose curve.

## 4. What we have found (consolidated)

1. **A well-set static ISF predicts realised drops as well as the loop and far better than v1/v2.**
   The exploitable individualisation is the **per-user level/scale** (≈ +7 mg/dL), learned online;
   it is *not* predictable cold-start from TDD/profile.
2. **The TDD law is √TDD** (empirical −0.5), not v1's 1/TDD or v2's 1/TDD²; this carries the chronic
   glucotoxic between-person resistance.
3. **No deployable steep glucose-ISF curve survives on oref.** Transplanted Diabeloop/power-law
   shapes degrade prediction; the best-fit shared steepness is k=0 once the level is individualised.
4. **The dominant glucose-dependent effect is near target (counterregulation)** — large, consistent,
   and a *safety* signal (ease off; do not fit the low realised ISF). This is common ground with
   Diabeloop's conservative-low tail and dynISF's low-BG inflation.
5. **The biggest residual error source is correction *magnitude*** (the loop over-trusts large
   IOB-driven corrections): affine, platform-invariant, not an IOB-conditioned ISF — belongs in the
   insulin-action model / level / basal.
6. **Mean reversion does not explain our rising high-BG curve.** The high-BG resistance question now
   rests entirely on **insulin-independent clearance (renal + glucose effectiveness)**, not yet
   removed.

---

## 5. Counterfactuals — what would change the conclusion

- **If, after subtracting glucose effectiveness / renal clearance, the clean high-BG curve falls
  (k>0):** genuine glucotoxic resistance is present overnight → a *gentle* glucose-ISF term is
  warranted (still far gentler than k=3.5) → Tim's position is vindicated at the top end. *(Next test
  #1.)*
- **If the marginal correction-dose response (drop per unit of a discrete correction bolus) falls
  with glucose:** the dose-attributed insulin effect *is* resistant, and our total-trajectory metric
  was inflated by non-insulin flux → adopt a glucose term. *(Next test #2 — this section's run.)*
- **If daytime/postprandial windows show resistance** (the power-law found daytime k→4; oref here is
  overnight only): a glucose term may be regime-specific (daytime) even if overnight is flat.
- **If the per-user adaptive-k survives nested CV** (the 36% is currently selection-inflated): a
  per-user *learned* steepness is justified for a real minority.
- **If bg is parameterised as the base instead of magnitude** and still loses out-of-user once
  mean-reversion and clearance are removed: confirms magnitude (not glucose) as the carrier despite
  the collinearity.

Conversely, the conclusion **stands** if: the clean curve stays flat after clearance subtraction;
the dose-response is flat; and the per-user k gains nothing under nested CV.

---

## 6. Known confounds (and how each was handled)

| confound | effect | handling |
|---|---|---|
| between-user level (chronic resistance) | inflates any pooled glucose slope | per-user centring / scale; √TDD |
| min-predicted-drop selection | flips the near-target sign | identified; compared filtered vs unfiltered |
| correction magnitude (`pred_drop`) | masquerades as a high-BG effect (r=0.59 with bg) | joint control; magnitude-only control; affine fit |
| mean reversion (high spikes falling) | inflates high-BG realised ISF | entry-trajectory stratification (excluded) |
| **renal / glucose-effectiveness clearance (>~180)** | **inflates high-BG realised ISF, insulin-independent** | **NOT yet removed — open** |
| dawn / EGP | adds baseline drift | partly via hour, start/pre-slope; intercept term |
| residual carb absorption | inflates drops | window carb-screen (imperfect) |
| bg↔magnitude collinearity | makes "glucose vs magnitude" attribution order-dependent | acknowledged; transfer + form tests used instead of attribution |

---

## 7. Recommendation as it stands (provisional)

> **ISF = `s_u · K/√TDD · (100/BG)^{k}`  + near-target easing clamp**, where
> - `K/√TDD` = population cold-start level at BG 100,
> - `s_u` = per-user effective-sensitivity scale + baseline, **adapted online from the user's own
>   outcomes** (the +7 mg/dL lever; not cold-start-predictable),
> - `k` = glucose steepness, **scheduled** ≈0.75 at cold-start → 0 as `s_u` is learned (a Bayesian
>   shrinkage: trust the population prior until you know the person),
> - near-target easing clamp = a one-sided safety guardrail (raise ISF / ease correction approaching
>   target), motivated by counterregulation and hypo-risk asymmetry, **not** fitted to realised ISF.

**Why no glucose curve, restated with the resolution:** not "there is no resistance" (there is), but
"resistance is **offset by clearance** in the net the loop acts on, so the loop wants the flat *net*
effective ISF; the power law's resistance term double-counts and over-doses highs." Where clearance
stops tracking resistance (impaired renal, dehydration, SGLT2 users, the minority), the net is no
longer flat — which is the argument for **online adaptation of the net** rather than a fixed curve.

**Caveats:** overnight-fasting scope only (daytime/postprandial untested — clearance may not offset
there); the clearance estimate is an upper bound; a per-user *learned* `k` for the ~36% is possible
but unvalidated (needs nested CV).

---

## 8. Path forward

The science has reached a natural resting point: the high-BG question is resolved (resistance real,
offset by clearance → flat net), test #2's failure shows observational closed-loop data weakly
identifies the causal ISF, and the deployable shape is settled (per-user net level + clamp, no glucose
curve). The remaining work is to turn that into a validated, deployable algorithm and close the two
genuine scientific gaps. Priority order:

**Phase A — lock the spec (immediate).**
- Write the v-next recommendation doc (Dynamic-ISF set) incorporating the clearance↔resistance
  cancellation, √TDD level, online net-scale adaptation, near-target clamp, cold-start gentle-`k` fade.
- This memo + memory done; commit the four new scripts (`bridge_diabeloop`, `magnitude_bias`,
  `gradient_isf_fit`, `insulin_dose_response`, `dose_response_db`, `clearance_corrected_isf`).

**Phase B — build the deployable mechanism (the 7 mg/dL lever).**
- Specify + prototype the **online per-user net-ISF adaptation**: estimate the affine
  `actual ≈ a_u + s_u·pred_drop` from the user's own outcomes, shrink toward the `√TDD` prior, with
  safety bounds and the near-target clamp. This is the concrete engineering deliverable.
- Forward-sim / shadow-mode validation harness (reuse the prior `*_forward_sim` / `*_shadow_eval`
  infrastructure) — because observational identification is weak, any change is validated
  prospectively, not from history.

**Phase C — close the two real scientific gaps.**
- **Daytime / postprandial replication** — the highest-value open question: clearance may *not*
  offset resistance when carb absorption dominates (the power law steepened to k→4 by day). If a
  glucose term is ever warranted, it is here, not overnight. Run the clearance-corrected analysis on
  daytime windows.
- **Nested-CV diagnostic** for the per-user adaptive-`k` minority (is the 36% real out-of-sample?).

**Phase D — deploy.**
- Port the validated net-ISF + online adaptation to an oref instantiation (AAPS/Trio), behind a
  shadow flag first; monitor TIR / hypo. Flag renal-status modifiers (SGLT2, CKD) where the net
  cancellation breaks.

---

## 9. Artifacts

**Scripts** (`inv008/`, run `python -m inv008.<name>`): `head_to_head`, `err_common`, `err_curve`,
`err_consistency`, `isf_pattern_ml`, `bridge_diabeloop`, `magnitude_bias`, `gradient_isf_fit`,
`insulin_dose_response`. **Outputs:** `results/*.{json,md}`, `charts/inv008/fig_*.png`. **Data:**
`results/head_to_head_windows.parquet`. **Repo:** `github.com/tim2000s/dynamic-isf-calculations`
(committed through 199862e; `insulin_dose_response` + this memo pending). **Prior g(BG) work:** other
Drive, `Dynamic ISF data/` (Diabeloop, power-law, hybrid).
