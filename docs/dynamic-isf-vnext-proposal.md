# Dynamic ISF — a proposal for v-next

**2026-06-07** · Tim Street / Claude

---

## Context

Dynamic ISF sets correction sensitivity from total daily dose (TDD) and current glucose. The
original equation (v1, Chris Wilson) makes the sensitivity anchor inversely proportional to
TDD. A later revision (v2) makes it inversely proportional to TDD squared, a steeper
dependence.

Tested against sensitivity calculated directly from 171 people's own glucose and insulin data,
the steeper exponent moves the wrong way: the observed TDD dependence is shallower than even
v1. v2 is the worst-fitting of every equation we tried, and in practice it gives a weaker
correction than v1 for almost everyone (a median of 3×, most markedly at low glucose) while
over-estimating sensitivity more severely than v1.

This document proposes the next version, and answers the question the data raises: one fixed
equation is not enough on its own. The exponent of the TDD law can be universal, but the level
of the curve has to be set per person, because the same TDD predicts ISFs that differ around
ninefold between people and most of that is a stable per-person offset. The proposal is a
universal shape with a per-user constant, recalibrated weekly from the person's own recent
data once dynamic ISF is enabled.

---

## 1. Proposal

Set sensitivity at normal target in inverse proportion to the **square root of TDD**, with a
**per-user constant K**:

```
ISF at normal target = K_user / √TDD
```

`K_user / √TDD` is the **TDD term** — it replaces v1's `1800/TDD` and v2's `115000/TDD²`. ISF
at other glucose levels comes from a glucose curve `g(BG)` that **replaces v1/v2's logarithmic
scaler** with the Diabeloop clinical curve (§4). In the same form as v1 and v2:

```
v1:      ISF(BG) = 1800    / ( TDD  · ln(BG_capped/75 + 1) )
v2:      ISF(BG) = 115000  / ( TDD² · ln(BG_floored/75)    )
v-next:  ISF(BG) = ( K_user / √TDD ) · g(BG)
```

The exponent (−½) is universal; `K_user` is calibrated per person (§5–6); `g(BG)` is the same
for everyone. The TDD *blend* and the high-glucose cap are kept from v1/v2; the *glucose
scaler* changes from log to `g(BG)`.

### The implementation equation

Substituting the safe-default (Tier-1) constant `K_user = profile_ISF · √(median 14-day TDD)`
and the Diabeloop quartic, the single closed-form ISF a device computes each cycle is:

```
                                        272 − 3.121·BG + 0.01511·BG² − 3.305e-5·BG³ + 2.69e-8·BG⁴
ISF(BG) = profile_ISF · √(TDD₁₄ / TDD) · ──────────────────────────────────────────────────────
                                                                81.63
```

where
- `profile_ISF` — the person's existing static profile ISF (mg/dL per U)
- `TDD₁₄` — median blended TDD over the last 14 days (updated weekly)
- `TDD` — the current blended TDD (the existing v1/v2 5-window blend, unchanged)
- `BG` — current glucose, clamped to `[54, cap(210, excess/3)]`
- `81.63 = q(99)` — the quartic at the normal target, so the glucose factor is 1.0 at target
- **level floor (§8.2):** clamp `profile_ISF · √(TDD₁₄/TDD) ≥ profile_ISF / 1.5` before applying
  the glucose factor (bounds how far the *level* can strengthen; does not touch `g(BG)`)

At the user's median TDD and target glucose this returns their **existing profile ISF
exactly** (`√(TDD₁₄/TDD)=1`, glucose factor `=1`) — a behaviour-preserving generalisation of
their current setting, adding a √TDD level response and the glucose curve. The Tier-2 variant
swaps `profile_ISF` for the person's measured sensitivity (§6).

- **Exponent: universal −½**, robust across the TDD construct used (§2).
- **K_user: per-user, recalibrated weekly** from the person's own recent data (§6). The safe
  default (Tier 1) anchors K to the user's existing profile ISF, preserving their dosing
  *level*; an optional stronger setting (Tier 2) anchors K to measured sensitivity (needs
  validation).
- **g(BG): the Diabeloop clinical curve** (§4), shared by all users — falls with glucose
  (firmer corrections when high, protective when low).
- **Evaluate in shadow first** (§8) with the level clamp before any live dosing.

The earlier question of a single global constant (the cohort values were ≈355 against tuned
profiles and ≈145 against measured sensitivity) is **superseded**: no global K is adequate
(§5), so the level is set per user rather than chosen once for everyone.

---

## 2. Why a square-root law

Every TDD~ISF relationship available was scored by leave-one-user-out cross-validation
against two independent targets: sensitivity calculated from each person's own data
(n = 114) and tuned-profile ISF (n = 138). The √TDD law wins on both.

![Candidate equations under leave-one-user-out cross-validation — median log error per target (green = fitted, blue/red/grey = fixed rules)](charts/inv008/fig_best_fit.png)

**Target: tuned-profile ISF (n = 138)**

| candidate | median abs err | median log err | within ±30% |
|---|---|---|---|
| **K/√TDD (K = 355)** | **12.8** | **0.256** | **51%** |
| Power law A·TDD^b (free exponent) | 13.5 | 0.268 | 49% |
| TDD-quartile bands | 13.8 | 0.292 | 48% |
| Multivariate (TDD, CR, basal, target) | 14.7 | 0.306 | 46% |
| **v1 (TDD⁻¹)** | 17.0 | 0.324 | 38% |
| 1700-rule | 17.9 | 0.454 | 34% |
| **v2 (TDD⁻²)** | 124.0 | 1.304 | 7% |

**Target: calculated sensitivity (n = 114)**

| candidate | median abs err | median log err | within ±30% |
|---|---|---|---|
| Power law + basal fraction | 6.0 | 0.294 | 44% |
| **K/√TDD (K = 145)** | 6.2 | 0.297 | 45% |
| Power law A·TDD^b | 6.1 | 0.304 | 45% |
| 1700-rule | 16.2 | 0.605 | 15% |
| **v1 (TDD⁻¹)** | 26.0 | 0.813 | 7% |
| **v2 (TDD⁻²)** | 171.1 | 2.323 | 2% |

(v2's anchor sits ~3× above v1 at target, so as a between-person level predictor it is far
the worst — its strength is low-glucose protection, not the level.)

Four points decide it:

1. **The exponent the data wants is ≈ −0.5.** Free-fitted exponents bracket it: on the
   cross-sectional treatments-per-day TDD they are −0.43 (profiles) and −0.38 (calculated
   sensitivity); re-fit on each person's median *blended per-tick* TDD — the windowed
   quantity the equation actually consumes, and an independent construct — they steepen to
   −0.62 (profiles) and −0.55 (calculated sensitivity). Either way the confidence intervals
   exclude −1, and −0.5 sits in the middle. Fixing the exponent at −½ costs nothing
   measurable; if anything the blended construct argues for a touch steeper, which a future
   revision could revisit.
2. **Extra inputs do not help.** Carb ratio and target add noise; basal fraction is
   marginal and not robust. The blend weight on the tuned profile value fitted to **zero**
   when predicting calculated sensitivity — a person's tuned ISF carries no information
   about their actual sensitivity beyond what TDD already provides.
3. **It beats both v1 and v2 decisively** — and v2 is last on both targets.
4. **It reproduces on an independent TDD construct.** The same √TDD law, re-fit on median
   blended per-tick TDD instead of the cross-sectional scalar, again scores best
   (median log-error 0.24 profiles / 0.29 sensitivity) — so the result is not an artefact
   of one TDD definition.

---

## 3. What it does, relative to v1 and v2

All three differ in the TDD term; v2 and v-next also differ from v1 in the glucose term:

```
v1:      ∝ 1 / TDD          glucose: ln(BG/div + 1)
v2:      ∝ 1 / TDD²         glucose: ln(BG/div)  (no +1; BG floored at div+1)
v-next:  = K_user / √TDD    glucose: Diabeloop quartic g(BG)  (§4)
```

ISF at normal target (where `g(BG) = 1`), mg/dL per U, divisor 75. The v-next column uses the
cohort-representative K = 355 purely to show the **level shape** — the shipped equation uses a
per-user K_user (§1, §6), not this single value:

| TDD (U/day) | v1 | v2 | **v-next (355/√TDD)** | v-next vs v1 |
|---|---|---|---|---|
| 15 | 143 | 1840 | 92 | ~1.6× stronger corrections |
| 25 | 86 | 663 | 71 | ~1.2× stronger |
| 36 | 59 | 320 | 59 | equal |
| 50 | 43 | 166 | 50 | ~15% weaker |
| 80 | 27 | 65 | 40 | ~1.5× weaker |
| 120 | 18 | 29 | 32 | ~1.8× weaker |

(v2's at-target ISF is high and rises steeply *below* target, where `ln(BG/div)→0` gives
near-zero corrections. So v2 is far gentler than v1 everywhere and strongly hypo-protective
when low; the v1/v2 ratio is glucose-dependent. See the v1-vs-v2 analysis doc.)

(ISF in mg/dL per U, at normal target.) v-next is gentler than v1 for heavy insulin users
and stronger for light ones — the opposite tilt to v2, and the direction the data supports.

---

## 4. The glucose curve g(BG)

v-next's glucose term is a **power-law / Diabeloop curve** — ISF falling with glucose (more
insulin per mg/dL when high; strongly protective at low). Its provenance is deliberately
different from the TDD level:

- **The shape and exponent come from the Diabeloop clinical population model** (ISF vs
  glucose from controlled clinical data, piecewise polynomial / quartic). That is the basis
  for the glucose curve; v-next adopts it rather than fitting an exponent from device data.
- **The absolute exponent cannot be determined from observational AID data.** Prediction-
  error backtests on closed-loop data estimate a candidate curve's curvature only *relative
  to the loop's existing DynISF* (the device already applies a glucose-dependent ISF), are
  horizon-dependent, and are underpowered — the glucose exponent moves prediction MAE by
  ≤2 mg/dL and a power-law sits within noise of the log scaler. So those backtests can
  confirm a power-law is *no worse* than log, but they cannot set or validate the exponent.
  (This is the glucose-axis counterpart of the level result: observational data, confounded
  by the loop's own behaviour and by glucose mean-reversion, cannot recover g(BG).)
- **Validation is prospective, not retrospective.** Establishing the exponent for this
  population requires a prospective / closed-loop trial; retrospective fitting on AID logs
  does not answer it.
- One robust observational note: at the full-action horizon the loop's own ISF predictions
  are essentially unbiased (≈ −0.8 mg/dL), so the much larger positive bias seen at a 2-hour
  horizon is a horizon artefact (insulin unfinished), not real over-aggression.

**Concrete instantiation.** The curve adopted is the Diabeloop population quartic, normalised
to 1.0 at the normal target so it composes with the level without moving the per-user anchor:

```
q(BG) = 272 − 3.121·BG + 0.01511·BG² − 3.305e-5·BG³ + 2.69e-8·BG⁴
g(BG) = q(BG) / q(target)        # BG high-capped at 210 (excess/3), low-floored at 54
```

It falls monotonically across the physiological range — more insulin per mg/dL when high,
strongly protective when low — and corresponds to a mild power law of exponent ≈ 1.3 over
BG 70–250, far gentler than the exponents retrospective device fits suggested (consistent
with those being artefacts). A one-parameter `(target/BG)^k` form is available as an
alternative of the same family.

So the glucose dimension is settled in **direction and form** (power-law, hypo-protective,
from the Diabeloop clinical model) and **open in exponent**, which is taken from that clinical
model and flagged for prospective validation — unlike the √TDD level and per-user K, which
the cohort data establishes directly.

---

## 5. Why a single constant is not enough — the level must be per-user

The √TDD *shape* fits everyone; the *level* does not. Decomposing ISF variance across the
cohort settles which parts of the equation can be ubiquitous and which cannot:

- **Exponent → universal.** Within one person, TDD barely moves — median p90/p10 ≈ 1.8×
  over their whole record (far less inside any 14-day window) — against ~15× between people.
  There is no lever arm to estimate a per-user exponent from one person's data, and the
  population evidence supports a single value. The exponent must be global.
- **Constant K → per-user, and necessarily so.** With the best global K, the typical person
  is still off by ~57% (between-user residual SD 0.45 in log; per-user constants span 35–325,
  a 9× range). Crucially, **~84–98% of that error is a stable per-person offset**, not noise:
  the scatter in *measuring* a person's own sensitivity is only ~6% (full history), rising to
  perhaps ~15–18% for a 14-day window — far smaller than the 57% a global constant leaves. So
  a short window of the person's own data removes most of the error a one-size equation makes.

This is why the earlier "which global K — 355 or 145?" framing was the wrong question. 355
(tuned-profile) and 145 (measured-sensitivity) are just the cohort medians of a quantity that
varies 9× between individuals. The level is not one number to choose; it is a per-user
constant to measure.

---

## 6. Personalisation: how K_user is set

`K_user` is recalibrated **weekly** from the person's own recent data; the universal √TDD
term then provides the within-week response as TDD moves. Two tiers, in increasing
aggressiveness:

**Tier 1 — profile-anchored (safe default).**
```
K_user = profile_ISF × √(median TDD over the last 14 days)
```
This re-expresses the user's *existing* static ISF as a TDD-responsive curve. At their typical
TDD **and target glucose** the ISF level is unchanged, so the dosing *level* is preserved; the
new behaviour is the glucose curve (§4), which firms corrections when high and eases them when
low, plus the √TDD adjustment as TDD drifts. The cohort shadow run (§7) confirms the level is
held — the §8.2 level clamp binds for ~0% of readings — while g(BG) reshapes dosing across the
glucose range. It needs nothing more than 14 days of TDD history.

**Tier 2 — sensitivity-anchored (stronger; needs validation).**
```
K_user = measured_ISF × √(median TDD over the last 14 days)
```
where `measured_ISF` is the person's empirically observed sensitivity (a ΔIOB regression of
glucose change on insulin absorbed over fasting windows). This doses to measured insulin
effect, which across the cohort is ~2.4× stronger than tuned profiles and carries the same
"is the measured level biased low / is it safe?" caveat as any move toward true sensitivity.
It is gated behind shadow evaluation and forward validation.

**Robustness requirements (both tiers):**
- **Fit-quality gate.** Per-user sensitivity fits are often weak over short windows (cohort
  median R² ≈ 0.22). When a 14-day estimate is too noisy (low R² / wide CI), fall back to
  Tier 1 rather than act on a noisy number.
- **Change clamp.** Bound the week-to-week move in K_user (and the deviation from profile ISF)
  so a single noisy window cannot lurch the curve.
- **Cadence.** Weekly is appropriate: K tracks slow sensitivity drift (season, activity,
  illness, life stage); the √TDD term already handles faster TDD swings within the week.

**Data the device needs.** Tier 1 needs only 14 days of TDD (delivery history AAPS already
holds). Tier 2 additionally needs IOB, glucose and treatment records retained over the window
to run the sensitivity regression — pulled from Nightscout, or persisted by AAPS if not
already stored.

---

## 7. Evidence status

**Settled:**
- The equation implementations are exact (18 unit tests) and reproduce device-logged ISF
  for all dynamic-ISF users, to within unmodelled per-person settings.
- The TDD exponent (≈ −0.5) is tested against sensitivity *calculated from people's own
  glucose and insulin data*, not against profile settings or the equations themselves; it is
  robust across the TDD construct (cross-sectional and blended per-tick) and to outlier and
  duration checks.
- The level varies ~9× between people and is dominated by a stable per-person offset, so
  per-user calibration is well-founded, not a tuning convenience.
- A cohort shadow evaluation (138 people, 9.1M readings; counterfactual replay of the ISF
  each device would have used) confirms Tier-1 preserves the dosing *level* — the §8.2 level
  clamp would bind for a median of 0% of readings (worst person ~19%) — and that the dosing
  change it introduces (median 38% per correction) is the glucose curve, not a level shift.
  Across the glucose range it tracks today's DynISF (v1) within ±30% (median ratio 1.01),
  a little firmer at high glucose where the Diabeloop curve is steeper than v1's log scaler.

**Provisional:**
- The **Tier-2 sensitivity anchor** rests on a regression estimator that may be biased low by
  unrecorded carbohydrate or endogenous-glucose effects, and is weak per-user over short
  windows (median R² ≈ 0.22). Its *shape* is trustworthy; its absolute *level* is not yet. A
  cohort shadow run (114 people with a usable fit) puts measured sensitivity at a median 0.39×
  the profile ISF, so unclamped Tier-2 would dose about 2.6× the correction insulin of the
  person's current setting. Under the §8.2 level clamp that level is pulled back to the
  profile/1.5 ceiling for essentially everyone (the clamp binds on ~100% of readings), which
  means the clamped form barely expresses the measured anchor at all — and the data-derived
  study shows that anchor is hypo-biased. Tier-2 is therefore not deployable without forward,
  outcome-based validation and, most likely, a different level bound than the Tier-1 clamp.
- The cohort is largely composed of users who never ran dynamic ISF (≈87%); the small
  dynamic-ISF subgroup hints at a steeper slope, with confidence intervals too wide to resolve.
- All evidence is retrospective and decision-level (counterfactual replay), single-cohort
  (open-source AID, mostly 2016–2023), with no closed-loop outcomes.

**Implication:** the universal √TDD shape and the Tier-1 (profile-anchored) personalisation
can go to shadow evaluation now; the Tier-2 (sensitivity-anchored) level needs forward
validation before it doses.

---

## 8. Path to deployment

1. **Shadow evaluation (done retrospectively; repeat live).** A counterfactual cohort sweep
   (§7) has computed the Tier-1 v-next ISF against each person's real per-tick history: the
   level is preserved and the change is the glucose curve, tracking v1 within ±30%. The live
   step is to compute it on-device alongside the running equation, log without acting, and
   confirm the weekly K_user recalibration is stable in production.
2. **Low-TDD safety clamp — on the level term.** Clamp the **level** `K_user/√TDD` so it is
   never more than ~1.5× stronger than the user's profile-ISF value (a floor on the level),
   applied to the level *only*, not the full ISF — otherwise the intended high-glucose
   aggression of g(BG) would be clipped ~40% of the time. The shadow run shows this level
   clamp rarely binds under Tier-1 (≈0% of readings); it is a guard against TDD spikes and
   the Tier-2 sensitivity anchor rather than a routine limiter.
3. **Live trial (after shadow review).** Enable for opt-in testers at **Tier 1**, clamp
   active, shadow comparison still logging.
4. **Tier 2 (later).** Only with forward-validated outcomes should the sensitivity-anchored
   K be offered, and then with the fit-quality gate and change clamp in force.

---

## 9. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Light insulin users dosed too strongly | Low-TDD clamp / threshold gate (§8.2); shadow review before live |
| Per-user K calibrated on a noisy 14-day window | Fit-quality gate → fall back to Tier 1; week-to-week change clamp |
| Sensitivity (Tier-2) anchor over-aggressive | Not deployed by default; Tier 1 leaves average dosing unchanged; Tier 2 gated behind validation |
| TDD reconstruction error feeds the equation | √TDD is *less* TDD-sensitive than v1 and far less than v2, so the same TDD error moves ISF less — a robustness gain |
| Single-cohort generalisation; few true DynISF users | Trial as opt-in; monitor across TDD bands; revisit exponent if the DynISF population shows steeper |

---

## 10. The change, concretely

The TDD blend, glucose cap, autosensitivity and temp-target handling are untouched. Three
pieces change: a weekly per-user K, the anchor formula, and the glucose scaler (log → the
Diabeloop quartic g(BG)).

```
# weekly, per user (Tier 1 shown; Tier 2 swaps profile_ISF for measured_ISF):
K_user = profile_ISF * sqrt(median_TDD_14d)        // clamp vs previous K_user

# every cycle:
level   = K_user / sqrt(TDD)             // the TDD term; was 1800/TDD (v1) or 115000/TDD² (v2)
level   = max(level, profile_ISF / 1.5)  // §8.2 clamp on the LEVEL only (not the full ISF)
ISF(BG) = level * g(BG)                   // g(BG) = Diabeloop quartic, normalised to 1 at target
```

Because `K_user = profile_ISF · √(median TDD)`, at the user's typical TDD this returns their
profile ISF exactly — Tier 1 is a behaviour-preserving generalisation of their current
setting, with the √TDD response added.

---

## Appendix — artefacts

- Methodology: companion methodology paper
- v1 vs v2 analysis: companion analysis document
- Equation search: `fit_best_isf.py`, `results/best_isf_fit_results.{json,md}`
- Personalisation analysis (blended-TDD refit, variance decomposition): `inv008/fit_personalisation.py`
- Glucose curve + v-next equation: `inv008/dynisf.py` (`g_quartic`, `isf_vnext`, `k_user_tier1`), tests in `inv008/tests/`
- Cohort shadow evaluation (Tier 1): `inv008/phase11_shadow_eval.py` → `results/phase11_shadow_eval.{json,md}`, `charts/inv008/fig_shadow_eval.png`
- Cohort shadow evaluation (Tier 2): `inv008/phase12_shadow_eval_tier2.py` → `results/phase12_shadow_eval_tier2.{json,md}`, `charts/inv008/fig_shadow_eval_tier2.png`
- Comparison figure: `charts/inv008/fig_best_fit.png`
- Device validation: `inv008/validate_device_isf.py`, `results/device_isf_validation.{json,md}`
- Repository: `github.com/tim2000s/dynamic-isf-calculations`
