# Dynamic ISF — a proposal for v-next

**2026-06-07** · Tim Street / Claude

---

## Context

Dynamic ISF sets correction sensitivity from total daily dose (TDD) and current glucose.
The original equation (**v1**, Chris Wilson) makes the sensitivity anchor inversely
proportional to TDD. A later revision of the maths (**v2**) makes it inversely proportional
to TDD squared — a steeper dependence of sensitivity on TDD.

Testing both against sensitivity calculated directly from 171 people's own glucose and
insulin data shows the steeper exponent moved the wrong way: the observed TDD dependence is
*shallower* than even v1. v2 is the worst-fitting of every equation tested against
calculated sensitivity, and in practice it weakens corrections for the 77% of people below
~64 U/day while over-estimating their sensitivity more severely than v1. This document
proposes the next version of the equation.

---

## 1. Proposal

Set the sensitivity at normal target in inverse proportion to the **square root of TDD**:

```
ISF at normal target = 355 / √TDD          ("the 355 rule")
```

This is the anchor. The existing glucose scaler then sets ISF at every other glucose
level, unchanged — so the full equation, with the TDD blend and glucose cap also unchanged,
is:

```
ISF(BG) = (355 / √TDD) · ln(target/divisor + 1) / ln(bg_capped/divisor + 1)
```

**The 355 constant carries no insulin divisor.** It was fit directly to ISF (mg/dL per U)
against √TDD, so it states the sensitivity at normal target for the cohort as a whole,
independent of insulin type. The divisor — which varies with insulin peak time (≈75 for
Lyumjev, 65 for Fiasp, 55 for a standard rapid analogue) — enters *only* through the
glucose scaler, exactly as it does in v1 and v2; at normal target the scaler is 1 for every
insulin type, so 355/√TDD is the ISF-at-target regardless of insulin. (Avoid folding the
constant and the divisor's normal-target log term into a single number such as a "300
rule": that product is only valid for one divisor and silently changes the implied ISF for
other insulins.)

- **Anchor constant K = 355**, matched to the sensitivity experienced users tune to.
- **Evaluate in shadow first** (§6); enable for live dosing only after real-world
  divergence data and a low-TDD safety clamp are in place.

A more physiologically-faithful constant (K = 145, matched to *measured* sensitivity) is
**not proposed for deployment yet**: it doses about 2.4× more strongly across the board and
rests on a provisional estimate (§5).

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
| **v2 (TDD⁻²)** | 33.1 | 0.566 | 28% |

**Target: calculated sensitivity (n = 114)**

| candidate | median abs err | median log err | within ±30% |
|---|---|---|---|
| Power law + basal fraction | 6.0 | 0.294 | 44% |
| **K/√TDD (K = 145)** | 6.2 | 0.297 | 45% |
| Power law A·TDD^b | 6.1 | 0.304 | 45% |
| 1700-rule | 16.2 | 0.605 | 15% |
| **v1 (TDD⁻¹)** | 26.0 | 0.813 | 7% |
| **v2 (TDD⁻²)** | 44.6 | 1.229 | 7% |

Four points decide it:

1. **The exponent the data wants is ≈ −0.5.** Free-fitted exponents land at −0.43
   (profiles) and −0.38 (calculated sensitivity); the log-log fit of calculated
   sensitivity against TDD gives −0.4 to −0.56, with confidence intervals excluding −1.
   Fixing the exponent at exactly −0.5 costs nothing measurable and gives a clean,
   memorable form.
2. **Extra inputs do not help.** Carb ratio and target add noise; basal fraction is
   marginal and not robust. The blend weight on the tuned profile value fitted to **zero**
   when predicting calculated sensitivity — a person's tuned ISF carries no information
   about their actual sensitivity beyond what TDD already provides.
3. **It beats both v1 and v2 decisively** — and v2 is last on both targets.
4. **It is corroborated** by an independent TDD-band analysis of the same cohort, whose
   band/log-linear family sits in the same performance cluster; √TDD is its continuous
   form.

---

## 3. What it does, relative to v1 and v2

All three equations share the identical glucose scaler; they differ only in how TDD sets
the anchor. With the divisor shown as a parameter (≈75 Lyumjev / 65 Fiasp / 55 rapid):

```
v1:      ISF(BG) = 1800   / ( TDD  · ln(bg_capped/divisor + 1) )
v2:      ISF(BG) = 115000 / ( TDD² · ln(bg_capped/divisor + 1) )
v-next:  ISF(BG) = (355 / √TDD) · ln(target/divisor + 1) / ln(bg_capped/divisor + 1)
```

Note a structural difference the divisor exposes: v1 and v2 fold the normal-target log term
into their leading constants, so their ISF *at normal target* shifts with insulin type;
v-next's anchor, 355/√TDD, does not — the divisor in v-next acts only on the glucose-scaling
shape above target.

The proposed curve rotates v1 about TDD ≈ 36 U/day (ISF at normal target, mg/dL per U,
shown for divisor 75):

| TDD (U/day) | v1 | v2 | **v-next (355/√TDD)** | v-next vs v1 |
|---|---|---|---|---|
| 15 | 143 | 607 | 92 | ~1.6× stronger corrections |
| 25 | 86 | 219 | 71 | ~1.2× stronger |
| 36 | 59 | 105 | 59 | equal |
| 50 | 43 | 55 | 50 | ~15% weaker |
| 80 | 27 | 21 | 40 | ~1.5× weaker |
| 120 | 18 | 9.5 | 32 | ~1.8× weaker |

(The v2 column tracks v1 by exactly the 63.9/TDD ratio — far higher than v1 below the
64 U/day crossover, far lower above it. At 15 U/day v2 implies an ISF of ~600 mg/dL per
unit, i.e. almost no correction, which is the extreme over-estimation the data rejects.)

(ISF in mg/dL per U, at normal target.) v-next is gentler than v1 for heavy insulin users
and stronger for light ones — the opposite tilt to v2, and the direction the data supports.

---

## 4. Shape versus level: the constant is a safety decision

The √TDD *shape* fits both targets. The two natural constants differ by 355/145 ≈ 2.45 —
exactly the ratio between tuned-profile ISF and calculated sensitivity measured separately.
So the equation splits into:

- a **statistical question** — what exponent? — answered by the data: −0.5; and
- a **safety question** — what constant? — which the data alone cannot settle:
  - **K = 355** reproduces the sensitivity experienced users have tuned themselves to. Safe
    to trial because it is anchored to dosing strengths people already run.
  - **K = 145** reproduces measured insulin effect — stronger, "truer" dosing, but in
    territory no one in the cohort actually operates at.

We propose **K = 355** for any first deployment.

---

## 5. Evidence status

**Settled:**
- The equation implementations are exact (18 unit tests) and reproduce device-logged ISF
  for all dynamic-ISF users, to within unmodelled per-person settings.
- The TDD exponent (≈ −0.5) is tested against sensitivity *calculated from people's own
  glucose and insulin data*, not against profile settings or the equations themselves, and
  is robust across cohorts and to outlier and duration sensitivity checks.

**Provisional:**
- The **calculated-sensitivity constant (K = 145)** comes from a regression estimator that
  may be biased low by unrecorded carbohydrate or endogenous-glucose effects. It has
  confidence intervals but no external ground truth. The *shape* it implies is trustworthy;
  the *level* is not yet.
- All evidence is retrospective and decision-level (counterfactual replay), single-cohort
  (open-source AID, mostly 2016–2023), with no closed-loop outcomes.

**Implication:** K = 355 can go to shadow evaluation now; K = 145 needs forward validation
before it is even a candidate for dosing.

---

## 6. Path to deployment

1. **Shadow evaluation (now).** Compute the v-next ISF alongside the live equation and log
   it, without acting on it, for several weeks. Measure how often, and by how much, the
   proposed ISF would have changed dosing, stratified by TDD band. No dosing impact.
2. **Low-TDD safety clamp.** Below ~36 U/day the proposed curve doses more strongly than v1
   (≈1.6× at 15 U/day) — relevant for children and very insulin-sensitive adults. Before any
   live dosing, clamp the result so it is never more than ~1.5× stronger than the v1 value,
   or gate the new law above a TDD threshold. The data says v1 *over-estimates* ISF for
   these people, so the clamp is a conservatism, not a correction.
3. **Live trial (after shadow review).** Enable for opt-in testers with the clamp active
   and the shadow comparison still logging, at K = 355.
4. **Revisit the constant (later).** Only with forward-validated outcome data should a
   constant between 355 and 145 be considered.

---

## 7. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Light insulin users dosed too strongly | Low-TDD clamp / threshold gate (§6.2); shadow review before live |
| Calculated-sensitivity constant over-aggressive | Not deployed; K = 355 only until forward validation |
| TDD reconstruction error feeds the equation | √TDD is *less* TDD-sensitive than v1 and far less than v2, so the same TDD error moves ISF less — a robustness gain |
| Single-cohort generalisation | Trial as opt-in; monitor across TDD bands; revisit with broader data |
| Insulin-type / divisor variation | The 355 anchor is divisor-free and applies at normal target for any insulin; the divisor enters only through the unchanged glucose scaler. Whether the anchor *should* also vary with insulin peak (as v1/v2 implicitly make it) is not resolvable from this cohort and is a question for the shadow evaluation |

---

## 8. The change, concretely

The only change is the sensitivity-anchor computation; the TDD blend, glucose cap, glucose
scaler, autosensitivity, and temp-target handling are untouched:

```
anchor = 355 / sqrt(TDD)        // ISF at normal target; no divisor term
                                // was: 1800 / (TDD · ln(target/divisor + 1))   for v1
                                //  or: 2300 / (ln(target/divisor + 1) · TDD² · 0.02)  for v2
ISF(BG) = anchor · scaler       // scaler is unchanged and equals 1 at normal target
```

The anchor is a plain `355 / sqrt(TDD)` — unlike v1 and v2, it does not multiply the
normal-target log term into the constant, so the insulin divisor never touches the anchor.
The divisor continues to act only inside the unchanged glucose `scaler`.

---

## Appendix — artefacts

- Methodology: companion methodology paper
- v1 vs v2 analysis: companion analysis document
- Equation search: `fit_best_isf.py`, `results/best_isf_fit_results.{json,md}`
- Comparison figure: `charts/inv008/fig_best_fit.png`
- Device validation: `inv008/validate_device_isf.py`, `results/device_isf_validation.{json,md}`
- Repository: `github.com/tim2000s/dynamic-isf-calculations`
