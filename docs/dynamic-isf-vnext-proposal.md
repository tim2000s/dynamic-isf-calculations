# Dynamic ISF: withdrawn v-next proposal and current evidence status

**2026-06-07; audited 2026-09-16** · Tim Street / Claude

> This document records a model proposed in June 2026 and the evidence that led to its
> withdrawal. It is not a dosing recommendation. The full technical review is in
> `DYNAMIC-ISF-AUDIT-2026-09.md`.

---

## Summary

The proposed v-next model combined three ideas:

1. an ISF anchor proportional to `1/√TDD`;
2. a per-person constant, `K_user`;
3. a glucose curve adapted from a Diabeloop quartic.

The model fitted two historical targets better than Dynamic ISF v1 or v2. One target was the
person's entered profile ISF. The other was a coefficient calculated from glucose change and
change in IOB during fasting windows.

The audit changes the interpretation. The change-in-IOB coefficient is an observational outcome
proxy, not an independent measurement of insulin action. New insulin is delivered during the
window, the controller chooses that delivery in response to glucose, and unrecorded carbohydrate
and endogenous glucose remain possible. The fitted `1/√TDD` relationship therefore does not
establish a physiological TDD law. The same problem affects constants derived from that proxy.

The proposal is withdrawn for live dosing. It can still be used as a named candidate in
simulation, shadow replay or a prospective study designed to identify insulin effect.

---

## 1. The historical model

At normal target, the proposed sensitivity was:

```
ISF = K_user / √TDD
```

At other glucose levels:

```
v1:      ISF(BG) = 1800   / (TDD  · ln(BG_capped/75 + 1))
v2:      ISF(BG) = 115000 / (TDD² · ln(BG_floored/75))
v-next:  ISF(BG) = (K_user / √TDD) · g(BG)
```

For the standard rapid-acting configuration, including NovoRapid, v2 uses divisor 75 and floors
glucose at 76 mg/dL. The author-confirmed v2 equation does not contain `+1` in its logarithm.

The proposed glucose function was:

```
q(BG) = 272 − 3.121·BG + 0.01511·BG² − 3.305e-5·BG³ + 2.69e-8·BG⁴
g(BG) = q(BG) / q(99)
```

The Tier 1 form anchored `K_user` to the existing profile ISF. The Tier 2 form replaced the
profile anchor with the change-in-IOB proxy. Tier 2 is invalid as a physiological calibration.
Tier 1 preserves the person's entered level at their reference TDD, but the added TDD and glucose
responses remain unvalidated.

---

## 2. What the historical fit showed

The candidate equations were scored by leave-one-user-out cross-validation against entered
profile ISF in 138 people and the change-in-IOB proxy in 114 people.

![Historical candidate comparison. These are predictive fits to entered settings and an observational proxy, not physiological validation.](charts/inv008/fig_best_fit.png)

| target | best simple candidate | median log error | interpretation |
|---|---|---:|---|
| entered profile ISF | `355/√TDD` | 0.256 | predicts how this cohort had configured ISF |
| change-in-IOB proxy | approximately `145/√TDD` | 0.297 | predicts an observational closed-loop coefficient |

These are valid predictive descriptions of those targets. They are not physiological constants.
The apparent exponent near `−0.5` may reflect controller behaviour, treatment patterns, selection
of fasting windows, residual carbohydrate and the construction of the proxy itself.

V2 performed poorly against both targets. That comparison remains useful for understanding how
different the equations are, but it does not prove that the square-root candidate is the correct
biological relationship.

---

## 3. What survives the audit

### Equation arithmetic

The implementations of v1, v2 and the candidate v-next equation are reproducible. The v2 replay
already used the author-confirmed no-`+1` form. Unit tests now include explicit checks that the
v2 floor follows the configured divisor, 75/76 for standard rapid-acting insulin and 55/56 for
the alternative configuration.

### Direct equation comparison

V2 generally produces a much higher ISF than v1 in this cohort, especially near its glucose
floor. The v2-to-v1 ratio depends on both glucose and TDD. This is a mathematical comparison of
equation outputs and remains valid.

### Predictive comparison within the loop model

In the selected overnight windows, a tuned static ISF and the loop's own ISF prediction produced
smaller prediction errors than v1 or v2. This is a comparison within the loop's linear IOB
prediction model. It is not an independent measurement of physiological sensitivity.

### Glucose-shape result

After changing the shape analysis from shuffled folds to contiguous time folds, the best fitted
power exponent remained `k = 0`. Adding a glucose-dependent multiplier did not improve median
absolute error. The analysed windows therefore provide no predictive support for adding the
proposed glucose curve. They do not prove that glucose can never affect insulin sensitivity.

---

## 4. Why the dosing proposal is withdrawn

Four identification problems matter.

1. **Change in IOB is not absorbed insulin.** IOB changes because insulin acts, because new
   insulin is delivered, and because time passes through the assumed action curve.
2. **Delivery is endogenous.** A closed-loop system gives more or less insulin in response to
   glucose and trend. Insulin delivery is therefore correlated with the outcome being explained.
3. **Carbohydrate screening is incomplete.** A future-rise filter uses information from the
   outcome interval and cannot prove that carbohydrate effect is absent. Separate post-COB work
   found evidence consistent with carbohydrate action outlasting the recorded COB model.
4. **The action-balance estimator is model dependent.** It uses model-derived IOB and suggested
   delivery. Suggested SMB and temporary basal are not always confirmed enacted delivery.

These limitations affect the claimed TDD exponent, the proposed constants, and the idea that a
short observational window can set a person's physiological ISF. They do not invalidate the
equation replay itself.

---

## 5. Appropriate future use

The v-next equation may be retained as a research candidate, with these restrictions:

- no live dosing recommendation from the retrospective fits;
- no claim that `1/√TDD` is a physiological law;
- no use of the historical Tier 2 constant as measured ISF;
- evaluation first in simulation or shadow mode;
- prospective validation using confirmed delivery, explicit carbohydrate handling and a design
  that separates insulin exposure from controller response.

A prospective study should compare a static profile, v1, v2 and the candidate v-next form on
pre-specified safety and glucose outcomes. It should use non-overlapping evaluation periods and
report results by TDD, glucose, insulin formulation and controller type.

---

## Reproducibility

- Equations and tests: `inv008/dynisf.py`, `inv008/tests/`
- Historical candidate search: `fit_best_isf.py`, `results/best_isf_fit_results.{json,md}`
- Predictive glucose-shape test: `inv008/gradient_isf_fit.py`, `results/gradient_isf_fit.{json,md}`
- Same-window comparison: `inv008/head_to_head.py`, `results/head_to_head.{json,md}`
- Audit: `DYNAMIC-ISF-AUDIT-2026-09.md`
