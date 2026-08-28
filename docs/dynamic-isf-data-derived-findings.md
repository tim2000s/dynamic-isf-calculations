# Dynamic ISF — can sensitivity be derived from data? (feasibility findings)

**2026-06-07** · Tim Street / Claude · companion to the methodology and v1/v2 analysis

---

## Question

Can a user's insulin sensitivity (ISF) be derived from device-observable data well enough to
remove the need to enter it — and, if so, how? This summarises a six-phase feasibility study
on 138 open-source AID users (Trio, AAPS, OpenAPS), reusing the empirical-ISF ΔIOB estimator
and the √TDD law established in the equation-fit analysis.

## The two separable parts of ISF

ISF has a shape, meaning how it scales with TDD and glucose, and a per-user level. The two
behave very differently.

The shape is derivable and settled. Sensitivity scales with total daily dose as roughly
1/√TDD, and that holds across users, with the glucose curve a separate dimension on top. The
shape needs no input from anyone.

The level is the hard part, and the subject of this study. A single global level leaves the
typical person about 57% off, so the question is whether the per-user level can be recovered
from data rather than entered.

## The feasibility arc (Phases 1–6)

| phase | question | finding |
|---|---|---|
| 1 | Is a short-window empirical ISF dosing-grade? | No. ±10% CI only at 60–90 d; a 14-day estimate swings **±34%** block-to-block, **memoryless** (lag-1 ACF ≈ 0). |
| 2 | Is the variation trackable structure? | Mostly not. Circadian ~2–4%, weekly ~1%, slow ~2%; **~90% residual, near-white**, partly measurement noise. |
| 3 | Can a stable baseline be recovered by shrinkage? | Yes, in self-consistency terms: cold-start ±42%; shrinking own-data toward the √TDD prior → **±22% (14 d) → ±13% (90 d)**. |
| 5 | Does the derived (sensitivity-anchored) level beat tuned ISF on real outcomes? | **No — and it's unsafe.** Users dosing far "weaker" than their measured sensitivity have **more** hypoglycaemia, not less (logR-vs-TBR **+0.38**). The empirical ISF is hypo-biased. |
| 6 | Does overnight data de-bias the sensitivity estimate? | **No.** Overnight is noisier (±50% test-retest) and carries the **same** bias (+0.39). The bias is intrinsic to observational sensitivity, not a daytime artefact. |
| 7 | Can outcome-based tuning give a safe per-user level? | **Yes, in principle.** Working-anchored cold start (355/√TDD) lands **±28% of the users' tuned ISF** with no entry; a bounded outcome nudge removes the residual hypo-signature (logR-vs-TBR +0.21 → **−0.04**) where the sensitivity route made it worse (+0.38). |

## What this settles

1. The √TDD shape is derivable and needs no entry. Anchored to working (tuned) profiles the
   constant is K ≈ 355 (ISF ≈ 355/√TDD); anchored to measured sensitivity it is ≈ 145. Since
   Phases 5 and 6 show the measured anchor is unsafe, the working-anchored 355 is the right
   cold start.

2. A zero-entry cold start is feasible and reasonably good. From observable TDD alone,
   ISF ≈ 355/√TDD places a brand-new user within ±28% of the value they would have tuned to
   (Phase 7). Against the measured empirical baseline the figure is ±42%, but that baseline is
   the hypo-biased one, so the working-ISF comparison is the number that matters for deployment.

3. The per-user level cannot be safely derived from observed sensitivity, in any window. The
   ΔIOB sensitivity estimate is entangled with the user's existing insulin excess or deficit:
   it reads as "very sensitive" (low ISF) precisely for people who run low, so dosing to it
   would give the most insulin to the most hypo-prone. Restricting to clean overnight data does
   not fix this. It is intrinsic to measuring sensitivity from observational glucose against
   insulin.

4. Beyond a stable baseline there is little left to track for correction ISF (Phase 2). The
   right architecture is therefore not an elaborate adaptive ISF but a stable baseline plus the
   √TDD demand-track, with autosens covering the bounded sub-daily and safety role it already
   serves.

## The constructive conclusion

A data-derived, zero-entry ISF is feasible **only in this form**:

> **√TDD shape (universal) + working-anchored cold start (≈ 355/√TDD, ±42%, no entry) +
> bounded outcome-based self-tuning of the level** — nudging ISF weaker when unexplained lows
> appear and stronger when sustained unexplained highs appear, damped and clamped.

It is **not** feasible as "measure the user's sensitivity and dose to it" — that route is
hypo-biased regardless of how clean the window is. The level must be tuned from **outcomes**
(the autotune/autosens philosophy), not from a sensitivity regression. That is the autotune
design, now empirically explained: it is built on the human-tuned anchor and outcome-attribution
precisely because the measurement route fails.

## Caveats

Decision-level and observational (no closed-loop simulation); single cohort, mostly users who
never ran dynamic ISF; the empirical estimator carries carb/endogenous-glucose confounds; "level"
errors in Phase 3 are self-consistency, not outcome-validated correctness. A closed-loop trial is
required before any of this drives dosing. Nothing here is medical advice.

## Reproducibility

`inv008/phase{1,2,3,5,6}_*.py` and `results/phase*_*.{json,md}` in
`github.com/tim2000s/dynamic-isf-calculations`. Inputs device-observable; the user's profile ISF
is used as a benchmark only.


## Addendum, 27 August 2026: replication on independent data

Phases 5 and 6 above concluded that observational sensitivity is biased low, that
dosing to it is unsafe, and that the bias is intrinsic rather than a daytime
artefact. Those conclusions rested on 138 users of open-source systems with
insulin on board taken from the loop's own record.

They have since been tested on 1,679 people across six Jaeb trial cohorts, ages 2
to 82, with insulin on board reconstructed from delivery records rather than
reported. The result does not carry over unchanged. Running the construction used
here, the fall over the horizon divided by insulin on board at the start, gives a
ratio of measured to entered sensitivity of 0.72 to 1.12 on those cohorts, and the
open-loop cohort exceeds its entered settings at six hours. Constructions that
instead divide by the action of all insulin above the programmed basal give 0.30
to 0.66 on the same records.

The 0.41 reported here does not survive re-examination on its own records either.
Its two anchors were fitted on different subsets, 114 people for the measured
constant and 138 for the entered one, and both were regression fits rather than
medians. Matched person by person across the same 112 records and taken as medians,
the constants are 207 and 333, a ratio of 0.62, and the median of per-person ratios
is 0.60. As a Walsh-style constant that is 1432 measured against 2219 entered.

So the ratio of 0.41 reported here is not reproduced at that value on trial data,
and the spread between constructions on identical records is wider than the spread
between cohorts. Why the oref archive sits at 0.41 where trial data sits near 1.0
under the same construction is open, and is a question about these records rather
than about the method. The conclusions of Phases 5 and 6, that the measured level
is unsafe to dose to and that overnight does not de-bias it, are untouched by this:
they rest on outcomes rather than on the ratio. See Dynamic-ISF-7.
