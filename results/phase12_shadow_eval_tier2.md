# Phase 12 — cohort shadow evaluation of Tier-2 v-next (sensitivity-anchored)

98 people with a usable measured-sensitivity fit (r² ≥ 0.1), 5,605,658 per-tick readings. Tier-2 anchors K_user = measured_ISF · √(median TDD); the glucose curve g(BG) is unchanged.

## Headline

- Measured sensitivity is a median **0.438× the profile ISF**, so at the level Tier-2 doses about **2.28× the insulin** of the person's current setting — much more aggressive than Tier-1, which preserves the level.
- **Without the clamp**, a correction dose changes by a median **214.3%** vs profile, and **97%** of readings dose more than 1.5× stronger than profile.
- **With the §8.2 level clamp** (floor the level at profile/1.5), the level clamp binds for a median **100%** of readings — i.e. for most people Tier-2 is pulled back to the same ceiling Tier-1 would allow; dose change vs profile falls to a median **96.1%**.
- vs Tier-1: median ISF ratio **0.438** unclamped, **0.667** clamped.

## Sensitivity ratio vs profile, by glucose band

| BG band | Tier-2 no clamp | Tier-2 clamped |
|---|---|---|
| 54-80 | 0.653 | 0.943 |
| 80-100 | 0.509 | 0.744 |
| 100-120 | 0.402 | 0.591 |
| 120-150 | 0.294 | 0.453 |
| 150-200 | 0.198 | 0.317 |
| 200-260 | 0.143 | 0.235 |

*<1 = more aggressive (more insulin than the current profile). Unclamped Tier-2 is well below 1 across the range; the clamp lifts the level back toward profile, leaving the glucose-shape behaviour intact.*

## By per-user median TDD

| TDD band | users | median measured/profile | median frac clamp binds |
|---|---|---|---|
| <15 | 1 | 0.492 | 0.999 |
| 15-25 | 11 | 0.433 | 1.0 |
| 25-40 | 39 | 0.497 | 0.993 |
| 40-65 | 36 | 0.411 | 1.0 |
| 65+ | 11 | 0.377 | 1.0 |

The chart draws one mark per person and is kept in the private data folder (REVIEW.md item 22).

**Reading.** Tier-2 is materially more aggressive than the person's current profile — about 2.28× the correction insulin at the level, before glucose scaling. The §8.2 level clamp is doing real work here: it binds for the majority of readings and converts Tier-2 into "no more than 1.5× stronger than profile", which is also the Tier-1 ceiling. So clamped Tier-2 and Tier-1 differ mainly where measured sensitivity is *weaker* than profile (a minority of users).

This is the dosing-magnitude side of the Tier-2 question only. The data-derived study (Phases 5–6) showed the measured-sensitivity anchor is hypo-biased — it reads most sensitive for the people who already run low — so even the clamped form needs forward, outcome-based validation before it doses. It is not a deployable default.

*Caveat: counterfactual decision-level replay, not closed-loop outcomes; measured ISF carries carb/endogenous-glucose confounds and per-user CI; single cohort; median-TDD anchor stands in for the weekly recalibration.*