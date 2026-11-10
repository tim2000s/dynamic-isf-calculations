---
title: "A profile-anchored dynamic sensitivity for open-source automated insulin delivery"
subtitle: "Initial proposal, with the TDD-based dynamic ISF in AndroidAPS and Trio as the comparator"
author: "Tim Street, Diabettech"
date: "10 November 2026"
---

## Status

This is a first draft for discussion. Nothing described here has run on a device or
in a simulator, and no parameter in it was fitted. The candidate is a modelling
choice assembled from the findings of documents 5, 6 and 7 in this series, and each
choice is marked where it is made.

The June 2026 proposal (document 3) was withdrawn for live dosing on 16 September
2026 and stays withdrawn. This document uses neither its square-root law nor its
learned level.

The earlier work did not leave one named model behind, so the reading taken of "the
updated model" is set out in the section headed "The candidate". It can be changed
there before anything is built on it.

Glucose is in mg/dL throughout because the code is. 100 mg/dL is 5.6 mmol/L and
200 mg/dL is 11.1 mmol/L.

## Summary

The existing equation, 1800 / (TDD × ln(glucose / divisor + 1)), does three jobs at
once. Each was examined separately in this series.

| Job | How the existing equation does it | What the analysis found | In the candidate |
|---|---|---|---|
| Sets the person's level | From total daily dose (TDD), a constant and the insulin type | Entered settings already follow TDD between people. The level implied changes by 22% with the insulin setting alone | The entered profile sensitivity |
| Moves the level with recent dosing | Blend of 8 hour, 24 hour and 7 day doses | Supported, with an exponent between about -0.65 and -1.0 within a person | Kept, as a bounded ratio to the person's own average |
| Makes a unit count for more as glucose rises | Logarithm of glucose | No support between 100 and 200 mg/dL away from meals. Open between 200 and 300 | Removed above 100 mg/dL, kept below it |

For one person with a profile sensitivity of 50 mg/dL per unit, 40 U a day and
rapid-acting insulin, at 200 mg/dL, stock AndroidAPS calculates 29.3, Trio's
logarithmic mode 41.7 and the candidate 50.0.

The equation is the small part of the change. The larger parts are a check of each
person's profile before switching (in AndroidAPS the profile sensitivity has not
been used for dosing while dynamic ISF was on) and a measurement of what the
removed glucose term was delivering after meals. Both come before a period in which
the two values are logged and only the existing one doses.

## The comparator

The comparator is the dynamic ISF in stock AndroidAPS (upstream master f77abd9b82,
31 October 2026) and stock Trio (upstream main 0b29cd534, 7 October 2026), read
from source on 10 November 2026. oref0 as published has autosens and no TDD-based
dynamic ISF, so it has nothing to replace. Loop and iAPS were not read.

### One person, worked through

Take the person above with glucose at 200 mg/dL and recent dosing equal to their
average.

AndroidAPS takes the 40 U, the glucose and a divisor of 55 for rapid-acting
insulin, and calculates 1800 / (40 × ln(200 / 55 + 1)) = 29.3 mg/dL per unit. The
profile's 50 does not enter. By that number alone, a correction to 100 mg/dL is
3.41 U.

Trio's logarithmic mode starts from the same equation and expresses it as a ratio
of the profile value, with a default adjustment factor of 0.8. The ratio comes out
above the upper autosens limit of 1.2 and is held there, which gives 50 / 1.2 =
41.7 and a correction of 2.40 U.

So the same equation returns 29.3 in one system and 41.7 in the other for this
person. The difference is Trio's bound around the profile and its default factor.

### Across glucose

Values in mg/dL per unit for the same person, from `inv013/candidate.py`
(`results/inv013_comparator.md`). AndroidAPS is at an adjustment factor of 100%,
Trio at its defaults.

| Glucose, mg/dL | AndroidAPS, divisor 55 | AndroidAPS, divisor 65 | AndroidAPS, divisor 75 | Trio logarithmic | Trio sigmoid | Candidate |
|---:|---:|---:|---:|---:|---:|---:|
| 70 | 54.8 | 61.6 | 68.3 | 68.5 | 55.7 | 63.1 |
| 80 | 50.1 | 56.1 | 62.0 | 62.6 | 53.7 | 57.7 |
| 100 | 43.4 | 48.3 | 53.1 | 54.3 | 50.0 | 50.0 |
| 120 | 38.9 | 43.0 | 47.1 | 48.6 | 47.1 | 50.0 |
| 150 | 34.2 | 37.6 | 41.0 | 42.8 | 44.3 | 50.0 |
| 200 | 29.3 | 32.0 | 34.6 | 41.7 | 42.4 | 50.0 |
| 250 | 26.3 | 28.5 | 30.7 | 41.7 | 41.8 | 50.0 |
| 300 | 24.1 | 26.1 | 28.0 | 41.7 | 41.7 | 50.0 |

Units to bring this person to 100 mg/dL, by the sensitivity alone (insulin on
board, predictions and delivery limits are left out):

| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Candidate |
|---:|---:|---:|---:|
| 150 | 1.46 | 1.17 | 1.00 |
| 200 | 3.41 | 2.40 | 2.00 |
| 250 | 5.71 | 3.60 | 3.00 |
| 300 | 8.29 | 4.80 | 4.00 |

This person's profile multiplied by their TDD is 2000. The AndroidAPS equation at
100 mg/dL implies 1737 with a divisor of 55, 1932 with 65 and 2124 with 75, so the
candidate's distance from the existing value at 100 mg/dL depends on the person and
on their insulin setting. In INV-010 the profile value ran from 0.56 to 1.80 of the
dynamic value across six people using Boost.

### What the code does

AndroidAPS (`OpenAPSSMBPlugin.calculateRawDynIsf`, `DetermineBasalSMB.kt`):

- The dose term is 0.33 of a weighted eight hour total scaled to a day, 0.34 of the
  seven day average and 0.33 of the most recent day, multiplied by the adjustment
  factor (default 100%, range 1 to 300%) and the profile percentage.
- The divisor is 55, 65 or 75 according to the insulin peak time.
- In dynamic mode the dose calculation uses this value in place of the profile
  sensitivity, so the entered value and its schedule through the day are unused.
- Each predicted glucose in the three prediction arrays gets its own sensitivity
  from the same equation. The dose is sized on a further value calculated at a
  glucose chosen by one of three branches.
- With any of the five dose totals missing, the plugin falls back to the profile
  and raises a notification.

Trio (`DynamicISF.swift`, a Swift port of the `trio-oref` JavaScript):

- Logarithmic mode multiplies the profile sensitivity out of the equation and puts
  it back, so the profile matters only through the autosens limits (defaults 0.7
  and 1.2), which bound the result between profile / 1.2 and profile / 0.7.
- Sigmoid mode is 1 at the profile target, so its level is the profile value. TDD
  enters only as the ratio of a weighted recent dose to the ten day average.
- The dose term is 0.65 of the last two hours' average and 0.35 of the ten day
  average by default. The insulin factor is 120 minus the peak time.

Trio's sigmoid mode therefore already has two of the candidate's properties, a
profile level and a bounded ratio. It differs from the candidate in what happens
above target.

## What the evidence supports

| Part | Finding | Source | Limit |
|---|---|---|---|
| Level between people | Entered sensitivity against TDD has slope -0.979 [-1.030, -0.934] in 581 people, median product 2067. Restricted to adults the product is 1915 | Document 6, constants update | The pooled record is 51% children. Entered settings record a decision and most are reached through the 1800 rule |
| Level for one person | Observational estimates run from 0.30 to 1.12 of the entered value according to construction | Document 7 | No construction was validated as a level |
| Level, by prediction | Overnight, a static 1800 rule missed by 34.1 mg/dL, v1 by 36.7 and the entered setting by 37.3 (1,626 people, trial cohorts). On 89 people using oref, the static profile missed by 20.3 and v1 by 24.6 | Document 6, INV-008 | The two cohorts disagree on whether the profile or a TDD rule is the better level. INV-008 windows were selected on the later outcome |
| TDD within a person | Exponent -0.645 (standard error 0.033, 1,313 people). The same estimator on simulated people returns -0.67 for a true -0.50 and -1.02 for a true -1.00 | Document 6 | Between about -0.65 and -1.0 is supported; the squared term of v2 was best for none of 1,626 people |
| Glucose, 100 to 200 mg/dL | Six shapes predicted within 0.05 mg/dL of each other against a bar of 0.5 set in advance, the flat one first. Apparent glucose dependence was +0.610 within four hours of a meal and +0.207 at six to nine hours | Document 6 | Observational, closed-loop records |
| Glucose, bias of prediction | On common windows the static profile's bias stayed between 5.0 and 7.8 mg/dL from 100 to 230 mg/dL. That of v1 ran from +3.7 to -33.2 | INV-008 | Rescaled loop predictions, outcome-selected windows |
| Glucose near target | Relative sensitivity measured 0.64 at 90 to 120 mg/dL against 1.00 at 120 to 150. The existing equation assumes 1.18 | Document 6 | A unit appears to do less near target, which is where the body's own defence against a fall acts |
| Glucose, 200 to 300 mg/dL | Fall per unit at 300 against 200 mg/dL was 0.66 [0.49, 0.73] in REPLACE-BG (24 people) and 0.73 [0.64, 0.93] under Loop (49 people). The existing equation's own ratio is 0.81 to 0.82 | Document 5 section 7, INV-012 | No simulated arm reached 300 mg/dL, so this estimator is unchecked there |

Two things follow for the design. The level evidence is mixed, and the profile is
chosen as the anchor on other grounds, given below. Above 200 mg/dL a decline of
about the size the existing equation assumes is not excluded.

## The candidate

sensitivity = profile sensitivity now ÷ clamp((recent dose ÷ own average dose)^b) × easing(glucose)

easing(glucose) = ln(P / D + 1) ÷ ln(min(glucose, P) / D + 1)

Worked for the same person at 80 mg/dL, with recent dosing 10% above their average:
the profile gives 50, the dose ratio of 1.10 is inside the limits and gives 45.5,
and the easing at 80 mg/dL is 1.15, so the value is 52.4. At any glucose from 100
mg/dL upward the easing is 1 and the value is 45.5.

Each part, and why it is this way:

- Where the level comes from: the profile sensitivity for the time of day. This is
  a modelling choice. The reasons are that the person and their clinic can see it
  and change it, that Trio's sigmoid mode already works this way, that the existing level moves by
  22% with the insulin setting for no physiological reason, and that document 7
  found no observational estimate fit to replace an entered value. The prediction
  evidence does not settle it, so a level of 1800 over the person's average TDD is
  carried as an evaluation arm.
- How much recent dosing moves it: the system's existing blended dose divided by
  the average it already holds (seven days in AndroidAPS, ten in Trio), raised to
  b. The default b = 1 is what the existing equation does, which leaves this half
  of the behaviour unchanged. b = 0.65 is an evaluation arm. The ratio is held
  inside the existing autosens limits. Across 960 people the blended dose ran from
  0.72 to 1.14 of the person's own average (tenth to ninetieth percentile), which
  moves the value by 1.39 and 0.88.
- What glucose does below the pivot P (default 100 mg/dL, a modelling choice): the
  existing logarithm, which raises the value by 15% at 80 mg/dL and 26% at 70. No
  correction is being given there, so this acts through the predictions: insulin
  on board is predicted to lower glucose further, which brings a reduced or zero
  basal rate forward. It runs against the measurement, which shows a unit doing
  less near target. It is kept because the measured effect is read as the body
  defending against a fall, and a prediction that relies on that defence errs
  towards more insulin. An arm without easing is carried.
- What glucose does above the pivot: nothing by default. The existing shape from
  200 mg/dL upward is an optional arm, since the 200 to 300 result is open.
- The divisor D stays in the easing only, where 55, 65 and 75 differ by under 0.03
  at 70 mg/dL. It could be fixed at one value.
- The adjustment factor has no role. Its job in the existing equation is to move
  the level, and the profile now does that.

With TDD history missing, the ratio is 1 and the value is the profile with easing.
The existing fallback and the candidate then agree above the pivot, so the first
week on a new pump does not switch mode.

## What has to be added

| Component | AndroidAPS today | Trio today | Needed for the candidate |
|---|---|---|---|
| Profile sensitivity in the dynamic path | Present in the profile object, unused in dynamic mode | Used | AndroidAPS: pass it into the calculation |
| Recent and average dose | Five totals, all required | Two hour and ten day averages | Nothing new collected. Expose the ratio |
| Equation | One function in the plugin | One function, two modes | One new function per system, with shared reference vectors from `inv013/candidate.py` |
| Sensitivity at predicted glucose | Three prediction lines | Forecast adjustment in logarithmic mode | Call the new function. Above the pivot the result no longer changes with predicted glucose |
| Sensitivity for the dose | Three branches | In the oref port | Call the new function at the same glucose each branch already chooses |
| Stored value | `variableSens` in the loop result | In the determination | Store the candidate beside it during shadow running |
| Settings | Adjustment factor, option to adjust sensitivity and target | Mode, two adjustment factors, weighting, limits | A mode selector. No new number is shown to the person in the first version |
| Profile check on switching | None needed today | None needed today | New. Compare the entered profile with what the existing equation has been calculating for that person, and show the ratio before the mode can be turned on |
| Tests | `ReplayApsResultsTest` with JavaScript parity assets | `DynamicISFTests`, `DynamicISFEnableTests` | Unit tests on the reference vectors, replay assets for the new mode, and a guard test that the existing mode's output is unchanged |
| After meals | Covered in part by the glucose term | Covered within the limits | See the next section |

Autosens needs a decision in AndroidAPS. INV-010 found under Boost that autosens
judged deviations against the profile value while the loop dosed on the dynamic
one, so part of its ratio reported the gap between the two. With the candidate the
two values coincide above the pivot when the dose ratio is 1. The stock autosens
path in dynamic mode was not re-read for this draft.

## What the replacement removes

By shape alone the existing equation gives 0.68 of its 100 mg/dL value at 200 mg/dL
and 0.56 at 300, so it asks for about 1.5 and 1.8 times the insulin for the same
distance from target. The candidate removes that.

The series' reading is that much of this extra insulin has been standing in for
carbohydrate still arriving after the absorption model stopped counting it. That
is the best-supported explanation and it is not proven. For people who announce
meals accurately the removal should cost little. For people who do not announce,
or whose meals run long, the glucose term may be doing useful work for a reason
other than the one it was designed on, and removing it without a replacement would
leave them with less insulin after meals.

The systems already have a route for uncounted carbohydrate: oref's deviation-based
unannounced-meal prediction and its dynamic carbohydrate absorption, which are
credited here as the right place for this work. Whether they are sufficient
without the glucose term is the first thing the evaluation has to answer. Two
outcomes are possible and both are acceptable: the candidate goes forward as
written, or it goes forward only with the above-200 option or a meal-side addition
for people who do not announce.

## Evaluation path

Arms: the existing equation as shipped; the profile alone with no dynamic term;
the candidate; the candidate without easing; the candidate with the option above
200 mg/dL; the candidate at b = 0.65; and the candidate with its level from 1800
over average TDD.

| Stage | What is done | What it decides |
|---|---|---|
| 1. Replay of the equation | Each arm's value calculated on recorded histories (the INV-008 pipeline), split by time since the last meal and by whether meals were announced | How much insulin the removed term was asking for, in whom and when. Descriptive, no pass rule |
| 2. Prediction | Each arm added to the same-window and overnight comparisons of INV-008 and INV-009 | The candidate is carried forward when its typical miss is within 0.5 mg/dL of the best arm, the bar document 6 used |
| 3. Simulation | Paired runs by seed in TimSim, with arms for unannounced meals and for a mis-set profile | A harm rule on time below 54 mg/dL, written before the runs, with thresholds taken from a null arm |
| 4. Shadow running | Both values calculated on the device each cycle, the existing one doses, both are stored | The distribution of the ratio per person, and the profile check above |
| 5. Use | The developer's own device first, then released lines | Stops on the harm rule |

Stages 1 and 2 are analysis and carry the analysis checklist. Stages 4 and 5
change code that can change a dose, so they carry unit and guard tests, the
repository's review checklist, a replay against real history and a fresh-context
review of the diff.

The way back at every stage is the mode selector: the existing equation stays in
the code unchanged and is the default until stage 5 has run. A build that includes
the new mode names the previous build.

## Open questions

- Whether the level should be the profile or a TDD rule. The two cohorts examined
  disagree, and a profile that has gone unused under dynamic ISF may be stale.
- Whether sensitivity falls between 200 and 300 mg/dL. Two cohorts suggest it does
  and the estimator has no simulated check in that range.
- Whether the easing below the pivot should exist, given that it runs against the
  measurement, and whether the pivot should be 100 or 120 mg/dL or follow the
  target.
- Whether TimSim's oref adapter runs the existing dynamic ISF closely enough to be
  the comparator in stage 3. This was not checked.
- How the bolus calculator in each system takes its sensitivity in dynamic mode.
  This was not read.
- Whether the maintainers of AndroidAPS and Trio want one shared definition. The
  two existing implementations of one equation differ by 42% for the worked person
  at 200 mg/dL, which argues for shared reference vectors whatever is decided here.

## Reproducing the numbers

`python3 -m inv013.candidate` writes `results/inv013_comparator.md`, which holds
every number in the comparator tables and the worked examples.
`inv013/tests/test_candidate.py` checks the stock equation against the replay code
of INV-008 and the candidate's shape. All other figures are quoted from the
documents named beside them.
