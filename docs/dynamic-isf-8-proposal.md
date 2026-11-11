---
title: "Adding a glucose-indexed sensitivity curve to open-source automated insulin delivery"
subtitle: "Initial proposal, with the TDD-based dynamic ISF in AndroidAPS and Trio as the comparator"
author: "Tim Street, Diabettech"
date: "11 November 2026"
---

## Status

This is a draft for discussion. Nothing described here has run on a device or in a
simulator, and nothing was fitted.

The model proposed is the curve of sensitivity against glucose published by
Romero-Ugalde and colleagues in 2025 [1], linked to total daily dose (TDD). The
shape is theirs. The link to TDD is derived in this document and appears nowhere in
the publication.

An earlier draft is in git (commit fd3662e). It proposed a sensitivity
anchored to the entered profile and flat above 100 mg/dL, and is kept in the code
as one arm of the evaluation.

Glucose is in mg/dL throughout because the code is. 100 mg/dL is 5.6 mmol/L and
200 mg/dL is 11.1 mmol/L.

## Summary

The curve model and the existing dynamic ISF have the same form. Each is a function
of glucose divided by a daily dose:

| | Sensitivity | Value of the glucose function at 100, 150, 200, 300 mg/dL |
|---|---|---|
| Existing | 1800 / ln(glucose / divisor + 1), over TDD | 1737, 1368, 1173, 965 (divisor 55) |
| Curve model | 1915 × Q(glucose) / Q(180), over TDD | 4378, 2492, 1674, 1148 |

Q is the published curve in mg/dL per unit. 1915 is the product of entered
sensitivity and daily dose measured in adults in this series (document 6), taken
to hold at 180 mg/dL, the reference glucose the publication uses. That choice sets
the level and is a modelling choice.

Replacing one equation with the other is a small change to the code: the glucose
function is swapped and the insulin-type divisor goes. The change in dosing is
large. At the same TDD the curve asks for less insulin than stock AndroidAPS at
every glucose up to 362 mg/dL with a divisor of 55, by a factor of 2.52 at 100
mg/dL, 1.43 at 200 and 1.19 at 300. It falls by 2.61 from 100 to 200 mg/dL where
the existing equation falls by 1.48.

A second way of scaling the curve is carried as an alternative: the entered
profile value in place of 1915 / TDD. It suits Trio, which already expresses
dynamic ISF as a ratio of the profile, and it does not move with recent dosing.

Most of what has to be added sits around the equation. The level is the least
settled part, because the publication gives no daily dose: three ways of setting
it give reference doses of 24, 54 and 71 to 78 U. A decision is needed on where
along the curve a correction is read. Trio's default limits would remove most of
the curve.

The evidence is set out in its own section. The curve is reproducible. The
analysis in this series does not support reading its decline between 150 and 250
mg/dL as a change in sensitivity, and whether dosing to it does better than the
existing equation has not been tested. The path proposed is therefore replay and
simulation first, then logging beside the existing equation, with dosing only
after those.

## The model

### The published curve

The publication [1] used records from 11,201 adults with type 1 diabetes on the DBLG1
closed-loop system. One event, as it defines them: glucose peaks at 240
mg/dL at least an hour after a meal, with under 2.5 g of carbohydrate still on
board, and falls to a minimum of 160 mg/dL more than 20 minutes later. Between
the two, insulin on board dropped and some insulin was delivered above the
reference basal, 2.0 U in all. The sensitivity for that event is 80 / 2.0 = 40
mg/dL per unit. The window ends where glucose stopped falling, so events are
selected on their outcome and the figure is a ratio of outcomes.

The reported result is a sensitivity that rises from 65 to 79 mg/dL per unit
between 50 and 100 mg/dL, falls to 35 at 180 and is below 20 above 300. Against a
dose sized at 180 mg/dL, the authors put the dose at about half for glucose up to
120 mg/dL and about double above 340.

The poster's figure fits the curve in two pieces, and the coefficients used here
were checked against that figure on 11 November 2026. For glucose above 100 and
up to 380 mg/dL it is a quartic (`inv008.dynisf.quartic_isf`), which gives 80.6 at
100 mg/dL, 35.3 at 180 and 21.1 at 300, and dose multiples of 0.56 at 120 and 1.95
at 340, in line with the abstract. For glucose up to 100 mg/dL it is a quadratic,
which gives 79.0 at 100 and 65.3 at 60. The two pieces differ by 2% where they
meet.

### Linking it to daily dose

The published values are absolute. A person on 20 U a day and a person on 80 U a
day would get the same 35 mg/dL per unit at 180 mg/dL, so the curve needs a
scale for each person before it can dose. Two links are set out.

The TDD link, which is the candidate:

sensitivity = 1915 × Q(glucose) ÷ Q(180) ÷ TDD

The profile link, carried as an alternative:

sensitivity = profile sensitivity × Q(glucose) ÷ Q(180)

For a person on 40 U a day with a profile of 50 mg/dL per unit, the TDD link
gives 109 mg/dL per unit at 100 mg/dL, 62 at 150, 48 at 180, 42 at 200 and 29 at
300. The profile link gives 114, 65, 50, 44 and 30. The two agree for anyone whose
profile times daily dose is 1915, and differ by that ratio for everyone else.

| | TDD link | Profile link |
|---|---|---|
| Fits | AndroidAPS, where the profile value is unused in dynamic mode | Trio, where dynamic ISF is a ratio of the profile |
| Entered sensitivity and its schedule through the day | Unused | Kept |
| Recent dosing | Moves the value, as now | No effect |
| What has to be chosen | The level (1915) and the definition of TDD | The glucose at which the profile is taken to hold |

The reference implementation (`inv013/candidate.py`) anchors both links at 180
mg/dL, a modelling choice made because that is the publication's own reference
point. Its treatment of the ends is a modelling choice as well. Below 100 mg/dL
the curve is held at its value at 100, which gives less insulin there than the
published quadratic, which falls towards 65. Work in this series in March 2026
declined the quadratic for dosing on the same ground. Above 380 mg/dL the curve is
held at its value at 380, where the published fit ends.

### The level

The publication reports no daily dose, so the level cannot be taken from it. The
reference dose is the daily dose at which the published values would apply
unchanged, and it follows from whichever level is assumed:

| Level assumed | Reference dose, U | Curve × TDD at 200 mg/dL |
|---|---:|---:|
| Adult entered product, 1915, at 180 mg/dL (this draft) | 54.3 | 1,674 |
| The 1800 rule at 105 mg/dL (this series, March 2026) | 23.7 | 731 |
| The same construction on Loop, Trio and oref0 records, at 200 mg/dL | 70.8 to 77.8 | 2,182 to 2,398 |

These differ by a factor of three. The first asks that a person's entered
sensitivity be right at 180 mg/dL. The second asks that it be right near target
and puts the whole curve at less than half the level. The third is the level
measured by the method on open-source records, 30 to 43% above this draft's, and
it comes from event rules that differ from the published ones. The existing
adjustment factor moves the level for one person in either direction. The default
is the open point.

## The comparator

The comparator is the dynamic ISF in stock AndroidAPS (upstream master f77abd9b82,
31 October 2026) and stock Trio (upstream main 0b29cd534, 7 October 2026), read
from source on 10 November 2026. oref0 as published has autosens and no TDD-based
dynamic ISF, so it has nothing to replace. Loop and iAPS were not read.

AndroidAPS (`OpenAPSSMBPlugin.calculateRawDynIsf`, `DetermineBasalSMB.kt`):

- The dose term is 0.33 of a weighted eight hour total scaled to a day, 0.34 of the
  seven day average and 0.33 of the most recent day, multiplied by the adjustment
  factor (default 100%, range 1 to 300%) and the profile percentage.
- The divisor is 55, 65 or 75 according to the insulin peak time, which moves the
  level at 100 mg/dL by 22% between the first and last.
- In dynamic mode the dose calculation uses this value in place of the profile
  sensitivity, so the entered value and its schedule through the day are unused.
- Each predicted glucose in the three prediction arrays gets its own sensitivity
  from the same equation. The dose is sized on a further value calculated at a
  glucose chosen by one of three branches.
- With any of the five dose totals missing, the plugin falls back to the profile
  and raises a notification.

Trio (`DynamicISF.swift`, a Swift port of the `trio-oref` JavaScript):

- Logarithmic mode expresses the same equation as a ratio of the profile
  sensitivity, with a default adjustment factor of 0.8, and holds the ratio inside
  the autosens limits (defaults 0.7 and 1.2). The result is bounded between
  profile / 1.2 and profile / 0.7.
- Sigmoid mode is 1 at the profile target, so its level is the profile value, and
  it is bounded by the same limits.
- The dose term is 0.65 of the last two hours' average and 0.35 of the ten day
  average by default.

## Side by side

One person: profile sensitivity 50 mg/dL per unit, 40 U a day, rapid-acting
insulin, recent dosing equal to their average. Values in mg/dL per unit, from
`results/inv013_comparator.md`.

| Glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Trio sigmoid | Curve, TDD link | Curve, profile link | Profile |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | 50.1 | 62.6 | 53.7 | 109.4 | 114.3 | 50.0 |
| 100 | 43.4 | 54.3 | 50.0 | 109.4 | 114.3 | 50.0 |
| 120 | 38.9 | 48.6 | 47.1 | 86.2 | 90.1 | 50.0 |
| 150 | 34.2 | 42.8 | 44.3 | 62.3 | 65.1 | 50.0 |
| 180 | 31.0 | 41.7 | 42.9 | 47.9 | 50.0 | 50.0 |
| 200 | 29.3 | 41.7 | 42.4 | 41.9 | 43.7 | 50.0 |
| 250 | 26.3 | 41.7 | 41.8 | 33.7 | 35.1 | 50.0 |
| 300 | 24.1 | 41.7 | 41.7 | 28.7 | 30.0 | 50.0 |
| 380 | 21.8 | 41.7 | 41.7 | 20.7 | 21.7 | 50.0 |

Against each comparator the TDD link moves dosing differently:

- Against stock AndroidAPS it gives less insulin up to 362 mg/dL with a divisor of
  55, 339 with 65 and 313 with 75, by most near target. Above those it gives more.
- Against Trio at its defaults it gives less insulin below about 200 mg/dL and more
  above, because Trio's upper limit holds its value at 41.7 by 180 mg/dL.
- Against this person's profile it gives less insulin below 174 mg/dL and more
  above.

## Reading the curve for a dose

The poster's figure plots sensitivity against a glucose labelled G(t), in bands of
20 mg/dL from 60 to 380. Neither the abstract nor the figure says which glucose of
an event that is, and the rest of the poster could not be retrieved. This draft
takes it to be the glucose halfway between the peak and the minimum, which is an
assumption. The existing code reads its equation at the current glucose, or at a
glucose chosen by one of the three branches. For a flat sensitivity the two are
the same. For a curve that falls this steeply they differ.

Units to bring the same person to 100 mg/dL, by the sensitivity alone (insulin on
board, predictions and delivery limits are left out):

| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Profile | Curve, read at start | Curve, read at midpoint |
|---:|---:|---:|---:|---:|---:|
| 150 | 1.46 | 1.17 | 1.00 | 0.80 | 0.61 |
| 200 | 3.41 | 2.40 | 2.00 | 2.39 | 1.61 |
| 250 | 5.71 | 3.60 | 3.00 | 4.46 | 3.01 |
| 300 | 8.29 | 4.80 | 4.00 | 6.97 | 4.78 |

Reading the curve at the starting glucose asks for 31 to 49% more insulin than
reading it at the midpoint. From 200 mg/dL the difference is 2.39 against 1.61 U,
on either side of the profile's 2.00. Integrating along the path gives 1.62 U,
close to the midpoint reading.

This is a decision the curve model forces and the existing equation did not. The
proposal's default is the midpoint between the glucose being corrected and the
target, a modelling choice. The reading at the starting glucose is carried as an
arm.

## What has to be added

| Component | Existing | Curve model | To add |
|---|---|---|---|
| Glucose function | 1800 / ln(glucose / divisor + 1) | The published quartic over its value at 180 mg/dL | One function per system with shared test vectors |
| Insulin-type divisor | 55, 65 or 75 | None | Removed in this mode |
| Level | 1800 | 1915 at 180 mg/dL | A default, chosen after stage 2 below. The three candidates span a factor of three |
| Daily dose | AndroidAPS blend of five totals. Trio weighted two hour and ten day averages | Not defined by the publication | Default: the existing blend, so that this half of the behaviour is unchanged. Arm: the plain average |
| Level for one person | Adjustment factor | The same | The factor stays. Its guidance has to be worked out again, since at 100% the value at 100 mg/dL is 2.5 times the existing one |
| Profile link | AndroidAPS none. Trio a ratio of the profile | Profile × Q(glucose) / Q(180) | A second mode, or the Trio form of the first. A decision for the maintainers |
| Where the curve is read | Current glucose and three branches | Not stated in the publication | A reading rule, tested in replay |
| Below 100 and above 380 mg/dL | The equation continues | Published quadratic falls below 100. No fit above 380 | Held at the end values. Arms below 100 mg/dL: the published quadratic, and the existing logarithm's rise |
| Limits | AndroidAPS none. Trio 0.7 to 1.2 around the profile | Runs from 0.41 to 2.19 of the worked person's profile | A decision for Trio. Inside its default limits the curve runs only from 137 to 201 mg/dL for the worked person |
| Predictions | Sensitivity at each predicted glucose | Rises 2.61 times as predicted glucose falls from 200 to 100 mg/dL | The same call sites. Replay of the effect on minimum predicted glucose and on time at zero basal |
| Stored value | `variableSens` in the loop result | | Both values stored while the existing one doses |
| Settings | Adjustment factor and mode options | | A mode selector. No new number shown to the person in the first version |
| Tests | `ReplayApsResultsTest` with JavaScript parity assets. `DynamicISFTests` | | Unit tests on the reference vectors, replay assets for the new mode, a guard test that the existing mode's output is unchanged |

Code change points. AndroidAPS: `calculateRawDynIsf`, the three per-prediction
lines and the three dose-sensitivity branches in `DetermineBasalSMB.kt`, and the
preference screen. Trio: `DynamicISF.swift`, the forecast adjustment used in
logarithmic mode, and the `trio-oref` JavaScript it was ported from.

## What is known about the curve

The reproduction in this series (INV-012) was built before the abstract was read
and its event rules differ: three hours clear of a meal where the publication uses
one hour with under 2.5 g on board, and a peak that stands 15 mg/dL clear where
the publication uses a 20 minute interval. It also multiplied each event by the
person's TDD. The rows below compare ratios, which do not depend on that.

| Question | Finding | Source |
|---|---|---|
| Is the curve reproducible? | A similar construction on 1,790 people in eight cohorts gave a ratio of 2.2 to 4.0 between 100 and 200 mg/dL, against 2.61 in the published curve. With no algorithm dosing (REPLACE-BG) it was 2.55 [2.25, 2.82] | Document 5 section 7, INV-012 |
| Does the level carry across systems? | At a midpoint of 200 mg/dL the construction gave 2,182 under Loop, 2,273 under Trio, 2,398 under oref0 and 3,283 in REPLACE-BG, against 1,674 from this draft's link. TDD definitions were not matched between sources | `results/inv012_peak_trough.md` |
| What does the method return for a constant sensitivity? | Simulated people with a true value of 2,000 at every glucose returned a falling curve reading 1.4 to 3.4 times the true value. Their ratio from 200 to 250 mg/dL was 0.87 [0.83, 0.90] and 0.92 [0.89, 0.95], against 0.84 to 1.05 in the six trial cohorts and 0.80 in the published curve | Document 5 section 7 |
| What produces the decline between 100 and 200 mg/dL? | The median fall is 30 to 43 mg/dL at 100 and 39 to 64 at 200. The median insulin used is 0.02 to 0.2 U at 100 and 0.3 to 1.5 U at 200. At 100 mg/dL, 9 to 58% of events had no net insulin used. A fitted fall of 27 to 44 mg/dL at zero net insulin reflects basal and earlier doses acting, and the choice of events that end at a minimum | Document 5 section 7, `results/inv012_peak_trough.md` |
| Above 200 mg/dL? | Fall per extra unit at 300 against 200 mg/dL was 0.66 [0.49, 0.73] in REPLACE-BG (24 people) and 0.73 [0.64, 0.93] under Loop (49 people), against 0.69 in the published curve. No simulated arm reached 300 mg/dL | Document 5 section 7 |
| Has the shape been tried as a predictor on oref records? | On 80 people and 62,344 overnight windows, with each shape anchored to the profile at 100 mg/dL, the typical miss was 24.2 mg/dL for the static profile and 29.8 for this shape. At 205 to 260 mg/dL it was 37.8 for the static profile and 73.7 for this shape | `oref-investigations/results/bridge_diabeloop.md` |

Four limits apply to that table. The simulation did not reach 100 or 300 mg/dL,
so the steep part of the curve near target has not been reproduced from a constant
sensitivity and the case there rests on the split between fall and insulin used.
The overnight windows in the last row were selected on the later outcome, as the
audit of 16 September 2026 recorded, and that test anchored the shape at 100 mg/dL
where this draft anchors at 180. Every row is observational. And the publication
is a conference abstract from the system's manufacturer, with no full paper found.

What follows for the proposal is narrower than a verdict on the curve. As a
description of net fall per net unit in events chosen from peak to minimum, it is
reproducible on other people's data. As a statement that a unit does less at
higher glucose, it is not supported between 150 and 250 mg/dL and is open above
that. Whether a loop that doses to it does better than one dosing to the existing
equation or to a flat profile is a separate question, and no evidence in this
series bears on it directly.

A June 2026 proposal that used this curve's shape (document 3) was withdrawn for
live dosing on 16 September 2026. The audit left such candidates open to
simulation, shadow replay and prospective study. The path below stays inside that.

## Evaluation path

Arms: the existing equation as shipped; the profile alone; the TDD link read at
the midpoint; the TDD link read at the starting glucose; the TDD link on plain
average TDD; the profile link; and the profile-anchored form of this document's
first draft, which is flat above 100 mg/dL and is the arm the curve has to beat
for its shape to matter.

| Stage | What is done | What it decides |
|---|---|---|
| 1. Replay of the equation | Each arm's value calculated on recorded histories (the INV-008 pipeline), split by glucose, by time since the last meal and by whether meals were announced | How much insulin each arm asks for, in whom and when. Descriptive, no pass rule |
| 2. Curve on each system | `inv012/peak_trough.py` run on AndroidAPS and Trio records under the published event rules, with the simulated constant-sensitivity arms extended to 100 and 300 mg/dL | The default level, and whether the decline near target and above 250 mg/dL appears without a real effect |
| 3. Simulation | Paired runs by seed in TimSim, with arms for unannounced meals and for a mis-set adjustment factor | A two-sided harm rule written before the runs: time below 54 mg/dL, and time above 250 mg/dL, since against stock AndroidAPS this model gives less insulin over most of the range |
| 4. Shadow running | Both values calculated on the device each cycle, the existing one doses, both are stored | The distribution of the ratio per person and the starting factor |
| 5. Use | The developer's own device first, then released lines | Follows only after stages 3 and 4 pass. Stops on the harm rule |

Stages 1 and 2 are analysis and carry the analysis checklist. Stages 4 and 5
change code that can change a dose, so they carry unit and guard tests, the
repository's review checklist, a replay against real history and a fresh-context
review of the diff.

The way back at every stage is the mode selector. The existing equation stays in
the code unchanged and remains the default until stage 5 has run, and a build
that includes the new mode names the previous build.

## Definitions to fix

- The level, and with it the reference dose.
- Which daily dose the curve is divided by, and over how many days.
- How the reference basal is defined, and which insulin action model gives the
  insulin on board. The abstract states neither and the level depends on both.
- Which glucose an event is filed under. The poster's axis is labelled G(t) without
  a definition, and this draft assumes the midpoint.
- Whether the published quadratic is used below 100 mg/dL. This draft holds the
  value at 100 instead.

## Open questions

- Whether the published curve may be built into open-source code. It was derived
  by a manufacturer from its own users' records, and its authors have not been
  asked.
- Whether a level set from entered settings at 180 mg/dL applies to any
  open-source system. The construction on those systems' records reads 30 to 43%
  higher at 200 mg/dL.
- Whether the steep part between 100 and 150 mg/dL should be used for dosing at
  all, given that a large share of events there had no net insulin used.
- How the bolus calculator in each system takes its sensitivity in dynamic mode,
  and what stock autosens does there. Neither was read.
- Whether TimSim's oref adapter runs the existing dynamic ISF closely enough to be
  the comparator in stage 3. This was not checked.
- Whether the maintainers of AndroidAPS and Trio want one shared definition. Their
  two implementations of the existing equation already differ by 42% for the worked
  person at 200 mg/dL.

## Reproducing the numbers

`python3 -m inv013.candidate` writes `results/inv013_comparator.md`, which holds
every number in the comparison tables and the reference doses.
`inv013/tests/test_candidate.py` checks the stock equation against the replay code
of INV-008, the curve against the values stated in the abstract, and the ways of
reading it. All other figures are quoted from the documents named beside them. The
worked event in "The published curve" is an illustration with round numbers and is
no one's record.

## Reference

1. Romero-Ugalde HM, Adenis A, Cartier-Michaud T, Gauthier P, Desir C, Le
   Roux-Mallouf T, Huneker E, Benhamou PY. 920-P: The Relationship between
   Glycemia and Insulin Sensitivity Factor. Diabetes 2025;74(Supplement 1):920-P.
   doi:10.2337/db25-920-P
