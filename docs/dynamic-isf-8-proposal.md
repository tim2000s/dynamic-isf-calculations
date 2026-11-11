---
title: "Adding a glucose-indexed sensitivity curve to open-source automated insulin delivery"
subtitle: "Initial proposal, with the TDD-based dynamic ISF in AndroidAPS and Trio as the comparator"
author: "Tim Street, Diabettech"
date: "11 November 2026"
---

## Summary

Take an adult who uses 40 U of insulin a day and has entered an insulin sensitivity
factor (ISF) of 50 mg/dL per unit. At a glucose of 200 mg/dL, the dynamic ISF in
stock AndroidAPS works to 29.3 mg/dL per unit for that person, and the curve
proposed here works to 41.9. At 100 mg/dL the two values are 43.4 and 109.4. A
correction from 200 down to 100 mg/dL comes to 3.41 U on the AndroidAPS value and
1.61 U on the proposal's. Those two doses are the arithmetic of the sensitivity
alone, before insulin on board and the loop's predictions change them.

The curve is the relationship between glucose and sensitivity published by
Romero-Ugalde and colleagues in 2025 [1], measured on the records of 11,201 adults
using a commercial closed-loop system. Its shape is theirs. The published values
take no account of how much insulin a person uses, so this document scales the
curve by total daily dose (TDD). That scaling is a modelling choice made here and
appears nowhere in the publication.

The existing equation and the proposed one have the same form, a function of
glucose divided by a daily dose. Replacing one with the other is a small change to
the code: the glucose function is swapped and the divisor that depends on insulin
type is dropped. The change in dosing is large. At the same daily dose the curve
gives a higher sensitivity than stock AndroidAPS, and so less insulin, at every
glucose up to 362 mg/dL with the divisor at 55. The ratio between the two is 2.52
at 100 mg/dL, 1.82 at 150, 1.43 at 200 and 1.19 at 300.

Two parts of the model are unsettled. The first is its level. The publication
reports no daily dose, and the choice made in this draft amounts to saying that
the published values fit a person on 54 U a day. An earlier choice in this series
put that figure at 24 U, and levels measured on other records in this series put
it between 71 and 107 U. The second is where on the curve a correction is read.
Reading it at the glucose being corrected asks for 31 to 49% more insulin than
reading it halfway to target, and the publication does not say which of those its
figure corresponds to.

The evidence in this series supports part of the curve. A similar construction on
1,790 people in eight cohorts reproduced the steep decline between 100 and 200
mg/dL. The same construction on simulated people whose sensitivity was constant
also returned a falling curve, so the decline between 150 and 250 mg/dL is not
evidence that a unit of insulin does less at higher glucose. Whether a loop dosing
to the curve does better than one dosing to the existing equation has not been
tested.

The evaluation proposed starts with replay and simulation. It then logs the new
value on the device while the existing equation doses, and it reaches dosing only
after both.

This is a draft for discussion. The proposed equation has not run on a device or
in a closed-loop simulator, and no parameter was fitted to data. Glucose is in
mg/dL throughout because the code is. 100 mg/dL is 5.6 mmol/L and 200 mg/dL is
11.1 mmol/L.

## The existing equation

Dynamic ISF replaces the single entered sensitivity with a value that falls as
glucose rises and as recent insulin use rises. The formulation is due to Chris
Wilson, and it is what people using AndroidAPS and Trio run today:

sensitivity = 1800 ÷ (TDD × ln(glucose ÷ divisor + 1))

For the person above, with a divisor of 55, that is 1800 ÷ (40 × ln(200 ÷ 55 + 1)),
or 29.3 mg/dL per unit at 200 mg/dL.

AndroidAPS was read at upstream master f77abd9b82 of 31 October 2026
(`OpenAPSSMBPlugin.calculateRawDynIsf` and `DetermineBasalSMB.kt`). The divisor
depends on the insulin's peak time and runs from 55 to 75, which moves the value at
100 mg/dL by 22%. The daily dose is a blend. A third of it is the last eight hours
of insulin, weighted and scaled up to a day, and the remaining two thirds are split
almost evenly between the seven day average and the most recent day. The blend is
then multiplied by the adjustment factor, which defaults to 100% and can be set
from 1 to 300%, and by the profile percentage. In dynamic mode this value takes
the place of the profile sensitivity in the dose calculation, so the entered value
and its schedule through the day go unused. Each predicted glucose in the loop's
forecasts gets its own sensitivity from the same equation, and the dose itself is
sized on a further value read at a glucose chosen by one of three branches in the
code. The plugin falls back to the profile and raises a notification whenever one
of the dose totals is missing.

Trio was read at upstream main 0b29cd534 of 7 October 2026 (`DynamicISF.swift`, a
Swift port of the `trio-oref` JavaScript). Its logarithmic mode expresses the same
equation as a ratio of the profile sensitivity, with a default adjustment factor
of 0.8, and holds that ratio inside the autosens limits. At the default limits of
0.7 and 1.2 the result stays between the profile divided by 1.2 and the profile
divided by 0.7, which for the person above is 41.7 to 71.4 mg/dL per unit. Its
sigmoid mode equals the profile value at the profile target and is held by the
same limits. The daily dose is 0.65 of the last two hours' average plus 0.35 of
the ten day average by default.

The two systems already give different values to the same person. At 200 mg/dL
Trio's logarithmic mode gives the person above 41.7 mg/dL per unit against 29.3
from AndroidAPS, a difference of 42%.

oref0 as published has autosens and no dynamic ISF based on daily dose, so there
is nothing in it to replace. Loop and iAPS were not read.

## The published curve

The publication [1] is a conference abstract and poster from the manufacturer of
the DBLG1 closed-loop system, using records from 11,201 adults with type 1
diabetes since March 2021. No full paper was found.

Its unit of analysis is a fall in glucose. On days with at least 75% of time in
closed loop, the authors took each glucose peak between 50 and 400 mg/dL that came
at least an hour after a meal, with less than 2.5 g of carbohydrate still on
board, and followed it to the next minimum more than 20 minutes later. The
sensitivity for the event is the fall divided by the insulin used over it. Insulin
used is the drop in insulin on board plus whatever was delivered above the
reference basal.

An illustration with round numbers, which is no one's record: glucose peaks at 240
mg/dL and reaches a minimum of 160 mg/dL, and over that time 2.0 U were used. The
event's sensitivity is 80 ÷ 2.0, or 40 mg/dL per unit.

Each event ends where glucose stopped falling, so events are selected on what
glucose did. The figure is therefore a ratio of two outcomes, the fall and the
insulin used during it, which bears on how the curve is read later in this
document.

The abstract reports a sensitivity that rises from 65 to 79 mg/dL per unit between
50 and 100 mg/dL. Above that it falls, to 35 at 180 mg/dL and below 20 above 300.
Relative to a dose sized at 180 mg/dL, the authors suggest the dose could be halved
at 120 mg/dL or below and may need doubling above 340.

The poster's results figure gives the fitted curve in two pieces, for glucose G in
mg/dL. Up to 100 mg/dL:

`ISF = 98.03 - 1.077 G + 0.008868 G^2`

Above 100 and up to 380 mg/dL:

`ISF = 272 - 3.121 G + 0.01511 G^2 - 3.305e-05 G^3 + 2.69e-08 G^4`

The coefficients in the code (`inv008.dynisf.quartic_isf`) were checked against
that figure on 11 November 2026 and match digit for digit. The quartic gives 80.6
mg/dL per unit at 100 mg/dL, 35.3 at 180, 21.1 at 300 and 15.3 at 380. It implies
dose multiples of 0.56 at 120 mg/dL and 1.95 at 340, in line with the abstract's
half and double. The quadratic gives 79.0 at 100 mg/dL and 65.3 at 60, so the two
pieces differ by 2% where they meet.

Several things the model needs are absent from the publication. It gives no daily
dose for the people measured. It does not define the reference basal or say which
insulin action model produced the insulin on board. The figure's horizontal axis
is labelled G(t) with no statement of which glucose in an event that is, and the
rest of the poster could not be retrieved.

## Scaling the curve to a person

The published values are the same for everyone. At 180 mg/dL the curve gives 35
mg/dL per unit to a person on 20 U a day and to a person on 80 U a day. Entered
settings do not look like that. Among 422 adults with entered pump settings
analysed earlier in this series, sensitivity multiplied by daily dose was 1915
(95% CI 1828 to 1986), which puts the first of those people near 96 mg/dL per unit
and the second near 24. The curve needs a scale for each person before it can be
used to dose.

The publication's own reference point is 180 mg/dL, the glucose against which it
states its halving and doubling of doses. The scaling used here takes a person's
sensitivity to be 1915 ÷ TDD at that glucose and keeps the published shape on
either side of it:

sensitivity = (1915 ÷ TDD) × (published value at this glucose ÷ published value at 180 mg/dL)

Both the anchor at 180 mg/dL and the adult figure of 1915 are modelling choices.
The figure pooled across all ages in the same analysis is higher, and 51% of that
group are children, so it was not used.

For the person on 40 U a day this gives 47.9 mg/dL per unit at 180 mg/dL, and the
published shape then carries it to 109.4 at 100 mg/dL, 62.3 at 150, 41.9 at 200
and 28.7 at 300. Put another way, the published values apply unchanged to someone
whose daily dose is 54.3 U. A person on less insulin gets proportionally higher
values and a person on more gets lower ones.

The alternative is to scale by the entered profile in place of the daily dose:

sensitivity = profile sensitivity × (published value at this glucose ÷ published value at 180 mg/dL)

For the same person this gives 114.3, 65.1, 43.7 and 30.0 mg/dL per unit at 100,
150, 200 and 300 mg/dL. The two scalings agree for anyone whose profile
sensitivity multiplied by daily dose is 1915. This person's product is 2,000,
which is why their two sets of values sit about 4% apart.

The choice between the two follows the system. AndroidAPS does not use the profile
sensitivity in dynamic mode, so scaling by daily dose fits its code as it stands,
and the value keeps moving with recent insulin use as it does now. The cost is the
one people already bear: the entered sensitivity and its schedule through the day
play no part. Trio already expresses dynamic ISF as a ratio of the profile.
Scaling by the profile fits that structure and keeps the schedule, and the value
then stops responding to recent dosing. The first scaling needs a level and a
definition of daily dose to be chosen. The second needs only the glucose at which
the profile is taken to hold. This document carries scaling by daily dose as the
candidate, because the comparator is defined that way, and scaling by the profile
as an arm of the evaluation.

The published fit stops at both ends, and the reference implementation
(`inv013/candidate.py`) makes a choice at each. Below 100 mg/dL it holds the value
reached at 100. The published quadratic falls there, towards 65 mg/dL per unit at
60 mg/dL, and following it would give more insulin as glucose drops below 100.
Work in this series in March 2026 declined the quadratic for dosing on that
ground. Above 380 mg/dL, where the fit ends, the value at 380 is held.

## The level

How high the whole curve sits is the least settled part of the model, because the
publication gives nothing to set it from. One way to compare the candidates is to
ask what daily dose a person would need for the published values to apply to them
unchanged. This document calls that the reference dose.

| Basis for the level | People | Sensitivity × daily dose at 200 mg/dL (95% interval) | Reference dose, U |
|---|---:|---:|---:|
| Adult entered settings, taken to hold at 180 mg/dL (this draft) | 422 | 1,674 (1,598 to 1,736) | 54.3 |
| The 1800 rule, taken to hold at 105 mg/dL (this series, March 2026) | | 731 | 23.7 |
| Measured under Loop | 842 | 2,182 (2,065 to 2,266) | 70.8 |
| Measured under Trio | 22 | 2,273 (1,569 to 3,218) | 73.7 |
| Measured under oref0 | 86 | 2,398 (1,779 to 5,809) | 77.8 |
| Measured with no algorithm dosing (REPLACE-BG) | 196 | 3,283 (2,898 to 3,854) | 106.5 |

The first row takes a person's entered sensitivity to be right at 180 mg/dL. The
second takes the conventional rule to be right near target, and it puts the whole
curve at less than half the height of the first. The last four rows are the level
that the publication's kind of measurement returned on other records in this
series (INV-012), read at 200 mg/dL. Trio's interval, from 22 people, includes
this draft's 1,674. The other three intervals lie above it.

The measured rows are not like for like with the published curve. Their event rules
differ from the published ones, as the section on evidence sets out, and daily dose
was not defined the same way in each source. They probably sit too high as well,
since the same measurement on simulated people with a true value of 2,000 read 1.4
to 3.4 times above it.

No default can be chosen from these figures. The reference doses run from 24 to
107 U, and narrowing that range is the purpose of the second stage of the
evaluation. For any one person the existing adjustment factor moves the level in
either direction, as it does now.

## Compared with the existing equation

The table gives the sensitivity each method assigns to the worked person, in mg/dL
per unit. Their profile value is 50 throughout, their daily dose is 40 U, their
insulin is rapid-acting and their recent dosing equals their average.

| Glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Trio sigmoid | Curve scaled by daily dose | Curve scaled by profile |
|---:|---:|---:|---:|---:|---:|
| 80 | 50.1 | 62.6 | 53.7 | 109.4 | 114.3 |
| 100 | 43.4 | 54.3 | 50.0 | 109.4 | 114.3 |
| 120 | 38.9 | 48.6 | 47.1 | 86.2 | 90.1 |
| 150 | 34.2 | 42.8 | 44.3 | 62.3 | 65.1 |
| 180 | 31.0 | 41.7 | 42.9 | 47.9 | 50.0 |
| 200 | 29.3 | 41.7 | 42.4 | 41.9 | 43.7 |
| 250 | 26.3 | 41.7 | 41.8 | 33.7 | 35.1 |
| 300 | 24.1 | 41.7 | 41.7 | 28.7 | 30.0 |
| 380 | 21.8 | 41.7 | 41.7 | 20.7 | 21.7 |

The curve is steeper as well as higher. Between 100 and 200 mg/dL its value falls
by a factor of 2.61, where the existing equation falls by 1.48. From 200 to 300
mg/dL the curve's value is multiplied by 0.69 and the existing equation's by 0.82.

What that does to dosing depends on the comparator. Against stock AndroidAPS the
curve scaled by daily dose gives less insulin at every glucose up to 362 mg/dL
with a divisor of 55, and up to 313 mg/dL with a divisor of 75. The gap is widest
near target. Above those crossing points the curve gives more. Against Trio at its
defaults the curve gives less insulin below about 200 mg/dL and more above,
because Trio's upper limit has pinned its value at 41.7 by 180 mg/dL. Against this
person's own profile the curve gives less insulin below 174 mg/dL and more above.

## Where on the curve a correction is read

A correction from 200 mg/dL to a target of 100 passes through every glucose in
between, and the curve gives a different sensitivity at each. The existing code
reads its equation at the current glucose, or at a glucose chosen by one of the
three branches mentioned earlier. The publication's figure plots sensitivity
against a glucose it does not define. This draft assumes an event was filed under
the glucose halfway between its peak and its minimum, and that assumption has not
been confirmed.

The table gives the units needed to bring the worked person to 100 mg/dL under
each reading. It is the arithmetic of the sensitivity alone.

| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Profile | Curve, read at the start | Curve, read at the midpoint | Curve, summed along the path |
|---:|---:|---:|---:|---:|---:|---:|
| 150 | 1.46 | 1.17 | 1.00 | 0.80 | 0.61 | 0.62 |
| 200 | 3.41 | 2.40 | 2.00 | 2.39 | 1.61 | 1.62 |
| 250 | 5.71 | 3.60 | 3.00 | 4.46 | 3.01 | 2.97 |
| 300 | 8.29 | 4.80 | 4.00 | 6.97 | 4.78 | 4.58 |

Reading the curve at the starting glucose asks for 31 to 49% more insulin than
reading it at the midpoint. From 200 mg/dL the two readings give 2.39 and 1.61 U,
on either side of the 2.00 U the flat profile gives. Adding up the curve along the
whole path gives 1.62 U. The midpoint reading stays within 0.04 U of that sum up
to a starting glucose of 250 mg/dL and is 0.20 U above it at 300.

The two readings coincide for a flat sensitivity and move apart as a curve
steepens. The default proposed is the midpoint between the glucose being corrected
and the target, which is a modelling choice. Reading at the starting glucose is
kept as an arm of the evaluation.

## What has to be added

The glucose function is the smallest part. One function per system returns the
published quartic divided by its value at 180 mg/dL, with test vectors shared
between the two codebases so that both return the same values. The divisor for
insulin type has no counterpart in the curve and is not used in this mode.

The level needs a default, and this document offers only the provisional 1915. The
daily dose it is divided by has to be defined as well, since the publication is
silent on it. The default proposed is each system's existing blend, so that this
half of the behaviour stays as people know it. With that blend the value moves
with recent dosing by the same multiplier as today. It is multiplied by 1.39 for a
person whose blended dose is 0.72 of their own average and by 0.88 at 1.14 of it,
which are the tenth and ninetieth percentiles among 960 people in document 6. A
plain average is an arm of the evaluation. The adjustment factor stays. Its
guidance has to be worked out again, because at 100% the new value at 100 mg/dL is
2.52 times the existing one.

A rule for where the curve is read has to be written into the dose calculation, as
the previous section describes. The treatment below 100 and above 380 mg/dL has to
be fixed too. Holding the end values is the default. Two alternatives below 100
mg/dL go into the evaluation: the published quadratic, and the rise of the
existing logarithm.

Trio's limits need a decision from its maintainers. For the worked person the
curve runs from 20.7 to 109.4 mg/dL per unit, and Trio's default limits would hold
it between 41.7 and 71.4. That leaves the curve free to move only between 137 and
201 mg/dL, so with the limits unchanged most of it would never reach the dose.

The forecasts are affected as well as the dose. Both systems calculate a
sensitivity at each predicted glucose, and the call sites can stay as they are.
The curve's value rises 2.61 times as predicted glucose falls from 200 to 100
mg/dL, against 1.48 for the existing equation, so a forecast heading towards
target is computed with a higher sensitivity than now. The effect on minimum
predicted glucose and on time spent at zero basal has not been measured and is
part of the replay.

The remaining additions are smaller. The loop result already stores the
sensitivity it used (`variableSens` in AndroidAPS), and during evaluation both the
existing and the new value are stored while the existing one doses. Settings gain
a mode selector and, in a first version, no new number for the person to set. Each
system gets unit tests on the shared reference vectors and replay assets for the
new mode, alongside its existing tests (`ReplayApsResultsTest` in AndroidAPS,
`DynamicISFTests` in Trio). A guard test confirms that the existing mode's output
is unchanged.

In AndroidAPS the code touched is `calculateRawDynIsf` and the preference screen,
together with the per-prediction lines and the dose-sensitivity branches in
`DetermineBasalSMB.kt`. In Trio it is `DynamicISF.swift` and the forecast
adjustment used in logarithmic mode. The `trio-oref` JavaScript that the Swift was
ported from would need the same change.

## What is known about the curve

The construction was repeated in this series as INV-012, on 1,790 people in eight
cohorts, and is reported in section 7 of document 5. It was built before the
abstract was read, and two of its event rules differ from the published ones. It
required three hours clear of a meal where the publication requires one hour with
under 2.5 g on board, and a peak standing 15 mg/dL clear of its surroundings where
the publication requires a 20 minute interval. It also multiplied each event by the
person's daily dose. The comparisons below are of ratios between glucose levels,
which that multiplication does not affect.

The shape is reproducible. Within person, the value at 100 mg/dL was 2.2 to 4.0
times the value at 200 mg/dL across the eight cohorts, against 2.61 in the
published curve. In REPLACE-BG, where no algorithm was dosing, the ratio was 2.55
(95% interval 2.25 to 2.82), which includes the published figure.

A constant sensitivity produces a falling curve by the same method. Simulated
people with a true value of 2,000 at every glucose returned values 1.4 to 3.4
times too high, declining as glucose rose. Their ratio of the value at 250 to the
value at 200 mg/dL was 0.87 (0.83 to 0.90) in one arm and 0.92 (0.89 to 0.95) in
the other. The six Jaeb cohorts gave 0.84 to 1.05, a range that contains both
simulated arms, so the decline in those cohorts between 150 and 250 mg/dL needs no
change in what a unit of insulin does to explain it. The published curve gives
0.80, which is steeper than either simulated arm by an amount this series cannot
account for. The simulation produced no usable events at 100 or at 300 mg/dL and
says nothing about either end.

Near target the decline comes from the insulin, with the fall almost unchanged.
Across cohorts the median fall was 30 to 43 mg/dL for events at 100 mg/dL and 39
to 64 for events at 200. The median insulin used was 0.02 to 0.2 U at 100 mg/dL
and 0.3 to 1.5 U at 200. A similar fall divided by much less insulin gives a much
higher figure. At 100 mg/dL between 9 and 58% of events had no net insulin used at
all, and a line fitted through each person's events put the fall at zero net
insulin at 27 to 44 mg/dL. In type 1 diabetes a raised glucose does not fall
without insulin, so that fall is basal and earlier doses acting, together with the
selection of events that end at a minimum. The steep part of the curve below 150
mg/dL therefore rests on very small net doses, and it has not been tested against
a known truth.

Above 200 mg/dL the question is open. Measured as the extra fall per extra unit,
the value at 300 against 200 mg/dL was 0.66 (0.49 to 0.73) in REPLACE-BG, from 24
people, and 0.73 (0.64 to 0.93) under Loop, from 49 people. The published curve
gives 0.69. Both intervals exclude 1, so a decline above 200 mg/dL is seen on that
measure. No simulated arm reached 300 mg/dL, and the measure therefore has no
known-truth check there either.

The shape has been tried once as a predictor on oref records, where it did worse
than the static profile. On 80 people and 62,344 overnight windows, with each
shape anchored to the profile at 100 mg/dL, the typical miss was 24.2 mg/dL for
the static profile and 29.8 for this shape. At 205 to 260 mg/dL it was 37.8
against 73.7. No interval was reported for those figures. Two things limit the
test. Its windows were selected on the later outcome, as the audit of 16 September
2026 recorded, and it anchored the shape at 100 mg/dL where this draft anchors at
180.

All of this is observational, and the publication is a conference abstract from
the system's manufacturer. As a description of net fall per net unit in events
chosen from peak to minimum, the curve is reproducible on other people's data. As
a statement that a unit of insulin does less at higher glucose, it is unsupported
between 150 and 250 mg/dL and open on either side of that range. Whether a loop
that doses to it does better than one dosing to the existing equation or to a flat
profile is a separate question. No evidence in this series bears on it directly.

A proposal from June 2026 that used this shape (document 3) was withdrawn for live
dosing on 16 September 2026. The audit of that date left the way open for such
candidates to be simulated and replayed in shadow, and to be studied
prospectively. The evaluation below stays inside that.

## Evaluation

Seven arms are compared throughout. Two are references: the existing equation as
shipped, and the profile alone. Four are forms of the curve. They are the curve
scaled by daily dose and read at the midpoint, the same read at the starting
glucose, the same on a plain average daily dose, and the curve scaled by the
profile. The seventh is the first draft of this document, a sensitivity anchored
to the profile and flat above 100 mg/dL. It is the arm the curve has to beat for
its shape to matter.

The first stage replays each arm's value on recorded histories using the INV-008
pipeline. Results are split by glucose and by time since the last meal, and
reported separately for people who announce meals and people who do not. The stage
shows how much insulin each arm asks for and in whom. It is descriptive and has no
pass rule.

The second stage runs the publication's measurement (`inv012/peak_trough.py`) on
AndroidAPS and Trio records under the published event rules, and extends the
simulated constant-sensitivity arms to reach 100 and 300 mg/dL. It decides the
default level. It also shows whether the decline near target and above 250 mg/dL
appears with no real effect behind it, as it does in between.

The third stage is simulation in TimSim, with runs paired by seed and arms for
unannounced meals and for a mis-set adjustment factor. A harm rule is written
before the runs and has two sides, time below 54 mg/dL and time above 250 mg/dL.
The second side is there because this model gives less insulin than stock
AndroidAPS over most of the range.

The fourth stage is shadow running. Both values are calculated on the device each
cycle and stored, and the existing one doses. That gives the distribution of the
ratio between them for each person, and from it a starting adjustment factor.

The fifth stage is use, on the developer's own device first and on released lines
after that. It follows only when the third and fourth stages have passed, and it
stops on the harm rule.

The first two stages are analysis and carry the analysis checklist. The last two
change code that can change a dose. They carry unit and guard tests, the
repository's review checklist, a replay against real history and a fresh-context
review of the diff.

The way back at every stage is the mode selector. The existing equation stays in
the code unchanged and remains the default until the fifth stage has run, and any
build that includes the new mode names the build before it.

## What is not settled

Some of what is missing can only come from the authors. The level depends on four
things the publication does not give: which glucose an event is filed under, how
the reference basal was defined, which insulin action model gave the insulin on
board, and what daily doses the people measured were using. The authors have not
been asked. They have also not been asked whether a curve derived by a
manufacturer from its own users' records may be built into open-source code.

The level is open, as described above. A level set from entered settings at 180
mg/dL may not apply to open-source systems, since the measurement on Loop and
oref0 records read higher at 200 mg/dL and only Trio's wide interval includes it.

It is an open question whether the steep part between 100 and 150 mg/dL should be
used for dosing at all, given how many events there involved no net insulin. The
question of whether to use the published quadratic below 100 mg/dL belongs with
it.

Parts of the existing systems were not read. The bolus calculator's sensitivity in
dynamic mode is one, and the behaviour of stock autosens in that mode is another.
TimSim's oref adapter has not been checked for whether it runs the existing
dynamic ISF closely enough to act as the comparator in the third stage.

The last question is for the maintainers of AndroidAPS and Trio, and it is whether
they want one shared definition. Their two implementations of the existing
equation already differ by 42% for the worked person at 200 mg/dL, and a reference
implementation with shared test vectors would be new to both.

## Reproducing the numbers

`python3 -m inv013.candidate` writes `results/inv013_comparator.md`, which holds
every number in this document's tables and worked examples.
`inv013/tests/test_candidate.py` checks the stock equation against the replay code
of INV-008, and the curve against the values in the abstract and on the poster's
figure. It also checks the ways of reading the curve.

Figures from other analyses are quoted from where they were produced. The
reproduction is in `results/inv012_peak_trough.md` and section 7 of document 5.
The adult figure of 1915 is from `docs/preprint-dosing-constants.md`, and the
overnight test is in `oref-investigations/results/bridge_diabeloop.md`.

Earlier drafts of this document are in git. The first (commit fd3662e) proposed
the profile-anchored, flat form now kept as an arm.

## Reference

1. Romero-Ugalde HM, Adenis A, Cartier-Michaud T, Gauthier P, Desir C, Le
   Roux-Mallouf T, Huneker E, Benhamou PY. 920-P: The Relationship between
   Glycemia and Insulin Sensitivity Factor. Diabetes 2025;74(Supplement 1):920-P.
   doi:10.2337/db25-920-P
