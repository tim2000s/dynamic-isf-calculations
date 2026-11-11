---
title: "An alternative proposal for dynamic ISF in open-source AID systems"
subtitle: "The Diabeloop glucose curve, tested on other cohorts and set against the logarithmic and sigmoid models"
author: "Tim Street, Diabettech"
date: "11 November 2026"
---

## Summary

The dynamic ISF models in AndroidAPS and Trio have two origins. The logarithmic
model began with Chris Wilson, who described the work in a talk to Loop and Learn
in February 2022 [2]. Through an online form he asked volunteers from the diabetes
online community what correction dose they would give at a stated glucose, with
glucose steady and no insulin on board. About 50 people answered, from a
two-year-old using 8 U a day to people in their seventies using more than 100 U,
and the usable answers gave about 300 data points. The doses and their effects were
reported by the respondents from experience, and neither was measured. He trained a neural network on those points and
gave the surface it found across glucose and TDD as roughly 277,700 ÷ (glucose ×
TDD). AndroidAPS first carried that equation, from February 2022, and the glucose
term was replaced by a logarithm in June 2022, which is the form in use now [3].
The talk predates that change and does not describe it.

The sigmoid model offered in Trio came from its author's observation of glucose
response, and no analysis of data has been published in support of it.

Diabeloop has since examined the same idea on the records of its own closed-loop
system and presented the result at the American Diabetes Association's Scientific
Sessions in 2025 [1]. From 11,201 adults it reported a sensitivity of 79 mg/dL per
unit at a glucose of 100 mg/dL, 35 at 180 and below 20 above 300, with a fitted
equation.

This document builds on that work. It sets out the equation and how it was
derived, and reports what the same measurement shows on other cohorts. It then
describes how the equation could be offered in AndroidAPS and Trio as an
alternative to both the logarithmic and the sigmoid models.

The measurement was repeated on 1,790 people in eight cohorts held in this series.
Within person, the value at 100 mg/dL was 2.2 to 4.0 times the value at 200,
against 2.61 in the published equation, so the steep fall towards 200 mg/dL is
reproducible. On records where no algorithm was dosing, the published equation
fits the measured values from 100 to 200 mg/dL at a single scale. On Loop records
no single scale fits. Refitting the equation's five coefficients to each cohort
agreed with both results, and returned a flatter curve than the published one
above 200 mg/dL. Simulated people whose sensitivity was constant also
returned a falling curve by this method, so the decline between 150 and 250 mg/dL
is not evidence that a unit of insulin does less at higher glucose.

The published equation has no term for the person, and it needs one before it can
dose. Two are described here: total daily dose (TDD), as the logarithmic model
uses, and the profile ISF. For an adult on 40 U a day with a profile ISF of 50
mg/dL per unit, the equation scaled by TDD gives 109.4 mg/dL per unit at 100 mg/dL
and 41.9 at 200. AndroidAPS at an adjustment factor of 100% gives that person 43.4
and 29.3. Both scalings are modelling choices made in this document and appear
nowhere in the publication.

Two parts of the model are unsettled. The first is its level, since the
publication reports no daily dose and the candidates in this series differ by
more than a factor of four. The second is where on the curve a correction is read, which moves
the dose by 31 to 49%. Whether a loop dosing to this equation does better than one
dosing to either existing model has not been tested.

This is a draft for discussion. The equation has not run on a device or in a
closed-loop simulator, and the candidate uses the published coefficients
unchanged. Readers are assumed to
know how dynamic ISF works in both systems. Glucose is in mg/dL throughout because
the code is.

## What is the Diabeloop equation?

The source [1] is a conference abstract and poster from Diabeloop, the maker of
the DBLG1 closed-loop system. It used the records of 11,201 adults with type 1
diabetes from March 2021 onwards. No full paper was found.

The abstract reports a sensitivity that rises from 65 to 79 mg/dL per unit between
50 and 100 mg/dL. Above that it falls, to 35 at 180 mg/dL and below 20 above 300.
Relative to a dose sized at 180 mg/dL, the authors suggest the dose could be halved
at 120 mg/dL or below and may need doubling above 340.

The poster's results figure gives the fit in two pieces, for glucose G in mg/dL.
Up to 100 mg/dL:

`ISF = 98.03 - 1.077 G + 0.008868 G^2`

Above 100 and up to 380 mg/dL:

`ISF = 272 - 3.121 G + 0.01511 G^2 - 3.305e-05 G^3 + 2.69e-08 G^4`

The coefficients in the code (`inv008.dynisf.quartic_isf`) were checked against
that figure on 11 November 2026 and match digit for digit. The quartic gives 80.6
mg/dL per unit at 100 mg/dL, 35.3 at 180, 21.1 at 300 and 15.3 at 380. It implies
dose multiples of 0.56 at 120 mg/dL and 1.95 at 340, in line with the abstract's
half and double. The quadratic gives 79.0 at 100 mg/dL and 65.3 at 60, so the two
pieces differ by 2% where they meet.

The equation is steeper than the logarithmic model. Between 100 and 200 mg/dL its
value falls by a factor of 2.61, where the logarithmic model with a divisor of 55
falls by 1.48. From 200 to 300 mg/dL its value is multiplied by 0.69 and the
logarithmic model's by 0.82.

It has no term for the person. Everyone gets 35.3 mg/dL per unit at 180 mg/dL,
whatever their daily dose or their profile.

## How was it derived?

The unit of analysis is a fall in glucose. On days with at least 75% of time in
closed loop, the authors took each glucose peak between 50 and 400 mg/dL that came
at least an hour after a meal, with less than 2.5 g of carbohydrate still on
board, and followed it to the next minimum more than 20 minutes later. The
sensitivity for the event is the fall divided by the insulin used over it. Insulin
used is the drop in insulin on board plus whatever was delivered above the
reference basal. The curve is the fit through those event values against glucose.

Each event ends where glucose stopped falling, so events are selected on what
glucose did. The figure is therefore a ratio of two outcomes, the fall and the
insulin used during it. That bears on how the curve is read in the next section.

Four things the derivation depends on are absent from the publication. It does not
define the reference basal or say which insulin action model produced the insulin
on board. The figure's horizontal axis is labelled G(t) with no statement of which
glucose in an event that is. It gives no daily dose for the people measured. The
rest of the poster could not be retrieved, and the authors have not been asked.

## What does it show on our data?

The same construction was run in this series as INV-012, on 1,790 people in eight
cohorts, and is reported in section 7 of document 5. It was built before the
abstract was read, and two of its event rules differ from the published ones. It
required three hours clear of a meal where the publication requires one hour with
under 2.5 g on board, and a peak standing 15 mg/dL clear of its surroundings where
the publication requires a 20 minute interval. Each event was filed under the
glucose halfway between its peak and its minimum, and multiplied by the person's
daily dose.

### The shape between 100 and 200 mg/dL

The steep fall is reproducible. Within person, the value at 100 mg/dL was 2.2 to
4.0 times the value at 200 mg/dL across the eight cohorts, against 2.61 in the
published equation. In REPLACE-BG, where no algorithm was dosing, the ratio was
2.55 (95% interval 2.25 to 2.82), which includes the published figure. Under
closed loop it was steeper: 3.26 (3.11 to 3.39) for Loop, 3.86 (2.72 to 6.04) for
oref0 and 3.21 (2.24 to 5.62) for Trio. The first two intervals lie above 2.61 and
Trio's includes it.

### Fitting the equation to the measured values

The table divides the measured value in each glucose bin by the published value at
that glucose. The result is the daily dose, in U, at which the published equation
would return what was measured. A row that stays level means the published shape
fits that cohort at one scale. Only the four cohorts with carbohydrate records are
shown, with the 95% interval over people.

| Cohort | People | 100 mg/dL | 150 | 200 | 250 | 300 |
|---|---:|---:|---:|---:|---:|---:|
| No algorithm (REPLACE-BG) | 196 | 103 (89 to 115) | 100 (94 to 108) | 106 (94 to 125) | 120 (103 to 137) | 174 (132 to 233) |
| Loop | 842 | 88 (84 to 91) | 58 (56 to 60) | 71 (67 to 73) | 82 (77 to 89) | 108 (95 to 119) |
| Trio | 22 | 90 (63 to 111) | 52 (40 to 64) | 74 (51 to 104) | | |
| oref0 | 86 | 109 (90 to 136) | 70 (59 to 86) | 78 (58 to 188) | | |

With no algorithm dosing, the published equation fits from 100 to 200 mg/dL at one
scale: 103, 100 and 106 U, each inside the other two intervals. Above 200 mg/dL
the row rises to 120 and then 174 U, so the published equation falls faster there
than these records do. Under Loop no single scale fits. The figure is 88 U at 100
mg/dL, 58 at 150 and 71 at 200, and those three intervals do not overlap. Trio and
oref0 show the same dip at 150 mg/dL with wider intervals.

Two cautions apply to the table. The columns hold different people, since fewer
reach the higher bins, and the within-person ratios in this section are the firmer
comparison. The published values were read at the centre of each bin on the
assumption that the publication files an event under its midpoint glucose, which
it does not state.

The scale itself is not a daily dose anyone was using. It is about 100 U in
REPLACE-BG because the measurement reads high, as the simulation below shows.

### Refitting the equation

The equation was also refitted. The quartic, with all five coefficients free, was
fitted to each cohort's values in 10 mg/dL bins, from 100 mg/dL to the highest bin
holding ten people. Each bin was weighted by the square root of the people in it,
which is a modelling choice, since the poster does not state its weighting. For
comparison the published fit was laid over the same bins with one free number, its
scale in U. The miss is the typical distance between a fit and the bin values, as
a percentage of each value.

| Cohort | People | Bins fitted, mg/dL | Best scale of the published fit, U | Miss, published fit | Miss, refit |
|---|---:|---|---:|---:|---:|
| No algorithm (REPLACE-BG) | 196 | 100 to 370 | 103 (96 to 113) | 14% | 7% |
| Loop | 842 | 100 to 380 | 71 (69 to 74) | 17% | 6% |
| Trio | 22 | 100 to 200 | 65 (48 to 80) | 18% | 7% |
| oref0 | 86 | 100 to 240 | 82 (71 to 99) | 23% | 14% |

The refit misses by less in every cohort, as a fit with four more free numbers
must. The second table shows where the two differ. It gives each fit relative to
its own value at 180 mg/dL.

| Fit | People | 100 mg/dL | 150 | 200 | 250 | 300 |
|---|---:|---:|---:|---:|---:|---:|
| Published | | 2.29 | 1.30 | 0.87 | 0.70 | 0.60 |
| Refit, no algorithm (REPLACE-BG) | 196 | 2.17 (1.91 to 2.40) | 1.28 (1.21 to 1.33) | 0.89 (0.86 to 0.92) | 0.81 (0.72 to 0.87) | 0.90 (0.77 to 1.02) |
| Refit, Loop | 842 | 3.19 (3.02 to 3.34) | 1.25 (1.22 to 1.27) | 1.00 (0.98 to 1.01) | 1.07 (1.03 to 1.11) | 0.94 (0.88 to 1.01) |
| Refit, Trio | 22 | 3.29 (2.15 to 4.76) | 1.19 (0.85 to 1.53) | 1.12 (0.66 to 1.61) | | |
| Refit, oref0 | 86 | 2.96 (2.22 to 3.83) | 1.00 (0.78 to 1.24) | 1.02 (0.82 to 1.38) | | |

With no algorithm dosing, the refit returns the published shape from 100 to 200
mg/dL, and all three intervals there include the published figure. At 250 and 300
mg/dL the refit sits at 0.81 and 0.90 of its value at 180, where the published fit
is at 0.70 and 0.60. Both of those intervals exclude the published figure. Under
Loop the refit differs from the published fit at every glucose in the table. It is
higher at 100 mg/dL, 3.19 against 2.29, and from 200 mg/dL upward it stays close
to its value at 180. Trio's intervals, from 22 people, include the published
figure at each of its three points. The oref0 values resemble Loop's, with an
interval that excludes the published figure only at 150 mg/dL.

Two limits apply. A quartic suits these records poorly, because the values fall
steeply to about 150 mg/dL and then level off. Under Loop the refit rises from
1.00 to 1.07 between 200 and 250 mg/dL, while the bin values themselves go from
2,182 to 2,025. Above 200 mg/dL the bins and the within-person ratios reported
below are the better guide. The second limit is that the refit describes the same
outcome ratio as the published curve, so the simulation in the next subsection
bears on it equally.

The refit gives no single equation to offer in place of the published one. Its
coefficients differ between cohorts, and its value at 180 mg/dL is 3,657 (3,431 to
4,055) with no algorithm dosing and 2,123 (2,048 to 2,214) under Loop. It confirms
the published shape between 100 and 200 mg/dL on records that no loop shaped. It
also agrees with the within-person ratios that, by this measure, the published
fall above 200 mg/dL is not reproduced here.

### What a constant sensitivity returns

A constant sensitivity produces a falling curve by the same method. Simulated
people with a true value of 2,000 at every glucose returned values 1.4 to 3.4
times too high, declining as glucose rose. Their ratio of the value at 250 to the
value at 200 mg/dL was 0.87 (0.83 to 0.90) in one arm and 0.92 (0.89 to 0.95) in
the other. The six Jaeb cohorts gave 0.84 to 1.05, a range that contains both
simulated arms, so the decline in those cohorts between 150 and 250 mg/dL needs no
change in what a unit of insulin does to explain it. The published equation gives
0.80, which is steeper than either simulated arm by an amount this series cannot
account for. The simulation produced no usable events at 100 or at 300 mg/dL and
says nothing about either end.

### Near target

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

### Above 200 mg/dL

The two measures disagree here. By the ratio of fall to insulin used, the value at
300 against 200 mg/dL within person was 0.90 (0.77 to 1.26) in REPLACE-BG, from 60
people, and 0.96 (0.84 to 1.08) under Loop, from 137 people. Both intervals
exclude the published 0.69. Measured as the extra fall per extra unit, the same
comparison gave 0.66 (0.49 to 0.73) in REPLACE-BG, from 24 people, and 0.73 (0.64
to 0.93) under Loop, from 49 people, and both of those include 0.69. No simulated
arm reached 300 mg/dL, so neither measure has a known-truth check there.

### As a predictor

The shape has been tried once as a predictor on oref records, where it did worse
than the static profile. On 80 people and 62,344 overnight windows, with each
shape anchored to the profile at 100 mg/dL, the typical miss was 24.2 mg/dL for
the static profile and 29.8 for this shape. At 205 to 260 mg/dL it was 37.8
against 73.7. No interval was reported for those figures. Two things limit the
test. Its windows were selected on the later outcome, as the audit of 16 September
2026 recorded, and it anchored the shape at 100 mg/dL where this document anchors
at 180.

### Where that leaves the equation

All of this is observational, and the publication is a conference abstract from
the system's manufacturer. As a description of net fall per net unit in events
chosen from peak to minimum, the equation is reproducible on other people's data
between 100 and 200 mg/dL. As a statement that a unit of insulin does less at
higher glucose, it is unsupported between 150 and 250 mg/dL and open on either
side of that range. Whether a loop that doses to it does better than
one dosing to either existing model or to a flat profile is a separate question. No
evidence in this series bears on it directly.

A proposal from June 2026 that used this shape (document 3) was withdrawn for live
dosing on 16 September 2026. The audit of that date left the way open for such
candidates to be simulated and replayed in shadow, and to be studied
prospectively. The evaluation proposed below stays inside that.

## How could it be implemented?

### What it would replace

AndroidAPS and Trio run the same logarithmic equation:

sensitivity = 1800 ÷ (adjustment factor × TDD × ln(glucose ÷ divisor + 1))

Trio computes a ratio, profile ISF × adjustment factor × TDD × ln(glucose ÷ divisor
+ 1) ÷ 1800, and divides the profile ISF by it, so the profile value cancels. Given
the same inputs the two systems return the same number. This was checked on 11 November 2026 against AndroidAPS dev 0c593a4b09 and
Trio dev 902be75e8, both of 10 November 2026. Neither calculation differs from the
released versions (AndroidAPS master f77abd9b82, Trio main 0b29cd534), and Trio's
`DynamicISF.swift` was last changed on 24 August 2026.

What separates the two in use is the defaults around the equation. Trio's
adjustment factor defaults to 0.8 where AndroidAPS defaults to 100%, which puts
Trio's value 25% higher. Trio then holds its ratio inside the autosens limits, 0.7
and 1.2 by default, and `calculateRawDynIsf` in AndroidAPS applies no limit. For
the worked person Trio's upper limit binds from 157 mg/dL. Above that glucose its
value stays at 41.7 mg/dL per unit while AndroidAPS continues down to 29.3 at 200
mg/dL. The daily dose is blended differently in each system, and the divisors
differ for faster insulins: 55, 65 or 75 by peak time in AndroidAPS, and 55 or 70
in Trio, or 120 less a custom peak time.

Trio's sigmoid mode returns the profile ISF at the profile target and moves it
within the same limits. oref0 as published has autosens and no dynamic ISF based
on daily dose. Loop and iAPS were not read.

### Linked to TDD, or to the profile ISF?

The published values need a scale for each person. At 180 mg/dL the equation gives
35.3 mg/dL per unit to a person on 20 U a day and to a person on 80 U a day. Among
422 adults with entered pump settings analysed earlier in this series, sensitivity
multiplied by daily dose was 1915 (95% CI 1828 to 1986), which puts the first of
those people near 96 mg/dL per unit and the second near 24.

The publication's own reference point is 180 mg/dL, the glucose against which it
states its halving and doubling of doses. Both links below hold a value at that
glucose and keep the published shape on either side of it.

Linked to TDD, a person's sensitivity at 180 mg/dL is taken to be 1915 ÷ TDD:

sensitivity = (1915 ÷ TDD) × (published value at this glucose ÷ published value at 180 mg/dL)

This is the logarithmic model with its glucose term swapped. The daily dose is
each system's existing blend and the adjustment factor multiplies it as now. The
divisor for insulin type has no counterpart and drops out. The value keeps moving
with recent insulin use by the same multiplier as today: 1.39 for a person whose
blended dose is 0.72 of their own average and 0.88 at 1.14 of it, which are the
tenth and ninetieth percentiles among 960 people in document 6. The profile ISF
and its schedule through the day play no part, as in AndroidAPS now. The adult
figure of 1915 is a modelling choice. The figure pooled across all ages in the
same analysis is higher, and 51% of that group are children, so it was not used.

Linked to the profile, the entered ISF is taken to hold at 180 mg/dL:

sensitivity = profile ISF × (published value at this glucose ÷ published value at 180 mg/dL)

This fits Trio's structure, which already expresses dynamic ISF as a ratio of the
profile. The schedule through the day is kept and no level has to be chosen. The
value stops responding to recent dosing. It is also only as good as the entered
number. It suits a profile ISF the person has
evaluated, and it carries the whole of any error in one they have not.

The two links agree for anyone whose profile ISF multiplied by daily dose is 1915.
The worked person has a profile ISF of 50 mg/dL per unit, a daily dose of 40 U,
rapid-acting insulin and recent dosing equal to their average. Their product is
2,000, which is why their two sets of values sit about 4% apart. The table gives
the sensitivity each method assigns them, in mg/dL per unit, with both systems at
their defaults.

| Glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Trio sigmoid | Equation linked to TDD | Equation linked to profile |
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

Against AndroidAPS at 100% the equation linked to TDD gives less insulin at every
glucose up to 362 mg/dL with a divisor of 55, and up to 313 mg/dL with a divisor
of 75. The ratio between the two is 2.52 at 100 mg/dL, 1.82 at 150, 1.43 at 200
and 1.19 at 300. Against Trio at its defaults, in either mode, the equation gives
less insulin below about 200 mg/dL and more above. Against this person's own
profile it gives less insulin below 174 mg/dL and more above.

The published fit stops at both ends, and the reference implementation
(`inv013/candidate.py`) makes a choice at each. Below 100 mg/dL it holds the value
reached at 100. The published quadratic falls there, towards 65 mg/dL per unit at
60 mg/dL, and following it would give more insulin as glucose drops below 100.
Work in this series in March 2026 declined the quadratic for dosing on that
ground. Above 380 mg/dL, where the fit ends, the value at 380 is held.

This document carries the link to TDD as the candidate, because the logarithmic
model is defined that way, and the link to the profile as an arm of the
evaluation.

### The level

The link to TDD needs a level, and it is the least settled part of the model
because the publication gives nothing to set it from. One way to compare the
candidates is the daily dose at which the published values would apply unchanged,
called the reference dose here.

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
curve at less than half the height of the first. The last four rows are the
measured values from the previous section, read at 200 mg/dL. Trio's interval,
from 22 people, includes this draft's 1,674. The other three intervals lie above
it. The measured rows probably sit too high, since the same measurement on
simulated people with a true value of 2,000 read 1.4 to 3.4 times above it, and
daily dose was not defined the same way in each source.

No default can be chosen from these figures. The reference doses run from 24 to
107 U, and narrowing that range is the purpose of the second stage of the
evaluation. For any one person the adjustment factor moves the level in either
direction, as it does now. Its guidance has to be worked out again, because at
100% the new value at 100 mg/dL is 2.52 times the existing one.

### Where on the curve a correction is read

A correction from 200 mg/dL to a target of 100 passes through every glucose in
between, and the equation gives a different sensitivity at each. The existing code
reads the logarithmic model at the current glucose, or at a glucose chosen by one
of three branches in `DetermineBasalSMB.kt`. This draft assumes the publication
filed each event under its midpoint glucose, and that assumption has not been
confirmed.

The table gives the units needed to bring the worked person to 100 mg/dL under
each reading. It is the arithmetic of the sensitivity alone, before insulin on
board and the loop's predictions change it.

| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Profile | Equation, read at the start | Equation, read at the midpoint | Equation, summed along the path |
|---:|---:|---:|---:|---:|---:|---:|
| 150 | 1.46 | 1.17 | 1.00 | 0.80 | 0.61 | 0.62 |
| 200 | 3.41 | 2.40 | 2.00 | 2.39 | 1.61 | 1.62 |
| 250 | 5.71 | 3.60 | 3.00 | 4.46 | 3.01 | 2.97 |
| 300 | 8.29 | 4.80 | 4.00 | 6.97 | 4.78 | 4.58 |

Reading the equation at the starting glucose asks for 31 to 49% more insulin than
reading it at the midpoint. From 200 mg/dL the two readings give 2.39 and 1.61 U,
on either side of the 2.00 U the flat profile gives. Adding up the equation along
the whole path gives 1.62 U. The midpoint reading stays within 0.04 U of that sum
up to a starting glucose of 250 mg/dL and is 0.20 U above it at 300.

The default proposed is the midpoint between the glucose being corrected and the
target, which is a modelling choice. Reading at the starting glucose is kept as an
arm of the evaluation.

### What changes in the code

The glucose function is the smallest part. One function per system returns the
published quartic divided by its value at 180 mg/dL, with test vectors shared
between the two codebases so that both return the same values.

Trio's limits need a decision from its maintainers. For the worked person the
equation runs from 20.7 to 109.4 mg/dL per unit, and Trio's default limits would
hold it between 41.7 and 71.4. That leaves it free to move only between 137 and
201 mg/dL, so with the limits unchanged most of the curve would never reach the
dose.

The forecasts are affected as well as the dose. Both systems calculate a
sensitivity at each predicted glucose, and the call sites can stay as they are.
The equation's value rises 2.61 times as predicted glucose falls from 200 to 100
mg/dL, against 1.48 for the logarithmic model, so a forecast heading towards
target is computed with a higher sensitivity than now. The effect on minimum
predicted glucose and on time spent at zero basal has not been measured and is
part of the replay.

The remaining additions are smaller. The loop result already stores the
sensitivity it used (`variableSens` in AndroidAPS), and during evaluation both the
existing and the new value are stored while the existing one doses. Settings gain
a mode selector and, in a first version, no new number for the person to set. Each
system gets unit tests on the shared reference vectors and replay assets for the
new mode, alongside its existing tests (`ReplayApsResultsTest` in AndroidAPS,
`DynamicISFTests` in Trio). A guard test confirms that the existing modes' output
is unchanged.

In AndroidAPS the code touched is `calculateRawDynIsf` and the preference screen,
together with the per-prediction lines and the dose-sensitivity branches in
`DetermineBasalSMB.kt`. In Trio it is `DynamicISF.swift` and the forecast
adjustment used in logarithmic mode.

## Evaluation

Eight arms are compared throughout. Three are references. The logarithmic model
as shipped and Trio's sigmoid are the two existing models, and the profile alone
is the baseline for both. Four arms are forms of the equation: linked to TDD and read at the midpoint, the same read at the
starting glucose, the same on a plain average daily dose, and linked to the
profile. The eighth is the first draft of this document, a sensitivity anchored to
the profile and flat above 100 mg/dL. It is the arm the equation has to beat for
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
The second side is there because this model gives less insulin than AndroidAPS at
100% over most of the range.

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

The way back at every stage is the mode selector. The existing modes stay in the
code unchanged and remain the default until the fifth stage has run, and any build
that includes the new mode names the build before it.

## What is not settled

Some of what is missing can only come from the authors. The level depends on four
things the publication does not give: which glucose an event is filed under, how
the reference basal was defined, which insulin action model gave the insulin on
board, and what daily doses the people measured were using. The authors have not
been asked. They have also not been asked whether an equation derived by a
manufacturer from its own users' records may be built into open-source code.

The level is open, as described above. A level set from entered settings at 180
mg/dL may not apply to open-source systems, since the measurement on Loop and
oref0 records read higher at 200 mg/dL and only Trio's wide interval includes it.

It is an open question whether the steep part between 100 and 150 mg/dL should be
used for dosing at all, given how many events there involved no net insulin. The
question of whether to use the published quadratic below 100 mg/dL belongs with
it.

The fall above 200 mg/dL is open in the same way. By the ratio of fall to insulin
used, neither REPLACE-BG nor Loop reproduces it, whichever way that ratio is
summarised. By the extra fall per extra unit, both do.

Parts of the existing systems were not read. The bolus calculator's sensitivity in
dynamic mode is one, and the behaviour of autosens in that mode is another.
TimSim's oref adapter has not been checked for whether it runs the logarithmic
model closely enough to act as the comparator in the third stage.

The last question is for the maintainers of AndroidAPS and Trio, and it is whether
they want one shared definition. The two systems run the same equation and, at
their defaults, return values 42% apart for the worked person at 200 mg/dL. A
reference implementation with shared test vectors would be new to both.

## Reproducing the numbers

`python3 -m inv013.candidate` writes `results/inv013_comparator.md`, which holds
the numbers in this document's worked examples and in every table except the
refit. `python3 -m inv013.refit` writes `results/inv013_refit.md`, with the refit
tables and the refitted coefficients. It reads the per-person bins that `python3
-m inv012.peak_trough` writes, a file too large to keep in git, and checks them
against the committed INV-012 table before fitting.
`inv013/tests/test_candidate.py` checks the logarithmic equation against the
replay code of INV-008, Trio's form against AndroidAPS's, and the published fit
against the values in the abstract and on the poster's figure. It also checks the
ways of reading the curve, and that the refit returns a quartic it was given.

Figures from other analyses are quoted from where they were produced. The
measurement on other cohorts is in `results/inv012_peak_trough.md` and section 7
of document 5. The adult figure of 1915 is from `docs/preprint-dosing-constants.md`,
and the overnight test is in `oref-investigations/results/bridge_diabeloop.md`.

Earlier drafts of this document are in git. The first (commit fd3662e) proposed
the profile-anchored, flat form now kept as an arm. The
second (commit b324528) opened on the comparison with AndroidAPS and described both
systems for a general reader.

## References

1. Romero-Ugalde HM, Adenis A, Cartier-Michaud T, Gauthier P, Desir C, Le
   Roux-Mallouf T, Huneker E, Benhamou PY. 920-P: The Relationship between
   Glycemia and Insulin Sensitivity Factor. Diabetes 2025;74(Supplement 1):920-P.
   doi:10.2337/db25-920-P
2. Wilson C. Chris Wilson on Insulin Sensitivity (Correction Factor) with Loop and
   Learn, 2/6/2022. Loop and Learn, YouTube.
   https://www.youtube.com/watch?v=oL49FhOts3c. The survey is described from 10 to
   16 minutes and the equation is given at 28 minutes.
3. AndroidAPS source history, github.com/nightscout/AndroidAPS. Commit f98175755e
   (12 February 2022) introduced 277700 / (TDD * bg), and commit c233f9ba12 (13
   June 2022) replaced it with 1800 / (TDD * ln(bg / divisor + 1)).
