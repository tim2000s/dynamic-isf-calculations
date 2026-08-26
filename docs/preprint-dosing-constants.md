---
title: "Reassessing the Walsh insulin dosing constants in type 1 diabetes"
subtitle: "Age composition in entered settings, and the limits of measurement from device data"
author: "Tim Street, Diabettech"
date: "26 August 2026"
abstract: |
  Insulin pump therapy is initialised from three rules that divide a constant by
  total daily dose: 1700 for the sensitivity factor, 500 for the carbohydrate
  ratio, and half the daily dose for basal. The rules predate continuous glucose
  monitoring and have been validated mostly in small cohorts of mixed age. This
  analysis recomputed them in 719 people with entered pump settings and 1,678
  people with device data, drawn from six Jaeb Center trials and the OpenAPS Data
  Commons, spanning ages 2 to 82 and daily doses of 7 to 107 U. Pooled across all
  ages the sensitivity constant was 2108 (95% CI 2048 to 2183), which excludes
  1700. Restricted to adults it was 1915 (1828 to 1986) and in people under 18 it
  was 2390 (2281 to 2526). The pooled log-log slope of sensitivity on daily dose
  was -0.98 (-1.03 to -0.93), apparently confirming the inverse rule, but the same
  people restricted to adults gave -0.73 (-0.82 to -0.65). Because 51% of the
  pooled sample was under 18, the apparent confirmation was an artefact of age
  composition. The carbohydrate ratio constant measured from 23,596 announced
  meals that ended within 20 mg/dL of their starting glucose was 415 (405 to 425),
  and entered divided by measured was 1.02, indicating no detectable bias. Basal
  share was 0.51 overall and depended on the automation in use rather than on the
  person, at 0.51 with no algorithm, 0.50 under one commercial system, 0.40 under
  a bionic pancreas and 0.63 under a do-it-yourself loop. The sensitivity constant
  could not be measured in the same way. In the single open-loop cohort a
  correction of 1 to 3.5 U produced 15 to 20 mg/dL of glucose fall per acting unit
  against 45 entered, and in every closed-loop cohort the estimate was near zero
  or negative because comparison nights had already been corrected through basal.
  Entered and measured sensitivity are different quantities, since a pump setting
  encodes the whole expected fall including spontaneous reversion of 21 to 30
  mg/dL over four hours. Measured directly from isolated overnight corrections in
  the open-loop cohort, with a six-hour window, the constant was 1483 (1291 to
  1654) in 136 people, an interval covering Walsh's 1700 and lying below the 1915
  people enter. Measured in the paediatric cohorts it was around 2744, so the
  age effect is confirmed independently of settings and is larger in measurement,
  at a ratio near 1.6 against 1.25 entered.
---

# Introduction

Initial insulin pump settings are conventionally derived from total daily dose
using constants that long predate continuous glucose monitoring. The sensitivity
factor, the fall in glucose expected from one unit of rapid-acting insulin, is
taken as 1700 divided by the daily dose in the formulation attributed to Walsh,
or 1800 in the variant attributed to Davidson [1,2]. The carbohydrate ratio is
taken as 500 divided by the daily dose, and basal insulin is taken as half of it.

These rules remain in wide use. They are the starting point in pump initiation
protocols, they are embedded in the onboarding of automated insulin delivery
systems, and the sensitivity factor in particular has acquired a second role: in
several automated systems it is not only a starting value but a live parameter
that the algorithm scales when it adjusts dosing.

Two features of the evidence behind them are worth stating. The cohorts used for
validation have generally been small and of mixed age, and the quantity validated
has generally been the setting people had entered rather than the response
measured from their glucose. Both matter. If sensitivity factors differ
systematically between children and adults, then a constant estimated by pooling
them describes neither group, and the apparent fit of an inverse-proportional
form can be produced by between-group differences in level rather than by any
within-person relationship. Separately, an entered setting records what somebody
was advised or has settled on, which need not be what a unit of insulin does.

Public release of several large trial datasets, together with the OpenAPS Data
Commons, makes it possible to address both points in the same analysis. This
paper recomputes the three constants in 719 people with entered settings and 1,678
people with device data, tests whether the inverse-proportional form survives
stratification by age, and asks how far each constant can be measured from
behaviour rather than read from a settings file.

An earlier analysis by the author, conducted in May 2026 on 138 open-source loop
users, found all three rules inconsistent with the data and reported a log-log
slope of -0.43 where the rule implies -1. That cohort is included here unchanged,
and one purpose of the present work is to establish whether its conclusions were
a property of that population.

# Methods

## Data sources

Device data came from six trials released by the Jaeb Center for Health Research
through its public diabetes dataset portal, and from three extractions of the
OpenAPS Data Commons.

| Cohort | System | People | CGM readings | Median days |
|---|---|---|---|---|
| Loop | DIY closed loop | 851 | 89,203,827 | 416 |
| REPLACE-BG | Pump and CGM, no algorithm | 196 | 12,034,638 | 238 |
| IOBP2 | Bionic pancreas | 343 | 7,891,499 | 102 |
| PEDAP | Control-IQ, ages 2–5 | 99 | 6,622,228 | 298 |
| DCLP3 | Control-IQ, adolescent and adult | 112 | 5,673,867 | 183 |
| DCLP5 | Control-IQ, ages 6–13 | 100 | 5,535,124 | 218 |
| v5 Trio | DIY closed loop | 22 | not applicable | not applicable |
| v6 AAPS | DIY closed loop | 19 | not applicable | not applicable |
| v7 oref0 | DIY closed loop | 97 | not applicable | not applicable |

The corpus holds 132,467,585 CGM readings, 68,995,020 basal records, 5,154,053
bolus records and 1,617,910 carbohydrate entries. The three OpenAPS Data Commons
extractions carry settings and summary statistics rather than the raw streams, so
they contribute to the entered-settings analysis alone.

REPLACE-BG occupies a distinct position in this corpus and carries much of the
argument below. It ran in 2015, before automated delivery was available to its
participants, so insulin was delivered by a pump following a schedule the person
had programmed, with corrections given by hand. It is the only cohort here in
which a night without a correction is a night on which nothing responded to
glucose.

Dates are de-identified in every Jaeb release, by a per-participant random shift
in some studies and by rebasing to a common epoch in others. Time of day survives
both transformations, which is what this analysis requires. Loop timestamps carry
a fixed per-participant offset from UTC, so daylight saving places time of day one
hour out for part of each year; the effect on window start hours between 23:00 and
03:00 is tolerable and is noted as a limitation.

## Entered settings

Sensitivity factors and carbohydrate ratios were taken as the per-person median
of the values recorded at the time of dosing. For Loop and REPLACE-BG these came
from the bolus calculator record, and for DCLP5 and PEDAP from the pump settings
captured on case report forms, which store up to ten daily segments each with a
start and end time; segment values were combined by time weighting. DCLP3 and
IOBP2 ship no settings file, and the bionic pancreas has no sensitivity factor to
record by design, so both are absent from the entered analysis.

Sensitivity factors in the Loop and REPLACE-BG calculator records are stored in
mmol/L per unit for almost all participants. A per-person rule converted values
whose median fell below 20 by a factor of 18.018, which reclassified 112,776 of
112,780 Loop rows. Conversion was verified by checking that the resulting
distribution fell within 10 to 400 mg/dL per unit.

The settings loader was validated against an independent stream: of 80,039 Loop
calculator rows carrying a carbohydrate entry, 100.00% matched a record in the
separately loaded carbohydrate table at the same timestamp and amount.

## Daily dose

Total daily dose was measured from each person's own pump record rather than
taken from a case report form, summing delivered basal and bolus insulin over
calendar days on which the pump recorded continuously, and requiring at least 30
such days. Measured daily dose correlated with the clinician-recorded value at
Spearman 0.92 in DCLP5 and 0.89 in PEDAP.

One field required correction. A column labelled as a seven-day basal total
proved on inspection to be a daily total, which was established by regressing it
on the measured daily basal delivery.

## Overnight windows

The sensitivity analyses used four-hour windows starting at 23:00, 00:00, 01:00,
02:00 and 03:00 local time. A window was retained when starting glucose lay
between 90 and 300 mg/dL, no carbohydrate had been recorded for six hours before
it, no bolus for three hours before it, glucose had not risen by more than 15
mg/dL in the preceding 30 minutes, and CGM coverage was complete. Screening
conditions were evaluated only on information available before the window opened,
so that no window was selected on the basis of its own outcome. In cohorts with no
carbohydrate stream, fasting was inferred from the absence of any bolus large
enough to be a meal dose, defined as 5% of that person's daily dose, which keeps
the criterion free of scale.

Of 2,179,131 candidate windows, 708,106 passed screening, from 1,678 people.

## Insulin action

Insulin action was computed by convolving delivery with the remaining-action
curve of the appropriate model. For the Loop cohort the model is LoopKit's
exponential form, with the peak and duration of the preset each participant was
using, recovered by comparing recomputed insulin on board against the value the
app recorded at each bolus; agreement reached r = 0.927. Elsewhere the oref
exponential form with a six-hour duration and 75-minute peak was used.

Action within a window was separated into two components by convolution. The
first is action arising from insulin delivered before the window opened, which is
predetermined with respect to anything that happens inside it. The second is
action from insulin delivered during the window, which is not.

Two unit errors in the underlying data required correction before any of this was
reliable. Bolus durations were recorded in milliseconds rather than seconds for
three studies, which placed insulin in 41,591 of one participant's 46,176 time
bins; durations above 86,400 were rescaled and those implying more than eight
hours were treated as instantaneous. Separately, timestamps returned at microsecond
resolution were being cast as though they were nanoseconds, dividing every basal
total by 1000 while leaving boluses correct. Both were caught by a mass
conservation check comparing grid totals against database sums, which now agrees
to 0.0000%.

## Exposure and matching

The exposure for the sensitivity analysis is insulin delivered above or below the
person's own scheduled basal profile, by whichever route it arrived. An algorithm
running temporary basal above the schedule has given a correction, and one
suspending delivery has given a negative correction. Defining the exposure as
boluses alone would discard most of the dosing variation in every automated
cohort. Where a study records a programmed schedule this was used directly; where
it does not, the reference is that person's own rolling median delivery for each
half hour of the day over the preceding 30 days.

Two estimators were applied. The first is a matched comparison: for each isolated
overnight correction, comparison windows were drawn from the same person at the
same hour with starting glucose within 15 mg/dL and no bolus, and sensitivity was
taken as the difference in four-hour fall divided by the difference in acting
insulin. Estimates were pooled as a ratio of sums across events rather than as a
median of per-event ratios, because dividing each event by its own small and noisy
denominator both inflates the tail and biases the centre toward zero.

The second is a stratified regression. Windows were grouped by person, hour,
starting glucose to 15 mg/dL, slope into the window to 1 mg/dL per 5 min, and
glucose one hour earlier to 25 mg/dL. Within each stratum the fall and the
exposure were centred, and the pooled within-stratum slope was taken as the
estimate. Strata with fewer than four windows or an exposure range below 0.4 U
were dropped.

## Measured carbohydrate ratio

For a meal that was announced and dosed for, and that left glucose where it
started, the ratio that person required is the grams divided by the units, with no
insulin model or basal assumption involved. Meals of at least 20 g with a dose within 20
minutes were selected, requiring no other carbohydrate for four hours before or
five hours after and no further insulin during. Of 78,705 such meals from 1,011
people, 23,596 ended within 20 mg/dL of their starting glucose. Meals that did not
end neutral were excluded rather than corrected, and people with fewer than eight
qualifying meals were dropped, leaving 742.

## Statistics

Intervals are percentile bootstrap over people, with 4,000 resamples for
constants and 1,500 for slopes, resampling people rather than observations
because windows from one person are not independent. Slopes are least squares on
logarithms. Where per-person estimates were pooled, DerSimonian and Laird random
effects pooling was used and heterogeneity is reported. Analysis ran locally in
Python against a PostgreSQL 17.9 instance with TimescaleDB 2.26.1.

# Results

## Entered constants

| Constant | Walsh | All 719 | 95% CI | Adults (422) | Under 18 (297) |
|---|---|---|---|---|---|
| Sensitivity x daily dose | 1700 | 2108 | 2048–2183 | 1915 | 2390 |
| Carb ratio x daily dose | 500 | 409 | 392–422 | 404 | 414 |
| Basal share of daily dose | 0.50 | 0.49 | 0.48–0.50 | 0.52 | 0.47 |

The pooled sensitivity constant sits 24% above Walsh's value and its interval
excludes it. Restricted to adults it falls to 1915 (1828 to 1986), which also
excludes 1700. Within the adult group the two sources disagree substantially, at
1799 in the Jaeb trials against 2381 in the OpenAPS Data Commons, whose members
set their own values over years rather than receiving them from a clinic. The one
subgroup whose interval covers 1700 is Jaeb participants aged 45 and over, at 1766
(1680 to 1839).

The carbohydrate ratio constant misses 500 by approximately a fifth in every
cohort and at every age. The basal rule holds, at 0.49 with an interval covering
0.50.

## Age composition

| Age band | n | Sensitivity x dose | 95% CI | Covers 1700 |
|---|---|---|---|---|
| 2 to 5 | 115 | 2163 | 2045–2389 | no |
| 6 to 12 | 128 | 2561 | 2406–2695 | no |
| 13 to 17 | 54 | 2350 | 2243–2624 | no |
| 18 to 44 | 162 | 1830 | 1738–1944 | no |
| 45 and over | 122 | 1766 | 1680–1839 | yes |

The constant falls from age six onward and reaches Walsh's value only in the
oldest group. A child of ten carries a sensitivity constant approximately 45%
above an adult of 45, which in dosing terms means the rule would correct that
child considerably harder than their own settings prescribe.

The carbohydrate ratio constant moves in the opposite direction and more weakly,
from 297 at ages 2 to 5, to 479 at 6 to 12, settling near 400 in adults. Only the
school-age bands are consistent with 500.

## The parametric form

Walsh's rule implies that sensitivity is inversely proportional to daily dose,
which is a log-log slope of exactly -1.

| Cohort | n | Slope | 95% CI | Excludes -1 |
|---|---|---|---|---|
| All Jaeb pooled | 581 | -0.98 | -1.03 to -0.93 | no |
| Jaeb adults only | 284 | -0.73 | -0.82 to -0.65 | yes |
| Jaeb under 18 | 297 | -0.91 | -0.96 to -0.85 | yes |
| OpenAPS Data Commons | 138 | -0.43 | -0.59 to -0.27 | yes |
| Everything | 719 | -0.87 | -0.93 to -0.82 | yes |

The pooled Jaeb slope of -0.98 appears to confirm the inverse rule. Restricting
the same people to adults moves it to -0.73, and the interval then excludes -1.
What produces the value near -1 is the mixing of young children, who use little
insulin and carry high sensitivity factors, with adults who do neither.
Between-group differences in level generate a gradient across the pooled sample
that need not reflect how sensitivity varies within any individual.

This is not attributable to the range of daily doses sampled. The v7 and Loop
cohorts span an identical spread of 1.43 in logs and return slopes of -0.46 and
-0.94 respectively.

A validation of the 1700 rule performed on a cohort of mixed age will therefore
tend to confirm it, and the confirmation will reflect composition. This applies to
the author's May 2026 analysis as much as to any other.

## Results by cohort

| Cohort | n | Sensitivity x dose | Carb ratio x dose | Basal share | Slope |
|---|---|---|---|---|---|
| Loop | 194 | 2043 (1913–2249) | 375 (357–407) | 0.58 | -0.94 |
| REPLACE-BG | 192 | 1853 (1791–1936) | 431 (410–460) | 0.51 | -0.72 |
| DCLP5 | 100 | 2429 (2306–2580) | 513 (479–542) | 0.46 | -1.00 |
| PEDAP | 95 | 2162 (2045–2317) | 291 (266–317) | 0.41 | -0.89 |
| v5 Trio | 22 | 2682 (2216–2962) | 373 (279–459) | 0.45 | -0.53 |
| v6 AAPS | 19 | 2953 (2201–5713) | 526 (409–607) | 0.39 | -0.03 |
| v7 oref0 | 97 | 2164 (2044–2566) | 387 (315–436) | 0.48 | -0.46 |

The paediatric cohorts differ from one another more than they differ from the
adult cohorts. PEDAP, covering ages 2 to 5, carries both the lowest carbohydrate
ratio constant at 291 and the lowest basal share at 0.41. DCLP5, covering ages 6
to 13, carries the highest paediatric sensitivity constant at 2429 and a
carbohydrate ratio constant of 513 whose interval covers 500.

The three OpenAPS Data Commons cohorts return shallow slopes, between -0.03 and
-0.53, against -0.72 to -1.00 in the Jaeb cohorts. These are people who have
adjusted their own settings over years, against people whose settings were set
within a trial protocol. The rule fits best where it is followed.

## Measured carbohydrate ratio

| Cohort | n | Carb ratio x dose | 95% CI | Covers 500 |
|---|---|---|---|---|
| Loop | 594 | 418 | 408–429 | no |
| REPLACE-BG | 148 | 390 | 362–427 | no |
| Everyone | 742 | 415 | 405–425 | no |

Among the 289 people who also had an entered carbohydrate ratio, entered divided
by measured was 1.02. The estimator therefore carries no detectable bias, and 415
can be read as the constant the rule should use rather than as a description of
what people have entered. It agrees closely with the entered figure of 409.

## Basal share

Delivered basal over delivered total requires no model, so it covers the two
cohorts with no settings file, adding 438 people omitted from the tables above.

| Cohort | System | n | Basal share | Covers 0.50 |
|---|---|---|---|---|
| Loop | DIY closed loop | 830 | 0.60 | no |
| REPLACE-BG | None | 196 | 0.51 | yes |
| DCLP3 | Control-IQ | 112 | 0.48 | yes |
| DCLP5 | Control-IQ | 100 | 0.46 | no |
| PEDAP | Control-IQ | 95 | 0.41 | no |
| IOBP2 | Bionic pancreas | 326 | 0.39 | no |
| Everyone | Mixed | 1,659 | 0.51 | no |

Restricting to adults removes age as an explanation and sharpens the pattern.
REPLACE-BG, where no algorithm was running, sits at 0.51 and Control-IQ at 0.50,
both consistent with Walsh. The bionic pancreas delivers 0.40 and the
do-it-yourself loop 0.63. Walsh's half therefore describes what people programme
when nothing is adjusting it, and survives one commercial system unchanged. The
remaining two move it substantially and in opposite directions, which is a
property of the algorithm rather than of the person using it.

## Measured sensitivity

| Cohort | Comparison | Events | Matched fall, 1–3.5 U | Per acting unit | Entered |
|---|---|---|---|---|---|
| REPLACE-BG | Open loop | 239 | 22 to 35 mg/dL | 15 to 20 | 45 |
| Loop | DIY closed loop | 269 | -4 to 5 | -2 to 5 | 57 |
| DCLP3 | Control-IQ | 79 | -8 to -4 | -5 to -4 | not recorded |
| DCLP5 | Control-IQ | 38 | -7 | -9 | 68 |
| PEDAP | Control-IQ | 50 | -6 | -11 | 161 |

Only the open-loop cohort returns a positive estimate. The explanation lies in
the comparison nights rather than in the treated ones. The four-hour fall on
nights when no correction was given, by starting glucose, was as follows.

| Cohort | 120–150 | 150–180 | 180–220 | 220–260 | 260+ |
|---|---|---|---|---|---|
| REPLACE-BG | 12 | 21 | 30 | 35 | 43 |
| Loop | 20 | 40 | 62 | 85 | 104 |
| DCLP5 | 11 | 37 | 68 | 102 | 126 |
| PEDAP | 11 | 38 | 72 | 111 | 146 |

Under an algorithm a night without a bolus is not an untreated night. Glucose
fell two to three times as far as in the open-loop cohort, because the algorithm
delivered the correction through basal and arrived at a comparable endpoint. An
added bolus therefore has little left to produce and little to measure.

This is not a consequence of the matching being too loose, and tightening it does
not help. Progressive conditioning in the Loop cohort moved the estimate from 11.6
mg/dL per unit with person and hour alone, to 4.3 adding glucose to 50 mg/dL, to
2.8 at 15 mg/dL, and then to a plateau: 2.6 adding slope, 2.9 adding glucose one
hour earlier, 2.6 tightening that further. A plateau argues against measurement
error in the reconstructed exposure, which would continue to drive the slope
toward zero as conditioning removed genuine variation. The large first step is
consistent with removal of mean reversion, since glucose that is high both
attracts insulin and falls without it.

What remains is that the algorithm's dose is tied to need the CGM trace does not
carry. Within an identical stratum, the correlation between insulin committed
before a window and insulin delivered during it was 0.054 in REPLACE-BG and 0.447
in Loop, with intermediate values of 0.062 in DCLP3, 0.140 in DCLP5, 0.151 in
PEDAP and 0.212 in IOBP2. Someone who gives a correction and sleeps generates
almost none of this association. A controller sampling every five minutes
continues to respond to information the trace does not contain, so conditioning on
the trace cannot break the association between its dosing and the outcome.

In the one cohort where the question can be posed, a unit produced 15 to 20 mg/dL
of fall overnight against 45 entered, which as a constant is 630 to 840 against
1886. The ratio of measured to entered lay between 0.39 and 0.50 across four
cohorts.

## The constant fitted from observed response

A settings file claims something narrow: that a correction of U units will lower
glucose by ISF times U. That claim can be measured directly. Isolated overnight
corrections of at least 1 U were selected at a starting glucose between 150 and
300 mg/dL, with a six-hour window so that the response of a six-hour insulin model
falls inside the observation, and the fall per unit given was taken as the median
across each person's nights and multiplied by their measured daily dose.

| Cohort | Ages | People | Nights | Fall per unit, mg/dL | Constant | 95% CI |
|---|---|---|---|---|---|---|
| REPLACE-BG, open loop | adult | 136 | 500 | 34.6 | 1483 | 1291 to 1654 |
| DCLP3, Control-IQ | 14+ | 52 | 79 | 36.8 | 1715 | 1407 to 2078 |
| Loop, DIY closed loop | mixed | 261 | 565 | 46.0 | 1896 | 1751 to 2096 |
| DCLP5, Control-IQ | 6 to 13 | 30 | 54 | 62.1 | 2749 | 1907 to 3309 |
| PEDAP, Control-IQ | 2 to 5 | 47 | 93 | 197.8 | 2744 | 2157 to 2962 |

The dose threshold is 2.5% of each person's own daily dose rather than a flat
number of units, which is 1.0 U at a daily dose of 40 and 0.34 U at 13.5. A flat
threshold of 1 U silently excluded the youngest cohorts, keeping 11 of PEDAP's 295
isolated corrections, and those cohorts carry the only measured check on the age
pattern. IOBP2 cannot be measured at all: the bionic pancreas doses through
automatic micro-boluses and there is no user correction to isolate, leaving 2
qualifying nights from 343 people.

REPLACE-BG is the cohort to read for the level. It is the only one in which
nothing but the person's own correction and their programmed basal is acting, so
the fall can be attributed to the dose. Its interval covers Walsh's 1700 and lies
below the 1915 those people and their peers have entered, which says entered
settings are slightly weak: a correction delivers a little more than the setting
predicts.

Every closed-loop cohort is biased upward, because the algorithm continues to act
through the window and part of the fall is the controller rather than the bolus.
The size of that bias can be approximated. DCLP3 and REPLACE-BG are both adult
cohorts, and DCLP3 sits 16% higher at 1715 against 1483, which is the only handle
this data offers on the magnitude. It rests on two cohorts that differ in more
than their automation, so it is an indication rather than a correction factor.

That approximation matters for the age question, because both paediatric cohorts
run Control-IQ and neither has an open-loop counterpart. Taking their figures at
face value gives 2749 and 2744 against an adult 1483, and discounting them by 16%
gives roughly 2370 against 1483. Either way the measured constant is around 1.6
times higher in children than in adults.

That is the same direction as the entered settings and a larger gap. Entered, the
ratio of the under-18 constant to the adult one is 2390 to 1915, or 1.25. Measured,
it is nearer 1.6. So paediatric settings reproduce the direction of the age effect
and understate its size, and children appear to be more sensitive than their own
settings record. This is the only place in this analysis where measurement and
entered settings disagree about something other than level, and it is the strongest
form of the age finding because it does not depend on what anybody typed into a
pump.

The horizon was tested rather than assumed. Rebuilding every window at six hours
instead of four moved the pooled fitted constant from 880 to 874, so truncation of
the insulin tail is not a material source of bias in this design.

## An estimator withdrawn, and a constant with it

Two earlier attempts on this question are withdrawn, and both failed in ways worth
recording because the failures are not obvious from their output.

The first regressed a per-person measured sensitivity on daily dose and reported a
slope of +0.107 with an interval of -0.024 to +0.240 as evidence that one constant
fits the whole dose range. It was fitted only to the 1,182 people of 1,660 whose
measured sensitivity came out above zero, which is selection on the outcome, and
it pooled cohorts whose bias differs. Within-cohort slopes from that estimator run
-0.053 in REPLACE-BG to +0.323 in DCLP5, with the largest cohort excluding zero at
+0.245. The defence offered was that a multiplier common to everyone cannot alter
a slope, and the multiplier is not common: the endogeneity measure above spans
0.054 to 0.447 between cohorts and cohort median daily dose spans 13.6 to 55.5 U,
so a cohort-specific bias varies along the axis the slope is measured on.

The second fitted ISF = K / TDD^b to predict the overnight fall out of sample and
reported K = 880 at Walsh's exponent. The exposure was the action within the
window from insulin given before it opened, chosen because pre-committed insulin
cannot respond to what happens inside the window. That choice discards the
correction. On isolated correction nights only 22 to 29% of the action is
pre-window, and the remainder is the tail of the basal rate, so the fitted
coefficient was the effect of the basal tail rather than of a dose. Refitting on
total action does not repair it, giving 355, because overnight most of that total
is also basal, and basal holds glucose level rather than lowering it. A per-person
intercept does not separate the two, since within a person both exposures vary
mostly through basal.

The direct measurement above avoids this by restricting to nights when the dose is
identifiable, at the cost of using 78 people rather than 1,511.

# Discussion

Four findings follow from this analysis.

The first concerns method, and it arose twice in this work. An
inverse-proportional relationship between sensitivity and daily dose can be
manufactured by pooling groups that differ in level, and in this corpus it was:
the pooled slope of -0.98 became -0.73 in the same people restricted to adults,
with 51% of the pooled sample under 18. Any validation of a dose-derived constant
performed on a cohort of mixed age is therefore weak evidence for the functional
form, and the direction of the artefact is toward confirmation.

The same artefact then appeared in the measured analysis, where the author
initially reported a pooled slope near zero as evidence that one constant fits the
whole dose range. It does not survive stratification. Within-cohort slopes from
that estimator run from -0.05 to +0.32 and the largest cohort excludes zero at
+0.245. The independent predictive fit reaches the same conclusion by a sounder
route, since the best-fitting exponent there spans 0.50 to 1.10 between cohorts
and the constant at a fixed exponent spans a factor of 3.39. Pooling estimates
whose bias differs by group, across groups that differ in the exposure, reproduces
the error in a second place. This applies to the author's earlier work and to an
earlier draft of the present paper.

The second concerns the constants themselves. The carbohydrate ratio rule is the
clearest departure from practice and the most consistent, at 415 measured and 409
entered against 500, with agreement between the two approaches and no detectable
bias in the measurement. The basal rule holds where nothing is adjusting delivery
and does not describe what an automated system delivers, which is unsurprising
once stated but is worth separating from any claim about the person. The
sensitivity rule is 13% out in adults and 41% out in children, in the same
direction, and covers 1700 only in adults aged 45 and over.

The third concerns which of the three original constants survives measurement,
and the answer inverts the picture from entered settings alone. Measured from
isolated corrections where nothing else was intervening, the sensitivity constant
is 1483 with an interval of 1291 to 1654, which covers Walsh's 1700. What people
enter, 1915 in adults, sits above both, so entered settings are slightly weak and
a correction delivers a little more than the setting predicts. In the paediatric
cohorts the measured constant is near 2744, so the age pattern found in entered
settings is reproduced by measurement and is larger there, at a ratio near 1.6
against 1.25. The carbohydrate
ratio goes the other way and further: 415 measured against 500, with entered
settings at 409 agreeing with the measurement rather than with the rule.

So Walsh's sensitivity constant is closer to what insulin does than to what people
have written down, and the carbohydrate rule is wrong by a fifth in both. That
distinction is only visible because both quantities were measured rather than read
off settings, and it is the reason this analysis was worth doing on device data at
all.

A caution attaches to the sensitivity figure that does not attach to the
carbohydrate one. The level rests on 136 people in a single cohort,
because attributing a fall to a dose requires no algorithm to be running, and the
age comparison rests on paediatric cohorts that are all closed loop and therefore
biased upward by an amount estimated from a single adult pair. The carbohydrate
constant rests on 742 people across cohorts. The
intervals reflect this and the sensitivity interval is correspondingly wide.

The fourth concerns study design. The absence of a measurable bolus effect in
four automated cohorts is a finding about those systems rather than a failure of
measurement, since the algorithm reached a comparable glucose endpoint through
basal adjustment alone. Observational estimation of
insulin sensitivity from automated insulin delivery data is obstructed by the
controller acting on information that is not recoverable from the CGM trace, and
the strength of that obstruction can be quantified by the association reported
above. This bears on the growing literature that estimates sensitivity from
closed-loop device data, and suggests such estimates require an open-loop
comparison or an experimental manipulation of dose.

# Limitations

The analysis is observational throughout, and no dose was assigned by the
investigator. Confounding by indication is the dominant threat to the sensitivity
estimates and is the explicit subject of the closed-loop results.

The window horizon is four hours against insulin models of six, so part of the
glucose response falls outside the observation. Several mechanisms could each produce a
measured-to-entered ratio near 0.4. Truncation at four hours is one. Selection into
correcting at times when glucose is expected to remain high is a second, and
optimism in entered settings a third. The present design does not separate them,
and extending the horizon to six hours would address the first. That extension is
planned.

The measured sensitivity estimator returns a value at or below zero for 29% of
people, 478 of 1,660, and any analysis restricted to positive values is selecting
on the outcome. Survival of that filter is close to independent of daily dose, at
Spearman +0.042 with p = 0.084, so it does not by itself generate a spurious
scaling. It does mean the surviving group is not a random sample of the cohort,
and the pooled estimates in the measured section should be read with that in mind.

Insulin action is reconstructed rather than observed. For the Loop cohort the
model was selected per person against the app's own recorded insulin on board,
reaching r = 0.927, but elsewhere a single model was assumed.

Entered settings were available for 719 people of the 1,678 with device data, and
two cohorts ship no settings file. The under-18 entered figures rest on Jaeb
participants alone, since the OpenAPS Data Commons cohort contains no children.

Daylight saving displaces time of day by one hour for part of the year in the Loop
cohort. Body weight is not available in the loaded data, so sensitivity per
kilogram could not be examined, and weight is a plausible mediator of the age
pattern reported here. The OpenAPS Data Commons participants are self-selected and
technically engaged, which is consistent with their distinct slopes and should
temper generalisation from them.

# Conclusions

Recomputed across 719 people with entered settings and 1,678 with device data, the
carbohydrate ratio constant is 415 rather than 500 and the basal share is 0.51
where no algorithm is adjusting delivery. The sensitivity constant is 1915 in
adults and 2390 in people under 18, and pooling those groups produces both a
constant and a functional form that describe neither. The inverse-proportional
form is not supported once age is held constant. Measured directly from isolated
overnight corrections where no algorithm was intervening, the sensitivity constant
is 1483 with an interval of 1291 to 1654, which covers Walsh's 1700 and sits below
the 1915 people enter, and it is around 1.6 times higher in children than adults
against 1.25 in entered settings. Of the three rules, sensitivity is therefore the one whose
original constant survives measurement, and the carbohydrate ratio is the one that
does not.

# Data availability

The Jaeb datasets are available from https://public.jaeb.org/datasets/diabetes
under the terms stated there. The OpenAPS Data Commons is available on
application to its curators. No participant-level data are redistributed with this
work.

# Code availability

Analysis code is at https://github.com/tim2000s/dynamic-isf-calculations. The
entered constants are produced by `inv009/walsh_constants.py`, the measured
carbohydrate ratio and basal share by `inv009/measured_basal_cr.py`, the
per-cohort sensitivity measurement by `inv009/measured_isf_by_study.py`, and the
net-dose estimator and its diagnostics by `inv009/net_dose_isf.py`.

# Acknowledgements and attribution

The source of the data is the Loop Study (sponsored by the Jaeb Center for Health
Research and funded by the Helmsley Charitable Trust), but the analyses, content
and conclusions presented herein are solely the responsibility of the authors and
have not been reviewed or approved by the study sponsor.

The source of the data is the Insulin Only Bionic Pancreas Pivotal Trial, but the
analyses, content and conclusions presented herein are solely the responsibility
of the authors and have not been reviewed or approved by the Bionic Pancreas
Research Group or Beta Bionics.

Data from REPLACE-BG, DCLP3, DCLP5 and PEDAP were supplied by the Jaeb Center for
Health Research. The analyses, content and conclusions presented herein are solely
the responsibility of the author and have not been reviewed or approved by the
study sponsors. The author thanks the participants of all seven cohorts, and the
curators of the OpenAPS Data Commons.

# References

1. Walsh J, Roberts R. Pumping Insulin. 6th ed. San Diego: Torrey Pines Press; 2016.
2. Davidson PC, Hebblewhite HR, Steed RD, Bode BW. Analysis of guidelines for
   basal-bolus insulin dosing: basal insulin, correction factor, and
   carbohydrate-to-insulin ratio. Endocrine Practice. 2008;14(9):1095-1101.
3. Aleppo G, Ruedy KJ, Riddlesworth TD, et al. REPLACE-BG: A Randomized Trial
   Comparing Continuous Glucose Monitoring With and Without Routine Blood Glucose
   Monitoring in Adults With Well-Controlled Type 1 Diabetes. Diabetes Care.
   2017;40(4):538-545.
4. Brown SA, Kovatchev BP, Raghinaru D, et al. Six-Month Randomized, Multicenter
   Trial of Closed-Loop Control in Type 1 Diabetes. New England Journal of
   Medicine. 2019;381(18):1707-1717.
5. Breton MD, Kanapka LG, Beck RW, et al. A Randomized Trial of Closed-Loop
   Control in Children with Type 1 Diabetes. New England Journal of Medicine.
   2020;383(9):836-845.
6. Wadwa RP, Reed ZW, Buckingham BA, et al. Trial of Hybrid Closed-Loop Control in
   Young Children with Type 1 Diabetes. New England Journal of Medicine.
   2023;388(11):991-1001.
7. Bionic Pancreas Research Group. Multicenter, Randomized Trial of a Bionic
   Pancreas in Type 1 Diabetes. New England Journal of Medicine.
   2022;387(13):1161-1172.
8. Lewis D, Leibrand S, #OpenAPS Community. Real-World Use of Open Source
   Artificial Pancreas Systems. Journal of Diabetes Science and Technology.
   2016;10(6):1411.
9. DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical
   Trials. 1986;7(3):177-188.
