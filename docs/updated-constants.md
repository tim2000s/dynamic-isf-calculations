---
title: "Updated constants from data analysis"
subtitle: "The three Walsh rules recomputed across 719 people in seven cohorts"
author: "Tim Street"
date: "26 August 2026"
---

## What this replaces

An analysis in May tested Walsh's three rules against 138 open-source loop users
and found all three wrong, including the parametric form: the slope of log
sensitivity on log daily dose came out at -0.43 where the rule implies -1.

This repeats that work on 719 people. The original 138 are included unchanged, so
the earlier result can be read beside the new ones rather than replaced by them,
and it reproduces exactly.

The main finding is one the smaller cohort could not have shown. **Walsh's
sensitivity constant is nearer right for adults than for children and wrong for
both**, and the pooled figures everyone quotes, the May analysis included, are
driven by which ages happen to be in the sample. In adults it is out by 13% and
in children by 41%, in the same direction; the carb ratio rule misses by a fifth
at every age; only the basal rule holds.

## The three constants

| Constant | Walsh | All 719 | 95% interval | Adults only | Under 18 |
|---|---|---|---|---|---|
| Sensitivity x daily dose | 1700 | 2108 | 2048 to 2183 | **1915** | 2390 |
| Carb ratio x daily dose | 500 | 409 | 392 to 422 | 404 | 414 |
| Basal share of daily dose | 0.50 | **0.49** | 0.48 to 0.50 | 0.52 | 0.47 |

Read across the sensitivity row rather than down it. Both the pooled and the adult
columns cover all 719 people; the under-18 column is the JAEB studies alone, because
the OpenAPS Commons cohort contains no children. Pooled, the constant sits 24% above
Walsh and the interval excludes it. Restricted to adults it lands at 1915, and within
that the two sources disagree: 1799 in the JAEB studies against 2381 in Commons, whose
members set their own numbers rather than being given them. For JAEB people over 45 it
falls to 1766 with an interval of 1680 to 1839, which is the one cell in this entire
analysis consistent with 1700.

Half of the pooled sample is under 18, so the headline 2108 is close to an average of
two different populations rather than a description of either. That is the reason to
read the row and not the corner, and it is the trap the May analysis fell into.

The carb ratio rule is the clearest failure and the most consistent. At 409 it
misses 500 by a fifth in every cohort and at every age, so people enter more
insulin per gram than the rule prescribes.

The basal rule survives. At 0.49 across everyone the interval covers Walsh's half
exactly. That is worth stating plainly, because it is the only one of the three
that holds, and it holds despite wide variation underneath it.

## Sensitivity by age, which is where the pooled figure comes from

| Age band | n | Sensitivity x dose | Interval | Consistent with 1700 |
|---|---|---|---|---|
| 2 to 5 | 115 | 2163 | 2045 to 2389 | no |
| 6 to 12 | 128 | 2561 | 2406 to 2695 | no |
| 13 to 17 | 54 | 2350 | 2243 to 2624 | no |
| 18 to 44 | 162 | 1830 | 1738 to 1944 | no |
| 45 and over | 122 | 1766 | 1680 to 1839 | **yes** |

The constant falls steadily from age six onward and only reaches Walsh's value in
the oldest group. A child of ten needs roughly half again the sensitivity factor
the rule would give them, which in dosing terms means the rule would correct a
child about 50% harder than their own settings say.

The carb ratio moves the other way and more weakly: 297 at ages 2 to 5, rising to
479 by 6 to 12, then settling near 400 in adults. Only the school-age bands are
consistent with the 500 rule.

## The parametric form, and a caution about how it is usually tested

Walsh's rule implies that sensitivity is inversely proportional to daily dose, so
a log-log slope of exactly -1.

| Cohort | n | Slope | Interval | Excludes -1 |
|---|---|---|---|---|
| All JAEB pooled | 581 | -0.98 | -1.03 to -0.93 | no |
| JAEB adults only | 284 | -0.73 | -0.82 to -0.65 | yes |
| JAEB under 18 | 297 | -0.91 | -0.96 to -0.85 | yes |
| OpenAPS Commons | 138 | -0.43 | -0.59 to -0.27 | yes |
| Everything | 719 | -0.87 | -0.93 to -0.82 | yes |

The pooled JAEB slope of -0.98 looks like a confirmation of Walsh and is not one.
Restricting the same people to adults moves it to -0.73. What produces the tidy
-1 is mixing young children, who use little insulin and have high sensitivity
factors, with adults who do neither. Between-group differences in level create a
gradient that has nothing to do with how sensitivity varies within anybody.

This is not a range effect. The v7 and Loop cohorts span an identical spread of
daily doses, 1.43 in logs, and return slopes of -0.46 and -0.94.

So a validation of the 1700 rule performed on a mixed-age cohort will tend to
confirm it, and the confirmation will be an artefact of composition. That is the
methodological point of this analysis and it applies to the May result as much as
to anyone else's.

## Results by cohort

| Cohort | n | Sensitivity x dose | Carb ratio x dose | Basal share | Slope |
|---|---|---|---|---|---|
| Loop | 194 | 2043 (1913 to 2249) | 375 (357 to 407) | 0.58 | -0.94 |
| REPLACE-BG | 192 | 1853 (1791 to 1936) | 431 (410 to 460) | 0.51 | -0.72 |
| DCLP5 | 100 | 2429 (2306 to 2580) | 513 (479 to 542) | 0.46 | -1.00 |
| PEDAP | 95 | 2162 (2045 to 2317) | 291 (266 to 317) | 0.41 | -0.89 |
| v5 Trio | 22 | 2682 (2216 to 2962) | 373 (279 to 459) | 0.45 | -0.53 |
| v6 AAPS classic | 19 | 2953 (2201 to 5713) | 526 (409 to 607) | 0.39 | -0.03 |
| v7 oref0 | 97 | 2164 (2044 to 2566) | 387 (315 to 436) | 0.48 | -0.46 |

DCLP3 and IOBP2 are absent because neither release ships a settings file. The
bionic pancreas has no sensitivity factor to record by design.

Three cohort-level observations.

The paediatric cohorts differ from each other more than from the adults. PEDAP,
ages 2 to 5, has the lowest carb ratio constant at 291 and the lowest basal share
at 0.41. DCLP5, ages 6 to 13, has the highest sensitivity constant among the
paediatric groups at 2429 and a carb ratio constant of 513 that is consistent
with Walsh.

Loop's basal share of 0.58 is the highest here and the furthest from Walsh. It is
also a delivered share rather than a programmed one, and under an automated
system those differ, so it is not directly comparable to the OpenAPS Commons
figures, which are programmed.

The three OpenAPS Commons cohorts all have shallow slopes, from -0.03 to -0.53,
while the JAEB cohorts sit between -0.72 and -1.00. These are people who have
tuned settings themselves over years against people whose settings came from a
clinic within a trial. A rule fits best where people follow it.

## Measured rather than entered

Everything above multiplies a setting somebody typed into a pump by their daily
dose, which answers what people run. It does not answer what the constant should
be. Two of the three rules can be measured from behaviour instead.

### Carb ratio: 415, and the measurement is trustworthy

For a meal that was announced, dosed for, and left glucose where it started, the
ratio that person needed was the grams divided by the units. No insulin model, no
basal assumption, no settings file.

Of 78,705 announced meals from 1,011 people, 23,596 ended within 20 mg/dL of their
starting glucose. Restricting to people with at least eight such meals gives 742.

| Cohort | n | Carb ratio x dose | 95% interval | Covers 500 |
|---|---|---|---|---|
| Loop | 594 | 418 | 408 to 429 | no |
| REPLACE-BG | 148 | 390 | 362 to 427 | no |
| Everyone | 742 | **415** | 405 to 425 | no |

The estimator passes its own honesty test. Across the 289 people who also have an
entered carb ratio, entered divided by measured is **1.02**. So this measurement
carries no meaningful bias, and 415 is what the rule should be rather than what
people happen to have typed. It agrees with the entered figure of 409 because in
this case people have already found the right answer.

### Basal share: 0.51, but it belongs to the system rather than the person

Delivered basal over delivered total needs no model at all, so it covers the two
cohorts with no settings file. That is 438 people every table above omits.

| Cohort | System | n | Basal share | Covers 0.50 |
|---|---|---|---|---|
| Loop | DIY loop | 830 | 0.60 | no |
| REPLACE-BG | none | 196 | 0.51 | yes |
| DCLP3 | Control-IQ | 112 | 0.48 | yes |
| DCLP5 | Control-IQ | 100 | 0.46 | no |
| PEDAP | Control-IQ | 95 | 0.41 | no |
| IOBP2 | bionic pancreas | 326 | 0.39 | no |
| Everyone | mixed | 1,659 | 0.51 | no |

Restricting to adults removes age as the explanation and the pattern sharpens.
REPLACE-BG, where no algorithm is running, sits at 0.51 and Control-IQ at 0.50,
both on Walsh. The bionic pancreas delivers 0.40 and the do-it-yourself loop 0.63.

So Walsh's half describes what people programme when nothing is adjusting it, and
survives one automated system unchanged. The other two move it substantially, in
opposite directions, and that is a property of the algorithm rather than of the
person using it.

### Sensitivity: the level cannot be measured here

The third rule resists the same treatment, and INV-009 established why before this
analysis began. Fitting the overnight fall against insulin recovers the shape of a
person's sensitivity and not its level, because insulin action is reconstructed
rather than observed, residual carbohydrate rides along with meal boluses, and a
loop delivers most insulin when glucose is refusing to move.

The size of that is now measured. Against 394 people with both, the regression
reads 15% of a working sensitivity. A matched-correction estimator needing neither
settings nor a basal model disagrees, reading between 50% in the open-loop cohort
and 19% in the youngest, and an attenuation that tracks how reactive each
controller is rules out a single correction factor.

One thing survives, because a common multiplier cannot change it. The measured
constant does not drift with daily dose: slope +0.106, interval -0.026 to +0.241.
A single sensitivity constant does fit across the dose range, which is the part of
Walsh's claim this data can settle. What that constant equals, it cannot.

Establishing the level would need a fixed correction, no carbohydrate, no
algorithm intervening, and glucose watched afterwards. REPLACE-BG is the only
cohort here with no algorithm.

## What to do with this

For an adult, the measured figure is 1915 divided by daily dose and the interval,
1828 to 1986, excludes 1700. Walsh survives only in the JAEB over-45s at 1766. So
1700 is a reasonable floor rather than the centre: it will read slightly strong for
most adults, which is the safer direction to be wrong in for a starting number, but
it should not be defended as the measured value. For a child it is too aggressive by
roughly a third, and the gap widens through the school-age years.

The 500 rule should be about 415. That comes from measurement rather than from
settings and carries no detectable bias, so it is the firmest number in this
document.

The half-of-daily-dose basal convention can stay for someone with no automation
and for Control-IQ. It does not describe what a do-it-yourself loop delivers, at
0.63 in adults, or a bionic pancreas at 0.40, and neither of those is a setting
anybody chose.

None of these are settings to hold to. INV-009 measured the spread of what people
actually run against the rule at 0.68 to 2.07, so any constant is a starting point
that then needs individual adjustment, and a quarter of people end up beyond what
an automated system's own adjustment range could reach from it.

## Method

Settings are those people had entered, taken as a per-person median: from the
bolus calculator for Loop and REPLACE-BG, and from the case report forms for
DCLP5 and PEDAP. Daily dose is measured from each person's own pump record over
days the pump was recording completely, requiring at least thirty. The OpenAPS
Commons figures are the May extraction unchanged.

Intervals are percentile bootstrap over people, 4,000 resamples. Slopes are least
squares on the logarithms with 1,500 bootstrap resamples.

Basal share is the programmed share where a study records a schedule and the
delivered share otherwise, and the two are labelled separately in the machine
readable output because they are not the same quantity.

Measured carb ratio takes announced meals of at least 20 g with a dose within 20
minutes, no other carbohydrate for four hours before or five hours after, and no
further insulin during, keeping those that ended within 20 mg/dL of their starting
glucose. Meals that did not end neutral are excluded rather than corrected.

Code is `inv009/walsh_constants.py` for the entered figures and
`inv009/measured_basal_cr.py` for the measured ones; `inv009/measured_constant.py`
holds the sensitivity attempt and documents why its level is not reported.
Results are in `results/inv009_walsh_by_cohort.parquet`,
`results/inv009_measured_basal_cr.json` and their companions.
