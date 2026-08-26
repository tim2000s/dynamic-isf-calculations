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
constants are close to right for adults and substantially wrong for children**,
and the pooled figures everyone quotes, the May analysis included, are driven by
which ages happen to be in the sample.

## The three constants

| Constant | Walsh | All 719 | 95% interval | Adults only | Under 18 |
|---|---|---|---|---|---|
| Sensitivity x daily dose | 1700 | 2108 | 2048 to 2183 | **1799** | 2390 |
| Carb ratio x daily dose | 500 | 409 | 392 to 422 | 408 | 414 |
| Basal share of daily dose | 0.50 | **0.49** | 0.48 to 0.50 | 0.54 | 0.47 |

Read across the sensitivity row rather than down it. Pooled, the constant sits 24%
above Walsh and the interval excludes it. Restricted to adults it lands at 1799,
and for people over 45 at 1766 with an interval of 1680 to 1839, which is the one
cell in this entire analysis consistent with 1700.

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

## What to do with this

For an adult, 1700 divided by daily dose remains a defensible starting sensitivity
factor and the data supports it more than the May analysis suggested. For a child
it is too aggressive by roughly a third, and the gap widens through the school-age
years.

The 500 rule should be about 400 for adults. The half-of-daily-dose basal
convention can stay as it is.

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

Code is `inv009/walsh_constants.py`; results are `results/inv009_walsh_by_cohort.parquet`
and `results/inv009_walsh_constants.json`.
