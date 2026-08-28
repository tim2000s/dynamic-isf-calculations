---
title: "How far observational insulin sensitivity sits below entered settings depends on how it is measured"
subtitle: "Six trial cohorts, 1,679 people, and a gap that ranges from 0.30 to 1.12 according to construction"
author: "Tim Street, Diabettech"
date: "27 August 2026"
---

## Summary

Sensitivity estimated from what glucose did after insulin was given generally comes
out below the sensitivity people have entered into their pumps. That was
established in 2026 on 138 users of open-source closed-loop systems, where the
constant anchored to measured sensitivity was 145 against 355 anchored to the
profiles people had tuned, a ratio of 0.41. The same work tested whether dosing to
the measured value would be safe and found it would not: people dosing well weaker
than their measured sensitivity had more hypoglycaemia rather than less, and
restricting the estimate to overnight periods left the bias unchanged.

This paper asks whether that ratio holds on other data. Across 1,679 people in six
trial cohorts it does not hold at a single value, and the reason is that the ratio
depends on how the estimate is constructed more than on which cohort is measured.
Holding the cohort fixed and varying the construction moves it from 0.30 to 1.12.
Running the original construction unchanged on the trial cohorts gives 0.72 to
1.12, with the open-loop cohort exceeding its entered settings at six hours.

So the direction of the earlier finding is not universal and the magnitude is not
transferable. What is transferable is the caution: an observational sensitivity
estimate is a statement about a method as much as about a person, and two
defensible methods on the same records differ by more than threefold.

## What was already known

The earlier analysis reached three conclusions relevant here. The shape of the
relationship between sensitivity and total daily dose is derivable from data and
scales approximately as one over the square root of daily dose. The level is not:
a short-window empirical estimate swings by a third between blocks and carries
almost no memory from one block to the next. And the level that comes out is
biased in a specific direction, low, such that dosing to it produces more
hypoglycaemia than dosing to a tuned profile.

That last point is the one that matters for interpretation. The gap between
measured and entered sensitivity is not evidence that entered settings are wrong.
It had already been tested against outcomes and the entered settings won.

## What this adds

The cohorts are the Jaeb Center public releases: REPLACE-BG, which ran in 2015
with no algorithm intervening, the Loop observational study, DCLP3, DCLP5 and
PEDAP on Control-IQ across ages 2 to adult, and IOBP2 on a bionic pancreas. Ages
run from 2 to 82 and daily doses from 7 to 107 U. Insulin on board is
reconstructed from the delivery record and an insulin action curve, which is the
weaker provenance: the earlier work had the loop's own figure.

The original construction is the fall over four hours divided by insulin on board
at the start, taken overnight, with carbohydrate screened by the shape of the trace
rather than from a log. Run unchanged on these cohorts it gives:

| Cohort | 4 h | 6 h | Entered | Ratio at 4 h | Ratio at 6 h |
|---|---|---|---|---|---|
| REPLACE-BG | 42.5 | 50.4 | 45.0 | 0.94 | 1.12 |
| Loop | 45.9 | 50.7 | 55.0 | 0.83 | 0.92 |
| PEDAP | 137.4 | 144.0 | 158.8 | 0.87 | 0.91 |
| DCLP3 | 31.2 | 33.6 | 41.0 | 0.76 | 0.82 |
| DCLP5 | 54.9 | 49.4 | 68.3 | 0.80 | 0.72 |

That is 0.72 to 1.12, not 0.41, and the open-loop cohort exceeds its entered
settings once the horizon covers the insulin duration.

Three other constructions were also run, differing in decisions that turn out to
matter as much as the cohort does. Whether a correction is counted as a bolus
alone or as anything delivered above the programmed basal, which matters because
56% of Loop's correcting insulin and 43% of PEDAP's arrives as raised temporary
basal. Whether the denominator is the units given, the insulin action across the
window, or insulin on board at the start alone. And whether episodes are required
to open with little insulin already on board.

| Construction | Cohorts | Ratio |
|---|---|---|
| Insulin on board at the start, 4 to 6 h | 5 | 0.72 to 1.12 |
| Bolus corrections, per unit given | 5 | 0.49 to 0.92 |
| Any route, insulin action across the window | 5 | 0.48 to 0.66 |
| Any route, action, low insulin on board at the start | 4 | 0.30 to 0.65 |
| oref archive, loop-recorded insulin on board | 1 | 0.41 |

The constructions bracket rather than agree, and they do so in an interpretable
direction. Dividing by insulin on board at the start ignores anything delivered
inside the window, which lowers glucose and enlarges the numerator without
appearing below the line, so that construction reads high and reads highest where
an algorithm is dosing throughout. Dividing by the action of all above-schedule
insulin counts insulin that was covering something else, so it reads low. The
honest statement is that the answer lies between them and this data does not place
it.

## Why the magnitude cannot be pinned here

Three mechanisms produce a low estimate and this design separates none of them.

Glucose falls on its own from a raised level. Overnight and carbohydrate-free,
between 150 and 250 mg/dL, it drops a median 21 to 30 mg/dL over four hours with
no correction given, so any estimator that removes that fall is measuring
something smaller than a settings file encodes, and any estimator that leaves it
in is measuring mean reversion.

Insulin is given when glucose is expected to stay up. Within an identical glucose
trace the correlation between insulin committed before a window and insulin
delivered during it is 0.054 in the open-loop cohort and 0.447 under a
do-it-yourself loop, so a controller keeps responding to something the trace does
not carry and its dose stays tied to the outcome.

And the insulin action curve is assumed rather than observed. Reconstructed
insulin on board agrees with the pump's own record at a median absolute difference
of 0.079 U and a correlation of 0.848 across 188 people in REPLACE-BG, which is
good but not exact.

## A limit on what a simulation can settle

The obvious way to choose between constructions is to run them against records
whose true sensitivity is known. That was done and it does not resolve the
question. Simulated open-loop records with unlogged carbohydrate return 51.6
against a true 50 for the bolus-and-units-given construction, but only about 22
qualifying episodes per person, and the simulated controller cannot be made to
correct through temporary basal at a realistic rate: at any gain that leaves
glucose in a plausible range it never delivers half a unit above schedule in
thirty minutes. So the constructions that matter for the closed-loop cohorts
cannot be tested against a known answer at all.

Reporting that limitation is the point of this section. An estimator validated on
the one case that can be simulated should not be presented as validated for the
cases that cannot.

## What follows

For anyone setting a pump, nothing changes. The entered constants remain the
number to use and the earlier outcome work is the reason.

For anyone deriving sensitivity from device data, which includes every dynamic
sensitivity equation in current use, the finding is that the quantity being derived
is defined by the construction as much as by the person. Two defensible choices of
denominator on the same overnight records differ by more than threefold, and the
choice that ignores insulin delivered during the window reads highest exactly where
an algorithm is delivering it. An equation calibrated against one construction
inherits its bias and there is no way to read that bias off the published number.

For anyone reporting such an estimate, the practical recommendation is to publish
the denominator in full and the ratio against entered settings alongside the
estimate. Every construction here looks defensible in isolation and they span 0.30
to 1.12, so neither the estimate nor the ratio means anything without the
construction attached.

The remaining question is why the oref archive returns 0.41 under a construction
that returns 0.72 to 1.12 on trial data. That is now a question about those
records rather than about method, and one difference is not yet held fixed: the
oref figure came from fitting a constant over the square root of daily dose across
users, where these are medians of per-person ratios.

## Method

Analysis is local, in Python, against the Jaeb releases loaded into PostgreSQL
17.9 with TimescaleDB 2.26.1. The corpus holds 132,467,585 CGM readings,
68,995,020 basal records and 5,154,053 bolus records.

Delivery is taken net of the programmed basal throughout, since basal offsets
hepatic glucose output rather than lowering glucose. The programmed rate comes
from the recorded schedule in Loop and REPLACE-BG, from the 48 half-hourly rates
on the case report form in DCLP3, DCLP5 and PEDAP, and from each person's own
median delivery at that half hour in IOBP2, which programmes no basal at all.

Entered settings are the per-person median of what was recorded at the time of
dosing or on the case report form. They were checked against the pump's own
arithmetic: the sensitivity implied by dividing the glucose excess by the
recommended correction agrees with the recorded column at a ratio of 1.000 in
Loop.

Code is at https://github.com/tim2000s/dynamic-isf-calculations. The assembled
ratios are produced by `inv009/ratio_evidence.py`, the estimators by
`inv009/correction_landing.py` and `inv009/correction_routes.py`, and the
simulation check by `inv009/audit_estimators.py`.

## Acknowledgements

The source of the data is the Loop Study (sponsored by the Jaeb Center for Health
Research and funded by the Helmsley Charitable Trust), but the analyses, content
and conclusions presented herein are solely the responsibility of the authors and
have not been reviewed or approved by the study sponsor.

The source of the data is the Insulin Only Bionic Pancreas Pivotal Trial, but the
analyses, content and conclusions presented herein are solely the responsibility
of the authors and have not been reviewed or approved by the Bionic Pancreas
Research Group or Beta Bionics.

Data from REPLACE-BG, DCLP3, DCLP5 and PEDAP were supplied by the Jaeb Center for
Health Research. The author thanks the participants of all six cohorts and the
curators of the OpenAPS Data Commons.
