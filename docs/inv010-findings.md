---
title: "OREF-INV-010: autosens against a real device"
subtitle: "Three phases run on the Boost cohort, and what the first one costs"
author: "Tim Street"
date: "26 August 2026"
---

## Summary

Phase one was the gate. The autosens reconstruction built for INV-009 does not
track what a device actually computed, so the detector comparison in that work
should not be relied on. That comparison never reached the published paper or its
summary, so nothing needs retracting, but the conclusion it reached about
widening the autosens clamp does not stand.

Phases two and three use only values the device wrote down, so they hold whatever
the reconstruction does. Both produced findings worth reporting.

Four of the ten users show a ratio that never varies from 1.000, and Boost's own
guidance is to switch autosens off, so those are switched off rather than finding
nothing. Among the six where the computation is demonstrably live, autosens is
neutral for a median of 39% of decisions rather than the 90% a naive average
across all ten gives.

Where it does move, it agrees only weakly with the dose-derived ratio Boost
applies alongside it, and when both are away from neutral they point in opposite
directions a third of the time.

## A correction to the scope

The scoping note reported saturation rates of 47.7% and 50.7% at the ceiling for
two users. Those figures came from `sensitivityRatio`, and
`boostAutosens_mode` reads 'tdd' for every user here, so that field is Boost's
own dose-derived ratio rather than autosens. The oref autosens value sits in
`boostAutosens_orefRatio`. Everything below uses the latter.

## Phase one: the reconstruction does not match the device

Ten users carry a recorded oref autosens ratio. Four of them return exactly
1.000 for every record, leaving six with a signal to compare against.

| User | Correlation | Sign agreement | Implied max basal | Observed |
|---|---|---|---|---|
| C | +0.670 | 76% | 3.86 | 10.82 |
| E | +0.516 | 43% | 0.76 | 3.00 |
| tim | +0.434 | 2% | -0.33 | 3.20 |
| H | +0.330 | 56% | 0.34 | 4.80 |
| F | +0.269 | 2% | 0.04 | 3.15 |
| I | +0.134 | 68% | 1.47 | 5.00 |

Median correlation is +0.382 and median sign agreement is 49%, which is a coin
toss. For two users the reconstruction points the opposite way to the device
almost every time, and the maximum daily basal implied by matching the device
comes out negative for one of them, which is not a quantity that can be negative.

The fault is not the sensitivity fed in. Substituting the profile sensitivity for
the dynamic one moves the correlation from +0.142 to +0.129 for one user and
from -0.059 to -0.049 for another. The reconstruction's deviations run
systematically positive, at a median of +1.4 to +2.7 mg/dL per five minutes for
the users it disagrees with most, while the device reads neutral over the same
periods.

What this costs: the INV-009 detector comparison, and with it the finding that
widening the autosens clamp from 0.7 to 1.2 out to 0.5 to 1.5 was worth about
0.9 mg/dL. That result had already required three corrections, the last of which
inverted it. It should now be treated as unsupported.

## Phase two: how much autosens moves, once the switched-off users are removed

An earlier version of this note reported that autosens is neutral for 90% of
decisions across all ten users. That figure is wrong to quote, because Boost's
guidance is to switch autosens off when using it, and four of the ten users show
a ratio that takes exactly one distinct value across their whole record. Those
are disabled, not quiet.

Separating them on whether the value ever varies:

| User | Distinct values | Neutral share | At ceiling | At floor | Status |
|---|---|---|---|---|---|
| H | 85 | 21.5% | 33.6% | 0% | live |
| I | 51 | 27.3% | 55.2% | 1.3% | live |
| E | 47 | 49.8% | 7.0% | 0.2% | live |
| C | 44 | 27.1% | 0% | 0% | live |
| tim | 31 | 84.7% | 0% | 2.2% | live |
| F | 18 | 96.1% | 0% | 0.2% | live |
| A, B, D, K | 1 | 100% | 0% | 0% | off |

Among the six live users the neutral share has a median of **39%**, so autosens
is adjusting something for roughly three fifths of decisions. That is a very
different picture from the one the pooled figure gave, and it is the one to
carry forward.

Two of the six are still mostly neutral. tim at 84.7% and F at 96.1% look more
like the disabled group than like their live peers, and the plugin fades its
ratio toward 1 when it has under an hour of valid deviations, so a record with
frequent carbohydrate entries can sit near neutral while running normally. Those
two should not be read as evidence either way.

Where autosens is active it reaches both bounds, which the INV-009 reconstruction
never managed: I sits at the ceiling 55.2% of the time and H 33.6%, while tim
reaches the floor 2.2% of the time.

## Phase three A: what autosens is measuring under dynamic ISF

Under Boost with dynamic ISF, autosens computes its deviations against the
profile sensitivity while the loop doses on the dynamic value. The two are
different numbers, often by a lot, so the reference autosens judges against is
not the one in use.

That has a consequence which can be checked. If a person's true sensitivity is
nearer the dynamic value, autosens predicts the glucose fall from the profile
value and is wrong by the gap between them. A profile sitting above the dynamic
value makes it over-predict the fall, deviations go positive, and the ratio
rises. So the ratio should track that gap upwards.

| User | Median profile ISF | Median dynamic ISF | Profile / dynamic | Correlation with autosens |
|---|---|---|---|---|
| C | 30.6 | 54.9 | 0.56 | +0.628 |
| tim | 122.4 | 116.4 | 0.73 | +0.391 |
| H | 36.6 | 46.6 | 0.77 | +0.212 |
| I | 85.0 | 73.6 | 1.15 | +0.001 |
| E | 70.0 | 49.5 | 1.41 | +0.330 |
| F | 97.2 | 54.0 | 1.80 | +0.232 |

The correlation is positive for all six, with a median of +0.281. The gap itself
is large: for one user the profile sits at 56% of the dynamic value and for
another at 180% of it.

So a meaningful share of what autosens reports under dynamic ISF is the distance
between two settings rather than anything about the person. It is not the whole
of it, since correlations of 0.28 leave most of the variance elsewhere, but it is
enough that the signal should not be read as a sensitivity measurement in this
configuration.

This is evidence for the guidance rather than against it. Boost recommends
switching autosens off when using it, and this is a mechanism by which leaving it
on would mislead.

## Phase three B: the two applied mechanisms are largely unrelated

Nine of these users run a dose-derived ratio and oref autosens at once, both
scaling the same sensitivity from different evidence.

| User | Applied ratio | oref ratio | Correlation | Both away from neutral | Opposed |
|---|---|---|---|---|---|
| C | 0.874 | 0.934 | +0.351 | 67.5% | 18% |
| E | 0.835 | 1.022 | +0.304 | 46.0% | 50% |
| H | 0.936 | 1.161 | +0.133 | 71.7% | 49% |
| F | 0.850 | 0.997 | +0.122 | 3.5% | 0% |
| tim | 0.971 | 0.980 | -0.091 | 1.5% | 0% |
| I | 0.849 | 1.117 | -0.035 | 71.4% | 82% |

Median correlation is +0.128. When both are away from neutral they point in
opposite directions for a median of 34% of the time, and for one user 82% of the
time.

So the two are close to independent, and for some people actively opposed. One
reads recent dose, the other reads recent deviations, and they reach different
conclusions about the same person at the same moment. Which is right is not
something this data settles. That they disagree is worth knowing before either
is trusted as a sensitivity signal.

## What follows

Phase three A also narrows what any future validation can use. A cohort running
dynamic ISF cannot validate an autosens reconstruction, because the device's own
autosens is partly reporting a settings gap rather than a sensitivity. The
INV-009 port failed against these users, and this gives a reason why that has
nothing to do with the port.

The INV-009 detector work needs redoing against a device before any of it is
quoted. That means a cohort running oref autosens as the applied mechanism, which
this one is not.

The phase two question resolved once the switched-off users were separated out.
Among users whose computation is live, autosens is adjusting for about three
fifths of decisions, which is a working mechanism rather than a dormant one.

What remains open is whether that generalises. Every user here runs Boost in its
dose-derived mode, so even the live oref ratio is a shadow rather than the number
doing the dosing. Whether autosens behaves the same way when it is the applied
mechanism, on profile sensitivity, is not something this cohort can show.

## Scope

Eleven people, self-selected, running a development branch, over five weeks. This
supports statements about mechanism and about whether a reconstruction is
faithful. It supports no population claim.
