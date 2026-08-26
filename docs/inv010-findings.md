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

In this cohort oref autosens sits at exactly 1.000 for a median of 90% of
decisions, and never moves at all for four of the ten users. Where it does move,
it agrees only weakly with the dose-derived ratio Boost applies alongside it, and
when both are away from neutral they point in opposite directions a third of the
time.

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

## Phase two: autosens is mostly not moving

| User | Exactly 1.000 | At ceiling | At floor | Active |
|---|---|---|---|---|
| A, B, D, K | 100% | 0% | 0% | 0% |
| F | 96.1% | 0% | 0.2% | 3.9% |
| tim | 84.7% | 0% | 2.2% | 15.3% |
| E | 49.8% | 7.0% | 0.2% | 50.2% |
| I | 27.3% | 55.2% | 1.3% | 72.7% |
| C | 27.1% | 0% | 0% | 72.9% |
| H | 21.5% | 33.6% | 0% | 78.5% |

Median across the ten: neutral for 90% of decisions.

Two readings are available and this data cannot separate them. Either autosens
is examining the record and finding nothing worth adjusting, or it is not being
fed in this configuration and returns neutral by default. Four users at exactly
100% neutral points at the second. Every user here runs Boost in its
dose-derived mode, where the oref ratio is computed alongside rather than
applied, so a plugin that is not receiving data would look exactly like this.

Distinguishing them needs a user running oref autosens as the applied mechanism.
None of these ten is.

Where autosens is active it does reach both bounds, which the INV-009
reconstruction never managed: I sits at the ceiling 55.2% of the time and H
33.6%, while tim reaches the floor 2.2% of the time.

## Phase three: the two mechanisms are largely unrelated

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

The INV-009 detector work needs redoing against a device before any of it is
quoted. That means a cohort running oref autosens as the applied mechanism, which
this one is not.

The phase two question is the more interesting one and is unresolved for a
mundane reason. If autosens really is neutral 90% of the time in normal use, the
mechanism is doing much less than its reputation suggests, and the INV-009
finding that per-person base error is large and mostly uncorrected would matter
more, not less. If instead it is a configuration artefact, the number means
nothing. One user on profile sensitivity with autosens applied would separate
them.

## Scope

Eleven people, self-selected, running a development branch, over five weeks. This
supports statements about mechanism and about whether a reconstruction is
faithful. It supports no population claim.
