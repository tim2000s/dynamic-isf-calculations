---
title: "OREF-INV-010: what the Boost cohort can and cannot answer about autosens"
subtitle: "Scoping note"
author: "Tim Street"
date: "26 August 2026"
---

## Why this was scoped

INV-009 established that a sensitivity factor derived from a seven day total
daily dose needs a real-time adjustment on top, and that per-person base error is
large: entered settings run from 0.68 to 2.07 of what the 1800 rule gives, and
24% of people sit beyond what even a widened autosens clamp could absorb.

The design question that followed was whether autosens output can drive a slow
loop that titrates the base. Testing it on Loop found nothing, but that test had
two holes. Loop has no autosens, so nobody there is acting on one, and the
reconstruction I used never reached the lower clamp, which left half the
hypothesis untestable.

The Boost cohort has what Loop lacks: autosens ratios recorded by real devices.
This note establishes what it can support before any work starts.

## What is in the data

Two tables. `boost_decisions` holds 2.1 million per-decision rows from 36 users
between August 2025 and August 2026. `boost_devicestatus_raw` holds the openaps
block, including the recorded `sensitivityRatio`.

| Quantity | Coverage |
|---|---|
| Recorded `sensitivityRatio` | 200,118 readings, 11 users, 2026-07-18 to 2026-08-23 |
| Profile sensitivity recorded | 8 users, spans of 82 to 175 days |
| Users with any profile sensitivity change above 10% | 3 (H with 12, tim with 9, B with 6) |
| Those users' time on dynamic ISF | 99 to 100% |
| Users on profile sensitivity during the ratio window | 2 (G at 0% dynamic, H at 27%) |

Saturation is common and both bounds are reached, which the Loop reconstruction
never managed:

| User | Mean ratio | At ceiling | At floor |
|---|---|---|---|
| G | 1.205 | 47.7% | 0.0% |
| H | 1.369 | 50.7% | 6.4% |
| F | 0.850 | 0.0% | 20.2% |
| E | 0.835 | 5.3% | 18.7% |
| I | 0.915 | 10.5% | 16.7% |
| tim | 0.971 | 0.0% | 2.3% |

## What this cohort cannot answer

The question that prompted the scoping is not viable here. Asking whether
autosens saturation precedes a change to the sensitivity factor requires people
whose dosing sensitivity IS that setting, and who change it.

Those two groups do not overlap. Every user with a recorded profile sensitivity
runs dynamic ISF for 99 to 100% of their decisions, so their profile value is a
fallback rather than the number doing the work. The two users on profile
sensitivity during the ratio window are G, who has no profile sensitivity
recorded at all, and H, at 27%. That is one usable person.

Recruiting for it would need users on profile sensitivity with autosens, over
enough months to catch settings changes. Nothing in the current cohort supplies
that, and no amount of analysis will change it.

## What this cohort can answer, in order of value

### One: does the autosens reconstruction work at all

INV-009 built an autosens port and compared detectors with it. That port needed
correcting three times: it used total delivered insulin where AAPS uses insulin
net of the scheduled basal, it applied no meal exclusion to the four cohorts that
record no carbohydrate, and it scored against an uncalibrated base. Each
correction changed the conclusion, and the last one inverted it.

It has never been checked against a device. Here there are 200,118 recorded
ratios from eleven people over five weeks, alongside the pump records the port
consumes.

The work is to run the port over the same people and periods and compare, per
person: correlation of the reconstructed ratio against the recorded one, agreement
on which bound is being approached, and agreement on the fraction of time spent
saturated. A port that tracks the device validates the INV-009 detector
comparison. One that does not invalidates it, which is worth knowing either way
and is cheap to establish.

This is the piece that decides whether the earlier work stands.

### Two: how often does autosens actually saturate

Nobody appears to have characterised this. The table above suggests it is common
and asymmetric, with two users pinned at a bound around half the time.

That matters directly for the clamp question INV-009 raised and could not settle.
If a user spends half their time at the ceiling, the clamp is not a safety margin
that occasionally binds, it is the operating point, and the loop is running a
sensitivity the detector disagrees with for half the day. Whether that is
protective or is suppressing a real signal is exactly what the outcome data
alongside it can address, since `boost_cgm` carries paired glucose.

The work is descriptive: saturation rates per user, by time of day, by whether
carbohydrate was recently entered, and paired against time in range and time
below range during and after saturation episodes.

### Three: what dynamic ISF and autosens do to each other

Nine of the eleven users run both at once. Dynamic ISF recomputes sensitivity
from recent total daily dose several times an hour. Autosens then multiplies it
from recent deviations. The two mechanisms are adjusting the same quantity from
different evidence, and nothing establishes whether they agree, are independent,
or oppose each other.

INV-009 found that the figure dynamic ISF uses swings by about 41% of a daily
total inside a single record. If autosens is partly correcting that swing, the
two are fighting, and the combination is doing something neither was designed for.

The work is to pair the recorded `sensitivityRatio` against the recorded
`variable_sens` and `dynamic_isf`, per decision, and measure whether they move
together or against each other.

This is exploratory rather than confirmatory. It has no preregistered prediction
and should be reported as description.

## Scope and cost

Phase one is the validation and is the gate. Roughly a day, reusing
`inv009/sensitivity_detectors.py` against `boost_devicestatus_raw`. Phases two
and three are each two to three days and should only start if phase one shows the
reconstruction is sound, since both lean on the same machinery.

Data is local and no new extraction is needed. The 2026-06-26 refresh note in the
extraction records says `extract_boost.py` is idempotent per user, so extending
the `sensitivityRatio` window forward is a re-run rather than new work, and would
widen phase one's overlap.

## What would need saying in any output

The cohort is eleven people, self-selected, running a development branch, over
five weeks. It supports statements about mechanism and about whether a
reconstruction is faithful. It does not support population claims, and it is not
the group to draw a clamp recommendation from on its own.
