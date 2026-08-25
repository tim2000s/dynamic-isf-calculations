---
title: "Does insulin sensitivity really scale with total daily dose, and by how much?"
subtitle: "The dynamic ISF equations tested against 1,684 people in seven public studies"
author: "Tim Street"
date: "25 August 2026"
---

## Key points

- Across 581 people with a recorded sensitivity factor, entered settings scale
  with total daily dose almost exactly as the 1800 rule says: a log-log slope of
  **-0.98** (95% CI -1.03 to -0.93), with sensitivity times daily dose sitting at
  2067 against the 2139 the v1 equation implies. That is v1's assumption, and it
  holds.
- Measured from how glucose actually responds, the relationship is real but
  **shallower than v1 and nothing like v2**. Between people it is -0.83, and
  within a single person, which is what the equations actually do several times an
  hour, pooling across 1,313 people gives **-0.65** (p = 5e-87).
- **No glucose term earns its place.** Six candidate shapes, including the log
  scalers from both equations, sit within 0.05 mg/dL of each other out of sample,
  and the flat one wins.
- Asked to predict the overnight fall, a static 1800 rule is the best of the lot
  at 34.1 mg/dL. v1 costs 2.6 mg/dL more. **v2 comes in at 140.9 mg/dL** and is
  never the best candidate for a single person out of 1,626.
- Simulated people whose sensitivity really did follow a squared law came back
  through the same pipeline reading -1.56 to -1.85, even with 45% of their meals
  never recorded. Nothing we measured is anywhere near that.

## Why ask this again

The two dynamic ISF equations both say that insulin sensitivity falls as total
daily dose rises, and disagree about how fast. The original, which I will call
v1 throughout, has sensitivity going as one over daily dose. The revised maths,
v2, has it going as one over the square of daily dose. Between two people whose
doses differ threefold, v1 predicts a threefold difference in sensitivity and v2
a ninefold one. That is not a detail, and the two cross at around 64 units a day,
so which one you are on determines whether the newer equation makes your
corrections stronger or weaker.

Earlier work on this question, my own included, ran on loop data: a few hundred
people, mostly on doses between 20 and 80 units a day, with sensitivity inferred
by leaning on the loop's own predictions. That work suggested the dose
relationship was real but shallower than either equation assumed. It could not
close the question, for four reasons that were structural rather than fixable.
The numbers were small. The range of doses was narrow, so there was little
leverage either side of the crossover. Every estimate depended on the loop's own
insulin model, which is the thing being questioned. And there was no open-loop
comparison at all, which matters more than it sounds and is the subject of a
whole section below.

The public study archives change all four. What follows uses seven of them,
covering roughly 1,700 people, ages 2 to 82, daily doses from 7 to 107 units, and
four different controllers, one of which is not a controller at all.

## The data

| Study | What people were using | People | Carbohydrate logged | Median days | Median dose |
|---|---|---|---|---|---|
| Loop | Do-it-yourself Loop, static sensitivity | 842 | yes | 397 | 38.6 U |
| REPLACE-BG | **Pump and CGM, no automation** | 196 | yes | 246 | 41.9 U |
| DCLP3 | Control-IQ, adults and adolescents | 112 | no | 183 | 50.0 U |
| DCLP5 | Control-IQ, ages 6 to 13 | 100 | no | 203 | 37.4 U |
| PEDAP | Control-IQ, ages 2 to 5 | 98 | no | 273 | 13.6 U |
| IOBP2 | Bionic pancreas | 336 | no | 93 | 49.0 U |

FLAIR is in the archive too but its release carries no pump file, so it
contributes daily totals and nothing this work can use.

Two things about that table do most of the work later. REPLACE-BG ran in 2015 on
a pump and a sensor with nobody's algorithm in between, which makes it the only
cohort where the insulin someone received was not chosen by a machine watching
their glucose. And Loop and REPLACE-BG are the only two where people recorded
what they ate, so they are the only two where a fasting window can be screened
rather than guessed at.

The rest is roughly 2.2 million candidate overnight windows, of which 708,106
pass the fasting screen, from 1,659 people with at least forty each.

## Measuring sensitivity when nobody wrote it down

None of these archives contains an insulin-on-board figure, a loop prediction, or
in most cases a sensitivity setting. Everything has to be rebuilt from delivered
basal and boluses.

A four-hour window starting on the hour between 11pm and 3am gives a fall in
glucose. Against that we need the insulin that acted during those four hours,
which comes from convolving the delivery record with an insulin model. Loop
subjects get the models Loop itself shipped, taken from the LoopKit source rather
than approximated, and everyone else gets the oref exponential. The fall divided
by the acting insulin, holding the starting glucose and where it had been, is a
number in milligrams per decilitre per unit, which is what a sensitivity factor
is.

Two details in that are load-bearing.

**Insulin is split by when it was decided.** Insulin already in the body when a
window opens was committed before anything in that window happened, so it cannot
be a response to the glucose the window goes on to measure. Insulin delivered
during the window can be, and under a controller usually is. The two are kept in
separate columns throughout.

**The screen only looks backwards.** A window is kept or dropped on its starting
glucose, on how complete the sensor trace is, and on how long since the person
last ate. It is not screened on what glucose then did. It is tempting to throw
out nights where glucose rose, on the grounds that something must have been
eaten, and it is wrong: those are disproportionately the nights insulin worked
least well, and dropping them does not clean the estimate, it inflates it. Those
filters exist in the code as a sensitivity arm and they move nothing important.

### The check that the arithmetic is right

Loop recorded its own insulin-on-board figure at every bolus, and that is the one
place in this entire dataset where an app's internal state was written down. It
can be recomputed from the pump record and compared with what the app actually
believed.

Across 158 people with enough of those records, the reconstruction tracks Loop's
own figure at a **median correlation of 0.927**, with a median absolute
difference of 0.32 units. That is the dose reconstruction validated against
ground truth rather than against my own assumptions.

Doing it needed something the archives do not obviously contain. Loop counts
insulin on board net of scheduled basal, so a temp basal contributes only the
difference between what it delivered and what the profile would have delivered
anyway. Rebuilding that needs the basal profile, which is not in any settings
file, but is recoverable: every temp basal record carries the rate it suppressed.
Six and a half gigabytes of raw text collapses to 1.6 million rate changes and
gives 842 people's basal schedules.

I should be straight about the limit of that check. Every candidate insulin model
shows a small positive level bias, so the fastest-decaying one wins by absorbing
it, and which model each person was running is therefore only weakly identified.
The correlation is the finding; the model choice is not, and nothing downstream
depends on it.

## What people had entered

The narrowest question first, and the cleanest. Of the 596 people with a
recorded sensitivity factor, 581 also have thirty or more complete days of pump
record to take a daily dose from.

Their entered settings scale with dose at a log-log slope of **-0.979**, with a
95% interval from -1.030 to -0.934. Sensitivity times daily dose has a median of
2067 against the 2139 that v1 implies at a normal target. The Spearman
correlation between that product and dose is **+0.024**, which is what a rule
holding looks like. Run the same test on v2's squared law and it comes out at
**+0.859**, which is what a rule failing looks like, consistently, across the
whole range.

![Entered sensitivity against daily dose](../charts/inv009/fig_entered_scatter.png)

This needs a caveat stated plainly rather than buried, because it is close to
circular. An entered sensitivity factor is a decision, not a measurement. People
and their clinics reach it partly through the 1800 rule, which is the same rule
v1's exponent came from. Finding that entered settings sit on 1800 over dose is
therefore partly a finding about the rule's grip on practice. What it is not is
circular in the other direction: nothing in clinical practice pushes people
towards a squared law, and nothing in the data looks like one.

Two smaller results are worth keeping. Within a single age band the slope is
shallower, between -0.68 and -1.10, so some of the tidy -0.98 comes from pooling
small children with adults. And carb ratios scale at -0.63 with a median product
of 411, so the 500 rule is a looser fit to practice than the 1800 rule is.

## What the glucose response says about dose

Now the harder question: not what people entered, but what a unit of insulin
actually did.

![Every exponent we measured](../charts/inv009/fig_exponents.png)

Between people, using only insulin that was already committed, the exponent is
**-0.83** with a 95% interval of -1.23 to -0.36. Fitting the power law directly
across everyone, including those whose estimate came out negative, gives -0.75.

Within a single person the exponent is **-0.645**, with a standard error of
0.033, from pooling 1,313 people. 73% of them share the sign. This is the version
that has never been tested and the version that matters most, because a law
holding across a population does not have to hold inside one person, and it is
inside one person that these equations run, several times an hour, off a dose
figure blended from the last few hours and days.

Both numbers sit between the square root and v1, and neither is near v2.

### The confound, and why one small cohort carries weight

There is a reason to distrust all of this, and it is worth showing rather than
arguing about.

Under a controller, insulin is close to a function of recent glucose. More
insulin is given exactly when glucose is refusing to come down. Regress the fall
on all the insulin that acted and the slope goes negative, which would mean
insulin raises glucose. It does not; the estimate is measuring the controller's
policy rather than the person's physiology.

![The confound, made visible](../charts/inv009/fig_endogeneity.png)

Reading across the cohorts, the share of people whose sensitivity comes out
positive collapses from 78% to 8% as the controller gets more reactive, and
recovers when only predetermined insulin is used. REPLACE-BG is the exception in
the way the design predicts: it is the only cohort where including all the
insulin makes the estimate better rather than worse, because in 2015 there was no
algorithm choosing it.

That is why a 196-person cohort from a decade ago carries weight out of
proportion to its size, and why it is reported separately everywhere. Its own
within-person exponent is **-1.38**, steeper than the pooled figure and the
closest any cohort comes to v1. It does not come close to v2 either.

## Does sensitivity change with glucose?

Both equations say it does, through a logarithmic term that makes a unit of
insulin do less when glucose is high.

Testing that has a trap in it. High glucose falls further than low glucose
whatever insulin is doing, through mass action, renal clearance and simple
regression to the mean, and near target the body actively defends against a
further fall. So a model that lets glucose explain the size of the fall, and then
asks whether the fall per unit of insulin also depends on glucose, will answer
the first question and report it as the second. Every fit here therefore carries
an additive glucose term, and what carries the claim is the interaction between
insulin action and glucose, which is the only thing that means sensitivity itself
moved.

Two answers come back.

The exponent, pooled across 1,528 people, is **-2.67** and comfortably
significant. The sign is the one neither equation expects: net sensitivity rises
with glucose rather than falling. Only 40% of people show it going the other way.
This is consistent with what I found on loop data last year and with the reasons
given there, so I will not relitigate it here beyond noting that the effect the
equations are reaching for is competing against clearance that rises with glucose
and against counter-regulation near target.

The more useful answer is the practical one.

![No glucose shape earns its place](../charts/inv009/fig_glucose.png)

Out of sample, across 1,616 people, six candidate shapes sit within 0.05 mg/dL of
each other, and the flat one wins. Both equations' log scalers are among the
losers. The decision rule set before running any of this asked for an improvement
of half a milligram per decilitre before a glucose term could be said to earn its
place. Nothing came within a tenth of that.

## Predicting the night

The exponents are the science. This is the practical question: if you had used
each equation's sensitivity to predict how far glucose would fall overnight, how
wrong would you have been.

Each candidate gets a per-person intercept fitted on the first 70% of that
person's nights and is scored on the rest. That favours the equations rather than
handicapping them. Overnight most insulin is basal and is there to offset
endogenous glucose production, and charging that to the sensitivity factor would
blame an equation for something it was never meant to supply.

![Predicting the overnight fall](../charts/inv009/fig_head_to_head.png)

| Candidate | Median error | Median bias |
|---|---|---|
| Static 1800 rule | **34.1 mg/dL** | +2.9 |
| Best single number for that person | 34.4 | -1.4 |
| v1 equation | 36.7 | +0.9 |
| The person's own entered setting | 37.3 | -0.2 |
| 355 over the square root of dose | 37.9 | +3.8 |
| v2 equation | **140.9** | +41.1 |

Of 1,626 people, the best candidate is the best single number for 638, the static
1800 rule for 587, the square root form for 195, v1 for 135 and their own entered
setting for 71. **It is v2 for nobody.**

v2 fails hardest where the consequences are worst. For people on under twenty
units a day its median error is 322 mg/dL, because a squared law hands a small
child a sensitivity of several thousand milligrams per decilitre per unit and
therefore a correction dose near zero. Its error falls as dose rises, to 83
mg/dL above 64 units a day, which is the crossover doing what the algebra says it
should. Even there it is more than twice as wrong as anything else on the list.

## What would have shown up if the equations were right

Every number above is an estimate of something nobody observed, using a
reconstruction of insulin that is certainly imperfect. The fitted sensitivities
come out well below what people have entered, so something is attenuating them.
The question that matters is whether that attenuation also bends the exponent,
and the way to find out is to build people whose sensitivity really does follow a
law we chose and run them through the identical pipeline.

![Recovering a known exponent](../charts/inv009/fig_synthetic.png)

| True exponent | Recovered, open loop | With a reactive controller | With 45% of meals never recorded |
|---|---|---|---|
| -0.50 | -0.50 | -0.42 | -0.72 |
| -1.00 | **-1.02** | -0.80 | -1.03 |
| -2.00 | **-1.57** | -1.85 | -1.56 |

The pipeline recovers shallow exponents almost exactly and compresses steep ones
somewhat. The number that matters is on the bottom row: a squared law comes back
reading between -1.56 and -1.85 under every condition tried, including one where
nearly half of all meals go unrecorded and there is a third of a unit of error on
what actually acted. What we measured is -0.65 to -0.88. A squared law cannot
produce that.

Attenuation in the simulation runs between 0.38 and 0.63 of the true
sensitivity, which is less severe than the real gap against entered settings. The
level is not something this method recovers well and I am not claiming it does.
The exponent is.

## What this does not settle

**The level of sensitivity is not measured here.** The fitted numbers are
attenuated, so the constant in any equation cannot be read off this work. Only
the shape of the relationship can.

**The closed-loop cohorts are weak evidence on their own.** The share of people
with a positive sensitivity is barely better than a coin flip in the Control-IQ
studies, and for the youngest cohort it is below one. Those cohorts replicate a
direction; they do not establish it. The open-loop arm and the entered settings
carry the argument.

**Everything here is overnight and fasting.** This says nothing about the hours
after a meal, which is where a glucose term might still earn its place and where
carbohydrate that the absorption model has not fully caught is a live
alternative explanation for anything that looks like changing sensitivity.

**Sensitivity that varies with the day rather than the dose is out of scope.**
Exercise, illness and a change in how someone eats all move sensitivity on
timescales of hours to days. Nothing here measures those, and nothing here argues
they are not real.

**Two problems in the underlying data were found and repaired**, and anyone
reusing these archives should know about them. Extended bolus durations reach the
database in milliseconds for three of the six studies and in seconds for the
other two, so a one-hour square wave arrives as a thousand hours; spreading a
bolus over that put insulin into 41,591 of one person's 46,176 five-minute
periods. And the database returns timestamps at microsecond resolution, which
silently divided every basal total by a thousand until the reconstruction was
checked against the database's own sums.

## Where this leaves the two equations

v1's assumption about dose is close to right, and it is close to right in the
place it was drawn from: it matches what people actually run their pumps at,
almost exactly. Measured against how glucose responds it is a little steep, with
the honest range somewhere between the square root and v1, but it is the right
shape and the right order of magnitude. Its glucose term does not earn its place,
though it costs little.

v2's squared law is not supported by anything measured here. Entered settings
contradict it, the glucose response contradicts it, it is never the best
predictor for any of 1,626 people, and the simulation says a squared law could
not have hidden as anything we saw. The failure is worst for children and for
anyone on a small daily dose, which is the group least able to absorb a
correction dose that is far too small.

None of this is a criticism of the people who built these equations, who were
working from the evidence available and from a rule of thumb that this work finds
holds up well. It is a straightforward statement that a much larger and more
varied group of people, including one that used no algorithm at all, does not
show the relationship the revised maths assumes.

If a dose term is wanted, the evidence here supports an exponent between about
-0.65 and -1.0, applied to the person's own level rather than used to set it.
Given that a static 1800 rule outpredicted every dynamic form tested, the more
interesting question is not which exponent to use but whether the sensitivity
factor is the right place to put this at all.

## Method and reproduction

Code is in `inv009/` in the `dynamic-isf-calculations` repository, with 28 unit
tests. The window cache rebuilds in about nine minutes on a laptop:

```
psql -d oref -f inv009/ingest/sql/05_settings.sql
python3 inv009/ingest/load_settings.py
python3 -m inv009.build_cache
python3 -m inv009.entered_isf && python3 -m inv009.effective_isf
python3 -m inv009.tdd_axis && python3 -m inv009.glucose_axis
python3 -m inv009.head_to_head && python3 -m inv009.synthetic
python3 -m inv009.loop_model_infer && python3 -m inv009.figures
```

Windows are four hours from starts between 11pm and 3am, needing 80% sensor
coverage, a starting glucose between 90 and 300 mg/dL, and either six hours clear
of logged carbohydrate or, where none is logged, four hours clear of a bolus
large relative to that person's own daily dose. Per-person fits use ordinary
least squares with a robust standard error; population figures pool per-person
estimates with DerSimonian and Laird, as in the earlier work, so the two can be
read against each other.

## Data attribution

The source of the data is the Loop Study (sponsored by the Jaeb Center for Health
Research and funded by the Helmsley Charitable Trust), but the analyses, content
and conclusions presented herein are solely the responsibility of the authors and
have not been reviewed or approved by the study sponsor.

The source of the data is the Insulin Only Bionic Pancreas Pivotal Trial, but the
analyses, content and conclusions presented herein are solely the responsibility
of the authors and have not been reviewed or approved by the Bionic Pancreas
Research Group or Beta Bionics.

DCLP3, DCLP5, PEDAP and REPLACE-BG are also public datasets from the Jaeb Center
for Health Research, whose releases give no attribution wording of their own. The
analyses, content and conclusions here are solely mine and have not been reviewed
or approved by any study sponsor.

Dates in these archives are shifted or rebased per participant, so time of day is
intact and calendar date is not. Nothing here depends on calendar date.
