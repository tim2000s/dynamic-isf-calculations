---
title: "Does insulin sensitivity vary with both glucose and total daily dose, and is there a clear relationship between the three?"
subtitle: "Testing both dynamic ISF equations against 1,684 people in seven public study archives"
author: "Tim Street"
date: "25 August 2026"
---

## What this is about

Somewhere in your pump settings is a number that says how far one unit of insulin
drops your glucose. Most of us arrived at it from the 1800 rule, which divides
1800 by your total daily dose, then adjusted it when corrections kept landing too
hard or too soft.

Dynamic ISF replaces that fixed number with a calculated one. Instead of using
what you entered, the algorithm works out a fresh sensitivity every few minutes
from two things: how much insulin you have been using lately, and where your
glucose is right now.

There are two versions of the calculation and they disagree about how much your
daily dose should matter. The original, v1, keeps the 1800 rule's shape: if your
daily dose doubles, your sensitivity halves. The revised version, v2, squares the
dose term, so if your daily dose doubles, your sensitivity drops to a quarter.

That difference is not academic. Take someone on 20 units a day and someone on 60.
v1 says the second person's corrections should be three times smaller per unit.
v2 says nine times smaller. The two versions cross at around 64 units a day, so
which one you are running decides whether the newer maths makes your corrections
stronger or weaker than the older one.

Both also reduce your sensitivity when your glucose is high, on the reasoning that
insulin works less well when you are running high.

This is an attempt to check all of that against a lot of people.

## The short version

Your daily dose does matter, and roughly in the way v1 says. Sensitivity falls as
daily dose rises, at an exponent close to what the 1800 rule assumes, though a
little shallower. Four separate methods agree on this, including one that was
free to find no relationship at all.

Your glucose level also matters, but in the opposite direction to what both
equations assume. They make insulin most effective when you are near target and
least effective when you are high. What I measured is the reverse. A unit does
least near target, which is where your body is actively pushing back against
going lower, and slightly more when you are high.

The two do not combine the way either equation combines them. Both multiply a
dose term by a glucose term and assume the two are independent. They are not: how
your sensitivity changes with glucose depends on how much insulin you use.

And a real part of what the equations are picking up is food rather than
sensitivity. In the four hours after a meal you logged, measured sensitivity goes
negative, meaning glucose rises while insulin is acting, because carbohydrate is
still arriving after the absorption model has finished counting it.

When I asked each version to predict how far glucose would actually fall
overnight, a plain static 1800 rule did better than either of them, and v2 did
very badly indeed.

## What that means if you are running one of these

If you use v1, the dose half of it is close to right and the glucose half is
pointed the wrong way but is small enough not to do much harm.

If you use v2, the squared dose term is not supported by anything I can measure
here. It was not the best predictor for a single person out of 1,626, and it fails
worst for people on small daily doses. For a young child on ten units a day,
squaring the dose gives a sensitivity of several thousand, which means a
correction dose of almost nothing.

If you use neither, the 1800 rule you probably already have is a reasonable
description of what these archives show, and it beat every dynamic form I tested
at predicting the overnight fall.

None of this is a criticism of the people who built these equations. They were
working from the same rule of thumb most of us set our pumps from, and one of the
findings here is that the rule holds up well.

## Where the data came from

Seven public study archives, about 1,700 people, ages 2 to 82, daily doses from 7
to 107 units.

| Study | What people were using | People | Food logged | Median days | Median dose |
|---|---|---|---|---|---|
| Loop | Do-it-yourself Loop, fixed sensitivity | 842 | yes | 397 | 38.6 U |
| REPLACE-BG | Pump and sensor, nothing automating | 196 | yes | 246 | 41.9 U |
| DCLP3 | Control-IQ, adults and teenagers | 112 | no | 183 | 50.0 U |
| DCLP5 | Control-IQ, ages 6 to 13 | 100 | no | 203 | 37.4 U |
| PEDAP | Control-IQ, ages 2 to 5 | 98 | no | 273 | 13.6 U |
| IOBP2 | Bionic pancreas | 336 | no | 93 | 49.0 U |

That gives 708,106 usable overnight stretches from 1,659 people, each
contributing at least forty.

Two of those rows matter more than the others.

REPLACE-BG ran in 2015 on a pump and a sensor with no algorithm in between. That
makes it the only group where the insulin someone received was not chosen by a
machine watching their glucose, which turns out to matter a great deal.

Loop and REPLACE-BG are also the only groups where people recorded what they ate.
That makes them the only two where I can genuinely exclude a meal rather than
guess at one.

## How I worked out what a unit did

None of these archives records an insulin on board figure, or what the algorithm
predicted, or in most cases even a sensitivity setting. All of it has to be
rebuilt from the basal and boluses that were delivered.

The approach is the one you would sketch on paper. Take a stretch of time
overnight. Note how far glucose fell. Work out how much insulin was actually
acting across that stretch, adding up every dose that contributed, using a
six hour insulin curve. Divide one by the other and you have milligrams per
decilitre per unit, which is what a sensitivity factor is.

Two details in that do most of the work, and both are worth understanding because
they are where this kind of analysis usually goes wrong.

### Insulin your loop chose cannot measure your sensitivity

If an algorithm is running, it gives you more insulin exactly when your glucose
is refusing to come down. So if you compare the fall against all the insulin that
acted, you find that more insulin goes with less falling, and the answer comes out
negative. That would say insulin puts your glucose up. It does not. You have
measured the algorithm's behaviour rather than your body's.

I handle this by separating insulin according to when it was decided. Insulin
already in you when a stretch begins was committed before anything in that stretch
happened, so it cannot be a reaction to the glucose I am about to measure. Insulin
delivered during the stretch can be.

The size of that problem is worth seeing.

![The confound, made visible](../charts/inv009/fig_endogeneity.png)

Reading across the groups, the share of people whose sensitivity comes out
positive collapses from 78% to 8% as the system gets more reactive, then recovers
once I use only the insulin that was already committed. REPLACE-BG is the
exception, because in 2015 nothing was choosing it.

That is why 196 people from a decade ago carry more weight here than their number
suggests.

### Basal is not free

Over any stretch of time, your glucose change is insulin action minus the glucose
your liver is putting out. Basal exists to cancel that output. So when your basal
is right, your glucose is flat, and a naive calculation of fall divided by insulin
gives zero. That is not a sensitivity of zero, it is a correct basal rate.

When I ran that naive calculation it returned 3 to 9 mg/dL per unit against
entered settings of 25 to 60, which is that effect and not a finding. Everything
below therefore compares each stretch against what that person's glucose usually
does at that time of day, and against how much insulin usually acts then. What
survives is the effect of the insulin that was not routine, which is what a
correction dose actually is.

### A check that the arithmetic is right

Loop wrote down its own insulin on board figure at every bolus. That is the one
place in all of these archives where an app's internal state was recorded, so I
can rebuild it from the pump data and see whether I get the same answer.

Across 158 people, my reconstruction tracks Loop's own figure at a correlation of
0.927, with a typical difference of a third of a unit. That is the whole dose
calculation checked against something other than my own assumptions.

Getting there needed the basal profile, which is in no settings file in the
archive, because Loop counts insulin on board net of scheduled basal. It turned
out to be recoverable: every temp basal record carries the rate it was overriding.

## Question one: does your daily dose change your sensitivity?

Yes, and by roughly the amount the 1800 rule assumes.

### What people had entered

Of the 596 people with a sensitivity factor on record, 581 also have enough pump
data to work out a daily dose.

![Entered sensitivity against daily dose](../charts/inv009/fig_entered_scatter.png)

Plotted against dose on a log scale, the slope is -0.979, with a 95% interval
running from -1.030 to -0.934. An exponent of exactly -1 is the 1800 rule, so
that interval sits on it and excludes anything much shallower.

Multiply each person's sensitivity by their daily dose and the median comes to
2067. The 1800 rule, evaluated the way v1 evaluates it, gives 2139. If the rule
holds, that product should not drift as dose changes, and it does not: the
correlation between them is +0.024, which is as close to nothing as you will see
in real data. Run the same test on v2's squared version and you get +0.859, which
is what a rule looks like when it does not hold.

There is a caveat here that needs saying out loud, because it is close to
circular. What you have entered is a decision, not a measurement, and you and your
team probably reached it through the 1800 rule in the first place. So this is
partly a finding about how much grip that rule has on practice.

It is not circular the other way round, though. Nothing in clinical practice
pushes anyone towards a squared law, and nothing in what people run looks like one.

Carb ratios follow the same pattern more loosely. Multiplying each person's carb
ratio by their daily dose gives a median of 411, so the 500 rule describes
practice less tightly than the 1800 rule does.

### What people's glucose actually did

Measuring from the glucose response rather than the settings gives a slightly
shallower relationship. Between people it comes to -0.83. Within a single person,
which is what these equations actually do several times an hour, it comes to
-0.645 across 1,313 people, with a standard error of 0.033 and about three
quarters sharing the direction.

![Exponents measured by each method](../charts/inv009/fig_exponents.png)

For a sense of what those numbers mean: an exponent of -1 is v1, where doubling
your dose halves your sensitivity. An exponent of -2 is v2, where doubling your
dose quarters it. An exponent of -0.65 means doubling your dose multiplies your
sensitivity by about 0.64, so a bit more than half.

Both sit between the square root and v1. Neither is anywhere near v2.

### The same question asked at every reading

The measurements above fit a line through each person's data. The more direct
thing is to work out a sensitivity at every single glucose reading and put it
beside what each equation calculates at that same moment.

At each reading I look back six hours, add up the insulin that acted in that time
from every dose that contributed, and compare the glucose change against what that
person's glucose usually does over the same six hours at that time of day. Only
stretches with no food recorded are used. That gives 1,408,861 readings from 1,004
people.

The level that comes out is too low to read directly, because dividing one noisy
number by another pulls the middle of the distribution towards zero. The
relationship survives that, and there is a way to check it: I ran the same
calculation on the numbers v1 and v2 themselves produce, where the answer is known
in advance.

| What was measured | Exponent recovered | What it should be |
|---|---|---|
| Real sensitivity | -0.843 (95% CI -0.970 to -0.720) | unknown, that is the question |
| v1's own output | -1.026 | -1 |
| v2's own output | -2.085 | -2 |

The method recovers a known answer to within about 0.03. On real data it gives
-0.843, which lands on the -0.83 from the completely different method above.

The same calculation shows how far off the two equations are in absolute terms.
Across all 1.4 million readings, v1 sits at about two and a half times the
measured value and lands within 30% of it 16% of the time. v2 sits at five and a
half times and lands within 30% of it 9% of the time.

So the dose relationship is real, it is close to what the 1800 rule assumes, and
it is nothing like the squared version.

## Question two: does your glucose level change your sensitivity?

Yes, and in the opposite direction to the one both equations assume.

This one has a trap in it. Glucose that starts high falls further than glucose
that starts low whatever insulin is doing, partly through your kidneys clearing it
and partly because things that are high tend to come down. And near target your
body pushes back against going lower. So if you let glucose explain the size of the
fall, and then ask whether the fall per unit also depends on glucose, you will
answer the first question and write down the second.

Every calculation here therefore separates the two, and what carries the claim is
only the part that says a unit of insulin itself did more or less.

![Sensitivity by glucose band against both equations](../charts/inv009/fig_glucose_profile.png)

| Where your glucose was | What I measured | v1 assumes | v2 assumes |
|---|---|---|---|
| 90 to 120 mg/dL | 0.64 | 1.18 | 1.75 |
| 120 to 150 | 1.00 | 1.00 | 1.00 |
| 150 to 190 | 1.10 | 0.87 | 0.72 |
| 190 to 300 | 1.16 | 0.75 | 0.54 |

Those are relative numbers, scaled so the 120 to 150 band is 1.00.

Read the first column downwards. A unit of insulin does about a third less near
target than it does in the middle of the range, then a little more as glucose
rises. Read the other two columns and they go the other way entirely.

The per-reading calculation says the same thing in real units, and it fits
nothing at all.

![Measured sensitivity against both equations at the same reading](../charts/inv009/fig_pointwise.png)

| Where your glucose was | What I measured | v1 calculates | v2 calculates |
|---|---|---|---|
| 70 to 100 mg/dL | -0.0 | 72.9 | 791.5 |
| 100 to 120 | -1.5 | 58.7 | 261.5 |
| 120 to 150 | 1.8 | 49.1 | 158.2 |
| 150 to 190 | 11.3 | 40.5 | 100.3 |
| 190 to 250 | 21.1 | 33.8 | 70.4 |
| 250 to 400 | 26.9 | 30.5 | 58.0 |

Three different methods return that reversal, and the third one fits no curve to
anything.

What it means in practice is that both equations correct hardest where insulin is
working best, and ease off where it is working least. The near-target end is the
one worth thinking about, because that is where your body is defending you against
a hypo, and both equations read that defence as high sensitivity and dose into it.

The good news is that the glucose term is small enough not to matter much. Out of
sample, six different glucose shapes including both equations' own came within
0.05 mg/dL of each other at predicting the overnight fall, and the flat one won. I
had set the bar at 0.5 mg/dL before a glucose term could be said to earn its
place. Nothing came within a tenth of that.

## Question three: do the two work together the way the equations assume?

No.

Both equations multiply a dose term by a glucose term, which quietly assumes the
two do not interact. In plain terms, they assume your glucose profile looks the
same whether you take 15 units a day or 80.

![Glucose profile at each dose band](../charts/inv009/fig_surface.png)

| Daily dose | 90 to 120 | 120 to 150 | 150 to 190 | 190 to 300 |
|---|---|---|---|---|
| Under 25 U | 0.81 | 1.04 | 1.15 | 1.21 |
| 25 to 40 U | 0.53 | 0.75 | 1.14 | 1.13 |
| 40 to 60 U | 0.58 | 0.94 | 1.18 | 1.26 |
| Over 60 U | 0.48 | 0.96 | 1.10 | 1.15 |

For the equations to be right, those rows would have to be identical. They are
not. The dip near target is mild for people on under 25 units a day and around
twice as deep for everyone else. Fitting the interaction directly gives -0.444,
which is about half the size of the dose term itself and comfortably real.

One thing does behave. Comparing across people, someone's own glucose profile
shows no relationship at all with how much insulin they use (a rank correlation
of -0.017 across 971 people). It is within a single
person that the two axes tangle.

## Is any of this really sensitivity, or is it food?

This is the question underneath all of it. Both explanations predict the same
thing. If your sensitivity genuinely falls when you are high and when your recent
dose is large, the equations describe something real. If instead carbohydrate from
a recent meal is still absorbing after the model has finished counting it, then
your glucose is high for that reason, insulin appears to be achieving less than it
should, and your recent dose is large because you ate. The measurement looks
identical either way.

Loop and REPLACE-BG recorded meals, so I can use time since eating to tell them
apart.

![Measured sensitivity by time since the last meal](../charts/inv009/fig_carb_tail.png)

Within four hours of a meal you logged, measured sensitivity is -1.62 mg/dL per
unit. A negative sensitivity is not a real thing. It says glucose went up while
insulin was working, which is what happens when carbohydrate is still arriving
from a meal the absorption model has already written off. It recovers to 6.06 by
six to nine hours out.

That is the absorption tail, measured directly, and it is not small.

![Apparent glucose dependence by time since the last meal](../charts/inv009/fig_carb_glucose.png)

The apparent glucose effect follows the same pattern. Close to a meal it is about
three times what it is in clean fasting. So a real part of what the glucose term
is responding to is food rather than sensitivity.

It does not go away entirely, though. Even well clear of a meal, some glucose
dependence remains, so food does not explain all of it.

The dose relationship is a different story. Holding recent carbohydrate constant,
83% of it survives.

| Group | Before | After holding food constant | Kept |
|---|---|---|---|
| Loop | -0.586 | -0.500 | 85% |
| REPLACE-BG | -1.443 | -0.978 | 68% |
| Together | -0.648 | -0.539 | 83% |

So the dose half is mostly describing something real. The glucose half is
substantially food.

## Putting it to work: predicting the night

The numbers above are the science. This is what happens if you actually use them.

Each candidate gets its own offset for each person, worked out from the first 70%
of that person's nights and scored on the rest. That is generous to the equations
rather than harsh on them, because most overnight insulin is basal, and asking a
sensitivity factor to account for your liver would be blaming it for something it
was never meant to do.

![Predicted overnight fall by candidate](../charts/inv009/fig_head_to_head.png)

| What you use | How far out it was, typically |
|---|---|
| A static 1800 rule | 34.1 mg/dL |
| The best single number for that person | 34.4 |
| v1 | 36.7 |
| Your own entered setting | 37.3 |
| 355 divided by the square root of your dose | 37.9 |
| v2 | 140.9 |

v2 is also biased, not merely imprecise: it overshoots the fall by 41 mg/dL at the
median, meaning it consistently expects insulin to do far more than it does.

Out of 1,626 people, the best option was a single well-chosen number for 638, the
static 1800 rule for 587, the square root version for 195, v1 for 135, and their
own setting for 71. There was nobody it was v2 for.

v2 is worst where it matters most. For people on under twenty units a day its
typical error is 322 mg/dL, for the reason given earlier. The error comes down as
dose goes up, to 83 mg/dL above 64 units a day, which is the crossover doing what
the algebra says it should. Even there it is more than twice as far out as
anything else.

## How do I know the method itself is not producing this?

Fair question, and the way to answer it is to build people whose sensitivity
follows a rule I choose, then run them through the identical machinery and see
what comes back.

![Recovery of a known relationship](../charts/inv009/fig_synthetic.png)

| What I built in | What came back |
|---|---|
| -0.50 | -0.50 |
| -1.00 | -1.02 |
| -2.00 | -1.57 |

Shallow relationships come back almost exactly. Steep ones come back slightly
flattened. The bottom row is the one that matters: a squared law comes back
reading between -1.56 and -1.85 under everything I tried, including a run where
nearly half of all meals went unrecorded. What I measured was between -0.65 and
-0.88. A squared law cannot look like that.

## What this does not tell you

It does not give you a number to set. The sensitivities I measure are shrunk by
the method, so you cannot read a constant off this work, only the shape of the
relationship.

It is all overnight and fasting. The four hours after a meal are measured here
only well enough to show sensitivity going negative. Working out what a meal is
really doing over those hours would need the carbohydrate absorption model itself
under the microscope, which is a different piece of work and probably a more
useful one.

It says nothing about sensitivity that moves with the day rather than the dose.
Exercise, illness and a change in how you eat all move it over hours and days.
None of that is measured here, and none of it is being argued against.

The groups running automated systems are weak evidence on their own. In the
Control-IQ studies the share of people with a positive sensitivity is close to a
coin toss, and for the youngest group it is worse than one. They back up a
direction rather than establishing it.

And two faults in the underlying data had to be repaired along the way. Extended
bolus durations arrive in milliseconds for three of the six studies and in seconds
for the other two, so a one hour square wave reads as a thousand hours. The
database also returns timestamps at microsecond resolution, which divided every
basal total by a thousand until I checked the reconstruction against the
database's own sums. Anyone else working from these archives will meet both.

## Where that leaves the two equations

The dose half of v1 is close to right, and it is closest in the place it came
from: it describes what people actually run their pumps at to within a few per
cent. Measured against what glucose does it looks a touch steep, with the
supported range somewhere between the square root and v1.

Its glucose half is pointed the wrong way, and so is v2's. Both make insulin most
effective near target when the measurement says it is least effective there. The
saving grace is that the term is small enough that it costs little.

For v2's squared dose term I could find no support here. What people run
contradicts it, what glucose does contradicts it, it was not the best predictor
for anyone, and the simulation says a squared law could not have looked like
anything I measured. Its worst behaviour is in children and in anyone on a small
daily dose, who are the least able to absorb a correction that turns out to be far
too small.

I want to be careful about the scope of that. It describes what a large and varied
group of people did. It says nothing about anyone's intent or competence, it is
one analysis of overnight fasting periods from archives never collected for this
purpose, and the equations under examination moved the conversation to where this
work could start. Somebody coming at it another way would be welcome.

If you want a dose term, what is here supports something between about -0.65 and
-1.0, applied to your own level rather than used to set it. But given that a plain
static 1800 rule out-predicted every dynamic version I tested, the question I keep
coming back to is not which exponent to use. It is whether sensitivity is the
right place to be making this adjustment at all, when the thing most obviously
needing attention is what a meal is still doing four hours after you logged it.

## Methods

Windows run four hours from starts between 23:00 and 03:00, requiring 80% sensor
coverage, a starting glucose between 90 and 300 mg/dL, and either six hours clear
of logged carbohydrate or, where none is logged, four hours clear of a bolus large
relative to that person's own daily dose. Screening looks only backwards, never at
what glucose subsequently did, because dropping the nights where glucose rose
removes disproportionately the nights insulin worked least well and inflates the
result.

Insulin action comes from convolving the delivery record with an insulin model.
Loop subjects use the models LoopKit itself ships, taken from its source rather
than approximated. Everyone else uses the oref exponential at six hours with a 75
minute peak.

Per-person fits are ordinary least squares with a heteroskedasticity-robust
standard error. Population figures pool the per-person estimates by the method of
DerSimonian and Laird, as in the earlier work in this series, so the two can be
read side by side. Confidence intervals are percentile bootstrap over people,
using the same weights as the estimate they wrap.

The per-reading calculation looks back six hours at every reading, takes both the
glucose change and the insulin action as departures from that person's own median
for that half hour of the day, and requires the excess insulin action to exceed
0.3 units so that the denominator is not noise.

Separability is tested by fitting insulin action against glucose, against dose,
and against the product of the two, and reading the third coefficient. The glucose
profile is measured both as a power law and, because the power law turned out to
describe it poorly, by interacting insulin action with glucose band indicators,
which imposes no shape.

Fitted as a power law rather than as bands, the glucose exponent comes to +0.238
with a standard error of 0.044. That value is not quoted in the body above because
it is unstable: restricting to a narrower glucose range moves it to about +1.0,
a shift larger than the number itself. A power law cannot describe a step, so
fitting one across the profile in the table averages the step into a slope of
almost nothing whose sign follows whichever range was chosen. The band table is
reported in its place.

An earlier draft of this document reported a glucose exponent of -2.67. That was
an indexing defect: the interaction term had been inserted ahead of the insulin
coefficient in the design matrix while the coefficient was still read from a fixed
position, so the interaction was divided by the starting glucose term. It was
found when a second implementation disagreed on the sign. Design matrix columns
are now addressed by name and a regression test covers it.

Code is in `inv009/` in the `dynamic-isf-calculations` repository, with 29 tests.

```
psql -d oref -f inv009/ingest/sql/05_settings.sql
python3 inv009/ingest/load_settings.py
python3 -m inv009.build_cache
python3 -m inv009.entered_isf && python3 -m inv009.effective_isf
python3 -m inv009.tdd_axis && python3 -m inv009.glucose_axis
python3 -m inv009.head_to_head && python3 -m inv009.synthetic
python3 -m inv009.loop_model_infer && python3 -m inv009.ml_shap
python3 -m inv009.carb_hypothesis && python3 -m inv009.joint_surface
python3 -m inv009.pointwise_isf && python3 -m inv009.figures
```

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
for Health Research, whose releases carry no attribution wording of their own. The
analyses, content and conclusions here are solely mine and have not been reviewed
or approved by any study sponsor.

Dates in these archives are shifted or rebased per participant, so time of day is
preserved and calendar date is not. Nothing here depends on calendar date.
