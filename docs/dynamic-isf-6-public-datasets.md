---
title: "Does your sensitivity really halve when your dose doubles?"
subtitle: "Testing both dynamic ISF equations against 1,684 people in seven public study archives"
author: "Tim Street"
date: "25 August 2026"
---

## Where this starts

If you loop, there's a number in your settings that decides how hard the system
corrects you: how far one unit moves your glucose. Most of us set it once, from a
rule of thumb, and then nudge it when something feels off.

The dynamic ISF equations take a different view. They say that number shouldn't
be fixed at all, because sensitivity depends on how much insulin you use, and
they recalculate it from your recent total daily dose several times an hour.

There are two versions. The original, which I'll call v1 throughout, has
sensitivity falling in step with your dose: double the dose, halve the
sensitivity. The revised maths, v2, has it falling with the square of the dose:
double the dose and sensitivity drops to a quarter. That's not a small
difference. For two people whose doses differ threefold, v1 says one is three
times more sensitive and v2 says nine times. The two cross at around 64 units a
day, so whether the newer version makes your corrections stronger or weaker
depends on which side of that you sit.

I've looked at this before, on data from people running loops, and came away
thinking the relationship was real but shallower than either version assumed. I
couldn't settle it. The group was a few hundred people, nearly all on doses
between 20 and 80 units a day, so there wasn't much to see either side of the
crossover. And every estimate I made leaned on the loop's own insulin model,
which is uncomfortably close to using the thing you're questioning to answer the
question.

The public study archives fix that. What follows uses seven of them: about 1,700
people, ages 2 to 82, daily doses from 7 to 107 units, four different systems,
and one group who weren't using any system at all.

I should say at the outset that this isn't an attack on anyone's work. These
equations came out of the same rule of thumb most of us set our pumps from, built
by people giving their time to make looping better, and one of the things I found
is that the rule holds up remarkably well. What I couldn't find was support for
the squared version.

## What I found

Sensitivity does fall as daily dose rises. That much is solid, and it shows up in
four different ways, including one that was free to find no relationship at all.

How steeply is the interesting part. Measured from what people had entered in
their pumps, it's almost exactly the 1800 rule, which is v1's assumption.
Measured from how glucose actually responded overnight, it's a bit shallower than
that. Nothing I looked at came anywhere near the squared law.

The other half of both equations, the bit that changes sensitivity according to
where your glucose is right now, didn't earn its place. Six different shapes
including the ones both equations use came out within 0.05 mg/dL of each other,
and the flat one was best. That matches what I found on looping data last year,
so I'll not go over the same ground again here.

And when I asked each version to predict how far glucose would actually fall
overnight, a plain static 1800 rule beat all of them.

## The archives

| Study | What people were using | People | Carbs logged | Median days | Median dose |
|---|---|---|---|---|---|
| Loop | Do-it-yourself Loop, fixed sensitivity | 842 | yes | 397 | 38.6 U |
| REPLACE-BG | **A pump and a sensor, nothing automating** | 196 | yes | 246 | 41.9 U |
| DCLP3 | Control-IQ, adults and teenagers | 112 | no | 183 | 50.0 U |
| DCLP5 | Control-IQ, ages 6 to 13 | 100 | no | 203 | 37.4 U |
| PEDAP | Control-IQ, ages 2 to 5 | 98 | no | 273 | 13.6 U |
| IOBP2 | Bionic pancreas | 336 | no | 93 | 49.0 U |

FLAIR is in the archive too, but its release ships no pump file, so it can tell
us about daily totals and nothing else.

Two rows there do most of the work later on. REPLACE-BG ran in 2015 on a pump and
a sensor with nobody's algorithm in between, which makes it the only group where
the insulin someone got wasn't chosen by a machine watching their glucose. And
Loop and REPLACE-BG are the only two where people logged what they ate, so
they're the only two where I can genuinely screen out a meal rather than guess.

That comes to about 2.2 million candidate overnight windows, of which 708,106
survive screening, from 1,659 people with at least forty each.

## Working out sensitivity when nobody wrote it down

None of these archives holds an insulin on board figure, or a loop's prediction,
or in most cases even a sensitivity setting. It all has to be rebuilt from the
basal and boluses that were delivered.

The idea is simple enough. Take four hours overnight, starting on the hour
somewhere between 11pm and 3am. Note how far glucose fell. Work out how much
insulin was actually acting across those four hours by running the delivery
record through an insulin model. Loop users get the models Loop itself shipped,
lifted from the LoopKit source rather than approximated; everyone else gets the
oref curve. The fall divided by the insulin that acted, allowing for where
glucose started and where it had been, gives milligrams per decilitre per unit,
which is what a sensitivity factor is.

Two details in there matter more than they look.

**I keep insulin separate depending on when it was decided.** Insulin already in
you when the four hours start was committed before anything in that window
happened, so it can't be a reaction to the glucose I'm about to measure. Insulin
delivered during the window can be, and if a system is running it usually is.
Those go in separate columns and stay there.

**The screen only looks backwards.** A window is kept or dropped on where glucose
started, on how complete the sensor trace is, and on how long since the person
last ate. It's never dropped on what glucose then did. It's tempting to bin the
nights where glucose rose, on the grounds that someone must have eaten, and it's
a mistake: those are exactly the nights insulin worked least well. Dropping them
doesn't clean the answer up, it flatters it.

### Checking the arithmetic against something real

Loop wrote down its own insulin on board figure at every bolus, and that's the
one place in all of this where an app's internal state got recorded. So I can
rebuild it from the pump record and see whether I get what the app got.

Across 158 people with enough of those records, my reconstruction tracks Loop's
own figure at a **median correlation of 0.927**, with a median difference of 0.32
units. That's the whole dose reconstruction checked against something other than
my own assumptions, which is not a luxury this kind of work usually gets.

Getting there needed something the archives don't obviously contain. Loop counts
insulin on board net of your scheduled basal, so a temp basal only contributes
the difference between what it gave and what your profile would have given
anyway. Rebuilding that needs the basal profile, which isn't in any settings file
in the archive. It turns out to be recoverable: every temp basal record carries
the rate it was overriding. Six and a half gigabytes of raw text boils down to
1.6 million rate changes, and gives basal schedules for 842 people.

I'll be straight about the limit of that check. Every candidate insulin model
comes out with a small positive bias, so the fastest-decaying one wins by
soaking it up, which means I can't really say which model each person was
running. The correlation is the finding. The model choice isn't, and nothing
downstream leans on it.

## What people had actually entered

The narrowest question first, and the cleanest one. Of the 596 people with a
sensitivity factor on record, 581 also have thirty or more complete days of pump
data to take a daily dose from.

Their settings track dose almost exactly as the 1800 rule says they should. The
slope is **-0.979**, with a 95% interval from -1.030 to -0.934. Multiply each
person's sensitivity by their daily dose and the median comes to 2067, against
the 2139 that v1 works out to at a normal target. If the rule holds, that product
shouldn't drift with dose, and it doesn't: the correlation between them is
**+0.024**. Do the same for the squared law and you get **+0.859**, which is what
it looks like when a rule doesn't hold, consistently, right across the range.

![Entered sensitivity against daily dose](../charts/inv009/fig_entered_scatter.png)

This needs a caveat said out loud rather than tucked away, because it's close to
circular. What you've entered isn't a measurement, it's a decision, and you and
your team probably arrived at it through the 1800 rule in the first place. That's
the same rule v1's maths came from. So finding that people sit on 1800 over dose
is partly a finding about how much grip that rule has on practice.

It isn't circular in the other direction, though. Nothing in clinical practice
nudges anyone towards a squared law, and nothing in what people run looks like
one.

Two smaller things worth keeping. Within a single age band the slope is
shallower, somewhere between -0.68 and -1.10, so a chunk of that tidy -0.98 comes
from pooling toddlers with adults. And carb ratios track dose at -0.63 with a
median product of 411, so the 500 rule is a looser fit to practice than the 1800
rule is.

## What the glucose actually did

Now the harder question. Not what people entered, but what a unit did.

![Every exponent I measured](../charts/inv009/fig_exponents.png)

Comparing between people, using only the insulin that was already committed, the
slope is **-0.83**, with a 95% interval of -1.23 to -0.36. Taking logs means
dropping anyone whose estimate came out negative, and those aren't a random
subset, so I also fitted the power law directly on the natural scale where
everybody survives. That gives -0.75, which is close enough to be reassuring.

Within a single person it's **-0.645**, with a standard error of 0.033, pooled
across 1,313 of them, and about three quarters share the direction. This is the version that's never been
tested, and it's the one that matters most in practice. A pattern that holds
across a population doesn't have to hold inside one person, and it's inside one
person that these equations run, several times an hour, off a dose figure blended
from your last few hours and days.

Both sit between the square root and v1. Neither is close to v2.

### A second opinion that assumes nothing

Everything above comes from fitting a shape and reading off a number, which makes
the answer only as good as the shape I picked. So I gave the whole thing to a
gradient boosted model instead, which assumes nothing. With 687,067 windows from
1,660 people it's free to decide sensitivity falls with dose, rises with it,
moves in steps, or doesn't depend on dose at all.

The thing to read off it is how insulin action and daily dose interact, not how
important dose is on its own. Dose has a big effect on how far glucose falls
overnight that has nothing to do with sensitivity, because people on more insulin
are different people. Sensitivity is the multiplier on insulin, so that's where
it lives.

| Daily dose | Shift in sensitivity attributable to dose |
|---|---|
| 15 U | +0.28 mg/dL per unit |
| 25 U | +0.06 |
| 34 U | +0.06 |
| 43 U | +0.02 |
| 55 U | 0.00 |
| 77 U | -0.04 |

It falls, steadily, without being asked to. And it flattens off above about forty
units a day, which is the shape a shallow relationship makes rather than a steep
one. I wouldn't read the sizes as sensitivity itself; they're shifts around the
model's own baseline and they carry the same shrinkage as everything else here.

### Why one small group matters more than its size

There's a reason to be suspicious of all of this, and it's better shown than
argued.

If a system is running, insulin is close to a function of your recent glucose.
You get more of it precisely when glucose isn't coming down. Regress the fall on
all the insulin that acted and the answer goes negative, which would mean insulin
puts glucose up. It doesn't. What you've measured is the algorithm's policy, not
anybody's physiology.

![The confound, made visible](../charts/inv009/fig_endogeneity.png)

Read across those groups and the share of people whose sensitivity comes out
positive collapses from 78% to 8% as the system gets more reactive, then recovers
once I use only the insulin that was already committed. REPLACE-BG behaves the
way the design predicts: it's the only group where counting all the insulin makes
the answer better rather than worse, because in 2015 there was no algorithm
choosing it.

That's why 196 people from a decade ago carry more weight than their number
suggests, and why I report them separately throughout. Their own within-person
slope is **-1.38**, steeper than the pooled figure and the closest anything here
gets to v1. It isn't close to v2 either.

## Does sensitivity change with where your glucose is?

Both equations say it does, through a term that makes a unit do less when you're
running high.

There's a trap in testing that. High glucose falls further than low glucose
whatever insulin is doing, through mass action, through what the kidneys clear,
and through plain regression to the mean. Near target your body actively pushes
back against a further fall. So if you let glucose explain the size of the drop,
and then ask whether the drop per unit also depends on glucose, you'll answer the
first question and write it down as the second. Every fit here carries a glucose
term of its own for that reason, and what carries the claim is the interaction
between insulin and glucose, because that's the only thing that means sensitivity
itself moved.

Two answers came back.

Pooled across 1,528 people the exponent comes out at **-2.67**, and the sign is
the opposite of what both equations expect: net sensitivity goes up with glucose,
not down. Only about 40% of people
show it going the other way. That agrees with what I found on looping data last
year, and with the reasoning there, so I won't rehearse it beyond noting that
whatever the equations are reaching for is competing with clearance that rises
with glucose and with your body defending you near target.

The more useful answer is the practical one.

![No glucose shape earns its place](../charts/inv009/fig_glucose.png)

Out of sample, across 1,616 people, six candidate shapes land within 0.05 mg/dL
of each other and the flat one wins. Both equations' curves are in the losing
half. Before running any of this I'd set the bar at half a milligram per
decilitre before a glucose term could be said to earn its keep. Nothing came
within a tenth of that.

## Predicting the night

The slopes are the science. This is the bit that matters if you're actually using
one of these: if you'd used each version's sensitivity to predict how far glucose
would fall overnight, how wrong would you have been?

Each candidate gets its own offset for each person, fitted on that person's first
70% of nights and scored on the rest. That's generous to the equations rather
than harsh on them. Overnight most of your insulin is basal, and it's there to
offset the glucose your liver is putting out, so charging that to the sensitivity
factor would be blaming an equation for something it was never asked to do.

![Predicting the overnight fall](../charts/inv009/fig_head_to_head.png)

| Candidate | Median error | Median bias |
|---|---|---|
| Static 1800 rule | **34.1 mg/dL** | +2.9 |
| Best single number for that person | 34.4 | -1.4 |
| v1 | 36.7 | +0.9 |
| Their own entered setting | 37.3 | -0.2 |
| 355 over the square root of dose | 37.9 | +3.8 |
| v2 | **140.9** | +41.1 |

Out of 1,626 people, the best option is a single well-chosen number for 638, the
static 1800 rule for 587, the square root form for 195, v1 for 135, and their own
setting for 71. There's nobody it's v2 for.

Where v2 struggles most is where I'd least want it to. For people on under twenty
units a day its median error is 322 mg/dL, because squaring the dose hands a
small child a sensitivity of several thousand and therefore a correction dose of
almost nothing. The error comes down as dose goes up, to 83 mg/dL above 64 units
a day, which is the crossover doing exactly what the algebra says it should. Even
there it's more than twice as far out as anything else on the list.

## What would have shown up if the steeper version were right

Everything above estimates something nobody observed, using a reconstruction of
insulin that's certainly imperfect. The sensitivities I get are well below what
people have entered, so something is shrinking them. The question that matters is
whether that shrinkage also bends the slope. The way to find out is to build
people whose sensitivity genuinely follows a law I chose, and run them through
the identical machinery.

![Recovering a known relationship](../charts/inv009/fig_synthetic.png)

| The truth I built in | What came back, no system | With a reactive system | With 45% of meals never logged |
|---|---|---|---|
| -0.50 | -0.50 | -0.42 | -0.72 |
| -1.00 | **-1.02** | -0.80 | -1.03 |
| -2.00 | **-1.57** | -1.85 | -1.56 |

Shallow relationships come back almost exactly. Steep ones come back a bit
flattened. The row that matters is the bottom one: a squared law comes back
reading between -1.56 and -1.85 under everything I tried, including a run where
nearly half of all meals go unlogged and there's a third of a unit of error on
what actually acted. What I measured was -0.65 to -0.88. I can't get a squared
law to look like that.

The shrinkage in the simulation runs between 0.38 and 0.63 of the truth, which is
less than the real gap against what people have entered. So the level isn't
something this method recovers well, and I'm not claiming it does. The shape is.

## What this doesn't settle

**It doesn't give you a number to set.** The sensitivities I fit are shrunk, so
you can't read a constant off this work, only the shape of the relationship.

**The groups using a system are weak evidence on their own.** In the Control-IQ
studies the share of people with a positive sensitivity is barely better than a
coin toss, and for the youngest group it's worse than one. They replicate a
direction. They don't establish it. The open-loop group and what people had
entered are carrying the argument.

**All of this is overnight and fasting.** It says nothing about the hours after a
meal, and that's exactly where a glucose term might still earn its place. It's
also where carbs the absorption model hasn't fully caught remain a perfectly good
alternative explanation for anything that looks like sensitivity moving. That
seems to me the more interesting open question, and it isn't one this work can
answer.

**Sensitivity that moves with the day rather than the dose is out of scope.**
Exercise, illness, a change in how you eat: all of them move sensitivity over
hours and days. Nothing here measures any of that, and nothing here says it isn't
real.

**Two faults in the underlying data turned up**, and anyone else working from
these archives will want to know. Extended bolus durations arrive in milliseconds
for three of the six studies and in seconds for the other two, so a one hour
square wave lands as a thousand hours; spreading a bolus across that put insulin
into 41,591 of one person's 46,176 five minute periods. And the database hands
back timestamps at microsecond resolution, which quietly divided every basal
total by a thousand until I checked my reconstruction against the database's own
sums. Both are fixed here at the point of reading rather than in the data itself.

## Where that leaves the two versions

v1's assumption about dose is close to right, and it's closest of all in the
place it came from: it matches what people actually run their pumps at, almost
exactly. Measured against how glucose behaves it looks a touch steep, with the
honest range somewhere between the square root and v1, but it's the right shape
and the right size. Its glucose term doesn't earn its place, though it costs
little either way.

For v2 I couldn't find support in any of this. What people run contradicts it,
the glucose response contradicts it, it wasn't the best predictor for a single
one of 1,626 people, and the simulation says a squared law couldn't have hidden
as anything I saw. What concerns me most is where it struggles hardest, which is
in children and in anyone on a small daily dose, the people least able to absorb
a correction that turns out to be far too small.

I want to be careful about what that does and doesn't mean. It's a statement
about what a large and varied group of people did, not about anybody's intent or
competence, and the work behind these equations has moved the conversation
forward in ways this analysis is standing on. It's also one analysis, overnight
and fasting, from archives that were never collected for this purpose. I'd be
glad to see someone come at it differently.

If a dose term is wanted, what's here supports something between about -0.65 and
-1.0, applied to a person's own level rather than used to set it. But given that
a plain static 1800 rule outpredicted every dynamic form I tested, the question I
keep coming back to isn't which exponent to use. It's whether the sensitivity
factor is the right place to be putting any of this, and whether the effort would
land better on the carbohydrate side, where the residual is bigger and nobody
thinks the problem is solved.

## Method, and how to repeat it

The code is in `inv009/` in the `dynamic-isf-calculations` repository, with 28
tests. The window cache rebuilds in about nine minutes on a laptop.

```
psql -d oref -f inv009/ingest/sql/05_settings.sql
python3 inv009/ingest/load_settings.py
python3 -m inv009.build_cache
python3 -m inv009.entered_isf && python3 -m inv009.effective_isf
python3 -m inv009.tdd_axis && python3 -m inv009.glucose_axis
python3 -m inv009.head_to_head && python3 -m inv009.synthetic
python3 -m inv009.loop_model_infer && python3 -m inv009.ml_shap
python3 -m inv009.figures
```

Windows are four hours from starts between 11pm and 3am, needing 80% sensor
coverage, a starting glucose between 90 and 300 mg/dL, and either six hours clear
of logged carbohydrate or, where none is logged, four hours clear of a bolus that
is large relative to that person's own daily dose. Per-person fits are ordinary
least squares with a heteroskedasticity-robust standard error; population figures
pool the per-person estimates using DerSimonian and Laird, as in the earlier
work, so the two can be read side by side.

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
for Health Research, whose releases carry no attribution wording of their own.
The analyses, content and conclusions here are solely mine and have not been
reviewed or approved by any study sponsor.

Dates in these archives are shifted or rebased per participant, so time of day is
intact and calendar date isn't. Nothing here depends on calendar date.
