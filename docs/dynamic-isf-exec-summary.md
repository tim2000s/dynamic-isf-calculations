# Does insulin sensitivity fall with glucose, or is it carbohydrate? The short read

**Tim Street, with analysis support from Claude (Anthropic) · June 2026.** *Dynamic ISF in open-source AID (oref0 / Trio)*

> Audit update, 16 September 2026: the outcome analyses support a predictive comparison and
> residual carbohydrate as a plausible explanation. They do not independently measure
> physiological ISF. The phrase "two independent ways" in the original summary was incorrect
> because both estimators use model-derived IOB or loop predictions. See
> `DYNAMIC-ISF-AUDIT-2026-09.md`.

Dynamic ISF lowers your insulin sensitivity factor as glucose rises, on the reasonable assumption that high glucose makes insulin work less well. It's a real biological effect, and plenty of people feel Dynamic ISF helps. I wanted to check something narrower: when you measure the glucose drop we actually get per unit of insulin, the number a loop really doses against, does that fall as glucose rises? On our own looping data, the answer is more interesting than a straight yes or no.

## What I found

- **A well-set static ISF already matches the loop prediction closely** and beats the dynamic equations on the selected windows. The earlier square-root TDD recommendation is withdrawn because its estimator did not identify physiological ISF.
- **In the selected overnight windows, apparent sensitivity did not fall with glucose.** Both outcome proxies were flat-to-rising with glucose, opposite to Dynamic ISF. Neither proxy is physiologically independent: one rescales the loop prediction and the other uses model-derived IOB plus suggested delivery.
- **The falling pattern lives around food, not high glucose itself.** For someone who doesn't announce carbs, the fall is a daytime thing that mostly vanishes overnight when they're genuinely not eating. And, more tentatively, even for a careful carb-announcer, the hour after carbs-on-board reaches zero tends to show a lower effective sensitivity that recovers over the following half hour, with glucose more likely to drift up against active insulin. Some absorption may run on a little after the loop has called the meal over. Two separate datasets show a similar shape.
- **The real, useful individual signal is each person's overall sensitivity level** (worth about 7 mg/dL), which has to be learned from their own results, not a glucose curve, whose best-fit steepness comes out at zero.

## What I think it might mean (carefully)

The effect Dynamic ISF was built for is real in the raw data, and our data show it too. What I'm questioning is the cause. One reasonable reading here is that some of what a glucose-driven ISF does may be compensating for carbohydrate the absorption model hasn't fully caught (unannounced meals, and possibly the tail of announced ones), rather than correcting a minute-to-minute loss of sensitivity. I offer it as a candidate explanation, not a settled one. It's not a knock on Dynamic ISF or on anyone's carb model. The meal tail is genuinely one of the hardest things to estimate from CGM, and some residual is expected of any such model.

Controlled studies show that insulin sensitivity varies over time. These observational data cannot
separate that variation from glucose effectiveness, controller delivery and residual carbohydrate.

## What I'd actually do

> **Get each person's sensitivity level right and keep adapting it from their own outcomes; keep a near-target easing of ISF as a hypo-safety measure; treat carbohydrate as carbohydrate. Don't reach for a steeper glucose curve by default.**

If you keep a glucose term, it's best understood as a rough proxy for carbohydrate the loop is missing, most useful for people who don't announce. It's worth remembering that when carbs are announced it can end up leaning on the same meal as the carb model.

**Scale.** ~115 people · millions of loop decisions · tens of thousands of carbohydrate-free overnight and post-meal windows · plus two external individuals (12 and 5 months) who, between them, show both halves of the carbohydrate story. *Full write-up and code: `github.com/tim2000s/dynamic-isf-calculations`.*
