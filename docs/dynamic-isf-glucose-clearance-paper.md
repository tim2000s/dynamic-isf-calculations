# Does insulin sensitivity really fall with glucose, or does carbohydrate absorption outlast its model? An observational look at open-source closed-loop data

**Tim Street, with analysis support from Claude (Anthropic)** · 2026-06-10

> Audit update, 16 September 2026: the action-balance calculation described below is not
> model-independent. IOB is model-derived, suggested delivery may differ from enacted delivery,
> and the carbohydrate screen selects windows using their later glucose trajectory. The paper's
> predictive findings remain useful, while physiological claims are withdrawn. See
> `DYNAMIC-ISF-AUDIT-2026-09.md`.

---

## The short version

If you run AAPS, Trio or oref0, you've probably met Dynamic ISF (or its cousins). The idea is simple and biologically reasonable: when your glucose is high, insulin doesn't work as well, so the loop should expect a smaller drop per unit and dose a bit harder. That's a real effect in the body, it's well described in the literature, and plenty of people feel that Dynamic ISF helps them. I'm not setting out to argue with any of that.

What I wanted to know was narrower and more practical. When you actually measure the glucose drop we get per unit of insulin on our own looping data, the number a loop is really dosing against, does it fall as glucose rises the way the dynamic equations assume?

The honest answer from this data is: not in the way you'd expect. Once you take carbohydrate out of the picture, the effective sensitivity doesn't fall with glucose at all. The falling-sensitivity pattern that Dynamic ISF is built around shows up specifically around food: both the meals people don't tell the loop about, and, more tentatively, the tail end of the meals they do announce, where absorption may carry on a little after the loop has decided carbs-on-board is back to zero.

So my reading, and I do want to stress that it's a reading rather than a proof, is that a glucose-driven ISF may be doing a useful job for a reason other than the one we usually give it. It may be quietly compensating for carbohydrate that the absorption model hasn't fully accounted for, rather than correcting a minute-to-minute loss of insulin sensitivity. That's not a criticism of the algorithm or of anyone's carb model. It's a suggestion about what the algorithm is actually fixing, and where the most useful work might sit next.

The rest of this is how I got there, with the caveats, because the caveats matter.

---

## 1. What Dynamic ISF assumes, and why it's reasonable

The insulin sensitivity factor, how far your glucose falls per unit of correction insulin, sits underneath every correction a loop makes. A static profile uses fixed values by time of day. Dynamic ISF varies it, lowering the ISF (so the loop expects less drop and gives more insulin) as glucose climbs and as total daily dose rises.

The reasoning behind the glucose part is sound. Sustained high glucose does blunt insulin's effect, through mechanisms that are well characterised in the cell biology: glucotoxicity, the knock-on effects on the insulin signalling cascade, oxidative stress and the rest [1–7]. Nobody sensible disputes that high glucose and insulin resistance travel together over the long run.

The question is whether that long-run, chronic relationship is the right thing to bake into a correction decision that plays out over the next hour or two, and whether what the loop sees as "insulin working less well when I'm high" is actually that, or something else wearing the same costume.

## 2. How I looked at it

The core idea is to ask, for each correction, what ISF would have explained the glucose drop that actually happened:

> effective ISF = the glucose drop we actually got ÷ the insulin that was actually working

Then you can lay any candidate (a static value, the v1 form ISF ∝ 1/TDD, the v2 form ∝ 1/TDD², or a steep glucose curve) against that same drop, on the same correction, for the same person. Comparing like-for-like on one window strips out the biggest confounder, which is that different people simply have different sensitivities.

I worked out the "insulin that was actually working" two ways, deliberately, so the answer doesn't lean on any one assumption:

- **From the loop's own prediction.** This rescales the loop's insulin-on-board forecast. Convenient, but it inherits the loop's insulin-action curve.
- **From conservation, model-free.** This is the insulin that genuinely acted, measured as the fall in insulin-on-board plus whatever was delivered during the window (SMBs and basal). It needs no assumption about how insulin acts over time; it just accounts for the insulin. Over a few hours most of a fast-acting dose has worked, so what's left over is small and the estimate is sturdy.

Importantly, the ISF the loop happened to be running cancels out of both, so this isn't circular. We're not grading the loop against itself.

To keep carbohydrate out where I wanted it out, I used fasting windows and threw away any window where glucose was rising in a way that smelt of absorption. For the overnight work I leaned on the small hours, when even someone who never announces a carb is genuinely not eating.

And then, because cohort averages can hide as much as they show, I ran the same method on two real individuals' live Nightscout data, entirely outside the main dataset: one who runs AAPS fully closed without announcing carbs, and one on oref who announces carefully. They turned out to illustrate the two halves of the story rather neatly.

The main dataset is per-tick logs from people on open-source AID (predominantly oref0, with some Trio), held locally in a database for the analysis. Numbers below are medians unless I say otherwise.

## 3. A well-set static ISF already does the job

Before anything subtle, the blunt result. On identical correction windows, the typical size of the prediction error (mg/dL) came out as:

| | loop as-run | static profile ISF | v1 (1/TDD) | v2 (1/TDD²) |
|---|---|---|---|---|
| typical error | 18.6 | 20.3 | 24.6 | 49.7 |

A well-set static ISF essentially matches the loop and comfortably beats both dynamic forms. The v2 form, with the steepest dependence, is far and away the worst. And when you fit the relationship between ISF and total daily dose directly, it scales roughly as 1/√TDD, much gentler than either the 1/TDD or 1/TDD² the equations assume.

That on its own doesn't tell you why. It just says the steep equations are reaching for something that isn't paying off on average. The rest is about what that something is.

## 4. Take carbohydrate out, and the fall disappears

Here's the part that made me sit up. When you compute the effective ISF in genuinely fasting, carbohydrate-free windows and plot it against glucose, it does not fall. If anything it's lowest near target (which makes sense, because as you come down towards normal the body starts defending against a hypo and the drop flattens off), and then it's flat to gently rising as you go higher. That's the opposite shape to Dynamic ISF.

**Effective ISF against glucose (each divided by the person's own profile ISF, so 1.0 = their profile).**

| glucose (mg/dL) | effective ISF (loop-based) | effective ISF (conservation) | v1 | v2 |
|---|---|---|---|---|
| 100–120 | 0.71 | 0.20 | 0.93 | 2.89 |
| 120–145 | 0.82 | 0.37 | 0.87 | 2.35 |
| 145–175 | 0.91 | 0.47 | 0.72 | 1.57 |
| 175–205 | 0.89 | 0.53 | 0.60 | 1.08 |
| 205–260 | 0.94 | 0.45 | 0.50 | 0.78 |

The two effective-ISF estimates sit at different levels (the conservation one reads low, for a reason I'll own up to in the limitations), but they agree on the thing that matters: the line trends gently up, while the dynamic forms trend steeply down. Restricting to windows where glucose was dead flat on entry barely moves it.

![Two model-dependent outcome ratios against glucose, both trending gently upward, set against the steeply falling dynamic form.](charts/inv008/fig_effective_isf_independent.png)

I also looked at where the genuinely useful individual signal lives. When you let a model learn a per-person intercept and scale but share the glucose shape across everyone, the gain comes almost entirely from getting each person's overall level right (worth about 7 mg/dL), and the best shared glucose steepness lands at essentially zero. That per-person level isn't predictable from total daily dose, so it has to be learned from the person's own results over time. A careful out-of-sample test found that learning a personal glucose curve only beat a flat line for about a quarter of people, below what you'd get by chance, which tells you the handful who looked like they wanted a glucose curve were mostly noise dressed up as signal.

This held in the small hours and in the middle of the day alike (the slope was flat in both), so it isn't a quirk of overnight physiology.

## 5. The falling pattern lives around food, not high glucose itself

If the effective ISF doesn't fall when you're fasting, why does it so obviously fall in the raw, all-day picture, and why do people feel Dynamic ISF earning its keep? Two things, and the second one is the one I didn't expect.

### 5a. The meals you don't announce

Take the first external individual (call them User A), running AAPS fully closed-loop and not announcing carbohydrate, relying on the loop to catch meals itself. Across all hours, their effective ISF falls steeply with glucose, and a glucose-driven ISF predicts their drops much better than a static one (typical error around 28 mg/dL for the Dynamic ISF they run, versus 49 static). On the face of it, a textbook case for Dynamic ISF.

But the fall is a daytime thing. Restrict to the overnight sleep window, when even a non-announcer is genuinely not eating, and the steepness collapses to a fraction of its daytime value. Roughly three-quarters of the apparent glucose dependence simply goes away once they're asleep and not eating.

Read plainly: a lot of what looks like "insulin working worse when I'm high" for this person is unannounced food sitting in the system and holding glucose up, depressing the apparent drop-per-unit during the day. Their Dynamic ISF is doing something genuinely useful here, leaning in harder to cover meals the loop wasn't told about. But it's covering carbs, not correcting a sensitivity that's actually changed.

### 5b. The tail of the meals you do announce

This second pattern is less clear-cut, so I'll flag the uncertainty as I go. I include it because it points in the same direction as the rest, not because I think it's settled.

Take the second external individual, User B, on oref, a careful carb-announcer. You might assume that once their carbs-on-board has counted down to zero, the meal is done and any window after that is clean. I looked at whether that holds, comparing the effective ISF in the hour straight after carbs-on-board reached zero against genuinely carbohydrate-free overnight windows for the same person.

The hour right after carbs reach zero does look a little different from true overnight fasting. The effective ISF in that hour tends to run lower, glucose is somewhat more likely to be drifting up even with active insulin on board, and the value climbs back towards the fasting level over the following half hour. The same general pattern showed up when I ran the test across the main cohort, one person at a time.

![After announced carbs reach zero, effective ISF tends to sit lower and recover over the following hour. Left: per person across the cohort. Right: the same person-by-person recovery.](charts/inv008/fig_post_cob_isf.png)

![The same test on User B's five months of data shows a similar shape.](charts/inv008/fig_post_cob_isf_dnzxy.png)

I'd be cautious about over-reading this, and there's a real confound to put on the table first. The hour after a meal is a busy one. There's more insulin still in flight, and some of its effect lands later than the glucose you can see inside the window, so part of that gap will be timing rather than carbohydrate. With that allowed for, timing alone wouldn't obviously push glucose up against active insulin, and it doesn't obviously explain why the effect eases off as the meal recedes. One explanation that's consistent with both the announcers and the non-announcers is that a little absorption can carry on past the point where the carbs-on-board count has finished. I'd put that no more firmly than as a possibility worth looking at properly, rather than something this analysis can confirm.

### What the cohort says when you pool it

Stepping back to the whole dataset and looking only at clean, carbohydrate-free overnight windows: the effective ISF sits a touch below each person's profile on average, but that offset doesn't track how much they announce carbs, so it reads as ordinary calibration (profiles set a little strong) rather than anything to do with food or glucose. And the glucose slope overnight is flat, not falling. Neither the level nor the slope, in clean fasting, looks like glucose-driven resistance.

![Across the cohort, in carbohydrate-free overnight windows: the small sub-profile offset doesn't track carb-announcing behaviour (it's calibration), and the glucose slope is flat, not the falling shape Dynamic ISF assumes.](charts/inv008/fig_cob_uam_cohort.png)

## 6. One more thing that exaggerates the picture

There's a second, more mundane reason the raw data looks like falling sensitivity, and it's worth naming so it isn't mistaken for physiology. The loop tends to over-predict big drops: when it expects a large fall, it gets roughly half of it. Across the data the realised drop works out at around 0.4 of the predicted drop plus a fairly fixed offset. Because the biggest corrections happen when you're highest, the size of the correction and the glucose level move together, and the loop's habit of over-promising large drops then masquerades as "insulin doesn't work when I'm high." That's really a property of the insulin-action model and basal, not of ISF, and it shouldn't be fixed by bending the ISF curve.

## 7. What I think this might mean

Putting it together, and keeping the hedges where they belong.

The effect Dynamic ISF was built to capture, glucose falling more slowly per unit of insulin when you're high, is real in the raw data. Our data show it too. What I'm questioning is the cause we usually attribute it to. When you remove carbohydrate, the effect largely goes with it. It reappears around unannounced meals, and possibly around the tail of announced ones. So one reasonable reading, on this data and this timescale, is that some of what a glucose-driven ISF does may be compensating for carbohydrate the loop hasn't fully accounted for, rather than tracking a minute-to-minute change in insulin sensitivity. I offer that as a candidate explanation, not a settled one.

I'd put it no more strongly than that. This is observational data, the effective-ISF estimates depend on a model, and I've only looked at the hour-to-hour correction timescale rather than weeks and months. I can't exclude a real physiological component, and chronic glucose-related resistance plainly exists. It's just that, on this timescale, it seems small, and where it does show up it shows up near target, in the opposite direction to what Dynamic ISF assumes. The chronic part is probably already sitting quietly inside each person's overall sensitivity level.

None of this means Dynamic ISF doesn't help people. Dosing harder when you're high and drifting will, often enough, bring you down, and if the reason you're high and drifting is uncovered carbohydrate, then a glucose-driven nudge is a reasonable, if indirect, way to catch it. The point is only that part of what it's catching may be carbs, and that's useful to know when you think about where to improve things next.

A word on the absorption model, because this is the part I most want to get right. The oref dynamic carbohydrate-absorption model was a genuine step forward. It already adapts absorption to what glucose is actually doing rather than assuming a fixed curve, and nothing here says it's wrong. The tail of a meal is one of the hardest things to recover from CGM alone, and it's very individual, so some residual beyond the modelled carbs-on-board is exactly what you'd expect from any model of this kind. At most, what this data gently hints at is that the loop may currently read part of that residual as reduced sensitivity. Whether better tail estimation, or better detection of carbohydrate the loop wasn't told about, could take on some of the work a glucose-driven ISF does today is, to me, a fair question for future work. I raise it as a question, nothing more.

## 8. So what would I actually change?

Modestly, and in this order:

- **Get each person's overall sensitivity level right, and keep adapting it from their own results.** This is where the real, recoverable gain sits (about 7 mg/dL), it isn't predictable from TDD, and it dwarfs anything a glucose curve adds. A 1/√TDD starting point for a cold start, refined online, is a sensible default.
- **Keep a near-target easing of ISF as a one-sided safety measure.** Not because the data is fitted to it, but because counterregulation near target is real and the cost of a hypo is asymmetric. This is the one place real physiology pokes through, and it argues for being more cautious near target, not less.
- **Treat carbohydrate as carbohydrate.** For people who don't announce, better catching of unannounced meals does the job a glucose-ISF is standing in for. For people who do, any residual is most likely to sit in the meal tail. Both feel like more honest levers than a glucose-indexed sensitivity.
- **If you keep a glucose term, understand what it's for.** It's a reasonable, crude proxy for carbohydrate the loop is missing, and it's most valuable for non-announcers. It's worth remembering that when carbs are announced, a glucose-driven ISF and the carb model can end up leaning on the same meal from two directions, so it shouldn't be left to run hot by default.

## 9. Limitations, plainly

This is observational data from people living their lives, not a controlled experiment, and it can only weakly pin down cause. A direct dose-to-drop test is no use here, because the loop doses reactively: it gives more correction exactly when glucose isn't budging, so more insulin can look like less effect. That's why I leaned on the loop's counterfactual prediction and the model-free conservation estimate, neither of which regresses outcome on dose.

The conservation estimate is mostly from oref0 users (I didn't assemble Trio basal profiles), and its absolute level reads low because the basal-deviation accounting over-counts a bit where someone's profile basal is set wrong, so I read its direction, which is sound, not its absolute numbers. The post-carb-clear effect is partly inflated by post-meal insulin still in flight, as noted above. The two external individuals are a person each, with their own quirks, and the AAPS one is noisy in the top glucose band overnight. And the recommendation, particularly anything touching the absorption tail, should be tried in shadow mode and watched before it goes anywhere near live dosing.

## 10. The resistance question stays open

I want to be clear that I haven't closed the biological question. An earlier draft of this work tried to resolve it, arguing that real high-glucose insulin resistance was being exactly cancelled by glucose-rising, insulin-independent clearance, leaving a flat net. It's a tidy idea and physiologically plausible, but when I went looking for the data to test it, it wasn't there: the windows you'd need (high glucose, essentially no insulin acting, no carbohydrate) barely exist, and the few that do are increasingly contaminated by suspected uncovered carbohydrate as glucose climbs. The estimate flipped sign depending on how I selected the windows, so I pulled the claim. I mention it because it's the right thing to do, and because it's a reminder of how thin the ground gets once you push past what this data can carry.

The practical conclusion doesn't depend on settling it, though. Whatever is going on underneath, the number a loop doses against, the realised net drop per unit, doesn't fall with glucose once carbohydrate is out of the way. That's the bit that matters for how we set up our loops.

## 11. Where that leaves us

For most of us, most of the time, a well-set, individually-adapted static ISF with a bit of near-target caution will do as well as a glucose curve, and it's easier to reason about. Dynamic ISF isn't doing nothing. But on this evidence, a fair part of what it's doing may be mopping up carbohydrate, announced and unannounced, that the loop hasn't fully accounted for, rather than tracking a sensitivity that genuinely changes minute to minute. If that's right, then the most useful place to push next probably isn't a steeper ISF curve. It's getting each person's baseline sensitivity right, and getting better at carbohydrate, which is, after all, the thing most of us were fighting in the first place.

I'd welcome people poking holes in this. The code's there to do it with.

---

## References

1. Rabbani N, Thornalley P (2024). *Front Endocrinol.* doi:10.3389/fendo.2023.1268308
2. Khalid M, Alkaabi J, Khan M, et al. (2021). *Int J Mol Sci.* doi:10.3390/ijms22168590
3. Zhao X, An X, Yang C, et al. (2023). *Front Endocrinol.* doi:10.3389/fendo.2023.1149239
4. Simon-Szabó L, Lizák B, Sturm G, et al. (2024). *Int J Mol Sci.* doi:10.3390/ijms25169113
5. Galicia-Garcia U, Benito-Vicente A, Jebari S, et al. (2020). *Int J Mol Sci.* doi:10.3390/ijms21176275
6. Beaupere C, Liboz A, Fève B, et al. (2021). *Int J Mol Sci.* doi:10.3390/ijms22020623
7. Młynarska E, Czarnik W, Dzięca N, et al. (2025). *Int J Mol Sci.* doi:10.3390/ijms26031094

## Reproducibility

Code: `github.com/tim2000s/dynamic-isf-calculations`, package `inv008/` (`python -m inv008.<name>`).
Key stages: `head_to_head` and `effective_isf_independent` (the two effective-ISF estimates),
`cob_uam_cohort` (the carbohydrate-free overnight decomposition), `post_cob_isf` and
`post_cob_isf_dnzxy` (the hour after carbs reach zero, cohort and external individual), `magnitude_bias`
(correction-size effect), `gradient_isf_fit` and `adaptive_k_nestedcv` (per-user level versus glucose
curve), `daytime_clearance` (daytime versus overnight).
