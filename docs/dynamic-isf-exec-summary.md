# Executive summary — Does ISF need to change with glucose?

**Tim Street, with Claude (Anthropic) · June 2026** — *Dynamic ISF in open-source automated insulin delivery (oref0 / Trio)*

**The question.** Dynamic-ISF algorithms (and the Diabeloop model) lower the insulin sensitivity factor — the glucose fall expected per unit of insulin — as glucose rises, assuming high glucose makes insulin less effective. In practice they rarely beat a well-set static ISF. We tested why, on real closed-loop data.

## What we saw

- A well-set **static ISF predicts the glucose drop as well as the loop itself**, and far better than the dynamic v1/v2 equations.
- The **effective ISF does not fall as glucose rises**. It is *lowest near target* and flat-to-higher when high — the opposite shape to dynamic ISF. Confirmed two independent ways, including one that makes no assumption about how insulin acts over time.
- The accuracy gain comes from **personalising each user's overall sensitivity level** (~7 mg/dL), *not* from adding a glucose curve.
- **Two real users, tested directly:** a carb-announcing user was flat (matching the cohort); a UAM user (no carb entries) fell steeply — **but that fall mostly disappeared overnight**, when they were genuinely not eating.

## Why we saw it

- **The "insulin works worse when high" effect is mostly carbohydrate, not physiology.** High glucose is usually carb-driven. For carb-announcing users the carb model handles it, so the fasting ISF is flat. For UAM users, *unannounced carbs leak in and look like resistance* — which is exactly why their curve falls by day and flattens overnight.
- **A second confound inflates it:** the loop over-trusts large corrections (~2×), and corrections are bigger at high glucose — so correction *size* masquerades as a glucose effect.
- **Real glucose physiology does appear — but near target** (the body defending against hypos), the *opposite* direction to what dynamic ISF assumes.

## What it means

> **Use a per-user-adapted static ISF (√TDD starting point, refined online from the person's own outcomes) + a near-target easing clamp. Don't add a glucose curve by default.**

A glucose-lowering ISF mainly compensates for *unannounced carbs* — genuinely useful for UAM users, but that is carb-handling, not resistance-correction. So: handle carbs as carbs, individualise the level, keep the static backbone.

**Honesty note.** Glucose really does cause insulin resistance physiologically — but on the minute-to-hour fasting timescale of closed-loop correction it is too small to see, and the chronic part already lives in the per-user level. An earlier draft claimed a "clearance cancels resistance" resolution; it was withdrawn after audit as untestable.

**Scale.** 119 individuals · ~9.6 M loop decisions · ~62,700 overnight + 64,300 daytime carb-screened windows · >560,000 ISF evaluations · plus 2 external Nightscout users (12 and 5 months). *Full paper and code: `github.com/tim2000s/dynamic-isf-calculations`.*
