# v1 vs v2 (UPDATED) comparison

Updated v2: `ISF = 115000/(TDD²·ln(BG_floored/divisor))` (anchor `2300/(ln(target/divisor)·TDD²·0.02)`, no +1; BG floored at divisor+1). 170 users, 9,517,761 ticks (cached replay).

## Headline

- updated v2 gives **3.04× the ISF of v1** (median) — i.e. much weaker corrections; weaker than v1 on **92%** of ticks.
- updated v2 is **1.96× the OLD v2** (median) — the dropped +1 raises ISF substantially.
- the ratio is now **BG-dependent** (old v2 was a flat 63.9/TDD):

| BG band | n | median ISF_v2updated/ISF_v1 |
|---|---|---|
| (40, 80] | 778,636 | 53.11 |
| (80, 100] | 1,671,653 | 6.24 |
| (100, 120] | 2,063,028 | 3.5 |
| (120, 150] | 2,265,314 | 2.53 |
| (150, 200] | 1,784,604 | 1.89 |
| (200, 360] | 906,116 | 1.48 |

## Reading

- The dropped +1 makes ISF blow up as BG approaches the divisor floor → **very high ISF (near-zero correction) at low BG** — strong hypo protection — tapering toward ~1.6× v1 at high BG.
- Net effect vs v1: updated v2 is *less* aggressive everywhere, dramatically so below ~100 mg/dL. Versus the old v2 it is uniformly higher-ISF (the +1 removal ≈3× at target).

*Caveat: counterfactual replay on plain-oref cohort; v1 unchanged (1800/(TDD·ln(BG/div+1))); high cap 210 retained, low floor divisor+1 added.*