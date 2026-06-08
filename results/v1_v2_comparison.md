# v1 vs v2 comparison

v1 `ISF = 1800/(TDD·ln(BG/75+1))`; v2 `ISF = 115000/(TDD²·ln(BG_floored/75))`. 170 users, 9,517,761 ticks (cached replay).

- v2 gives a median **3.04× the ISF of v1** — weaker corrections — and is weaker on **92%** of readings.
- The margin depends on glucose (the equations use different glucose terms):

| glucose band | n | median ISF_v2 / ISF_v1 |
|---|---|---|
| (40, 80] | 778,636 | 53.11 |
| (80, 100] | 1,671,653 | 6.24 |
| (100, 120] | 2,063,028 | 3.5 |
| (120, 150] | 2,265,314 | 2.53 |
| (150, 200] | 1,784,604 | 1.89 |
| (200, 360] | 906,116 | 1.48 |

The v2 log approaches zero as glucose nears its floor, so v2's ISF climbs steeply below ~100 mg/dL (near-zero correction when low) and settles to roughly 1.5× v1 when high. v2 is the gentler equation everywhere, markedly so at low glucose.

*Counterfactual replay on the open-source cohort; high cap 210, low floor 76.*