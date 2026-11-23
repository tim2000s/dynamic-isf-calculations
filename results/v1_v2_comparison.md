# v1 vs v2 comparison

v1 `ISF = 1800/(TDD·ln(BG/75+1))`; v2 `ISF = 115000/(TDD²·ln(BG_floored/75))`. 147 users, 7,481,069 ticks (cached replay).

- v2 gives a median **3.14× the ISF of v1** — weaker corrections — and is weaker on **98%** of readings.
- The margin depends on glucose (the equations use different glucose terms):

| glucose band | n | median ISF_v2 / ISF_v1 |
|---|---|---|
| (40, 80] | 573,532 | 56.25 |
| (80, 100] | 1,268,448 | 6.54 |
| (100, 120] | 1,606,508 | 3.71 |
| (120, 150] | 1,782,558 | 2.66 |
| (150, 200] | 1,440,451 | 2.0 |
| (200, 360] | 766,403 | 1.63 |

The v2 log approaches zero as glucose nears its floor, so v2's ISF climbs steeply below ~100 mg/dL (near-zero correction when low) and settles to roughly 1.5× v1 when high. v2 is the gentler equation everywhere, markedly so at low glucose.

*Counterfactual replay on the open-source cohort; high cap 210, low floor 76.*