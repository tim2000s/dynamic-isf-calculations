# v1 vs v2 comparison

v1 `ISF = 1800/(TDD·ln(BG/75+1))`; v2 `ISF = 115000/(TDD²·ln(BG_floored/75))`. 170 users, 9,461,361 ticks (cached replay).

- v2 gives a median **3.32× the ISF of v1** — weaker corrections — and is weaker on **98%** of readings.
- The margin depends on glucose (the equations use different glucose terms):

| glucose band | n | median ISF_v2 / ISF_v1 |
|---|---|---|
| (40, 80] | 775,966 | 59.74 |
| (80, 100] | 1,664,035 | 6.84 |
| (100, 120] | 2,050,609 | 3.82 |
| (120, 150] | 2,250,045 | 2.75 |
| (150, 200] | 1,771,259 | 2.06 |
| (200, 360] | 901,054 | 1.63 |

The v2 log approaches zero as glucose nears its floor, so v2's ISF climbs steeply below ~100 mg/dL (near-zero correction when low) and settles to roughly 1.5× v1 when high. v2 is the gentler equation everywhere, markedly so at low glucose.

*Counterfactual replay on the open-source cohort; high cap 210, low floor 76.*