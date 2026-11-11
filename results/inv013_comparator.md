# INV-013: the curve against the stock dynamic ISF

## Sensitivity x daily dose by glucose

The stock equation is 1800 / ln(glucose / divisor + 1). The curve is
1915 x quartic(glucose) / quartic(180), held below 100 and above 380 mg/dL.
The published column is the poster's fit in mg/dL per unit: the quadratic up to
100 mg/dL and the quartic above. The curve uses the quartic's value at 100 (80.6).

| Glucose, mg/dL | Stock, divisor 55 | Stock, divisor 65 | Stock, divisor 75 | Published | Curve | Curve over stock 55 | Curve over stock 75 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 70 | 2193 | 2463 | 2730 | 66.1 | 4378 | 2.00 | 1.60 |
| 80 | 2005 | 2243 | 2480 | 68.6 | 4378 | 2.18 | 1.77 |
| 100 | 1737 | 1932 | 2124 | 79.0 | 4378 | 2.52 | 2.06 |
| 120 | 1555 | 1721 | 1884 | 63.5 | 3449 | 2.22 | 1.83 |
| 150 | 1368 | 1505 | 1638 | 45.9 | 2492 | 1.82 | 1.52 |
| 180 | 1239 | 1357 | 1471 | 35.3 | 1915 | 1.55 | 1.30 |
| 200 | 1173 | 1281 | 1385 | 30.8 | 1674 | 1.43 | 1.21 |
| 250 | 1051 | 1141 | 1228 | 24.8 | 1346 | 1.28 | 1.10 |
| 300 | 965 | 1043 | 1118 | 21.1 | 1148 | 1.19 | 1.03 |
| 380 | 870 | 936 | 998 | 15.3 | 830 | 0.95 | 0.83 |

Ratio of the value at 100 to the value at 200 mg/dL: curve 2.61, stock 1.48. At 300 over 200: curve 0.69, stock 0.82.
Stock level at 100 mg/dL, divisor 75 over divisor 55: 1.22.

Glucose at which the curve equals the stock equation: 362 mg/dL (divisor 55), 339 mg/dL (divisor 65), 313 mg/dL (divisor 75).

## Reference dose

The published curve carries no daily dose. Linking it to one needs the dose at
which the published values apply, which is a level divided by the quartic.

| Level taken | Reference dose, U | Curve x daily dose at 200 mg/dL |
|---|---:|---:|
| Adult entered product 1915 at 180 mg/dL (this draft) | 54.3 | 1674 |
| Pooled entered product 2067 at 180 mg/dL | 58.6 | 1807 |
| 1800 rule at 105 mg/dL (March 2026 work) | 23.7 | 731 |
| Own peak-to-trough level under Loop at 200 mg/dL | 70.8 | 2182 |
| Own peak-to-trough level under Trio at 200 mg/dL | 73.7 | 2273 |
| Own peak-to-trough level under oref0 at 200 mg/dL | 77.8 | 2398 |
| Own peak-to-trough level under no algorithm (REPLACE-BG) at 200 mg/dL | 106.5 | 3283 |

The adult product carries a 95% interval of 1828 to 1986, which at 200 mg/dL is 1598 to 1736.
At 180 mg/dL the TDD link gives 95.8 mg/dL per unit on 20 U a day and 23.9 on 80 U a day. The published value there is 35.3 for both.

Dose multiple against the reference glucose: 0.56 at 120 mg/dL and 1.95 at 340 mg/dL. The abstract gives about a half and about double.

## One person

Profile sensitivity 50 mg/dL per unit, total daily dose 40 U, recent dose
equal to the person's own average. AndroidAPS at an adjustment factor of 100%, Trio at
its defaults (0.8 logarithmic, 0.5 sigmoid, limits 0.7 and 1.2, insulin factor 55).
Values are mg/dL per unit.

| Glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Trio sigmoid | Curve, TDD link | Curve, profile link | Profile, eased below 100 |
|---:|---:|---:|---:|---:|---:|---:|
| 70 | 54.8 | 68.5 | 55.7 | 109.4 | 114.3 | 63.1 |
| 80 | 50.1 | 62.6 | 53.7 | 109.4 | 114.3 | 57.7 |
| 100 | 43.4 | 54.3 | 50.0 | 109.4 | 114.3 | 50.0 |
| 120 | 38.9 | 48.6 | 47.1 | 86.2 | 90.1 | 50.0 |
| 150 | 34.2 | 42.8 | 44.3 | 62.3 | 65.1 | 50.0 |
| 180 | 31.0 | 41.7 | 42.9 | 47.9 | 50.0 | 50.0 |
| 200 | 29.3 | 41.7 | 42.4 | 41.9 | 43.7 | 50.0 |
| 250 | 26.3 | 41.7 | 41.8 | 33.7 | 35.1 | 50.0 |
| 300 | 24.1 | 41.7 | 41.7 | 28.7 | 30.0 | 50.0 |
| 380 | 21.8 | 41.7 | 41.7 | 20.7 | 21.7 | 50.0 |

Inside Trio's default limits (41.7 to 71.4) the curve runs only from 137 to 201 mg/dL for this person.

## Units to bring that person to 100 mg/dL

Distance from 100 mg/dL over the sensitivity. This is the arithmetic of the sensitivity
alone and leaves out insulin on board, predictions and delivery limits. An event is
taken to be filed under the glucose halfway down its fall, which the abstract does not
state. On that reading, using the curve at the starting glucose asks for more insulin
than its construction implies.

| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Profile | Curve, read at start | Curve, read at midpoint | Curve, along the path |
|---:|---:|---:|---:|---:|---:|---:|
| 150 | 1.46 | 1.17 | 1.00 | 0.80 | 0.61 | 0.62 |
| 200 | 3.41 | 2.40 | 2.00 | 2.39 | 1.61 | 1.62 |
| 250 | 5.71 | 3.60 | 3.00 | 4.46 | 3.01 | 2.97 |
| 300 | 8.29 | 4.80 | 4.00 | 6.97 | 4.78 | 4.58 |

Reading at the start asks for 31 to 49% more than reading at the midpoint.
For this person the TDD link equals the profile value at 174 mg/dL; the profile link does so at 180.

## Movement with recent dose

The daily dose in the link is this draft's to define. Dividing by
the existing blend moves the value by the same multiplier as the stock equation:
1.39 at a blend of 0.72 of the person's average, 0.88 at a blend of 1.14 of the person's average.
