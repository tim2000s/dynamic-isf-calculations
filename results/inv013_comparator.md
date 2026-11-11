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

## The published curve against the measured levels

Measured sensitivity x daily dose in each glucose bin (INV-012, median over events with
positive insulin used, 95% interval over people) divided by the published value at that
glucose. The result is the daily dose, in U, at which the published curve returns the
measurement. A row that stays level means the published shape fits that cohort.

| Cohort | People | 100 | 150 | 200 | 250 | 300 |
|---|---:|---:|---:|---:|---:|---:|
| Loop | 842 | 88 (84 to 91) | 58 (56 to 60) | 71 (67 to 73) | 82 (77 to 89) | 108 (95 to 119) |
| Trio | 22 | 90 (63 to 111) | 52 (40 to 64) | 74 (51 to 104) | - | - |
| oref0 | 86 | 109 (90 to 136) | 70 (59 to 86) | 78 (58 to 188) | - | - |
| no algorithm (REPLACE-BG) | 196 | 103 (89 to 115) | 100 (94 to 108) | 106 (94 to 125) | 120 (103 to 137) | 174 (132 to 233) |

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

The profile link over AndroidAPS at 100% for this person: 2.63 at 100, 1.90 at 150, 1.49 at 200, 1.24 at 300 mg/dL. The two are equal at 378 mg/dL with a divisor of 55 and at 330 with 75. The profile link equals Trio's logarithmic mode at 209 mg/dL and its sigmoid at 206.

Trio's logarithmic mode is the AndroidAPS equation divided by its adjustment factor, held inside
its limits. At the default factor of 0.8 it is 25% above AndroidAPS at 100% until the upper limit
binds, which for this person is from 157 mg/dL. With the factor at 1.0 and the limits wide it
returns the AndroidAPS value: 29.3 against 29.3 mg/dL per unit at 200 mg/dL.

Inside Trio's default limits (41.7 to 71.4) the curve on the profile link runs only from 141 to 208 mg/dL for this person.

## Units to bring that person to 100 mg/dL

Distance from 100 mg/dL over the sensitivity, with the curve on the profile link. This is
the arithmetic of the sensitivity alone and leaves out insulin on board, predictions and
delivery limits. An event is taken to be filed under the glucose halfway down its fall,
which the abstract does not state. On that reading, using the curve at the starting glucose asks for more insulin
than its construction implies.

| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Profile | Curve, read at start | Curve, read at midpoint | Curve, along the path |
|---:|---:|---:|---:|---:|---:|---:|
| 150 | 1.46 | 1.17 | 1.00 | 0.77 | 0.59 | 0.59 |
| 200 | 3.41 | 2.40 | 2.00 | 2.29 | 1.54 | 1.55 |
| 250 | 5.71 | 3.60 | 3.00 | 4.27 | 2.89 | 2.85 |
| 300 | 8.29 | 4.80 | 4.00 | 6.67 | 4.58 | 4.39 |

Reading at the start asks for 31 to 49% more than reading at the midpoint.
For this person the TDD link equals the profile value at 174 mg/dL; the profile link does so at 180.

## Movement with recent dose

The daily dose in the link is this draft's to define. Dividing by
the existing blend moves the value by the same multiplier as the stock equation:
1.39 at a blend of 0.72 of the person's average, 0.88 at a blend of 1.14 of the person's average.
