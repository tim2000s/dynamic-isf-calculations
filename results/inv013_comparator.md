# INV-013: what each equation calculates for one person

Profile sensitivity 50 mg/dL per unit, total daily dose 40 U, recent dose
equal to the person's own average. AndroidAPS at an adjustment factor of 100%, Trio at
its defaults (0.8 logarithmic, 0.5 sigmoid, limits 0.7 and 1.2, insulin factor 55).
Values are mg/dL per unit.

| Glucose, mg/dL | AndroidAPS, divisor 55 | AndroidAPS, divisor 65 | AndroidAPS, divisor 75 | Trio logarithmic | Trio sigmoid | Candidate | Candidate, option above 200 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 70 | 54.8 | 61.6 | 68.3 | 68.5 | 55.7 | 63.1 | 63.1 |
| 80 | 50.1 | 56.1 | 62.0 | 62.6 | 53.7 | 57.7 | 57.7 |
| 100 | 43.4 | 48.3 | 53.1 | 54.3 | 50.0 | 50.0 | 50.0 |
| 120 | 38.9 | 43.0 | 47.1 | 48.6 | 47.1 | 50.0 | 50.0 |
| 150 | 34.2 | 37.6 | 41.0 | 42.8 | 44.3 | 50.0 | 50.0 |
| 200 | 29.3 | 32.0 | 34.6 | 41.7 | 42.4 | 50.0 | 50.0 |
| 250 | 26.3 | 28.5 | 30.7 | 41.7 | 41.8 | 50.0 | 44.8 |
| 300 | 24.1 | 26.1 | 28.0 | 41.7 | 41.7 | 50.0 | 41.1 |

## Level the stock equation implies at 100 mg/dL

| Insulin setting | Divisor | Sensitivity x daily dose |
|---|---:|---:|
| rapid, 75 min peak | 55 | 1737 |
| ultra-rapid, 55 min peak | 65 | 1932 |
| 45 min peak | 75 | 2124 |

## Correction size against the stock equation

Stock AndroidAPS value at divisor 55 divided into the candidate's: above 1 the
candidate asks for less insulin for the same distance from target.

| Glucose, mg/dL | Candidate | Candidate, option above 200 |
|---|---:|---:|
| 70 | 1.15 | 1.15 |
| 80 | 1.15 | 1.15 |
| 100 | 1.15 | 1.15 |
| 120 | 1.29 | 1.29 |
| 150 | 1.46 | 1.46 |
| 200 | 1.70 | 1.70 |
| 250 | 1.90 | 1.70 |
| 300 | 2.07 | 1.70 |

Stock ratio of the value at 300 to the value at 200 mg/dL: 0.82 (divisor 55), 0.81 (divisor 65), 0.81 (divisor 75).

## Movement with recent dose

Multiplier on the sensitivity when the blended dose sits at the tenth and ninetieth
percentile of its ratio to the person's own average.

| Blended over average | Stock | Candidate, exponent 1 | Candidate, exponent 0.65 |
|---:|---:|---:|---:|
| 0.72 | 1.39 | 1.39 | 1.24 |
| 1.14 | 0.88 | 0.88 | 0.92 |

## Units to bring the worked person to 100 mg/dL

Distance from 100 mg/dL divided by the sensitivity at the starting glucose. This is
the arithmetic of the sensitivity alone and leaves out insulin on board, predictions
and delivery limits.

| Starting glucose, mg/dL | AndroidAPS, divisor 55 | Trio logarithmic | Candidate |
|---:|---:|---:|---:|
| 120 | 0.51 | 0.41 | 0.40 |
| 150 | 1.46 | 1.17 | 1.00 |
| 200 | 3.41 | 2.40 | 2.00 |
| 250 | 5.71 | 3.60 | 3.00 |
| 300 | 8.29 | 4.80 | 4.00 |

## Shape alone

Each value divided by the same equation's value at 100 mg/dL, divisor 55. The last two
columns are the candidate at the other two divisors.

| Glucose, mg/dL | Stock | Candidate | Candidate, divisor 65 | Candidate, divisor 75 |
|---:|---:|---:|---:|---:|
| 70 | 1.26 | 1.26 | 1.27 | 1.29 |
| 80 | 1.15 | 1.15 | 1.16 | 1.17 |
| 100 | 1.00 | 1.00 | 1.00 | 1.00 |
| 120 | 0.90 | 1.00 | 1.00 | 1.00 |
| 150 | 0.79 | 1.00 | 1.00 | 1.00 |
| 200 | 0.68 | 1.00 | 1.00 | 1.00 |
| 250 | 0.60 | 1.00 | 1.00 | 1.00 |
| 300 | 0.56 | 1.00 | 1.00 | 1.00 |

## Figures quoted in the proposal

- Level at 100 mg/dL, divisor 75 over divisor 55: 1.22.
- Trio logarithmic over AndroidAPS (divisor 55) at 200 mg/dL: 1.42.
- Candidate at 80 mg/dL with recent dose 1.10 of average: 52.4; from 100 mg/dL upward: 45.5.
- Spread of the easing at 70 mg/dL across divisors 55 to 75: 0.023.
