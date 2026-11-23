# Static-ISF prediction error vs glucose (full curve)

67 users, 54,635 overnight windows. `err_static > 0` means the person's static profile ISF *over-predicted the drop* (the real fall was smaller). Population value is the median across users; CI is bootstrapped over users.

## Absolute over-prediction (mg/dL)

| BG band | median | 90% CI | users over-pred | users |
|---|---|---|---|---|
| 80-100 | 16.79 | [13.0, 29.8] | 90% | 29 |
| 100-120 | 5.54 | [2.83, 12.49] | 73% | 51 |
| 120-145 | 6.45 | [2.33, 10.06] | 65% | 63 |
| 145-175 | 9.34 | [4.45, 13.21] | 65% | 65 |
| 175-205 | 12.77 | [0.84, 16.61] | 63% | 63 |
| 205-260 | 14.71 | [2.66, 23.41] | 64% | 58 |

## Fractional over-prediction (effective-sensitivity gap, %)

| BG band | median | 90% CI | users |
|---|---|---|---|
| 80-100 | 88.0% | [76.0%, 100.0%] | 29 |
| 100-120 | 24.0% | [14.0%, 46.0%] | 51 |
| 120-145 | 16.0% | [9.0%, 22.0%] | 63 |
| 145-175 | 14.0% | [5.0%, 19.0%] | 65 |
| 175-205 | 14.0% | [1.0%, 18.0%] | 63 |
| 205-260 | 12.0% | [2.0%, 17.0%] | 58 |

## Absolute curve, raw vs confound-adjusted (mg/dL)

| BG band | raw | adjusted |
|---|---|---|
| 80-100 | 16.79 | 33.9 |
| 100-120 | 5.54 | 21.37 |
| 120-145 | 6.45 | 14.75 |
| 145-175 | 9.34 | 9.23 |
| 175-205 | 12.77 | 0.26 |
| 205-260 | 14.71 | -2.08 |

![Error curve](charts/inv008/fig_err_curve.png)

*Absolute and fractional disagree by construction: corrections are bigger at high glucose, so a flat proportional mismatch shows up as a larger mg/dL error there. The fractional panel is the physiological view; the confound-adjusted panel shows how much of the curve survives removing dawn/IOB/hour structure.*