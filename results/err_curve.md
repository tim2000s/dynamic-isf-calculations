# Static-ISF prediction error vs glucose (full curve)

> Audit status: this is a prediction-error analysis within the loop model. It is not a
> physiological view of insulin sensitivity.

80 users, 62,344 overnight windows. `err_static > 0` means the person's static profile ISF *over-predicted the drop* (the real fall was smaller). Population value is the median across users; CI is bootstrapped over users.

## Absolute over-prediction (mg/dL)

| BG band | median | 90% CI | users over-pred | users |
|---|---|---|---|---|
| 80-100 | 16.79 | [13.0, 29.8] | 91% | 33 |
| 100-120 | 7.37 | [4.34, 12.49] | 78% | 64 |
| 120-145 | 8.59 | [5.99, 11.66] | 70% | 76 |
| 145-175 | 9.37 | [4.8, 13.21] | 65% | 78 |
| 175-205 | 15.03 | [7.67, 18.29] | 68% | 72 |
| 205-260 | 14.93 | [2.36, 24.1] | 63% | 67 |

## Fractional over-prediction (effective-sensitivity gap, %)

| BG band | median | 90% CI | users |
|---|---|---|---|
| 80-100 | 83.0% | [75.0%, 94.0%] | 33 |
| 100-120 | 30.0% | [20.0%, 44.0%] | 64 |
| 120-145 | 21.0% | [15.0%, 27.0%] | 76 |
| 145-175 | 14.0% | [6.0%, 19.0%] | 78 |
| 175-205 | 15.0% | [7.0%, 18.0%] | 72 |
| 205-260 | 13.0% | [1.0%, 18.0%] | 67 |

## Absolute curve, raw vs confound-adjusted (mg/dL)

| BG band | raw | adjusted |
|---|---|---|
| 80-100 | 16.79 | 35.96 |
| 100-120 | 7.37 | 23.26 |
| 120-145 | 8.59 | 17.57 |
| 145-175 | 9.37 | 9.32 |
| 175-205 | 15.03 | 7.0 |
| 205-260 | 14.93 | 0.78 |

![Error curve](charts/inv008/fig_err_curve.png)

*Absolute and fractional disagree by construction: corrections are bigger at high glucose, so a flat proportional mismatch shows up as a larger mg/dL error there. The fractional panel is the physiological view; the confound-adjusted panel shows how much of the curve survives removing dawn/IOB/hour structure.*
