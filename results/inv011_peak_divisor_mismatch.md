# INV-011: Peak and divisor mismatch

This is a software-model sensitivity analysis. It does not estimate clinical outcomes.

## Inverse-ISF effect of divisor 55

| Glucose | V1, 55 versus 75 | V2, 55 versus 75 |
|---:|---:|---:|
| 80 mg/dL | 1.24 times | 5.81 times |
| 99 mg/dL | 1.22 times | 2.12 times |
| 120 mg/dL | 1.21 times | 1.66 times |
| 150 mg/dL | 1.20 times | 1.45 times |
| 180 mg/dL | 1.19 times | 1.35 times |
| 210 mg/dL | 1.18 times | 1.30 times |

The multiplier is the insulin assigned per mg/dL of modelled effect. It does not mean that a positive correction is delivered at or below target.

## Effect still attributed to one unit at glucose 99 mg/dL

| Time | Reference 45-minute IOB | Assumed 70-minute IOB | V1 combined ratio | V2 combined ratio |
|---:|---:|---:|---:|---:|
| 15 min | 0.954 | 0.978 | 0.84 | 0.48 |
| 30 min | 0.850 | 0.922 | 0.89 | 0.51 |
| 45 min | 0.724 | 0.846 | 0.95 | 0.55 |
| 60 min | 0.598 | 0.760 | 1.04 | 0.60 |
| 70 min | 0.520 | 0.701 | 1.10 | 0.64 |
| 90 min | 0.381 | 0.582 | 1.25 | 0.72 |
| 120 min | 0.226 | 0.421 | 1.52 | 0.88 |
| 180 min | 0.067 | 0.186 | 2.27 | 1.31 |
| 240 min | 0.015 | 0.062 | 3.38 | 1.95 |
| 300 min | 0.002 | 0.011 | 4.92 | 2.85 |

V1 crosses from underestimating to overestimating effect still to come at 53 minutes. V2 crosses at 139 minutes in this reference case.

The mismatch increases calculated correction at delivery, understates early action, and retains too much modelled effect later. A full controller can moderate or amplify those errors through IOB limits and low-glucose prediction. This analysis does not predict the net delivered insulin or clinical outcome of a live controller.
