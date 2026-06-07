# Phase 10 — multi-patient prediction backtest for the default glucose exponent k

13 patients, 16,432 windows, end-of-insulin-action horizon. Global k shared, per-patient level α fit individually; prediction-error scored.

- **default k (equal patient weight): 0.5**; window-weighted: 1.6
- **leave-one-patient-out k: median 0.5** [range 0.5–0.6]
- mean per-patient MAE: power-law@best-k **17.33** vs log 17.35 vs loop 21.7

## How much does k matter? (mean per-patient MAE by k)

| k | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 4.0 |
|---|---|---|---|---|---|---|
| MAE | 17.46 | 17.72 | 18.02 | 18.36 | 18.68 | 19.28 |

![Multi-patient MAE vs k](charts/inv008/fig_multipatient_k.png)

## Per-patient best k

| patient | TDD | n | best k |
|---|---|---|---|
| henny425 | 6 | 704 | 0.5 |
| kelseyhuss | 17 | 658 | 2.3 |
| svns | 18 | 727 | 1.7 |
| boost(N=1) | 22 | 6663 | 2.5 |
| fuxchr | 27 | 1017 | 4.0 |
| aadiabetes | 27 | 1031 | 1.0 |
| mikens | 29 | 1177 | 0.5 |
| diajesse | 42 | 297 | 0.5 |
| andycgm | 43 | 272 | 2.2 |
| ns_rot6 | 46 | 2819 | 3.0 |
| noahr | 54 | 55 | 0.5 |
| eli | 77 | 40 | 0.5 |
| nightscout1 | 84 | 972 | 0.5 |

**Reading:** the population objective is minimised at k≈0.5 (LOPO median 0.5); but MAE varies only ~1.8 mg/dL across k 1–4, so the exponent is a weak lever — any moderate k in ~1.5–3 is near-optimal. Power-law beats both log and the loop on the mean. Use the LOPO k as the default; treat the curve shape as more important than the precise exponent.