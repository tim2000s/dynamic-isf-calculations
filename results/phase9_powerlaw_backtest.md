# Phase 9 — consolidated-equation backtest (N=1 overnight prediction data)

9031 overnight fasting cycles; tdd_7day range 7.7–26.7 U/day. Prediction-error design (scale the loop's 2h drop by isf_cand/isf_loop, compare to actual).

## Overall MAE (mg/dL, lower better)

| formula | MAE | bias |
|---|---|---|
| loop (its own ISF) | 26.26 | – |
| log scaler /TDD | 21.47 | +7.68 |
| **power-law (target/BG)^4.25, /TDD** | **20.15** | +3.49 |
| power-law, /√TDD (best k=4.0) | 20.3 | – |

- power-law beats log by **6.2%**, and the loop by 23.3%.
- best glucose exponent **k = 4.25** (prior work found ≈3.5).
- √TDD vs 1/TDD: MAE 20.3 vs 20.15 — within-patient TDD range is small, so this N=1 data barely distinguishes the TDD exponent (that is the cohort's job; cohort = √TDD).

## Per-BG-band MAE

| BG band | n | loop | log | power-law |
|---|---|---|---|---|
| (70, 90] | 2722 | 26.2 | 26.7 | 23.4 |
| (90, 105] | 2318 | 26.8 | 19.7 | 19.5 |
| (105, 120] | 2640 | 24.1 | 17.0 | 16.2 |
| (120, 150] | 1341 | 29.1 | 22.5 | 22.1 |
| (150, 200] | 10 | 109.8 | 74.3 | 57.3 |

*N=1 (author's own closed-loop). Validates the glucose-curve shape and this patient's constant; the population TDD exponent (√TDD) comes from the cohort.*