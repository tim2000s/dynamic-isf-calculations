# Phase 9b — consolidated-equation backtest at end-of-insulin-action (~3.17h)

6663 overnight cycles, horizon median 3.17h (end of the loop's IOB prediction). Compare to Phase 9 at +2h.

## Overall (MAE mg/dL; bias = actual − predicted, + = formula too aggressive)

| formula | MAE | bias |
|---|---|---|
| loop (its own ISF) | 21.1 | -0.76 |
| log scaler /TDD | 19.03 | +7.76 |
| **power-law (target/BG)^2.25 /TDD** | **18.37** | +4.1 |
| power-law /√TDD (k=2.5) | 18.33 | – |

- power-law beats log by **3.5%**, loop by 12.9%; best k = **2.25**.
- bias at ~3.17h: loop -0.76, log +7.76, power-law +4.1 (vs the larger +bias seen at 2h).

## Per-BG-band MAE

| BG band | n | loop | log | power-law |
|---|---|---|---|---|
| (50, 90] | 2330 | 26.4 | 26.1 | 24.7 |
| (90, 105] | 1793 | 17.5 | 14.8 | 14.8 |
| (105, 120] | 1721 | 17.5 | 13.9 | 13.2 |
| (120, 150] | 819 | 21.4 | 19.1 | 19.2 |

*N=1, end-of-insulin-action horizon captures more of the total ISF than +2h. Validates the glucose-curve shape; TDD exponent is the cohort's (√TDD).*