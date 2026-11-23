# Same-window head-to-head: static vs dynamic ISF as a drop predictor

89 v5/v7 users, 62,824 overnight future-rise-selected windows. Each ISF form is tested on the same windows by rescaling the loop's IOB prediction. Error = actual end BG − predicted end BG; lower MAE = better. The activity proxy comes from the loop prediction, and the historical rise filter selects on the later outcome.

## Overall (median |error|, mg/dL)

| form | MAE | bias | n |
|---|---|---|---|
| static | 20.5 | 8.1 | 62,824 |
| v1 | 23.3 | -1.3 | 57,151 |
| v2 | 60.1 | 52.5 | 57,151 |
| loop | 18.6 | 7.4 | 62,824 |

## Bias by glucose (mg/dL; >0 = over-predicts the drop)

| BG band | static | v1 | v2 | loop |
|---|---|---|---|---|
| 80-100 | 14.9 | 27.9 | 134.3 | 15.0 |
| 100-120 | 6.7 | 12.3 | 77.6 | 5.4 |
| 120-145 | 7.0 | 7.3 | 70.1 | 5.8 |
| 145-175 | 7.9 | -1.7 | 54.7 | 8.2 |
| 175-230 | 10.6 | -24.6 | 15.8 | 10.6 |

## Per-user best predictor (lowest MAE), 80 users

| form | users where best |
|---|---|
| static | 22 |
| v1 | 19 |
| v2 | 4 |
| loop | 35 |

![Head-to-head](charts/inv008/fig_head_to_head.png)

*Per-window dataset (features + realised ISF) saved to `results/head_to_head_windows.parquet` for the pattern/ML step. v6 excluded (iob/isf accounting did not reconcile). Trio ISF under 20 converted from mmol/L, reading by reading.*