# Same-window head-to-head: static vs dynamic ISF as a drop predictor

89 v5/v7 users, 62,751 overnight future-rise-selected windows. Each ISF form is tested on the same windows by rescaling the loop's IOB prediction. Error = actual end BG − predicted end BG; lower MAE = better. The activity proxy comes from the loop prediction, and the historical rise filter selects on the later outcome.

## Overall (median |error|, mg/dL)

| form | MAE | bias | n |
|---|---|---|---|
| static | 20.3 | 6.6 | 62,751 |
| v1 | 24.6 | -10.4 | 57,078 |
| v2 | 49.7 | 25.2 | 57,078 |
| loop | 18.6 | 7.4 | 62,751 |

## Bias by glucose (mg/dL; >0 = over-predicts the drop)

| BG band | static | v1 | v2 | loop |
|---|---|---|---|---|
| 80-100 | 14.7 | 14.9 | 64.2 | 15.0 |
| 100-120 | 6.5 | 3.7 | 43.4 | 5.4 |
| 120-145 | 6.0 | -1.8 | 43.2 | 5.8 |
| 145-175 | 5.0 | -12.2 | 29.9 | 8.0 |
| 175-230 | 7.8 | -33.2 | 1.4 | 10.6 |

## Per-user best predictor (lowest MAE), 80 users

| form | users where best |
|---|---|
| static | 24 |
| v1 | 17 |
| v2 | 4 |
| loop | 35 |

![Head-to-head](charts/inv008/fig_head_to_head.png)

*Per-window dataset (features + realised ISF) saved to `results/head_to_head_windows.parquet` for the pattern/ML step. v6 excluded (iob/isf accounting did not reconcile). Units cleaned (mmol-scale ISF ×18.018).*
