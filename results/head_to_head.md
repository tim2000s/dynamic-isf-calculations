# Same-window head-to-head: static vs dynamic ISF as a drop predictor

75 oref0 users, 55,001 overnight future-rise-selected windows. Each ISF form is tested on the same windows by rescaling the loop's IOB prediction. Error = actual end BG − predicted end BG; lower MAE = better. The activity proxy comes from the loop prediction, and the historical rise filter selects on the later outcome.

## Overall (median |error|, mg/dL)

| form | MAE | bias | n |
|---|---|---|---|
| static | 20.0 | 6.3 | 55,001 |
| v1 | 23.5 | -2.5 | 51,673 |
| v2 | 56.9 | 48.2 | 51,673 |
| loop | 18.0 | 6.0 | 55,001 |

## Bias by glucose (mg/dL; >0 = over-predicts the drop)

| BG band | static | v1 | v2 | loop |
|---|---|---|---|---|
| 80-100 | 15.8 | 29.3 | 139.3 | 16.0 |
| 100-120 | 6.0 | 12.9 | 77.6 | 5.0 |
| 120-145 | 5.2 | 5.8 | 64.1 | 4.0 |
| 145-175 | 6.0 | -3.4 | 48.8 | 7.0 |
| 175-230 | 6.4 | -26.9 | 11.8 | 8.2 |

## Per-user best predictor (lowest MAE), 67 users

| form | users where best |
|---|---|
| static | 19 |
| v1 | 12 |
| v2 | 3 |
| loop | 33 |

![Head-to-head](charts/inv008/fig_head_to_head.png)

*Per-window dataset (features + realised ISF) saved to `results/head_to_head_windows.parquet` for the pattern/ML step. v6 excluded (iob/isf accounting did not reconcile). OpenAPS Data Commons only.*