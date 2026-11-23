# Overnight glucose-drop per starting IOB at a 4-hour horizon

99 people, 78,668 future-rise-selected overnight windows (11pm–3am start). Outcome ratio = (BG(T) − BG(T+4h)) / IOB(T), mg/dL per U. The denominator omits insulin delivered during the window, including algorithmic temporary basal and boluses. This is not an identified physiological ISF.

## Population

- Per-person median sensitivity: median **31.6** mg/dL/U [IQR 21.4–43.6].
- Pooled per-window: median 29.6 [IQR 18.0–46.3].
- Observed sensitivity-vs-glucose exponent: **k ≈ -1.12** (v-next quartic is +1.3; positive = falls with glucose).

## Sensitivity vs starting glucose (normalised to target)

| BG(T) | observed | v-next g(BG) | v1 | v2 |
|---|---|---|---|---|
| 106 | 1.00 | 0.92 | 0.96 | 0.81 |
| 118 | 1.75 | 0.79 | 0.89 | 0.61 |
| 135 | 2.08 | 0.66 | 0.82 | 0.47 |
| 160 | 2.16 | 0.51 | 0.74 | 0.37 |
| 198 | 2.29 | 0.38 | 0.65 | 0.29 |

![Overnight sensitivity](charts/inv008/fig_overnight_sensitivity.png)

## Per-person (first 30 by median sensitivity)


*Caveat: the denominator omits insulin delivered during the four-hour window, including algorithmic temporary basal and boluses. The future-rise filter uses the endpoint trajectory for selection. Treat the result as a descriptive ratio.*