# Observed sensitivity vs glucose — per-window, against the calculated curves

98 people, 1,114,690 fasting windows. Within each glucose band, observed ISF = −b from ΔBG = a + b·ΔIOB + c·trend (the per-window measurement, pooled only within the band so it is not collapsed across glucose).

## Observed vs calculated (normalised to 1.0 at target)

| BG | n windows | observed ISF | observed (norm) | v-next g(BG) | v1 | v2 |
|---|---|---|---|---|---|---|
| 78 | 98059 | 10.0 | 0.70 | 1.30 | 1.19 | 8.47 |
| 90 | 133811 | 12.4 | 0.86 | 1.12 | 1.07 | 1.52 |
| 100 | 163360 | 14.3 | 1.00 | 0.99 | 0.99 | 0.97 |
| 110 | 162591 | 15.4 | 1.07 | 0.88 | 0.93 | 0.72 |
| 122 | 196124 | 15.7 | 1.10 | 0.76 | 0.87 | 0.57 |
| 140 | 163383 | 16.5 | 1.15 | 0.62 | 0.80 | 0.44 |
| 162 | 103053 | 15.6 | 1.09 | 0.50 | 0.73 | 0.36 |
| 192 | 62240 | 15.8 | 1.10 | 0.40 | 0.66 | 0.29 |

![Observed vs calculated](charts/inv008/fig_window_shape.png)

## Reading

- Observed sensitivity **falls with glucose**, fitted power-law exponent **k ≈ -0.4** over this range; the v-next quartic's exponent over the same range is 1.3.
- This is the like-for-like the averaged anchor could not give: the measured value compared to the calculated value *at each glucose*, rather than one pooled slope.

*Caveat: 30-min fasting windows, ΔIOB ∈ (0,2] U, trend-adjusted; observational, so counter-regulation at low BG and unrecorded carbs/EGP still confound the band estimates, most at the extremes where windows are fewer.*