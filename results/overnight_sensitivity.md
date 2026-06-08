# Overnight insulin sensitivity at a 4-hour horizon

112 people, 81,814 carb-screened overnight windows (11pm–3am start). Sensitivity = (BG(T) − BG(T+4h)) / IOB(T), mg/dL per U.

## Population

- Per-person median sensitivity: median **31.6** mg/dL/U [IQR 22.0–42.6, range 5.9–107.9].
- Pooled per-window: median 30.5 [IQR 18.8–46.9].
- Observed sensitivity-vs-glucose exponent: **k ≈ -1.13** (v-next quartic is +1.3; positive = falls with glucose).

## Sensitivity vs starting glucose (normalised to target)

| BG(T) | observed | v-next g(BG) | v1 | v2 |
|---|---|---|---|---|
| 106 | 1.00 | 0.92 | 0.96 | 0.81 |
| 118 | 1.69 | 0.79 | 0.89 | 0.61 |
| 135 | 2.05 | 0.66 | 0.82 | 0.47 |
| 160 | 2.20 | 0.51 | 0.74 | 0.37 |
| 198 | 2.26 | 0.38 | 0.65 | 0.29 |

![Overnight sensitivity](charts/inv008/fig_overnight_sensitivity.png)

## Per-person (first 30 by median sensitivity)

| user | n win | median BG(T) | median IOB | median drop | median sens | IQR |
|---|---|---|---|---|---|---|
| U125 | 37 | 170 | 4.84 | 33 | 6 | -13–13 |
| U180 | 56 | 153 | 3.84 | 15 | 6 | -13–11 |
| U170 | 172 | 156 | 4.60 | 54 | 8 | 5–20 |
| U033 | 488 | 116 | 2.74 | 22 | 8 | 6–12 |
| U076 | 36 | 117 | 3.66 | 30 | 9 | 7–12 |
| U089 | 58 | 138 | 0.79 | 31 | 12 | 5–118 |
| U116 | 33 | 141 | 3.59 | 52 | 13 | 3–23 |
| U166 | 2560 | 122 | 1.88 | 30 | 13 | 9–21 |
| U174 | 322 | 190 | 5.45 | 62 | 14 | 0–23 |
| U165 | 248 | 170 | 5.50 | 62 | 14 | 5–20 |
| U068 | 334 | 142 | 2.74 | 41 | 15 | 11–21 |
| U145 | 200 | 130 | 3.19 | 44 | 15 | 6–22 |
| U022 | 412 | 142 | 3.54 | 44 | 15 | 10–22 |
| U042 | 34 | 150 | 3.96 | 76 | 15 | 14–30 |
| U106 | 44 | 178 | 5.18 | 74 | 16 | 7–24 |
| U118 | 1809 | 143 | 2.45 | 44 | 16 | 11–23 |
| U147 | 1695 | 143 | 2.78 | 45 | 17 | 13–22 |
| U009 | 56 | 400 | 16.26 | 221 | 19 | 11–20 |
| U171 | 82 | 145 | 4.24 | 82 | 19 | 6–51 |
| U055 | 1055 | 151 | 2.45 | 47 | 19 | 15–23 |
| U065 | 758 | 161 | 2.26 | 48 | 19 | 13–26 |
| U153 | 725 | 134 | 1.91 | 42 | 20 | 16–26 |
| U182 | 198 | 146 | 2.71 | 57 | 20 | 13–26 |
| U111 | 59 | 171 | 5.25 | 96 | 20 | 10–28 |
| U088 | 128 | 187 | 3.24 | 81 | 21 | 15–39 |
| U142 | 632 | 126 | 1.63 | 34 | 21 | 8–37 |
| U083 | 6454 | 141 | 1.29 | 29 | 22 | 14–33 |
| U090 | 163 | 208 | 4.10 | 92 | 22 | 15–40 |
| U119 | 174 | 240 | 4.27 | 92 | 22 | 13–38 |
| U126 | 122 | 136 | 1.69 | 33 | 22 | 7–30 |

*Caveat: basal continues over the 4 h and roughly offsets endogenous glucose in a fasting state, so drop/IOB(T) is an approximation of sensitivity; residual dawn effect, basal mis-set, and counter-regulation still bias the tails.*