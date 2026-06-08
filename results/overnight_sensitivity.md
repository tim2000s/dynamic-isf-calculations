# Overnight insulin sensitivity at a 4-hour horizon

113 people, 83,653 carb-screened overnight windows (11pm–3am start). Sensitivity = (BG(T) − BG(T+4h)) / IOB(T), mg/dL per U.

## Population

- Per-person median sensitivity: median **31.4** mg/dL/U [IQR 20.9–41.7, range -15.5–107.9].
- Pooled per-window: median 30.1 [IQR 18.3–46.4].
- Observed sensitivity-vs-glucose exponent: **k ≈ -1.41** (v-next quartic is +1.3; positive = falls with glucose).

## Sensitivity vs starting glucose (normalised to target)

| BG(T) | observed | v-next g(BG) | v1 | v2 |
|---|---|---|---|---|
| 80 | -0.40 | 1.26 | 1.16 | 4.30 |
| 98 | 1.00 | 1.02 | 1.01 | 1.06 |
| 112 | 1.69 | 0.85 | 0.92 | 0.68 |
| 130 | 2.55 | 0.69 | 0.84 | 0.50 |
| 155 | 2.68 | 0.54 | 0.75 | 0.38 |
| 195 | 2.83 | 0.39 | 0.66 | 0.29 |

![Overnight sensitivity](charts/inv008/fig_overnight_sensitivity.png)

## Per-person (first 30 by median sensitivity)

| user | n win | median BG(T) | median IOB | median drop | median sens | IQR |
|---|---|---|---|---|---|---|
| U120 | 50 | 77 | 3.73 | -62 | -16 | -26–10 |
| U125 | 37 | 170 | 4.84 | 33 | 6 | -13–13 |
| U180 | 58 | 149 | 3.87 | 14 | 6 | -14–11 |
| U171 | 117 | 111 | 3.69 | 40 | 7 | 4–43 |
| U170 | 172 | 156 | 4.60 | 54 | 8 | 5–20 |
| U033 | 514 | 115 | 2.61 | 21 | 8 | 5–12 |
| U076 | 48 | 110 | 3.11 | 28 | 8 | 2–11 |
| U145 | 250 | 124 | 2.52 | 38 | 10 | 3–21 |
| U089 | 58 | 138 | 0.79 | 31 | 12 | 5–118 |
| U116 | 36 | 136 | 3.75 | 52 | 12 | 3–22 |
| U165 | 255 | 168 | 5.35 | 61 | 13 | 5–20 |
| U166 | 3068 | 116 | 1.61 | 26 | 13 | 9–22 |
| U174 | 337 | 189 | 5.34 | 57 | 14 | 1–23 |
| U068 | 334 | 142 | 2.74 | 41 | 15 | 11–21 |
| U042 | 42 | 140 | 3.24 | 74 | 15 | 13–24 |
| U022 | 421 | 139 | 3.47 | 43 | 15 | 10–24 |
| U118 | 1973 | 141 | 2.35 | 40 | 15 | 10–22 |
| U106 | 46 | 177 | 5.19 | 69 | 16 | 7–22 |
| U126 | 139 | 124 | 1.64 | 26 | 16 | -1–29 |
| U147 | 1727 | 141 | 2.70 | 44 | 17 | 13–22 |
| U111 | 64 | 169 | 4.94 | 94 | 18 | 9–27 |
| U127 | 528 | 148 | 2.84 | 51 | 18 | 6–45 |
| U142 | 702 | 124 | 1.62 | 30 | 19 | 7–36 |
| U009 | 56 | 400 | 16.26 | 221 | 19 | 11–20 |
| U055 | 1072 | 151 | 2.42 | 47 | 19 | 15–23 |
| U065 | 758 | 161 | 2.26 | 48 | 19 | 13–26 |
| U182 | 207 | 144 | 2.75 | 52 | 19 | 12–25 |
| U153 | 725 | 134 | 1.91 | 42 | 20 | 16–26 |
| U088 | 136 | 186 | 3.07 | 79 | 21 | 14–39 |
| U083 | 6552 | 141 | 1.28 | 29 | 21 | 13–32 |

*Caveat: basal continues over the 4 h and roughly offsets endogenous glucose in a fasting state, so drop/IOB(T) is an approximation of sensitivity; residual dawn effect, basal mis-set, and counter-regulation still bias the tails.*