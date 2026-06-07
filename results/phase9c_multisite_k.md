# Phase 9c — glucose exponent k across 12 sites (full insulin action)

Per-site prediction-error fit of ISF = α·(target/BG)^k at the end-of-insulin-action horizon (~150–225 min). TDD spans 6–84 U/day.

- **median best k = 0.75** [IQR 0.5–2.25, range 0.5–4.0], n-weighted 1.9
- power-law ≥ log for **12/12** sites
- between-site TDD exponent (ISF∝1/TDD^p): **p = 0.739** → closer to **0.5 (√TDD)**

| site | model | n | TDD | best k | MAE pl | MAE log | MAE loop |
|---|---|---|---|---|---|---|---|
| henny425 | sigmoid | 704 | 6 | **0.50** | 16.63 | 16.73 | 20.52 |
| kelseyhuss | log | 658 | 17 | **2.25** | 17.27 | 17.51 | 23.9 |
| svns | sigmoid | 727 | 18 | **1.75** | 21.98 | 22.27 | 29.21 |
| fuxchr | sigmoid | 1017 | 27 | **4.00** | 12.77 | 13.08 | 13.81 |
| aadiabetes | sigmoid | 1031 | 27 | **1.00** | 13.62 | 13.62 | 15.31 |
| mikens | sigmoid | 1177 | 29 | **0.50** | 16.89 | 17.0 | 21.15 |
| diajesse | sigmoid | 297 | 42 | **0.50** | 20.57 | 20.82 | 22.41 |
| andycgm | log | 272 | 43 | **2.25** | 17.41 | 18.79 | 18.64 |
| ns_rot6 | log | 2819 | 46 | **3.00** | 13.03 | 13.23 | 13.77 |
| noahr | log | 55 | 54 | **0.50** | 24.98 | 25.39 | 35.62 |
| eli | log | 40 | 77 | **0.50** | 13.14 | 13.14 | 18.52 |
| nightscout1 | log | 972 | 84 | **0.50** | 14.62 | 14.67 | 28.1 |

- Single-patient (boost cache) gave k≈2.25 at ~3.17h; the multisite median (0.75) differs.

*N small per site (40–2819 windows); k is noisy per site but the central tendency firms up the exponent. Prediction-error design; mixed sigmoid/log loop formulas (the scaling is formula-agnostic).*