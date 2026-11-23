# Phase 2 — sensitivity variance decomposition by timescale

115 users with ≥150 clean windows (ΔIOB ≥ 0.3 U).
Per-window log local-ISF, centred per user (baseline removed); variance attributed to each timescale (one-way η²); daily-series autocorrelation separates persistent (trackable) from white.

## Variance budget (median across users, fraction of within-user log-sensitivity variance)

- **circadian (hour-of-day)**: 2%  [Q1 1%, Q3 5%]
- **weekly (day-of-week)**: 1%  [Q1 1%, Q3 3%]
- **slow (~monthly)**: 2%  [Q1 1%, Q3 5%]
- **residual (day-to-day stochastic + measurement noise)**: ~92%

## Is the day-to-day residual trackable or white?

- median daily lag-1 autocorrelation: **+0.09** (weak/none → largely white at daily scale)
- median daily lag-7 autocorrelation: +0.02

## Circadian shape (population mean dev by bucket, log units)

| bucket | mean log-dev | ≈ relative ISF |
|---|---|---|
| night | -0.011 | ×0.99 |
| dawn | -0.035 | ×0.97 |
| midday | -0.001 | ×1.00 |
| afternoon | +0.019 | ×1.02 |
| evening | +0.030 | ×1.03 |
| late | +0.041 | ×1.04 |

## Real vs artefact (stricter ΔIOB ≥ 0.6 U)

- circadian η²: 2% (base) → 4% (strict)
- residual: 92% (base) → 88% (strict)
- → residual shrinks under stricter filter ⇒ part of it is estimator artefact; circadian survives ⇒ real physiology.