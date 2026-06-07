# Phase 1 — empirical-ISF convergence & reproducibility

138 users with a fit (of 138 canonical).

## (1) Growing trailing window — within-window precision

| window | n users fit | median rel. 95% CI half-width | ≤±10% | ≤±15% | ≤±20% | ≤±25% |
|---|---|---|---|---|---|---|
| 7d | 133 | ±30% | 5% | 14% | 24% | 35% |
| 14d | 136 | ±24% | 10% | 26% | 41% | 53% |
| 30d | 137 | ±17% | 21% | 44% | 64% | 74% |
| 60d | 137 | ±11% | 39% | 64% | 76% | 83% |
| 90d | 137 | ±10% | 48% | 72% | 81% | 88% |

## (2) Test-retest — reproducibility of a 14-day estimate

- users with ≥3 fourteen-day blocks: **129**
- median block-to-block CV (true reproducibility): **±34%** [Q1 26%, Q3 44%]
- median within-block regression SE (the model's own claim): ±10%
- → the regression SE understates true variability by ~3.3×
- fraction of users reproducible within ±10% (whole-history CV): 2%
- fraction of users reproducible within ±15% (whole-history CV): 4%
- fraction of users reproducible within ±20% (whole-history CV): 9%
- fraction of users reproducible within ±25% (whole-history CV): 20%

### Drift vs jitter (adjacent fortnight-to-fortnight)

- median adjacent-block change: **±34%** (vs whole-history ±34%)
- median lag-1 autocorrelation of block estimates: -0.00 (mostly uncorrelated noise — not trackable at this resolution)
- fraction with adjacent change ≤ ±10%: 3%
- fraction with adjacent change ≤ ±15%: 7%
- fraction with adjacent change ≤ ±20%: 16%
- fraction with adjacent change ≤ ±25%: 26%