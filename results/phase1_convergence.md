# Phase 1 — empirical-ISF convergence & reproducibility

141 users with a fit (of 141 canonical).

## (1) Growing trailing window — within-window precision

| window | n users fit | median rel. 95% CI half-width | ≤±10% | ≤±15% | ≤±20% | ≤±25% |
|---|---|---|---|---|---|---|
| 7d | 135 | ±30% | 5% | 13% | 24% | 36% |
| 14d | 139 | ±24% | 9% | 26% | 41% | 53% |
| 30d | 140 | ±17% | 21% | 44% | 64% | 74% |
| 60d | 140 | ±11% | 39% | 64% | 76% | 83% |
| 90d | 140 | ±10% | 48% | 72% | 81% | 88% |

## (2) Test-retest — reproducibility of a 14-day estimate

- users with ≥3 fourteen-day blocks: **132**
- median block-to-block CV (true reproducibility): **±34%** [Q1 27%, Q3 45%]
- median within-block regression SE (the model's own claim): ±11%
- → the regression SE understates true variability by ~3.3×
- fraction of users reproducible within ±10% (whole-history CV): 2%
- fraction of users reproducible within ±15% (whole-history CV): 4%
- fraction of users reproducible within ±20% (whole-history CV): 8%
- fraction of users reproducible within ±25% (whole-history CV): 20%

### Drift vs jitter (adjacent fortnight-to-fortnight)

- median adjacent-block change: **±33%** (vs whole-history ±34%)
- median lag-1 autocorrelation of block estimates: -0.00 (mostly uncorrelated noise — not trackable at this resolution)
- fraction with adjacent change ≤ ±10%: 3%
- fraction with adjacent change ≤ ±15%: 8%
- fraction with adjacent change ≤ ±20%: 17%
- fraction with adjacent change ≤ ±25%: 27%