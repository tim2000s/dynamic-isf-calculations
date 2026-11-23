# Phase 1 — empirical-ISF convergence & reproducibility

119 users with a fit (of 119 canonical).

## (1) Growing trailing window — within-window precision

| window | n users fit | median rel. 95% CI half-width | ≤±10% | ≤±15% | ≤±20% | ≤±25% |
|---|---|---|---|---|---|---|
| 7d | 113 | ±29% | 6% | 15% | 27% | 38% |
| 14d | 117 | ±21% | 11% | 29% | 45% | 56% |
| 30d | 118 | ±17% | 23% | 44% | 64% | 73% |
| 60d | 118 | ±11% | 38% | 64% | 72% | 81% |
| 90d | 118 | ±11% | 45% | 69% | 78% | 86% |

## (2) Test-retest — reproducibility of a 14-day estimate

- users with ≥3 fourteen-day blocks: **110**
- median block-to-block CV (true reproducibility): **±33%** [Q1 25%, Q3 45%]
- median within-block regression SE (the model's own claim): ±10%
- → the regression SE understates true variability by ~3.2×
- fraction of users reproducible within ±10% (whole-history CV): 3%
- fraction of users reproducible within ±15% (whole-history CV): 5%
- fraction of users reproducible within ±20% (whole-history CV): 10%
- fraction of users reproducible within ±25% (whole-history CV): 22%

### Drift vs jitter (adjacent fortnight-to-fortnight)

- median adjacent-block change: **±32%** (vs whole-history ±33%)
- median lag-1 autocorrelation of block estimates: +0.02 (mostly uncorrelated noise — not trackable at this resolution)
- fraction with adjacent change ≤ ±10%: 3%
- fraction with adjacent change ≤ ±15%: 8%
- fraction with adjacent change ≤ ±20%: 19%
- fraction with adjacent change ≤ ±25%: 30%