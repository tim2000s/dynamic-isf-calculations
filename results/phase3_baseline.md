# Phase 3 — baseline empirical-Bayes shrinkage estimator

113 users. Population prior K_pop = **130.3** (ISF·√TDD), between-user spread τ = 0.545 in log (prior-only baseline error **±47%**).

Error = median |log(estimate) − log(full-history baseline)|, as ± percent. Own-data error uses MEASURED window-to-baseline spread (not the optimistic regression SE). Shrinkage weight w = τ²/(τ²+σ²_W).

| trailing window | n users | shrink weight w | prior-only | own-only | **shrinkage** |
|---|---|---|---|---|---|
| 7d | 113 | 0.65 | ±47% | ±36% | **±27%** |
| 14d | 112 | 0.72 | ±47% | ±29% | **±24%** |
| 30d | 105 | 0.81 | ±47% | ±22% | **±19%** |
| 60d | 95 | 0.84 | ±47% | ±18% | **±17%** |
| 90d | 85 | 0.90 | ±47% | ±17% | **±13%** |

- Own-data overtakes the prior at **~7 days**.
- Best achievable baseline error here: **±13%** at 90-day windows (shrinkage).
- Prior-only (zero own-data, cold start) baseline error: **±47%** — the floor a brand-new user starts at.