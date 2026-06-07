# Phase 3 — baseline empirical-Bayes shrinkage estimator

132 users. Population prior K_pop = **129.2** (ISF·√TDD), between-user spread τ = 0.522 in log (prior-only baseline error **±42%**).

Error = median |log(estimate) − log(full-history baseline)|, as ± percent. Own-data error uses MEASURED window-to-baseline spread (not the optimistic regression SE). Shrinkage weight w = τ²/(τ²+σ²_W).

| trailing window | n users | shrink weight w | prior-only | own-only | **shrinkage** |
|---|---|---|---|---|---|
| 7d | 132 | 0.61 | ±42% | ±37% | **±26%** |
| 14d | 131 | 0.69 | ±42% | ±30% | **±22%** |
| 30d | 124 | 0.79 | ±42% | ±23% | **±19%** |
| 60d | 114 | 0.84 | ±42% | ±19% | **±16%** |
| 90d | 95 | 0.88 | ±42% | ±16% | **±13%** |

- Own-data overtakes the prior at **~7 days**.
- Best achievable baseline error here: **±13%** at 90-day windows (shrinkage).
- Prior-only (zero own-data, cold start) baseline error: **±42%** — the floor a brand-new user starts at.