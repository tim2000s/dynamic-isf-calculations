# Phase 3 — baseline empirical-Bayes shrinkage estimator

135 users. Population prior K_pop = **127.4** (ISF·√TDD), between-user spread τ = 0.514 in log (prior-only baseline error **±43%**).

Error = median |log(estimate) − log(full-history baseline)|, as ± percent. Own-data error uses MEASURED window-to-baseline spread (not the optimistic regression SE). Shrinkage weight w = τ²/(τ²+σ²_W).

| trailing window | n users | shrink weight w | prior-only | own-only | **shrinkage** |
|---|---|---|---|---|---|
| 7d | 135 | 0.60 | ±43% | ±37% | **±26%** |
| 14d | 134 | 0.68 | ±43% | ±30% | **±23%** |
| 30d | 127 | 0.79 | ±43% | ±23% | **±19%** |
| 60d | 116 | 0.84 | ±43% | ±19% | **±16%** |
| 90d | 97 | 0.88 | ±43% | ±16% | **±13%** |

- Own-data overtakes the prior at **~7 days**.
- Best achievable baseline error here: **±13%** at 90-day windows (shrinkage).
- Prior-only (zero own-data, cold start) baseline error: **±43%** — the floor a brand-new user starts at.