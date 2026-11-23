# Best-fit individualised glucose-ISF (gradient/separable-NLLS over a shared shape)

67 users, 54,635 windows. Model: actual_drop ≈ a_u + s_u·(pred_drop·(100/BG)^k). Per-user (a_u,s_u) fit by contiguous within-user 5-fold CV; shared k searched. **k=0 = flat magnitude model; does k>0 help?**

## Headline

| model | within-user out-of-fold MAE |
|---|---|
| static (no fit) | 24.44 |
| magnitude (k=0, per-user scale) | 18.85 |
| **best-fit glucose k*=0.0 (per-user scale)** | **18.85** |
| diabeloop shape, per-user scaled | 24.32 |
| cold-start best (no adaptation, k=0.25) | 22.98 |

**Glucose shape adds 0.0 mg/dL beyond magnitude** (k*=0.0). Per-user optimal k: median 0.25, 48% want k=0, 42% want k≥1.

Per-user scale s predictable from TDD: ρ=-0.07 (p=0.562); from profile ISF: ρ=0.01 (p=0.949).

![gradient fit](charts/inv008/fig_gradient_isf_fit.png)

*k=0 is flat (=magnitude model). best_fit_k>0 with lower MAE ⇒ an individualised glucose shape helps BEYOND magnitude. cold-start = population shape with no per-user adaptation. per-user optimal k shows whether people want different steepness.*