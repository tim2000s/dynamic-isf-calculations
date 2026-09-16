# Best-fit individualised glucose-ISF (gradient/separable-NLLS over a shared shape)

80 users, 62,344 windows. Model: actual_drop ≈ a_u + s_u·(pred_drop·(100/BG)^k). Per-user (a_u,s_u) fit by contiguous within-user 5-fold CV; shared k searched. **k=0 = flat magnitude model; does k>0 help?**

## Headline

| model | within-user out-of-fold MAE |
|---|---|
| static (no fit) | 24.16 |
| magnitude (k=0, per-user scale) | 18.66 |
| **best-fit glucose k*=0.0 (per-user scale)** | **18.66** |
| diabeloop shape, per-user scaled | 24.03 |
| cold-start best (no adaptation, k=0.75) | 22.44 |

**Glucose shape adds 0.0 mg/dL beyond magnitude** (k*=0.0). Per-user optimal k: median 0.25, 49% want k=0, 41% want k≥1.

Per-user scale s predictable from TDD: ρ=0.03 (p=0.819); from profile ISF: ρ=-0.17 (p=0.135).

![gradient fit](charts/inv008/fig_gradient_isf_fit.png)

*k=0 is flat (=magnitude model). best_fit_k>0 with lower MAE ⇒ an individualised glucose shape helps BEYOND magnitude. cold-start = population shape with no per-user adaptation. per-user optimal k shows whether people want different steepness.*