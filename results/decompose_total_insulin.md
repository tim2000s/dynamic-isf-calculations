# Overnight sensitivity on TOTAL delivered insulin (basal + temp + bolus)

89 v6/v7 users, 92,918 overnight carb-screened windows. Dose = insulin delivered over the window from the reconstructed grid (includes scheduled basal). Model: drop = a + b·INS + c·(BG−100) + d·INS·(BG−100).

## Result

- **ISF at 100 (pooled): -0.2 mg/dL per U** — still off.
- EGP rise with zero insulin: **-1.4 mg/dL/h** (glucose rises without insulin — as it should).
- Residual glucose term c = 0.814 (was the 'reversion'; smaller now that scheduled basal is counted).
- Glucose interaction d = 0.001 → -0.75%/mg/dL (rises with glucose).
- At the same glucose, high- vs low-insulin windows now differ by **1 mg/dL** (was ~7 with iob_iob) — insulin is starting to separate from the glucose level.

## Q2 — per-user ISF at 100

- Median **-1.3** mg/dL per U [IQR -5.5–-0.0], 24% positive (80 people).

![Total-insulin decomposition](charts/inv008/fig_decompose_total_insulin.png)

*v6/v7 only (they have the reconstructed delivery grid; v5/Trio excluded). The grid's total_u reflects actual delivery including temp reductions and suspensions. CGM aligned to the grid by the stage-2 recovered anchor.*