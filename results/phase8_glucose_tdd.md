# Phase 8 — glucose term + TDD exponent, jointly on the oref cohort

129 users, overnight clean windows. Reconciles the prior power-law glucose term with the cohort √TDD level.

## (b) Glucose term — power-law (target/BG)^k vs log scaler

- power-law fits per-window ISF better than log for **41%** of users
- median within-user R²: power-law 0.048 vs log 0.049
- fitted glucose exponent k: median **-0.6** [IQR -0.74–-0.36], positive (ISF falls as BG rises) for 7% of users

## (a) TDD exponent (with the power-law glucose term)

- between-user ISF-at-target ∝ 1/TDD^**0.288** → closer to **0.5 (√TDD)**

## Reading

- **(b) is not clearly supported by the observational oref data** (power-law wins only 41%; k often non-positive) — most likely the unannounced-carb confound biases the within-user glucose slope. The prior clean N=1 backtest remains the stronger evidence for the power-law; oref data neither confirms nor refutes it.
- **(a):** with a power-law glucose term, the between-user level scales as 1/TDD^0.288 — consistent with √TDD; the glucose-term choice does not overturn the TDD-exponent finding (the two axes are largely separable, glucose term ≈ 1 at target).

*Caveat: overnight clean windows reduce but do not remove the carb confound; an observed-ISF-vs-BG regression cannot recover g(BG) at all (regression-to-mean inverts the sign); the glucose curve is a control/safety construct validated by prediction-error, as in the prior backtest. Per-window local-ISF is also hypo-biased in level (Phase 5/6), affecting the constant more than the exponents.*