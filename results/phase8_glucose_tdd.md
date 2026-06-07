# Phase 8 — glucose term + TDD exponent, jointly on the oref cohort

129 users, overnight clean windows. Reconciles the prior power-law glucose term with the cohort √TDD level.

## (b) Glucose term — power-law (target/BG)^k vs log scaler

- power-law fits per-window ISF better than log for **41%** of users
- median within-user R²: power-law 0.048 vs log 0.049
- fitted glucose exponent k: median **-0.6** [IQR -0.74–-0.36], positive (ISF falls as BG rises) for 7% of users

## (a) TDD exponent (with the power-law glucose term)

- between-user ISF-at-target ∝ 1/TDD^**0.288** → closer to **0.5 (√TDD)**

## Reading

- **(b) is not testable on the observational oref data.** The fitted glucose exponent k is
  *negative* for ~93% of users — local-ISF measured as RISING with BG (0.70→1.23 normalised).
  That is the **opposite** of the established physiology: **ISF *falls* with BG — more insulin
  per mg/dL at high glucose** — which the log, power-law and Diabeloop curves all encode. The
  rising measured slope is therefore an **artefact**: regression-to-the-mean and
  counter-regulation (high BG is already trending down and that fall is mis-credited to insulin;
  lows are defended), **not** real sensitivity, and **not** the carb confound. The formula's
  falling-ISF direction is physiologically correct; an observed-ISF-vs-BG regression is
  wrong-signed by artefact and cannot recover g(BG). The prior clean N=1 prediction-error
  backtest remains the stronger (and right-signed) evidence for the power-law.
- **(a):** with a power-law glucose term, the between-user level scales as 1/TDD^0.288 —
  consistent with √TDD; the glucose-term choice does not overturn the TDD-exponent finding (the
  two axes are largely separable, glucose term ≈ 1 at target).

*Caveat: per-window local-ISF cannot measure g(BG) — its BG-slope is dominated by glucose
dynamics (mean-reversion/counter-regulation), giving the wrong sign. The glucose curve must be
established by prediction-error/outcome validation (the prior N=1 design), where ISF correctly
falls with BG. Local-ISF is also hypo-biased in level (Phase 5/6).*