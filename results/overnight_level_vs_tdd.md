# Does the clean overnight sensitivity level follow √TDD?

99 people. Level = per-person median overnight sensitivity (4-hour horizon, BG ≥ target, carb-screened); TDD = total daily dose. Fit log–log across people.

## Result

- Fitted exponent: **-0.445** [95% CI -0.701, -0.194] (n-weighted -0.338).
- √TDD (−0.5) is inside the CI; 1/TDD (−1) is excluded.
- √TDD constant: **K = 197** (ISF ≈ 197/√TDD).

## Fit comparison (median |log error|, lower = better)

| form | constant / exponent | median log err |
|---|---|---|
| √TDD | K=197 | 0.303 |
| free power | A=160, p=-0.445 | 0.316 |
| 1/TDD | K=1334 | 0.361 |

The chart draws one mark per person and is kept in the private data folder (REVIEW.md item 22).

## Reading

Measured overnight sensitivity falls with TDD at an exponent of -0.445 — the same negative direction as the earlier cross-sectional fits, and 1/TDD (v1) is firmly rejected. √TDD sits inside the confidence interval, so the clean data is consistent with — or a touch shallower than — a square-root law. This is the level result re-derived on measured (not fitted) sensitivity, independent of profile settings and of the equations themselves.

*Caveat: the level is a per-person median over supra-target overnight windows, so it still carries some of the glucose mean-reversion confound (people who run higher overnight read a touch more sensitive); this can bias the exponent if overnight glucose correlates with TDD. The direction and the rejection of 1/TDD are robust; the exact exponent is approximate.*