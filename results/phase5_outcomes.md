# Phase 5 — does the derived (sensitivity-anchored) ISF beat tuned ISF?

114 users with realised CGM outcomes + valid empirical ISF.
Median R = entered/empirical = **2.38** (users dose, on average, well weaker than their measured sensitivity).

## Directional test (Spearman across users)

If the derived/empirical ISF is *right*, weaker-than-sensitivity dosing (high R) should mean more time-high and less time-low.

| association | ρ | expected sign |
|---|---|---|
| logR_vs_TAR (expect +) | -0.01 | + |
| logR_vs_TBR (expect -) | +0.39 | − |
| logR_vs_TIR | -0.10 | ? |
| mismatch_vs_TIR (expect -) | -0.10 | − |
| mismatch_vs_TBR | +0.39 | ? |

## Outcomes by entered/empirical ratio tertile

| R band | n | median R | TIR | TAR | TBR |
|---|---|---|---|---|---|
| low (≈ doses to sensitivity) | 38 | 1.6 | 85% | 14% | 1.6% |
| mid | 38 | 2.4 | 82% | 16% | 3.0% |
| high (doses weak) | 38 | 4.5 | 82% | 12% | 4.4% |

## Reading

- Associations are weak / not in the predicted direction → the large entered-vs-empirical gap does **not** translate into the expected glycaemic signal. Most likely the empirical level is biased (over-estimates sensitivity), so a sensitivity-anchored ISF would dose too strongly. **Favours the profile-anchored (Tier-1) design**, which preserves the user's working level and only adds the √TDD shape.

*Caveat:* closed-loop basal/SMB/autosens partially compensate a mis-set ISF, attenuating these associations; this is decision-level/observational, single cohort.