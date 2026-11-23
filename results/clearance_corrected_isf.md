# Test #1 (bounding): high-BG resistance after removing insulin-independent clearance

75 users, 54,649 insulin-active windows. Non-insulin flux estimated from windows where the loop expects ~no insulin action (|loop predicted drop|<5).

## Data-derived non-insulin flux (mg/dL over 4h)

| BG band | flux | n |
|---|---|---|
| 100-120 | 10.0 | 1,948 |
| 120-145 | 22.0 | 793 |
| 145-175 | 44.0 | 187 |
| 175-205 | 67.5 | 98 |
| 205-260 | 80.0 | 152 |

## Realised ÷ profile ISF ratio: raw vs clearance-corrected

| BG band | raw | corrected | clean+corrected |
|---|---|---|---|
| 100-120 | 0.58 | 0.29 | 0.34 |
| 120-145 | 0.87 | 0.33 | 0.3 |
| 145-175 | 0.85 | 0.19 | 0.2 |
| 175-205 | 0.89 | 0.15 | 0.1 |
| 205-260 | 0.88 | 0.24 | 0.24 |

**Implied k**: raw -0.04, corrected **0.62**, clean+corrected 0.68 (k>0 = resistance).

## Power-law bound (clearance required vs data-derived)

| BG band | required clearance | data-derived | median drop / pred |
|---|---|---|---|
| 175-205 | 52.5 | 67.5 | 77.0 / 81.6 |
| 205-260 | 75.3 | 80.0 | 109.5 / 114.0 |

**Verdict: RESISTANCE appears after clearance removal (corrected ratio falls with BG).**

![clearance corrected](charts/inv008/fig_clearance_corrected_isf.png)

*PLAUSIBILITY BOUND, not a clean number: clearance↔resistance partially entangled; low-insulin high-BG windows sparse; EGP not perfectly stable. nonInsulin includes EGP (can be negative near target = glucose rising).*