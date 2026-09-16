# Test #1 (bounding): high-BG resistance after removing insulin-independent clearance

> Audit status: withdrawn as a causal result. Low-action closed-loop windows do not identify
> insulin-independent clearance, so subtraction does not identify insulin-only sensitivity.

88 users, 62,084 insulin-active windows. Non-insulin flux estimated from windows where the loop expects ~no insulin action (|loop predicted drop|<5).

## Data-derived non-insulin flux (mg/dL over 4h)

| BG band | flux | n |
|---|---|---|
| 100-120 | 10.0 | 2,017 |
| 120-145 | 23.0 | 848 |
| 145-175 | 44.0 | 220 |
| 175-205 | 67.0 | 106 |
| 205-260 | 81.0 | 158 |

## Realised ÷ profile ISF ratio: raw vs clearance-corrected

| BG band | raw | corrected | clean+corrected |
|---|---|---|---|
| 100-120 | 0.59 | 0.3 | 0.35 |
| 120-145 | 0.84 | 0.31 | 0.27 |
| 145-175 | 0.86 | 0.21 | 0.24 |
| 175-205 | 0.87 | 0.17 | 0.14 |
| 205-260 | 0.89 | 0.25 | 0.24 |

**Implied k**: raw -0.1, corrected **0.44**, clean+corrected 0.45 (k>0 = resistance).

## Power-law bound (clearance required vs data-derived)

| BG band | required clearance | data-derived | median drop / pred |
|---|---|---|---|
| 175-205 | 50.9 | 67.0 | 76.0 / 83.8 |
| 205-260 | 74.5 | 81.0 | 109.0 / 115.0 |

**Verdict: RESISTANCE appears after clearance removal (corrected ratio falls with BG).**

![clearance corrected](charts/inv008/fig_clearance_corrected_isf.png)

*PLAUSIBILITY BOUND, not a clean number: clearance↔resistance partially entangled; low-insulin high-BG windows sparse; EGP not perfectly stable. nonInsulin includes EGP (can be negative near target = glucose rising).*
