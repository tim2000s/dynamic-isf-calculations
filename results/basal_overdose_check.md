# Is the overnight 'reversion' actually scheduled-basal insulin?

82 people, 42,028 overnight carb-screened windows with **|iob_iob| < 0.2** (net IOB ≈ 0, so insulin delivery ≈ scheduled basal).

## Result

- BG drift: pooled median **-5.0 mg/dL/h**; per-user median **-5.5** [IQR -8.5, -3.0], **93% of users falling**.
- Drift by starting glucose (mg/dL/h):

| BG band | drift | n |
|---|---|---|
| 100-115 | -5.0 | 15,776 |
| 115-130 | -7.0 | 7,895 |
| 130-150 | -11.0 | 2,619 |
| 150-175 | -18.5 | 1,198 |
| 175-220 | -19.5 | 1,372 |

![Basal over-dose](charts/inv008/fig_basal_overdose.png)

**Conclusion.** At neutral net IOB overnight, glucose falls (median ~-5.5 mg/dL/h, 91% of users) — scheduled basal is delivering net insulin over EGP. This is the insulin that iob_iob does not capture; the 'reversion' is over-basalisation, not spontaneous glucose decline. Sensitivity from drop/iob_iob is therefore biased; total insulin (scheduled basal + temp + bolus) is the right denominator.

*IOB note: for v5/v7, iob_iob = bolusiob + basaliob exactly (internally consistent); basaliob is the temp-basal deviation and is typically negative, so iob_iob is net of scheduled basal. v6 (AAPS) iob_iob does not reconcile with its component columns and needs separate handling. A gross IOB calculation error is ruled out for v5/v7; the issue is that scheduled basal is not in iob_iob.*