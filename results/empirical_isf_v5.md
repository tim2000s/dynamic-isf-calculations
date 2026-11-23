# Delta-IOB outcome proxy on the canonical cohort

Historical filename: `empirical_isf_v5`. Per-user OLS: `ΔBG = a + b·ΔIOB + c·BG_pre_trend`; the reported proxy is `-b`. Delta IOB is model-derived and is changed by insulin delivered during the 30-minute window. Closed-loop delivery also responds to glucose. The coefficient does not identify physiological ISF. Per-user intervals describe regression uncertainty conditional on that model.

## Coverage

- Canonical cohort users in the analysis: 119
- Users with valid empirical_isf in [5,500] AND R² ≥ 0.10: **98**

## Cohort-level results

- Median empirical ISF: **22** [Q1 16, Q3 30] mg/dL/U
- Median entered ISF: **50** mg/dL/U
- Median empirical/entered ratio: **0.44** [Q1 0.28, Q3 0.60]
- Median R²: 0.23
- Median windows per user: 5092

## Stratified by DynISF group

| Group | n | Median entered | Median empirical | Median ratio | Median R² | Median per-user 95% CI half-width |
|---|---|---|---|---|---|---|
| **no_dynisf** | 98 | 50 | 22 | 0.44 | 0.23 | ±1 mg/dL/U |
| **dynisf_sigmoid** | 0 | (n<5) | | | | |
| **dynisf_log** | 0 | (n<5) | | | | |

**Per-user uncertainty caveat**: median per-user 95% CI half-width shows the regression noise per individual estimate. Cohort-level medians are robust because errors cancel; individual user numbers have wide CIs and should not be reported as 'this user's empirical ISF is X' without the band.
