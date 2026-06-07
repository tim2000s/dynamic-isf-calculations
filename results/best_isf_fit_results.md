# Best-fit ISF equation search — LOUO-CV comparison

2026-06-07 19:13 · empirical target n=114 · entered target n=138

Fitted power law (empirical): ISF = 86.0·TDD^-0.376
Fitted power law (entered):  ISF = 285.1·TDD^-0.433
Blend weight on entered ISF: 0.0

## Target: empirical ISF (observed sensitivity)

| candidate | fitted_cv | median_abs_err | p75_abs_err | median_log_err | frac_within_30pct |
|---|---|---|---|---|---|
| Power law + basal_frac | True | 6.0 | 10.9 | 0.294 | 0.439 |
| K/sqrt(TDD) | True | 6.2 | 10.7 | 0.297 | 0.447 |
| Fitted C/TDD | True | 6.7 | 12.0 | 0.298 | 0.386 |
| TDD-quartile bands | True | 6.6 | 10.5 | 0.303 | 0.421 |
| Power law A·TDD^b | True | 6.1 | 11.1 | 0.304 | 0.447 |
| Blend entered×power law | True | 6.1 | 11.1 | 0.304 | 0.447 |
| Multivariate (TDD,CR,basal,target) | True | 6.2 | 10.2 | 0.306 | 0.395 |
| Power law + ln(CR) | True | 6.5 | 11.2 | 0.33 | 0.412 |
| 1700-rule | False | 16.2 | 31.2 | 0.605 | 0.149 |
| v1 (TDD^-1) | False | 26.0 | 48.0 | 0.813 | 0.07 |
| Entered profile ISF | False | 32.0 | 50.3 | 0.896 | 0.044 |
| v2 (TDD^-2) | False | 171.1 | 464.1 | 2.323 | 0.018 |

## Target: entered ISF (user-tuned profile)

| candidate | fitted_cv | median_abs_err | p75_abs_err | median_log_err | frac_within_30pct |
|---|---|---|---|---|---|
| K/sqrt(TDD) | True | 12.8 | 30.8 | 0.256 | 0.514 |
| Power law A·TDD^b | True | 13.5 | 29.4 | 0.268 | 0.493 |
| Power law + basal_frac | True | 14.4 | 24.8 | 0.285 | 0.442 |
| TDD-quartile bands | True | 13.8 | 28.9 | 0.292 | 0.478 |
| Multivariate (TDD,CR,basal,target) | True | 14.7 | 28.5 | 0.306 | 0.457 |
| Power law + ln(CR) | True | 15.4 | 28.4 | 0.318 | 0.428 |
| v1 (TDD^-1) | False | 17.0 | 37.8 | 0.324 | 0.377 |
| Fitted C/TDD | True | 19.3 | 37.0 | 0.336 | 0.435 |
| 1700-rule | False | 17.9 | 36.7 | 0.454 | 0.341 |
| v2 (TDD^-2) | False | 124.0 | 375.7 | 1.304 | 0.065 |