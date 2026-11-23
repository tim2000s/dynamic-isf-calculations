# Best-fit ISF equation search — LOUO-CV comparison

2026-11-23 12:41 · empirical target n=116 · entered target n=141

Fitted power law (empirical): ISF = 101.4·TDD^-0.431
Fitted power law (entered):  ISF = 349.6·TDD^-0.5
Blend weight on entered ISF: 0.1

## Target: delta-IOB observational outcome proxy

| candidate | fitted_cv | median_abs_err | p75_abs_err | median_log_err | frac_within_30pct |
|---|---|---|---|---|---|
| TDD-quartile bands | True | 6.2 | 9.3 | 0.282 | 0.44 |
| Blend entered×power law | True | 5.9 | 9.8 | 0.283 | 0.466 |
| Fitted C/TDD | True | 6.8 | 10.9 | 0.284 | 0.448 |
| K/sqrt(TDD) | True | 6.0 | 10.5 | 0.295 | 0.414 |
| Power law A·TDD^b | True | 6.4 | 11.1 | 0.297 | 0.405 |
| Power law + basal_frac | True | 6.3 | 11.1 | 0.299 | 0.388 |
| Multivariate (TDD,CR,basal,target) | True | 6.2 | 10.0 | 0.307 | 0.405 |
| Power law + ln(CR) | True | 6.3 | 10.2 | 0.315 | 0.431 |
| 1700-rule | False | 18.6 | 33.8 | 0.627 | 0.138 |
| v1 (TDD^-1) | False | 29.0 | 50.2 | 0.856 | 0.052 |
| Entered profile ISF | False | 31.0 | 49.9 | 0.896 | 0.043 |
| v2 (TDD^-2) | False | 207.0 | 512.0 | 2.427 | 0.0 |

## Target: entered ISF (user-tuned profile)

| candidate | fitted_cv | median_abs_err | p75_abs_err | median_log_err | frac_within_30pct |
|---|---|---|---|---|---|
| Power law A·TDD^b | True | 12.1 | 30.0 | 0.243 | 0.539 |
| K/sqrt(TDD) | True | 12.6 | 27.7 | 0.256 | 0.511 |
| Power law + basal_frac | True | 13.6 | 26.1 | 0.256 | 0.518 |
| TDD-quartile bands | True | 13.8 | 29.0 | 0.26 | 0.511 |
| Multivariate (TDD,CR,basal,target) | True | 13.7 | 28.7 | 0.289 | 0.468 |
| Power law + ln(CR) | True | 13.9 | 29.5 | 0.297 | 0.454 |
| v1 (TDD^-1) | False | 16.7 | 39.5 | 0.322 | 0.397 |
| Fitted C/TDD | True | 18.1 | 38.5 | 0.326 | 0.369 |
| 1700-rule | False | 16.1 | 35.1 | 0.412 | 0.369 |
| v2 (TDD^-2) | False | 158.8 | 385.8 | 1.508 | 0.043 |