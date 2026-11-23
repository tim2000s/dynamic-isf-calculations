# Best-fit ISF equation search — LOUO-CV comparison

2026-11-23 15:34 · empirical target n=98 · entered target n=119

Fitted power law (empirical): ISF = 101.3·TDD^-0.43
Fitted power law (entered):  ISF = 311.7·TDD^-0.478
Blend weight on entered ISF: 0.2

## Target: delta-IOB observational outcome proxy

| candidate | fitted_cv | median_abs_err | p75_abs_err | median_log_err | frac_within_30pct |
|---|---|---|---|---|---|
| Blend entered×power law | True | 6.0 | 10.2 | 0.26 | 0.5 |
| TDD-quartile bands | True | 6.7 | 9.7 | 0.283 | 0.459 |
| Fitted C/TDD | True | 7.0 | 11.5 | 0.293 | 0.429 |
| K/sqrt(TDD) | True | 6.5 | 10.2 | 0.294 | 0.429 |
| Power law + basal_frac | True | 6.9 | 11.2 | 0.315 | 0.388 |
| Multivariate (TDD,CR,basal,target) | True | 6.5 | 10.6 | 0.318 | 0.367 |
| Power law A·TDD^b | True | 6.8 | 11.4 | 0.325 | 0.398 |
| Power law + ln(CR) | True | 6.9 | 10.8 | 0.334 | 0.408 |
| 1700-rule | False | 18.0 | 31.3 | 0.592 | 0.163 |
| v1 (TDD^-1) | False | 28.9 | 46.5 | 0.822 | 0.061 |
| Entered profile ISF | False | 26.6 | 48.7 | 0.826 | 0.051 |
| v2 (TDD^-2) | False | 207.0 | 409.5 | 2.407 | 0.0 |

## Target: entered ISF (user-tuned profile)

| candidate | fitted_cv | median_abs_err | p75_abs_err | median_log_err | frac_within_30pct |
|---|---|---|---|---|---|
| Power law A·TDD^b | True | 13.1 | 27.1 | 0.267 | 0.496 |
| Power law + basal_frac | True | 14.0 | 23.6 | 0.277 | 0.479 |
| K/sqrt(TDD) | True | 12.8 | 25.6 | 0.281 | 0.454 |
| Power law + ln(CR) | True | 15.2 | 28.5 | 0.315 | 0.412 |
| Multivariate (TDD,CR,basal,target) | True | 13.5 | 28.5 | 0.319 | 0.445 |
| v1 (TDD^-1) | False | 15.8 | 40.3 | 0.331 | 0.387 |
| TDD-quartile bands | True | 15.0 | 27.0 | 0.336 | 0.429 |
| Fitted C/TDD | True | 16.1 | 40.5 | 0.341 | 0.395 |
| 1700-rule | False | 14.9 | 33.8 | 0.403 | 0.387 |
| v2 (TDD^-2) | False | 158.8 | 359.3 | 1.508 | 0.05 |