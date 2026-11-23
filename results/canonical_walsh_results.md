# Canonical Walsh-constants analysis (n = 141)

Single cohort, single TDD definition, used by every downstream artefact.

- TDD = `basal + Σ treatments / span` for v6/v7 (treatments-derived); `mean_tdd` from live extraction for v5 (mathematically equivalent — same delivery integration, different upstream).
- Quality filter: ISF [10,300], CR [2,50], target [70,130], TDD [5,200], n_days ≥ 14.
- Bootstrap B = 2000.

## Cohort composition

| Source | n | TDD method | Median ISF | Median CR | Median TDD |
|---|---|---|---|---|---|
| v5_trio | 22 | live_extracted_mean_tdd | 61 | 10.0 | 40.1 |
| v6_aaps_classic | 22 | treatments_plus_basal | 42 | 7.5 | 61.2 |
| v7_oref0 | 97 | treatments_plus_basal | 50 | 9.0 | 40.0 |

## Walsh constants (95 % bootstrap CI)

| Constant | n | Median | 95 % CI | Walsh | CI excludes Walsh? |
|---|---|---|---|---|---|
| ISF × TDD | 141 | 2236 | [2030, 2482] | 1700 | **yes** |
| CR × TDD | 141 | 376 | [315, 402] | 500 | **yes** |
| basal / TDD | 141 | 0.47 | [0.45, 0.52] | 0.5 | no |

## log(ISF) ~ log(TDD) slope (95 % bootstrap CI on slope)

| Group | n | Slope b | 95 % CI on slope | Intercept a |
|---|---|---|---|---|
| **no_dynisf** | 123 | -0.46 | [-0.63, -0.30] | 5.67 |
| **dynisf_sigmoid** | 9 | -0.89 | [-1.44, -0.35] | 7.44 |
| **dynisf_log** | 9 | -0.84 | [-1.35, -0.25] | 7.36 |
| ALL | 141 | **-0.50** | [-0.67, -0.34] | 5.86 |

Walsh's slope = −1. The slope interval contains −1 for: dynisf_sigmoid, dynisf_log.

## Outlier sensitivity (drop |z(ISF×TDD)| > 2)

- Dropped 3 people
- ISF × TDD trimmed: 2155 [1986, 2469]
- CR × TDD trimmed: 366 [315, 401]
- basal / TDD trimmed: 0.47 [0.45, 0.53]

## Duration sensitivity

| Threshold | n | ISF × TDD median [95 % CI] |
|---|---|---|
| ≥ 14 | 141 | 2236 [2030, 2482] |
| ≥ 30 | 136 | 2198 [2032, 2524] |
| ≥ 60 | 124 | 2140 [1960, 2464] |
| ≥ 90 | 106 | 2108 [1792, 2459] |
| ≥ 180 | 77 | 2023 [1703, 2369] |
