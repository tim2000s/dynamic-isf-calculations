# Canonical Walsh-constants analysis (n = 119)

Single cohort, single TDD definition, used by every downstream artefact.

- TDD = `basal + Σ treatments / span` for v6/v7 (treatments-derived); `mean_tdd` from live extraction for v5 (mathematically equivalent — same delivery integration, different upstream).
- Quality filter: ISF [10,300], CR [2,50], target [70,130], TDD [5,200], n_days ≥ 14.
- Bootstrap B = 2000.

## Cohort composition

| Source | n | TDD method | Median ISF | Median CR | Median TDD |
|---|---|---|---|---|---|
| v6_aaps_classic | 22 | treatments_plus_basal | 42 | 7.5 | 61.2 |
| v7_oref0 | 97 | treatments_plus_basal | 50 | 9.0 | 40.0 |

## Walsh constants (95 % bootstrap CI)

| Constant | n | Median | 95 % CI | Walsh | CI excludes Walsh? |
|---|---|---|---|---|---|
| ISF × TDD | 119 | 2121 | [1875, 2464] | 1700 | **yes** |
| CR × TDD | 119 | 376 | [314, 403] | 500 | **yes** |
| basal / TDD | 119 | 0.48 | [0.46, 0.53] | 0.5 | no |

## log(ISF) ~ log(TDD) slope (95 % bootstrap CI on slope)

| Group | n | Slope b | 95 % CI on slope | Intercept a |
|---|---|---|---|---|
| **no_dynisf** | 119 | -0.48 | [-0.65, -0.30] | 5.74 |
| dynisf_sigmoid | 0 | (n<5) | | |
| dynisf_log | 0 | (n<5) | | |
| ALL | 119 | **-0.48** | [-0.65, -0.30] | 5.74 |

Walsh's slope = −1. No group's slope interval contains −1.

## Outlier sensitivity (drop |z(ISF×TDD)| > 2)

- Dropped 3 people
- ISF × TDD trimmed: 2108 [1826, 2415]
- CR × TDD trimmed: 366 [313, 402]
- basal / TDD trimmed: 0.48 [0.46, 0.53]

## Duration sensitivity

| Threshold | n | ISF × TDD median [95 % CI] |
|---|---|---|
| ≥ 14 | 119 | 2121 [1875, 2464] |
| ≥ 30 | 114 | 2116 [1850, 2454] |
| ≥ 60 | 103 | 2057 [1813, 2382] |
| ≥ 90 | 94 | 2050 [1769, 2376] |
| ≥ 180 | 67 | 1980 [1680, 2369] |
