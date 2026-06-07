# dynamic-isf-calculations

Analysis of equations that set a dynamic insulin sensitivity factor (ISF) from total daily
dose (TDD) and current glucose, evaluated against real-world data from 171 people using
open-source automated insulin delivery (AID) systems.

Dynamic ISF makes correction sensitivity a function of TDD. The original equation
(**v1**, Chris Wilson) makes the sensitivity anchor inversely proportional to TDD; a later
revision of the maths (**v2**) makes it inversely proportional to TDD squared. This work
generates each person's dynamic ISF under both equations from their own glucose and
insulin history, tests both against sensitivity calculated directly from the data, and
proposes a next version.

## The equations

Both share the same TDD-blending step (five windows → a weighted TDD) and the same
glucose scaler. They differ only in the TDD term:

| | sensitivity anchor at normal target | implied law |
|---|---|---|
| **v1** | `1800 / (TDD · ln(target/divisor + 1))` | ISF ∝ 1/TDD |
| **v2** | `2300 / (ln(target/divisor + 1) · TDD² · 0.02)` | ISF ∝ 1/TDD² |

Compared at the same glucose and TDD, the shared terms cancel **in the ratio between the
two equations** (within each equation the glucose scaler applies fully), giving
`ISF_v2 / ISF_v1 = 63.9 / TDD` — a pure function of TDD with a crossover at ~64 U/day.

## Headline results

- 77% of the 171-user cohort sits below the crossover: v2 computes weaker corrections for
  most people (≈3× weaker at 20 U/day), stronger only for the heaviest insulin users.
- Sensitivity calculated from each person's own data follows **ISF ∝ TDD^−0.4…−0.56** —
  shallower than v1's −1 and far from v2's −2. v2 is the worst-fitting of every equation
  tested against calculated sensitivity.
- Best simple equation across all candidates (leave-one-user-out cross-validation):
  **ISF ≈ K/√TDD** — K=355 anchored to tuned-profile ISF, K=145 anchored to calculated
  sensitivity (the anchor choice is a safety decision; see the proposal document).
- The v1 implementation reproduces device-logged ISF for all dynamic-ISF users, to within
  unmodelled per-person adjustment factor / divisor / autosensitivity.

## Documents (`docs/`)

| document | what it is |
|---|---|
| `dynamic-isf-methodology.md` | step-by-step methodology and reasoning (start here) |
| `dynamic-isf-v1-v2-analysis.md` | the v1-vs-v2 comparison results |
| `dynamic-isf-vnext-proposal.md` | the K/√TDD proposal for v-next + path to deployment |

Figures in `charts/inv008/`; candidate-equation and device-validation results in `results/`.

## Pipeline (`inv008/` package)

Per-user parallel replay (resumable, atomic per-person outputs):

```
python -m inv008.runner --stage 1 --platforms v6 v7    # delivery records → windowed TDD
python -m inv008.runner --stage 2 --platforms v5 v6 v7 # glucose readings → ISF under v1 & v2
python -m inv008.stage3_plots                          # per-person pages + cohort figures
python fit_best_isf.py                                 # cross-validated equation comparison
python -m inv008.validate_device_isf                   # implementation vs device-logged ISF
python -m pytest inv008/tests/                         # 18 unit tests
```

| module | role |
|---|---|
| `inv008/dynisf.py` | v1/v2 equations + TDD blend, vectorised, unit-tested |
| `inv008/tdd_windows.py` | delivery records → 5-min grid → the five TDD windows |
| `inv008/sources.py` | raw delivery/profile adapters + absolute-time anchor recovery |
| `inv008/stage1_tdd.py` | per-person windowed-TDD reconstruction worker |
| `inv008/stage2_replay.py` | per-person ISF replay worker (anchor validation, flat-TDD arm) |
| `inv008/stage3_plots.py` | per-person pages + cohort figures |
| `inv008/runner.py` | multiprocessing orchestrator (resume, logging, manifests) |
| `inv008/validate_device_isf.py` | replayed v1 ISF vs device-logged ISF (implementation check) |
| `fit_best_isf.py` | cross-validated comparison of candidate ISF equations |

Supporting extraction scripts (produce the inputs the pipeline expects):
`extract_treatments_tdd.py`, `extract_hourly_basal.py`, `canonical_cohort.py`,
`canonical_walsh.py`, `empirical_isf_v5.py`.

## Data

**No participant data is included in this repository.** The pipeline expects:

- a local time-series database with per-person decision tables extracted from public
  Nightscout samples and the OpenAPS Data Commons (anonymised),
- the raw delivery/profile archives on disk for TDD reconstruction,
- derived JSON inputs (cohort, calculated sensitivity, basal profiles, delivery totals,
  user mappings).

Point the pipeline at your local data tree with the `DYNISF_ROOT` environment variable
(it defaults to the current working directory):

```
export DYNISF_ROOT=/path/to/your/data
```

Only cohort-level figures and aggregate results are committed; one example per-person page
is included (anonymised, from the public OpenAPS Data Commons). The committed result JSONs
contain only anonymised participant IDs and derived statistics — no glucose timeseries or
identifiers.

## Caveats

This is a retrospective, decision-level analysis of equations — not a closed-loop outcome
study, and nothing here is dosing advice. See the caveats sections in the analysis and
proposal documents.

## Licence

- **Code** (everything except `docs/`) — MIT, see [`LICENSE`](LICENSE).
- **Documentation** (`docs/`) — Creative Commons Attribution 4.0 (CC BY 4.0), see
  [`docs/LICENSE`](docs/LICENSE).

Data provenance: derived from the OpenAPS Data Commons and publicly shared, de-identified
Nightscout datasets; only aggregate/derived results are redistributed here.
