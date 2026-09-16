# dynamic-isf-calculations

Analysis of equations that set a dynamic insulin sensitivity factor (ISF) from total daily
dose (TDD) and current glucose, first evaluated in 171 people using open-source automated
insulin delivery (AID) systems and subsequently tested in public JAEB cohorts.

Dynamic ISF makes correction sensitivity a function of TDD and current glucose. The original
equation (**v1**, Chris Wilson) makes the sensitivity anchor inversely proportional to TDD; a
later revision (**v2**) makes it inversely proportional to TDD squared. The repository replays
both equations over recorded histories and compares their predictions with static and
loop-calibrated alternatives.

An audit completed on 16 September 2026 found that the historical change-in-IOB and
action-balance estimators are observational outcome proxies, not independent physiological
ISF measurements. Claims that physiological ISF follows `1/√TDD`, and the associated live
v-next dosing proposal, are therefore withdrawn. The equation replay and within-model
prediction comparisons remain reproducible. See
[`docs/DYNAMIC-ISF-AUDIT-2026-09.md`](docs/DYNAMIC-ISF-AUDIT-2026-09.md).

## The equations

Both use the same TDD-blending step (five windows → a weighted TDD), but they differ in both
the TDD term and the glucose logarithm:

| | sensitivity anchor at normal target | implied law |
|---|---|---|
| **v1** | `1800 / (TDD · ln(target/divisor + 1))` | ISF ∝ 1/TDD |
| **v2** | `2300 / (ln(target/divisor) · TDD² · 0.02)` | ISF ∝ 1/TDD² |

v1 keeps a `+1` in its glucose log and v2 does not. V2 floors glucose at `divisor+1` so the
log remains positive. For the standard rapid-acting configuration, including NovoRapid, the
divisor/floor is 75/76 mg/dL. The alternative configuration uses 55/56. The ratio between the
equations depends on glucose as well as TDD.

## Audited headline results

- v2 computes a weaker correction than v1 for almost everyone — on 92% of readings, a median
  of 3× weaker, most markedly at low glucose. At target glucose the two equations would only
  meet near 194 U/day, beyond anyone in the cohort.
- The historical change-in-IOB proxy scales as approximately **TDD^−0.4…−0.56**. That is a
  reproducible description of the proxy, not evidence that physiological ISF follows the same
  law.
- In selected overnight windows, a tuned static ISF and the loop's own calibrated prediction
  have lower error than v1 or v2. This comparison remains inside the loop's linear IOB model.
- With contiguous time-fold validation, the best glucose-shape exponent is `k=0`; adding a
  glucose-dependent sensitivity multiplier does not improve median prediction error.
- Post-COB analyses are consistent with carbohydrate absorption sometimes continuing after the
  controller's recorded COB reaches zero. This is a credible explanation for part of the
  apparent benefit attributed to glucose-dependent ISF, but it is not proof of a universal
  mechanism.
- The v1 implementation reproduces device-logged ISF for all dynamic-ISF users, to within
  unmodelled per-person adjustment factor / divisor / autosensitivity.
- The historical `K/√TDD` candidate and its constants remain available for research comparison,
  but are not physiological constants or dosing recommendations.

## Documents (`docs/`)

| document | what it is |
|---|---|
| `dynamic-isf-methodology.md` | step-by-step methodology and reasoning (start here) |
| `dynamic-isf-v1-v2-analysis.md` | the v1-vs-v2 comparison results |
| `dynamic-isf-data-derived-findings.md` | can sensitivity be derived from data? feasibility findings |
| `dynamic-isf-vnext-proposal.md` | withdrawn v-next proposal and current evidence status |
| `DYNAMIC-ISF-AUDIT-2026-09.md` | code, estimator and physiological-evidence audit |

Figures in `charts/inv008/`; candidate-equation and device-validation results in `results/`.

## Pipeline (`inv008/` package)

Per-user parallel replay (resumable, atomic per-person outputs):

```
python -m inv008.runner --stage 1 --platforms v6 v7    # delivery records → windowed TDD
python -m inv008.runner --stage 2 --platforms v5 v6 v7 # glucose readings → ISF under v1 & v2
python -m inv008.stage3_plots                          # per-person pages + cohort figures
python fit_best_isf.py                                 # cross-validated equation comparison
python -m inv008.validate_device_isf                   # implementation vs device-logged ISF
python -m pytest inv008/tests/                         # 27 unit tests
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

This is a retrospective, decision-level analysis of equations, not a prospective closed-loop
outcome study. Closed-loop delivery is chosen in response to glucose, while carbohydrate,
basal need and endogenous glucose are incompletely observed. The fitted observational
coefficients must therefore not be treated as measured physiological ISF or used as dosing
defaults. Nothing here is dosing advice.

## Licence

- **Code** (everything except `docs/`) — MIT, see [`LICENSE`](LICENSE).
- **Documentation** (`docs/`) — Creative Commons Attribution 4.0 (CC BY 4.0), see
  [`docs/LICENSE`](docs/LICENSE).

Data provenance: derived from the OpenAPS Data Commons and publicly shared, de-identified
Nightscout datasets; only aggregate/derived results are redistributed here.
