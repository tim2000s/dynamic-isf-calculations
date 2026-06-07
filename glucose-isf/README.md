# glucose-isf — the glucose dimension g(BG)

Prior work (Mar–Apr 2026) on **glucose-dependent ISF** — the g(BG) half of
`ISF = f(TDD) × g(BG)`. The rest of this repository covers the f(TDD) half (the √TDD law);
this folder covers how ISF should vary with **current glucose**.

## What it establishes

- **Glucose-dependent ISF is real**, from the **Diabeloop** DBLG1 clinical *population* model
  (ADA scientific poster): ISF vs glucose as a piecewise polynomial (quadratic ≤100 mg/dL,
  quartic >100), with proportional IQR (~×0.5/×1.6 of median) — supporting a single
  multiplicative per-patient scaling factor. `Variable ISF.png` / `Diabeloop.png` are the model.
- **A power-law glucose term beats the logarithmic scaler.** `ISF = (C/TDD)·(target/BG)^k`,
  k≈3.5, cut 2-hour glucose-prediction error 12–18% vs `ln(BG/D+1)`, and is hypo-protective by
  shape (much higher ISF at low BG → less insulin where it's dangerous). See
  `DynamicISF_PowerLaw_Analysis.md`, `Hybrid_ISF_Analysis.md`.
- **A no-TDD hybrid** (Diabeloop quartic above 105 + power-law tail below) and a per-patient
  **calibration** strategy (scaling factor from TDD/profile day-1, refined from overnight
  prediction errors). `DynamicISF_TDD_Effective_Analysis.md`.

Evidence base: **N=1** (10 months of one patient's overnight closed-loop data) + synthetic
multi-patient validation + the Diabeloop population curve. Strong on the glucose **shape**; the
TDD **exponent** was not tested across real patients (that is the cohort work elsewhere in this
repo — which finds ≈√TDD, not the 1/TDD these papers assumed).

## How it reconciles with the cohort work (see `results/phase8_glucose_tdd.md`)

- **TDD exponent:** the oref cohort rejects 1/TDD in favour of ≈√TDD; use the cohort value.
- **Glucose term:** the cohort's *observational* per-window ISF regression **cannot** validate the
  glucose shape — the unannounced-carb confound flips the slope (apparent ISF rises with BG). So
  the power-law rests on the **clean prediction-error backtest** here, not the cohort regression.
  Validating it at scale needs a forward-prediction backtest (this design), not an ISF regression.

## Key files

`DynamicISF_PowerLaw_Analysis.md`, `Hybrid_ISF_Analysis.md`, `DynamicISF_TDD_Effective_Analysis.md`
(papers); `ns_bg_scaling*.py`, `ns_polynomial_backtest.py`, `ns_combined_powerlaw_tdd.py`,
`ns_overnight_backtest.py`, `ns_prediction_performance.py`, `empirical_isf_curves.py`,
`adaptive_isf_model.py`, `replay_simulation.py` (analysis).

## Data

**No raw participant data is included** (consistent with the rest of the repo): the N=1 CGM/insulin
CSVs, caches, and generated Word documents live only on the author's Drive. The Diabeloop model is
a published population curve. Scripts expect the local data to be supplied.

Provenance: originally in `Dynamic ISF data/` on the author's Drive; merged here so the f(TDD) and
g(BG) halves of the dynamic-ISF work live together.
