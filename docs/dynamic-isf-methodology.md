# Dynamic ISF from total daily dose — evaluation methodology

**2026-06-07** · Tim Street / Claude

---

## Abstract

This paper describes a reproducible computational method for evaluating equations that set
a dynamic insulin sensitivity factor (ISF) as a function of total daily dose (TDD). The
method (1) implements each candidate equation exactly, with a unit-test contract; (2)
reconstructs, from each person's recorded insulin-delivery history, the precise TDD signal
the equations consume; (3) replays every glucose reading through each equation to generate
the ISF it would produce; (4) compares the resulting ISF distributions between equations;
(5) searches a family of candidate equations for the best fit to independently-calculated
sensitivity using leave-one-user-out cross-validation; and (6) validates the equation
implementations against ISF values that devices themselves logged. We applied it to 171
people using open-source automated insulin delivery (AID) systems and roughly 9 million
glucose readings. This paper documents each step, the script that performs it, and the
reasoning that connects the outputs to the conclusions.

---

## 1. Background and objective

Dynamic ISF makes a person's correction sensitivity a function of two things: their total
daily dose of insulin, and their current glucose. The original formulation, due to Chris
Wilson, sets the sensitivity anchor in inverse proportion to TDD (the classic "1800 rule"
shape) and then scales it logarithmically with glucose. We refer to this as **v1**. A
later revision of the maths, also by Chris Wilson, makes the anchor inversely proportional
to **TDD squared**. We refer to this as **v2**.

Both forms share the same TDD-blending step and the same glucose scaler; they differ only
in how strongly TDD drives the sensitivity anchor. The questions this method answers:

1. How differently do v1 and v2 dose, across a real population, as a function of TDD?
2. Does either TDD power law match sensitivity as it is actually observed in user data?
3. Is there a better-fitting equation among a wider candidate family?

The scope is deliberately narrow: the **between-person sensitivity anchor** as a function
of TDD. The within-person glucose scaler is identical across the equations and is held
fixed throughout. This is a retrospective, decision-level analysis of equations; it is not
a closed-loop outcome study and provides no dosing advice.

---

## 2. The equations

Both v1 and v2 first derive a blended TDD, then a sensitivity anchor at "normal target",
then scale it by current glucose.

**TDD blend** (shared). From five windows of recorded delivery:

```
W8H = (1.4·TDD_4h + 0.6·TDD_8-4h) · 3
if W8H < 0.75·TDD_7d:
    adj7D = W8H + (W8H/TDD_7d)·(TDD_7d − W8H)
    TDD   = 0.34·adj7D + 0.33·TDD_1d + 0.33·W8H
else:
    TDD   = 0.33·W8H + 0.34·TDD_7d + 0.33·TDD_1d
```

**Glucose scaling** (shared). With glucose capped at 210 mg/dL (excess at one-third
weight), divisor 75 (for a rapid analogue such as Lyumjev), and normal target 99 mg/dL:

```
scaler = ln(target/divisor + 1) / ln(bg_capped/divisor + 1)
```

**Sensitivity anchor** (the only difference):

| | anchor at normal target | implied law |
|---|---|---|
| **v1** | `1800 / (TDD · ln(target/divisor + 1))` | ISF ∝ 1/TDD |
| **v2** | `2300 / (ln(target/divisor + 1) · TDD² · 0.02)` | ISF ∝ 1/TDD² |

The ISF used at any glucose is `anchor · scaler`. Writing both out in long form, with the
shared logarithmic terms collapsed:

```
v1:   ISF(BG) = 1800    / ( TDD   · ln(bg_capped/75 + 1) )
v2:   ISF(BG) = 115 000 / ( TDD²  · ln(bg_capped/75 + 1) )
```

A decisive analytic observation drives the whole method. **Comparing the two equations at
the same glucose and the same TDD, the glucose term and the divisor cancel in the ratio
between them** (not within either equation, where the scaler applies fully):

```
ISF_v2 / ISF_v1 = 2300 / (0.02 · 1800 · TDD) = 63.9 / TDD
```

The entire behavioural difference between v1 and v2 is therefore a pure function of TDD,
with a crossover at TDD ≈ 64 U/day. This is precisely why the method centres on
reconstructing TDD accurately and replaying it.

---

## 3. Data

Three anonymised cohorts of open-source AID users, held in a local time-series database:

| Platform | Users | Readings | TDD source |
|---|---|---|---|
| Trio | 23 | 3.03 M | device-logged TDD |
| AAPS | 39 | 1.31 M | reconstructed from treatments |
| OpenAPS | 109 | 6.58 M | reconstructed from treatments |

Each row is one 5-minute decision cycle carrying the glucose reading, the hour of day, and
an anonymised relative timestamp (seconds since that person's first record). Raw delivery
and profile archives (boluses, temp basals, basal schedules) remain on disk for the TDD
reconstruction. No participant-level data is published; only cohort-level aggregates and
figures are released.

---

## 4. Method, step by step

The pipeline is a Python package run in three stages plus two analyses. It is parallelised
per person (one person = one task), resumable, with atomic per-person outputs and
timestamped run manifests, and completes end to end in under ten minutes on a desktop
workstation.

### 4.1 Equation implementation and its test contract

The v1 and v2 equations and the TDD blend are implemented as vectorised array operations,
preserving every constant, the glucose cap, and the blend's branch logic. Correctness is
fixed by **18 unit tests** against hand-computed fixtures: the closed-form 63.9/TDD
identity, the glucose-cap compression, both branches of the TDD blend, and the
missing-data gates. *Why:* every downstream conclusion rests on these being exactly the
equations as defined; the tests are the contract that guarantees it.

### 4.2 TDD reconstruction

Where the device did not log TDD, we rebuild the five windows the blend consumes from raw
delivery records:

1. **Parse delivery events:** boluses (any record carrying insulin) and temp-basal
   segments (absolute rate or percent-of-profile, with supersession and cancellation),
   from each person's archive. Where an export has no temp-basal records, basal delivery
   uses the person's profile schedule (flagged).
2. **Build a 5-minute delivery grid:** lay the basal schedule across time, overlay
   temp-basal segments, and add boluses into their bins — units delivered per bin.
3. **Compute the windows:** trailing 4 h, the 8–4 h window, trailing 24 h, yesterday's
   calendar total, and the 7-day average. A calendar day counts toward the 7-day average
   only if it is complete on the grid and carries at least one bolus, so upload gaps do
   not silently deflate it.
4. **Blend:** apply the exact W8H formula; emit a missing value (the real-world fallback
   to profile ISF) wherever a component is absent.

*Why this matters:* the equations are functions of the *blended* TDD, not a flat daily
average. Reconstructing the actual windowed signal each device would have computed is what
makes the replay faithful rather than approximate. A flat-TDD arm is carried alongside, to
measure how much the blend itself contributes.

### 4.3 Absolute-time recovery and ISF replay

Reconstructed TDD lives in absolute time; decision rows carry only relative time. To join
them we recover each person's absolute-time anchor and **validate it against an independent
signal — the recorded hour-of-day**. Two candidate anchors are each refined by a whole-day
shift and then a fine sweep, scored by how well the implied hour-of-day matches the record;
the candidate giving the best join coverage wins, and residual mismatch is flagged when it
exceeds tolerance.

Then, for every glucose reading, the replay computes the ISF under both equations using
that reading's TDD (device-logged where available, reconstructed otherwise). The output is
a per-person table of (glucose, TDD, anchor, ISF) under each equation. Median join coverage
was 99.4%; 5 of 148 reconstructed users retain an uncertain anchor. The closed-form ratio
is reproduced in the replayed output (predicted 1.885 vs observed 1.882 for one
spot-check), confirming the join is sound.

### 4.4 Figures

Each person gets a page: the ISF–glucose curves under each equation at their median TDD
(with a TDD-interquartile band), a two-week sample of both dynamic-ISF traces over their
real glucose, and the per-reading ratio distribution. Cohort figures: the observed
ratio-versus-TDD crossover against the 63.9/TDD theory curve; per-person median ISF under
each equation; and the log-log ISF–TDD relationship with the observed points and fitted
slope.

### 4.5 Best-fit equation search

To ask whether a better equation exists, candidates are scored against two independent
targets:

- **calculated sensitivity** — ISF derived directly from each person's own data by a
  regression on fasting windows (`ΔBG = a + b·ΔIOB + c·BG_trend`; sensitivity = −b, the
  observed glucose drop per unit of insulin absorbed; n = 114 after quality gates). This
  is the ground truth for "what insulin actually does".
- **tuned-profile ISF** — each person's own profile setting (n = 138): "what experienced
  users converge to".

Candidates: the two equations above; the historical 1700-rule; the tuned profile value;
and fitted forms — a re-fitted constant/TDD, a free-exponent power law A·TDD^b, a
fixed-exponent square-root rule K/√TDD, TDD-quartile bands, multivariate log-linear models
adding carb ratio / basal fraction / target, and a geometric blend of the profile value
with the power law.

**Every fitted candidate is scored by leave-one-user-out cross-validation** — refit on all
but one person, predict the held-out person, repeat — so fitted forms are judged
out-of-sample and are directly comparable with the fixed rules. Metrics: median absolute
error, median absolute log-error (scale-free), and fraction of people predicted within
±30%. *Why cross-validation:* without it, fitted equations would carry an unfair in-sample
advantage; this measures genuine generalisation to a new person.

### 4.6 Implementation validation against device-logged ISF

Independently of the equation comparison, we check that the v1 implementation reproduces
what devices actually computed. Trio logs its own per-cycle ISF, giving ground truth for
its dynamic-ISF users. We correlate replayed v1 ISF against the device value, expecting,
per person, a positive log-log correlation (the curve shape tracks) and a tight, roughly
constant multiplicative offset — the device additionally applies an adjustment factor, an
insulin divisor, and an autosensitivity ratio, none of which the replay models — rather
than a ratio of exactly one.

This step also surfaced a data-quality issue: four users switched the *units* of their
logged ISF mid-history (mmol/L per unit vs mg/dL per unit). A per-reading correction
(values below 20 scaled by 18.018) resolves all four. It affects only this validation read;
the main pipeline never consumes that field.

---

## 5. How the outputs support the conclusions

**The v1/v2 difference is purely TDD, and v2 is weaker for most people.** The 63.9/TDD
ratio is exact algebra (§2) and is reproduced in the replayed data (§4.3). The crossover
figure shows per-person median ratios lying on that curve across the full TDD range, with
the majority of users below the 64 U/day crossover. The conclusion — that moving from v1
to v2 weakens corrections for most people and strengthens them only for the heaviest — is
therefore not a modelling artefact but a direct consequence of exact algebra confirmed in
data.

**Both TDD power laws are too steep.** The log-log fit of independently-calculated
sensitivity against reconstructed TDD has a slope near −0.5 (bootstrap interval excluding
−1). v1 assumes −1; v2 assumes −2. Because the target is sensitivity *calculated from each
person's own glucose and insulin data* — not profile settings, not the equations
themselves — this is an independent test of the TDD exponent, and it rejects both.

**A square-root law fits best.** Under cross-validation the K/√TDD form is best against
tuned-profile ISF and tied-best against calculated sensitivity, beating both equations and
every alternative. Added inputs (carb ratio, target, basal fraction) do not robustly
improve it. The two natural anchors differ by a constant factor, establishing that the
*shape* of the relationship is settled by the data while the *level* is a separate choice.

**The implementation is faithful.** The unit tests fix the implementation to the defined
equations, and the device-ISF validation shows the dynamic-ISF users tracking their
devices' logged ISF with stable per-person offsets once the units issue is corrected, with
no unexplained outliers. This licenses treating the replayed ISF as the equation's true
output.

---

## 6. Reproducibility

| Step | Script | Key output |
|---|---|---|
| Equation implementation + tests | `inv008/dynisf.py`, `inv008/tests/` | 18 passing tests |
| TDD reconstruction | `inv008/stage1_tdd.py` (+ `tdd_windows.py`, `sources.py`) | per-person TDD tables |
| ISF replay | `inv008/stage2_replay.py` | per-person ISF tables |
| Orchestration | `inv008/runner.py` | run logs + manifests |
| Figures | `inv008/stage3_plots.py` | cohort + per-person figures |
| Equation search | `fit_best_isf.py` | candidate-comparison tables |
| Device validation | `inv008/validate_device_isf.py` | validation tables |
| Delivery/basal inputs | `extract_treatments_tdd.py`, `extract_hourly_basal.py` | input data |
| Reference cohort + calculated sensitivity | `canonical_cohort.py`, `empirical_isf_v5.py` | targets for the search |

```
python -m pytest inv008/tests/
python -m inv008.runner --stage 1 --platforms v6 v7
python -m inv008.runner --stage 2 --platforms v5 v6 v7
python -m inv008.stage3_plots
python fit_best_isf.py
python -m inv008.validate_device_isf
```

---

## 7. Limitations

1. **Counterfactual replay.** We compare the ISF each equation *would have computed* on
   real histories, not closed-loop glycaemic outcomes.
2. **Basal approximation for one platform.** Where exports lack temp-basal records, basal
   TDD uses the profile schedule.
3. **Calculated-sensitivity level.** The regression estimator may be biased low by
   unrecorded carbohydrate or endogenous-glucose effects; it has confidence intervals but
   no external ground truth. The *shape* it implies is robust; its absolute *level* is
   provisional.
4. **Single cohort.** Open-source AID users, mostly 2016–2023; no commercial-system or
   demographic data. n = 114/138 for the calculated-sensitivity analyses.
5. **Anchor uncertainty.** Five reconstructed users retain time-anchor uncertainty and are
   flagged.
