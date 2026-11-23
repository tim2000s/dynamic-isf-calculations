# Review checklist

Read before committing. Each item is a way this repository has gone wrong, written as a
check a diff can be tested against. Findings are labelled important (fix before commit),
nit (optional) or pre-existing (record it, do not let it block).

## Insulin exposure

1. Net of schedule. A dose or a correction is insulin above the programmed basal, by any
   route. An estimator keyed on `bolus_u` alone misses algorithm basal and is wrong on
   closed-loop records.
2. No double count. A dose landing in bin T is inside IOB(T) already, so the delivered term
   starts after T. `inv009/tests/test_action.py` and `inv012/tests/test_peak_trough.py`
   hold this.
3. Suggested is not delivered. `sug_rate` and `sug_smb_units` in the oref tables are
   suggestions, and the tables hold no manual boluses. State which is used and why.

## Estimators

4. Known truth first. A new sensitivity estimator is run on `inv009/synthetic.py` records
   before its real-data number is read, and the write-up says which glucose range the
   simulation covered.
   The simulator's own arms are checked too: until 10 November 2026 its unlogged meals
   were written into `carbs_g` and absorbed at twice their size, so every arm described
   as hiding meals hid none and no arm reached 150 mg/dL away from a meal.
   `inv009/tests/test_synthetic.py` now catches this. A simulated arm's glucose
   percentiles are looked at before its result is read.
5. Selection on the outcome. A window that starts at a high reading, or ends where glucose
   stopped falling, selects episodes by what glucose did. Its ratio is an outcome ratio
   and is described as one (audit of 16 September 2026).
6. Ratios over a denominator near zero. Report the share of events with no net insulin
   used alongside any median of fall per unit.
7. The expected sign and the modelling choices are written in the module docstring before
   the run, and every reported effect carries an interval and a verdict.

## Data

8. Repeated records. `oref_v5` holds 3,028,426 rows and 2,783,877 distinct
   (user, timestamp) pairs. Collapse before counting or integrating.
9. Units. `sug_ISF` in `oref_v5` switches between mmol/L and mg/dL within some users.
   Bolus durations are in milliseconds for three JAEB studies.
   The profile sensitivity in `canonical_cohort.json` is already mg/dL per U. A rule that
   converts values below 20 is checked against the person's daily dose before it is applied
   (INV-016, 12 November 2026: three values of 10 to 14 at 62 to 115 U a day were mg/dL).
10. Cohorts are reported separately and never pooled in absolute units: each JAEB study,
    Trio, oref0 and AndroidAPS. A constant is stated with its population, since the pooled
    JAEB record is 51% children.
11. Withdrawn figures stay withdrawn: 880, 355, the square-root TDD law and the 0.41 ratio.
17. A statement that a source lacks a record type is checked by listing the files. Until
    18 November 2026 `inv008/sources.py`, the audit and the paper said the AAPS v6 export
    holds no temporary-basal records; it holds `TemporaryBasals.json` for all 44 people and
    the reader did not open it. Since 23 November 2026 it does, and `inv008/tests/test_v6_temps.py`
    tests the reading; nothing yet catches an unread record type in general.
18. The same-windows claim. A sentence that says every candidate was scored on the same
    windows is checked against the window count of each arm. `results/head_to_head.md`,
    `results/divisor_sensitivity.md` and the title of `fig_head_to_head.png` say so, while
    V1 and V2 were scored on 57,078 of the 62,751 (57,151 of 62,824 since the rerun of
    23 November 2026). `results/common_window.md` holds the
    comparison on identical windows and is the one to quote.
19. The source tables move. On 18 November 2026 `inv008/overnight_sensitivity.py`, run
    unchanged, returned 81,769 windows against the 81,814 it had committed, with three Trio
    people's `oref_v5` rows changed in between. A variant run reproduces the committed
    count first and states both where they differ (`inv008/overnight_variants.py` does).

20. One record in several files is counted once. The AndroidAPS uploads and the oref0 dump
    files overlap in time, and `inv008/sources.py` read every copy of a bolus, and
    `extract_treatments_tdd.py`, behind canonical_cohort.json, did the same and also summed
    deleted AndroidAPS records (found 23 November 2026; `results/reader_correction.md`
    gives the size). Both now count each (time, units) once. The same day `parse_ts` was
    found cutting times at 25 characters, which dropped every record written with
    milliseconds and an offset. Tested in `inv008/tests/test_v6_temps.py`. Every result built
    on the daily doses, the cohort file or the replay cache was rerun the same day, together
    with a unit fix: ISF values under 20 had been read as mmol/L where they are mg/dL
    (`inv008/tests/test_units.py`). A new reader of these archives checks the share of unique records
    and of unparsed times before summing.

21. A bootstrap that resamples people numbers them in a fixed order. `inv008/overnight_variants.py`
    took people in the order its worker pool returned them, so the same seed drew different
    people on each run and its intervals moved by up to 0.3 mg/dL per U between identical runs
    (found 23 November 2026, fixed by sorting). A rerun whose point estimates match and whose
    intervals do not is this fault.

22. Per-person results are never published (decided 23 November 2026). Scripts write their
    full output, per-person parts included, to the private data folder or a git-ignored cache;
    `inv008/publish.py` copies a result into results/ with per-person records and rows removed
    and fails if a participant ID remains. A diff that adds a participant ID (U and three
    digits) outside a cohort's ID range, a per-person table, a per-person chart or a per-person
    page to the repository is refused. Boost data is never in this repository.

## Running

12. `DYNISF_ROOT` moves the INV-009 results directory as well as the INV-008 inputs. Set
    for an INV-009 or INV-012 run, it writes results into another repository and the
    per-person Loop insulin model falls back without a message (10 November 2026). Pass
    paths as arguments instead.
13. Seven workers, a `__main__` guard, and every worker input in the job tuple. Time one
    worker before a cohort run. A stage that finishes far faster than expected, or two
    arms that agree to every decimal, is a fault until shown otherwise.

16. A long run writes its raw output to disk before any reduction or formatting, and the
    whole script is run end to end on two people and two days before the cohort run. On
    11 November 2026 the INV-014 arms ran for 2 h 47 min and were lost when the final
    parquet write met a run that had ended early, whose reason TimSim records as a dict.
    The one-person pilot had no such run. `inv014/known_truth.py` now pickles the pool's
    output first.

## Writing

14. Person-first language. In type 1 diabetes a raised glucose does not fall without
    insulin, so a fall with no correction is described as basal and earlier doses acting
    (correction of 27 September 2026).
15. Third-party slides and figures whose redistribution status is unknown are not
    committed.
