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
10. Cohorts are reported separately and never pooled in absolute units: each JAEB study,
    Trio, oref0 and AndroidAPS. A constant is stated with its population, since the pooled
    JAEB record is 51% children.
11. Withdrawn figures stay withdrawn: 880, 355, the square-root TDD law and the 0.41 ratio.

## Running

12. `DYNISF_ROOT` moves the INV-009 results directory as well as the INV-008 inputs. Set
    for an INV-009 or INV-012 run, it writes results into another repository and the
    per-person Loop insulin model falls back without a message (10 November 2026). Pass
    paths as arguments instead.
13. Seven workers, a `__main__` guard, and every worker input in the job tuple. Time one
    worker before a cohort run. A stage that finishes far faster than expected, or two
    arms that agree to every decimal, is a fault until shown otherwise.

## Writing

14. Person-first language. In type 1 diabetes a raised glucose does not fall without
    insulin, so a fall with no correction is described as basal and earlier doses acting
    (correction of 27 September 2026).
15. Third-party slides and figures whose redistribution status is unknown are not
    committed.
