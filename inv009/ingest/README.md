# Pump settings and bolus calculator records

BabelBetes emits five event streams (cgm, basal, bolus, carbs, age) and stops
there. That is the right scope for it, but it leaves out the number a study of
insulin sensitivity most wants: the sensitivity factor the person actually had
entered. Those numbers sit unloaded in the raw JAEB archives.

This loader adds three tables to the `studies` schema of the `oref` database. It
does not change anything BabelBetes wrote.

```
psql -d oref -v ON_ERROR_STOP=1 -f inv009/ingest/sql/05_settings.sql
python3 inv009/ingest/load_settings.py            # about 40 seconds end to end
```

It is idempotent. Every source stages through a temp table and upserts, so a
re-run replaces rather than collides.

## What is loaded

| Table | Source | Subjects | Rows |
|---|---|---|---|
| `studies.wizard` | Loop `LOOPDeviceWizard.txt` | 201 | 103,990 |
| `studies.wizard` | REPLACE-BG `HDeviceWizard.txt` | 192 | 218,900 |
| `studies.pump_settings` | DCLP5 `InsulinPumpSettings.txt` | 101 | 1,553 visits |
| `studies.pump_settings` | PEDAP `PEDAPInsulinDeliveryDetails.txt` | 102 | 1,569 visits |
| `studies.basal_sched` | Loop `LOOPDeviceBasal{1,2,3}.txt` | 842 | 1,627,707 changes |

596 people with a recorded sensitivity factor. Median entered sensitivity is
65 mg/dL/U in Loop, 50 in REPLACE-BG, 68 in DCLP5 and 165 in PEDAP, which is the
order the age of those cohorts would lead you to expect.

DCLP3 and IOBP2 ship no settings file at all, so neither appears here. IOBP2 is
the bionic pancreas, which by design has no sensitivity factor to record.

## Four things worth knowing before using it

**Units.** Loop and REPLACE-BG store the calculator fields in mmol/L whatever
the person saw on screen, and the archives give no unit column. The loader
decides per subject on the median entered sensitivity, not per row, so a single
stray value cannot flip a subject. 414 of the 419 subjects are mmol; the other
five recorded no sensitivity at all. DCLP5 and PEDAP come straight off US case
report forms and are already mg/dL. Carb ratio is grams per unit in both systems
and is never converted.

**Timestamps line up, and this was checked.** Wizard rows are put on the same
clock as the event streams by the same rule BabelBetes used: for Loop, UTC plus
the roster's fixed per-subject offset; for REPLACE-BG, days from a 2015-01-01
epoch. BabelBetes built `studies.carbs` for Loop out of this same wizard file, so
the two must agree, and they do: **all 80,039 Loop wizard rows carrying carbs
match a `studies.carbs` row exactly.** That is the check to re-run if the
timestamp handling is ever touched.

Rows for subjects BabelBetes excluded, or from outside the glucose record it
loaded, are dropped at the end of the load (28,793 rows). That enforces each
study's own inclusion window without re-deriving its exclusion rules.

**Visit dates do not align with device data.** `pump_settings.therapy_dt` is the
real clinic date from the case report form, while the device streams are
date-shifted. The two cannot be joined. Treat `pump_settings` as a per-subject
property and not as a time series.

**`basal_reported_u` is a daily total.** It comes from a raw column named
`TotBasalPreced7Days`, which reads as a weekly figure and is not one. Summing the
half-hourly schedule over 24 hours agrees with it (r = 0.92 in DCLP5, 0.89 in
PEDAP, median ratio 0.99), so the column is daily and the name is misleading.
That agreement is also what validates the schedule parser, because the CRF writes
a half hour only when the rate changes and leaves the rest blank: read the blanks
as zero and everyone's basal halves. They are forward filled, wrapping from the
end of the day when midnight itself is blank.

## Why `studies.basal_sched` exists

Loop computes insulin on board net of scheduled basal: a temp basal contributes
only the difference between what it delivered and what the profile would have
delivered. Reproducing what the app believed therefore needs the schedule, and
the delivered stream alone cannot give it.

The only record of the schedule is the rate each temp basal suppressed, in the
`SuprRate` column. It is observed only while a temp was running, and it is a step
function that changes a handful of times a day, so 6.7 GB of raw text collapses
to 1.6 million changes. **Forward fill it.**

This is what makes `inv009/loop_model_infer.py` possible: with the schedule in
hand, insulin on board can be recomputed under each candidate insulin model and
compared against the value Loop itself recorded, which both identifies the model
each person was using and checks the dose reconstruction against ground truth.
