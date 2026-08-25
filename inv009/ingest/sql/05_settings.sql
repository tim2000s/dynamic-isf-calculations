-- Pump settings and bolus-calculator records: what people had their pumps set to.
--
-- Not a BabelBetes output. BabelBetes emits five event streams and stops there,
-- which is the right scope for it but leaves out the one thing a study of
-- insulin sensitivity most wants: the sensitivity factor the person had entered.
-- Those numbers are in the raw archives, unloaded, in two different shapes.
--
-- The bolus calculator files (Loop, REPLACE-BG) record what the calculator was
-- working from at each bolus: the sensitivity factor, carb ratio, targets, the
-- glucose and carbs the person keyed in, and in Loop's case the insulin on
-- board the app itself computed. That last column is the useful one twice over:
-- it is a settings record, and it is an independent check that a dose
-- reconstruction from the pump streams reproduces what the app believed.
--
-- The CRF settings files (DCLP5, PEDAP) record the whole pump programme at each
-- study visit: a half-hourly basal schedule and up to ten time segments of carb
-- ratio and correction factor.
--
-- UNITS. Loop and REPLACE-BG both store the calculator fields in mmol/L, which
-- is what Tidepool normalised to regardless of what the person saw on screen.
-- The loader converts per subject on the median, rather than per row, so a
-- single stray reading cannot flip a subject's units. DCLP5 and PEDAP come
-- straight off US CRFs and are already mg/dL.
--
-- DATES. ts_local follows the same convention as the event streams, so wizard
-- rows line up with studies.cgm and studies.bolus. The CRF visit dates do NOT:
-- they are the real clinic dates and the device data is shifted, so
-- pump_settings.therapy_dt cannot be joined to a device timestamp. Treat
-- pump_settings as a per-subject property, not a time series.
--
-- Run:  psql -d oref -v ON_ERROR_STOP=1 -f inv009/ingest/sql/05_settings.sql

BEGIN;

-- One row per bolus-calculator event. Loop and REPLACE-BG only.
CREATE TABLE IF NOT EXISTS studies.wizard (
    subject_id         text      NOT NULL,
    ts_local           timestamp NOT NULL,
    bolus_rec_id       text      NOT NULL,   -- '' when the file records no link
    bg_input_mgdl      real,                 -- glucose the person keyed in
    carb_input_g       real,
    iob_u              real,                 -- Loop only, and only ~60% of rows
    cr_g_per_u         real,
    isf_mgdl_per_u     real,
    target_low_mgdl    real,
    target_high_mgdl   real,
    rec_correction_u   real,
    rec_net_u          real,
    PRIMARY KEY (subject_id, ts_local, bolus_rec_id)
);

COMMENT ON TABLE studies.wizard IS
    'Bolus calculator records from study raw files (Loop LOOPDeviceWizard.txt, '
    'REPLACE-BG HDeviceWizard.txt). Glucose and sensitivity converted from '
    'mmol/L per subject. iob_u is the app''s own insulin-on-board and excludes '
    'the bolus being recommended. Loaded by inv009/ingest/load_settings.py.';

-- One row per pump-settings record per study visit. DCLP5 and PEDAP only.
CREATE TABLE IF NOT EXISTS studies.pump_settings (
    subject_id           text NOT NULL,
    visit_seq            int  NOT NULL,   -- order within subject as the file presents it
    visit_label          text,
    therapy_dt           date,            -- clinic date: NOT aligned to device timestamps
    tdd_reported_u       real,            -- total daily insulin as reported at the visit
    basal_reported_u     real,           -- DAILY basal total, see note below
    basal_hh             real[],          -- 48 half-hourly rates, forward filled
    basal_total_u        real,            -- sum(basal_hh)/2, a cross-check on the schedule
    isf_segments         jsonb,           -- [{start_min, end_min, value}, ...]
    cr_segments          jsonb,
    isf_tw_mgdl_per_u    real,            -- time weighted over the 24 hours
    cr_tw_g_per_u        real,
    PRIMARY KEY (subject_id, visit_seq)
);

COMMENT ON TABLE studies.pump_settings IS
    'Programmed pump settings per study visit, from CRF files (DCLP5 '
    'InsulinPumpSettings.txt, PEDAP PEDAPInsulinDeliveryDetails.txt). Already '
    'mg/dL. therapy_dt is a real clinic date and does not align with the '
    'date-shifted device streams: use this per subject, not as a time series. '
    'basal_reported_u comes from a raw column named TotBasalPreced7Days, but it '
    'is a DAILY total, not a weekly one: it agrees with the half-hourly schedule '
    'summed over 24 hours (r=0.92 DCLP5, 0.89 PEDAP, median ratio 0.99).';

-- Loop's scheduled basal rate, recovered from what each temp basal suppressed.
CREATE TABLE IF NOT EXISTS studies.basal_sched (
    subject_id      text      NOT NULL,
    ts_local        timestamp NOT NULL,
    sched_rate_u_hr real
);

COMMENT ON TABLE studies.basal_sched IS
    'Scheduled (profile) basal rate for Loop subjects, taken from the SuprRate '
    'column of temp basal rows in LOOPDeviceBasal*.txt. Loop''s insulin on '
    'board is net of scheduled basal, so reproducing it needs this rather than '
    'total delivery. Observed only while a temp basal was running, so it is a '
    'sampled schedule: forward fill it.';

CREATE INDEX IF NOT EXISTS wizard_subject_ts_idx
    ON studies.wizard (subject_id, ts_local DESC);
CREATE INDEX IF NOT EXISTS basal_sched_subject_ts_idx
    ON studies.basal_sched (subject_id, ts_local DESC);

COMMIT;

\echo 'studies.wizard, studies.pump_settings, studies.basal_sched created. Load with inv009/ingest/load_settings.py.'
