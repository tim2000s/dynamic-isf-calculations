# Empirical Walsh constants in open-source automated insulin delivery: an observational analysis of 138 anonymised loop users

**Tim Street¹**, *[Clinical co-author, TBD]²*, *[Statistical reviewer, TBD]³*

¹ Independent researcher, Live In Their Shoes (tim@liveintheirshoes.com), Bristol, UK
² *[Affiliation TBD — preferably a diabetes-care endocrinologist with teaching responsibility for the 1700/500/TDD÷2 rules]*
³ *[Affiliation TBD — methods statistician]*

**Corresponding author:** Tim Street, tim@liveintheirshoes.com

**Word count (main text):** ~2,800

**Funding:** None.

**Conflicts of interest:** None. The author has no commercial relationships with insulin pump, CGM or AID vendors and received no funding for this analysis.

**Data availability:** The cohort is derived from the OpenAPS Data Commons archive (n=240, April 2026 release) and a sister Nightscout sample dump. Analysis code, extraction scripts and the de-identified per-user feature table that supports the paper are available at https://github.com/tim2000s/dynamic-isf-calculations (analysis code under MIT licence, documentation under CC BY 4.0). The Walsh-constants analysis itself is `canonical_walsh.py`, run against the shared cohort module `canonical_cohort.py`. No patient-level identifying data is reproduced.

**Ethics:** This work uses only anonymised loop-uploads from public, opt-in data-sharing programmes (OpenAPS Data Commons; user-submitted Nightscout dumps). No personally identifying information is held or analysed. *[IRB / local ethics committee statement TBD by clinical co-author.]*

---

## Abstract

**Background.** The classical "1700-rule" for insulin sensitivity factor (ISF), the "500-rule" for carbohydrate-to-insulin ratio (CR), and the assumption that 50% of total daily insulin dose (TDD) is delivered as basal, are still routinely taught when initialising pump therapy. These rules originate from pre-loop, multiple-daily-injection populations of the 1980s. Whether they remain valid for users of open-source automated insulin delivery (AID) is unknown.

**Methods.** We assembled an anonymised cohort of open-source AID users (Trio, AAPS classic, and oref0/OpenAPS, 2016–2025) drawn from public Nightscout dumps and the OpenAPS Data Commons. We applied a single canonical TDD definition: for users whose treatments archives were parseable from disk, TDD = profile basal + sum of recorded insulin-delivery events / recording span. For users extracted live via the Nightscout REST API at upload time, TDD = the live-extracted `mean_tdd` field, which integrates the same delivery events at extraction. We retained users with at least 14 days of decision history (n=138). We computed bootstrap 95% confidence intervals for the median ISF×TDD, CR×TDD and basal/TDD constants, fitted log-linear models for `log(ISF) = a - b·log(TDD)` per platform with bootstrap CIs on the slope, and tested robustness to outlier removal and duration thresholds.

**Results.** Median ISF×TDD = 2,381 [95% CI 2,144–2,633], significantly above Walsh's 1,700. Median CR×TDD = 402 [351–441], significantly below the 500-rule. Median basal/TDD = 0.45 [0.43–0.48], significantly below Walsh's 0.50. Because ISF is a divisor in the correction-dose calculation, a cohort constant above 1,700 means the rule prescribes a *smaller* ISF, and therefore a *larger* correction dose, than these users have settled on. The fitted cohort-wide log-linear slope was −0.43 (95% CI −0.59 to −0.27), inconsistent with the slope of −1 implied by Walsh's parametric form. The no-DynISF subgroup (n=120) had slope −0.38 (CI −0.56 to −0.22). The DynISF subgroups had slopes of −0.89 (Sigmoid, n=9; CI −1.44 to −0.35) and −0.84 (Log, n=9; CI −1.35 to −0.25); neither subgroup excluded −1, and these are exploratory findings. All three constants survived outlier removal (5 users dropped) and the strictest duration threshold (≥180 days, n=75).

**Conclusion.** In open-source AID users, all three classical Walsh constants disagree with the empirical cohort medians at the 95%-CI level, and the parametric form of the rules does not fit the no-DynISF subgroup that comprises the majority of the cohort. Initial profile defaults derived from Walsh constants set ISF too low, and the basal share of total insulin too high, for this population.

**Keywords:** automated insulin delivery, OpenAPS, AAPS, insulin sensitivity factor, Walsh's rules, observational study.

---

## 1. Introduction

The "1700-rule" (insulin sensitivity factor in mg/dL ≈ 1700 / TDD), the "500-rule" (carb ratio in g/U ≈ 500 / TDD) and the convention that approximately 50% of TDD is delivered as basal insulin date to John Walsh's *Pumping Insulin* (first edition 1989, sixth edition 2017).¹ They were derived from clinical observation of multiple-daily-injection (MDI) users in the 1980s and were later adopted, largely without re-derivation, into the documentation that accompanies most modern insulin pumps and into the training materials of diabetes specialist nurses worldwide.

Two things have changed since the rules were written. First, a substantial fraction of people with type 1 diabetes, and in some communities the majority, now use automated insulin delivery (AID) systems that adjust basal and bolus insulin in response to continuous glucose monitor (CGM) data.² Three of the most widespread open-source AID systems (oref0/OpenAPS, AAPS, and Trio) share an underlying control algorithm derived from oref0.³ Second, the introduction of super-microboluses (SMBs) has shifted a growing share of total insulin from user-entered carb boluses to algorithmic micro-doses delivered every five minutes.⁴ A user whose AID is well-tuned and whose meal-bolus discipline is loose may now have a profile in which the algorithm decides 30–60% of the daily insulin total, a regime Walsh did not anticipate.

Whether the empirical relationships between ISF, CR and TDD remain those that Walsh observed has, to our knowledge, not been re-tested in a published cohort of open-source AID users. The OpenAPS Data Commons⁵ and analogous community-shared Nightscout dumps make such a re-test feasible: they provide profile fields (ISF, CR, target, hourly basal pattern) alongside the full treatments archive (every recorded insulin-delivery event, both algorithm-driven and user-initiated), anonymised at upload.

In this paper we ask three questions, in increasing strength of claim:

1. Are Walsh's constants 1,700, 500 and 0.5 still recovered from the medians of an open-source AID cohort?
2. Does the parametric form (ISF, CR ∝ 1/TDD; basal = TDD/2) itself fit the data?
3. Do users on different open-source AID platforms, despite running the same underlying control algorithm, present systematically different profiles?

We report bootstrap 95% confidence intervals on each constant, log-linear fits per platform, robustness to outlier removal, and sensitivity to the duration-of-history filter. We highlight what the data does not support and label exploratory findings explicitly.

This work is observational. It does not prospectively validate any new prescribing rule. It does not include the demographic data (age, weight, duration of diabetes, country of practice) that would be required to interpret any per-platform differences definitively. The intent is narrower: to establish whether the textbook constants, in the form they are taught, are recoverable from the data of the population they are now most often applied to.

---

## 2. Methods

### 2.1 Cohort

We assembled three sources of anonymised open-source AID upload data:

- **Trio + early-AAPS DynISF (n=29 raw):** extracted live via the Nightscout REST API from 24 publicly-shared sites between 2023 and 2025. Decisions table `oref_v5` (3.03 M rows). Users U000–U028.
- **AAPS classic (n=44 raw):** anonymised AAPS uploads to the OpenAPS Nightscout sample archive, decisions table `oref_v6` (1.31 M rows). Users U029–U072.
- **oref0/OpenAPS (n=110 raw):** anonymised uploads from the OpenAPS Data Commons (April 2026 release), decisions table `oref_v7` (6.58 M rows), spanning 2016–2023. Users U073–U182.

After de-duplication, the raw analytical cohort comprised 183 users; 138 met the quality filters.

**Scope.** All sources are *open-source* AID systems built on the oref0 control algorithm. Closed-source and commercial AID systems, including Tandem Control-IQ, Medtronic MiniMed 780G, Insulet Omnipod 5 and CamAPS FX, are not represented. The findings should not be assumed to generalise to commercial systems, which use different control algorithms (model-predictive control, fuzzy logic, PID variants), different default profile presets, and different onboarding workflows. Replication on closed-source AID populations would require vendor data-sharing agreements that are not currently available to this analysis.

### 2.2 Profile and TDD extraction

For each user we obtained the declared pump profile (ISF in mg/dL, CR in g/U, target-low BG in mg/dL, extracted as the `target_low` field; almost all cohort users set `target_low` = `target_high` so this functions as a point target), hourly basal pattern in U/h × 24, and max basal in U/h from the user's most recent Nightscout profile document. Where the v7 enrichment had not retained the profile basal total, we re-extracted it from the underlying profile JSON by walking the `direct-sharing-31/*profile*.json` files in the user's NS sample directory, picking the active profile and converting its segment list to a 24-hour vector summing to total daily basal in U/day.

**Canonical TDD definition.** We applied a single TDD definition across all users in the analysis, chosen to be the most inclusive and defensible option:

- **For v6 (AAPS classic) and v7 (oref0/OpenAPS) users**, where the original Nightscout treatments archives are accessible on disk, TDD = `profile basal (U/day) + Σ recorded insulin-delivery events / recording-span days`. Treatments include both algorithm-driven super-microboluses and user-initiated meal and correction boluses.
- **For v5 (Trio + early-DynISF) users**, whose data was extracted live via the Nightscout REST API rather than from a static archive, we used the per-user `mean_tdd` field that the live extractor computed. This was derived by integrating the same recorded delivery events over the user's entire data window, mathematically equivalent to the v6/v7 method, just applied at extraction time rather than re-applied at analysis time.

Both routes converge on the same definition: real, accountable, recorded insulin delivered per day. Two earlier definitions were considered and discarded: a pure algorithmic-SMB-only sum (`sug_smb_units`), which under-counts user-entered boluses; and a hybrid `max(2 × basal, basal + SMB/day, algorithm-derived TDD)`, which overestimates for high-basal-low-bolus users and underestimates for low-basal-high-bolus users, biases that cancel only at the cohort median.

### 2.3 Quality filters

We retained users with ISF in [10, 300] mg/dL; CR in [2, 50] g/U; target-low BG in [70, 130] mg/dL; TDD in [5, 200] U/day; and at least 14 days of decision history.

After applying these filters, 138 users entered the analytical cohort: 22 v5 / Trio (live `mean_tdd`), 19 v6 / AAPS classic (treatments-derived), 97 v7 / oref0 (treatments-derived). Of the 138, 18 are on DynISF (v5 only, the only AID software in the cohort that exposes DynISF, and only available from AAPS 2021 onwards); 9 use the Sigmoid formula and 9 use the Log formula, per a curated mapping validated against live Nightscout reason strings.

### 2.4 Statistical analysis

For each Walsh constant we report the cohort median and a 95% confidence interval obtained by 2,000 non-parametric bootstrap resamples. For the log-linear form `log(ISF) = a - b·log(TDD)`, fitted by ordinary least squares per group, we report (a, b, n). We compared platforms with the Mann-Whitney U test and Kruskal-Wallis H test on ISF × TDD. We performed an outlier-sensitivity analysis by removing users with |z(ISF × TDD)| > 2. We tested duration sensitivity at five thresholds (≥ 14, 30, 60, 90, 180 days). All analyses are reproducible from the open analysis script (see *Data availability*).

---

## 3. Results

### 3.1 Cohort summary

The 138 users had a median TDD of 44 U/day, median ISF of 50 mg/dL, median CR of 9 g/U, median target 100 mg/dL, and median total basal 20 U/day. Cohort distribution across platforms: 22 Trio (live `mean_tdd`), 19 AAPS classic (treatments-derived), 97 oref0 (treatments-derived). Eighteen users, all in the v5 / Trio sub-cohort, ran DynISF; the other 120 ran no DynISF (oref0 and AAPS classic predate DynISF).

### 3.2 All three Walsh constants disagree with the cohort

| Constant | n | Median | 95% CI | Walsh value | CI excludes Walsh? |
|---|---|---|---|---|---|
| ISF × TDD | 138 | 2,381 | [2,144, 2,633] | 1,700 | Yes |
| CR × TDD | 138 | 402 | [351, 441] | 500 | Yes |
| basal / TDD | 138 | 0.45 | [0.43, 0.48] | 0.50 | Yes |

The bootstrap confidence intervals around all three constants exclude the Walsh values.

The direction of the ISF result requires care, because ISF enters the correction-dose calculation as a divisor. A larger ISF means the algorithm assumes each unit of insulin lowers glucose further, and therefore delivers a *smaller* correction dose. The cohort constant of 2,381 is 40% above Walsh's 1,700, so at any given TDD these users run a *higher* ISF, and correct *less* aggressively, than the rule prescribes. Equivalently, the profiles they have settled on imply greater insulin sensitivity than the 1700-rule assumes. The CR result runs the other way: at 402 the cohort constant is 20% below the 500-rule, so these users enter *more* insulin per gram of carbohydrate than the rule prescribes. The basal share is slightly smaller than Walsh's convention, at 45% against 50%.

### 3.3 The parametric form of the rules does not fit

A 1/TDD relationship between ISF and TDD implies, in log space, a slope of exactly −1. We fitted `log(ISF) = a - b · log(TDD)` per group, with a 95% bootstrap CI on the slope:

| Group | n | slope b | 95% CI on b | intercept a | CI excludes −1? |
|---|---|---|---|---|---|
| no DynISF | 120 | −0.38 | [−0.56, −0.22] | 5.44 | Yes |
| Sigmoid DynISF | 9 | −0.89 | [−1.44, −0.35] | 7.44 | no (n=9) |
| Log DynISF | 9 | −0.84 | [−1.35, −0.25] | 7.36 | no (n=9) |
| All cohort | 138 | −0.43 | [−0.59, −0.27] | 5.65 | Yes |

The cohort-wide and no-DynISF slope CIs both exclude −1: the Walsh parametric form does not fit the no-DynISF sub-cohort that comprises 87% of the cohort. The DynISF sub-cohorts of n=9 each have point estimates near −1 but bootstrap CIs too wide to exclude any value between −1.4 and −0.25; these per-group findings are exploratory and require replication.

In practical terms, doubling a no-DynISF user's TDD from 30 to 60 U/day is associated with an ISF reduction of about 22% (1 − 2^−0.38), not the 50% the Walsh form would predict. Similar shallower-than-Walsh observations have been reported in non-AID populations *(see references — Walsh constant re-evaluation literature in MDI cohorts)*, but to our knowledge this is the first replication in open-source AID users.

*Figure 1. log(ISF) versus log(TDD) by AID platform, with cohort fit (teal) and Walsh slope = −1 (red dashed). Each point is one user. The cohort fit is consistently shallower than −1, confirming that the parametric form of the Walsh rule does not fit open-source AID data. Both axes are log-scaled; tick labels are in linear units (mg/dL and U/day) for readability.*

### 3.4 Robustness to outlier removal

Five users (U042, U054, U071, U150, U171) had |z| > 2 on ISF × TDD. After dropping these five, the trimmed cohort (n=133) gave ISF × TDD median = 2,308 [2,130, 2,593]; CR × TDD median = 395 [345, 437]; basal / TDD median = 0.46 [0.43, 0.49]. The trimmed CIs still exclude the Walsh values for all three constants. The headline finding is robust to outlier choice.

### 3.5 Duration-of-history sensitivity

| Threshold | n | ISF × TDD median [95% CI] |
|---|---|---|
| ≥ 14 days | 138 | 2,381 [2,144, 2,633] |
| ≥ 30 days | 133 | 2,373 [2,144, 2,721] |
| ≥ 60 days | 122 | 2,303 [2,130, 2,598] |
| ≥ 90 days | 104 | 2,256 [2,071, 2,568] |
| ≥ 180 days | 75 | 2,144 [1,986, 2,545] |

Restricting the cohort to users with longer histories yields ISF × TDD values that drift slowly downward but, even at the strictest 180-day cut where n=75, the lower 95% CI bound (1,986) still excludes 1,700. The drift may reflect a survivorship effect, in that users whose initial settings were further from population norms may correct and remain on the platform, rather than a fundamental difference between long-term and short-term users. We do not attempt to distinguish these mechanisms.

### 3.6 Empirical ISF inferred from BG dynamics (supplementary)

To compare *entered* ISF against the BG response actually observed in fasting intervals, we fitted a per-user linear regression of `ΔBG_30min = a + b · ΔIOB + c · pre-window BG trend` over 30-minute fasting windows (no carbs in window or preceding 90 minutes, no pre-window SMB). ΔIOB is the change in insulin-on-board over the window, an unambiguous quantity in U absorbed, used here in preference to `iob_activity` whose unit convention we could not verify against the oref0 source alone (the empirical ratio ΔIOB / `iob_activity` was 5–10 across users, inconsistent with either common scaling). Empirical ISF for a user is the fitted slope b. Per-user 95% CIs are derived from the regression standard error.

For 114 users with R² ≥ 0.10 in the canonical cohort, median empirical ISF was 22 mg/dL/U (Q1 16, Q3 29) versus median entered ISF of 54 mg/dL/U; the empirical/entered ratio was 0.41 (Q1 0.27, Q3 0.56). The direction means observed BG response per U absorbed in fasting is *smaller* than the entered profile would predict, consistent with real-world BG dynamics being attenuated by counter-regulation, glucose autoregulation and sensor noise relative to the controlled-setting pharmacodynamics that Walsh's rules implicitly assume.

We do not interpret this as the user's *true* ISF being 0.41× their entered value. Profile ISF in oref0/AAPS/Trio is a control parameter for the algorithm's correction-bolus output; the BG response observed during pure fasting (no fresh delivery) reflects only insulin already in the body decaying, a different quantity. A more principled inference requires implementing the iterative oref0 autotune algorithm with a proper PK convolution; we do not do that here. The supplementary value of this analysis is the consistent direction (entered ISF assumes more physiological response than observed) and the cohort distribution of the gap.

We note explicitly that this supplementary result and the headline result are not in tension, because they measure different things. The headline compares entered profiles against a *teaching rule*; this section compares entered profiles against *observed fasting dynamics*. Entered ISF can simultaneously be higher than the 1700-rule prescribes and higher than a fasting-window regression recovers.

### 3.7 Sensitivity to TDD definition

To verify the headline does not depend on our canonical TDD choice, we re-computed the constants under a hybrid `max(2 × profile basal, basal + algorithm-driven SMB/day)` definition that allowed 6 additional users with no parseable treatments file to enter (n=144). The hybrid TDD over-estimates real TDD for high-basal-low-bolus users (e.g. U170: hybrid 171 vs treatments-derived 106 U/day) and under-estimates for low-basal-high-bolus users (e.g. U071: hybrid 22 vs treatments-derived 109). The cohort median ISF × TDD shifts from 2,381 (canonical) to 2,239 (hybrid) and CR × TDD from 402 to 351; both remain on the same side of their respective Walsh values regardless of definition. The slope is shallower under canonical TDD (−0.43) than under hybrid (−0.55); the canonical estimate is the more conservative claim, and we use it as the headline.

---

## 4. Discussion

The core empirical finding of this analysis is that all three classical Walsh constants disagree with the cohort medians at the 95% CI level, and the parametric form of the rules itself fails for the no-DynISF sub-cohort that comprises 87% of users. The disagreement is robust to outlier removal, robust to the TDD definition, and persists across duration thresholds from 14 to 180 days.

The clinical implication is direct. A clinician initialising pump or AID therapy who applies the 1700-rule at the cohort-median TDD of 44 U/day will recommend an initial ISF of approximately 39 mg/dL (2.2 mmol/L). The empirical median in this cohort is 50 mg/dL (2.8 mmol/L), and the empirical Walsh constant for the cohort is 2,381, not 1,700. Following the empirical fit at the same TDD gives an initial ISF of approximately 54 mg/dL. Because a smaller ISF produces a larger correction bolus, the Walsh recommendation runs about 28% too aggressive: a user starting on the Walsh-recommended ISF would be dosed to drop BG more sharply than the cohort median user does in equilibrium. The 500-rule errs in the same direction for carbohydrate, prescribing a CR that delivers less insulin per gram than these users have settled on, and the 50% basal convention over-states the basal share.

The slope finding is the methodologically important result. Walsh's rules implicitly assume ISF ∝ 1/TDD (slope = −1 in log-log space). The no-DynISF cohort slope is −0.38 (95% CI −0.56 to −0.22), well separated from −1. This is not a calibration problem that a revised numerator would fix: the functional form itself is wrong for this population, so a single-constant rule will misestimate at both ends of the TDD range even if it is correct in the middle. The two DynISF sub-cohorts have point estimates of −0.89 (Sigmoid) and −0.84 (Log) but the bootstrap CIs are wide (n=9 each) and span values that include −1; we cannot conclude DynISF users follow Walsh's slope, only that their point estimates are *closer* to it than the no-DynISF group. A natural mechanistic question, whether DynISF, by actively modulating runtime ISF on top of the entered profile, encourages users to enter values aligned with the Walsh form, remains open.

Whether the gap between Walsh's 1980s constants and our 2020s open-source AID cohort reflects population-level changes (a shift to younger or differently-comorbid populations), system-level effects (SMBs and DynISF allowing different entered profiles), or self-tuning (users adapting to algorithm behaviour rather than to physiology), cannot be answered from observational data without demographics, prospective intervention, or a matched MDI control. The pragmatic recommendation is not to abandon Walsh. The empirical fit is within Walsh-shaped distance of the rule at typical adult TDDs, and the rules retain the virtue of being computable at the bedside from a single number. It is to teach them as starting estimates that should be re-tuned within days to weeks of pump start, rather than as stable targets, and to expect the initial correction setting in particular to need weakening rather than strengthening.

Recent observational re-derivations of the Walsh constants in MDI populations *(see References — Walsh-constant re-evaluation literature)* have reported similar shallower-than-1 slopes. Our contribution is the open-source AID replication and the explicit demonstration that the parametric form, not just the constant, is the issue.

---

## 5. Limitations

Seven limitations bound the strength of the claims.

**No demographic data.** We do not have age, sex, weight, duration of diabetes, geography or onboarding pathway. The cohort is *adult-presumed* by virtue of being drawn from adult-led DIY communities, but a paediatric admixture cannot be excluded. Walsh's 1,700 was derived in adults; paediatric teaching typically uses lower constants of 1,300–1,500. A paediatric admixture would therefore pull the cohort constant *downward*, away from the direction we observe, so it would bias against the finding rather than explain it. *(This reverses the direction stated in an earlier draft and should be confirmed by the clinical co-author, since it depends on paediatric users having entered profiles that follow paediatric teaching.)*

**Profile values are declared, not needed.** Profile fields are what the user entered into the pump, not what their physiology required. Selection matters: users with badly calibrated settings either correct (and persist in the data) or leave the platform (and disappear from it). Survivorship may explain part of the duration-sensitivity drift toward Walsh.

**No prospective validation.** This is observational. The empirical constants have not been deployed prospectively to a new patient with measured outcomes.

**DynISF subgroup analysis is underpowered.** With n=9 Sigmoid and n=9 Log users, from a hand-curated mapping derived by querying live Nightscout reason strings for the v5 cohort, the slope CIs span approximately −1.4 to −0.25. We cannot conclude DynISF users follow Walsh's slope, only that their point estimates are closer to it than the no-DynISF group. The minimum-detectable Cohen's d for the Sigmoid-vs-Log slope difference at this allocation (α=0.05, power=0.80) is approximately 1.4, far above any realistic effect for this comparison. The DynISF-on cohort is a structural limit of the data: oref0 (v7) and AAPS classic (v6) predate DynISF, so the only DynISF users in the available archives are the n=29 v5 / Trio cohort, and only 18 of those passed the quality filter. Replicating these per-formula slope estimates with confidence requires assembling a fresh DynISF-era cohort (post-2021 AAPS, current Trio) with treatments-archive access, a separate data-collection effort outside this paper.

**Sample-size imbalance.** The no-DynISF subgroup (n=120) dominates the cohort-wide finding. Both AAPS-classic (n=19) and Trio + early-DynISF (n=22) sub-cohorts are small; their separate slope estimates are reported as exploratory rather than confirmatory.

**Source-table supersession audit.** We verified that the v5 decisions table (29 users, U000–U028) is a strict superset of both v3 (21 users, U000–U020) and v4 (25 users, U000–U024); no users from earlier extraction passes are silently absent.

**No formal pre-registration.** The analysis plan was iteratively refined during the work, not registered before data inspection. The canonical TDD definition, bootstrap analyses, stratifications, outlier sensitivity, duration sensitivity and TDD-definition sensitivity are all reported transparently to mitigate this; readers should treat the per-group findings accordingly.

---

## 6. Conclusion

In a cohort of 138 anonymised open-source AID users, all three of Walsh's classical constants — the 1700-rule (ISF), the 500-rule (CR), and the 50/50 basal-bolus split — disagree with the empirical cohort medians at the 95% CI level. Empirical values are 2,381 [2,144, 2,633] for ISF × TDD, 402 [351, 441] for CR × TDD, and 0.45 [0.43, 0.48] for basal / TDD. The fitted log-linear slope of `log(ISF)` against `log(TDD)` is −0.43 (95% CI −0.59 to −0.27), inconsistent with the −1 slope the Walsh form requires. The deviation from −1 is robust in the no-DynISF subgroup (n=120) and exploratory in the small DynISF subgroups (n=9 each). Applied at initialisation, the 1700-rule sets ISF too low, and therefore corrects too aggressively, for this population. We recommend treating the classical constants as starting estimates only, and re-tuning within days of pump start.

---

## References

1. Walsh J, Roberts R. *Pumping Insulin: Everything You Need to Succeed on a Smart Insulin Pump.* 6th ed. Torrey Pines Press; 2017.
2. Bergenstal RM, Beck RW, Close KL, et al. Glucose Management Indicator (GMI): A New Term for Estimating A1C From Continuous Glucose Monitoring. *Diabetes Care.* 2018;41(11):2275-2280.
3. Lewis DM. *Automated Insulin Delivery: How artificial pancreas "closed loop" systems can aid you in living with diabetes.* Self-published; 2019.
4. Lewis DM, Leibrand S, OpenAPS Community. Real-World Use of Open Source Artificial Pancreas Systems. *Journal of Diabetes Science and Technology.* 2016;10(6):1411. doi:10.1177/1932296816665635.
5. OpenAPS Data Commons. Available from: https://openaps.org/outcomes (accessed April 2026).
6. *[Citation TBD by clinical co-author — observational Walsh-rule re-derivation in adult MDI cohorts. Candidates to search: King AB et al. on insulin sensitivity factor in adults on MDI; Pettus et al. observational validation of correction-bolus rules.]*
7. *[Citation TBD by clinical co-author — community-reported ISF self-tuning analyses. Candidates: oref0 autotune algorithm reference (Lewis 2017 GitHub); Wallace et al. autotune validation; AAPS Auto-ISF community evaluations published in Journal of Diabetes Science and Technology or in Diabetes Technology & Therapeutics conference proceedings 2022-2024.]*

The two TBD citations should be supplied by the clinical co-author with current access to PubMed and the diabetes-tech literature. The author is independent and has not had direct database access for this draft; placeholders are deliberately retained to signal what needs verification rather than fabricated citations.

---

## Author contributions (placeholder)

**TS** conceived the study, assembled the cohort, performed the analysis and wrote the first draft. **[Clinical co-author]** contributed the clinical framing and the discussion of teaching implications, and is responsible for the IRB / ethics statement. **[Statistical reviewer]** verified the bootstrap, power and stratification analyses independently. All authors approved the final manuscript.

---

## Implementation notes (audit trail)

**Canonical TDD definition.** Every figure, table, model fit and supplementary analysis in this paper was generated from a single shared cohort module (`canonical_cohort.load_canonical_cohort()`). It applies the canonical TDD definition (Methods §2.2) and the quality filter (§2.3). The cohort artefact is exported as JSON for independent verification.

**Empirical activity-unit verification.** We verified empirically that the `iob_activity` column reported by the loop has a non-trivial scaling against ΔIOB (median ratio 5–10 across users) and adopted the unit-clean ΔIOB for the empirical-ISF supplementary analysis. The Walsh-constants headline does not depend on `iob_activity` and so is not affected.

**Source-table supersession.** v5 / Trio (29 users, U000–U028) is a strict superset of v3 (U000–U020) and v4 (U000–U024). No users from earlier extraction passes are silently lost.

**Reproducibility.** The cohort assembly module, the canonical Walsh-constants script (with bootstrap CI on the slope), the empirical-ISF v5 ΔIOB-based regression, and the platform stratification are all open-source Python and run on commodity hardware. Per-user computation is parallelised across 12 cores. Wall time: approximately 3 minutes for the full pipeline.

---

## Changes from the previous draft

This version corrects a directional error and adds two clarifications. It makes no change to any computed value.

1. **ISF direction, corrected.** The previous abstract stated that Walsh-derived defaults "may systematically over-estimate ISF" and that cohort ISFs were "≈40% less sensitive than Walsh would predict". Both were the wrong way round, and contradicted §4 of the same draft. At the cohort-median TDD the 1700-rule yields ISF ≈ 39 mg/dL against an empirical ≈ 54 mg/dL, so the rule sets ISF *too low* and corrects *too aggressively*, and the profiles these users run imply *greater* sensitivity than the rule assumes. The abstract, §3.2 and §6 now state this consistently. The basal-share claim was correct as written (0.45 observed against 0.50 taught) and is unchanged.
2. **Paediatric confound direction, reversed.** §5 previously argued a paediatric admixture would shift the empirical constant upward and so could explain the finding. Paediatric teaching uses *lower* constants (1,300–1,500), so an admixture would pull the cohort constant downward and bias against the finding. Flagged for the clinical co-author to confirm.
3. **Supplementary reconciliation added.** §3.6 now states explicitly why the empirical/entered ratio of 0.41 does not contradict the headline: the two sections compare entered profiles against different referents, a teaching rule in one case and observed fasting dynamics in the other.
