# Tier-2 anchor — quantifying the glucose double-count

98 people. The measured ISF is a ΔIOB²-weighted average of the true ISF over each person's fasting windows; with no glucose-level term in the regression, it carries the average g(BG) of those windows. Tier-2 then re-applies g(BG), so its ISF is off by the mean-g factor at every glucose.

## Result

- Fasting windows sit at a median BG of **112 mg/dL** (target 99); a median **70%** of each person's windows are above target.
- Mean g(BG) over the windows (the bias factor): median **0.782** [IQR 0.677–0.849].
- So the corrected at-target anchor = measured_ISF / 0.782 → a median **1.278×** the current anchor (IQR 1.18–1.477×).

![Tier-2 anchor de-bias](charts/inv008/fig_tier2_anchor_debias.png)

## Reading

The double-count has a median factor of 0.782, and the factor was below 1 for 96% of people. Where windows sat above target (g<1) the current anchor is too low and Tier-2 over-doses; where below (g>1) it under-doses. The clean fix is to divide the measured slope by this mean-g factor before using it as the at-target anchor, or to fit the slope on near-target windows only.

## Per-person
