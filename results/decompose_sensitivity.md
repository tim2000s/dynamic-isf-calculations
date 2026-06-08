# Can the glucose↔sensitivity relation and the at-100 anchor be measured from logs?

100 people, 81,182 overnight carb-screened windows (BG ≥ target). Model: drop = a + b·IOB + c·(BG−100) + d·IOB·(BG−100); c is meant to absorb glucose reversion so b is the insulin effect at 100 and d its glucose dependence.

## The attempt fails to identify insulin — and that is the answer

- The reversion term alone predicts a **46 mg/dL** drop from BG 160 with **zero insulin** — overnight glucose reverts toward target on its own.
- IOB is **collinear with glucose** (corr IOB vs BG−100 = 0.555): the loop doses by glucose, so the insulin and reversion terms cannot be separated. The fitted ISF at 100 collapses to 1.94 mg/dL per U — not physiological, the signature of an unidentified model.
- The cleanest test: at the same glucose, high-IOB and low-IOB windows drop almost the same (median gap ≈ **7 mg/dL**). The 4-hour drop is set by where glucose started, not by how much insulin was on board.

## Q1 — sensitivity vs glucose: not recoverable

The insulin effect is not separable from reversion in closed-loop data, so the glucose shape g(BG) cannot be measured from these logs. It must come from controlled (clamp / clinical) data.

## Q2 — per-user ISF at 100: not from drop/IOB

At target the loop holds glucose there, so there is essentially no excursion to measure — only a handful of users have enough near-target 4h windows, and those read a loop-suppressed value. The usable per-user anchor at 100 is the **tuned profile ISF** (cohort median ≈ 50 mg/dL per U), or a **K/√TDD cold-start** where no profile exists.

## Q3 — value at 100 vs the changes: moot

With the glucose changes unobservable and the anchor coming from profiles, the link cannot be derived from logs; the multiplicative form ISF(BG) = ISF₁₀₀ · g(BG) is a modelling choice for clinical data to validate.

![Diagnostic](charts/inv008/fig_decompose_sensitivity.png)

**Conclusion.** Closed-loop observational data cannot identify the insulin sensitivity here. The homeostasis the loop enforces creates two killers: (1) glucose reverts toward target regardless of insulin (~46 mg/dL drop from BG 160 with zero insulin), and that reversion is collinear with IOB (corr 0.56, the loop doses by glucose), so the insulin term collapses when reversion is controlled; (2) at target the loop holds BG there, so almost no excursion exists to measure (few near-target 4h windows, and those that exist read a loop-suppressed value). The best per-user ISF at 100 therefore comes from tuned profiles or the cross-sectional K/√TDD cold-start, not from drop/IOB.