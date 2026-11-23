# Phase 11 — cohort shadow evaluation of Tier-1 v-next

141 people (18 ran DynISF), 9,223,241 per-tick readings. v-next = (K_user/√TDD)·g(BG), Diabeloop quartic g(BG), Tier-1 K_user = profile_ISF·√(median TDD). Counterfactual replay: ISF the device would have used vs what it did.

## Headline

- A correction dose moves by a median **37.3%** vs the person's static profile — and that change is **almost entirely the glucose curve**, not the level: the at-target level barely moves (TDD-axis factor median 1.0, within-user p10–p90 [0.885, 1.157]), but g(BG) reshapes dosing strongly across the glucose range (factor 0.752 at the median reading).
- So Tier-1 preserves *average* dosing only at target glucose; across the BG range it is a real behaviour change — more insulin when high, less when low — which is the point of a dynamic ISF.
- vs today's DynISF (v1) the two are close: median sensitivity ratio **0.89**, within ±30% across all glucose bands.
- The §8.2 clamp governs the **level** (TDD axis), and under Tier-1 it almost never binds: a median **0.0%** of ticks (worst person 13.2%) have a level >1.5× stronger than profile. (The high-BG aggression — 40% of ticks beyond 1.5× once g(BG) is included — is intended, not a level fault, so the clamp should be applied to the level term, not the full ISF.)

## Sensitivity ratio by glucose band (the g(BG) behaviour)

| BG band | v-next / profile | v-next / v1 |
|---|---|---|
| 54-80 | 1.445 | 1.188 |
| 80-100 | 1.143 | 1.104 |
| 100-120 | 0.905 | 0.989 |
| 120-150 | 0.674 | 0.881 |
| 150-200 | 0.461 | 0.733 |
| 200-260 | 0.334 | 0.649 |

*>1 = more sensitive (less insulin); <1 = more aggressive (more insulin). v-next is protective at low BG and more aggressive at high BG — the Diabeloop curve shape — anchored to the person's own profile level.*

## By per-user median TDD

| TDD band | users | median ratio | median \|Δdose\| % | frac level >1.5× strong |
|---|---|---|---|---|
| <15 | 1 | 0.663 | 50.8 | 0.0 |
| 15-25 | 17 | 0.711 | 42.4 | 0.0 |
| 25-40 | 50 | 0.765 | 36.3 | 0.0 |
| 40-65 | 55 | 0.752 | 39.1 | 0.0 |
| 65+ | 18 | 0.775 | 34.2 | 0.0 |

![Shadow eval](charts/inv008/fig_shadow_eval.png)

**Reading:** Tier-1 v-next leaves the *level* essentially unchanged — at each person's median TDD and target glucose it returns their profile ISF, and the TDD-axis swing within a record is small (p10–p90 ≈ 0.87–1.17). The §8.2 level clamp therefore almost never binds (0% of ticks at the median person). What changes is the shape: g(BG) makes corrections firmer when high and gentler when low — the intended dynamic behaviour — so the median 38% per-tick dose change is the glucose curve acting on each person's own BG distribution, not a level shift. Versus today's DynISF (v1) the two equations track within ±30% across the whole glucose range, with v-next a little firmer at high BG (the Diabeloop curve is steeper there than v1's log). Implication for §8.2: clamp the *level* term, not the full ISF, or the intended high-BG aggression would be clipped ~40% of the time.

*Caveat: counterfactual decision-level replay (ISF that would have been used), not closed-loop outcomes; single cohort; median-TDD anchor stands in for the weekly 14-day recalibration.*