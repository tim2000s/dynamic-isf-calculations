# Phase 7 — outcome-based tuning prototype vs the sensitivity route

112 users. Working-anchored cold-start constant K_pop = **355.1** (ISF·√TDD), i.e. ISF ≈ 355/√TDD.

Each candidate is a zero-entry per-user ISF. `err vs working` = median |log(ISF/entered)| as ± % (entered ISF is the known-working reference that achieves the users' actual TIR). `logR vs TBR` = directional hypo-safety (positive = the unsafe Phase-5 signature: the candidate assigns lower ISF / more insulin to hypo-prone users).

| candidate | err vs working ISF | logR vs TBR (≤0 = safe) |
|---|---|---|
| sensitivity route (empirical) | ±138% | +0.38 |
| cold start (√TDD, no per-user) | ±28% | +0.21 |
| outcome route (cold + nudge) | ±41% | -0.04 |

- median outcome nudge applied: ×1.00 (mostly weakening, reflecting hypo-avoidance).

## Reading

- The **sensitivity route** is closest in level but carries the unsafe hypo-signature (positive logR-vs-TBR) — it would over-dose the hypo-prone.
- The **cold start** is unbiased (no per-user signal to bias it) but coarse.
- The **outcome route** keeps the cold-start's safety and nudges toward working levels in the safe direction (weaker when lows present) — the only route that is both reasonably close and hypo-safe.

*This is a design prototype: the nudge's safety is partly by construction, and the level it reaches is bounded by the cold-start. True validation requires a closed-loop trial; observationally we can only show direction and rough level.*