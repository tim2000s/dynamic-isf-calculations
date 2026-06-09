# Executive summary — Should the insulin sensitivity factor change with glucose?

**Tim Street, with Claude (Anthropic) · 2026-06-09**
*Dynamic ISF in open-source automated insulin delivery (oref0 / Trio)*

---

**The question.** Dynamic-ISF algorithms (and the Diabeloop population model) lower the insulin
sensitivity factor — the expected glucose fall per unit of insulin — as glucose rises, on the premise
that high glucose makes insulin less effective. In practice, dynamic equations rarely beat a well-set
static ISF. We tested why, on real closed-loop data.

**Bottom line.** The ISF an automated system should actually use — the factor relating its insulin to
the *realised net glucose drop* — **does not fall with glucose.** It is suppressed near target (the
body defends against hypos) and flat-to-mildly-rising above. A per-user-adapted static ISF is the
right design; a glucose-lowering curve predicts drops that do not happen and over-doses highs.

## Key findings

- **Static beats dynamic.** On identical fasting windows, a well-set static ISF predicts the glucose
  drop as well as the loop itself (median error 20 vs 19 mg/dL) and far better than the dynamic v1/v2
  equations (25 / 50 mg/dL). The dynamic forms over-steepen the dose dependence.
- **The net effective ISF does not decline with glucose** — confirmed two independent ways, including a
  method that uses no assumptions about how insulin acts over time. The Diabeloop/dynamic curves fall
  steeply and diverge from the data; transplanted onto this cohort they degrade prediction.
- **What helps is personalising the *level*, not adding a curve.** Fitting each person's own
  sensitivity scale recovers the largest gain (~7 mg/dL); the best-fit glucose steepness is zero, and a
  per-user glucose curve fails out-of-sample.
- **The apparent glucose effect is mostly correction *size*, not glucose.** The loop over-predicts
  large corrections by ~2×; this masquerades as a high-glucose effect and belongs in the insulin/basal
  model, not the ISF.

## Recommendation

> **ISF = a per-user-adapted sensitivity level (K/√TDD starting point, refined online from the
> person's own outcomes) + a near-target easing clamp for hypo safety. No glucose-dependent
> correction term.** Validate prospectively (shadow mode) before deployment; flag kidney-function /
> SGLT2 users, where the picture may differ.

## What we did *not* establish (honesty note)

Molecular physiology genuinely shows high glucose causes insulin resistance. A plausible reconciliation
— that resistance is real but cancelled by glucose-driven, insulin-independent clearance — could *not*
be confirmed: the clean data needed to separate the two does not exist in fasting records, and an
earlier claim to have resolved it was withdrawn after audit. The actionable conclusion is unaffected,
because it rests on the *net* response, which is what the controller doses against.

## Scale

**119** eligible individuals (29 Trio + 110 oref0; 73–100 per analysis), **~9.6 million** raw loop
decisions, **~62,700 overnight** + **~64,300 daytime** carb-screened correction windows, **>560,000**
candidate-ISF evaluations.

*Full paper, methods, audit and code: `github.com/tim2000s/dynamic-isf-calculations`.*
