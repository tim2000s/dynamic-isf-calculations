"""The sensitivity constant fitted from observed glucose response, by cohort.

The matched-correction estimator elsewhere in this package answers whether a
single extra unit can be shown to have done anything, and under an algorithm it
cannot, because the comparison nights were corrected through basal. That is a
statement about one estimator and not about the data. A constant can still be
fitted, by asking which law of the form ISF = K / TDD^b best predicts the
overnight fall on periods the fit never saw.

The fit gives every person their own intercept, because overnight most insulin is
basal offsetting hepatic glucose output, and no sensitivity factor should be asked
to carry that. K is therefore the marginal fall per acting unit relative to each
person's own baseline trajectory. This is why it is smaller than the constant in a
settings file, which encodes the whole expected fall including the part that
happens without insulin, and the two should not be compared without saying so.

Scoring is a per-person 70/30 split in time, median absolute error on the held-out
tail, pooled as the median across people.
"""
from __future__ import annotations

import json

import numpy as np

from . import config
from .recommend import _subject_frames, fit_constant, score

SOURCE = "tdd_7d"
GRID = [round(x, 2) for x in np.arange(0.0, 1.65, 0.05)]
MIN_PEOPLE = 25
N_BOOT = 400


def curve(frames, source: str = SOURCE) -> list[tuple[float, float, float]]:
    out = []
    for b in GRID:
        k = fit_constant(frames, source, b)
        m, _ = score(frames, source, b, k)
        out.append((float(b), float(k), float(m)))
    return out


def at(c, b: float) -> tuple[float, float]:
    row = min(c, key=lambda t: abs(t[0] - b))
    return row[1], row[2]


def main() -> int:
    config.ensure_dirs()
    frames = list(_subject_frames())
    by: dict[str, list] = {}
    for t in frames:
        by.setdefault(t[1], []).append(t)

    pooled = curve(frames)
    b_star, k_star, mae_star = min(pooled, key=lambda t: t[2])
    k1, mae1 = at(pooled, 1.0)

    rng = np.random.default_rng(20260826)
    boot = [fit_constant([frames[i] for i in rng.integers(0, len(frames), len(frames))],
                         SOURCE, 1.0) for _ in range(N_BOOT)]
    ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))

    rows = []
    for st, fr in sorted(by.items()):
        if len(fr) < MIN_PEOPLE:
            continue
        c = curve(fr)
        b_c, k_c, m_c = min(c, key=lambda t: t[2])
        k_c1, m_c1 = at(c, 1.0)
        rows.append(dict(cohort=st, n=len(fr), best_b=b_c, k_at_best=k_c,
                         mae_at_best=m_c, k_at_1=k_c1, mae_at_1=m_c1))

    out = dict(n_people=len(frames), source=SOURCE,
               pooled=dict(best_b=b_star, k_at_best=k_star, mae_at_best=mae_star,
                           k_at_1=k1, mae_at_1=mae1, k_at_1_ci=list(ci)),
               by_cohort=rows)
    (config.RESULTS / "inv009_constant_by_cohort.json").write_text(json.dumps(out, indent=1))

    print(f"Sensitivity constant fitted to observed overnight fall, {len(frames)} people\n")
    print(f"  pooled, at Walsh's exponent of 1.0 : K = {k1:.0f} "
          f"[{ci[0]:.0f}, {ci[1]:.0f}]   MAE {mae1:.2f} mg/dL")
    print(f"  pooled, at the best-fitting {b_star:.2f}   : K = {k_star:.0f}"
          f"                MAE {mae_star:.2f} mg/dL")
    print(f"  the two differ by {mae1 - mae_star:.2f} mg/dL, so Walsh's form costs "
          f"almost nothing here\n")
    print(f"{'cohort':<12s}{'people':>7s}{'best b':>8s}{'K there':>9s}"
          f"{'K at b=1':>10s}{'MAE':>8s}")
    for r in rows:
        print(f"{r['cohort']:<12s}{r['n']:7d}{r['best_b']:8.2f}{r['k_at_best']:9.0f}"
              f"{r['k_at_1']:10.0f}{r['mae_at_1']:8.2f}")
    bs = [r["best_b"] for r in rows]
    ks = [r["k_at_1"] for r in rows]
    print(f"\n  best exponent spans {min(bs):.2f} to {max(bs):.2f}")
    print(f"  constant at b=1 spans {min(ks):.0f} to {max(ks):.0f}, "
          f"a factor of {max(ks) / min(ks):.2f}")
    print("\n  So a constant is obtainable and a single universal one is not "
          "supported by\n  these cohorts. Both statements are needed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
