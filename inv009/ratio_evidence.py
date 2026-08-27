"""Every independent estimate of measured sensitivity against entered sensitivity.

The estimates in this package disagree with each other on the level. What none of
them disagree about is the direction: in every cohort, by every construction, the
sensitivity a correction achieves is below the sensitivity that was entered. This
assembles the ratios so the spread is visible rather than a single one being
picked, and so any claim made about it is made against the whole set.

Sources differ in three ways that matter. Whether insulin on board was recorded by
the loop itself or reconstructed from delivery records; whether the correcting dose
was counted as a bolus alone or as anything above the programmed basal; and whether
the denominator was the dose given or the insulin action across the window.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd

from . import config

WALSH = 1700.0


def rows() -> list[dict]:
    out = []

    # INV-008, oref archives. Insulin on board is the loop's own recorded value,
    # not a reconstruction, which is the strongest provenance in the series.
    p = config.RESULTS / "best_isf_fit_results.json"
    if p.exists():
        d = json.loads(p.read_text())
        ke, kn = d["sqrt_rule_empirical_K"], d["sqrt_rule_entered_K"]
        out.append(dict(source="INV-008 oref (loop-recorded IOB)", cohort="138 users",
                        n=int(d.get("n_empirical", 0)), measured=ke / math.sqrt(40) * 40,
                        entered=kn / math.sqrt(40) * 40, ratio=ke / kn,
                        note="K/sqrt(TDD) fits, expressed at a daily dose of 40 U"))

    # JAEB, corrections isolated as boluses, fall per unit given.
    p = config.RESULTS / "inv009_correction_landing.json"
    if p.exists():
        for r in json.loads(p.read_text()):
            if not np.isfinite(r.get("entered", np.nan)):
                continue
            out.append(dict(source="JAEB landing (bolus, per unit given)",
                            cohort=r["cohort"], n=r["n_people"],
                            measured=r["achieved"], entered=r["entered"],
                            ratio=r["achieved"] / r["entered"],
                            note="sensitivity in mg/dL/U, not a constant"))

    # JAEB, corrections by any route above the programmed basal.
    p = config.RESULTS / "inv009_correction_routes.json"
    if p.exists():
        for r in json.loads(p.read_text()):
            e = r.get("entered")
            if e is None or not np.isfinite(e):
                continue
            for key, lab in (("isf_settled", "any route, action"),
                             ("isf_clean", "any route, action, clean slate")):
                v = r.get(key)
                if v is None or not np.isfinite(v):
                    continue
                out.append(dict(source=f"JAEB {lab}", cohort=r["cohort"],
                                n=r["n_people"], measured=v, entered=e,
                                ratio=v / e, note="sensitivity in mg/dL/U"))
    return out


def main() -> int:
    config.ensure_dirs()
    r = pd.DataFrame(rows())
    if r.empty:
        print("no sources found")
        return 1
    r.to_parquet(config.RESULTS / "inv009_ratio_evidence.parquet", index=False)

    print("Measured sensitivity against entered sensitivity, every construction\n")
    print(f"{'source':<42s}{'cohort':<12s}{'n':>5s}{'measured':>10s}"
          f"{'entered':>9s}{'ratio':>7s}")
    for src, d in r.groupby("source", sort=False):
        for _, x in d.iterrows():
            print(f"{src:<42s}{x.cohort:<12s}{int(x.n):5d}{x.measured:10.1f}"
                  f"{x.entered:9.1f}{x.ratio:7.2f}")
        src = ""
    print()
    q = r.ratio
    print(f"  {len(q)} estimates across {r.cohort.nunique()} cohorts and "
          f"{r.source.nunique()} constructions")
    print(f"  ratio: min {q.min():.2f}  median {q.median():.2f}  max {q.max():.2f}")
    print(f"  estimates at or above 1.0: {int((q >= 1.0).sum())} of {len(q)}")
    print()
    print("  The direction is unanimous and the magnitude is not. Any claim from")
    print("  this work has to carry the range rather than the median.")

    (config.RESULTS / "inv009_ratio_evidence.json").write_text(json.dumps(
        dict(n_estimates=int(len(q)), ratio_min=float(q.min()),
             ratio_median=float(q.median()), ratio_max=float(q.max()),
             n_at_or_above_one=int((q >= 1.0).sum()),
             estimates=r.to_dict("records")), indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
