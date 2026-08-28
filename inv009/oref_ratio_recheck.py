"""Is the oref archive's 0.41 in the data, or in how it was aggregated?

The oref work reported a constant of 145 anchored to measured sensitivity against
355 anchored to tuned profiles, a ratio of 0.41. The Jaeb cohorts under the same
construction return 0.72 to 1.12. One difference has never been held fixed: the
oref figure comes from fitting a constant over the square root of daily dose across
users, where the Jaeb figures are medians of per-person ratios. A fit across users
and a median of within-person ratios are not the same statistic, and they can
differ by a lot when the underlying distribution is skewed.

The per-person overnight sensitivities from that same run are on disk. Joining them
to each person's entered sensitivity gives the ratio computed the Jaeb way on the
oref records, which isolates the aggregation from everything else.
"""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd

from . import config
from .walsh_constants import load as load_walsh


def main() -> int:
    config.ensure_dirs()
    p = config.RESULTS / "overnight_sensitivity.json"
    if not p.exists():
        print("overnight_sensitivity.json not found")
        return 1
    d = json.loads(p.read_text())
    per = pd.DataFrame(d["per_person"])
    print(f"oref overnight run: {d['n_patients']} people, "
          f"{d['total_windows']:,} windows")
    print(f"method: {d['method']}")
    print(f"population median measured sensitivity: "
          f"{d['population_median_sensitivity']:.1f} mg/dL/U\n")

    w = load_walsh()
    w = w[w.source == "OpenAPS Commons"][["user_id", "isf", "tdd", "isf_x_tdd"]]
    per["user_id"] = per.user.astype(str)
    w["user_id"] = w.user_id.astype(str)
    m = per.merge(w, on="user_id", how="inner")
    if m.empty:
        # user keys may differ in prefix; try a suffix match
        w["k"] = w.user_id.str.extract(r"(\d+)$")[0]
        per["k"] = per.user_id.str.extract(r"(\d+)$")[0]
        m = per.merge(w, on="k", how="inner", suffixes=("", "_w"))
    print(f"matched to entered settings: {len(m)} people\n")
    if m.empty:
        print("no key overlap between the two sources; cannot join")
        return 1

    m = m[(m.isf > 0) & (m.median_sens > 0) & (m.tdd > 0)].copy()
    m["ratio"] = m.median_sens / m.isf
    m["k_measured"] = m.median_sens * m.tdd

    q = m.ratio
    print("Ratio of measured to entered sensitivity, computed per person\n")
    print(f"  n {len(m)}")
    print(f"  median {q.median():.2f}   quartiles {q.quantile(.25):.2f} "
          f"to {q.quantile(.75):.2f}   range {q.min():.2f} to {q.max():.2f}")
    print(f"  people at or above 1.0: {int((q >= 1).sum())} of {len(q)} "
          f"({100*(q >= 1).mean():.0f}%)\n")

    print("The same records, aggregated the two ways\n")
    print(f"{'statistic':<44s}{'value':>8s}")
    print(f"{'median of per-person ratios':<44s}{q.median():8.2f}")
    kfit_m = float(np.median(m.median_sens * np.sqrt(m.tdd)))
    kfit_e = float(np.median(m.isf * np.sqrt(m.tdd)))
    print(f"{'K over sqrt(TDD), measured':<44s}{kfit_m:8.1f}")
    print(f"{'K over sqrt(TDD), entered':<44s}{kfit_e:8.1f}")
    print(f"{'their ratio':<44s}{kfit_m/kfit_e:8.2f}")
    print(f"{'reported in the earlier work':<44s}{145/355.1:8.2f}")
    print()
    print(f"{'measured constant, median ISF x TDD':<44s}{float(m.k_measured.median()):8.0f}")
    print(f"{'entered constant, median ISF x TDD':<44s}{float(m.isf_x_tdd.median()):8.0f}")

    (config.RESULTS / "inv009_oref_ratio_recheck.json").write_text(json.dumps(
        dict(n=int(len(m)), ratio_median=float(q.median()),
             ratio_q1=float(q.quantile(.25)), ratio_q3=float(q.quantile(.75)),
             frac_at_or_above_one=float((q >= 1).mean()),
             k_sqrt_measured=kfit_m, k_sqrt_entered=kfit_e,
             k_sqrt_ratio=kfit_m / kfit_e,
             k_measured=float(m.k_measured.median()),
             k_entered=float(m.isf_x_tdd.median())), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
