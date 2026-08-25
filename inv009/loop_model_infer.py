"""Which insulin model was each Loop user running, and does our arithmetic match theirs?

Loop recorded its own insulin on board at every bolus, and that number is a gift.
It is the one place in this entire dataset where an app's internal state was
written down, so it can be recomputed from the pump record and checked against
what the app actually believed.

It does two jobs at once. It identifies each person's insulin model, because
Loop offered several and the archives do not say which was chosen, and the one
whose recomputed insulin on board tracks the recorded value is the one they were
running. And it validates the dose reconstruction itself: if basal, temp basals,
extended boluses and the scheduled profile are being assembled correctly, the
numbers agree; if any of that is wrong, they do not.

The important subtlety is that Loop's insulin on board is NET of scheduled basal.
A temp basal contributes only the difference between what it delivered and what
the profile would have delivered anyway, and a person sitting on their profile
has no basal contribution at all. Reconstructing it from total delivery gives a
number that climbs all night and matches nothing, which is why the scheduled
profile had to be recovered from the raw archives in the first place.

    python3 -m inv009.loop_model_infer
"""
from __future__ import annotations

import json
import multiprocessing as mp
import sys

import numpy as np
import pandas as pd

from . import config, db, grid as gridmod, insulin_models as M

MIN_RECORDS = 40
ACCEPT_REL_MAD = 0.15
ACCEPT_MAD_U = 0.30


def assess(job) -> dict | None:
    subject_id = job
    streams = db.streams(subject_id)
    wiz = streams["wizard"]
    if wiz.empty or wiz.iob_u.notna().sum() < MIN_RECORDS:
        return None
    g = gridmod.build_grid(streams)
    if g is None or g.sched_u.isna().all():
        return None

    # Net of the profile, which is what Loop counts.
    net = (g.total_u - g.sched_u.fillna(0.0)).to_numpy(float)
    bolus = g.bolus_u.to_numpy(float)
    ts = g.ts.values

    w = wiz.dropna(subset=["iob_u"])
    idx = np.searchsorted(ts, w.ts_local.values) - 1
    ok = (idx >= 0) & (idx < len(g))
    idx, recorded = idx[ok], w.iob_u.to_numpy(float)[ok]
    if len(idx) < MIN_RECORDS:
        return None

    out = dict(subject_id=subject_id, n=int(len(idx)),
               recorded_median=float(np.median(recorded)))
    best, best_mad = None, np.inf
    for m in config.LOOP_MODELS:
        kern = M.kernel(m)
        iob = np.convolve(net, kern)[:len(net)]
        # The recorded value precedes the bolus it was recommending.
        mine = iob[idx] - bolus[idx]
        mad = float(np.median(np.abs(mine - recorded)))
        out[f"mad_{m}"] = mad
        out[f"bias_{m}"] = float(np.median(mine - recorded))
        out[f"corr_{m}"] = float(np.corrcoef(mine, recorded)[0, 1]) if len(idx) > 3 else np.nan
        if mad < best_mad:
            best, best_mad = m, mad
    # Scale on the spread rather than the level: plenty of these subjects sit at
    # essentially no insulin on board most of the night, and dividing by that
    # turns a third of a unit into a relative error of three million.
    scale = max(float(np.percentile(np.abs(recorded), 75)), 0.2)
    out["model"] = best
    out["mad"] = best_mad
    out["rel_mad"] = best_mad / scale
    out["corr_best"] = out[f"corr_{best}"]   # not "corr": that shadows DataFrame.corr
    out["bias_best"] = out[f"bias_{best}"]
    out["accepted"] = bool(out["rel_mad"] <= ACCEPT_REL_MAD or best_mad <= ACCEPT_MAD_U)
    return out


def main() -> int:
    config.ensure_dirs()
    subs = db.subjects("Loop").subject_id.tolist()
    with mp.Pool(config.WORKERS, maxtasksperchild=8) as pool:
        rows = [r for r in pool.imap_unordered(assess, subs, chunksize=4) if r]
    R = pd.DataFrame(rows)
    if R.empty:
        print("no Loop subject had enough recorded insulin on board")
        return 0
    R.to_parquet(config.RESULTS / "inv009_loop_model_choice.parquet", index=False)

    res = dict(n_subjects=int(len(R)),
               accepted=int(R.accepted.sum()),
               accept_rate=float(R.accepted.mean()),
               median_rel_mad=float(R.rel_mad.median()),
               median_mad_u=float(R.mad.median()),
               median_corr=float(R["corr_best"].median()),
               median_bias_u=float(R["bias_best"].median()),
               bias_by_model={m: float(R[f"bias_{m}"].median()) for m in config.LOOP_MODELS},
               corr_by_model={m: float(R[f"corr_{m}"].median()) for m in config.LOOP_MODELS},
               model_counts={k: int(v) for k, v in R.model.value_counts().items()},
               model_counts_accepted={k: int(v) for k, v
                                      in R[R.accepted].model.value_counts().items()})
    (config.RESULTS / "inv009_loop_model_choice.json").write_text(json.dumps(res, indent=2))

    print(f"Loop subjects with recorded insulin on board: {len(R)}")
    print(f"  agreement with what the app recorded: median correlation {res['median_corr']:.3f}, "
          f"median absolute difference {res['median_mad_u']:.3f} U "
          f"({100 * res['median_rel_mad']:.1f}% of a typical value)")
    print(f"  accepted at 15% or 0.3 U: {res['accepted']} of {len(R)} "
          f"({100 * res['accept_rate']:.0f}%)")
    print("\n  model chosen (count / accepted / median bias U / median correlation)")
    for m in config.LOOP_MODELS:
        print(f"    {m:11s} {res['model_counts'].get(m, 0):4d}  "
              f"{res['model_counts_accepted'].get(m, 0):4d}  "
              f"{res['bias_by_model'][m]:+7.3f}  {res['corr_by_model'][m]:.3f}")
    print("\n  The correlation is the finding: the shape of insulin on board over a")
    print("  night is reproduced from the pump record alone. WHICH model wins is a")
    print("  weaker claim, because a small positive level bias in the reconstruction")
    print("  trades off against how fast a model decays, and the fastest model")
    print("  absorbs it. The analyses do not depend on the choice.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
