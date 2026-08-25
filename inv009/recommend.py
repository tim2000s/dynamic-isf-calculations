"""What sensitivity law the archives actually support, and how to adjust it in real time.

Two questions, both scored the same way the equations were scored: predict the
overnight fall, fit only what a deployed system could fit, and measure on periods
the fit never saw.

Part one asks what to derive a sensitivity from. Three candidate inputs, all
available to a pump: a seven day average of total daily dose, a longer average
over the whole record, and the blended figure both dynamic ISF equations use.
Each is tried across a range of exponents, so the answer is a shape rather than
an opinion.

Part two asks whether a static number should then be adjusted, and by what. A
sensitivity that never moves ignores illness, exercise and everything else that
shifts week to week. The comparison here is against a mechanism of the same
family as autosens: watch how far recent nights departed from what the
sensitivity predicted, and carry that forward. It is not oref's autosens, which
works on a rolling window of five minute deviations and is bounded by
preferences. It tests the idea autosens embodies, which is adaptation from recent
error, at the cadence this data supports.

    python3 -m inv009.recommend
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from . import config, data, db

TRAIN_FRAC = 0.7
MIN_TEST = 20
MIN_WINDOWS = 80
# ISF = K / TDD**b, so b is POSITIVE when sensitivity falls as TDD rises. The
# 1800 rule is b = 1 and the square root form is b = 0.5. Reported elsewhere in
# this series as a slope of -0.83, which is the same thing with the sign of the
# log-log gradient rather than of the divisor.
EXPONENTS = np.round(np.arange(0.0, 1.61, 0.05), 2)
TDD_SOURCES = {"tdd_7d": "seven day average",
               "tdd_u": "whole record average",
               "tdd_blend": "the blend both equations use"}


def _subject_frames():
    for study in config.COHORTS:
        w = data.load(study)
        if w.empty:
            continue
        w = w.dropna(subset=["a_pre", "drop", "tdd_u", "tdd_7d", "tdd_blend"])
        w = w[(w.tdd_u > 0) & (w.tdd_7d > 0) & (w.tdd_blend > 0)]
        for sid, d in w.groupby("subject_id", sort=False):
            if len(d) < MIN_WINDOWS:
                continue
            d = d.sort_values("t0")
            n_tr = int(len(d) * TRAIN_FRAC)
            if len(d) - n_tr < MIN_TEST:
                continue
            yield sid, study, d, n_tr


def fit_constant(frames, source: str, b: float) -> float:
    """The single constant K in ISF = K / TDD^b, fitted across everyone's training half.

    One number for the whole population, which is what a rule of thumb is. Each
    person still gets their own intercept, because overnight most insulin is
    basal and offsets hepatic output, which no sensitivity factor should be asked
    to carry.
    """
    num = den = 0.0
    for _, _, d, n_tr in frames:
        tr = d.iloc[:n_tr]
        a = tr.a_pre.to_numpy(float)
        y = tr["drop"].to_numpy(float)
        x = a / np.power(tr[source].to_numpy(float), b)
        x = x - x.mean()                      # the intercept is per person
        y = y - y.mean()
        num += float(np.sum(x * y))
        den += float(np.sum(x * x))
    return num / den if den > 0 else np.nan


def score(frames, source: str, b: float, k: float) -> tuple[float, float]:
    errs = []
    for _, _, d, n_tr in frames:
        tr, te = d.iloc[:n_tr], d.iloc[n_tr:]
        isf_tr = k / np.power(tr[source].to_numpy(float), b)
        isf_te = k / np.power(te[source].to_numpy(float), b)
        pred_tr = isf_tr * tr.a_pre.to_numpy(float)
        pred_te = isf_te * te.a_pre.to_numpy(float)
        off = float(np.mean(tr["drop"].to_numpy(float) - pred_tr))
        e = te["drop"].to_numpy(float) - (pred_te + off)
        errs.append(np.median(np.abs(e)))
    return float(np.median(errs)), float(len(errs))


def adaptive_score(frames, source: str, b: float, k: float,
                   half_life: float, bound: tuple[float, float]) -> float:
    """Carry recent error forward, in the spirit of autosens.

    After each night, compare what happened with what the sensitivity predicted,
    and update a running ratio. The ratio multiplies the sensitivity on the next
    night. `half_life` is in nights and `bound` clamps the ratio, both matching
    the shape of the real mechanism rather than its exact parameters.
    """
    decay = 0.5 ** (1.0 / half_life)
    errs = []
    for _, _, d, n_tr in frames:
        tr, te = d.iloc[:n_tr], d.iloc[n_tr:]
        isf_tr = k / np.power(tr[source].to_numpy(float), b)
        off = float(np.mean(tr["drop"].to_numpy(float) - isf_tr * tr.a_pre.to_numpy(float)))
        # Warm the ratio on the training half so the test half does not pay for it.
        ratio = 1.0
        for src in (tr, te):
            a = src.a_pre.to_numpy(float)
            y = src["drop"].to_numpy(float)
            isf = k / np.power(src[source].to_numpy(float), b)
            out = []
            for i in range(len(src)):
                pred = ratio * isf[i] * a[i] + off
                out.append(y[i] - pred)
                if a[i] > 0.3:
                    # What ratio would have made tonight right, damped and clamped.
                    want = (y[i] - off) / (isf[i] * a[i])
                    if np.isfinite(want):
                        ratio = decay * ratio + (1 - decay) * float(np.clip(want, *bound))
                        ratio = float(np.clip(ratio, *bound))
            if src is te:
                errs.append(np.median(np.abs(np.asarray(out))))
    return float(np.median(errs))


def main() -> int:
    config.ensure_dirs()
    frames = list(_subject_frames())
    print(f"{len(frames)} people with at least {MIN_WINDOWS} nights\n")
    res: dict = {"n_people": len(frames)}

    print("PART ONE  what to derive the sensitivity from")
    print(f"  {'input':>28s} {'best exponent':>14s} {'constant':>10s} {'error':>9s}")
    best = None
    res["part_one"] = []
    for source, label in TDD_SOURCES.items():
        curve = []
        for b in EXPONENTS:
            k = fit_constant(frames, source, b)
            if not np.isfinite(k):
                continue
            mae, _ = score(frames, source, b, k)
            curve.append((float(b), float(k), mae))
        b_best, k_best, mae_best = min(curve, key=lambda t: t[2])
        res["part_one"].append(dict(source=source, label=label, exponent=b_best,
                                    constant=k_best, mae=mae_best, curve=curve))
        print(f"  {label:>28s} {b_best:>14.2f} {k_best:>10.0f} {mae_best:>9.2f}")
        if best is None or mae_best < best[3]:
            best = (source, b_best, k_best, mae_best)

    # The blend usually wins by a hair, and the hair is not physiology. Within a
    # person it moves with the very insulin the sensitivity multiplies, so
    # dividing by it shrinks the regressor toward its own mean. A stable average
    # cannot do that. Report the correlation and prefer the stable source unless
    # the blend wins by more than the margin below.
    from scipy import stats as sps
    corr = {c: [] for c in TDD_SOURCES}
    for _, _, d, _ in frames:
        for c in TDD_SOURCES:
            v = d[c].to_numpy(float)
            if np.std(v) > 1e-9:
                corr[c].append(sps.spearmanr(d.a_pre, v).statistic)
    res["regressor_leak"] = {c: (float(np.median(v)) if v else None) for c, v in corr.items()}
    print("\n  how far each input moves with the insulin it is dividing:")
    for c, v in res["regressor_leak"].items():
        note = "  cannot leak, constant per person" if v is None else ""
        print(f"    {TDD_SOURCES[c]:>28s} {('%+.3f' % v) if v is not None else '   n/a':>8s}{note}")
    stable = min((r for r in res["part_one"] if r["source"] == "tdd_7d"),
                 key=lambda r: r["mae"])
    if best[0] == "tdd_blend" and (stable["mae"] - best[3]) < 0.5:
        print(f"    the blend leads by only {stable['mae'] - best[3]:.2f} mg/dL and leaks, "
              f"so the seven day average is preferred")
        best = ("tdd_7d", stable["exponent"], stable["constant"], stable["mae"])
    src, b, k, mae = best
    print(f"\n  best: ISF = {k:.0f} / TDD^{abs(b):.2f} using the {TDD_SOURCES[src]}")
    print(f"  for somebody on 40 U/day that is {k / 40 ** b:.0f} mg/dL per unit, "
          f"and on 15 U/day {k / 15 ** b:.0f}")
    # How much worse are the round-number alternatives?
    for label, bb in (("1/TDD, the 1800 rule", 1.0), ("square root", 0.5),
                      ("1/TDD squared, v2's law", 2.0), ("no TDD term at all", 0.0)):
        kk = fit_constant(frames, src, bb)
        mm, _ = score(frames, src, bb, kk)
        print(f"    against {label:24s} K={kk:8.1f}, b={bb:.2f}  error {mm:.2f} "
              f"({mm - mae:+.2f})")
    res["best"] = dict(source=src, exponent=b, constant=k, mae=mae)

    print("\nPART TWO  whether to adjust it, and how fast")
    print(f"  {'mechanism':>34s} {'error':>9s} {'change':>8s}")
    print(f"  {'no adjustment at all':>34s} {mae:>9.2f} {'':>8s}")
    res["part_two"] = [dict(mechanism="none", mae=mae, delta=0.0)]
    for hl in (1.0, 3.0, 7.0, 14.0, 30.0):
        for bound in ((0.7, 1.2), (0.5, 1.5)):
            m = adaptive_score(frames, src, b, k, hl, bound)
            tag = f"adapt, half-life {hl:g} nights, {bound[0]}-{bound[1]}"
            res["part_two"].append(dict(mechanism=tag, mae=m, delta=m - mae))
            print(f"  {tag:>34s} {m:>9.2f} {m - mae:>+8.2f}")
    (config.RESULTS / "inv009_recommend.json").write_text(json.dumps(res, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
