"""Effective sensitivity from the overnight carb-free windows, by study.

This answers the question the entered-settings analysis cannot: not what people
typed into a pump, but what a unit of insulin did to glucose overnight when no
carbohydrate was involved. It is the same window machinery the dynamic ISF work
used, reported per study and converted to a Walsh-style constant.

The headline result is that the level does not survive the estimation, and the
reason is measurable rather than mysterious. Overnight glucose falls on its own,
faster from a higher start, so an estimator that does not remove that fall
measures mean reversion, and one that does remove it returns the marginal effect
of an extra unit. A settings file encodes neither: it encodes the total fall a
person expects. The three quantities are reported side by side.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import config, data
from inv008.err_common import boot_median_ci

STUDIES = ["ReplaceBG", "Loop", "DCLP3", "DCLP5", "PEDAP", "IOBP2"]
LABEL = {"ReplaceBG": "REPLACE-BG (open loop)", "Loop": "Loop (DIY closed loop)",
         "DCLP3": "DCLP3 (Control-IQ, adult)", "DCLP5": "DCLP5 (Control-IQ, 6-13y)",
         "PEDAP": "PEDAP (Control-IQ, 2-5y)", "IOBP2": "IOBP2 (bionic pancreas)"}
BANDS = [(90, 120), (120, 150), (150, 180), (180, 220), (220, 260), (260, 400)]


def spontaneous(w: pd.DataFrame) -> list[dict]:
    """Median four-hour fall on nights with no bolus at all, by starting glucose.

    This is the counterfactual every sensitivity estimate is measured against, and
    it is large enough that whether it is subtracted changes the answer fivefold.
    """
    q = w[(w.quiet_before_h >= config.EVENT_QUIET_H) & (w.bolus_in_u <= 1e-9)]
    rows = []
    for lo, hi in BANDS:
        d = q[(q.bg0 >= lo) & (q.bg0 < hi)]
        if len(d) < 30:
            continue
        rows.append(dict(lo=lo, hi=hi, n=int(len(d)), fall=float(d["drop"].median())))
    return rows


def decompose(ev: pd.DataFrame) -> dict:
    """Split the fall on a correction night into what insulin did and what did not.

    Pooled as a ratio of sums rather than a median of per-event ratios. Dividing
    each event by its own small and noisy denominator both explodes the tail and
    biases the centre toward zero; summing first does neither.

    `total` is the whole fall per acting unit above what the control nights
    received, which is the quantity a settings file claims to encode. `spont` is
    the part the matched control nights delivered without a correction, and
    `causal` is the remainder. In an open loop the control nights carry basal
    only, so `spont` is genuine mean reversion; under an algorithm they carry
    whatever it decided to give, which is why `causal` collapses there.
    """
    if len(ev) < 20 or ev.extra_action_u.sum() <= 0:
        return dict(n=int(len(ev)), n_subj=int(ev.subject_id.nunique()),
                    total=np.nan, spont=np.nan, causal=np.nan,
                    causal_lo=np.nan, causal_hi=np.nan)
    ex = ev.extra_action_u.to_numpy()
    dt, dc = ev.drop_t.to_numpy(), ev.drop_c.to_numpy()
    rng = np.random.default_rng(20260826)
    b = [((dt - dc)[i].sum() / ex[i].sum())
         for i in (rng.integers(0, len(ex), len(ex)) for _ in range(2000))]
    return dict(n=int(len(ev)), n_subj=int(ev.subject_id.nunique()),
                total=float(dt.sum() / ex.sum()), spont=float(dc.sum() / ex.sum()),
                causal=float((dt - dc).sum() / ex.sum()),
                causal_lo=float(np.percentile(b, 2.5)),
                causal_hi=float(np.percentile(b, 97.5)))


def main() -> int:
    config.ensure_dirs()
    ev_all = pd.read_parquet(config.RESULTS / "inv009_events_full.parquet")
    ui = pd.read_parquet(config.RESULTS / "inv009_user_isf.parquet")
    entered = pd.read_parquet(config.RESULTS / "inv009_entered_isf.parquet")
    entered = entered[(entered.isf > 0) & (entered.tdd_u > 0)]
    ent_by = entered.set_index("subject_id").isf

    rows, spon = [], {}
    for s in STUDIES:
        w = data.load(s)
        if w.empty:
            continue
        spon[s] = spontaneous(w)
        ev = ev_all[ev_all.study == s]
        dec = decompose(ev)
        u = ui[ui.study == s]
        reg = float(u.s_pre.median()) if len(u) else np.nan
        tdd = float(u.tdd_u.median()) if len(u) else np.nan
        ent = float(ent_by.reindex(u.subject_id).dropna().median()) if len(u) else np.nan
        rows.append(dict(
            study=s, label=LABEL[s],
            closed_loop=bool(u.closed_loop.iloc[0]) if len(u) else None,
            n_people=int(len(u)), n_events=dec["n"], n_event_subjects=dec["n_subj"],
            tdd_median=tdd, isf_regression=reg,
            isf_total=dec["total"], isf_spont=dec["spont"], isf_causal=dec["causal"],
            causal_lo=dec["causal_lo"], causal_hi=dec["causal_hi"], isf_entered=ent,
            k_regression=reg * tdd, k_total=dec["total"] * tdd,
            k_causal=dec["causal"] * tdd, k_entered=ent * tdd,
            total_over_entered=dec["total"] / ent if ent else np.nan,
            causal_over_entered=dec["causal"] / ent if ent else np.nan,
        ))

    t = pd.DataFrame(rows)
    t.to_parquet(config.RESULTS / "inv009_measured_isf_by_study.parquet", index=False)

    out = dict(studies=rows, spontaneous=spon,
               note="isf_* are mg/dL per acting unit; k_* are that times median TDD")
    (config.RESULTS / "inv009_measured_isf_by_study.json").write_text(json.dumps(out, indent=1))

    print("Overnight, carb-free, isolated correction nights: what the fall is made of")
    print("(mg/dL per acting unit above what the matched control nights received)\n")
    print(f"{'study':<26s}{'ppl':>5s}{'ev':>5s}{'TOTAL':>7s}{'spont':>7s}"
          f"{'causal':>8s}{'95% CI':>17s}{'entered':>8s}")
    for r in rows:
        f = lambda v, w=7, d=1: (f"%{w}.{d}f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        ci = (f"{r['causal_lo']:7.1f} to{r['causal_hi']:6.1f}"
              if np.isfinite(r['causal_lo']) else f"{'-':>17}")
        print(f"{r['label']:<26s}{r['n_people']:5d}{r['n_events']:5d}"
              f"{f(r['isf_total'])}{f(r['isf_spont'])}{f(r['isf_causal'], 8)}{ci}"
              f"{f(r['isf_entered'], 8, 0)}")
    print("\nAs a Walsh-style constant (sensitivity x median daily dose)\n")
    print(f"{'study':<26s}{'total':>8s}{'causal':>8s}{'regression':>12s}{'entered':>9s}")
    for r in rows:
        f = lambda v, w=8: (f"%{w}.0f" % v) if np.isfinite(v) else f"{'-':>{w}}"
        print(f"{r['label']:<26s}{f(r['k_total'])}{f(r['k_causal'])}"
              f"{f(r['k_regression'], 12)}{f(r['k_entered'], 9)}")
    print("\nSpontaneous four-hour fall with no bolus given (mg/dL)")
    print("Under an algorithm this is not an untreated night: it is a night the")
    print("algorithm handled, which is why the causal column above collapses.\n")
    hdr = "".join(f"{lo}-{hi}".rjust(9) for lo, hi in BANDS)
    print(f"{'study':<26s}{hdr}")
    for s_ in STUDIES:
        if s_ not in spon:
            continue
        m = {(r['lo'], r['hi']): r['fall'] for r in spon[s_]}
        print(f"{LABEL[s_]:<26s}" + "".join(
            (f"{m[b]:9.0f}" if b in m else f"{'-':>9s}") for b in BANDS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
