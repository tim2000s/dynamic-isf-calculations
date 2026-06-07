#!/usr/bin/env python3
"""Is there a better ISF~TDD equation than v1 or v2?

Pulls together every ISF~TDD relationship examined (the historical 1700-rule, a
ΔIOB-based estimate of observed sensitivity, a TDD-band lookup, and the v1/v2
equations) and evaluates candidate equations side by side with leave-one-user-out
cross-validation, so fitted forms are scored out-of-sample and comparable with the
fixed rules.

Targets:
  * empirical ISF (ΔIOB-derived observed sensitivity; n≈114) — "what actually happens"
  * entered ISF (user-tuned profile value; n=138)            — "what users converge to"

Candidates:
  fixed:   v1 (TDD^-1), v2 (TDD^-2), 1700-rule, entered profile ISF
  fitted:  C/TDD (re-fitted constant), power law A·TDD^b, power law + basal fraction,
           power law + ln CR, full multivariate, TDD-quartile band lookup,
           geometric blend of entered ISF with the fitted power law (weight CV-fitted)

Output: best_isf_fit_results.json / .md, charts/inv008/fig_best_fit.png
"""
from __future__ import annotations

import json
import math
from datetime import datetime
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("DYNISF_ROOT", Path.cwd()))
LOG_TERM = math.log(99.0 / 75 + 1.0)  # normal target 99, divisor 75 (Lyumjev)

# ---------------------------------------------------------------- data

def load_cohort() -> pd.DataFrame:
    coh = pd.DataFrame(json.loads((ROOT / "canonical_cohort.json").read_text()))
    emp = pd.DataFrame(json.loads((ROOT / "empirical_isf_v5.json").read_text()))
    emp = emp[["user_id", "empirical_isf", "r2", "n_windows"]]
    df = coh.merge(emp, on="user_id", how="left")
    df["basal_frac"] = df["basal"] / df["tdd"]
    df["emp_valid"] = (df["r2"] >= 0.10) & df["empirical_isf"].between(5, 500)
    return df


# ---------------------------------------------------------------- candidates
# Each candidate: fit(train_df, ycol) -> predict(test_df) -> np.ndarray

def cand_v1(train, ycol):
    return lambda t: 1800.0 / (t["tdd"].to_numpy() * LOG_TERM)

def cand_v2(train, ycol):
    return lambda t: 2300.0 / (LOG_TERM * t["tdd"].to_numpy() ** 2 * 0.02)

def cand_walsh(train, ycol):
    return lambda t: 1700.0 / t["tdd"].to_numpy()

def cand_entered(train, ycol):
    return lambda t: t["isf"].to_numpy()

def cand_c_over_tdd(train, ycol):
    c = float(np.median(train[ycol] * train["tdd"]))
    return lambda t: c / t["tdd"].to_numpy()

def _ols_log(train, ycol, cols):
    X = np.column_stack([np.ones(len(train))] +
                        [train[c].to_numpy(dtype=float) for c in cols])
    y = np.log(train[ycol].to_numpy(dtype=float))
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    def predict(t):
        Xt = np.column_stack([np.ones(len(t))] +
                             [t[c].to_numpy(dtype=float) for c in cols])
        return np.exp(Xt @ beta)
    return predict, beta

def cand_power(train, ycol):
    train = train.assign(ln_tdd=np.log(train["tdd"]))
    pred, _ = _ols_log(train, ycol, ["ln_tdd"])
    return lambda t: pred(t.assign(ln_tdd=np.log(t["tdd"])))

def cand_power_basal(train, ycol):
    tr = train.assign(ln_tdd=np.log(train["tdd"]))
    pred, _ = _ols_log(tr, ycol, ["ln_tdd", "basal_frac"])
    return lambda t: pred(t.assign(ln_tdd=np.log(t["tdd"])))

def cand_power_cr(train, ycol):
    tr = train.assign(ln_tdd=np.log(train["tdd"]), ln_cr=np.log(train["cr"]))
    pred, _ = _ols_log(tr, ycol, ["ln_tdd", "ln_cr"])
    return lambda t: pred(t.assign(ln_tdd=np.log(t["tdd"]), ln_cr=np.log(t["cr"])))

def cand_full(train, ycol):
    tr = train.assign(ln_tdd=np.log(train["tdd"]), ln_cr=np.log(train["cr"]),
                      ln_tgt=np.log(train["target_low"]))
    pred, _ = _ols_log(tr, ycol, ["ln_tdd", "ln_cr", "basal_frac", "ln_tgt"])
    return lambda t: pred(t.assign(ln_tdd=np.log(t["tdd"]), ln_cr=np.log(t["cr"]),
                                   ln_tgt=np.log(t["target_low"])))

def cand_bands(train, ycol):
    qs = np.quantile(train["tdd"], [0.25, 0.5, 0.75])
    meds = []
    edges = [-np.inf, *qs, np.inf]
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = train[(train["tdd"] > lo) & (train["tdd"] <= hi)][ycol]
        meds.append(float(np.median(m)) if len(m) else float(np.median(train[ycol])))
    def predict(t):
        idx = np.searchsorted(qs, t["tdd"].to_numpy(), side="left")
        return np.array([meds[i] for i in idx])
    return predict

def cand_sqrt(train, ycol):
    """Fixed exponent −0.5; only the constant K is fitted: ISF = K/√TDD."""
    k = float(np.median(train[ycol] * np.sqrt(train["tdd"])))
    def predict(t):
        return k / np.sqrt(t["tdd"].to_numpy())
    predict.k = k
    return predict


def cand_blend(train, ycol):
    """Geometric blend entered^w · powerlaw^(1−w); w chosen by inner LOUO on train."""
    tr = train.assign(ln_tdd=np.log(train["tdd"]))
    pred_pl, _ = _ols_log(tr, ycol, ["ln_tdd"])
    pl_train = pred_pl(tr)
    best_w, best_err = 0.0, np.inf
    for w in np.linspace(0, 1, 11):
        p = train["isf"].to_numpy() ** w * pl_train ** (1 - w)
        err = float(np.median(np.abs(p - train[ycol].to_numpy())))
        if err < best_err:
            best_w, best_err = float(w), err
    def predict(t):
        pl = pred_pl(t.assign(ln_tdd=np.log(t["tdd"])))
        return t["isf"].to_numpy() ** best_w * pl ** (1 - best_w)
    predict.w = best_w
    return predict


CANDIDATES = {
    "v1 (TDD^-1)":                    (cand_v1, False),
    "v2 (TDD^-2)":                  (cand_v2, False),
    "1700-rule":                     (cand_walsh, False),
    "Entered profile ISF":           (cand_entered, False),
    "Fitted C/TDD":                  (cand_c_over_tdd, True),
    "K/sqrt(TDD)":                   (cand_sqrt, True),
    "Power law A·TDD^b":             (cand_power, True),
    "Power law + basal_frac":        (cand_power_basal, True),
    "Power law + ln(CR)":            (cand_power_cr, True),
    "Multivariate (TDD,CR,basal,target)": (cand_full, True),
    "TDD-quartile bands":            (cand_bands, True),
    "Blend entered×power law":       (cand_blend, True),
}


def louo_eval(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
    rows = []
    y = df[ycol].to_numpy()
    for name, (factory, fitted) in CANDIDATES.items():
        if ycol == "isf" and name in ("Entered profile ISF", "Blend entered×power law"):
            continue  # circular against the entered target
        preds = np.full(len(df), np.nan)
        if fitted:
            for i in range(len(df)):
                train = df.drop(df.index[i])
                predict = factory(train, ycol)
                preds[i] = float(predict(df.iloc[[i]])[0])
        else:
            preds = factory(df, ycol)(df)
        ae = np.abs(preds - y)
        loge = np.abs(np.log(preds / y))
        rows.append({
            "candidate": name, "fitted_cv": fitted,
            "median_abs_err": round(float(np.median(ae)), 1),
            "p75_abs_err": round(float(np.percentile(ae, 75)), 1),
            "median_log_err": round(float(np.median(loge)), 3),
            "frac_within_30pct": round(float((loge < math.log(1.3)).mean()), 3),
        })
    return pd.DataFrame(rows).sort_values("median_log_err")


def _md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def main() -> None:
    df = load_cohort()
    emp = df[df["emp_valid"]].copy()
    print(f"cohort n={len(df)}, empirical-valid n={len(emp)}")

    res_emp = louo_eval(emp, "empirical_isf")
    res_ent = louo_eval(df.dropna(subset=["isf"]), "isf")

    # fitted parameters on the full data, for the writeup (guard: TDD must be > 0)
    tr = emp[emp["tdd"] > 0].assign(ln_tdd=np.log(emp[emp["tdd"] > 0]["tdd"]))
    _, beta = _ols_log(tr, "empirical_isf", ["ln_tdd"])
    power_emp = {"A": round(float(np.exp(beta[0])), 1), "b": round(float(beta[1]), 3)}
    ent = df.dropna(subset=["isf"])
    ent = ent[ent["tdd"] > 0]
    tr2 = ent.assign(ln_tdd=np.log(ent["tdd"]))
    _, beta2 = _ols_log(tr2, "isf", ["ln_tdd"])
    power_ent = {"A": round(float(np.exp(beta2[0])), 1), "b": round(float(beta2[1]), 3)}
    blend = cand_blend(emp, "empirical_isf")
    sqrt_emp = cand_sqrt(emp, "empirical_isf")
    sqrt_ent = cand_sqrt(ent, "isf")

    out = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "n_cohort": int(len(df)), "n_empirical": int(len(emp)),
        "power_law_empirical": power_emp,
        "power_law_entered": power_ent,
        "blend_weight_on_entered": blend.w,
        "sqrt_rule_empirical_K": round(sqrt_emp.k, 1),
        "sqrt_rule_entered_K": round(sqrt_ent.k, 1),
        "target_empirical_isf": res_emp.to_dict("records"),
        "target_entered_isf": res_ent.to_dict("records"),
    }
    (ROOT / "best_isf_fit_results.json").write_text(json.dumps(out, indent=1))

    md = [f"# Best-fit ISF equation search — LOUO-CV comparison",
          f"\n{out['generated']} · empirical target n={len(emp)} · entered target n={len(df)}",
          f"\nFitted power law (empirical): ISF = {power_emp['A']}·TDD^{power_emp['b']}",
          f"Fitted power law (entered):  ISF = {power_ent['A']}·TDD^{power_ent['b']}",
          f"Blend weight on entered ISF: {blend.w:.1f}",
          "\n## Target: empirical ISF (observed sensitivity)\n",
          _md_table(res_emp),
          "\n## Target: entered ISF (user-tuned profile)\n",
          _md_table(res_ent)]
    (ROOT / "best_isf_fit_results.md").write_text("\n".join(md))

    # figure: median log error per candidate, both targets
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, res, title in ((axes[0], res_emp, f"target: empirical ISF (n={len(emp)})"),
                           (axes[1], res_ent, f"target: entered ISF (n={len(df)})")):
        r = res.sort_values("median_log_err", ascending=True)
        colors = ["#d62728" if "V2" in c else "#1f77b4" if "V1" in c
                  else "#7f7f7f" if not f else "#2ca02c"
                  for c, f in zip(r["candidate"], r["fitted_cv"])]
        ax.barh(r["candidate"], r["median_log_err"], color=colors)
        ax.set_xlabel("median |log error| (lower = better)")
        ax.set_title(title)
        ax.invert_yaxis()
        ax.grid(alpha=0.3, axis="x")
    fig.suptitle("Candidate ISF equations, leave-one-user-out CV "
                 "(green = fitted, grey/blue/red = fixed rules)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(ROOT / "charts/inv008/fig_best_fit.png", dpi=150)

    print(res_emp.to_string(index=False))
    print()
    print(res_ent.to_string(index=False))


if __name__ == "__main__":
    main()
