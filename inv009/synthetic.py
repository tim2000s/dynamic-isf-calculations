"""Simulated people with a sensitivity law we chose, run through the real pipeline.

Everything in this investigation is an estimate of a quantity nobody observed.
The only way to know what the estimator does to a known truth is to build people
whose sensitivity really does follow a chosen law, put them through the same
window extraction and the same fits, and see what comes back.

Three things are being measured here, and only the third is a nuisance check.

Does the exponent survive? The fitted sensitivities are heavily attenuated
against entered settings, because insulin action is reconstructed rather than
observed and because meal tails contaminate the regressor. Attenuation that is
the same for everyone cancels out of a log-log slope. Attenuation that grows
with dose does not, and would bend the exponent. This says which happens.

What does a reactive controller do to it? Insulin that is chosen by watching
glucose is not independent of the thing it is supposed to explain. The simulation
can be run with the controller on and off, holding the true law fixed, so the
size of that bias is measured rather than argued about.

Does the machinery work at all? A recovered exponent near the truth on the easy
arm is the check that nothing is wired backwards.

    python3 -m inv009.synthetic
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

from . import config, insulin_models as M, tdd as tddmod, windows as windowmod

BIN = config.GRID_MIN
DAY_BINS = int(24 * 60 / BIN)


def simulate_subject(rng, tdd: float, isf: float, days: int = 200,
                     reactive: float = 0.0, noise: float = 2.5,
                     unannounced: float = 0.0, action_noise: float = 0.0) -> pd.DataFrame:
    """One person's record, on the same grid shape the real pipeline produces.

    Glucose moves by insulin action times their sensitivity, against carbohydrate
    appearing from meals and endogenous production, with a mild pull back toward
    baseline standing in for glucose effectiveness and renal clearance. The basal
    schedule is set so that its action at this person's sensitivity just offsets
    endogenous production, which is what makes a basal rate correct and is why
    the fitted intercept has to absorb it.

    Meals move around, vary in size, are dosed for imperfectly and are sometimes
    skipped, and corrections happen at bedtime. Without that variation every
    night looks the same, insulin on board at midnight barely moves, and there is
    nothing for a regression to read. Real records vary; so must these.

    `reactive` is a controller's gain, in units per hour per mg/dL above target.
    At zero, delivery is a fixed schedule plus the person's own boluses, which is
    open loop.

    `unannounced` is the share of meals eaten without being logged or dosed for,
    and `action_noise` is multiplicative error on how much of a dose actually
    acted. Together they are the harsh arm: they reproduce the heavy attenuation
    the real fits show against entered settings, and the question they answer is
    whether an exponent still survives it.
    """
    n = days * DAY_BINS
    ts = pd.date_range("2019-01-01", periods=n, freq=f"{int(BIN)}min")
    hod = (np.arange(n) % DAY_BINS) * BIN / 60.0
    day_of = np.arange(n) // DAY_BINS

    cr = 500.0 / tdd                       # grams per unit, the usual rule of thumb
    mg_per_gram = isf / cr                 # what a gram does, by definition of the two
    basal_u_day = tdd * 0.45
    shape = 1.0 + 0.35 * np.exp(-0.5 * ((hod - 5.0) / 2.0) ** 2)
    shape = shape / shape.mean()
    basal = basal_u_day / DAY_BINS * shape

    carbs = np.zeros(n)
    carbs_hidden = np.zeros(n)
    bolus = np.zeros(n)
    for base_h, size in ((7.5, 45.0), (12.5, 55.0), (19.0, 70.0)):
        for d in range(days):
            if rng.random() < 0.12:                       # skipped often enough to matter
                continue
            h = base_h + rng.normal(0, 1.1)               # meals move around
            i = d * DAY_BINS + int(np.clip(h, 0, 23.9) * 60 / BIN)
            if i >= n:
                continue
            g = float(np.clip(rng.lognormal(np.log(size), 0.42), 8, 190))
            carbs[i] += g
            if rng.random() >= unannounced:
                bolus[i] += g / cr * rng.normal(1.0, 0.18)   # dosed for, imperfectly
            else:
                carbs_hidden[i] += g                          # eaten, never recorded
    # Bedtime corrections, which is where most of the night-to-night variation in
    # insulin on board at midnight actually comes from.
    for d in range(days):
        if rng.random() < 0.35:
            i = d * DAY_BINS + int(rng.uniform(21.0, 23.5) * 60 / BIN)
            if i < n:
                bolus[i] += float(np.clip(rng.lognormal(np.log(0.05 * tdd), 0.5), 0.2, 0.3 * tdd))

    carb_kern = np.exp(-np.arange(60) * BIN / 50.0)
    carb_kern /= carb_kern.sum()
    ra = np.convolve(carbs + carbs_hidden, carb_kern)[:n] * mg_per_gram

    kern = M.kernel("oref_6h75")
    act_kern = -np.diff(np.append(kern, 0.0))             # fraction acting per bin
    egp = isf * basal                                     # what basal is there to offset
    sg = 0.014                                            # pull back to baseline per bin

    delivered = basal + bolus
    # What actually acted, which is not exactly what the model says acted.
    effective = delivered * (1.0 + action_noise * rng.normal(0, 1, n)) if action_noise > 0 \
        else delivered
    act = np.convolve(np.clip(effective, 0, None), act_kern)[:n]
    g = np.empty(n)
    g[0] = 130.0
    drift = rng.normal(0, noise, n)
    corr = np.zeros(n)
    K = len(act_kern)
    for i in range(1, n):
        g[i] = (g[i - 1] - isf * act[i - 1] + ra[i - 1] + egp[i - 1]
                - sg * (g[i - 1] - 110.0) + drift[i - 1])
        g[i] = min(max(g[i], 40.0), 400.0)
        if reactive > 0.0:
            c = max(0.0, reactive * (g[i] - 120.0)) * BIN / 60.0
            if c > 1e-9:
                corr[i] = c
                j = min(i + K, n)
                act[i:j] += c * act_kern[:j - i]          # its action, scattered forward
    total = delivered + corr

    # Unannounced carbohydrate is eaten but never appears in the record, which is
    # the point: the screen cannot see it.
    return pd.DataFrame({"ts": ts, "cgm": g, "basal_u": basal + corr, "bolus_u": bolus,
                         "sched_u": basal, "carbs_g": carbs, "total_u": total})


def run_cohort(beta_true: float, n_subjects: int = 60, reactive: float = 0.0,
               seed: int = 0, isf_at_40: float = 50.0, log_carbs: bool = True,
               unannounced: float = 0.0, action_noise: float = 0.0,
               tag: str = "") -> dict:
    """Build a cohort on a known law and recover the exponent from it.

    The constant is anchored so that a person on forty units a day has the same
    sensitivity whatever exponent is being simulated. Without that, steep and
    shallow laws are not comparable: holding the constant fixed instead would
    give the shallow arm a sensitivity of 664 mg/dL/U at ten units a day, which
    is not a person, and the arm would fail for that reason rather than because
    the exponent was hard to recover.
    """
    rng = np.random.default_rng(seed)
    k_true = isf_at_40 / (40.0 ** beta_true)
    tdds = np.exp(rng.uniform(np.log(10), np.log(100), n_subjects))
    rows, truth = [], []
    for i, tdd in enumerate(tdds):
        isf = k_true * tdd ** beta_true
        grid = simulate_subject(rng, tdd, isf, reactive=reactive,
                                unannounced=unannounced, action_noise=action_noise)
        if not log_carbs:
            grid = grid.assign(carbs_g=0.0)
        subj = tddmod.subject_level(grid)
        tw = tddmod.windowed(grid)
        w = windowmod.extract_windows(grid, tw, f"SIM:{i}", "ReplaceBG", 40.0,
                                      subj, ("oref_6h75",))
        if w.empty:
            continue
        w = w[windowmod.screen(w)]
        if len(w) < config.MIN_WINDOWS:
            continue
        from .effective_isf import _fit_subject
        from .data import attach_action
        w = attach_action(w.copy())
        f = _fit_subject(w, "a_pre")
        if not f:
            continue
        rows.append(dict(subject_id=f"SIM:{i}", tdd_u=subj["tdd_u"], s=f["s"],
                         se=f["se"], n=f["n"]))
        truth.append(dict(tdd=tdd, isf=isf))
    R = pd.DataFrame(rows)
    T = pd.DataFrame(truth)
    from . import stats as st
    pos = R[R.s > 0]
    r = st.loglog(pos.tdd_u, pos.s)
    ci = st.boot_slope_ci(pos.tdd_u, pos.s)
    return dict(beta_true=beta_true, reactive=reactive, log_carbs=log_carbs,
                tag=tag, unannounced=unannounced, action_noise=action_noise,
                n_fitted=int(len(R)), n_positive=int(len(pos)),
                frac_positive=float((R.s > 0).mean()) if len(R) else np.nan,
                recovered=r["slope"], ci=[ci[0], ci[1]],
                bias=float(r["slope"] - beta_true) if np.isfinite(r["slope"]) else np.nan,
                attenuation=float(np.median(pos.s) / np.median(T.isf)) if len(pos) else np.nan,
                truth_isf_median=float(T.isf.median()) if len(T) else np.nan,
                fitted_s_median=float(pos.s.median()) if len(pos) else np.nan)


def main() -> int:
    config.ensure_dirs()
    out = []
    print("Recovering a known exponent through the real pipeline\n")
    print(f"{'truth':>7s} {'controller':>11s} {'carbs':>6s} {'n+':>7s} {'recovered':>22s} "
          f"{'bias':>7s} {'atten':>6s}")
    for beta in (-0.5, -1.0, -2.0):
        for reactive, tag in ((0.0, "open"), (0.02, "reactive")):
            r = run_cohort(beta, reactive=reactive)
            out.append(r)
            print(f"{beta:+7.2f} {tag:>11s} {'yes':>6s} {r['n_positive']:3d}/{r['n_fitted']:3d} "
                  f"{r['recovered']:+8.3f} [{r['ci'][0]:+.2f},{r['ci'][1]:+.2f}] "
                  f"{r['bias']:+7.3f} {r['attenuation']:6.2f}")
    # And once with carbohydrate hidden, which is what the cohorts without a carb
    # record are actually working with.
    r = run_cohort(-1.0, reactive=0.0, log_carbs=False, tag="carbs hidden")
    out.append(r)
    print(f"{-1.0:+7.2f} {'open':>11s} {'no':>6s} {r['n_positive']:3d}/{r['n_fitted']:3d} "
          f"{r['recovered']:+8.3f} [{r['ci'][0]:+.2f},{r['ci'][1]:+.2f}] "
          f"{r['bias']:+7.3f} {r['attenuation']:6.2f}")

    # The harsh arm. The real fits are attenuated far harder than anything above,
    # so the question is whether an exponent survives that much contamination.
    print("\n  harsh arm: 45% of meals never recorded, 35% error on what acted")
    for beta in (-0.5, -1.0, -2.0):
        r = run_cohort(beta, reactive=0.0, unannounced=0.45, action_noise=0.35,
                       tag="harsh")
        out.append(r)
        print(f"{beta:+7.2f} {'open':>11s} {'yes':>6s} {r['n_positive']:3d}/{r['n_fitted']:3d} "
              f"{r['recovered']:+8.3f} [{r['ci'][0]:+.2f},{r['ci'][1]:+.2f}] "
              f"{r['bias']:+7.3f} {r['attenuation']:6.2f}")
    (config.RESULTS / "inv009_synthetic.json").write_text(json.dumps(out, indent=2, default=float))
    print("\n  attenuation is the fitted sensitivity as a fraction of the true one;")
    print("  bias is what the pipeline adds to the exponent, and is the number that")
    print("  decides whether the measured exponent can be read at face value")
    return 0


if __name__ == "__main__":
    sys.exit(main())
