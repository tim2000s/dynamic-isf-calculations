"""INV-020: the glucose curve against the logarithmic model at every recorded loop decision.

At each decision the loop reads a sensitivity at the current glucose, predicts where glucose
will end up and asks for insulin from that prediction. This replays every decision of every
person in the OpenAPS Data Commons with a daily-dose reconstruction (AndroidAPS, oref_v6, and
oref0, oref_v7) twice: once with the logarithmic model and once with the curve of
docs/dynamic-isf-8-proposal.md (profile sensitivity x Q(g)/Q(140),
`inv013.candidate.curve_profile`). It reports, across people and by glucose band, whether the
curve's sensitivity is higher or lower than the logarithmic model's, and how the prediction and
the insulin asked for move with it. Per-person results are computed for the summaries and
written only to the git-ignored inv020_cache/, never to the repository.
Nothing is simulated forward. Each decision keeps its recorded state, so this says what each
model would have asked for at that moment and not where glucose would have gone next.

How the prediction and the request move. In oref (determine-basal.js) the eventual glucose is
    naive eventual glucose + deviation = bg - IOB x S_iob + 6 x minDelta + 30 x activity x ISF
where S_iob is the ISF when IOB is above zero and min(ISF, profile ISF) otherwise. Everything
but the two ISF terms is the same for both models, so each model's eventual glucose is the
logged one with the logged ISF terms taken out and its own put in. minPredBG is shifted by the
same amount, a modelling choice: its curves are built from the same insulin activity. The
choice between minDelta and the average deltas depends on the sign of the deviation, which
depends on the ISF, so near zero deviation the shift is approximate. The request follows oref:
(min(eventual glucose, minPredBG) - target) / ISF with eventual glucose at or above target, and
2 x min(0, (eventual glucose - target) / ISF) below it, taking the logged target as the low
target, a modelling choice. That second form is oref's low temporary basal and is never logged,
so only the first can be checked: at the logged ISF it reproduces the logged request on the
share of decisions inv020/checks.py reports. The branches where oref holds the scheduled basal
in range are not modelled. A positive request is what oref answers with an SMB (half of it by
default) or a high temporary basal; a negative one cuts the temporary basal. No dosing limit
is applied.

Inputs and modelling choices, fixed before the run:
  Commons. Decisions from oref_v6 and oref_v7. A timestamp whose rows disagree is dropped
  (REVIEW.md item 8); identical repeats count once. The sensitivity the prediction was made with
  is the logged sug_isf (mg/dL per U; values under 20 were checked and are mg/dL: the request
  identity holds in mg/dL and fails in mmol/L). oref sets sens = profile sens / sensitivityRatio,
  so the profile sensitivity is sug_isf x sensitivityRatio. Daily dose is the INV-008
  reconstruction, rebuilt with the readers corrected on 23 November 2026 (AndroidAPS temp basals
  read, each bolus counted once; inv008/sources.py), joined on ts_relative_sec. The logarithmic arm is
  the AndroidAPS equation at 100% with divisor 55, as in INV-016, without autosens or its
  limits, and inherits the reconstruction's remaining errors. Nobody in the Commons ran it, so
  both arms are counterfactual there.
  Both arms omit autosens. Decisions whose logged prediction is at oref's clamp (39 or 401
  mg/dL) are dropped, since the clamped value hides the prediction. Glucose 40 to 400 mg/dL. A
  person enters with 288 decisions (one day). Intervals are 95%: per person by resampling days,
  across people by resampling people, 2,000 draws. A verdict is given where the interval
  excludes no change, and the result is called unresolved where it does not.

Expected before the run, from the curve's shape and INV-016:
  Commons. The curve's sensitivity is above the logarithmic model's at almost every decision, by
  about 2.4 times near 100 mg/dL and 1.1 times at 300 mg/dL for a person whose profile x daily
  dose is the cohort median (about 2,600), because Q(g) x ln(g/55 + 1) falls with glucose. So the
  curve asks for less insulin at most decisions, the gap largest near target and narrowing with
  glucose; people whose product is under about 1,600 should cross over and get more from the
  curve above about 250 mg/dL. The prediction: with IOB above 30 x activity the curve's higher
  sensitivity predicts a lower eventual glucose.
  The figures above were written from INV-016, whose daily doses carry the repeated boluses.
  The curve-to-log ratio of sensitivity is proportional to the daily dose, so the corrected
  doses should lower it in proportion.

The run of 23 November 2026 (148 people; numbers from results/inv020_replay.md and
results/inv020_checks.md). The median profile x daily dose was
2,217 (AndroidAPS) and 2,282 (oref0), not 2,600, so the expected ratios become about 2.0 near
100 mg/dL and 1.0 at 300. As expected, the curve read the higher sensitivity at most decisions
in every cohort, 1.67 to 1.75 times the log model's from 80 to 150 mg/dL and less above, and it
asked for less insulin: the positive request was 0.50 and 0.57 of the log model's, with
intervals below 1. Above 260 mg/dL the ratio was 1.05 and 1.15 rather than about 1.0, and the
below-80 band sat lower than 80 to 100. More people crossed over than the 7 and 21 with a
product under 1,600 (16 of 39 and 48 of 109, at a median of 170 and 210 mg/dL). Neither
departure was looked into further.

    python3 -m inv020.replay [--pilot]       writes results/inv020_replay.md and charts/inv020/
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

from inv008.config import DB
from inv013.candidate import curve_profile, stock_v1

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CHARTS = ROOT / "charts" / "inv020"
REPLAY = (
    Path.home() / "oref-investigations" / "inv008_cache" / "replay"
)  # rebuilt 23 Nov 2026, corrected readers
DIVISOR = 55.0
MIN_DECISIONS = 288
N_BOOT = 2000
WORKERS = 7
G_EDGES = (40, 80, 100, 150, 200, 260, 401)
G_LABELS = (
    "below 80",
    "80 to 100",
    "100 to 150",
    "150 to 200",
    "200 to 260",
    "above 260",
)
COHORTS = {
    "v6": "AndroidAPS (Commons)",
    "v7": "oref0 (Commons)",
}
PRIVATE = ROOT / "inv020_cache"  # git-ignored: per-person output stays here


def isf_terms(s, profile, iob, act):
    """The part of oref's eventual glucose that depends on the ISF: 30 x activity x ISF from the
    deviation, minus IOB x ISF, the ISF capped at the profile's when IOB is not above zero."""
    return 30.0 * act * s - iob * np.where(iob > 0, s, np.minimum(s, profile))


def changes(target, ev, mp, isf_pred, profile, iob, act, isf_log, isf_curve) -> dict:
    """Prediction and request under each arm, from one decision's logged values (arrays)."""
    logged = isf_terms(isf_pred, profile, iob, act)
    out = {}
    for name, s in (("log", isf_log), ("curve", isf_curve)):
        shift = isf_terms(s, profile, iob, act) - logged
        e, m = ev + shift, mp + shift
        out["p_" + name] = np.fmin(e, m)
        out["r_" + name] = np.where(
            e < target,
            2 * np.minimum(0, (e - target) / s),
            (np.fmin(e, m) - target) / s,
        )
    return out


def one_per_timestamp(d: pd.DataFrame, ts: str) -> pd.DataFrame:
    """Identical repeats once; a timestamp whose rows disagree is dropped (REVIEW.md item 8)."""
    d = d.drop_duplicates()
    return d[~d[ts].duplicated(keep=False)]


def commons(job):
    import psycopg2

    user, table, path = job
    tdd = pd.read_parquet(path, columns=["ts_relative_sec", "tdd"]).drop_duplicates(
        "ts_relative_sec"
    )
    with psycopg2.connect(**DB) as conn:
        d = pd.read_sql(
            f"""SELECT ts_relative_sec, sug_bg AS bg, sug_current_target AS target,
                   sug_isf AS isf_pred, sug_sensitivityratio AS sr, sug_eventualbg AS ev, reason_minpredbg AS mp,
                   iob_iob AS iob, iob_activity AS act
                FROM {table} WHERE user_id = %(u)s""",
            conn,
            params={"u": user},
        )
    d = one_per_timestamp(d, "ts_relative_sec").merge(
        tdd, on="ts_relative_sec", how="inner"
    )
    d["day"] = (d.ts_relative_sec // 86400).astype(int)
    d["profile"] = d.isf_pred * d.sr.fillna(1.0)
    d["isf_log"] = stock_v1(d.bg.to_numpy(float), d.tdd.to_numpy(float), DIVISOR)
    product = float((d.profile * d.tdd).median())
    return finish(d, user, "v6" if table == "oref_v6" else "v7", product)


def finish(d: pd.DataFrame, user: str, cohort: str, product: float):
    d = d[
        d.bg.between(40, 400)
        & (d.isf_pred > 0)
        & (d.isf_log > 0)
        & (d.profile > 0)
        & d[["target", "ev", "iob", "act"]].notna().all(axis=1)
    ].copy()
    d["mp"] = d.mp.fillna(d.ev)
    p = np.fmin(d.ev, d.mp)
    d = d[(p > 39) & (p < 401)]  # oref clamps its predictions there
    if len(d) < MIN_DECISIONS:
        return None
    d["isf_curve"] = curve_profile(d.bg.to_numpy(float), d.profile.to_numpy(float))
    d = d.assign(
        **changes(
            *(
                d[c].to_numpy(float)
                for c in (
                    "target",
                    "ev",
                    "mp",
                    "isf_pred",
                    "profile",
                    "iob",
                    "act",
                    "isf_log",
                    "isf_curve",
                )
            )
        )
    )
    d["g"] = np.digitize(d.bg, G_EDGES) - 1
    keep = [
        "day",
        "g",
        "bg",
        "profile",
        "isf_log",
        "isf_curve",
        "p_log",
        "p_curve",
        "r_log",
        "r_curve",
    ]
    return user, cohort, product, d[keep].astype({"day": "int32", "g": "int8"})


def day_ratio(
    d: pd.DataFrame, num: pd.Series, den: pd.Series, rng
) -> tuple[float, float, float]:
    """sum(num) / sum(den) with a 95% interval from resampling whole days."""
    days = (
        pd.DataFrame({"n": num, "d": den, "day": d.day}).groupby("day").sum().to_numpy()
    )
    k = rng.integers(0, len(days), (N_BOOT, len(days)))
    draws = days[k, 0].sum(axis=1) / days[k, 1].sum(axis=1)
    return days[:, 0].sum() / days[:, 1].sum(), *np.percentile(draws, [2.5, 97.5])


def verdict(ci, none: float, low: str, high: str) -> str:
    return high if ci[1] > none else low if ci[2] < none else "unresolved"


def per_person(user, cohort, product, d, rng) -> dict:
    ratio = d.isf_curve / d.isf_log
    one = pd.Series(1.0, index=d.index)
    pos = day_ratio(d, np.maximum(d.r_curve, 0), np.maximum(d.r_log, 0), rng)
    return {
        "person": user,
        "cohort": cohort,
        "decisions": len(d),
        "days": d.day.nunique(),
        "product": product,
        "profile": d.profile.median(),
        "isf_log": d.isf_log.median(),
        "isf_curve": d.isf_curve.median(),
        "isf_ratio": ratio.median(),
        "curve_higher": (ratio > 1).mean(),
        "dp": day_ratio(d, d.p_curve - d.p_log, one, rng),
        "less_insulin": (d.r_curve < d.r_log).mean(),
        "dr": day_ratio(d, d.r_curve - d.r_log, one, rng),
        "smb_log": np.maximum(d.r_log, 0).mean() / 2,
        "smb_curve": np.maximum(d.r_curve, 0).mean() / 2,
        "pos_ratio": pos,
        "pos_verdict": verdict(pos, 1.0, "less insulin", "more insulin"),
        "on_log": (d.r_log > 0).mean(),
        "on_curve": (d.r_curve > 0).mean(),
        "cross": cross_glucose(d),
    }


def cross_glucose(d) -> float:
    """Lowest 20 mg/dL band at or above 100 mg/dL from which the curve's sensitivity is below the
    log model's at most decisions in that band and every band above it with 20 decisions or more."""
    hi = d[d.bg >= 100]
    b = (hi.isf_curve < hi.isf_log).groupby((hi.bg // 20) * 20).agg(["mean", "size"])
    above = b[b["size"] >= 20]["mean"] > 0.5
    tail = (
        above[::-1].cumprod()[::-1].astype(bool)
    )  # true where this band and all above it are
    return float(tail[tail].index.min()) if tail.any() else np.nan


def med_ci(x, rng):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    draws = np.median(x[rng.integers(0, len(x), (N_BOOT, len(x)))], axis=1)
    return float(np.median(x)), *np.percentile(draws, [2.5, 97.5])


def fmt(t, p=2):
    return f"{t[0]:.{p}f} ({t[1]:.{p}f} to {t[2]:.{p}f})"


def report(people: pd.DataFrame, bands: pd.DataFrame) -> str:
    rng = np.random.default_rng(0)
    L = [
        "# INV-020: the glucose curve against the logarithmic model at every decision",
        "",
        "Generated by `inv020/replay.py`, whose docstring holds the method, the modelling choices and",
        "what was expected before the run. Every recorded decision is replayed with each model's",
        "sensitivity at the glucose of the decision, all else as recorded. Request is oref's insulin",
        "requirement before any limit. Nothing is simulated forward. Checks on the method are in",
        "results/inv020_checks.md.",
        "",
        "## Across people",
        "",
        "Median over people with a 95% interval from resampling people. A higher sensitivity lowers",
        "a positive request and shrinks a negative one (a cut to the temporary basal), so a lower request",
        "under the curve is less insulin where the request is positive and more where it is negative.",
        "The verdict is on the positive request, the insulin oref would add.",
        "",
        "| Cohort | People | Decisions | Profile ISF x daily dose | Curve ISF / log ISF | Decisions with curve ISF higher | Change in predicted glucose, mg/dL | Positive request, curve / log | Verdict | People: less / more / unresolved | People with a crossover, median mg/dL |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|",
    ]
    for c, name in COHORTS.items():
        p = people[people.cohort == c]
        if p.empty:
            continue
        pos = med_ci([r[0] for r in p.pos_ratio], rng)
        prod = "-" if p["product"].isna().all() else f"{p['product'].median():.0f}"
        L.append(
            f"| {name} | {len(p)} | {p.decisions.sum():,} | {prod} | {fmt(med_ci(p.isf_ratio, rng))} | "
            f"{fmt(med_ci(100 * p.curve_higher, rng), 0)}% | {fmt(med_ci([r[0] for r in p.dp], rng), 1)} | "
            f"{fmt(pos)} | {verdict(pos, 1.0, 'less insulin', 'more insulin')} | "
            + " / ".join(
                str((p.pos_verdict == v).sum())
                for v in ("less insulin", "more insulin", "unresolved")
            )
            + f" | {p.cross.notna().sum()}, {p.cross.median():.0f} |"
        )
    L += [
        "",
        "## By glucose at the decision",
        "",
        "Per person, the median of curve ISF / log ISF and the mean change in request (U) in the band; "
        "then the median over people with a 95% interval.",
        "",
    ]
    for c, name in COHORTS.items():
        sub = bands[bands.cohort == c]
        if sub.empty:
            continue
        L += [
            f"### {name}",
            "",
            "| Glucose, mg/dL | People | Curve ISF / log ISF | Change in request per decision, U | Verdict on the request |",
            "|---|---:|---:|---:|---|",
        ]
        for i, lab in enumerate(G_LABELS):
            s = sub[(sub.g == i) & (sub.n >= 20)]
            if len(s) >= 3:
                dr = med_ci(s.dr, rng)
                L.append(
                    f"| {lab} | {len(s)} | {fmt(med_ci(s.ratio, rng))} | {fmt(dr)} | "
                    f"{verdict(dr, 0.0, 'lower under the curve', 'higher under the curve')} |"
                )
        L.append("")
    return "\n".join(L) + "\n"


def work(job):
    out = commons(job)
    if out is None:
        return None
    user, cohort, product, d = out
    rng = np.random.default_rng(zlib.crc32(user.encode()))
    band = (
        d.assign(ratio=d.isf_curve / d.isf_log, dr=d.r_curve - d.r_log)
        .groupby("g")
        .agg(ratio=("ratio", "median"), dr=("dr", "mean"), n=("dr", "size"))
        .reset_index()
        .assign(person=user, cohort=cohort)
    )
    return per_person(user, cohort, product, d, rng), band


def jobs() -> list:
    out = []
    for f in sorted(REPLAY.glob("*.parquet")):
        meta = json.loads(f.with_suffix(".meta.json").read_text())
        if meta["platform"] in ("v6", "v7"):
            out.append((f.stem, "oref_" + meta["platform"], str(f)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--pilot",
        action="store_true",
        help="one person from each cohort, report printed and nothing written",
    )
    a = ap.parse_args()
    js = jobs()
    if a.pilot:
        js = [
            next(j for j in js if j[1] == "oref_v6"),
            next(j for j in js if j[1] == "oref_v7"),
        ]
    t = time.time()
    first = work(js[0])  # one worker timed before the pool starts
    print(f"{len(js)} people; the first took {time.time() - t:.1f} s", flush=True)
    with mp.Pool(WORKERS) as pool:
        rest = pool.map(work, js[1:])
    res = [r for r in [first, *rest] if r is not None]
    people = pd.DataFrame([r[0] for r in res])
    bands = pd.concat([r[1] for r in res])
    print(f"{len(people)} people kept of {len(js)}; {time.time() - t:.0f} s")
    if a.pilot:
        print(report(people, bands))
        return
    # the per-person output goes to disk before the report is formatted (REVIEW.md item 16), in the
    # git-ignored cache: per-person results are not published
    PRIVATE.mkdir(exist_ok=True)
    flat = people.copy()
    for col in ("pos_ratio", "dp", "dr"):
        flat[col], flat[col + "_lo"], flat[col + "_hi"] = zip(*people[col])
    flat.to_csv(PRIVATE / "inv020_people.csv", index=False)
    bands.to_csv(PRIVATE / "inv020_bands.csv", index=False)
    (RESULTS / "inv020_replay.md").write_text(report(people, bands))
    from inv020.figures import main as figures

    figures()
    print(RESULTS / "inv020_replay.md")


if __name__ == "__main__":
    main()
