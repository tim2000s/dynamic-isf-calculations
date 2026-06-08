"""Stage 3 — per-user pages and cohort figures for the V1 vs V2 ISF comparison.

Usage:
    python -m inv008.stage3_plots                  # all users + cohort figs, 12 workers
    python -m inv008.stage3_plots --users U073 ... # subset
    python -m inv008.stage3_plots --cohort-only

Output: charts/inv008/users/<uid>.png, charts/inv008/fig_*.png, cohort_summary.json
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from inv008 import config
from inv008.dynisf import isf_v1, isf_v2

CHART_DIR = config.ROOT / "charts" / "inv008"
USER_DIR = CHART_DIR / "users"

C_V1 = "#1f77b4"   # blue
C_V2 = "#d62728"   # red
C_BG = "#bbbbbb"
C_EMP = "#2ca02c"  # green

EMPIRICAL = {r["user_id"]: r for r in
             json.loads((config.ROOT / "empirical_isf_v5.json").read_text())} \
    if (config.ROOT / "empirical_isf_v5.json").exists() else {}


def _load(user_id: str):
    df = pd.read_parquet(config.REPLAY_DIR / f"{user_id}.parquet")
    meta = json.loads((config.REPLAY_DIR / f"{user_id}.meta.json").read_text())
    # recompute v2 from the current equation so figures match the live implementation
    df["isf_v2"] = isf_v2(df["bg"].to_numpy(), df["tdd"].to_numpy())
    return df[np.isfinite(df["isf_v1"]) & np.isfinite(df["isf_v2"])].copy(), meta


def _densest_window(df: pd.DataFrame, days: int = 14) -> pd.DataFrame:
    ts = df["ts_relative_sec"].to_numpy()
    span = days * 86400
    if ts[-1] - ts[0] <= span:
        return df
    # slide in 1-day steps; pick the window with the most ticks
    starts = np.arange(ts[0], ts[-1] - span, 86400)
    counts = np.searchsorted(ts, starts + span) - np.searchsorted(ts, starts)
    s = starts[np.argmax(counts)]
    return df[(df["ts_relative_sec"] >= s) & (df["ts_relative_sec"] < s + span)]


def render_user_page(args: tuple[str,]) -> dict:
    user_id, = args
    try:
        df, meta = _load(user_id)
        if len(df) < 500:
            return {"user": user_id, "status": "skip", "reason": f"{len(df)} valid ticks"}
        tdd_med = float(np.nanmedian(df["tdd"]))
        tdd_p25, tdd_p75 = np.nanpercentile(df["tdd"], [25, 75])

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
        fig.suptitle(
            f"{user_id} ({meta['platform']}) — dynamic ISF, v1 vs v2   |   "
            f"median TDD {tdd_med:.0f} U/day   |   V2/V1 ratio {63.89/tdd_med:.2f}",
            fontsize=12)

        # --- A: ISF vs BG at the user's TDD ---
        ax = axes[0]
        bg = np.linspace(70, 300, 200)
        ax.plot(bg, isf_v1(bg, tdd_med), color=C_V1, lw=2, label="v1 (TDD$^{-1}$)")
        ax.plot(bg, isf_v2(bg, tdd_med), color=C_V2, lw=2, label="v2 (TDD$^{-2}$)")
        ax.fill_between(bg, isf_v1(bg, tdd_p75), isf_v1(bg, tdd_p25), color=C_V1, alpha=0.15)
        ax.fill_between(bg, isf_v2(bg, tdd_p75), isf_v2(bg, tdd_p25), color=C_V2, alpha=0.15)
        emp = EMPIRICAL.get(user_id)
        if emp:
            ax.axhline(emp["empirical_isf"], color=C_EMP, ls="--", lw=1.5,
                       label=f"empirical ISF ({emp['empirical_isf']:.0f})")
            ax.axhspan(emp["ci_low_isf"], emp["ci_high_isf"], color=C_EMP, alpha=0.10)
        ax.set_xlabel("glucose (mg/dL)")
        ax.set_ylabel("ISF (mg/dL per U)")
        ax.set_title(f"ISF–glucose curve at TDD {tdd_med:.0f} (band: TDD IQR)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        # --- B: densest 14-day window ---
        ax = axes[1]
        w = _densest_window(df)
        t_days = (w["ts_relative_sec"] - w["ts_relative_sec"].iloc[0]) / 86400.0
        ax2 = ax.twinx()
        ax2.plot(t_days, w["bg"], color=C_BG, lw=0.4, alpha=0.7)
        ax2.set_ylabel("glucose (mg/dL)", color="#888888")
        ax2.tick_params(axis="y", colors="#888888")
        ax.plot(t_days, w["isf_v1"], color=C_V1, lw=0.6, label="ISF V1")
        ax.plot(t_days, w["isf_v2"], color=C_V2, lw=0.6, label="ISF V2")
        ax.set_xlabel("days")
        ax.set_ylabel("ISF (mg/dL per U)")
        ax.set_title("14-day sample (densest data window)")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(alpha=0.3)

        # --- C: per-tick ratio distribution ---
        ax = axes[2]
        ratio = (df["isf_v2"] / df["isf_v1"]).to_numpy()
        ax.hist(ratio, bins=60, color="#9467bd", alpha=0.85)
        ax.axvline(1.0, color="k", lw=1, ls=":")
        ax.axvline(float(np.nanmedian(ratio)), color="#9467bd", lw=1.5, ls="--",
                   label=f"median = {float(np.nanmedian(ratio)):.2f}")
        ax.set_xlabel("ISF V2 / ISF V1 (per tick)")
        ax.set_ylabel("ticks")
        ax.set_title("Correction strength shift (ratio > 1 → V2 doses less)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        fig.tight_layout(rect=(0, 0, 1, 0.93))
        out = USER_DIR / f"{user_id}.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        return {"user": user_id, "status": "ok"}
    except Exception as e:
        return {"user": user_id, "status": "error", "reason": f"{type(e).__name__}: {e}"}


def per_user_summary() -> pd.DataFrame:
    rows = []
    for f in sorted(config.REPLAY_DIR.glob("*.meta.json")):
        m = json.loads(f.read_text())
        if m.get("median_tdd") is None:
            continue
        df, _ = _load(m["user"])
        if len(df) < 500:
            continue
        ratio = float(np.nanmedian(df["isf_v2"] / df["isf_v1"]))
        rows.append({
            "user": m["user"], "platform": m["platform"],
            "n_ticks": int(len(df)),
            "median_tdd": m["median_tdd"],
            "median_isf_v1": m["median_isf_v1"], "median_isf_v2": m["median_isf_v2"],
            "median_ratio": round(ratio, 3),
            "empirical_isf": EMPIRICAL.get(m["user"], {}).get("empirical_isf"),
            "anchor_uncertain": bool(m.get("anchor_uncertain", False)),
            "basal_source": m.get("basal_source", "device"),
        })
    return pd.DataFrame(rows)


def render_cohort_figs(summary: pd.DataFrame) -> None:
    colors = {"v5": "#1f77b4", "v6": "#ff7f0e", "v7": "#2ca02c"}
    labels = {"v5": "Trio (v5)", "v6": "AAPS classic (v6)", "v7": "OpenAPS (v7)"}

    # --- crossover: observed median ratio vs TDD ---
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for plat, g in summary.groupby("platform"):
        ax.scatter(g["median_tdd"], g["median_ratio"], s=22, alpha=0.75,
                   color=colors[plat], label=f"{labels[plat]} (n={len(g)})")
    ax.axhline(1.0, color="k", lw=0.8, alpha=0.5)
    ax.annotate("ratio = 1 (v1 = v2)", (summary["median_tdd"].max(), 1.0),
                fontsize=8, va="bottom", ha="right")
    ax.set_xscale("log")
    ax.set_xlabel("median TDD (U/day, log scale)")
    ax.set_ylabel("median ISF V2 / ISF V1")
    ax.set_title(f"v2 vs v1 correction-strength shift — {len(summary)} users\n"
                 "ratio > 1: V2 estimates weaker corrections; < 1: stronger")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(CHART_DIR / "fig_crossover.png", dpi=150)
    plt.close(fig)

    # --- V1 vs V2 median ISF per user ---
    fig, ax = plt.subplots(figsize=(6.5, 6))
    for plat, g in summary.groupby("platform"):
        ax.scatter(g["median_isf_v1"], g["median_isf_v2"], s=22, alpha=0.75,
                   color=colors[plat], label=labels[plat])
    lim = (0, min(400.0, max(summary["median_isf_v1"].max(),
                             summary["median_isf_v2"].max()) * 1.05))
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("median ISF under V1 (mg/dL per U)")
    ax.set_ylabel("median ISF under V2 (mg/dL per U)")
    ax.set_title("Per-user median dynamic ISF: V1 vs V2")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(CHART_DIR / "fig_isf_scatter.png", dpi=150)
    plt.close(fig)

    # --- the TDD relationship itself: log-log ISF vs TDD, slope -1 vs slope -2 ---
    e_all = summary.dropna(subset=["empirical_isf"])
    if len(e_all) >= 10:
        log_term = np.log(config.NORMAL_TARGET / config.INSULIN_DIVISOR + 1.0)
        fig, ax = plt.subplots(figsize=(8, 6))
        for plat, g in e_all.groupby("platform"):
            ax.scatter(g["median_tdd"], g["empirical_isf"], s=26, alpha=0.8,
                       color=colors[plat], label=f"{labels[plat]} empirical (n={len(g)})")
        tt = np.geomspace(e_all["median_tdd"].min() * 0.8,
                          e_all["median_tdd"].max() * 1.2, 100)
        ax.plot(tt, 1800.0 / (tt * log_term), color=C_V1, lw=2,
                label="V1: 1800/(TDD·logTerm) — slope −1")
        ax.plot(tt, 2300.0 / (log_term * tt ** 2 * 0.02), color=C_V2, lw=2,
                label="V2: 2300/(logTerm·TDD²·0.02) — slope −2")
        # fitted power law of the empirical data
        x = np.log(e_all["median_tdd"].to_numpy(dtype=float))
        y = np.log(e_all["empirical_isf"].to_numpy(dtype=float))
        ok = np.isfinite(x) & np.isfinite(y)
        slope, intercept = np.polyfit(x[ok], y[ok], 1)
        ax.plot(tt, np.exp(intercept) * tt ** slope, "k-", lw=1.6, alpha=0.8,
                label=f"empirical fit: slope {slope:.2f}")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("median TDD (U/day)")
        ax.set_ylabel("ISF at normal target (mg/dL per U)")
        ax.set_title("Which TDD power law does observed sensitivity follow?\n"
                     "V1 assumes ISF ∝ 1/TDD; V2 assumes ISF ∝ 1/TDD²")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3, which="both")
        fig.tight_layout()
        fig.savefig(CHART_DIR / "fig_tdd_loglog.png", dpi=150)
        plt.close(fig)
        print(f"empirical ISF~TDD power-law slope: {slope:.3f} "
              f"(V1 implies -1, V2 implies -2)")

    # --- agreement with empirical ISF ---
    e = summary.dropna(subset=["empirical_isf"])
    if len(e) >= 10:
        fig, axs = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True, sharex=True)
        for ax, col, name, c in ((axs[0], "median_isf_v1", "v1 (TDD$^{-1}$)", C_V1),
                                 (axs[1], "median_isf_v2", "v2 (TDD$^{-2}$)", C_V2)):
            ax.scatter(e["empirical_isf"], e[col], s=24, alpha=0.75, color=c)
            lim = (0, float(np.nanpercentile(
                np.concatenate([e["empirical_isf"], e[col]]), 99)) * 1.1)
            ax.plot(lim, lim, "k--", lw=1)
            ax.set_xlim(lim); ax.set_ylim(lim)
            err = np.abs(np.log(e[col] / e["empirical_isf"]))
            mae = float(np.median(np.abs(e[col] - e["empirical_isf"])))
            ax.set_title(f"{name}\nmedian |error| {mae:.1f} mg/dL/U; "
                         f"median log-error {float(np.median(err)):.2f}")
            ax.set_xlabel("empirical ISF (observed)")
            ax.grid(alpha=0.3)
        axs[0].set_ylabel("formula median ISF")
        fig.suptitle(f"Formula vs empirically observed ISF — {len(e)} users", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        fig.savefig(CHART_DIR / "fig_empirical.png", dpi=150)
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--users", nargs="+", default=None)
    ap.add_argument("--workers", type=int, default=config.DEFAULT_WORKERS)
    ap.add_argument("--cohort-only", action="store_true")
    args = ap.parse_args()

    USER_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if not args.cohort_only:
        users = args.users or sorted(p.stem for p in config.REPLAY_DIR.glob("*.parquet"))
        tasks = [(u,) for u in users]
        print(f"{stamp} rendering {len(tasks)} user pages with {args.workers} workers")
        ctx = mp.get_context("spawn")
        with ctx.Pool(args.workers, maxtasksperchild=config.MAXTASKSPERCHILD) as pool:
            results = list(pool.imap_unordered(render_user_page, tasks))
        ok = sum(1 for r in results if r["status"] == "ok")
        bad = [r for r in results if r["status"] == "error"]
        print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} pages: {ok}/{len(tasks)} ok")
        for r in bad:
            print("  ERROR", r["user"], r["reason"])

    summary = per_user_summary()
    render_cohort_figs(summary)
    summary.to_json(CHART_DIR / "cohort_summary.json", orient="records", indent=1)
    n_weak = int((summary["median_ratio"] > 1).sum())
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} cohort figs done: "
          f"{len(summary)} users, {n_weak} get weaker corrections under V2 "
          f"({100*n_weak/len(summary):.0f}%)")


if __name__ == "__main__":
    main()
