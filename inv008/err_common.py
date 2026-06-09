#!/usr/bin/env python3
"""Shared definitions for the static-ISF error-curve analysis (err_curve / err_consistency).

The same-window dataset (head_to_head_windows.parquet) gives, per overnight window, the
glucose error of the person's *static profile* ISF as a drop predictor:

    err_static = actual_end_BG − predicted_end_BG(profile ISF)      [mg/dL]
                 > 0  ⇒ the profile ISF OVER-predicted the drop (real drop was smaller)

Two views of the same quantity, and the gap between them is the point of the analysis:

  absolute  (mg/dL)        err_static — what the controller actually feels at the wheel.
  fractional (unit-free)   frac_overpred = 1 − realised_isf / profile_isf
                           = err_static / predicted_drop(profile)
                           — the effective-sensitivity gap. Removes the mechanical fact that
                             corrections (and so absolute errors) are simply bigger at high BG.

realised_isf is a noise-amplifying ratio (sug_isf · drop / predicted_drop); its tails are
meaningless, so the fractional metric is winsorised. err_static is glucose-scale and bounded,
used as-is with medians.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inv008 import config

WINDOWS = config.ROOT / "results" / "head_to_head_windows.parquet"

# reporting bands across the actionable overnight range
BG_BANDS = [(80, 100), (100, 120), (120, 145), (145, 175), (175, 205), (205, 260)]
BAND_CTR = [(a + b) / 2 for a, b in BG_BANDS]
BAND_LBL = [f"{a}-{b}" for a, b in BG_BANDS]

FRAC_CLIP = 3.0          # winsorise fractional over-prediction to ±300%
MIN_USER_WINDOWS = 60    # a user must clear this to enter per-user statistics
MIN_BAND_WINDOWS = 8     # min windows for a per-user-per-band median


def load_windows() -> pd.DataFrame:
    """Load the window dataset with the derived error metrics attached."""
    d = pd.read_parquet(WINDOWS)
    for c in ["bg", "tdd", "iob", "hour", "profile_isf", "realised_isf",
              "err_static", "start_slope", "bg_end"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d.profile_isf > 0) & d.err_static.notna() & d.bg.notna()].copy()
    d["act_drop"] = d.bg - d.bg_end
    d["eff_ratio"] = d.realised_isf / d.profile_isf          # <1 ⇒ less sensitive than profile
    d["frac_overpred"] = (1.0 - d.eff_ratio).clip(-FRAC_CLIP, FRAC_CLIP)
    d["late_hour"] = d.hour.isin([1, 2]).astype(float)       # 1-2am vs 23-0
    # restrict to users with enough windows for stable per-user curves
    keep = d.groupby("user").bg.transform("size") >= MIN_USER_WINDOWS
    return d[keep].reset_index(drop=True)


def band_of(bg: np.ndarray) -> np.ndarray:
    """Index of the BG band each reading falls in, or -1 if outside all bands."""
    out = np.full(len(bg), -1)
    for k, (a, b) in enumerate(BG_BANDS):
        out[(bg >= a) & (bg < b)] = k
    return out


def boot_median_ci(x: np.ndarray, n_boot: int = 2000, seed: int = 0, lo=5, hi=95):
    """Percentile bootstrap CI for a median (resampling the values)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return (np.nan, np.nan, np.nan)
    rng = np.random.default_rng(seed)
    meds = np.median(x[rng.integers(0, len(x), size=(n_boot, len(x)))], axis=1)
    return (float(np.median(x)), float(np.percentile(meds, lo)), float(np.percentile(meds, hi)))


def ols(y: np.ndarray, X: np.ndarray):
    """Plain OLS. X already includes an intercept column. Returns (beta, se, n)."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    y, X = y[ok], X[ok]
    n, p = X.shape
    if n <= p:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(n - p, 1)
    sigma2 = float(resid @ resid) / dof
    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return None
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    return beta, se, n


def dersimonian_laird(betas: np.ndarray, ses: np.ndarray):
    """Random-effects meta-analysis of per-user coefficients (DerSimonian-Laird).

    Pools many noisy per-user estimates into a population effect while estimating the genuine
    between-user spread τ. Answers 'is this effect real *and* shared across people?'.
    """
    from scipy import stats
    b = np.asarray(betas, float)
    s = np.asarray(ses, float)
    ok = np.isfinite(b) & np.isfinite(s) & (s > 0)
    b, s = b[ok], s[ok]
    k = len(b)
    if k < 3:
        return None
    w = 1.0 / s**2
    b_fixed = float(np.sum(w * b) / np.sum(w))
    Q = float(np.sum(w * (b - b_fixed) ** 2))
    C = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1.0 / (s**2 + tau2)
    b_re = float(np.sum(wr * b) / np.sum(wr))
    se_re = float(np.sqrt(1.0 / np.sum(wr)))
    z = b_re / se_re if se_re > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return {"k": k, "b_re": b_re, "se_re": se_re, "z": float(z), "p": p,
            "tau": float(np.sqrt(tau2)), "I2_pct": float(I2),
            "frac_same_sign": float(np.mean(np.sign(b) == np.sign(b_re)))}
