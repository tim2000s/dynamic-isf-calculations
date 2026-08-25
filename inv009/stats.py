"""Shared statistics. Small, and deliberately explicit about what each one assumes.

statsmodels is not installed on this machine, and a mixed model is not needed
here anyway: every question is either a slope across subjects or a slope within
one, and the second is handled the way INV-008 handled it, by fitting each
subject separately and pooling with DerSimonian and Laird.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inv008.err_common import dersimonian_laird, ols  # noqa: F401  (re-exported)


def loglog(x, y, w=None):
    """Slope of log y on log x, with a standard error and an R squared.

    Weights are optional and are the inverse variance of each point where the
    y values are themselves estimates with known precision.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if w is not None:
        w = np.asarray(w, dtype=float)
        m &= np.isfinite(w) & (w > 0)
    lx, ly = np.log(x[m]), np.log(y[m])
    n = len(lx)
    if n < 5:
        return dict(slope=np.nan, se=np.nan, intercept=np.nan, n=n, r2=np.nan)
    ww = np.ones(n) if w is None else w[m]
    ww = ww / ww.mean()
    X = np.column_stack([np.ones(n), lx])
    W = np.diag(ww)
    xtwx_inv = np.linalg.pinv(X.T @ W @ X)
    beta = xtwx_inv @ (X.T @ W @ ly)
    resid = ly - X @ beta
    dof = max(n - 2, 1)
    s2 = float((ww * resid ** 2).sum() / dof)
    cov = s2 * xtwx_inv
    ss_tot = float((ww * (ly - np.average(ly, weights=ww)) ** 2).sum())
    return dict(slope=float(beta[1]), se=float(np.sqrt(max(cov[1, 1], 0))),
                intercept=float(beta[0]), n=int(n),
                r2=float(1 - (ww * resid ** 2).sum() / ss_tot) if ss_tot > 0 else np.nan)


def deming(x, y, lam: float = 1.0):
    """Slope when BOTH variables carry measurement error.

    Ordinary least squares assumes the x axis is known exactly. Daily dose is
    not: it is an average over a finite number of days, from a pump record with
    gaps. Error in x drags an ordinary slope toward zero, which here would make
    every candidate exponent look shallower than it is, so the errors-in-
    variables slope is reported beside it. `lam` is the ratio of error variances,
    and one means the two are equally noisy in logs.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    lx, ly = np.log(x[m]), np.log(y[m])
    if len(lx) < 5:
        return np.nan
    sxx = float(np.var(lx, ddof=1))
    syy = float(np.var(ly, ddof=1))
    sxy = float(np.cov(lx, ly, ddof=1)[0, 1])
    if abs(sxy) < 1e-12:
        return np.nan
    return float((syy - lam * sxx + np.sqrt((syy - lam * sxx) ** 2 + 4 * lam * sxy ** 2))
                 / (2 * sxy))


def boot_slope_ci(x, y, w=None, n_boot: int = 2000, seed: int = 0, lo=2.5, hi=97.5):
    """Percentile interval for a log-log slope, resampling subjects.

    Takes the same weights as the point estimate. An unweighted interval around
    a weighted slope is not an interval for that slope, and the two can disagree
    enough to put the estimate outside its own confidence bounds.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if w is not None:
        w = np.asarray(w, dtype=float)
        m &= np.isfinite(w) & (w > 0)
    lx, ly = np.log(x[m]), np.log(y[m])
    ww = np.ones(len(lx)) if w is None else w[m]
    n = len(lx)
    if n < 5:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    out = np.full(n_boot, np.nan)
    for i in range(n_boot):
        j = rng.integers(0, n, n)
        sx, sy, sw = lx[j], ly[j], ww[j]
        sw = sw / sw.mean()
        mx = np.average(sx, weights=sw)
        v = float(np.sum(sw * (sx - mx) ** 2))
        if v > 1e-12:
            out[i] = float(np.sum(sw * (sx - mx) * (sy - np.average(sy, weights=sw))) / v)
    return tuple(float(v) for v in np.nanpercentile(out, [lo, hi]))


def band_of(v, bands):
    """Label a value by the band it falls in."""
    v = np.asarray(v, dtype=float)
    out = np.full(len(v), "", dtype=object)
    for b in bands:
        lo, hi = b[0], b[1]
        label = b[2] if len(b) > 2 else f"{lo}-{hi}"
        out[(v >= lo) & (v < hi)] = label
    return out


def describe_slope(res: dict, ci: tuple) -> str:
    """A slope with its interval, and which candidate laws it is consistent with."""
    s, (lo, hi) = res["slope"], ci
    laws = []
    for name, val in (("v1 (1/TDD)", -1.0), ("v2 (1/TDD^2)", -2.0),
                      ("root TDD", -0.5), ("no relationship", 0.0)):
        if lo <= val <= hi:
            laws.append(name)
    return (f"{s:+.3f} (95% CI {lo:+.3f} to {hi:+.3f}, n={res['n']}); "
            f"consistent with: {', '.join(laws) if laws else 'none of the candidates'}")
