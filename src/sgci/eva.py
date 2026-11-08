"""Extreme value analysis: GEV on block maxima, GPD on declustered peaks-over-threshold.

Shape convention: xi > 0 heavy (Frechet) tail, xi < 0 bounded (Weibull) tail.
(scipy.stats.genextreme uses c = -xi; we keep our own likelihood to make the sign explicit.)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

EULER = 0.5772156649015329
XI_EPS = 1e-6


def gev_nll(params, x) -> float:
    mu, log_sigma, xi = params
    sigma = np.exp(log_sigma)
    z = (np.asarray(x) - mu) / sigma
    if abs(xi) < XI_EPS:
        return float(len(z) * log_sigma + np.sum(z + np.exp(-z)))
    t = 1.0 + xi * z
    if np.any(t <= 0):
        return 1e10
    return float(len(z) * log_sigma + (1.0 + 1.0 / xi) * np.sum(np.log(t)) + np.sum(t ** (-1.0 / xi)))


@dataclass
class GEV:
    mu: float
    sigma: float
    xi: float
    nll: float
    n: int


def fit_gev(x, fix_xi: float | None = None, start: GEV | None = None) -> GEV:
    x = np.asarray(x, dtype=float)
    s0 = np.std(x) * np.sqrt(6) / np.pi            # Gumbel method-of-moments start
    m0 = np.mean(x) - EULER * s0
    if fix_xi is not None:
        f = lambda p: gev_nll((p[0], p[1], fix_xi), x)
        r = optimize.minimize(f, [m0, np.log(s0)], method="Nelder-Mead",
                              options={"xatol": 1e-9, "fatol": 1e-11, "maxiter": 20000})
        return GEV(float(r.x[0]), float(np.exp(r.x[1])), float(fix_xi), float(r.fun), len(x))
    best = None
    starts = [[m0, np.log(s0), xi0] for xi0 in (-0.2, -0.05, 0.05, 0.2)]
    if start is not None:  # warm start (bootstrap): previous optimum plus one generic start
        starts = [[start.mu, np.log(start.sigma), start.xi], starts[2]]
    for p0 in starts:
        r = optimize.minimize(lambda p: gev_nll(p, x), p0, method="Nelder-Mead",
                              options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 20000})
        if best is None or r.fun < best.fun:
            best = r
    mu, ls, xi = best.x
    return GEV(float(mu), float(np.exp(ls)), float(xi), float(best.fun), len(x))


def gev_return_level(g: GEV, T) -> np.ndarray:
    """Level exceeded on average once every T blocks (annual exceedance probability 1/T)."""
    T = np.asarray(T, dtype=float)
    y = -np.log(1.0 - 1.0 / T)
    if abs(g.xi) < XI_EPS:
        return g.mu - g.sigma * np.log(y)
    return g.mu + g.sigma / g.xi * (y ** (-g.xi) - 1.0)


def gev_cdf(g: GEV, z) -> np.ndarray:
    return stats.genextreme.cdf(z, -g.xi, loc=g.mu, scale=g.sigma)


def bootstrap_gev(x, T, n_boot: int = 2000, seed: int = 20261106, ci=(2.5, 97.5)):
    """Non-parametric bootstrap of block maxima -> percentile CI of return levels and xi."""
    x = np.asarray(x, dtype=float)
    rng = np.random.default_rng(seed)
    rl = np.empty((n_boot, len(T)))
    xis = np.empty(n_boot)
    g0 = fit_gev(x)
    for b in range(n_boot):
        g = fit_gev(rng.choice(x, size=len(x), replace=True), start=g0)
        rl[b] = gev_return_level(g, T)
        xis[b] = g.xi
    return {"rl_lo": np.percentile(rl, ci[0], axis=0), "rl_hi": np.percentile(rl, ci[1], axis=0),
            "xi_lo": float(np.percentile(xis, ci[0])), "xi_hi": float(np.percentile(xis, ci[1])),
            "samples": rl}


def lr_test_gumbel(x) -> dict:
    """Likelihood-ratio test of H0: xi = 0 (Gumbel) against the full GEV."""
    full = fit_gev(x)
    gum = fit_gev(x, fix_xi=0.0)
    d = 2.0 * (gum.nll - full.nll)
    return {"deviance": float(d), "p_value": float(stats.chi2.sf(max(d, 0.0), df=1))}


# --------------------------------------------------------------------------- POT / GPD

def decluster_runs(values: np.ndarray, times_h: np.ndarray, u: float, run_h: float = 72.0) -> np.ndarray:
    """Cluster maxima of exceedances of u; a cluster ends after > run_h hours without exceedance."""
    idx = np.flatnonzero(values > u)
    if idx.size == 0:
        return np.array([])
    peaks = []
    cur = [idx[0]]
    for i in idx[1:]:
        if times_h[i] - times_h[cur[-1]] > run_h:
            peaks.append(values[cur].max())
            cur = [i]
        else:
            cur.append(i)
    peaks.append(values[cur].max())
    return np.asarray(peaks)


def gpd_nll(params, y) -> float:
    log_sigma, xi = params
    sigma = np.exp(log_sigma)
    y = np.asarray(y)
    if abs(xi) < XI_EPS:
        return float(len(y) * log_sigma + np.sum(y) / sigma)
    t = 1.0 + xi * y / sigma
    if np.any(t <= 0):
        return 1e10
    return float(len(y) * log_sigma + (1.0 + 1.0 / xi) * np.sum(np.log(t)))


@dataclass
class GPD:
    u: float
    sigma: float
    xi: float
    rate: float   # cluster exceedances per year
    n: int


def fit_gpd(peaks, u: float, years: float) -> GPD:
    y = np.asarray(peaks) - u
    best = None
    for xi0 in (-0.2, 0.0, 0.2):
        r = optimize.minimize(lambda p: gpd_nll(p, y), [np.log(np.mean(y)), xi0], method="Nelder-Mead",
                              options={"xatol": 1e-9, "fatol": 1e-11, "maxiter": 20000})
        if best is None or r.fun < best.fun:
            best = r
    return GPD(u, float(np.exp(best.x[0])), float(best.x[1]), len(y) / years, len(y))


def gpd_return_level(g: GPD, T) -> np.ndarray:
    T = np.asarray(T, dtype=float)
    m = g.rate * T
    if abs(g.xi) < XI_EPS:
        return g.u + g.sigma * np.log(m)
    return g.u + g.sigma / g.xi * (m ** g.xi - 1.0)
