"""Tide-gauge handling: UHSLC hourly record -> detrended season maxima relative to MSL(1995-2014).

Vertical reference chain (documented in the README):
    raw value      : millimetres above the Tanjong Pagar chart datum (gauge zero)
    MSL(1995-2014) : mean of the annual means 1995-2014 (years with >= MIN_COVERAGE data)
    all extreme-value work below is done in metres above MSL(1995-2014), which is also the
    baseline of the IPCC AR6 projections, so the two can be added without a datum shift.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

MISSING = -32767
MIN_COVERAGE = 0.80          # fraction of hourly values required for a year/season to count
BASELINE = (1995, 2014)      # AR6 reference period
BASELINE_MID = 2005.0        # midpoint used to anchor the linear trend


def read_uhslc_csv(path) -> pd.Series:
    """Hourly sea level in metres above chart datum, UTC index, missing values dropped."""
    d = pd.read_csv(path, header=None, names=["year", "month", "day", "hour", "mm"])
    d = d[d["mm"] != MISSING]
    t = pd.to_datetime(d[["year", "month", "day"]]) + pd.to_timedelta(d["hour"], unit="h")
    return pd.Series(d["mm"].to_numpy() / 1000.0, index=pd.DatetimeIndex(t), name="wl_cd_m")


def hours_in_year(year: int) -> int:
    return 8784 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 8760


def annual_means(s: pd.Series, min_coverage: float = MIN_COVERAGE) -> pd.DataFrame:
    g = s.groupby(s.index.year)
    out = pd.DataFrame({"n": g.size(), "mean": g.mean()})
    out["coverage"] = [n / hours_in_year(y) for y, n in out["n"].items()]
    out["valid"] = out["coverage"] >= min_coverage
    return out


def baseline_msl(am: pd.DataFrame, period=BASELINE) -> float:
    sel = am[(am.index >= period[0]) & (am.index <= period[1]) & am["valid"]]
    if len(sel) == 0:
        raise ValueError("no valid years in baseline period")
    return float(sel["mean"].mean())


@dataclass
class Trend:
    slope_m_per_yr: float
    intercept_at_mid: float   # fitted annual mean at BASELINE_MID (m above chart datum)
    n_years: int


def linear_trend(am: pd.DataFrame) -> Trend:
    v = am[am["valid"]]
    x = v.index.to_numpy(dtype=float) + 0.5 - BASELINE_MID   # annual mean sits mid-year
    slope, icept = np.polyfit(x, v["mean"].to_numpy(), 1)
    return Trend(float(slope), float(icept), int(len(v)))


def season_year(idx: pd.DatetimeIndex) -> np.ndarray:
    """July-June 'season year' labelled by the year in which it starts.

    Singapore's highest water levels occur in the Northeast Monsoon (Nov-Feb); a calendar year
    would split one monsoon season into two blocks and could take two maxima from one season.
    """
    return np.where(idx.month >= 7, idx.year, idx.year - 1)


def detrend(s: pd.Series, msl0: float, trend: Trend | None = None, am: pd.DataFrame | None = None,
            method: str = "linear") -> pd.Series:
    """Water level in metres above MSL(1995-2014) with the long-term MSL change removed.

    method="linear": subtract slope * (t - 2005) (keeps interannual MSL variability, e.g. ENSO).
    method="annual_mean": subtract each calendar year's own mean and add back msl0 (removes
    interannual variability and any undocumented gauge level shifts; a sensitivity case).
    """
    rel = s - msl0
    if method == "linear":
        if trend is None:
            raise ValueError("trend required")
        t = s.index.year + (s.index.dayofyear - 0.5) / 365.25
        return rel - trend.slope_m_per_yr * (t - BASELINE_MID)
    if method == "annual_mean":
        if am is None:
            raise ValueError("annual means required")
        m = am["mean"].reindex(s.index.year).to_numpy()
        return s - m
    raise ValueError(method)


def season_maxima(x: pd.Series, min_coverage: float = MIN_COVERAGE) -> pd.DataFrame:
    sy = season_year(x.index)
    df = pd.DataFrame({"v": x.to_numpy(), "sy": sy}, index=x.index)
    rows = []
    for y, grp in df.groupby("sy"):
        hours = 8784 if hours_in_year(y + 1) == 8784 else 8760  # season spans Feb of y+1
        cov = len(grp) / hours
        i = grp["v"].idxmax()
        rows.append({"season": int(y), "coverage": cov, "max": float(grp["v"].max()), "time": i})
    out = pd.DataFrame(rows).set_index("season")
    out["valid"] = out["coverage"] >= min_coverage
    return out
