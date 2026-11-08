import numpy as np
import pandas as pd
import pytest
from scipy import stats

from sgci import eva, tide


def _hourly(start, end, fn):
    idx = pd.date_range(start, end, freq="h")
    t = idx.year + (idx.dayofyear - 0.5) / 365.25
    return pd.Series(fn(t, idx), index=idx)


def test_read_uhslc_drops_missing_keeps_negative(tmp_path):
    p = tmp_path / "h.csv"
    p.write_text("2000,1,1,0,1500\n2000,1,1,1,-32767\n2000,1,1,2,-30\n")
    s = tide.read_uhslc_csv(p)
    assert len(s) == 2
    assert s.iloc[0] == pytest.approx(1.5)
    assert s.iloc[1] == pytest.approx(-0.03)
    assert s.index[1] == pd.Timestamp("2000-01-01 02:00")


def test_hours_in_year():
    assert tide.hours_in_year(2000) == 8784
    assert tide.hours_in_year(1900) == 8760
    assert tide.hours_in_year(2023) == 8760
    assert tide.hours_in_year(2024) == 8784


def test_season_year_july_to_june():
    idx = pd.DatetimeIndex(["2000-06-30 23:00", "2000-07-01 00:00", "2001-01-15"])
    assert list(tide.season_year(idx)) == [1999, 2000, 2000]


def test_annual_means_and_baseline_skip_invalid_years():
    s = _hourly("1994-01-01", "2015-12-31 23:00", lambda t, i: np.where(i.year == 2000, 9.0, 1.0))
    s = s[~((s.index.year == 2000) & (s.index.month > 2))]   # 2000 only Jan-Feb -> invalid
    am = tide.annual_means(s)
    assert not am.loc[2000, "valid"]
    assert am.loc[2001, "valid"]
    assert tide.baseline_msl(am) == pytest.approx(1.0)


def test_linear_trend_recovers_slope_and_anchor():
    s = _hourly("1990-01-01", "2019-12-31 23:00", lambda t, i: 1.7 + 0.003 * (t - 2005))
    tr = tide.linear_trend(tide.annual_means(s))
    assert tr.slope_m_per_yr == pytest.approx(0.003, rel=1e-3)
    assert tr.intercept_at_mid == pytest.approx(1.7, abs=1e-3)


def test_detrend_linear_removes_trend_relative_to_baseline():
    s = _hourly("1990-01-01", "2019-12-31 23:00",
                lambda t, i: 1.7 + 0.003 * (t - 2005) + 0.5 * np.sin(2 * np.pi * np.arange(len(t)) / 12.42))
    am = tide.annual_means(s)
    msl0 = tide.baseline_msl(am)
    x = tide.detrend(s, msl0, tide.linear_trend(am), am, method="linear")
    m = x.groupby(x.index.year).mean()
    assert abs(np.polyfit(m.index, m.values, 1)[0]) < 1e-5
    assert abs(m.mean()) < 0.01


def test_detrend_annual_mean_zero_mean_each_year():
    s = _hourly("2001-01-01", "2004-12-31 23:00", lambda t, i: 1.0 + 0.2 * (i.year - 2001))
    am = tide.annual_means(s)
    x = tide.detrend(s, 1.0, None, am, method="annual_mean")
    assert np.allclose(x.groupby(x.index.year).mean(), 0.0)


def test_season_maxima_picks_max_and_flags_coverage():
    s = _hourly("2000-07-01", "2002-06-30 23:00", lambda t, i: np.zeros(len(t)))
    s.loc["2001-01-10 05:00"] = 2.0
    s.loc["2001-12-01 03:00"] = 1.5
    s = s.drop(s.loc["2001-07-01":"2002-03-01"].index[::2])   # thin season 2001 below 80%
    sm = tide.season_maxima(s)
    assert sm.loc[2000, "max"] == 2.0
    assert sm.loc[2000, "time"] == pd.Timestamp("2001-01-10 05:00")
    assert sm.loc[2000, "valid"]
    assert not sm.loc[2001, "valid"]


# ------------------------------------------------------------------ GEV

@pytest.mark.parametrize("xi", [-0.2, 0.0, 0.15])
def test_return_level_matches_scipy(xi):
    g = eva.GEV(mu=1.5, sigma=0.07, xi=xi, nll=0, n=0)
    T = np.array([2, 10, 100, 1000])
    expect = stats.genextreme.ppf(1 - 1 / T, -xi, loc=1.5, scale=0.07)
    assert np.allclose(eva.gev_return_level(g, T), expect, rtol=1e-9, atol=1e-9)


@pytest.mark.parametrize("xi", [-0.15, 0.2])
def test_gev_nll_matches_scipy_logpdf(xi):
    x = stats.genextreme.rvs(-xi, loc=1.6, scale=0.08, size=50, random_state=1)
    nll = eva.gev_nll((1.6, np.log(0.08), xi), x)
    assert nll == pytest.approx(-stats.genextreme.logpdf(x, -xi, loc=1.6, scale=0.08).sum(), rel=1e-9)


def test_gev_nll_gumbel_limit_continuous():
    x = np.linspace(1.4, 1.9, 30)
    a = eva.gev_nll((1.6, np.log(0.07), 0.0), x)
    b = eva.gev_nll((1.6, np.log(0.07), 1e-4), x)
    assert a == pytest.approx(b, rel=1e-3)


def test_fit_gev_recovers_parameters():
    x = stats.genextreme.rvs(-0.1, loc=1.6, scale=0.07, size=4000, random_state=7)
    g = eva.fit_gev(x)
    assert g.mu == pytest.approx(1.6, abs=0.01)
    assert g.sigma == pytest.approx(0.07, rel=0.05)
    assert g.xi == pytest.approx(0.1, abs=0.04)


def test_fit_gev_agrees_with_scipy_small_sample():
    x = stats.genextreme.rvs(0.1, loc=1.6, scale=0.06, size=40, random_state=3)
    g = eva.fit_gev(x)
    c, loc, sc = stats.genextreme.fit(x)
    assert g.xi == pytest.approx(-c, abs=0.01)
    assert g.mu == pytest.approx(loc, abs=0.002)
    assert g.nll <= eva.gev_nll((loc, np.log(sc), -c), x) + 1e-6


def test_fixed_xi_fit_is_gumbel():
    x = stats.gumbel_r.rvs(loc=1.5, scale=0.05, size=3000, random_state=2)
    g = eva.fit_gev(x, fix_xi=0.0)
    assert g.xi == 0.0
    assert g.mu == pytest.approx(1.5, abs=0.005)


def test_bootstrap_ci_brackets_and_is_reproducible():
    x = stats.genextreme.rvs(0.0, loc=1.6, scale=0.07, size=37, random_state=5)
    T = [10, 100]
    a = eva.bootstrap_gev(x, T, n_boot=40, seed=1)
    b = eva.bootstrap_gev(x, T, n_boot=40, seed=1)
    rl = eva.gev_return_level(eva.fit_gev(x), T)
    assert np.allclose(a["rl_lo"], b["rl_lo"])
    assert np.all(a["rl_lo"] < rl) and np.all(rl < a["rl_hi"])
    assert a["rl_hi"][1] - a["rl_lo"][1] > a["rl_hi"][0] - a["rl_lo"][0]   # wider for rarer events


def test_lr_test_detects_heavy_tail():
    heavy = stats.genextreme.rvs(-0.4, loc=0, scale=1, size=800, random_state=4)
    gum = stats.gumbel_r.rvs(size=800, random_state=4)
    assert eva.lr_test_gumbel(heavy)["p_value"] < 0.01
    assert eva.lr_test_gumbel(gum)["p_value"] > 0.01


# ------------------------------------------------------------------ POT / GPD

def test_decluster_runs_separates_clusters():
    v = np.array([0, 2, 3, 0, 0, 0, 5, 0, 0, 0, 0, 0, 4], float)
    t = np.arange(len(v)) * 24.0
    peaks = eva.decluster_runs(v, t, 1.0, run_h=72.0)
    assert list(peaks) == [3.0, 5.0, 4.0]
    peaks = eva.decluster_runs(v, t, 1.0, run_h=100.0)
    assert list(peaks) == [5.0, 4.0]


def test_gpd_nll_matches_scipy():
    y = stats.genpareto.rvs(0.1, scale=0.05, size=60, random_state=1)
    assert eva.gpd_nll((np.log(0.05), 0.1), y) == pytest.approx(
        -stats.genpareto.logpdf(y, 0.1, scale=0.05).sum(), rel=1e-9)


def test_gpd_return_level_formula():
    g = eva.GPD(u=1.6, sigma=0.05, xi=-0.1, rate=3.0, n=0)
    T = np.array([10.0, 100.0])
    expect = 1.6 + stats.genpareto.ppf(1 - 1 / (3.0 * T), -0.1, scale=0.05)
    assert np.allclose(eva.gpd_return_level(g, T), expect)


def test_fit_gpd_recovers():
    y = stats.genpareto.rvs(-0.1, scale=0.05, size=4000, random_state=9)
    g = eva.fit_gpd(y + 1.0, 1.0, years=1000.0)
    assert g.xi == pytest.approx(-0.1, abs=0.04)
    assert g.sigma == pytest.approx(0.05, rel=0.06)
    assert g.rate == pytest.approx(4.0)
