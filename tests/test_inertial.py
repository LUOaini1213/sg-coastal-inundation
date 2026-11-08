import numpy as np
import pytest

from sgci import inertial as M


def test_hunter_analytical_front():
    """Bates et al. (2010) test 1: front over a flat plane, n = 0.03, u = 1 m/s, t = 1 h."""
    u, n, dx = 1.0, 0.03, 10.0
    nx = int(5000 / dx)
    z = np.zeros((3, nx))
    b = np.zeros_like(z, bool)
    b[:, 0] = True
    r = M.run(z, dx, n, b, lambda t: float(M.hunter_front(-dx / 2, t, u, n)), 3600.0)
    x = np.arange(nx) * dx - dx / 2
    ha = M.hunter_front(x, 3600.0, u, n)
    m = ha > 0
    rmse = np.sqrt(np.mean((r.h[1][m] - ha[m]) ** 2))
    assert rmse / ha.max() < 0.05
    front = x[r.h[1] > 0.01].max()
    assert abs(front - 3600.0) < 200.0


def test_lake_at_rest_is_preserved():
    rng = np.random.default_rng(0)
    z = rng.uniform(-1.0, 0.5, size=(20, 20))
    h0 = 1.0 - z
    r = M.run(z, 30.0, 0.035, np.zeros_like(z, bool), None, 3600.0, h0=h0)
    assert np.abs(z + r.h - 1.0).max() < 1e-9


def test_closed_basin_conserves_mass():
    z = np.zeros((30, 30))
    z[10:20, 10:20] = -0.5
    h0 = np.where(z < 0, 1.0, 0.2)
    h0[5, 5] = 2.0
    r = M.run(z, 20.0, 0.03, np.zeros_like(z, bool), None, 1800.0, h0=h0)
    assert r.h.sum() == pytest.approx(h0.sum(), rel=1e-10)
    assert abs(r.volume_error) < 1e-6
    assert r.h.min() >= 0.0


def test_steady_stage_fills_low_land_to_boundary_level():
    z = np.full((5, 40), 0.3)
    z[:, 0] = -2.0
    b = np.zeros_like(z, bool)
    b[:, 0] = True
    r = M.run(z, 30.0, 0.03, b, lambda t: 1.0, 6 * 3600.0)
    assert np.allclose((z + r.h)[:, 1:], 1.0, atol=0.02)
    assert r.volume_in == pytest.approx(r.h[:, 1:].sum() * 900.0, rel=1e-6)


def test_ridge_above_stage_stays_dry_behind():
    z = np.full((5, 30), 0.2)
    z[:, 0] = -2.0
    z[:, 10] = 1.5
    b = np.zeros_like(z, bool)
    b[:, 0] = True
    r = M.run(z, 30.0, 0.03, b, lambda t: 1.0, 3 * 3600.0)
    assert r.hmax[:, 11:].max() == 0.0
    assert r.hmax[:, 5].min() > 0.5


def test_short_pulse_floods_less_than_steady():
    z = np.zeros((5, 120))
    z[:, 0] = -2.0
    b = np.zeros_like(z, bool)
    b[:, 0] = True
    pulse = lambda t: 0.5 + 0.5 * np.exp(-((t - 3600.0) / 1200.0) ** 2)
    rp = M.run(z, 30.0, 0.05, b, pulse, 4 * 3600.0)
    rs = M.run(z, 30.0, 0.05, b, lambda t: 1.0, 4 * 3600.0)
    assert (rp.hmax > 0.05).sum() < (rs.hmax > 0.05).sum()


def test_coarsen():
    a = np.arange(16, dtype=float).reshape(4, 4)
    assert M.coarsen(a, 2, "mean")[0, 0] == pytest.approx(2.5)
    assert M.coarsen(a, 2, "min")[1, 1] == 10
    assert M.coarsen(a > 4, 2, "any")[0, 0] and not M.coarsen(a > 4, 2, "all")[0, 0]


def test_first_step_flux_uses_higher_bed_for_flow_depth():
    """Water on a raised cell spilling onto a lower one: h_f = max(eta) - max(z) = 0.5 m."""
    dx, dt = 10.0, 0.05
    z = np.array([[1.0, 0.0]])
    h0 = np.array([[0.5, 0.1]])
    r = M.run(z, dx, 0.03, np.zeros_like(z, bool), None, dt, h0=h0, dt_max=dt)
    q = M.G * 0.5 * dt * (1.5 - 0.1) / dx            # friction term vanishes because q0 = 0
    assert h0[0, 0] - r.h[0, 0] == pytest.approx(q * dt / dx, rel=1e-6)


def _rough():
    rng = np.random.default_rng(1)
    z = rng.uniform(0, 1.5, (60, 60))
    z[rng.random((60, 60)) < 0.15] += rng.uniform(3, 15)      # building-like blocks in a DSM
    z[:, :3] = -1.0
    b = np.zeros_like(z, bool)
    b[:, :3] = True
    return z, b


def test_ramped_cold_start_reproduces_bathtub_without_overshoot():
    from sgci import inundation
    z, b = _rough()
    fill = inundation.fill_level(z, b)
    h0 = np.where((fill < 0.5) & ~b, np.maximum(0.5 - z, 0), 0)
    r = M.run(z, 30.0, 0.035, b, lambda t: 0.5 + 1.5 * min(1, t / 3600), 4 * 3600.0, h0=h0)
    eta = np.where(r.hmax > 0.01, z + r.hmax, -9)
    assert eta[:, 3:].max() < 2.0 + 0.02
    bath = (fill < 2.0) & (2.0 - z >= 0.05) & ~b
    assert ((r.hmax > 0.05) & ~b).sum() == bath.sum()


def test_theta_blending_damps_dam_break_overshoot():
    z, b = _rough()
    plain = M.run(z, 30.0, 0.035, b, lambda t: 2.0, 2 * 3600.0, theta=1.0)
    damped = M.run(z, 30.0, 0.035, b, lambda t: 2.0, 2 * 3600.0)
    top = lambda r: np.where(r.hmax > 0.01, z + r.hmax, -9)[:, 3:].max()
    assert top(damped) < top(plain) - 0.3
