"""2-D local-inertial overland flow model (Bates, Horritt & Fewtrell 2010; de Almeida et al. 2012).

The shallow-water momentum equation without the convective acceleration term, solved on a
regular grid with a staggered (face) discharge q [m^2/s] and cell-centred depth h [m]:

    q_new = (q - g*h_f*dt*dS) / (1 + g*dt*n^2*|q| / h_f^(7/3))
    h_new = h + dt*(sum of face inflows - outflows)/dx

h_f is the flow depth at a face (max water surface minus max bed). The time step follows the
CFL-type condition dt = alpha*dx/sqrt(g*h_max) (alpha <= 0.7 recommended by Bates et al.).
Cells flagged as ``boundary`` (the open sea) have their water surface prescribed by a
stage hydrograph. No rainfall, drainage, infiltration or coastal defences are represented.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

G = 9.81


@dataclass
class Result:
    h: np.ndarray
    hmax: np.ndarray
    t_end: float
    steps: int
    volume_in: float     # net volume that entered from boundary cells (m^3)
    volume_error: float  # stored volume - initial - inflow (m^3), should be ~0


def run(z: np.ndarray, dx: float, n_manning, boundary: np.ndarray, stage, t_end: float,
        alpha: float = 0.7, h0: np.ndarray | None = None, h_dry: float = 1e-3,
        dt_max: float = 10.0, closed: bool = True, theta: float = 0.7) -> Result:
    """Integrate the local-inertial equations.

    z        : bed elevation (m), 2-D
    boundary : bool mask of cells whose water level is prescribed: eta = stage(t)
    stage    : callable t[s] -> water level (m), or None for no boundary forcing
    closed   : outer edge of the grid is a wall (no flux)
    theta    : de Almeida et al. (2012) weighting; theta < 1 blends each face discharge with its
               two neighbours along the flow direction, which damps the checkerboard oscillations
               the original scheme develops at low friction and on steep, rough (DSM) beds
    """
    z = np.asarray(z, dtype=float)
    ny, nx = z.shape
    n2 = np.broadcast_to(np.asarray(n_manning, dtype=float) ** 2, z.shape)
    n2x = np.maximum(n2[:, 1:], n2[:, :-1])
    n2y = np.maximum(n2[1:, :], n2[:-1, :])
    h = np.zeros_like(z) if h0 is None else np.array(h0, dtype=float)
    qx = np.zeros((ny, nx - 1))
    qy = np.zeros((ny - 1, nx))
    zx = np.maximum(z[:, 1:], z[:, :-1])
    zy = np.maximum(z[1:, :], z[:-1, :])
    interior = ~boundary
    v0 = float(h[interior].sum() * dx * dx)
    vin = 0.0
    hmax = h.copy()
    t = 0.0
    steps = 0
    while t < t_end - 1e-9:
        if stage is not None and boundary.any():
            eta_b = stage(t)
            h[boundary] = np.maximum(eta_b - z[boundary], 0.0)
        hm = float(h[interior].max()) if interior.any() else 0.0
        if stage is not None and boundary.any():
            hm = max(hm, float(np.max(np.maximum(stage(t) - z[boundary], 0.0))))
        dt = min(dt_max, alpha * dx / np.sqrt(G * max(hm, h_dry)), t_end - t)
        eta = z + h
        # x faces
        hf = np.maximum(eta[:, 1:], eta[:, :-1]) - zx
        wet = hf > h_dry
        s = (eta[:, 1:] - eta[:, :-1]) / dx
        hf_s = np.where(wet, hf, 1.0)
        qx = _blend(qx, theta, axis=1)
        qx = np.where(wet, (qx - G * hf_s * dt * s) / (1.0 + G * dt * n2x * np.abs(qx) / hf_s ** (7.0 / 3.0)), 0.0)
        # y faces
        hf = np.maximum(eta[1:, :], eta[:-1, :]) - zy
        wet = hf > h_dry
        s = (eta[1:, :] - eta[:-1, :]) / dx
        hf_s = np.where(wet, hf, 1.0)
        qy = _blend(qy, theta, axis=0)
        qy = np.where(wet, (qy - G * hf_s * dt * s) / (1.0 + G * dt * n2y * np.abs(qy) / hf_s ** (7.0 / 3.0)), 0.0)
        # limit outflow so that no cell is drained below zero depth in one step
        out = np.zeros_like(h)
        out[:, :-1] += np.clip(qx, 0, None)
        out[:, 1:] += np.clip(-qx, 0, None)
        out[:-1, :] += np.clip(qy, 0, None)
        out[1:, :] += np.clip(-qy, 0, None)
        avail = h * dx / dt
        f = np.where(out > avail, avail / np.where(out > 0, out, 1.0), 1.0)
        f[boundary] = 1.0
        qx = np.where(qx > 0, qx * f[:, :-1], qx * f[:, 1:])
        qy = np.where(qy > 0, qy * f[:-1, :], qy * f[1:, :])
        dh = np.zeros_like(h)
        dh[:, :-1] -= qx
        dh[:, 1:] += qx
        dh[:-1, :] -= qy
        dh[1:, :] += qy
        dh *= dt / dx
        vin -= float(dh[boundary].sum() * dx * dx)  # what boundary cells lost went into the domain
        h = np.maximum(h + dh, 0.0)
        hmax = np.maximum(hmax, np.where(interior, h, 0.0))
        t += dt
        steps += 1
    v1 = float(h[interior].sum() * dx * dx)
    return Result(h=h, hmax=hmax, t_end=t, steps=steps, volume_in=vin, volume_error=v1 - v0 - vin)


def _blend(q: np.ndarray, theta: float, axis: int) -> np.ndarray:
    """theta*q_i + (1-theta)/2*(q_{i-1} + q_{i+1}) along ``axis`` (end faces unchanged)."""
    if theta >= 1.0 or q.shape[axis] < 3:
        return q
    out = q.copy()
    if axis == 1:
        out[:, 1:-1] = theta * q[:, 1:-1] + 0.5 * (1.0 - theta) * (q[:, :-2] + q[:, 2:])
    else:
        out[1:-1, :] = theta * q[1:-1, :] + 0.5 * (1.0 - theta) * (q[:-2, :] + q[2:, :])
    return out


def coarsen(z: np.ndarray, k: int, how: str = "mean") -> np.ndarray:
    """Aggregate a grid by an integer factor k (trailing rows/cols dropped)."""
    ny, nx = (z.shape[0] // k) * k, (z.shape[1] // k) * k
    b = z[:ny, :nx].reshape(ny // k, k, nx // k, k)
    if how == "mean":
        return b.mean(axis=(1, 3))
    if how == "min":
        return b.min(axis=(1, 3))
    if how == "any":
        return b.any(axis=(1, 3))
    if how == "all":
        return b.all(axis=(1, 3))
    raise ValueError(how)


def hunter_front(x, t, u, n):
    """Analytical depth for a front advancing at constant speed u over a flat frictional plane
    (Hunter et al. 2005, used as Test 1 in Bates et al. 2010):
        h(x,t) = [ (7/3) * (C - n^2 u^2 (x - u t)) ]^(3/7),  with C = 0 here (front at x = u t)."""
    arg = (7.0 / 3.0) * (-(n ** 2) * u ** 2 * (np.asarray(x) - u * t))
    return np.where(arg > 0, arg, 0.0) ** (3.0 / 7.0)
