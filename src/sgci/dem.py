"""Copernicus GLO-30 DEM mosaic for Singapore and the vertical datum shift to local MSL.

Copernicus DEM heights are orthometric heights above the EGM2008 geoid. Local mean sea level
is not the geoid: it sits above it by the mean dynamic topography (MDT). We convert

    z_msl = z_egm2008 - MDT

with the MDT taken from DTU22MDT (Knudsen et al.) averaged over the ocean cells around
Singapore. This shift is the single most influential assumption in the study, so results are
also reported for MDT +/- 0.3 m and for MDT = 0 (no shift, a common but incorrect default).
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

BBOX = (103.59, 1.15, 104.10, 1.48)  # lon_min, lat_min, lon_max, lat_max
MDT_BOX = (103.55, 1.05, 104.25, 1.45)  # ocean cells used to average DTU22MDT
PIXEL_DEG = 1.0 / 3600.0


def load_dem(tile_paths, bbox=BBOX):
    import rasterio
    from rasterio.merge import merge

    srcs = [rasterio.open(p) for p in tile_paths]
    try:
        arr, transform = merge(srcs, bounds=bbox)
    finally:
        for s in srcs:
            s.close()
    return arr[0].astype(np.float64), transform


def cell_area_m2(transform, rows: int) -> np.ndarray:
    """Area of each row's cells (m^2) on the WGS84 sphere approximation (varies <0.1% here)."""
    R = 6371008.8
    lat = transform.f + (np.arange(rows) + 0.5) * transform.e
    dy = abs(transform.e) * np.pi / 180 * R
    dx = transform.a * np.pi / 180 * R * np.cos(np.deg2rad(lat))
    return dx * dy


def read_dtu_mdt(path: Path):
    """Parse a gzipped GRAVSOFT ASCII grid -> (lats N->S, lons 0..360, grid)."""
    tok = gzip.open(path, "rt").read().split()
    la1, la2, lo1, lo2, dla, dlo = map(float, tok[:6])
    nla = int(round((la2 - la1) / dla)) + 1
    nlo = int(round((lo2 - lo1) / dlo)) + 1
    g = np.asarray(tok[6:], dtype=float).reshape(nla, nlo)
    lats = la2 - np.arange(nla) * dla
    lons = lo1 + np.arange(nlo) * dlo
    return lats, lons, g


def mdt_offset(path: Path, box=MDT_BOX) -> dict:
    lats, lons, g = read_dtu_mdt(path)
    i = (lats >= box[1]) & (lats <= box[3])
    j = (lons >= box[0]) & (lons <= box[2])
    sub = g[np.ix_(i, j)]
    sub = sub[np.isfinite(sub) & (np.abs(sub) < 9)]
    return {"mean": float(sub.mean()), "min": float(sub.min()), "max": float(sub.max()), "n": int(sub.size)}
