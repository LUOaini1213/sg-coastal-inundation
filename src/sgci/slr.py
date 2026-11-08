"""IPCC AR6 relative sea-level projections (Fox-Kemper et al. 2021; Garner et al. 2021 data).

Values are total relative sea-level change (including vertical land motion) in metres
relative to the 1995-2014 mean, medium-confidence workflow, at the PSMSL location of the
Tanjong Pagar tide gauge (PSMSL id 1746, 1.27N 103.85E).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

TANJONG_PAGAR_PSMSL = 1746
QUANTILES = {"p05": 0.05, "p17": 0.17, "p50": 0.50, "p83": 0.83, "p95": 0.95}


def read_location(nc_path: Path, location_id: int = TANJONG_PAGAR_PSMSL) -> dict:
    import h5py

    with h5py.File(nc_path, "r") as f:
        locs = f["locations"][:]
        j = int(np.flatnonzero(locs == location_id)[0])
        q = f["quantiles"][:]
        years = f["years"][:].astype(int)
        slc = f["sea_level_change"][:, :, j].astype(float)  # (quantile, year) in mm
        lat, lon = float(f["lat"][j]), float(f["lon"][j])
    out = {"lat": lat, "lon": lon, "years": years.tolist(), "quantiles": {}}
    for name, qq in QUANTILES.items():
        i = int(np.argmin(np.abs(q - qq)))
        if abs(q[i] - qq) > 1e-6:
            raise ValueError(f"quantile {qq} not in file")
        out["quantiles"][name] = (slc[i] / 1000.0).round(4).tolist()
    return out


def value(proj: dict, year: int, q: str = "p50") -> float:
    return float(proj["quantiles"][q][proj["years"].index(year)])
