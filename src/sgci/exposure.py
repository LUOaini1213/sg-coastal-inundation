"""Overlay flood extents with buildings, roads, MRT exits and Census 2020 population.

DSM caveat handled here: Copernicus GLO-30 is a surface model, so the cell under a building
centroid (or a covered MRT exit) carries the roof height and would essentially never flood.
Point assets are therefore counted as exposed when a flooded cell lies within ``NEAR_CELLS``
cells (~60 m) of the point. Roads are sampled on the road line itself (no buffer).
"""
from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import ndimage

NEAR_CELLS = 2
ROAD_STEP_DEG = 10.0 / 111_320.0   # ~10 m sampling interval along road lines


def rowcol(transform, lon, lat):
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    c = np.floor((lon - transform.c) / transform.a).astype(int)
    r = np.floor((lat - transform.f) / transform.e).astype(int)
    return r, c


def in_grid(r, c, shape):
    return (r >= 0) & (r < shape[0]) & (c >= 0) & (c < shape[1])


def near(mask: np.ndarray, cells: int = NEAR_CELLS) -> np.ndarray:
    if cells <= 0:
        return mask
    st = np.ones((2 * cells + 1, 2 * cells + 1), dtype=bool)
    return ndimage.binary_dilation(mask, structure=st)


def sample(mask: np.ndarray, r, c) -> np.ndarray:
    ok = in_grid(r, c, mask.shape)
    out = np.zeros(len(r), dtype=bool)
    out[ok] = mask[r[ok], c[ok]]
    return out


def load_buildings(path) -> pd.DataFrame:
    els = json.load(open(path, encoding="utf8"))["elements"]
    rows = [(e["center"]["lon"], e["center"]["lat"]) for e in els if "center" in e]
    return pd.DataFrame(rows, columns=["lon", "lat"])


def load_road_samples(path) -> pd.DataFrame:
    """Densify OSM road ways into points, each carrying the length (km) it represents."""
    els = json.load(open(path, encoding="utf8"))["elements"]
    lons, lats, lens = [], [], []
    for e in els:
        g = e.get("geometry")
        if not g or len(g) < 2:
            continue
        xy = np.array([(p["lon"], p["lat"]) for p in g])
        for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:]):
            dx = (x1 - x0) * np.cos(np.deg2rad(y0))
            dy = y1 - y0
            seg_deg = float(np.hypot(dx, dy))
            n = max(1, int(np.ceil(seg_deg / ROAD_STEP_DEG)))
            t = (np.arange(n) + 0.5) / n
            lons.append(x0 + t * (x1 - x0))
            lats.append(y0 + t * (y1 - y0))
            lens.append(np.full(n, seg_deg * 111.32 / n))
    return pd.DataFrame({"lon": np.concatenate(lons), "lat": np.concatenate(lats),
                         "km": np.concatenate(lens)})


def subzone_population(csv_path) -> pd.DataFrame:
    """Census 2020 residents by subzone (rows that are not '<PA> - Total' and not 'Total')."""
    d = pd.read_csv(csv_path)
    d = d[["Number", "Total_Total"]].rename(columns={"Number": "name", "Total_Total": "pop"})
    d["pop"] = pd.to_numeric(d["pop"].astype(str).str.replace(",", "").replace("-", "0"), errors="coerce").fillna(0)
    d = d[(d["name"].str.strip() != "Total") & ~d["name"].str.contains(r"-\s*Total\s*$", regex=True)]
    d["key"] = d["name"].str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    return d


@dataclass
class Assets:
    b_r: np.ndarray
    b_c: np.ndarray
    b_zone: np.ndarray        # subzone index per building (-1 outside)
    r_r: np.ndarray
    r_c: np.ndarray
    r_km: np.ndarray
    e_r: np.ndarray
    e_c: np.ndarray
    e_station: np.ndarray     # station name per exit
    zone_pop: np.ndarray      # Census population per subzone index
    zone_pa: np.ndarray       # planning-area name per subzone index


def exposure(flooded: np.ndarray, a: Assets, cell_area: np.ndarray, sg: np.ndarray) -> dict:
    fl = flooded & sg
    nearby = near(fl)
    b_hit = sample(nearby, a.b_r, a.b_c)
    r_hit = sample(fl, a.r_r, a.r_c)
    e_hit = sample(nearby, a.e_r, a.e_c)
    stations = sorted(set(a.e_station[e_hit]))
    nz = len(a.zone_pop)
    b_tot = np.bincount(a.b_zone[a.b_zone >= 0], minlength=nz)
    b_fl = np.bincount(a.b_zone[(a.b_zone >= 0) & b_hit], minlength=nz)
    share = np.divide(b_fl, b_tot, out=np.zeros(nz), where=b_tot > 0)
    pop = a.zone_pop * share
    area = float((fl * cell_area[:, None]).sum() / 1e6)
    return {"area_km2": area, "buildings": int(b_hit.sum()), "road_km": float(a.r_km[r_hit].sum()),
            "mrt_stations": len(stations), "station_names": stations, "population": float(pop.sum()),
            "pop_by_zone": pop, "bldg_by_zone": b_fl}
