"""2-D local-inertial storm-tide run for one coastal sub-domain, compared with the bathtub.

Called from run_pipeline.py (needs the DEM grids it builds). Writes outputs/inertial.json and
docs/data/model_*.json. Experiments:
  base        : observed 1999-12-23 storm-tide hydrograph shape, raised so its peak equals the
                scenario still-water level; Manning n = 0.035; 30 m grid; alpha = 0.7
  steady      : stage held at the peak level for the whole run (bathtub assumption) - the model
                should then approach the bathtub extent
  alpha_0.35  : half the time step (temporal convergence)
  dx_60/dx_90 : 2x / 3x coarser grid (spatial convergence)
  n_0.02/0.06 : roughness sensitivity
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from sgci import inertial, inundation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
WEB = ROOT / "docs" / "data"
HEADLINE = "rp100_ssp585_2100"
BOX_CELLS = 180          # sub-domain size (cells); ~5.6 km at 30 m
DEPTH_MIN = 0.05         # a cell counts as flooded if max depth >= 5 cm
N_BASE = 0.035


def choose_box(G, fill, wl):
    """Box centred on the densest cluster of bathtub-flooded mainland Singapore.

    Offshore-island hot spots (Pulau Tekong polder, Jurong Island) are excluded: much of that
    land was reclaimed or still under construction when the DEM was acquired."""
    fl = (fill < wl) & G["sg"] & ~np.isin(G["pa_idx"], G["island_idx"])
    from scipy import ndimage
    dens = ndimage.uniform_filter(fl.astype(float), size=BOX_CELLS)
    r, c = np.unravel_index(np.argmax(dens), dens.shape)
    r0 = int(np.clip(r - BOX_CELLS // 2, 0, fl.shape[0] - BOX_CELLS))
    c0 = int(np.clip(c - BOX_CELLS // 2, 0, fl.shape[1] - BOX_CELLS))
    return r0, c0


def stage_fn(hours, wl):
    t = hours * 3600.0
    return lambda s: float(np.interp(s, t - t[0], wl))


def run_case(z, sea, inland, sg, dx, n, hyd_t, hyd_wl, alpha=0.7, steady=False):
    if steady:
        st = (lambda s, v=float(hyd_wl.max()): v)
    else:
        st = stage_fn(hyd_t, hyd_wl)
    h0 = None
    t0 = time.time()
    r = inertial.run(z, dx, n, sea, st, t_end=(hyd_t[-1] - hyd_t[0]) * 3600.0, alpha=alpha, h0=h0)
    wet = (r.hmax >= DEPTH_MIN) & sg
    return {"area_km2": float(wet.sum() * dx * dx / 1e6), "volume_m3": float((r.hmax * sg).sum() * dx * dx),
            "steps": r.steps, "seconds": round(time.time() - t0, 1),
            "mass_error_rel": float(abs(r.volume_error) / max(abs(r.volume_in), 1.0))}, r


def main(G, fill, scen, tide_res):
    sc = next(s for s in scen if s["id"] == HEADLINE)
    wl = sc["wl"]
    r0, c0 = choose_box(G, fill, wl)
    sl = (slice(r0, r0 + BOX_CELLS), slice(c0, c0 + BOX_CELLS))
    z = G["z"][sl].copy()
    z[~np.isfinite(z)] = 50.0      # inland water bodies act as walls (held at regulated level)
    sea = G["sea"][sl].copy()
    inland = G["inland_zero"][sl].copy()
    sg = G["sg"][sl].copy()
    tr = G["tr"]
    dx = 30.0
    hyd = pd.read_csv(OUT / "hydrograph_1999_12_23.csv")
    hyd_t = hyd["hour"].to_numpy(float)
    hyd_wl = hyd["wl_m"].to_numpy() + (wl - hyd["wl_m"].max())
    # bathtub reference inside the box, same depth threshold
    bt = (fill[sl] < wl - DEPTH_MIN) & sg
    bath_km2 = float(bt.sum() * dx * dx / 1e6)
    res = {"scenario": HEADLINE, "wl": wl, "box_rows_cols": [r0, c0, BOX_CELLS],
           "box_lonlat": [round(tr.c + c0 * tr.a, 4), round(tr.f + (r0 + BOX_CELLS) * tr.e, 4),
                          round(tr.c + (c0 + BOX_CELLS) * tr.a, 4), round(tr.f + r0 * tr.e, 4)],
           "hydrograph_hours": float(hyd_t[-1] - hyd_t[0]), "depth_min_m": DEPTH_MIN,
           "bathtub_area_km2": round(bath_km2, 3), "runs": {}}
    cases = {"base": dict(n=N_BASE), "steady": dict(n=N_BASE, steady=True),
             "alpha_0.35": dict(n=N_BASE, alpha=0.35), "n_0.02": dict(n=0.02), "n_0.06": dict(n=0.06)}
    rbase = None
    for name, kw in cases.items():
        out, r = run_case(z, sea, inland, sg, dx, kw["n"], hyd_t, hyd_wl, alpha=kw.get("alpha", 0.7),
                          steady=kw.get("steady", False))
        res["runs"][name] = out
        print(name, out)
        if name == "base":
            rbase = r
    for k in (2, 3):
        zc = inertial.coarsen(z, k, "mean")
        seac = inertial.coarsen(sea.astype(float), k, "mean") > 0.5
        inlc = inertial.coarsen(inland.astype(float), k, "mean") > 0.5
        sgc = inertial.coarsen(sg.astype(float), k, "mean") > 0.5
        out, _ = run_case(zc, seac, inlc & ~seac, sgc & ~seac, dx * k, N_BASE, hyd_t, hyd_wl)
        res["runs"][f"dx_{int(dx * k)}"] = out
        print(k, out)
    b = res["runs"]["base"]["area_km2"]
    res["base_vs_bathtub_ratio"] = round(b / bath_km2, 3) if bath_km2 else None
    res["steady_vs_bathtub_ratio"] = round(res["runs"]["steady"]["area_km2"] / bath_km2, 3) if bath_km2 else None
    # agreement maps (cell counts) for the base run
    mb = (rbase.hmax >= DEPTH_MIN) & sg
    res["cells_both"] = int((mb & bt).sum())
    res["cells_bathtub_only"] = int((bt & ~mb).sum())
    res["cells_model_only"] = int((mb & ~bt).sum())
    (OUT / "inertial.json").write_text(json.dumps(res, indent=1) + "\n", encoding="utf8")
    np.save(ROOT / "data_cache" / "model_hmax.npy", rbase.hmax)
    # web layer: polygons of model-flooded vs bathtub-only cells
    import rasterio.features
    from affine import Affine
    from shapely.geometry import mapping, shape

    tsub = tr * Affine.translation(c0, r0)
    cls = np.zeros(z.shape, dtype=np.int16)
    cls[bt & ~mb] = 1
    cls[mb] = 2
    feats = []
    for geom, val in rasterio.features.shapes(cls, mask=cls > 0, transform=tsub, connectivity=4):
        g = shape(geom).simplify(0.00012, preserve_topology=True)
        gj = json.loads(json.dumps(mapping(g)), parse_float=lambda s: round(float(s), 5))
        feats.append({"type": "Feature", "properties": {"c": int(val)}, "geometry": gj})
    web = json.loads((WEB / "scenarios.json").read_text(encoding="utf8"))
    web["model_box"] = res["box_lonlat"]
    web["model"] = {k: res[k] for k in ("bathtub_area_km2", "base_vs_bathtub_ratio", "steady_vs_bathtub_ratio")}
    web["model"]["base_area_km2"] = res["runs"]["base"]["area_km2"]
    (WEB / "scenarios.json").write_text(json.dumps(web, indent=1) + "\n", encoding="utf8")
    (WEB / "model_compare.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                          separators=(",", ":")), encoding="utf8")
    return res
