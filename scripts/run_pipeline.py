"""End-to-end pipeline: open data -> extremes -> scenarios -> inundation -> exposure -> web data.

    python scripts/run_pipeline.py            # everything (downloads into data_cache/ first)
    python scripts/run_pipeline.py --skip-model

Writes small derived outputs to outputs/ and docs/data/ (both committed). Raw inputs stay in
data_cache/ (git-ignored). Runtime on a laptop CPU: ~10 minutes, most of it the 2-D model.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sgci import dem as demmod  # noqa: E402
from sgci import eva, exposure, fetch, inundation, slr, tide  # noqa: E402

OUT = ROOT / "outputs"
WEB = ROOT / "docs" / "data"
CACHE = ROOT / "data_cache"
RETURN_PERIODS = [2, 10, 50, 100, 1000]
MAP_RPS = [10, 100]
HORIZONS = [2050, 2100]
SSPS = ["ssp126", "ssp245", "ssp585"]
SSP_LABEL = {"ssp126": "SSP1-2.6", "ssp245": "SSP2-4.5", "ssp585": "SSP5-8.5"}
N_BOOT = 1000
DATUM_SHIFTS = [0.0, -0.30, 0.30, None]   # offsets added to the MDT; None = no MDT (naive)
AREA_CURVE_MAX_MM = 7000
# planning areas that are offshore islands (Pulau Tekong & Ubin; Jurong Island & Bukom; Sentosa etc.)
OFFSHORE = ["NORTH-EASTERN ISLANDS", "WESTERN ISLANDS", "SOUTHERN ISLANDS"]


def jdump(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, default=float) + "\n", encoding="utf8")


# ----------------------------------------------------------------------------- 1. extremes
def stage_tide():
    s = tide.read_uhslc_csv(CACHE / "uhslc" / "h699a.csv")
    am = tide.annual_means(s)
    msl0 = tide.baseline_msl(am)
    tr = tide.linear_trend(am)
    res = {"record_start": str(s.index.min()), "record_end": str(s.index.max()),
           "n_hours": int(len(s)), "msl_1995_2014_above_cd_m": round(msl0, 4),
           "trend_mm_per_yr": round(tr.slope_m_per_yr * 1000, 2), "trend_n_years": tr.n_years}
    # trend standard error (OLS, no autocorrelation correction - stated in README)
    v = am[am["valid"]]
    x = v.index.to_numpy(float) + 0.5 - tide.BASELINE_MID
    resid = v["mean"].to_numpy() - (tr.intercept_at_mid + tr.slope_m_per_yr * x)
    se = np.sqrt(resid.var(ddof=2) / ((x - x.mean()) ** 2).sum())
    res["trend_se_mm_per_yr"] = round(se * 1000, 2)
    res["observed_rise_2005_2020_m"] = round(tr.slope_m_per_yr * 15, 3)
    fits = {}
    for method in ["linear", "annual_mean"]:
        x_d = tide.detrend(s, msl0, tr, am, method=method)
        sm = tide.season_maxima(x_d)
        good = sm[sm["valid"]]
        g = eva.fit_gev(good["max"].to_numpy())
        rl = eva.gev_return_level(g, RETURN_PERIODS)
        boot = eva.bootstrap_gev(good["max"].to_numpy(), RETURN_PERIODS, n_boot=N_BOOT)
        lr = eva.lr_test_gumbel(good["max"].to_numpy())
        fits[method] = {
            "n_seasons": int(len(good)), "excluded_seasons": [int(i) for i in sm.index[~sm["valid"]]],
            "mu": round(g.mu, 4), "sigma": round(g.sigma, 4), "xi": round(g.xi, 4),
            "xi_ci95": [round(boot["xi_lo"], 3), round(boot["xi_hi"], 3)],
            "return_levels": {str(T): {"rl": round(float(r), 3), "lo": round(float(lo), 3), "hi": round(float(hi), 3)}
                              for T, r, lo, hi in zip(RETURN_PERIODS, rl, boot["rl_lo"], boot["rl_hi"])},
            "gumbel_lr_p": round(lr["p_value"], 3),
            "max_observed": round(float(good["max"].max()), 3),
            "max_observed_time": str(good.loc[good["max"].idxmax(), "time"]),
        }
        if method == "linear":
            out = good[["max", "time", "coverage"]].copy()
            out["max"] = out["max"].round(4)
            out["coverage"] = out["coverage"].round(4)
            out.to_csv(OUT / "season_maxima.csv")
            months = pd.to_datetime(good["time"]).dt.month
            fits[method]["share_maxima_nov_to_feb"] = round(float(months.isin([11, 12, 1, 2]).mean()), 3)
            # POT / GPD cross-check on the same detrended hourly series
            vals = x_d.to_numpy()
            th = (x_d.index - x_d.index[0]).total_seconds().to_numpy() / 3600.0
            u = float(np.quantile(vals, 0.999))
            peaks = eva.decluster_runs(vals, th, u, run_h=72.0)
            years = len(x_d) / 8766.0
            gp = eva.fit_gpd(peaks, u, years)
            fits["pot_gpd"] = {"threshold_m": round(u, 3), "threshold_quantile": 0.999, "n_clusters": gp.n,
                               "rate_per_yr": round(gp.rate, 3), "sigma": round(gp.sigma, 4),
                               "xi": round(gp.xi, 4),
                               "return_levels": {str(T): round(float(r), 3) for T, r in
                                                 zip(RETURN_PERIODS, eva.gpd_return_level(gp, RETURN_PERIODS))}}
            # storm-tide hydrograph shape for the 2-D model: +-36 h around the largest event
            t0 = pd.Timestamp(fits[method]["max_observed_time"])
            win = x_d[(x_d.index >= t0 - pd.Timedelta(hours=36)) & (x_d.index <= t0 + pd.Timedelta(hours=36))]
            hyd = pd.DataFrame({"hour": ((win.index - t0).total_seconds() / 3600).astype(int),
                                "wl_m": win.round(3).to_numpy()})
            hyd.to_csv(OUT / "hydrograph_1999_12_23.csv", index=False)
    res["fits"] = fits
    jdump(res, OUT / "tide.json")
    return res


# ----------------------------------------------------------------------------- 2. AR6
def stage_slr():
    proj = {s: slr.read_location(CACHE / "ar6" / f"total_{s}_medium_confidence_values.nc") for s in SSPS}
    out = {"psmsl_id": slr.TANJONG_PAGAR_PSMSL, "baseline": "1995-2014", "units": "m",
           "workflow": "medium confidence, total relative sea-level change incl. vertical land motion",
           "scenarios": {}}
    for s, p in proj.items():
        out["lat"], out["lon"] = round(p["lat"], 2), round(p["lon"], 2)
        out["scenarios"][s] = {str(y): {q: round(slr.value(p, y, q), 3) for q in slr.QUANTILES}
                               for y in (2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100)}
    jdump(out, OUT / "ar6_tanjong_pagar.json")
    return out


# ----------------------------------------------------------------------------- 3. scenarios
def build_scenarios(tide_res, slr_res):
    fit = tide_res["fits"]["linear"]
    scen = []
    states = [("2020", "ssp245", 2020)] + [(f"{SSP_LABEL[s]} {y}", s, y) for y in HORIZONS for s in SSPS]
    for T in MAP_RPS:
        rl = fit["return_levels"][str(T)]
        for label, s, y in states:
            q = slr_res["scenarios"][s][str(y)]
            sid = f"rp{T}_{'2020' if y == 2020 else s + '_' + str(y)}"
            scen.append({
                "id": sid, "label": f"{T}-yr event, {label}", "rp": T, "ssp": None if y == 2020 else s,
                "year": y, "rl": rl["rl"], "slr_p50": q["p50"], "slr_p17": q["p17"], "slr_p83": q["p83"],
                "wl": round(rl["rl"] + q["p50"], 3),
                "wl_slr_lo": round(rl["rl"] + q["p17"], 3), "wl_slr_hi": round(rl["rl"] + q["p83"], 3),
                "wl_rl_lo": round(rl["lo"] + q["p50"], 3), "wl_rl_hi": round(rl["hi"] + q["p50"], 3),
            })
    return scen


# ----------------------------------------------------------------------------- 4/5. grids
def load_grids():
    import geopandas as gpd
    from rasterio.features import rasterize

    tiles = [CACHE / "dem" / f"{t}.tif" for t in fetch.DEM_TILES]
    raw, tr = demmod.load_dem(tiles)
    mdt = demmod.mdt_offset(CACHE / "mdt" / "dtu22mdt.grd.gz")
    delta = round(mdt["mean"], 2)
    pa = gpd.read_file(CACHE / "datagov" / "planning_areas.geojson")
    sg_all = rasterize([(g, 1) for g in pa.geometry], out_shape=raw.shape, transform=tr, fill=0).astype(bool)
    sea, inland_zero = inundation.split_sea(raw, sg_all)
    plateau = inundation.flat_plateaus(raw, sg_all) & ~sea
    inland_zero = inland_zero | plateau          # all inland water bodies
    # inland water is held at its regulated level: it neither floods nor conveys sea water
    z = np.where(inland_zero, np.inf, raw - delta)
    pa_idx = rasterize([(g, i + 1) for i, g in enumerate(pa.geometry)], out_shape=raw.shape,
                       transform=tr, fill=0).astype(np.int32) - 1
    sg = sg_all & ~inland_zero & ~sea
    island_idx = np.flatnonzero(pa["PLN_AREA_N"].isin(OFFSHORE).to_numpy())
    return dict(raw=raw, tr=tr, mdt=mdt, delta=delta, sea=sea, inland_zero=inland_zero, z=z, sg=sg,
                pa=pa, pa_idx=pa_idx, island_idx=island_idx)


def load_assets(G):
    import geopandas as gpd

    tr, shape = G["tr"], G["raw"].shape
    b = exposure.load_buildings(CACHE / "osm" / "buildings.json")
    zones = gpd.read_file(CACHE / "datagov" / "subzones.geojson")
    pts = gpd.GeoDataFrame(b, geometry=gpd.points_from_xy(b.lon, b.lat), crs="EPSG:4326")
    j = gpd.sjoin(pts, zones[["geometry"]].reset_index(names="zi"), how="left", predicate="within")
    j = j[~j.index.duplicated()]
    b_zone = j["zi"].fillna(-1).astype(int).to_numpy()
    pop = exposure.subzone_population(CACHE / "datagov" / "population_census2020.csv")
    key = zones["SUBZONE_N"].str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    popmap = dict(zip(pop["key"], pop["pop"]))
    zone_pop = key.map(popmap).fillna(0).to_numpy(float)
    unmatched = sorted(set(pop["key"]) - set(key))
    roads = exposure.load_road_samples(CACHE / "osm" / "roads.json")
    ex = gpd.read_file(CACHE / "datagov" / "mrt_exits.geojson")
    br, bc = exposure.rowcol(tr, b.lon, b.lat)
    rr, rc = exposure.rowcol(tr, roads.lon, roads.lat)
    er, ec = exposure.rowcol(tr, ex.geometry.x, ex.geometry.y)
    sg = G["sg"] | G["inland_zero"] | G["sea"]  # keep roads inside Singapore polygons only
    road_in = exposure.sample(G["pa_idx"] >= 0, rr, rc)
    assets = exposure.Assets(b_r=br, b_c=bc, b_zone=b_zone, r_r=rr[road_in], r_c=rc[road_in],
                             r_km=roads.km.to_numpy()[road_in], e_r=er, e_c=ec,
                             e_station=ex["STATION_NA"].to_numpy(), zone_pop=zone_pop,
                             zone_pa=zones["PLN_AREA_N"].to_numpy())
    meta = {"osm_buildings": int(len(b)), "buildings_in_subzones": int((b_zone >= 0).sum()),
            "road_km_in_sg": round(float(assets.r_km.sum()), 1), "mrt_exits": int(len(ex)),
            "mrt_stations": int(ex["STATION_NA"].nunique()), "census_pop_subzone_sum": int(pop["pop"].sum()),
            "census_pop_matched": int(zone_pop.sum()), "unmatched_subzones": unmatched}
    del sg
    return assets, meta


def stage_inundation(scen):
    G = load_grids()
    cache_fill = CACHE / f"fill_{G['delta']:.2f}.npy"
    if cache_fill.exists():
        fill = np.load(cache_fill)
    else:
        fill = inundation.fill_level(G["z"], G["sea"], cap=AREA_CURVE_MAX_MM / 1000 + 1.0)
        np.save(cache_fill, fill)
    area = demmod.cell_area_m2(G["tr"], G["raw"].shape[0])[:, None] * np.ones(G["raw"].shape)
    sg = G["sg"]
    # area curve (exact on a 1 mm grid): connected and naive bathtub
    levels_mm = np.arange(0, AREA_CURVE_MAX_MM + 1)
    def curve(vals):
        v = vals[sg]
        a = area[sg]
        o = np.argsort(v)
        vs, cs = v[o], np.cumsum(a[o])
        k = np.searchsorted(vs, levels_mm / 1000.0, side="left")   # count of v < level
        return np.where(k > 0, cs[np.maximum(k - 1, 0)], 0.0) / 1e6
    ac = pd.DataFrame({"level_mm": levels_mm, "connected_km2": curve(fill).round(4),
                       "naive_km2": curve(G["z"]).round(4)})
    ac.to_csv(OUT / "area_curve.csv", index=False)
    sg_land_km2 = float(area[sg].sum() / 1e6)

    assets, meta = load_assets(G)
    rows, pa_rows = [], []
    for sc in scen:
        fl = inundation.flooded_from_fill(fill, G["sea"], sc["wl"])
        e = exposure.exposure(fl, assets, area[:, 0], sg)
        sc_row = {"id": sc["id"], "wl": sc["wl"], "area_km2": round(e["area_km2"], 3),
                  "buildings": e["buildings"], "population": int(round(e["population"], -1)),
                  "road_km": round(e["road_km"], 2), "mrt_stations": e["mrt_stations"],
                  "station_names": ";".join(s.title() for s in e["station_names"])}
        for k in ("wl_slr_lo", "wl_slr_hi"):
            fk = inundation.flooded_from_fill(fill, G["sea"], sc[k])
            ek = exposure.exposure(fk, assets, area[:, 0], sg)
            tag = k.replace("wl_", "")
            sc_row[f"area_km2_{tag}"] = round(ek["area_km2"], 3)
            sc_row[f"buildings_{tag}"] = ek["buildings"]
            sc_row[f"population_{tag}"] = int(round(ek["population"], -1))
        isl = np.isin(G["pa_idx"], G["island_idx"])
        sc_row["area_km2_mainland"] = round(float(area[fl & sg & ~isl].sum() / 1e6), 3)
        rows.append(sc_row)
        # per planning area (central estimate)
        flsg = fl & sg
        pa_area = np.bincount(G["pa_idx"][flsg], weights=area[flsg], minlength=len(G["pa"])) / 1e6
        df = pd.DataFrame({"zone_pa": assets.zone_pa, "pop": e["pop_by_zone"], "bldg": e["bldg_by_zone"]})
        agg = df.groupby("zone_pa").sum()
        for i, name in enumerate(G["pa"]["PLN_AREA_N"]):
            pa_rows.append({"id": sc["id"], "planning_area": name.title(), "area_km2": round(float(pa_area[i]), 3),
                            "buildings": int(agg["bldg"].get(name, 0)),
                            "population": int(round(float(agg["pop"].get(name, 0.0)), -1))})
    ex_df = pd.DataFrame(rows)
    ex_df.to_csv(OUT / "exposure_scenarios.csv", index=False)
    pa_df = pd.DataFrame(pa_rows)
    pa_df.to_csv(OUT / "exposure_by_planning_area.csv", index=False)

    # datum sensitivity: fill shifts by a constant when the datum offset changes
    sens = []
    for sc in scen:
        r = {"id": sc["id"]}
        for d in DATUM_SHIFTS:
            delta_alt = 0.0 if d is None else G["delta"] + d
            wl_eff = sc["wl"] + delta_alt - G["delta"]        # fill_alt = fill + (delta - delta_alt)
            key = "naive_datum" if d is None else f"mdt{d:+.2f}"
            r[key] = round(float(area[(fill < wl_eff) & sg].sum() / 1e6), 3)
        r["no_connectivity"] = round(float(area[(G["z"] < sc["wl"]) & sg].sum() / 1e6), 3)
        r["eight_connected"] = None
        sens.append(r)
    # 8-connectivity for the headline scenarios only (one extra priority flood)
    fill8 = inundation.fill_level(G["z"], G["sea"], cap=4.5, eight=True)
    for r, sc in zip(sens, scen):
        r["eight_connected"] = round(float(area[(fill8 < sc["wl"]) & sg].sum() / 1e6), 3)
    pd.DataFrame(sens).to_csv(OUT / "sensitivity_area.csv", index=False)

    # simple plausibility check of the datum choice: share of land below 5 m above MSL
    zl = G["z"][sg]
    al = area[sg]
    below5 = {f"{lab}": round(float(al[(zl + G['delta'] - dd) < 5.0].sum() / al.sum()), 3)
              for lab, dd in (("mdt", G["delta"]), ("no_shift", 0.0))}
    grid_meta = {"mdt": G["mdt"], "datum_offset_m": G["delta"], "sg_land_km2": round(sg_land_km2, 1),
                 "share_land_below_5m": below5, "cell_size_m": round(float(np.sqrt(area[600, 900])), 1),
                 "grid_shape": list(G["raw"].shape), **meta}
    jdump(grid_meta, OUT / "grid_meta.json")
    return G, fill, area, ex_df, pa_df


# ----------------------------------------------------------------------------- 6. web data
def stage_web(G, fill, scen, ex_df, pa_df, tide_res, slr_res):
    import rasterio.features
    from shapely.geometry import mapping, shape

    levels = sorted({sc["wl"] for sc in scen})
    band = np.full(fill.shape, 0, dtype=np.int16)
    prev = -np.inf
    for k, L in enumerate(levels, start=1):
        band[(fill >= prev) & (fill < L) & G["sg"]] = k
        prev = L
    feats = []
    for geom, val in rasterio.features.shapes(band, mask=band > 0, transform=G["tr"], connectivity=4):
        g = shape(geom).simplify(0.00012, preserve_topology=True)
        if g.is_empty:
            continue
        gj = mapping(g)
        gj = json.loads(json.dumps(gj), parse_float=lambda s: round(float(s), 5))
        feats.append({"type": "Feature", "properties": {"b": int(val)}, "geometry": gj})
    (WEB / "flood_bands.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                        separators=(",", ":")), encoding="utf8")
    pa = G["pa"][["PLN_AREA_N", "geometry"]].copy()
    pa["geometry"] = pa.geometry.simplify(0.0003, preserve_topology=True)
    pa["name"] = pa["PLN_AREA_N"].str.title()
    pa[["name", "geometry"]].to_file(WEB / "planning_areas.geojson", driver="GeoJSON", COORDINATE_PRECISION=5)
    web = {"levels": levels, "scenarios": [], "tide": {"rl": tide_res["fits"]["linear"]["return_levels"]},
           "ar6": slr_res["scenarios"]}
    exi = ex_df.set_index("id")
    for sc in scen:
        e = exi.loc[sc["id"]]
        pas = pa_df[(pa_df["id"] == sc["id"]) & (pa_df["area_km2"] > 0)].sort_values("area_km2", ascending=False)
        web["scenarios"].append({**{k: sc[k] for k in ("id", "label", "rp", "ssp", "year", "wl", "rl", "slr_p50",
                                                         "slr_p17", "slr_p83")},
                                 "band": levels.index(sc["wl"]) + 1,
                                 "area_km2": float(e["area_km2"]), "area_lo": float(e["area_km2_slr_lo"]),
                                 "area_hi": float(e["area_km2_slr_hi"]), "buildings": int(e["buildings"]),
                                 "population": int(e["population"]), "road_km": float(e["road_km"]),
                                 "mrt_stations": int(e["mrt_stations"]),
                                 "station_names": [s for s in str(e["station_names"]).split(";") if s and s != "nan"],
                                 "by_pa": pas[["planning_area", "area_km2", "buildings", "population"]].head(12)
                                 .to_dict(orient="records")})
    inert = OUT / "inertial.json"
    if inert.exists():
        res = json.loads(inert.read_text())
        web["model_box"] = res["box_lonlat"]
        web["model"] = {k: res[k] for k in ("bathtub_area_km2", "base_vs_bathtub_ratio", "steady_vs_bathtub_ratio")}
        web["model"]["base_area_km2"] = res["runs"]["base"]["area_km2"]
    jdump(web, WEB / "scenarios.json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-model", action="store_true")
    ap.add_argument("--no-fetch", action="store_true")
    ap.add_argument("--model-only", action="store_true", help="re-run only the 2-D model from cached grids")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    if args.model_only:
        import run_model  # noqa: E402
        G = load_grids()
        fill = np.load(CACHE / f"fill_{G['delta']:.2f}.npy")
        scen = json.loads((OUT / "scenarios.json").read_text())
        run_model.main(G, fill, scen, json.loads((OUT / "tide.json").read_text()))
        return
    WEB.mkdir(parents=True, exist_ok=True)
    if not args.no_fetch:
        fetch.main()
    t = stage_tide()
    print("tide done", t["fits"]["linear"]["return_levels"]["100"])
    s = stage_slr()
    scen = build_scenarios(t, s)
    jdump(scen, OUT / "scenarios.json")
    G, fill, area, ex_df, pa_df = stage_inundation(scen)
    print(ex_df[["id", "wl", "area_km2", "buildings", "population", "road_km", "mrt_stations"]].to_string())
    stage_web(G, fill, scen, ex_df, pa_df, t, s)
    if not args.skip_model:
        import run_model  # noqa: E402
        run_model.main(G, fill, scen, t)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    main()
