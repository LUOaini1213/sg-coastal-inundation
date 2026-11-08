"""Recompute every number quoted in README.md from committed outputs and compare.

README.md is rendered from README.template.md: each ``{{key}}`` is replaced by a value this
script derives from outputs/*.json|csv. Where a number can be recomputed independently of the
pipeline (return levels from the GEV parameters, water levels from RL + SLR, flooded areas
from the 1 mm area curve, model ratios, the analytical-solution test) it is, and the
recomputation must agree with the stored output before the value is accepted.

    python tools/check_readme_numbers.py          # exit 1 on any mismatch
    python tools/check_readme_numbers.py --write  # regenerate README.md
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from sgci import eva, inertial  # noqa: E402

OUT = ROOT / "outputs"
HEADLINE = ["rp100_2020", "rp100_ssp245_2050", "rp100_ssp585_2100", "rp10_2020", "rp100_ssp126_2100",
            "rp100_ssp245_2100", "rp10_ssp585_2100"]


class Mismatch(Exception):
    pass


def agree(a, b, tol, what):
    if not np.isfinite(a) or not np.isfinite(b) or abs(a - b) > tol:
        raise Mismatch(f"{what}: recomputed {a} vs stored {b} (tol {tol})")


def f(v, d):
    return f"{v:,.{d}f}"


def pct(v, d=0):
    return f"{100 * v:.{d}f}%"


def hunter_rmse():
    u, n, dx = 1.0, 0.03, 10.0
    nx = int(5000 / dx)
    z = np.zeros((3, nx))
    b = np.zeros_like(z, bool)
    b[:, 0] = True
    r = inertial.run(z, dx, n, b, lambda t: float(inertial.hunter_front(-dx / 2, t, u, n)), 3600.0)
    x = np.arange(nx) * dx - dx / 2
    ha = inertial.hunter_front(x, 3600.0, u, n)
    m = ha > 0
    return float(np.sqrt(np.mean((r.h[1][m] - ha[m]) ** 2))), float(ha.max())


def values() -> dict:
    V = {}
    t = json.loads((OUT / "tide.json").read_text())
    lin = t["fits"]["linear"]
    am = t["fits"]["annual_mean"]
    gp = t["fits"]["pot_gpd"]
    sm = pd.read_csv(OUT / "season_maxima.csv")
    a6 = json.loads((OUT / "ar6_tanjong_pagar.json").read_text())["scenarios"]
    scen = {s["id"]: s for s in json.loads((OUT / "scenarios.json").read_text())}
    ex = pd.read_csv(OUT / "exposure_scenarios.csv").set_index("id")
    pa = pd.read_csv(OUT / "exposure_by_planning_area.csv")
    sens = pd.read_csv(OUT / "sensitivity_area.csv").set_index("id")
    ac = pd.read_csv(OUT / "area_curve.csv")
    gm = json.loads((OUT / "grid_meta.json").read_text())
    mo = json.loads((OUT / "inertial.json").read_text())

    # ---- tide gauge
    agree(len(sm), lin["n_seasons"], 0, "season count")
    agree(sm["max"].max(), lin["max_observed"], 0.0006, "max observed")
    V["n_seasons"] = str(lin["n_seasons"])
    V["n_hours"] = f(t["n_hours"], 0)
    V["msl_cd"] = f(t["msl_1995_2014_above_cd_m"], 2)
    V["trend"] = f(t["trend_mm_per_yr"], 1)
    V["trend_se"] = f(t["trend_se_mm_per_yr"], 1)
    agree(t["trend_mm_per_yr"] * 15 / 1000, t["observed_rise_2005_2020_m"], 0.0006, "observed rise")
    V["obs_rise_2020"] = f(t["observed_rise_2005_2020_m"] * 100, 1)
    V["ar6_2020"] = f(a6["ssp245"]["2020"]["p50"] * 100, 1)
    V["max_obs"] = f(lin["max_observed"], 2)
    V["max_obs_date"] = lin["max_observed_time"][:10]
    V["share_ne_monsoon"] = pct(lin["share_maxima_nov_to_feb"])
    g = eva.GEV(lin["mu"], lin["sigma"], lin["xi"], 0, 0)
    for T in ("10", "100", "1000"):
        rl = float(eva.gev_return_level(g, int(T)))
        agree(rl, lin["return_levels"][T]["rl"], 0.002, f"RL{T} from GEV parameters")
        V[f"rl{T}"] = f(lin["return_levels"][T]["rl"], 2)
        V[f"rl{T}_lo"] = f(lin["return_levels"][T]["lo"], 2)
        V[f"rl{T}_hi"] = f(lin["return_levels"][T]["hi"], 2)
    V["gev_mu"], V["gev_sigma"], V["gev_xi"] = f(lin["mu"], 3), f(lin["sigma"], 3), f(lin["xi"], 2)
    V["xi_lo"], V["xi_hi"] = f(lin["xi_ci95"][0], 2), f(lin["xi_ci95"][1], 2)
    V["gumbel_p"] = f(lin["gumbel_lr_p"], 2)
    V["rl100_am"] = f(am["return_levels"]["100"]["rl"], 2)
    V["rl100_am_hi"] = f(am["return_levels"]["100"]["hi"], 2)
    V["rl100_gpd"] = f(gp["return_levels"]["100"], 2)
    V["gpd_n"] = str(gp["n_clusters"])
    V["gpd_rate"] = f(gp["rate_per_yr"], 1)

    # ---- AR6
    for s in ("ssp126", "ssp245", "ssp585"):
        for y in ("2050", "2100"):
            q = a6[s][y]
            V[f"slr_{s}_{y}"] = f(q["p50"], 2)
            V[f"slr_{s}_{y}_lo"] = f(q["p17"], 2)
            V[f"slr_{s}_{y}_hi"] = f(q["p83"], 2)

    # ---- scenarios, areas (recomputed from the 1 mm area curve)
    curve = ac.set_index("level_mm")
    def area_at(wl, col="connected_km2"):
        return float(curve.loc[int(round(wl * 1000)), col])
    for sid, sc in scen.items():
        y = "2020" if sc["year"] == 2020 else str(sc["year"])
        s = sc["ssp"] or "ssp245"
        agree(sc["rl"] + a6[s][y]["p50"], sc["wl"], 0.0011, f"{sid} water level")
        agree(area_at(sc["wl"]), ex.loc[sid, "area_km2"], 0.0006, f"{sid} area vs curve")
        for d, col in ((0.3, "mdt+0.30"), (-0.3, "mdt-0.30")):
            agree(area_at(sc["wl"] + d), sens.loc[sid, col], 0.0006, f"{sid} datum {col}")
        agree(area_at(sc["wl"], "naive_km2"), sens.loc[sid, "no_connectivity"], 0.0006, f"{sid} naive")
        p = pa[pa["id"] == sid]
        agree(p["area_km2"].sum(), ex.loc[sid, "area_km2"], 0.03, f"{sid} PA area sum")
        agree(p["population"].sum(), ex.loc[sid, "population"], 300, f"{sid} PA population sum")
    for sid in HEADLINE:
        k = sid.replace("rp", "")
        e = ex.loc[sid]
        V[f"wl_{k}"] = f(scen[sid]["wl"], 2)
        V[f"area_{k}"] = f(area_at(scen[sid]["wl"]), 1)
        V[f"area_{k}_lo"] = f(e["area_km2_slr_lo"], 1)
        V[f"area_{k}_hi"] = f(e["area_km2_slr_hi"], 1)
        V[f"main_{k}"] = f(e["area_km2_mainland"], 1)
        V[f"bldg_{k}"] = f(e["buildings"], 0)
        V[f"pop_{k}"] = f(e["population"], 0)
        V[f"road_{k}"] = f(e["road_km"], 0)
        V[f"mrt_{k}"] = str(int(e["mrt_stations"]))
        V[f"naive_{k}"] = f(sens.loc[sid, "no_connectivity"], 1)
        V[f"d0_{k}"] = f(sens.loc[sid, "naive_datum"], 1)
        V[f"dlo_{k}"] = f(sens.loc[sid, "mdt-0.30"], 1)
        V[f"dhi_{k}"] = f(sens.loc[sid, "mdt+0.30"], 1)
        V[f"e8_{k}"] = f(sens.loc[sid, "eight_connected"], 1)
    h = "100_ssp585_2100"
    V["naive_ratio"] = f(sens.loc["rp" + h, "no_connectivity"] / ex.loc["rp" + h, "area_km2"], 1)
    V["mrt_names_" + h] = ", ".join(n.replace(" Mrt Station", "") for n in str(ex.loc["rp" + h, "station_names"]).split(";"))
    top = pa[pa["id"] == "rp" + h].sort_values("area_km2", ascending=False).head(3)
    V["top_pa"] = ", ".join(f"{r.planning_area} ({r.area_km2:.1f} km²)" for r in top.itertuples())
    isl_share = 1 - ex.loc["rp" + h, "area_km2_mainland"] / ex.loc["rp" + h, "area_km2"]
    V["island_share"] = pct(isl_share)

    # ---- grid / datum
    V["mdt"] = f(gm["datum_offset_m"], 2)
    V["below5_mdt"] = pct(gm["share_land_below_5m"]["mdt"])
    V["below5_raw"] = pct(gm["share_land_below_5m"]["no_shift"])
    V["sg_land_km2"] = f(gm["sg_land_km2"], 0)
    V["cell_size"] = f(gm["cell_size_m"], 0)
    V["n_buildings"] = f(gm["osm_buildings"], 0)
    V["n_exits"] = str(gm["mrt_exits"])
    V["n_stations"] = str(gm["mrt_stations"])
    V["road_km_total"] = f(gm["road_km_in_sg"], 0)

    # ---- 2-D model
    runs = mo["runs"]
    agree(runs["base"]["area_km2"] / mo["bathtub_area_km2"], mo["base_vs_bathtub_ratio"], 0.0006, "model ratio")
    agree(runs["steady"]["area_km2"] / mo["bathtub_area_km2"], mo["steady_vs_bathtub_ratio"], 0.0006, "steady ratio")
    V["m_box"] = "{:.3f}–{:.3f}°E, {:.3f}–{:.3f}°N".format(*[mo["box_lonlat"][i] for i in (0, 2, 1, 3)])
    V["m_bath"] = f(mo["bathtub_area_km2"], 2)
    for k in ("base", "steady", "alpha_0.35", "n_0.02", "n_0.06", "dx_60", "dx_90"):
        V[f"m_{k}"] = f(runs[k]["area_km2"], 2)
    V["m_ratio"] = pct(mo["base_vs_bathtub_ratio"])
    V["m_steady_ratio"] = pct(mo["steady_vs_bathtub_ratio"])
    V["m_reduction"] = pct(1 - mo["base_vs_bathtub_ratio"])
    V["m_alpha_diff"] = pct(abs(runs["alpha_0.35"]["area_km2"] / runs["base"]["area_km2"] - 1), 1)
    V["m_dx90_diff"] = pct(abs(runs["dx_90"]["area_km2"] / runs["base"]["area_km2"] - 1), 1)
    V["m_mass"] = f"{max(r['mass_error_rel'] for r in runs.values()):.0e}"
    rmse, hmx = hunter_rmse()
    V["hunter_rmse"] = f(rmse * 100, 1)
    V["hunter_rel"] = pct(rmse / hmx, 1)
    return V


def render(V: dict) -> str:
    tpl = (ROOT / "README.template.md").read_text(encoding="utf8")
    missing = sorted(set(re.findall(r"\{\{(\w[\w.]*)\}\}", tpl)) - set(V))
    if missing:
        raise Mismatch(f"template keys without a computed value: {missing}")
    return re.sub(r"\{\{(\w[\w.]*)\}\}", lambda m: V[m.group(1)], tpl)


def main():
    try:
        text = render(values())
    except Mismatch as e:
        print(f"FAIL: {e}")
        return 1
    readme = ROOT / "README.md"
    if "--write" in sys.argv:
        readme.write_text(text, encoding="utf8")
        print("README.md written")
        return 0
    cur = readme.read_text(encoding="utf8")
    if cur != text:
        import difflib
        for line in list(difflib.unified_diff(cur.splitlines(), text.splitlines(), "README.md", "recomputed",
                                              lineterm=""))[:40]:
            print(line)
        print("FAIL: README.md numbers differ from recomputed values")
        return 1
    n = len(re.findall(r"\{\{", (ROOT / "README.template.md").read_text(encoding="utf8")))
    print(f"OK: {n} README numbers recomputed and matched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
