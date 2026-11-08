"""Consistency of the committed outputs (no raw data needed)."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
sys.path.insert(0, str(ROOT / "tools"))


def test_readme_values_recompute_without_mismatch():
    import check_readme_numbers as c
    v = c.values()
    assert len(v) > 50


def test_area_curve_monotone_and_connected_le_naive():
    ac = pd.read_csv(OUT / "area_curve.csv")
    assert (np.diff(ac["connected_km2"]) >= -1e-9).all()
    assert (np.diff(ac["naive_km2"]) >= -1e-9).all()
    assert (ac["connected_km2"] <= ac["naive_km2"] + 1e-9).all()


def test_scenario_ordering():
    sc = {s["id"]: s for s in json.loads((OUT / "scenarios.json").read_text())}
    for k in ("2020", "ssp126_2050", "ssp585_2100"):
        assert sc[f"rp100_{k}"]["wl"] > sc[f"rp10_{k}"]["wl"]
    assert sc["rp100_ssp585_2100"]["wl"] > sc["rp100_ssp245_2100"]["wl"] > sc["rp100_ssp126_2100"]["wl"]
    for s in sc.values():
        assert s["wl_slr_lo"] <= s["wl"] <= s["wl_slr_hi"]
        assert s["wl_rl_lo"] <= s["wl"] <= s["wl_rl_hi"]


def test_web_bands_match_scenario_levels():
    w = json.loads((ROOT / "docs" / "data" / "scenarios.json").read_text())
    for s in w["scenarios"]:
        assert w["levels"][s["band"] - 1] == s["wl"]
    assert w["levels"] == sorted(w["levels"])
    gj = json.loads((ROOT / "docs" / "data" / "flood_bands.geojson").read_text())
    bands = {f["properties"]["b"] for f in gj["features"]}
    assert bands <= set(range(1, len(w["levels"]) + 1))


def test_exposure_monotone_in_water_level():
    ex = pd.read_csv(OUT / "exposure_scenarios.csv").sort_values("wl")
    assert (np.diff(ex["area_km2"]) >= -1e-9).all()
    assert (np.diff(ex["buildings"]) >= 0).all()
