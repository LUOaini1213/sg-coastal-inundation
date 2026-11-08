import gzip

import numpy as np
import pandas as pd
import pytest
from affine import Affine

from sgci import dem, exposure, inundation, slr


# ------------------------------------------------------------------ inundation

def _basin():
    """Sea on the left (0), a ridge at 2.0 m, a closed basin at 0.5 m behind it."""
    z = np.full((7, 9), 3.0)
    z[:, 0:2] = 0.0          # sea
    z[:, 2] = 1.0            # low coastal strip
    z[1:6, 3] = 2.0          # ridge
    z[2:5, 4:7] = 0.5        # basin
    return z


def test_sea_mask_edge_connected_only():
    z = np.ones((6, 6))
    z[:, 0] = 0.0
    z[3, 3] = 0.0            # enclosed zero cell (reservoir)
    sea = inundation.sea_mask(z)
    assert sea[:, 0].all()
    assert not sea[3, 3]


def test_connectivity_blocks_basin_until_ridge_overtopped():
    z = _basin()
    sea = inundation.sea_mask(z)
    f = inundation.flood_label(z, sea, 1.5)
    assert f[:, 2].all()
    assert not f[3, 5]                           # basin protected by the 2.0 m ridge
    assert inundation.flood_naive(z, sea, 1.5)[3, 5]   # naive bathtub floods it
    assert inundation.flood_label(z, sea, 2.01)[3, 5]
    assert not inundation.flood_label(z, sea, 2.0)[3, 4]   # strict '<'


def test_fill_level_equals_minimax_path():
    z = _basin()
    sea = inundation.sea_mask(z)
    fill = inundation.fill_level(z, sea)
    assert fill[3, 2] == 1.0
    assert fill[3, 5] == 2.0        # lowest way in is over the 2.0 m ridge
    assert np.isneginf(fill[0, 0])


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("eight", [False, True])
def test_fill_level_matches_label_flood_on_random_terrain(seed, eight):
    rng = np.random.default_rng(seed)
    z = rng.uniform(0.1, 3.0, size=(40, 50))
    z[:, :3] = 0.0
    sea = inundation.sea_mask(z)
    fill = inundation.fill_level(z, sea, eight=eight)
    for wl in (0.5, 1.0, 1.7, 2.3):
        a = inundation.flooded_from_fill(fill, sea, wl)
        b = inundation.flood_label(z, sea, wl, eight=eight)
        assert (a == b).all()


def test_fill_cap_marks_high_cells_inf():
    z = _basin()
    sea = inundation.sea_mask(z)
    fill = inundation.fill_level(z, sea, cap=1.5)
    assert fill[3, 2] == 1.0
    assert np.isinf(fill[3, 5]) and fill[3, 5] > 0


def test_eight_connectivity_leaks_through_diagonal():
    z = np.full((5, 5), 5.0)
    z[:, 0] = 0.0
    z[1, 1] = 0.5
    z[2, 2] = 0.5            # only diagonally connected to (1,1)
    sea = inundation.sea_mask(z)
    assert not inundation.flood_label(z, sea, 1.0)[2, 2]
    assert inundation.flood_label(z, sea, 1.0, eight=True)[2, 2]


def test_depth():
    z = np.array([[0.2, 1.0]])
    d = inundation.depth(z, np.array([[True, False]]), 0.7)
    assert d[0, 0] == pytest.approx(0.5) and d[0, 1] == 0.0


# ------------------------------------------------------------------ dem / datum

def test_read_dtu_mdt_orientation(tmp_path):
    # GRAVSOFT: header lat1 lat2 lon1 lon2 dlat dlon, rows from north to south
    p = tmp_path / "m.grd.gz"
    vals = np.array([[3.0, 4.0, 5.0], [1.0, 1.1, 1.2], [0.0, 0.1, 0.2]])   # north row first
    with gzip.open(p, "wt") as f:
        f.write("0.0 2.0 103.0 105.0 1.0 1.0\n")
        f.write("\n".join(" ".join(map(str, r)) for r in vals))
    lats, lons, g = dem.read_dtu_mdt(p)
    assert lats[0] == 2.0 and lats[-1] == 0.0
    assert g[1, 1] == pytest.approx(1.1)
    off = dem.mdt_offset(p, box=(103.5, 0.5, 105.0, 1.5))
    assert off["mean"] == pytest.approx(1.15) and off["n"] == 2


def test_cell_area_near_equator():
    tr = Affine(1 / 3600, 0, 103.6, 0, -1 / 3600, 1.48)
    a = dem.cell_area_m2(tr, 3)
    assert a[0] == pytest.approx(30.87 ** 2, rel=0.01)


# ------------------------------------------------------------------ AR6 extraction

def test_slr_read_location(tmp_path):
    h5py = pytest.importorskip("h5py")
    p = tmp_path / "t.nc"
    q = np.round(np.linspace(0, 1, 101), 2)
    q = np.unique(np.concatenate([q, [0.17, 0.83]]))
    years = np.array([2020, 2050, 2100])
    locs = np.array([12, 1746, 99])
    data = np.zeros((len(q), len(years), len(locs)), dtype=np.int16)
    for i, qq in enumerate(q):
        data[i, :, 1] = (np.array([40, 200, 700]) * (0.5 + qq)).astype(np.int16)
        data[i, :, 0] = 9999
    with h5py.File(p, "w") as f:
        f["locations"] = locs
        f["quantiles"] = q
        f["years"] = years
        f["sea_level_change"] = data
        f["lat"] = np.array([0, 1.27, 0], np.float32)
        f["lon"] = np.array([0, 103.85, 0], np.float32)
    out = slr.read_location(p)
    assert slr.value(out, 2100, "p50") == pytest.approx(0.7)
    assert slr.value(out, 2050, "p83") == pytest.approx(0.266)
    assert slr.value(out, 2020, "p05") == pytest.approx(0.022)


# ------------------------------------------------------------------ exposure

def test_rowcol_and_sample():
    tr = Affine(0.1, 0, 100.0, 0, -0.1, 2.0)
    r, c = exposure.rowcol(tr, [100.05, 100.25, 99.0], [1.95, 1.75, 1.0])
    assert list(r[:2]) == [0, 2] and list(c[:2]) == [0, 2]
    m = np.zeros((5, 5), bool)
    m[2, 2] = True
    assert list(exposure.sample(m, r, c)) == [False, True, False]


def test_near_radius():
    m = np.zeros((9, 9), bool)
    m[4, 4] = True
    n = exposure.near(m)
    assert n[4, 6] and n[2, 2] and not n[4, 7] and not n[1, 4]
    assert (exposure.near(m, 0) == m).all()


def test_subzone_population_parsing(tmp_path):
    p = tmp_path / "pop.csv"
    p.write_text("Number,Total_Total,Total_0_4\nTotal,300,1\nAng Mo Kio - Total,250,1\n"
                 "Ang Mo Kio Town Centre,200,1\nCheng San,\"50\",1\nChangi- Total,7,0\nWestern Islands,-,0\n")
    d = exposure.subzone_population(p)
    assert list(d["key"]) == ["ANGMOKIOTOWNCENTRE", "CHENGSAN", "WESTERNISLANDS"]
    assert list(d["pop"]) == [200, 50, 0]


def test_road_samples_length():
    import json
    import tempfile
    from pathlib import Path

    way = {"elements": [{"type": "way", "geometry": [{"lon": 103.8, "lat": 1.3}, {"lon": 103.81, "lat": 1.3}]}]}
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "r.json"
        p.write_text(json.dumps(way))
        rs = exposure.load_road_samples(p)
    assert rs["km"].sum() == pytest.approx(1.113, rel=0.005)
    assert len(rs) >= 100


def test_exposure_counts_and_population_share():
    shape = (6, 6)
    fl = np.zeros(shape, bool)
    fl[0, 0] = True
    sg = np.ones(shape, bool)
    a = exposure.Assets(
        b_r=np.array([0, 2, 5, 5]), b_c=np.array([2, 0, 5, 4]), b_zone=np.array([0, 0, 1, 1]),
        r_r=np.array([0, 0, 3]), r_c=np.array([0, 1, 3]), r_km=np.array([0.5, 0.25, 9.0]),
        e_r=np.array([1, 5]), e_c=np.array([1, 5]), e_station=np.array(["A", "B"]),
        zone_pop=np.array([1000.0, 500.0]), zone_pa=np.array(["P", "Q"]))
    e = exposure.exposure(fl, a, np.full(shape[0], 1e6), sg)
    assert e["area_km2"] == pytest.approx(1.0)
    assert e["buildings"] == 2                      # both zone-0 buildings within 2 cells
    assert e["population"] == pytest.approx(1000.0)
    assert e["road_km"] == pytest.approx(0.5)       # only the sample on the flooded cell
    assert e["mrt_stations"] == 1 and e["station_names"] == ["A"]
    e2 = exposure.exposure(fl, a, np.full(shape[0], 1e6), np.zeros(shape, bool))
    assert e2["area_km2"] == 0 and e2["buildings"] == 0


def test_split_sea_reservoir_inside_land_polygon_is_inland():
    dem = np.full((10, 12), 5.0)
    dem[:, :3] = 0.0                 # open sea
    dem[3:7, 6:10] = 0.0             # reservoir...
    dem[5, 3:6] = 0.0                # ...joined to the sea by a zero-valued channel (barrage gap)
    land = np.zeros_like(dem, bool)
    land[:, 3:] = True
    sea, inland = inundation.split_sea(dem, land, near_cells=1)
    assert sea[:, :3].all()
    assert sea[5, 3]                 # within 1 cell of the outside sea: shoreline mismatch
    assert not sea[4, 8] and inland[4, 8]
    assert not inland[0, 0]
    # with the plain edge-connected mask the reservoir would wrongly count as sea
    assert inundation.sea_mask(dem)[4, 8]


def test_flat_plateaus_detects_only_large_exactly_flat_patches():
    rng = np.random.default_rng(0)
    dem = rng.uniform(2, 3, size=(40, 40))
    dem[5:25, 5:25] = 1.0            # 400 identical cells -> water body
    dem[30:35, 30:35] = 1.0          # 25 identical cells -> too small
    within = np.ones_like(dem, bool)
    m = inundation.flat_plateaus(dem, within, min_cells=300)
    assert m[10, 10] and not m[32, 32] and not m[0, 0]
    assert not inundation.flat_plateaus(dem, np.zeros_like(within), min_cells=300).any()
