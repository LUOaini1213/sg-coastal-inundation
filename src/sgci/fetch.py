"""Download open input datasets into ``data_cache/`` (git-ignored).

Every source here is public and needs no account or API key. Licences are
listed in the README. Re-running is cheap: files that already exist are kept.
"""
from __future__ import annotations

import io
import json
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data_cache"
UA = {"User-Agent": "sg-coastal-inundation research script (github.com/LUOaini1213)"}

UHSLC_HOURLY = "https://uhslc.soest.hawaii.edu/data/csv/rqds/pacific/hourly/h699a.csv"
UHSLC_META = "https://uhslc.soest.hawaii.edu/rqds/pacific/doc/qa699a.dmt"
AR6_ZIP = "https://zenodo.org/api/records/6382554/files/ar6-regional-confidence.zip/content"
AR6_LOCS = "https://zenodo.org/api/records/6382554/files/location_list.lst/content"
AR6_MEMBER = ("ar6-regional-confidence/regional/confidence_output_files/medium_confidence/"
              "{s}/total_{s}_medium_confidence_values.nc")
COP_DEM = ("https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{t}_DEM/"
           "Copernicus_DSM_COG_10_{t}_DEM.tif")
DEM_TILES = ["N01_00_E103_00", "N01_00_E104_00"]
OVERPASS = "https://overpass-api.de/api/interpreter"
DATAGOV_POLL = "https://api-open.data.gov.sg/v1/public/api/datasets/{d}/poll-download"
SCENARIOS = ["ssp126", "ssp245", "ssp585"]


def _get(url: str, headers: dict | None = None, data: bytes | None = None, timeout=300) -> bytes:
    h = dict(UA)
    h.update(headers or {})
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers=h, data=data)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # network hiccups: retry with backoff
            if attempt == 3:
                raise
            print(f"  retry {url[:80]} ({exc})")
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("unreachable")


def _save(path: Path, url: str, **kw) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"download {url[:100]}")
    blob = _get(url, **kw)
    path.write_bytes(blob)
    return path


class _HttpRange(io.RawIOBase):
    """Seekable read-only file over HTTP range requests (to pull one member out of a 9 GB zip)."""

    def __init__(self, url: str):
        self.url = url
        req = urllib.request.Request(url, method="HEAD", headers=UA)
        with urllib.request.urlopen(req, timeout=120) as r:
            self.size = int(r.headers["Content-Length"])
        self.pos = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        d = _get(self.url, headers={"Range": f"bytes={self.pos}-{self.pos + n - 1}"})
        b[: len(d)] = d
        self.pos += len(d)
        return len(d)


def fetch_tide_gauge() -> Path:
    _save(CACHE / "uhslc" / "qa699a.dmt", UHSLC_META)
    return _save(CACHE / "uhslc" / "h699a.csv", UHSLC_HOURLY)


def fetch_ar6() -> list[Path]:
    _save(CACHE / "ar6" / "location_list.lst", AR6_LOCS)
    out = []
    z = None
    for s in SCENARIOS:
        p = CACHE / "ar6" / f"total_{s}_medium_confidence_values.nc"
        if not p.exists():
            if z is None:
                z = zipfile.ZipFile(io.BufferedReader(_HttpRange(AR6_ZIP), buffer_size=8 << 20))
            print(f"extract {s} from remote AR6 zip (~72 MB)")
            p.write_bytes(z.read(AR6_MEMBER.format(s=s)))
        out.append(p)
    return out


def fetch_dem() -> list[Path]:
    return [_save(CACHE / "dem" / f"{t}.tif", COP_DEM.format(t=t)) for t in DEM_TILES]


def fetch_datagov(dataset_id: str, name: str) -> Path:
    p = CACHE / "datagov" / name
    if p.exists():
        return p
    meta = json.loads(_get(DATAGOV_POLL.format(d=dataset_id)))
    url = meta["data"]["url"]
    return _save(p, url)


OSM_QUERIES = {
    "buildings": """[out:json][timeout:600];
area["ISO3166-1"="SG"][admin_level=2]->.sg;
(way["building"](area.sg);relation["building"](area.sg););
out center tags;""",
    "roads": """[out:json][timeout:600];
area["ISO3166-1"="SG"][admin_level=2]->.sg;
way["highway"~"^(motorway|trunk|primary|secondary|tertiary|residential|unclassified|motorway_link|trunk_link|primary_link|secondary_link|tertiary_link|service|living_street)$"](area.sg);
out geom tags;""",
}


def fetch_osm(name: str) -> Path:
    p = CACHE / "osm" / f"{name}.json"
    if p.exists():
        return p
    body = urllib.parse.urlencode({"data": OSM_QUERIES[name]}).encode()
    return _save(p, OVERPASS, data=body, timeout=900)


# data.gov.sg dataset ids (Singapore Open Data Licence v1.0)
DATAGOV = {
    "planning_areas.geojson": "d_4765db0e87b9c86336792efe8a1f7a66",  # URA MP2019 Planning Area Boundary (No Sea)
    "subzones.geojson": "d_8594ae9ff96d0c708bc2af633048edfb",        # URA MP2019 Subzone Boundary (No Sea)
    "mrt_exits.geojson": "d_b39d3a0871985372d7e1637193335da5",       # LTA MRT Station Exit (GeoJSON)
    "population_census2020.csv": "d_d95ae740c0f8961a0b10435836660ce0",  # SingStat Census 2020 residents by PA/subzone
}


MDT_URL = "https://ftp.space.dtu.dk/pub/DTU22/MDT/dtu22mdt.grd.gz"


def fetch_mdt() -> Path:
    return _save(CACHE / "mdt" / "dtu22mdt.grd.gz", MDT_URL)


def main():
    fetch_tide_gauge()
    fetch_mdt()
    fetch_dem()
    for name, did in DATAGOV.items():
        fetch_datagov(did, name)
    for name in OSM_QUERIES:
        fetch_osm(name)
    fetch_ar6()


if __name__ == "__main__":
    main()
