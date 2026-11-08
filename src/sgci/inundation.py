"""Static ("bathtub") inundation with hydraulic connectivity.

A land cell floods for a still-water level ``wl`` (m above MSL) only if a path of cells that
are all below ``wl`` connects it to the open sea. Instead of re-running a flood fill for every
scenario we compute, once, the *fill level* of every cell: the lowest still-water level at which
it becomes connected to the sea (a minimax path problem, solved with a priority flood,
Barnes et al. 2014). Then  flooded(wl) == fill_level < wl  for any wl.
``flood_label`` is the textbook connected-component version used to cross-check it in tests.
"""
from __future__ import annotations

import heapq

import numpy as np
from scipy import ndimage

FOUR = ndimage.generate_binary_structure(2, 1)
EIGHT = ndimage.generate_binary_structure(2, 2)


def sea_mask(dem: np.ndarray, sea_value: float = 0.0) -> np.ndarray:
    """Open sea = cells at the DEM's flattened sea value that connect to the raster edge.

    Copernicus GLO-30 sets ocean cells to exactly 0 m. Zero-valued cells not connected to the
    edge (e.g. reservoirs behind dams) are *not* sea.
    """
    zero = dem <= sea_value
    lab, _ = ndimage.label(zero, structure=FOUR)
    edge = np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
    edge = edge[edge > 0]
    return np.isin(lab, edge)


def flat_plateaus(dem: np.ndarray, within: np.ndarray, min_cells: int = 300) -> np.ndarray:
    """Water bodies flattened by DEM editing: connected patches of one identical non-zero value.

    Copernicus GLO-30 flattens lakes and reservoirs to a single height. Natural or built land
    surfaces are never exactly flat over >= min_cells cells (~0.29 km^2 at 30 m).
    """
    out = np.zeros(dem.shape, dtype=bool)
    q = np.where(within & (dem != 0), np.round(dem, 3), np.nan)
    vals, counts = np.unique(q[np.isfinite(q)], return_counts=True)
    for v in vals[counts >= min_cells]:
        lab, n = ndimage.label(q == v, structure=FOUR)
        sizes = np.bincount(lab.ravel())
        big = np.flatnonzero(sizes >= min_cells)
        out |= np.isin(lab, big[big > 0])
    return out


def split_sea(dem: np.ndarray, land_poly: np.ndarray, near_cells: int = 3):
    """Open sea vs. inland water for a DEM whose water cells are flattened to 0.

    Zero-valued cells inside the land polygons (reservoirs behind barrages and dams, ponds,
    areas still water when the DEM was acquired) are *inland water*, unless they lie within
    ``near_cells`` of open sea outside the polygons (shoreline mismatch between the polygon
    and the DEM coastline). Returns (sea, inland_water).
    """
    edge_sea = sea_mask(dem)
    outside = edge_sea & ~land_poly
    st = np.ones((2 * near_cells + 1, 2 * near_cells + 1), dtype=bool)
    near_out = ndimage.binary_dilation(outside, structure=st)
    sea = outside | (edge_sea & land_poly & near_out)
    inland = (dem <= 0) & ~sea
    return sea, inland


def flood_label(z: np.ndarray, sea: np.ndarray, wl: float, eight: bool = False) -> np.ndarray:
    """Land cells flooded at still-water level wl (connected-component version)."""
    wet = sea | (z < wl)
    lab, _ = ndimage.label(wet, structure=EIGHT if eight else FOUR)
    keep = np.unique(lab[sea])
    keep = keep[keep > 0]
    return np.isin(lab, keep) & ~sea


def flood_naive(z: np.ndarray, sea: np.ndarray, wl: float) -> np.ndarray:
    """No-connectivity bathtub: every land cell below wl, used to show why connectivity matters."""
    return (z < wl) & ~sea


def fill_level(z: np.ndarray, sea: np.ndarray, cap: float = np.inf, eight: bool = False) -> np.ndarray:
    """Minimum still-water level that floods each cell via a sea-connected path (inf if > cap).

    Strictly, cell c floods at wl iff there is a path from the sea on which every land cell has
    z < wl, i.e. iff max(z along the best path) < wl. The priority flood returns that minimax
    value; the strict '<' is applied by the caller.
    """
    rows, cols = z.shape
    out = np.full(z.shape, np.inf)
    done = np.zeros(z.shape, dtype=bool)
    heap: list[tuple[float, int]] = []
    sea_idx = np.flatnonzero(sea)
    done.flat[sea_idx] = True
    out.flat[sea_idx] = -np.inf
    if eight:
        nbrs = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    else:
        nbrs = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    # seed: land cells adjacent to the sea
    pad = np.pad(sea, 1)
    adj = np.zeros_like(sea)
    for dr, dc in nbrs:
        adj |= pad[1 + dr:rows + 1 + dr, 1 + dc:cols + 1 + dc]
    adj &= ~sea
    for i in np.flatnonzero(adj):
        v = z.flat[i]
        if v <= cap:
            heapq.heappush(heap, (float(v), int(i)))
            done.flat[i] = True
    zf = z.ravel()
    while heap:
        lvl, i = heapq.heappop(heap)
        out.flat[i] = lvl
        r, c = divmod(i, cols)
        for dr, dc in nbrs:
            rr, cc = r + dr, c + dc
            if 0 <= rr < rows and 0 <= cc < cols:
                j = rr * cols + cc
                if not done.flat[j]:
                    v = max(lvl, zf[j])
                    if v <= cap:
                        done.flat[j] = True
                        heapq.heappush(heap, (v, j))
    return out


def flooded_from_fill(fill: np.ndarray, sea: np.ndarray, wl: float) -> np.ndarray:
    return (fill < wl) & ~sea


def depth(z: np.ndarray, flooded: np.ndarray, wl: float) -> np.ndarray:
    d = np.where(flooded, wl - z, 0.0)
    return np.clip(d, 0.0, None)
