"""Static figures for the README and report (outputs/figures/*.png).

Charts use committed outputs only; the two maps also need the cached grids (data_cache/).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from sgci import eva  # noqa: E402

OUT = ROOT / "outputs"
FIG = OUT / "figures"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e3e2de"
INLAND = "#c9c4b8"
SSP_COL = {"ssp126": "#2a78d6", "ssp245": "#1baf7a", "ssp585": "#eb6834"}
SSP_LABEL = {"ssp126": "SSP1-2.6", "ssp245": "SSP2-4.5", "ssp585": "SSP5-8.5"}
plt.rcParams.update({"font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150,
                     "savefig.bbox": "tight", "lines.linewidth": 2})


def fig_return_levels():
    t = json.loads((OUT / "tide.json").read_text())
    f = t["fits"]["linear"]
    sm = pd.read_csv(OUT / "season_maxima.csv")
    x = np.sort(sm["max"].to_numpy())
    n = len(x)
    p = (np.arange(1, n + 1) - 0.44) / (n + 0.12)          # Gringorten
    T_obs = 1 / (1 - p)
    g = eva.GEV(f["mu"], f["sigma"], f["xi"], 0, n)
    T = np.logspace(np.log10(1.1), 3, 200)
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(T, eva.gev_return_level(g, T), color="#2a78d6", label=f"GEV fit (ξ = {f['xi']:.2f})")
    Ts = [int(k) for k in f["return_levels"]]
    lo = [f["return_levels"][str(k)]["lo"] for k in Ts]
    hi = [f["return_levels"][str(k)]["hi"] for k in Ts]
    ax.fill_between(Ts, lo, hi, color="#2a78d6", alpha=0.15, lw=0, label="95% bootstrap CI")
    gp = t["fits"]["pot_gpd"]
    ax.plot(Ts, [gp["return_levels"][str(k)] for k in Ts], "--", color="#eb6834", lw=1.5,
            label="GPD peaks-over-threshold check")
    ax.plot(T_obs, x, "o", ms=4.5, color=INK, mfc="white", mew=1.2, label=f"Season maxima 1984–2023 (n = {n})")
    ax.set_xscale("log")
    ax.set_xlabel("Return period (years)")
    ax.set_ylabel("Still-water level (m above MSL 1995–2014)")
    ax.set_title("Tanjong Pagar extreme sea levels (detrended)", loc="left", fontsize=10.5)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    fig.savefig(FIG / "return_levels.png")
    plt.close(fig)


def fig_slr():
    a = json.loads((OUT / "ar6_tanjong_pagar.json").read_text())
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for s, d in a["scenarios"].items():
        yrs = sorted(int(y) for y in d)
        med = [d[str(y)]["p50"] for y in yrs]
        ax.fill_between(yrs, [d[str(y)]["p17"] for y in yrs], [d[str(y)]["p83"] for y in yrs],
                        color=SSP_COL[s], alpha=0.14, lw=0)
        ax.plot(yrs, med, color=SSP_COL[s], label=SSP_LABEL[s])
        ax.annotate(f"{med[-1]:.2f} m", (yrs[-1], med[-1]), xytext=(4, 0), textcoords="offset points",
                    va="center", fontsize=8.5, color=INK2)
    ax.set_xlim(2020, 2110)
    ax.set_ylabel("Relative sea-level rise (m, vs 1995–2014)")
    ax.set_title("IPCC AR6 projections at Tanjong Pagar (median, shaded = likely range)", loc="left", fontsize=10.5)
    ax.legend(frameon=False, loc="upper left")
    fig.savefig(FIG / "ar6_slr.png")
    plt.close(fig)


def fig_area_curve():
    ac = pd.read_csv(OUT / "area_curve.csv")
    sc = pd.read_csv(OUT / "exposure_scenarios.csv")
    meta = json.loads((OUT / "grid_meta.json").read_text())
    lv = ac["level_mm"] / 1000
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(lv, ac["naive_km2"], color="#8a8984", lw=1.5, label="No connectivity check (naive bathtub)")
    ax.plot(lv, ac["connected_km2"], color="#2a78d6", label=f"Connected to sea, datum offset {meta['datum_offset_m']:.2f} m")
    for d, ls in ((-0.3, ":"), (0.3, "--")):
        k = np.round((lv + d) * 1000).astype(int)
        ok = (k >= 0) & (k < len(ac))
        ax.plot(lv[ok], ac["connected_km2"].to_numpy()[k[ok]], ls, color="#2a78d6", lw=1.1,
                label=f"Connected, offset {meta['datum_offset_m'] + d:.2f} m")
    rp100 = sc[sc["id"].str.startswith("rp100")]
    ax.plot(rp100["wl"], rp100["area_km2"], "o", ms=5, color="#eb6834", label="100-yr scenarios")
    ax.set_xlim(1.0, 3.5)
    ax.set_ylim(0, max(ac["naive_km2"][lv <= 3.5].max(), 10) * 1.05)
    ax.set_xlabel("Still-water level (m above MSL 1995–2014)")
    ax.set_ylabel("Flooded land in Singapore (km²)")
    ax.set_title("Flooded area vs water level, and the datum assumption", loc="left", fontsize=10.5)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(FIG / "area_vs_level.png")
    plt.close(fig)


def fig_maps():
    cache = ROOT / "data_cache"
    sys.path.insert(0, str(ROOT / "scripts"))
    import run_pipeline as rp

    G = rp.load_grids()
    fill = np.load(cache / f"fill_{G['delta']:.2f}.npy")
    sc = {s["id"]: s for s in json.loads((OUT / "scenarios.json").read_text())}
    tr = G["tr"]
    ext = [tr.c, tr.c + tr.a * fill.shape[1], tr.f + tr.e * fill.shape[0], tr.f]
    fig, ax = plt.subplots(figsize=(8, 5.2))
    land = np.where(G["sea"], np.nan, 1.0)
    ax.imshow(land, extent=ext, cmap=matplotlib.colors.ListedColormap(["#ecebe7"]), interpolation="nearest")
    water = np.where(G["inland_zero"], 1.0, np.nan)
    ax.imshow(water, extent=ext, cmap=matplotlib.colors.ListedColormap([INLAND]), interpolation="nearest")
    ax.plot([], [], "s", color=INLAND, ms=8, label="Inland water (reservoirs; held at regulated level)")
    order = [("rp100_ssp585_2100", "#86b6ef"), ("rp100_ssp245_2100", "#3987e5"), ("rp100_2020", "#184f95")]
    for sid, col in order:
        m = (fill < sc[sid]["wl"]) & G["sg"]
        ax.imshow(np.where(m, 1.0, np.nan), extent=ext, cmap=matplotlib.colors.ListedColormap([col]),
                  interpolation="nearest")
        ax.plot([], [], "s", color=col, ms=8, label=f"{sc[sid]['label']} ({sc[sid]['wl']:.2f} m)")
    G["pa"].boundary.plot(ax=ax, color="#7d7c78", lw=0.3)
    m = json.loads((OUT / "inertial.json").read_text())
    b = m["box_lonlat"]
    ax.plot([b[0], b[2], b[2], b[0], b[0]], [b[1], b[1], b[3], b[3], b[1]], color="#eb6834", lw=1.2,
            label="2-D model sub-domain")
    ax.set_xlim(103.6, 104.09)
    ax.set_ylim(1.16, 1.475)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.legend(frameon=True, fontsize=8, loc="lower left", facecolor="white", edgecolor=GRID)
    ax.set_title("Connected still-water flooding, 100-yr event (screening; not an official flood map)",
                 loc="left", fontsize=10.5)
    fig.savefig(FIG / "flood_map.png")
    plt.close(fig)

    # model vs bathtub in the sub-domain
    hmax = np.load(cache / "model_hmax.npy")
    r0, c0, n = m["box_rows_cols"]
    sl = (slice(r0, r0 + n), slice(c0, c0 + n))
    bt = (fill[sl] < m["wl"]) & (m["wl"] - G["z"][sl] >= m["depth_min_m"]) & G["sg"][sl]
    mb = (hmax >= m["depth_min_m"]) & G["sg"][sl]
    fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.6), layout="constrained")
    sub_ext = [tr.c + c0 * tr.a, tr.c + (c0 + n) * tr.a, tr.f + (r0 + n) * tr.e, tr.f + r0 * tr.e]
    base = np.where(G["sea"][sl], np.nan, 1.0)
    for a in axs:
        a.imshow(base, extent=sub_ext, cmap=matplotlib.colors.ListedColormap(["#ecebe7"]))
        a.imshow(np.where(G["inland_zero"][sl], 1.0, np.nan), extent=sub_ext,
                 cmap=matplotlib.colors.ListedColormap([INLAND]))
        a.ticklabel_format(useOffset=False)
        a.set_aspect("equal")
        a.grid(False)
        a.tick_params(labelsize=7)
    d = np.where(mb, hmax, np.nan)
    im = axs[0].imshow(d, extent=sub_ext, cmap="Blues", vmin=0, vmax=2.0)
    fig.colorbar(im, ax=axs[0], fraction=0.046, pad=0.02, label="Max depth (m)")
    axs[0].set_title(f"2-D local-inertial model: {m['runs']['base']['area_km2']:.2f} km²", fontsize=9.5, loc="left")
    cls = np.full(bt.shape, np.nan)
    cls[bt & mb] = 0
    cls[bt & ~mb] = 1
    cls[mb & ~bt] = 2
    axs[1].imshow(cls, extent=sub_ext, cmap=matplotlib.colors.ListedColormap(["#1c5cab", "#eb6834", "#1baf7a"]),
                  vmin=0, vmax=2, interpolation="nearest")
    for lab, col in (("both", "#1c5cab"), ("bathtub only", "#eb6834"), ("model only", "#1baf7a")):
        axs[1].plot([], [], "s", color=col, label=lab)
    axs[1].legend(fontsize=7.5, loc="lower left")
    axs[1].set_title(f"Bathtub: {m['bathtub_area_km2']:.2f} km²", fontsize=9.5, loc="left")
    fig.suptitle("100-yr storm tide + SSP5-8.5 2100 median, observed 1999-12-23 hydrograph shape", fontsize=10, x=0.02,
                 ha="left")
    fig.savefig(FIG / "model_vs_bathtub.png")
    plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    fig_return_levels()
    fig_slr()
    fig_area_curve()
    if "--charts-only" not in sys.argv:
        fig_maps()


if __name__ == "__main__":
    main()
