"""Mutation check: plant known bugs one at a time; each must make the test suite fail.

A green test suite is only evidence if the tests can go red. Every mutant below is a
plausible bug (sign error, off-by-one, wrong exponent, dropped mask...). The script applies
one mutant, runs the relevant test module, restores the file, and fails if any mutant
survives (tests still pass) or if a mutant's anchor text is not found exactly once.

    python tools/mutation_check.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S = "src/sgci/"
TE, TS, TI = "tests/test_tide_eva.py", "tests/test_spatial.py", "tests/test_inertial.py"

MUTANTS = [
    # (file, original, mutated, tests, description)
    (S + "tide.py", "MISSING = -32767", "MISSING = -32768", TE, "wrong missing-value flag"),
    (S + "tide.py", "year % 400 == 0", "year % 400 == 1", TE, "leap-year rule"),
    (S + "tide.py", "idx.month >= 7", "idx.month >= 6", TE, "season starts in June"),
    (S + "tide.py", "+ 0.5 - BASELINE_MID   # annual", "- BASELINE_MID   # annual", TE, "half-year trend offset"),
    (S + "tide.py", "return rel - trend.slope_m_per_yr", "return rel + trend.slope_m_per_yr", TE, "detrend sign"),
    (S + "tide.py", "& (am.index <= period[1]) & am[\"valid\"]]", "& (am.index <= period[1])]", TE,
     "baseline uses incomplete years"),
    (S + "eva.py", "y ** (-g.xi) - 1.0", "y ** (g.xi) - 1.0", TE, "GEV return-level exponent sign"),
    (S + "eva.py", "(1.0 + 1.0 / xi) * np.sum(np.log(t)) + np.sum(", "(1.0 - 1.0 / xi) * np.sum(np.log(t)) + np.sum(",
     TE, "GEV likelihood term"),
    (S + "eva.py", "return g.mu - g.sigma * np.log(y)", "return g.mu + g.sigma * np.log(y)", TE, "Gumbel return level"),
    (S + "eva.py", "np.percentile(rl, ci[0], axis=0)", "np.percentile(rl, ci[1], axis=0)", TE, "bootstrap CI bound"),
    (S + "eva.py", "times_h[cur[-1]] > run_h", "times_h[cur[-1]] > 2 * run_h", TE, "declustering run length"),
    (S + "eva.py", "m = g.rate * T", "m = T", TE, "GPD ignores exceedance rate"),
    (S + "inundation.py", "edge = edge[edge > 0]", "edge = edge[edge >= 0]", TS, "background label counted as sea"),
    (S + "inundation.py", "v = max(lvl, zf[j])", "v = zf[j]", TS, "fill level not minimax"),
    (S + "inundation.py", "wet = sea | (z < wl)", "wet = sea | (z <= wl)", TS, "non-strict flooding threshold"),
    (S + "inundation.py", "sea = outside | (edge_sea & land_poly & near_out)", "sea = edge_sea", TS,
     "reservoirs treated as open sea"),
    (S + "inundation.py", "big = np.flatnonzero(sizes >= min_cells)", "big = np.flatnonzero(sizes >= 1)", TS,
     "plateau size threshold dropped"),
    (S + "exposure.py", "st = np.ones((2 * cells + 1, 2 * cells + 1)", "st = np.ones((2 * cells - 1, 2 * cells - 1)",
     TS, "proximity radius off by one"),
    (S + "exposure.py", "c = np.floor((lon - transform.c) / transform.a)", "c = np.floor((lon - transform.c) / transform.e)",
     TS, "column from wrong pixel size"),
    (S + "exposure.py", "pop = a.zone_pop * share", "pop = a.zone_pop * share * 0.5", TS, "population share"),
    (S + "exposure.py", 'r"-\\s*Total\\s*$"', 'r" - Total$"', TS, "planning-area total rows kept"),
    (S + "slr.py", "slc[i] / 1000.0", "slc[i] / 100.0", TS, "mm to m conversion"),
    (S + "dem.py", "lats = la2 - np.arange(nla) * dla", "lats = la1 + np.arange(nla) * dla", TS, "MDT grid flipped"),
    (S + "inertial.py", "np.asarray(n_manning, dtype=float) ** 2", "np.asarray(n_manning, dtype=float) ** 1", TI,
     "Manning n not squared"),
    (S + "inertial.py", "dh[:, 1:] += qx", "dh[:, 1:] -= qx", TI, "continuity sign"),
    (S + "inertial.py", "hf = np.maximum(eta[:, 1:], eta[:, :-1]) - zx",
     "hf = np.maximum(eta[:, 1:], eta[:, :-1]) - np.minimum(z[:, 1:], z[:, :-1])", TI, "face flow depth uses min bed"),
    (S + "inertial.py", "(3.0 / 7.0)", "(3.0 / 5.0)", TI, "analytical solution exponent"),
    (S + "inertial.py", "alpha * dx / np.sqrt", "4 * alpha * dx / np.sqrt", TI, "time step violates CFL"),
    (S + "inertial.py", "0.5 * (1.0 - theta) * (q[:, :-2] + q[:, 2:])", "0.0 * (q[:, :-2] + q[:, 2:])", TI,
     "theta blending drops neighbour discharges"),
]


def run_tests(test: str) -> bool:
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", test],
                       cwd=ROOT, capture_output=True, text=True, timeout=900)
    return r.returncode == 0


def main():
    for t in (TE, TS, TI):
        if not run_tests(t):
            print(f"unmutated tests fail in {t}; fix them first")
            return 1
    survivors, broken = [], []
    for i, (path, old, new, test, desc) in enumerate(MUTANTS, 1):
        p = ROOT / path
        src = p.read_text(encoding="utf8")
        if src.count(old) != 1:
            broken.append(f"{path}: anchor found {src.count(old)}x: {old!r}")
            continue
        try:
            p.write_text(src.replace(old, new), encoding="utf8")
            passed = run_tests(test)
        finally:
            p.write_text(src, encoding="utf8")
        print(f"[{i:2d}/{len(MUTANTS)}] {'SURVIVED' if passed else 'killed  '}  {desc}  ({path})", flush=True)
        if passed:
            survivors.append(desc)
    print(f"\n{len(MUTANTS) - len(survivors) - len(broken)}/{len(MUTANTS)} mutants killed")
    for b in broken:
        print("BROKEN ANCHOR:", b)
    for s in survivors:
        print("SURVIVOR:", s)
    return 1 if (survivors or broken or len(MUTANTS) < 15) else 0


if __name__ == "__main__":
    sys.exit(main())
