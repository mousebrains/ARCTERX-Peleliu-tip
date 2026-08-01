# ARCTERX 2023 — Peleliu tip vortex

Analysis of an island-wake eddy shed from the Peleliu tip (Palau), from four
GPS wave-buoy drifters and a twelve-gauge bottom-pressure array.

**The repository is self-contained.** Clone it, install four Python packages,
and every result in the two analysis documents reproduces with no external
data volume mounted.

---

## Start here

```bash
pip install -r requirements.txt

# 1. the parser reproduces the delivered CORDC files bit-for-bit
python3 tests/test_mwb_dat.py

# 2. the vorticity analysis (~30 s)
python3 src/eddy_analysis.py --output eddy_out

# 3. figures
python3 src/eddy_summary_figure.py --input eddy_out
python3 src/eddy_frame_scatter.py  --input eddy_out
```

Then read **`DRIFTER_ANALYSIS.md`** and **`PRESSURE_ANALYSIS.md`**. Both open
with the traps rather than the results, because most of the wrong answers here
looked like successes while they were wrong.

## What was found

**The eddy.** Vorticity **−1.19 × 10⁻³ s⁻¹**, Rossby number **−67**,
anticyclonic, rotation period 2.9 h, coherent for the entire 25.4 h record
(Okubo–Weiss < 0 in 100 % of windows). Divergence under 2 % of |ζ| by three
independent estimators. The vortex is **not solid-body** — vorticity falls off
with radius, core radius ~1.1–1.2 km.

**The array.** Tides resolved to 99 % of variance, M2 = 0.502 ± 0.004 m across
all twelve gauges with a 2.3° phase spread over 11 km. A co-tidal chart
showing M2 propagating 314° at 11 m/s — a shallow-water speed set by the bank
tops, not the 1500 m channel.

**Two things that do not work**, established with numbers so nobody repeats
them: drifter GPS height cannot measure sea level (off by ~100×), and bottom
pressure cannot recover the current here (biased 2–7×, because 83 % of the
flow is not barotropic).

## Layout

```
src/          analysis modules; run any of them with --help
matlab/       MATLAB port of the vorticity analysis, plus the original approach
tests/        parser regression tests against the delivered CORDC files
tools/        make_data_subset.py — rebuilds data/ from the SeaChest archive
data/         the vendored subset; see data/README.md for provenance
papers/       citations, DOIs and BibTeX (no PDFs — see LICENSE-DATA)
figures/      a handful that explain something; the rest are regenerable
```

### The modules

| module | what it does |
|---|---|
| `mwb_dat.py` | decodes the raw 85-byte instrument records |
| `mwb_nc.py` | writes CF-1.13 netCDF from them; recovers ~38 % more data than delivered |
| `pressure_array.py` | loads the pressure gauges, harmonic analysis |
| `eddy_kinematics.py` | vorticity, divergence, strain; three cross-checked estimators |
| `eddy_analysis.py` | driver — runs the analysis and writes figures |
| `paths.py` | resolves the vendored subset vs the full archive |

## Working with the full archive

The vendored pressure records are 1-minute block means, which is everything
the tidal and residual work needs. For finer-than-1-minute analysis at gauges
other than the two full-rate segments, point at the archive:

```bash
export ARCTERX_ARCHIVE=/Volumes/SeaChest/ARCTERX/2023/Wake
```

Without it, `pressure_array.load` says so rather than failing obscurely.

## MATLAB

`matlab/eddy_kinematics_drifters.m` reproduces the Python analysis and adds
bootstrap confidence intervals on the vortex-structure fit. Needs R2019b+ with
the Mapping, Statistics and Curve Fitting toolboxes. Point `cfg.dataDir` at
`data/drifters`.

`matlab/explore_drifter_paths.m` is the earlier geometric approach, kept
because `DRIFTER_ANALYSIS.md` §3 explains what it does and why the velocity
method replaced it.

## Related work

The X-band radar pipeline lives in a separate repository (`wamos_tpw`), which
also carries the tidal-harmonic module (`wamos_tpw.tides`) used here for the
constituent fits and its Rayleigh/conditioning guards.

## Licences

Code is MIT (`LICENSE`). Data is **not uniformly licensed** — read
`LICENSE-DATA` before redistributing anything under `data/`. Journal PDFs are
deliberately not included.
