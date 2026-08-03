# ARCTERX 2023 — Peleliu tip vortex

Analysis of an island-wake eddy shed from the Peleliu tip (Palau), from four
GPS wave-buoy drifters and a twelve-gauge bottom-pressure array.

**The repository is self-contained** in its data: clone it, install four Python
packages, and the analyses run with no external volume mounted.

**Reproduction is not yet complete, and one number is contested.** The drifter
half of `DRIFTER_ANALYSIS.md` reproduces end to end. The tidal results in
`PRESSURE_ANALYSIS.md` reproduce via `src/pressure_analysis.py`; its residual
spectra, noise-floor split and gradient-inversion sections do not yet have a
driver. That script also reports a **contradiction in the published M2 phase
gradient** — see the flag in `PRESSURE_ANALYSIS.md` §8.1.

---

## Start here

```bash
pip install -r requirements.txt

# 1. the vorticity analysis (~30 s)
python3 src/eddy_analysis.py --output eddy_out

# 2. figures
python3 src/eddy_summary_figure.py --input eddy_out
python3 src/eddy_frame_scatter.py  --input eddy_out

# 3. round-trip: the regenerated NetCDF matches the raw bursts
python3 tests/test_mwb_nc.py data/drifters/mwb458d02_gps_timeseries.nc data/raw_dat/458

# 4. null tests: does the pipeline recover a vortex whose answer we know?
python3 tests/test_synthetic_recovery.py

# 5. the tidal half: harmonic fits, co-tidal chart, consistency check
python3 src/pressure_analysis.py
```

`tests/test_mwb_dat.py` additionally checks the parser against the
CORDC-delivered NetCDF bit-for-bit. Those files are not redistributed here, so
it skips unless you point it at your own copy:

```bash
python3 tests/test_mwb_dat.py data/raw_dat/458 /path/to/mwb458d02_gps_timeseries.nc
```

Then read **`DRIFTER_ANALYSIS.md`** and **`PRESSURE_ANALYSIS.md`**. Both open
with the traps rather than the results, because most of the wrong answers here
looked like successes while they were wrong.

**`docs/`** derives the mathematics behind every method — the velocity-gradient
tensor, Stokes' theorem, Lamb-Oseen, harmonic analysis, the pressure response
function — and says where each one breaks. Start at `docs/README.md`.

## What was found

**The eddy.** Vorticity **−1.19 × 10⁻³ s⁻¹**, Rossby number **−67**,
anticyclonic, rotation period 2.9 h, coherent for the entire 25.4 h record
(Okubo–Weiss < 0 in 100 % of windows). Divergence is 2.3 % of |ζ| by the
contour estimator and 0.1–0.2 % by the other two. The vortex is **not solid-body** — vorticity falls off
with radius, core radius ~1.1–1.2 km.

**The array.** Tides resolved to 99 % of variance, M2 = 0.502 ± 0.004 m across
all twelve gauges with a 2.3° phase spread over 11 km. A co-tidal chart showing
M2 propagating **297° at 35 m/s**, resolved at only 3.2σ — a speed between the
bank-top and channel shallow-water limits, matching neither. An earlier
**314° at 11 m/s** is withdrawn; see `PRESSURE_ANALYSIS.md` §8.1.

**Two things that do not work**, established with numbers so nobody repeats
them: drifter GPS height cannot measure sea level (off by ~100×), and bottom
pressure cannot recover the current here (biased 2–7×, because 83 % of the
flow is not barotropic).

## Layout

```
src/          analysis modules; run any of them with --help
docs/         the mathematics and reasoning behind every method
matlab/       MATLAB port of the vorticity analysis, plus the original approach
tests/        parser, round-trip and synthetic-recovery tests
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
| `pressure_analysis.py` | driver — regenerates the tidal numbers and checks them |
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

The vorticity analysis exists in both languages. Paths resolve relative to the
script, so both run from a clone with nothing to edit:

```matlab
run matlab/eddy_kinematics_drifters.m    % the analysis, ~1 min
run matlab/explore_drifter_paths.m       % the earlier geometric approach
```

| script | toolboxes beyond base MATLAB |
|---|---|
| `eddy_kinematics_drifters.m` | Mapping, Statistics and Machine Learning, Curve Fitting |
| `explore_drifter_paths.m` | Mapping, Parallel Computing (`parfor`; runs serially without it) |

Written against R2026a. Output goes to `eddy_out_matlab/` (git-ignored).

### Which one is authoritative

**`src/` is the reference implementation; the MATLAB is a port.** Two
implementations of one method will drift apart, so when they disagree the
Python is right by definition and the MATLAB has the bug — that convention is
the whole defense against silently diverging answers.

Both were run from a clean path on 2026-08-01. Every median statistic agrees
to the digits printed:

| | Python | MATLAB |
|---|---|---|
| ζ median | −1.187 × 10⁻³ s⁻¹ | −1.187 × 10⁻³ s⁻¹ |
| Rossby | −66.9 | −66.9 |
| rotation period | 2.94 h | 2.94 h |
| leave-one-out 1σ | 1.622 × 10⁻⁴ (16 %) | 1.623 × 10⁻⁴ (16 %) |
| constellation revolutions | −6.86 | −6.86 |
| Okubo–Weiss < 0 | 100 % of windows | 100 % of windows |
| **Lamb–Oseen R** | **1206 m** | **1116 m**, 95 % CI [1005, 1605] |
| **circulation Γ** | **−6114 m² s⁻¹** | **−5439**, 95 % CI [−9051, −4591] |

**The last two rows disagree, and that is expected** — it is the single
sensitive step in the analysis, not a porting bug. The two implementations
differ in percentile convention and moving-average edge handling, which is
enough to move R by ~90 m; a grid over bin counts and radius cutoffs spans
1161–1286 m. `DRIFTER_ANALYSIS.md` §5 documents this, and it is why the
MATLAB adds a 1000-replicate bootstrap the Python does not have. Each
implementation's point estimate falls inside the other's interval.

So: quote the interval for R and Γ, never a point estimate. Treat a
disagreement in the *other* rows as a regression — those are medians over
hundreds of windows and should match exactly.

If you change the method, change `src/` first and re-run both.

### What has no MATLAB equivalent

**The entire pressure-array analysis** — `src/pressure_array.py` and
everything in `PRESSURE_ANALYSIS.md`: harmonic fits, the co-tidal chart, the
C05 rotation, the instrument noise floor. A MATLAB-only user gets the drifter
half of the project and needs Python for the other half. Porting it is not
planned; the tidal fitting leans on `wamos_tpw.tides` and its Rayleigh and
conditioning guards, which have no MATLAB counterpart.

## Related work

The X-band radar pipeline lives in a separate repository (`wamos_tpw`), which
also carries the tidal-harmonic module (`wamos_tpw.tides`) used here for the
constituent fits and its Rayleigh/conditioning guards.

## Licenses

Code is MIT (`LICENSE`). Data is **not uniformly licensed** — read
`LICENSE-DATA` before redistributing anything under `data/`. Journal PDFs are
deliberately not included.
