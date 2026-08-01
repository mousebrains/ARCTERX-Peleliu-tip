# data/ — provenance and what was trimmed

Everything here is derived from the ARCTERX 2023 Wake archive
(`/Volumes/SeaChest/ARCTERX/2023/Wake`, ~2.6 GB for the parts used). This
directory is the **self-contained subset**: every analysis in
`PRESSURE_ANALYSIS.md` and `DRIFTER_ANALYSIS.md` runs from it with no external
volume mounted.

Rebuild it with `python3 tools/make_data_subset.py --archive <path>`. That
script *is* the provenance record — read it if you need to know exactly what
was kept.

Licensing is **not** uniform here. Read `../LICENSE-DATA` before
redistributing: `bathy/` is Coral Reef Research Foundation (CRRF) data.

---

## What is here

| path | size | content |
|---|---|---|
| `raw_dat/458/` | 19 MB | 117 raw `.dat` bursts, mwb458 deployment 2 |
| `drifters/` | 28 MB | 4 regenerated CF-1.13 timeseries |
| `pressure/1min/` | 5 MB | 12 gauges, 1-minute block means |
| `pressure/highrate/` | 1.6 MB | 2 short full-rate segments |
| `adcp/` | 268 KB | C05 depth-averaged current + the event list |
| `bathy/` | 1.7 MB | 25 m DEM, cropped, Z only, float32 |
| `eddy_kinematics.npz` | 100 KB | output of `eddy_analysis.py`, for the figure scripts |

## What was trimmed, and why it costs nothing

**Pressure: 2.35 GB → 5 MB.** The gauges run at 1–16 Hz for a month. Every
tidal, residual, co-tidal and gradient result comes from **1-minute block
means**, so that is what ships — a ~180× reduction losing nothing above the
wave band.

Decimation is by **block average, never striding**. Striding aliases swell
into the tidal band: it inflates the post-fit residual from ~4 cm to 7.5–14 cm
and produces a physically impossible 149° M2 phase spread across the array.
See `PRESSURE_ANALYSIS.md` §3.2.

`pressure/highrate/` keeps two windows at full rate so the lessons that need
raw data stay runnable without the archive:

- `Pe3_2023-05-22_fullrate.npz` — 16 Hz × 6 h, for the `cosh(kh)` wave
  attenuation demo
- `HBN_2023-05-22_fullrate.npz` — 2 Hz × 24 h, for the stride-vs-average demo

**Bathymetry: 135 MB → 1.7 MB.** Cropped to the working box
(6.88–7.05 °N, 134.05–134.28 °E), `Z` only, float32. The `xUTM`/`yUTM` grids
are dropped. Note the DEM is unreliable over steep reef — within 150 m of the
Peleliu and Angaur gauges the relief is 18–59 m (`PRESSURE_ANALYSIS.md` §6).

By the **Coral Reef Research Foundation (CRRF)**, Koror, Palau —
<https://coralreefpalau.org/>, redistributed here with their permission.
Attribute CRRF in anything that depends on it.

**C05: 56 MB → 268 KB.** Depth-averaged east/north velocity, depth,
temperature and the range bins. The per-bin velocity profiles are not kept.

⚠ **The rotation is NOT applied** to the stored velocities. C05's compass is
skewed — suspected hard-iron from the iron anchor beneath it — and the
correction is `W · exp(−i·14.4°)` absolute. Left to the analysis so the raw
instrument frame stays visible. See `PRESSURE_ANALYSIS.md` §4.

**Excluded entirely**: `fits.mat` (72 MB, regenerable), the CORDC-delivered
NetCDF files (not ours to redistribute — `tests/test_mwb_dat.py` skips without
them and takes a path to your own copy), and the journal PDFs (copyright —
`papers/README.md` has the DOIs).

## Reading it

Pressure and ADCP records are `.npz` with `time_ms` as int64 milliseconds
since the Unix epoch:

```python
import numpy as np
z = np.load("data/pressure/1min/Pe1_1min.npz")
t = z["time_ms"].astype("datetime64[ms]")
dep, temp = z["dep"], z["temp"]
lat, lon = float(z["lat"]), float(z["lon"])
rate = float(z["source_rate_hz"])   # the ORIGINAL instrument rate
```

Or, preferably, through the loader, which handles the vendored/archive choice
and the block-averaging for you:

```python
import sys; sys.path.insert(0, "src")
import pressure_array as pa
t, dep, temp, lat, lon, rate = pa.load("Pe1", step_s=60.0)
```

`step_s` must be a multiple of 60 s against the vendored subset. Anything
finer needs the archive, and the loader will say so.

## Gauge inventory

| group | sites | rate | duration |
|---|---|---|---|
| Peleliu tip | Pe1, **Pe2**, Pe3 | 2, **1**, 16 Hz | 28–31 d |
| Hydrographer Bank | HBN, HBS, HBE, HBW, HBM | 2 Hz | 19–26 d |
| Angaur | An1, An2, An3, An4 | 2 Hz | 29–31 d |

**Pe2 declares `seconds since` in its time units; the other eleven declare
`milliseconds since`.** All are correct and self-describing. Code that
hardcodes milliseconds turns Pe2's 31-day 1 Hz record into a phantom
45-minute record at 1000 Hz — this happened. `xarray.open_dataset` decodes CF
time correctly; raw `netCDF4` readers must parse the unit string.
