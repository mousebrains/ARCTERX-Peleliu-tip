# 7. Data formats

The raw instrument record, and the conventions the generated netCDF files
satisfy.

Code: `src/mwb_dat.py` (parser), `src/mwb_nc.py` (writer),
`tests/test_mwb_dat.py`, `tests/test_mwb_nc.py`.

---

## 7.1 The raw `.dat` record

The Miniature Wave Buoys write **85-byte fixed records**, little-endian, with
no file header and no framing. The layout is a **truncated `UBX-NAV-PVT`**,
u-blox's standard position-velocity-time message, with one important
modification: `lon`, `lat`, `height`, `velN`, `velE`, `velD`, `gSpeed` and
`headMot` are stored as **float32 rather than u-blox's scaled integers**.

That substitution has a real consequence. u-blox stores longitude as an int32
in units of 1e-7 degrees — about 1 cm. A float32 has 24 bits of mantissa, so
near 134° the spacing is

$$\mathrm{ulp}(134°) = 2^{-24}\times 2^{8} \approx 1.53\times10^{-5}\ \text{degrees} \approx 1.69\ \mathrm{m}$$

while latitude near 7° has $\mathrm{ulp} \approx 4.8\times10^{-7}$ degrees
$\approx 0.05$ m. **Longitude is quantized 34× more coarsely than latitude**,
purely because of where the exponent boundary falls. It is below the GPS noise
floor so it does not affect any result here, but it is a real asymmetry and it
would matter for a precise-positioning application.

`src/mwb_dat.py` documents the layout byte by byte and — importantly —
separates what was **verified bit-for-bit against the delivered netCDF** from
what is **inferred from the u-blox convention**. Those are different epistemic
categories and the file marks which is which.

### There is no temperature, and no hidden field

The record ends with **7 bytes that carry no named field** (offsets 78–84).
In a full `UBX-NAV-PVT` those offsets hold `reserved1[6]` and the start of
`headVeh`; in this truncated 85-byte variant they could in principle hold
anything the buoy chose to append — a thermistor reading being the obvious
candidate.

They do not. Across all **239,616 records** of the vendored deployment every
one of the seven bytes is `0x00`, with exactly one distinct value each. It is
padding, not data.

So the `.dat` files are **purely GPS**: position, velocity, time and the
receiver's own quality estimates. There is **no temperature channel**, and no
other undocumented measurement hiding in the record. Sea-surface temperature
for this experiment has to come from the pressure gauges (which record it), the
CTD, or the ship.

### Time reconstruction

UTC comes from the `year/month/day/hour/min/sec` fields plus a signed `nano`
correction. Converting to epoch milliseconds uses Howard Hinnant's civil-days
algorithm, which is exact integer arithmetic with no floating-point date
handling and no library dependency.

## 7.2 What regeneration recovers

Reprocessing the raw `.dat` rather than using the delivered files gains three
things:

**1. Bit-for-bit agreement on the shared fields.** `time`, `lat`, `lon`, `u`,
`v`, `w` and `numsats` reproduce with **max |difference| exactly 0** over
172,032 samples, and `sog`, `cog`, `pdop` to the delivered files' own
quantization. That is what makes the rest trustworthy.

**2. About 38 % more data.** The delivered files trim bursts at both ends of
each deployment.

**3. Fields the delivered files drop entirely**: `height`, `hmsl`, `h_acc`,
`v_acc`, `s_acc`, `head_acc`, `t_acc`, `nano`, `itow`, `fix_type`, `flags`.

`v_acc` is the one that mattered — it is the only usable quality filter for
GPS height (§6.5) and it is absent from the delivered product.

Also `pdop` at full 0.01 resolution: the delivered files store
`round(pDOP/100)` in a byte, so theirs is only ever 1 or 2.

> The bit-for-bit test lives in `tests/test_mwb_dat.py`. The CORDC-delivered
> files are not redistributed here, so it skips unless you pass a path to your
> own copy. `tests/test_mwb_nc.py` needs no delivered file and checks the
> round trip raw → netCDF → raw.

## 7.3 CF conventions

The generated files declare `Conventions = "CF-1.13, ACDD-1.3"`.

CF (Climate and Forecast) governs the physical content: units, standard names,
coordinates, cell methods, the calendar. ACDD (Attribute Convention for Data
Discovery) governs discovery metadata: title, summary, creator, extents. They
are complementary and both are worth satisfying.

### Choices worth explaining

**`standard_name` only where one genuinely exists.** GPS accuracy estimates,
DOP, satellite counts and receiver bitfields have no CF standard name (checked
against table v94). Inventing one to satisfy an ACDD recommendation would
*break* CF compliance, which is the more important of the two. Those variables
get a `long_name` and nothing else.

**`calendar = "standard"`, not `"utc"`.** `"utc"` is not a CF calendar value.
This is an easy mistake because the data genuinely is UTC.

**A separate `sampling_period` dimension.** CF forbids a multidimensional
variable sharing a name with a dimension, so the natural naming collides and
has to be broken.

**`positive = "up"` on the height variables.** CF §4.3 requires any variable
carrying a vertical-coordinate standard name to declare `positive`. `hmsl`
carries `standard_name = "altitude"`, which is a vertical coordinate, so
omitting it is a hard error. This was missed initially — `cfchecker` only
implements up to CF-1.8 and does not flag it; the IOOS `compliance-checker` at
`cf:1.11` does.

**Compression**: zlib level 4 by default. Levels above ~4 buy very little on
float32 geophysical data while costing noticeably more time.

### Verifying compliance

```bash
pip install compliance-checker
compliance-checker --test cf:1.11 data/drifters/mwb458d02_gps_timeseries.nc
```

All four files return **0 errors**. Two warnings remain, and both are artifacts
of checking a CF-1.13 file against CF-1.11 rules:

- `Conventions does not contain "CF-1.11"` — it contains CF-1.13, which is
  correct and newer.
- `units_metadata is recommended` — a CF-1.11 recommendation that CF-1.13
  withdrew. Adding it would satisfy the 1.11 checker while being wrong for the
  version the file declares.

Do not "fix" either one. Checking against an older version than the file
targets will always produce this kind of noise; read the messages rather than
chasing a zero.

## 7.4 The pressure and ADCP subset

The vendored records are `.npz` with `time_ms` as int64 milliseconds since the
Unix epoch:

```python
import numpy as np
z = np.load("data/pressure/1min/Pe1_1min.npz")
t = z["time_ms"].astype("datetime64[ms]")
dep, temp = z["dep"], z["temp"]
rate = float(z["source_rate_hz"])    # the ORIGINAL instrument rate
```

Preferably go through the loader, which handles the vendored-versus-archive
choice and the block averaging:

```python
import sys; sys.path.insert(0, "src")
import pressure_array as pa
t, dep, temp, lat, lon, rate = pa.load("Pe1", step_s=60.0)
```

### The time-units trap

**Pe2 declares `seconds since` in its time units; the other eleven gauges
declare `milliseconds since`.** All twelve are correct and self-describing.

Code that hardcodes milliseconds turns Pe2's 31-day 1 Hz record into a phantom
45-minute record at 1000 Hz — and this happened, and was misdiagnosed twice
(first as a dead sensor, then as a bug in the file) before the header was read
properly. Pe2 is in fact the **longest** record in the array and its M2 phase
matches the rest to 0.3°.

`xarray.open_dataset` decodes CF time correctly. Raw `netCDF4` readers must
parse the unit string — `_time_scale` in `tools/make_data_subset.py` shows the
pattern.

The general rule this illustrates: **a header is a claim that must be read, not
a convention that can be assumed.** Never carry a convention across
instruments, cruises or file versions without verifying it in the data at hand.

### A stability note

`pressure_array` deliberately never closes its netCDF `Dataset` objects and
holds them in a module-level list. Pe2 segfaults inside HDF5 on close after a
large read. Leaking the handle for the life of the process is the lesser evil.

## References

- u-blox (2020). *u-blox 8 / u-blox M8 Receiver Description — Protocol
  Specification*, UBX-13003221. §32.17.15.1 `UBX-NAV-PVT`.
- Eaton, B. et al. (2024). *NetCDF Climate and Forecast (CF) Metadata
  Conventions*, version 1.11+. <https://cfconventions.org/>
- NOAA/IOOS (2015). *Attribute Convention for Data Discovery (ACDD) 1.3*.
- Hinnant, H. (2021). `chrono`-compatible low-level date algorithms.
  <https://howardhinnant.github.io/date_algorithms.html>
