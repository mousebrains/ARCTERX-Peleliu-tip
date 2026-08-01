# Bottom-pressure array, ARCTERX 2023 Peleliu Wake

What the pressure sensors can and cannot tell us, how to read them without
being misled, and what has been established so far.

Written for someone joining the project. Read the **Traps** section before you
write any code — every one of them produced a confidently wrong answer first
time round, and most of them look like success while they are doing it.

Code: `pressure_array.py` (loader + harmonic analysis).
Tidal machinery lives in `wamos_tpw.tides` (the WAMOS repo), branch
`feature/tidal-harmonics`.

---

## 1. The instruments

Twelve bottom-mounted wave gauges, all in 17–21 m of water, spanning ~11 km.

| group | sites | location | notes |
|-------|-------|----------|-------|
| Peleliu tip | Pe1, Pe2, Pe3 | 6.971–6.983 N, 134.218–134.233 E | forms a usable gradient triangle |
| Hydrographer Bank | HBN, HBS, HBE, HBW, HBM | 6.921–6.935 N, 134.193–134.203 E | HBM is on the SeaSpider |
| Angaur | An1, An2, An3, An4 | 6.906–6.923 N, 134.135–134.155 E | |

`/Volumes/SeaChest/ARCTERX/2023/Wake/Pressure Sensors/` holds eleven;
HBM is in `Bank Seaspider/060466_20230602_1919_HBM_clipped.nc`.

Variables: `pdbar` (total pressure, dbar), `dep` (water level, m), `temp`,
`time`, scalar `lat`/`lon`. HBM adds `cond` and `salt`.

**Sample rates differ and are not documented anywhere but the file.**

| site | rate | duration |
|------|------|----------|
| Pe3 | 16 Hz | 28.2 d |
| Pe2 | **1 Hz** | **31.1 d** (longest) |
| all others | 2 Hz | 19.2–30.9 d |

Near the middle of Hydrographer Bank, within ~3 m of each other:

- **C05**, an RDI Workhorse ADCP — 12-min ensembles, Dec 2022 – Jun 2023 (182.5 d)
- the **SeaSpider**, carrying Sig1000 (16 Hz), a CTD (2 Hz), and two
  MicroRiders (1 kHz) — all with pressure, colocated within ~1 m

That near-colocation is scientifically valuable: see §5.

**Rates, verified against the files** (not from memory — three of the four
recollections were right, one was not):

| instrument | measured | evidence |
|---|---|---|
| CTD | 2 Hz | `hotel.mat`, dt = 0.5 s; time is Unix-epoch **seconds** |
| Sig1000 | 16 Hz | `Burst_SamplingRate = 16`; 43 cells x 0.5 m, 4 beams, ENU |
| MicroRiders | **1024 Hz** | `rate = 1024` in the ODAS header, s/n 330 and 429 |
| C05 Workhorse | **12 min** | 720.0 s uniformly, raw and processed |

The MicroRiders are the only genuine kHz instruments in the experiment.

---

## 2. The physics, in plain terms

A bottom pressure sensor weighs the water column above it. Three things change
that weight:

1. **Tide** — ±1 m here, and it dominates everything (99 % of the variance).
2. **Surface waves** — but they are *attenuated with depth*. A wave of length
   L is felt at the bottom by a factor `1/cosh(kh)`. Short waves die before
   they reach 19 m; long swell mostly survives. This is why the gauges see the
   12 s swell and are blind to the 4.4 s wind sea.
3. **Everything else** — weather, internal waves, currents. This is the ~4 cm
   residual left after the tide is removed, and most of §5 is about what it is.

The practical consequence for wake work: a bottom pressure sensor **is** a sea
level gauge, so in principle it sees the surface depression of a passing eddy.
In practice the numbers do not work — see §6.

---

## 3. Traps

### 3.1 Read the declared time units. Do not assume.

Ten gauges declare `milliseconds since ...`. **Pe2 declares `seconds since ...`.**
Both are correct and self-describing.

Hardcoding milliseconds turns Pe2's 31-day 1 Hz record into a phantom
45-minute record at 1000 Hz — which is what happened here, and it was reported
as a dead sensor before the error was found. The physics gives it away: the
record showed 2.7 m of range with 0.46 m sd, and **you cannot get a full tidal
swing in 45 minutes.**

`xarray.open_dataset` decodes CF time correctly and is immune. Raw `netCDF4`
readers are not — `pressure_array.load` now parses the unit string explicitly.

### 3.2 Block-average, never stride.

At 19 m depth an 8 s wave still reaches the bottom at ~58 % amplitude, so 0.5 m
of swell puts ~25 cm into bottom pressure. Decimating a 2 Hz record to 1 min by
taking every 120th sample folds all of that into the tidal band.

Measured cost of getting this wrong:

| | strided | block-averaged |
|---|---|---|
| tidal variance explained | 90–97 % | **98.9–99.3 %** |
| post-fit residual | 7.5–14.4 cm | **3.8–4.5 cm** |
| M2 phase spread over 11 km | 149° | **1.7°** |

A 149° spread implies the tide crossed 11 km in 5.1 hours — a 0.6 m/s wave
speed, when shallow water gives 13 m/s. **Physically impossible answers are
the cheapest bug detector you have.** Compute the implied quantity and sanity
check it.

### 3.3 Reference all phases to one epoch.

The records start at different times. Phase referenced to each record's own
start is meaningless between instruments. `wamos_tpw.tides` carries the epoch
in the result and defaults to a fixed one.

### 3.4 Do not fit constituents the record cannot separate.

Two constituents `ds` degrees/hour apart need `360/ds` hours to separate
(Rayleigh). K1/P1 and S2/K2 need 183 days; Mf/MSf likewise; M2/N2 needs 27.6 d,
so the 26-day Hydrographer records cannot carry N2.

But **pairwise Rayleigh is necessary, not sufficient**, and this is subtle:

- K1 and P1 *alone* in a 29-day record recover fine — they still drift 57°
  apart. Rayleigh is conservative.
- What actually breaks it is stacking several close lines. Adding **S1
  (24.000 h)** beside K1 (23.934) and P1 (24.066) makes a three-way degeneracy.
  The design-matrix condition number goes from ~2.9e3 to ~7.7e4 and P1 comes
  out at 2–9× truth.

And the trap inside the trap: **the corrupted fit has the *lower* residual.**
On real data, forcing the bad set gave 3.08 cm residual against 4.06 cm for
the clean fit. Residual and R² both say "better". Only the condition number
exposes it. `TidalFit.condition` reports it; a warning fires above 1e4.

The right fix for a short record is **inference**: take the amplitude ratio and
phase lag of the minor constituent from a long reference record and impose it,
adding no free parameter. See §7.

### 3.5 Match footprints before comparing anything.

Two examples that both produced fake discrepancies here:

- Comparing the SeaSpider gauge's 2-min means *sampled* at 12 min against the
  ADCP's 12-min *averages* inflated their difference. Block-averaging both to
  12 min fixed it.
- Comparing buoy Hs (3–25 s band) against bottom-pressure Hs (5–32 s band,
  all the sensor can resolve at 19 m) gave 0.70 vs 0.40 m and looked like a
  measurement problem. It was two different bands.

### 3.6 File-specific gotchas

- **`C05_2023.mat` `pressure`** is in **deca-pascals** and contains
  `4294967286` fill values (uint32 overflow). `depth` is in metres but carries
  `0.000` for 181 out-of-water ensembles. Use `depth` with a `>5 m` mask.
- **Pe2 segfaults inside HDF5 on `Dataset.close()`** after a large read
  ("There are 1 HDF5 objects open!"). The reads are correct. `pressure_array`
  deliberately never closes datasets — do not reintroduce a context manager.
- **C05 is 12-minute ensembles**, uniformly 720.0 s, in both the raw and
  processed files.

---

## 4. Tidal analysis

Mixed, mainly semidiurnal (form factor ≈ 0.61). Constituents from
`wamos_tpw.tides.harmonic_fit`, phases on a common epoch:

| | amplitude | notes |
|---|---|---|
| M2 | 0.494–0.506 m | spread of only 2.3° in phase across the whole 11 km array |
| S2 | 0.178–0.183 m | |
| K1 | 0.233–0.245 m | |
| O1 | 0.161–0.173 m | |

Fit quality depends on record length, because the 26-day Hydrographer records
must drop N2:

| records | constituents | variance explained | residual |
|---|---|---|---|
| 28–31 d (Pe, An) | 9 | 99.05–99.17 % | 4.0–4.2 cm |
| 26 d (HB) | 8 (N2 dropped) | 96.8–97.2 % | 7.1–7.6 cm |

The HB penalty is entirely the missing N2. Inferring it (§7) should recover it.

### Tidal currents from pressure gradients — WITHDRAWN

An earlier version of this document reported M2 tidal currents of 63.6 cm/s
(Peleliu tip), 45.9 (Angaur) and 20.3 (Hydrographer) from the cluster
sea-level gradients. **Those numbers are not supported. Do not use them.**

Three independent reasons, in increasing order of how fundamental they are:

1. **The gradient is not resolvable.** The M2 sea-level difference across a
   cluster, against the harmonic-amplitude uncertainty:

   | cluster | span | signal | uncertainty | SNR |
   |---|---|---|---|---|
   | Hydrographer | 1.9 km | 4.2 mm | 4.4 mm | **1.0** |
   | Peleliu tip | 2.1 km | 16.3 mm | 3.8 mm | 4.2 |
   | Angaur | 2.9 km | 10.0 mm | 3.7 mm | 2.7 |
   | whole array | 11.4 km | 28.2 mm | 3.8 mm | 7.4 |

   At SNR 1.0 the Hydrographer gradient is pure noise.

2. **It fails validation.** Reconstructing the bank current from the HB
   gradient and comparing with the C05 ADCP gives r = −0.18 (east) and −0.22
   (north). The 11 km regional reconstruction, where SNR is 7.4, gives
   r ≈ 0.00.

3. **The reason it can never work here.** Only **17 %** of the measured east
   current is phase-locked tide (10 % north). The other 83 % is subtidal and
   baroclinic — invisible to bottom pressure, which senses only the
   barotropic mode. No improvement in gradient precision recovers it.

This is a structural limit of the instrument, not a precision problem. The
pressure array measures **sea level** well and **current** not at all.

### Use the measured current instead — C05

The RDI Workhorse on the bank measures current directly: 181 days
(Dec 2022 – Jun 2023), 12-min ensembles, 43 bins.

⚠ **Apply the compass rotation.** The WAMOS `thompson2023/adcp_validation`
work found the Workhorse heading is skewed — suspected hard-iron from the
iron anchor under C05. Independently reproduced here from the complex
transfer against the Sig1000 (4–14 m band, |v| > 0.15 m/s):

| | WAMOS notes | reproduced here |
|---|---|---|
| gain | 0.985 | 0.983 |
| angle | +21.09° | **+21.00°** (n = 2958) |

Corrections are complex multipliers: `W·exp(−i·21.09°)` maps the Workhorse
into the Sig1000 frame, and since SIG is itself +6.6° off true, the absolute
correction is **`W_corr = W · exp(−i·14.4°)`**. The skew is stable — +20.2 to
+21.7° in every 5-day chunk of the deployment.

After correction:

| | |
|---|---|
| principal axis | **113.4° / 293.4° true**, ellipticity 0.35 |
| depth-averaged sd | 55.5 (east) / 21.3 (north) cm/s, max speed 220 cm/s |
| phase-locked tide | 17 % east, 10 % north |
| **direction set by** | **subtidal flow** — sign(total) = sign(subtidal) 73 % of the time |

The tide supplies 61 % of the along-axis *variance* but the subtidal flow sets
the *sign*. That matters for shedding direction.

**Event list**: `c05_high_flow_events.csv` — 141 events longer than 1 h above
the 85th percentile (81 m/s threshold): **81 eastward, 60 westward**.

The east/west labels are **insensitive to the rotation** (0.0 % of ensembles
change sign under it), because the flow is strongly rectilinear. The rotation
matters for the axis bearing, not the classification.

For May 2023, the subtidal daily means organise into clean regimes:

| dates | regime |
|---|---|
| May 14–23 | **WEST** (−36 to −74 cm/s) |
| May 24–27 | **EAST** (+22 to +71) |
| May 28–30 | WEST (−13 to −83) |
| May 31 | EAST (+23) |

The drifter-observed eddy sits inside a 7-hour westward event,
**2023-05-22 04:12–11:12, peak −115 cm/s**, during a westward subtidal
regime — consistent with its observed translation west at bearing 266°.

**Testable prediction for the radar**: eddies detaching **eastward** around
**May 24–27** and **May 31**, westward through May 14–23.

⚠ Caveat: C05 is on the bank, 5–6 km from the tip, in 18 m of water while the
tip drops to 600–1900 m. Since 83 % of its signal is non-tidal it cannot be
extrapolated to the tip by any simple rule. Treat it as a **regime indicator
whose sign is meaningful**, not as the tip current. A geometric check against
the eddy translation bearing was tried and is *not* valid for this reason —
the bank's principal axis has no obligation to match a deep-water eddy's
heading.

## 5. What the residual is

After removing the tide, ~4 cm remains. It is **red** — 83 % of its variance is
at periods of 6 h to 5 days — so it does *not* average down (1 min → 2 h
smoothing takes 4.17 → 3.82 cm only).

Things it is **not**:

- **Not temperature/steric.** Regressing residual on bottom temperature removes
  ~1 % at every site.
- **Not simple sensor noise.** Two gauges 600 m apart differ by 1.5 cm, so most
  of it is coherent at that scale.

The colocated pair settles the split. HBM (SeaSpider) and C05 (ADCP) sit 3 m
apart, so any difference between them is instrumental by construction:

| | |
|---|---|
| raw difference, matched 12-min footprints | 3.63 cm |
| after removing a **−308 s clock offset** and a **−1.25 % gain difference** | 3.08 cm |
| **per instrument** | **2.18 cm** |
| after 1 h / 3 h / 12 h smoothing | 1.72 / 1.13 / 0.66 cm |

So of the ~4.0 cm residual, roughly **2.2 cm is instrument and 3.3 cm is real
ocean** (in quadrature). Both are red, which is why neither averages away.

The 5-minute clock offset between two instruments on the same mooring is worth
remembering — see §7 on the multi-year record.

---

## 6. Established negatives

Recording these so nobody spends a month rediscovering them.

### The eddy is not detectable in bottom pressure

From the drifter-measured velocity profile, gradient-wind balance predicts a
**4.1 cm depression at the eddy centre** (centrifugal beats Coriolis 24:1, so
it is a low regardless of rotation sense), falling to 2.5 cm at 1 km and
0.85 cm at 2 km.

The eddy's closest approach to any sensor was **1.68 km (Pe1)**, where the
signature is ~1 cm — against a **3.3 cm oceanographic** floor that no averaging
removes. It is below the noise, and the noise is irreducible.

A tidal-phase composite of the Peleliu-local anomaly *looked* like a 5.4 cm
locked signal, but its amplitude scales with baseline length to the reference
cluster (5.44 / 3.34 / 2.31 / 0.82 cm at 11.4 / 11.8 / 6 / 1 km) — it is the
barotropic tidal gradient, not a wake signature. **Always test a "local" signal
against baseline length.**

### Drifter GPS height cannot see it either — by ~100×

Burst-mean height scatter is 3.3–6.4 m; differencing between drifters barely
helps (2.0–4.4 m) and 3 h smoothing leaves 1.4–3.9 m, against a 1–3 cm
expected signature. The failure of common-mode cancellation is diagnostic: the
error is **platform-specific** (hull multipath, antenna motion), not
atmospheric.

Decomposed against the drifters' own Doppler velocity (which is far more
precise than GPS position), height is usable only near the wave peak and is
essentially all instrument below ~0.02 Hz — the band that matters for sea
level. Use `v_acc` as a quality filter; it is the best predictor of height
error (r = 0.46–0.99) while pdop and numsats are not.

### Wave spectra are only trustworthy on the flat bank

The `cosh(kh)` correction assumes a locally flat bottom over a wavelength
(~150 m for 12 s waves). Within 150 m of each sensor the DEM relief is:

| site | relief | verdict |
|---|---|---|
| HBN, HBS, HBE | 0.8–1.6 m | flat — trustworthy |
| HBW | 32 m | bank edge — suspect |
| Pe1, Pe2, Pe3, An1–4 | 18–59 m | steep reef — **not defensible** |

The three flat sites agree on Hs (0.49–0.55 m, Tp 12.0–12.3 s). Pe3 gives
0.39 m and should not be used.

Incidental but useful: the DEM bias on the flat bank comes out **+0.73 m** from
the wave-gauge `dep`, independently reproducing the **+1.2 m** `--depth-adjust`
calibrated in `wamos_tpw` from the Seaspider pressure. At Peleliu and Angaur
the DEM is wrong by +22 to +58 m and should not be used at all.

---

---

## 7. Open work

**Multi-year C05 for tidal admittance.** The single 182.5-day deployment is
*just* short of the 182.6 days K1/P1 needs. A multi-year record resolves K1/P1,
S2/K2, and the long-period constituents (Mm, Mf, Ssa, Sa) that no 30-day
deployment can touch. Being on the bank itself, it beats transferring ratios
from Malakal across the reef. Feed the result to
`harmonic_fit(..., infer={"P1": ("K1", ratio, lag)})`.

⚠ **Time skew and drift between C05 recoveries and redeployments must be
compensated before concatenating deployments.** This is not hypothetical: a
−308 s offset was measured between two instruments on the *same mooring* in a
single deployment (§5). Concatenating without correction will smear the
diurnal band, which is exactly the band the exercise is meant to resolve.
Estimate the offset per deployment by regressing the difference against
`d(elevation)/dt` — a clock offset τ appears as `τ · dH/dt`.

**Shedding direction vs subtidal flow.** The strongest lead in this work.
C05 (rotation applied) gives 181 days of measured current whose *sign* is set
by the subtidal flow, and the May 2023 regimes make a specific prediction
(§4). The radar can test it directly, and if it holds the mechanism extends to
the full multi-year C05 record. Note this supersedes the earlier plan to use
pressure-derived tidal phase, which does not work (§4).

**Circulation budget.** The observed eddy circulation Γ = −6114 m² s⁻¹ needs
roughly 2–3× more forcing than one M2 half-cycle at the measured 0.64 m/s tip
current would supply. Either the flow past the tip exceeds the array-averaged
value or circulation accumulates over several cycles. The Thompson ADCP spot
measurements could discriminate.

**Unverified instrument rates.** CTD ~2 Hz, Sig1000 16 Hz, MicroRiders 1 kHz
are from memory and have not been checked against the files. The MicroRiders
are the only genuine kHz instruments in the experiment.

---

## 8. Co-tidal chart, and what one mooring can validate

Two products the array *does* support, and the limits they imply for the 2025
radar work, where **C05 is the only current mooring**.

### 8.1 The co-tidal chart (resolved)

A plane fit of the complex tidal amplitude across all twelve gauges resolves
the **phase** gradient cleanly, though not the amplitude gradient:

| constituent | phase gradient | propagation | significance |
|---|---|---|---|
| **M2** | **0.761 ± 0.102 °/km** | **314° at 11 m/s** | 7.5σ |
| S2 | 0.877 ± 0.283 | 311° at 9 m/s | 3.1σ |
| O1 | 0.610 ± 0.308 | 328° at 6 m/s | 2.0σ |
| K1 | 0.061 ± 0.215 | — | not resolved |
| M2 amplitude | 1.4 ± 1.8 mm/km | — | not resolved |

Three constituents independently agreeing on 311–328° is real. The amplitude
being flat across 11 km is consistent with the near-identical M2 amplitudes in
§4.

**The phase speed is the interesting part.** 11 m/s implies an effective depth
of ~11 m — close to √(g·19) = 14 m/s for the bank tops, and nowhere near
√(g·1500) = 121 m/s for the channel. The array is sensing a *shallow-water
controlled* wave over the reef complex, not the deep-channel barotropic tide.

That may be the physical reason the gradient-to-current inversion in §4 fails:
the momentum balance used there assumes the deep barotropic mode, which is not
what these gauges are measuring.

### 8.2 The legitimate route from sea level to current

Not an inversion — a **forward model**. Use the co-tidal chart to constrain a
barotropic tidal model of the channel with real bathymetry, and take the
current from the model. That respects friction and topography, which the
analytic momentum balance ignores and which set the local response. The numbers
in §8.1 are exactly the constraint such a model needs.

### 8.3 How far one mooring reaches

From the 7-platform structure function in
`thompson2023/adcp_validation/NOTES.md` (session 1, step 3), against a field
standard deviation of 0.42 m/s per component:

| separation | correlation | variance explained | uncertainty |
|---|---|---|---|
| 1.5 km | 0.90 | 81 % | ±0.19 m/s |
| 4 km | 0.84 | 70 % | ±0.24 |
| 6.5 km | 0.79 | 63 % | ±0.27 |
| 10 km | 0.69 | 48 % | ±0.33 |

A single mooring constrains the field usefully to ~10 km — **enough for regime
and sign, marginal for magnitude**. That is exactly the precision the
shedding-direction question in §4 needs, and it is why C05 alone supports a
coarse channel-wide picture across its whole record.

### 8.4 Validation footprint for the 2025 Angaur radar

Geometry (radar tower at 6.91677 N, 134.14840 E):

| | from tower | from C05 |
|---|---|---|
| C05 | 5.83 km | — |
| HBM / SeaSpider | 5.78 km | 0.05 km |
| Angaur An2 | 0.59 km | 5.98 km |
| Peleliu tip (Pe2) | 10.26 km | 5.23 km |

**A perfect radar cannot agree with C05 better than the ocean decorrelates.**
Expected scatter from spatial decorrelation alone, against an assumed 0.15 m/s
radar error:

| tile–C05 separation | decorrelation floor | combined | dominated by |
|---|---|---|---|
| 0.12 km | 0.033 | 0.154 | radar error |
| 0.75 km | 0.138 | 0.204 | radar error |
| **0.9 km** | **0.15** | 0.21 | **crossover** |
| 2.5 km | 0.204 | 0.253 | decorrelation |
| 6.5 km | 0.271 | 0.310 | decorrelation |
| 10 km | 0.331 | 0.363 | decorrelation |

Consequences for 2025:

1. **Validate within ~1 km of C05.** Beyond that a radar-vs-mooring scatter
   plot measures the ocean, not the instrument. A 0.3 m/s scatter at 6 km is
   the *expected* result for a perfect radar and must not be read as radar
   error.
2. **Quote the decorrelation floor alongside every comparison**, matched to the
   tile–mooring separation. Without it, skill is systematically understated at
   range.
3. **C05 is 5.8 km from the tower**, so the tiles that can validate it sit at
   mid-range, not near the radar. Tiles near the tower (where the radar is
   best) are ~6 km from C05 and cannot be validated tightly against it.
4. The 2023 seven-platform ensemble is the only configuration that can
   validate the radar *across* the channel. Use it to characterise radar skill
   as a function of range, then carry that characterisation into 2025 where
   only C05 remains.
