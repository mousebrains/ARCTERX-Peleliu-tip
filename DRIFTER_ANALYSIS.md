# Drifters and the Peleliu tip vortex

How four GPS wave buoys were turned into a vorticity measurement, what the
answer is, and which approaches were tried and found not to work.

Companion to `PRESSURE_ANALYSIS.md`. Read the **Traps** section there too —
several lessons apply to both halves of the project.

Code: `src/mwb_dat.py`, `src/mwb_nc.py`, `src/eddy_kinematics.py`,
`src/eddy_analysis.py`. MATLAB port in `matlab/eddy_kinematics_drifters.m`.

---

## 1. What happened

Four CORDC Miniature Wave Buoys (mwb458, 788, 790, 793) were seeded into a
vortex shed from the Peleliu tip. Between **2023-05-22 05:11 and 2023-05-23
06:34 UTC** all four traveled with it. The eddy started 1.68 km from the Pe1
pressure gauge and tracked **west, bearing 266°, for 12 km** over 25.4 h.

The buoys sample GPS position and Doppler velocity at 2 Hz in bursts of 2048
samples — that block is the chunk the buoy uses for its onboard directional
wave spectrum (§2), not a storage artifact.

## 2. The instrument records

Raw `.dat` files are **85-byte fixed records**, little-endian, no header: a
truncated `UBX-NAV-PVT` with lon/lat/height/velN/velE/velD/gSpeed/headMot
stored as float32 rather than u-blox's scaled integers. `src/mwb_dat.py`
documents the layout byte by byte and separates what was **verified
bit-for-bit against the delivered netCDF** from what is inferred from the
u-blox convention.

One burst = 2048 samples = 1023.5 s at 2 Hz. Bursts are internally gap-free;
consecutive bursts are separated by **35–202 s** while the buoy processes and
telemeters, giving a 93.8 % duty cycle.

**The 2048 is a processing block, not a buffer.** The evidence: 2048 × 85 =
174,080 bytes is a multiple of 2048 B but *not* 4096 B, and an 85-byte record
(5 × 17) is storage-hostile — a write-driven design would pad to 88 or 96. The
power of two lives in the *sample* domain. Burst duration is rigid (1023.5 s,
MAD exactly 0) while the *interval* varies (1059–1225 s), which is
"collect N, then process", not a clock scheduler.

### Regenerating the netCDF

`src/mwb_nc.py` decodes the raw bursts to CF-1.13 files. Against the delivered
CORDC files this reproduces `time`, `lat`, `lon`, `u`, `v`, `w` and `numsats`
**bit-for-bit** (max |difference| exactly 0 over 172,032 samples), and `sog`,
`cog`, `pdop` to the delivered files' own quantization.

It also recovers **~38 % more data** — the delivered files trim bursts at both
ends of each deployment — and keeps fields the delivered files drop entirely:
`height`, `hmsl`, `h_acc`, `v_acc`, `s_acc`, `head_acc`, `t_acc`, `nano`,
`itow`, `fix_type`, `flags`, plus **`pdop` at full 0.01 resolution** (the
delivered files store `round(pDOP/100)` in a byte, so theirs is only ever 1
or 2).

## 3. The method

The earlier approach (`matlab/explore_drifter_paths.m`) fits each drifter's
*position* as a shared center plus a per-drifter radius, taking the orbital
direction from the drifter's own velocity. Vorticity is then only reachable
indirectly, and the model has a ±90° sign degeneracy: nothing stops the fit
from returning a negative radius and flipping the phase to compensate.

How often that happens is itself diagnostic. Over 434 windows × 4 drifters
(1736 fitted radii, `matlab/explore_drifter_paths.m` on the regenerated
netCDF):

| interval | fitted radii | negative |
|---|---|---|
| the eddy-coherent window (05-22 05:11 → 05-23 06:34) | 1216 | **11.9 %** |
| outside it | 520 | **30.4 %** |
| whole record | 1736 | 17.5 % |

The degeneracy nearly triples where there is no coherent vortex to fit, which
is the point: the model cannot tell you it has failed, it just returns a
signed radius. The velocity method below has no such mode.

Vorticity is a property of the **velocity** field, and these buoys measure
velocity directly. So estimate it from the velocities, three ways, and
cross-check.

**Primary — circulation / Stokes.** Walk the closed loop joining the four
drifters, integrate the along-track velocity; that is the circulation Γ, and
Γ/A is the area-averaged vorticity. Exact regardless of how complicated the
interior flow is, which matters because the cluster is comparable in size to
the vortex.

**Cross-check — least-squares velocity gradient.** Fit an affine field over a
sliding window; the gradient tensor gives vorticity, divergence and both
strain components, plus Okubo–Weiss. This is Okubo & Ebbesmeyer (1976), whose
Eqs. (1)–(2) are exactly the model implemented — worth more than any
re-derivation.

**Independent check — constellation rotation**, from positions only, using no
velocity data at all.

Divergence likewise gets three estimators: contour flux, least squares, and
`d(ln A)/dt` from the cluster area (positions only).

### Preprocessing

2 Hz velocity is dominated by waves (spectral peak 4.4 s). Each burst is
block-averaged; **bursts are gap-free so averaging never straddles a gap**.
A smoothing window that spans one is built from one-sided data and can be
wrong by 0.62 m/s at burst edges — comparable to the orbital velocity being
measured.

The answer is invariant to this choice: ζ median varies **0.7 %** across a 16×
range of averaging length (64 s → 1024 s).

## 4. The answer

| | |
|---|---|
| **vorticity ζ** | **−1.19 × 10⁻³ s⁻¹** (median, circulation/Stokes) |
| **Rossby number** | **−67** (f = 1.78 × 10⁻⁵ s⁻¹ at 6.99 °N) |
| rotation period | 2.9 h |
| uncertainty | **±17 %** (leave-one-out), *not* the 5 % formal error |
| divergence | < 2 % of \|ζ\| by three independent estimators |
| Okubo–Weiss | **< 0 in 100 % of windows** — coherent for the whole 25.4 h |
| \|strain\|/\|ζ\| | 0.36 |
| cluster | 0.41 km² median (~640 m), vortex core ~1.2 km |

**The vortex is not solid-body.** \|ζ\| falls monotonically with cluster scale
(r = −0.59): 1.58 × 10⁻³ at ~350 m down to 0.98 × 10⁻³ at ~950 m. A
Lamb–Oseen fit gives core radius **R ≈ 1.1–1.2 km** (quote the bootstrap
interval, not a point estimate — §5) and Γ ≈ −5400 m² s⁻¹.

That also resolves an apparent discrepancy: the constellation turned **−6.86**
revolutions where ζ/2 integrated predicts −8.95. The 23 % shortfall *is* the
non-solid-body profile — outer drifters orbit more slowly than the core
vorticity implies.

Radial velocity is 4 % of azimuthal, confirming a coherent, non-dispersing
vortex and validating the fitted center.

Three of the four drifters roughly **double their orbital radius** over the
record (+111 %, +122 %, +178 %); mwb793d02 started near the core edge and
stayed (−4 %). A slow outward spiral, consistent with the small persistent
positive v_r.

## 5. Traps

**Quote the resampling error, not the formal one.** The leave-one-out spread
(17 % of \|ζ\|) is three times the least-squares formal error (5 %), because
the formal error cannot know the flow is curved across the cluster. Molinari &
Kirwan hit the same thing in 1975 — their series were "ragged with frequent
changes in sign" wherever shear was small relative to observational error.

**The Lamb–Oseen fit is the sensitive part.** Two correct implementations
differing only in percentile convention and moving-average edge handling
returned R = 1040 m and R = 1206 m. Repeating the fit over a grid of bin counts
and radius cutoffs spans 1161–1286 m. The medians are robust; this single
global fit to ~12 binned points is not. Hence the bootstrap CI.

**Each ζ is an area average** over whatever the polygon spanned at that
moment. Because the vortex is not solid-body, part of the apparent time
variation is the cluster sampling different radii, not the vortex changing.
The radial profile separates them.

**Gate on polygon shape.** Leave-one-out triangles can go near-collinear, and
Γ/A then diverges — spikes to ζ = 0.27 s⁻¹ (Rossby 15,000) before gating on
the isoperimetric quotient 4πA/P².

**Condition the design matrix by hand.** Its columns carry different units
(1, m, m, s); unscaled, `cond` reports the unit mismatch, ~2000, rather than
the drifter geometry, ~2.

**The center needs an assumption, unavoidably.** At a single instant a uniform
background flow and a displacement of the vortex center are *exactly*
degenerate: ω ẑ×(x−c) contains the constant −ω ẑ×c. We break it by low-passing
the drifter-mean velocity over ~2 orbital periods, so every center estimate is
conditional on that. If the center track looks wrong, suspect this first.

## 6. Established negatives

**Drifter GPS height cannot measure sea level — by ~100×.** Burst-mean height
scatter is 3.3–6.4 m; differencing between drifters barely helps (2.0–4.4 m)
and 3 h smoothing leaves 1.4–3.9 m, against a 1–3 cm expected eddy signature.

The failure of common-mode cancellation is diagnostic: the error is
**platform-specific** (hull multipath, antenna motion), not atmospheric.
Decomposed against the buoys' own Doppler velocity — far more precise than GPS
position — height is usable near the wave peak and is essentially all
instrument below ~0.02 Hz, the band that matters for sea level. In the 5–32 s
band, Hs from height is 1.39 m against 0.65 m from velocity and 0.49–0.55 m
from flat-bank bottom pressure.

Use `v_acc` as the quality filter: it predicts height error (r = 0.46–0.99)
while `pdop` and `numsats` do not. It is in the regenerated netCDF; the
delivered files drop it.

**The eddy is invisible in bottom pressure.** Gradient-wind balance predicts a
4.1 cm depression at the center, ~1 cm at the nearest gauge's 1.68 km
approach, against a 3.3 cm *oceanographic* residual floor that no averaging
removes. See `PRESSURE_ANALYSIS.md` §6.

## 7. Open work

**Circulation budget.** Γ = −5400 to −6100 m² s⁻¹ needs roughly 2–3× more
forcing than one M2 half-cycle at the measured 0.64 m/s tip current supplies.
Either flow past the tip exceeds the array-averaged value, or circulation
accumulates over several cycles. The Thompson ADCP spot measurements could
discriminate.

**Shedding direction.** C05 shows the *sign* of the along-axis flow is set by
the subtidal current, with clean regimes in May 2023 (west 14–23, east 24–27).
The drifter eddy sits inside a westward event, consistent with its westward
translation. X-band radar can test the prediction out of sample —
`PRESSURE_ANALYSIS.md` §4 and §8.

**MATLAB port** (`matlab/eddy_kinematics_drifters.m`) reproduces the Python to
the digits shown in §4, and adds bootstrap CIs on the Lamb–Oseen fit. It needs
the Mapping, Statistics and Machine Learning, and Curve Fitting toolboxes.
The Python is the reference implementation — see the parity table in
`README.md` for what "reproduces" means and what to do when they disagree.
