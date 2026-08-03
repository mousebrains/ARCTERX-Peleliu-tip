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
| uncertainty | **±16 %** (leave-one-out), *not* the 5 % formal error |
| divergence | 2.3 % of \|ζ\| (contour); 0.1 % (d(ln A)/dt), 0.2 % (LSQ) |
| Okubo–Weiss | **< 0 in 100 % of windows** — coherent for the whole 25.4 h |
| \|strain\|/\|ζ\| | 0.36 |
| cluster | 0.41 km² median (~640 m), vortex core ~1.2 km |

**The vortex is not solid-body.** \|ζ\| falls with cluster scale
(r = −0.56, 95 % CI [−0.63, −0.48], n = 343): 1.58 × 10⁻³ in the smallest
quartile of cluster scale (median 419 m) down to 0.98 × 10⁻³ in the largest
(903 m). Against radius from the fitted centre instead it is −0.37 — the
conclusion is insensitive to that choice but the value is not, so quote the
definition alongside the number. A
Lamb–Oseen fit gives core radius **R ≈ 1.1–1.2 km** (quote the bootstrap
interval, not a point estimate — §5) and Γ ≈ −5400 m² s⁻¹.

That also resolves an apparent discrepancy: the constellation turned **−6.86**
revolutions where ζ/2 integrated predicts −8.83. The 22 % shortfall *is* the
non-solid-body profile — outer drifters orbit more slowly than the core
vorticity implies.

Radial velocity is 4 % of azimuthal, confirming a coherent, non-dispersing
vortex and validating the fitted center.

**Corroborated at this headland by a different instrument class.** Johnston et
al. (2019) measured wake eddies at the Peleliu tip in 2016 from moored ADCPs
and shipboard survey: ~1–2 km diameter, anticyclonic, shed every tidal cycle,
with Rossby number "reaching 80 and 65" — same convention, ζ/f. Our −67 sits
inside their range. Their scaling-based "Ro ~ 30" uses a one-sided shear rather
than ζ, so their own numbers span a factor of two; the circulation/Stokes
estimate here is the better-defined one.

Three of the four drifters roughly **double their orbital radius** over the
record (+111 %, +122 %, +178 %); mwb793d02 started near the core edge and
stayed (−4 %). A slow outward spiral, consistent with the small persistent
positive v_r.

## 5. Traps

**Quote the resampling error, not the formal one.** The leave-one-out spread
(16 % of \|ζ\|) is three times the least-squares formal error (5 %), because
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

Poulain et al. (2023) hit the identical bias in the Cyprus Gyre — a vortex
four orders of magnitude away in Rossby number — and say the drifter estimate
"is an overestimate because the vorticity, in absolute value, always decreases
with increasing distance from the gyre center." It is a property of the
estimator, not of Palau. The cost here: inverting the observed constellation
rotation (−6.86 rev / 25.4 h → 3.70 h) under a solid-body assumption gives
ζ = 4π/P = 9.4 × 10⁻⁴ s⁻¹ against 1.19 × 10⁻³ from four-drifter circulation,
**21 % low** — the §4 shortfall again. Neither is wrong; they average over
different radii. `docs/02-vortex-structure.md` §2.3 works it through.

**Single-trajectory methods inherit that 21 %.** Wavelet ridge analysis,
rotary Fourier and complex demodulation all reach ζ only through ζ = 4π/P, so
none of them can beat the four-drifter circulation estimate on this dataset —
and at 6.9 orbital cycles the record is too short for the one thing wavelets
buy over Fourier. Evaluated and declined; `papers/README.md` entry 7 has the
arithmetic.

Zeiden et al. (2022) ran this exact comparison on Palau drifter clusters and
found a **factor of 2** between a single-drifter wavelet estimate (Ro ~ 6) and
an all-drifter least-squares fit (Ro ~ 3). Same effect, opposite geometry:
their 5 km cluster is larger than the eddy core so their *cluster* estimate
reads low, where our orbital loop is larger than the cluster so our
*constellation* estimate reads low. The rule is the same either way — whichever
estimator integrates over more area reads lower.

**Gate threshold is uncalibrated and looser than the literature.** Huntley et
al. (2022) derive, against a model where truth is known, a cluster-shape
cutoff of Λ = 0.20 (Λ = 1 equilateral, 0 collinear) at scales of 5 km and
below. Our `min_quality = 0.10` on the isoperimetric quotient
(`src/eddy_kinematics.py:336`) is equivalent to Λ = 0.165 — about 20 % more
permissive, against a recommendation the authors describe as already less
stringent than common practice. We also apply one absolute threshold to both
the four-drifter polygon (quotient max 0.785) and the leave-one-out triangles
(max 0.605), so it is relatively looser where there are fewer constraints.
Left as it is: it is uncalibrated but demonstrably not load-bearing (§1.6 of
`docs/`).

**The quotient was also the wrong variable, and that one is now fixed.**
Spydell et al. (2019) state it directly: "cluster area or ellipticity were used
as criteria to distinguish error. We show that the drifter cluster **minor
axis** (narrowness) is a key time-dependent factor." Their Eq. (16) gives
σ_ζ ~ σ_u/λ_a, so a narrow cluster is noisy however well-proportioned its
polygon looks. Our quotient is an area/ellipticity measure. Checking our own
cluster against their criterion:

| | |
|---|---|
| minor axis λ_a, median | 241 m (p5 = 98 m after gating; min 12 m before) |
| epochs below their 50 m danger threshold | 12 of 355 |
| caught by the quotient gate alone | 9 |
| **leaked through** | **3**, carrying median \|ζ\| = 3.2 × 10⁻³ (2.7× the record median) |
| correlation, quotient vs λ_a | r = 0.82 |

**A second gate, `λ_a ≥ 50 m`, is now applied** — their threshold, the point at
which the error exceeds 5f even for velocity errors under 0.004 m s⁻¹. The two
gates correlate at r = 0.82, so the quotient was doing much of this by proxy,
but not all of it. Effect on the published numbers:

| | before | after |
|---|---|---|
| epochs used | 343 | 340 |
| median ζ | −1.1901 × 10⁻³ | **−1.1874 × 10⁻³** (+0.23 %) |
| Rossby | −67.0 | **−66.9** |
| leave-one-out 1σ | 1.710 × 10⁻⁴ (17 %) | **1.622 × 10⁻⁴ (16 %)** |
| 5th percentile \|ζ\| | 2.421 × 10⁻³ | **2.312 × 10⁻³** (−4.5 %) |

A cleaner tail, not a different answer — which is the only honest reason to add
a gate. Sensitivity across λ_a,min ∈ [0, 150] m is 1.3 %, the same insensitivity
the quotient shows; both sweeps are in `docs/01-velocity-gradient-kinematics.md`
§1.6. The MATLAB port carries the same gate and reproduces the new medians.

**Instrument error is not what the 16 % is measuring.** Spydell's Eq. (16)
propagates instrument error alone. On our geometry it gives σ_ζ of 0.4–1.9 % of
|ζ| across plausible velocity errors (0.002–0.010 m s⁻¹), against our 5 % formal
and 16 % jackknife. GPS noise is roughly a tenth of the jackknife, so the
leave-one-out spread really is dominated by flow curvature and unresolved
scales — exactly what it was claimed to measure. Both caveats run the same way:
the calculation set the error correlation to zero, and our Doppler velocities
are block-averaged over 2048 samples, so the true instrument term is smaller
still.

**A better estimator exists and is not implemented.** Zeiden et al. use a
solid-body-*constrained* least-squares fit alongside the plane fit, solving for
`[u0, v0, ζ/2]` in one system. It returns the eddy centre as a fit parameter
rather than by inverting `A⁻¹(c − U₀)`, and its error carries **no
aspect-ratio term** — so it does not degrade as the polygon goes collinear,
which is the failure our gate exists to suppress. It also filters
tides: for their cluster C7 the semidiurnal band held 48 % of plane-fit
vorticity variance but only 18 % of solid-body-fit variance. See
`papers/README.md` entry 15.

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

*Both horns now have supporting evidence, and they are not exclusive.*

For the second horn — circulation accumulating over several cycles — MacKinnon
et al. (2019) measured at Palau's separation point where tides and sheared
low-frequency flow act together, and conclude that including high-frequency
oscillatory currents boosts the net flux of vorticity into the interior by a
depth-dependent factor of **2 to 25**, warning that models omitting them
misestimate momentum and energy loss. A 2–3× shortfall sits at the bottom of
that range. Neither horn is confirmed for our event.

For the first horn: Johnston et al.
(2019) report that flow "intensifies" at the south point as it is constrained
around the topography, and that the velocity difference across the wake eddy
**exceeds 1 m/s** — against our array-averaged 0.64 m/s, a factor ~1.6 on
Γ ~ UL. They also give an intrinsic shedding timescale of ~6 h from a Strouhal
scaling (St ~ 0.2, L = 2 km, U = 0.5 m/s), close enough to the 6.2 h M2
half-cycle that they suggest the tide "may effectively generate eddies." A
1.6× stronger tip flow over a comparable interval plausibly closes a 2–3× gap.
**Still a hypothesis** — it rests on their numbers at their moorings, not ours,
and the Thompson ADCP remains the discriminating measurement.

Two dimensionless numbers from that paper are worth adopting and are not
currently computed anywhere here: an effective Reynolds number
Ref = H/(Cd L) ~ 100, which places this flow in the vortex-street shedding
regime, and the Strouhal number above. Both are cheap to evaluate from data
already in `data/`.

*Attribution corrected.* This document previously credited Ref = H/(Cd L) to
Wolanski et al. (1984). It is not their formula. Their island wake parameter
(Eq. 8) is **P = U H²/(K_z W)**, with K_z the vertical eddy diffusivity and W
the island width; H/(Cd L) is the effective-Reynolds-number form, which
Johnston et al. present alongside a citation to both Wolanski and Signell &
Geyer (1991). Substituting Wolanski's own closure K_z = 0.067 H U∗ makes the
two differ by a factor of ~3.9 at the drag coefficient Johnston measured
(P ≈ 26 against Ref ≈ 100), converging only near Cd = 4.5 × 10⁻³. Both are
≫ 1, so the shedding-regime conclusion is unchanged — but quote the form you
actually used and cite it correctly. `papers/README.md` entry 29 has the
comparison.

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
