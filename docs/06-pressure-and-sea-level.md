# 6. Pressure and sea level

Turning bottom pressure into sea level, establishing what the instruments can
actually resolve, and two negative results that are worth more than they look.

Code: `src/pressure_array.py`. Results: `PRESSURE_ANALYSIS.md`.

---

## 6.1 The hydrostatic relation

For a motionless column of water the pressure at the bottom is

$$p = p_{\text{atm}} + \rho g h$$

so a change in bottom pressure maps to a change in sea-surface height:

$$\Delta \eta = \frac{\Delta p}{\rho g}$$

With $\rho \approx 1025\ \mathrm{kg\,m^{-3}}$ and $g = 9.81\ \mathrm{m\,s^{-2}}$,
1 dbar ≈ 0.995 m. At tidal frequencies this is excellent: the flow is slow, the
column is close to hydrostatic, and the gauges give sea level directly.

Three corrections matter in principle:

1. **Atmospheric pressure.** Over a 29-day record the inverse-barometer effect
   is real (~1 cm/mbar) but at tidal frequencies it is small and largely common
   across an 11 km array, so it cancels in the phase gradients of §5.5.
2. **Density.** A temperature-driven density change alters $\rho g h$ without
   any change in sea level. The gauges record temperature for this reason.
3. **Dynamic pressure.** Where the flow is fast the column is *not*
   hydrostatic. This is what §6.5 turns out to be about.

## 6.2 Depth attenuation: why a bottom gauge is a low-pass filter

Surface waves do not reach the bottom equally. Linear wave theory gives the
pressure perturbation at depth $z$ below the surface, in water of depth $h$:

$$p'(z) = \rho g\, a \frac{\cosh\!\left(k(h+z)\right)}{\cosh(kh)}\, e^{i(kx-\omega t)}$$

At the bottom ($z = -h$) this becomes $p' = \rho g a / \cosh(kh)$, so the
pressure response function is

$$K_p = \frac{1}{\cosh(kh)}$$

with $k$ from the dispersion relation $\omega^2 = gk\tanh(kh)$.

$\cosh$ grows exponentially, so $K_p$ **falls off exponentially with $kh$**.
Numbers for 18 m water:

| period | wavelength | $kh$ | $K_p$ |
|---|---|---|---|
| 20 s | ~250 m | 0.45 | 0.90 |
| 10 s | ~150 m | 0.75 | 0.75 |
| 5 s | ~39 m | 2.9 | 0.11 |
| 3 s | ~14 m | 8.0 | 0.0007 |

A bottom gauge in 18 m sees swell but is essentially blind to wind sea. This is
not a defect — it is free anti-aliasing, and it is why 1-minute block averages
of bottom pressure are clean. It also means **you cannot recover $H_s$ from
bottom pressure without dividing by $K_p$**, and dividing by 0.0007 amplifies
noise catastrophically, so a high-frequency cutoff is mandatory.

### A footprint trap

Comparing $H_s$ from the buoys against $H_s$ from bottom pressure initially
gave 0.70 m versus 0.40 m and looked like a real discrepancy. It was not: the
buoy figure integrated 3–25 s and the pressure figure 5–32 s. **Different
integration bands are different quantities.** Recomputed over a common band the
disagreement vanished.

The same class of error appeared twice in this project (see also 2-minute means
sampled at 12 min versus 12-minute averages). Compare only commensurable
measurements — same footprint, same averaging, same band.

## 6.3 Establishing the noise floor from colocated instruments

Two independent instruments 3 m apart on the Hydrographer Bank — an RDI
Workhorse ADCP (C05) and a SeaSpider carrying a Sig1000, a CTD and two
MicroRiders — measure the same sea level. Their **difference** contains no
ocean signal, only instrument error, which is the cleanest possible noise
estimate.

Raw difference: **3.63 cm**. Two corrections then apply:

- a **−308 s clock offset** between the instruments, and
- a **−1.25 % gain** difference,

after which the difference falls to **3.08 cm**. Assuming the two instruments
have comparable and independent errors,

$$\sigma_{\text{single}} = \frac{3.08}{\sqrt{2}} = 2.18\ \mathrm{cm}$$

So the post-tidal-fit residual of 3.8–4.5 cm splits roughly into **2.2 cm
instrument** and **3.3 cm ocean** — real non-tidal sea-level variability that no
amount of averaging removes because it is signal, not noise.

That 3.3 cm floor is what makes §6.4 decisive.

## 6.4 Negative result: the eddy is invisible in bottom pressure

From gradient-wind balance (§2.6) the vortex depresses the surface by
**4.1 cm at its center**. The eddy's closest approach to the Pe1 gauge is
1.68 km, where the depression has fallen to about **1 cm**.

Against a 3.3 cm irreducible oceanographic residual, a 1 cm signal is not
detectable. This is not a failure of effort or method — it is a statement about
signal versus floor that could have been made in advance, and now can be.

**Why it is worth writing down**: it tells the next person not to spend a month
looking, and it quantifies what would be needed to succeed — a gauge inside
~500 m of the track, or a much longer record with many eddy passages to
average.

## 6.5 Negative result: drifter GPS height cannot measure sea level

Off by about **100×**. The numbers:

| | scatter |
|---|---|
| burst-mean height, single drifter | 3.3–6.4 m |
| differenced between drifters | 2.0–4.4 m |
| after 3 h smoothing | 1.4–3.9 m |
| **required for the eddy signature** | **1–3 cm** |

### The diagnostic that identifies the cause

Differencing between drifters barely helps. That is the key observation. If the
error were **atmospheric** — ionospheric or tropospheric delay — it would be
strongly common between four buoys a few hundred meters apart, and differencing
would cancel most of it. It does not cancel, so the error is
**platform-specific**: hull multipath and antenna motion on a small buoy in a
seaway.

This is a good example of using the *structure* of an error to identify its
origin rather than guessing.

### Decomposing against the buoys' own velocity

The buoys measure Doppler velocity far more precisely than GPS position.
Comparing the two in the frequency domain shows height is usable near the wave
peak and is **essentially all instrument noise below ~0.02 Hz** — precisely the
band that matters for sea level.

In the 5–32 s band, $H_s$ from height is 1.39 m against 0.65 m from velocity and
0.49–0.55 m from flat-bank bottom pressure. Height is high by more than a
factor of two.

### The usable quality filter

`v_acc`, the receiver's own vertical accuracy estimate, predicts height error
with $r = 0.46$–$0.99$. `pdop` and `numsats` do **not**. Use `v_acc`.

Note that `v_acc` is one of the fields the CORDC-delivered files drop and the
regenerated files keep (§7) — an example of why regenerating from raw was
worth doing.

## 6.6 The ADCP compass rotation

C05's compass is skewed, suspected hard-iron from the iron anchor beneath it.
The correction was established by reproducing the complex transfer function
between the Workhorse and the Sig1000:

$$W = G\, e^{i\phi} S$$

Fitting complex $G e^{i\phi}$ by least squares over $n = 2958$ concurrent
ensembles gives **gain 0.983, angle +21.00°**, against the WAMOS project's
independently derived 0.985 / +21.09° ($n = 2854$). Two analyses agreeing to
0.09° is strong.

The correct absolute correction is $W\cdot e^{-i\,14.4°}$, giving a principal
axis of 113.4°/293.4°. Applying it changes east/west flow labels by 0.0 %,
which is why the shedding-direction result of `PRESSURE_ANALYSIS.md` §8 is
unaffected by it.

### A sign check that was not valid

The rotation sign was once "validated" against the eddy's 266° bearing and
declared backwards. That check was meaningless: C05 sits on an 18 m bank 5–6 km
from an eddy in 600–1900 m of water, so there is no reason the two should
align. **Only the instrument-to-instrument transfer function is a valid
calibration**, because it compares commensurable measurements of the same flow
at the same place.

The stored velocities in `data/adcp/` have the rotation **not** applied,
deliberately, so the raw instrument frame stays visible and the correction is
an explicit step in the analysis.

## References

- Gill, A. E. (1982). *Atmosphere–Ocean Dynamics*. Academic Press. — hydrostatic
  balance and linear wave theory.
- Dean, R. G. and R. A. Dalrymple (1991). *Water Wave Mechanics for Engineers
  and Scientists*. World Scientific. — the pressure response factor $K_p$.
- Bendat, J. S. and A. G. Piersol (2010). *Random Data: Analysis and
  Measurement Procedures*, 4th ed. Wiley. — complex transfer functions and
  coherence.

Full citations with DOIs in `papers/README.md`.
