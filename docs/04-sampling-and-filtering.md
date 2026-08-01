# 4. Sampling and filtering

Waves are the dominant signal in almost every instrument here, and none of the
science is about waves. How you remove them decides whether the answer is
right.

Code: `load_blocks` in `src/eddy_kinematics.py`; `pressure_1min` in
`tools/make_data_subset.py`; `src/pressure_array.py`.

---

## 4.1 The separation of scales that makes this possible

| process | period | where it appears |
|---|---|---|
| wind sea and swell | 3–30 s (peak 4.4 s) | everything |
| infragravity | 30–300 s | pressure gauges |
| **eddy rotation** | **2.9 h** | drifters |
| **tides** | **12.4, 25.8 h** | pressure gauges |
| subtidal / synoptic | days | C05, shedding regime |

The eddy and the tides are two to four decades away from the wave band. That
gap is what makes averaging work: a filter can remove essentially all wave
variance while leaving the signal untouched. Without such a gap none of this
would be separable and the project would need a fundamentally different design.

## 4.2 Block averaging, and why the burst structure matters

The drifters sample at 2 Hz in bursts of 2048 samples. A burst is
**1023.5 s and internally gap-free**; consecutive bursts are separated by
35–202 s while the buoy processes and telemeters, giving a 93.8 % duty cycle.

The code averages **strictly inside a burst**, splitting each 2048-sample burst
into `blocks_per_burst` equal blocks (default 4, so 512 samples = 256 s).

This is not a cosmetic choice. A smoothing window that spans an inter-burst gap
is built from one-sided data, and at burst edges that is wrong by up to
**0.62 m s⁻¹** — comparable to the orbital velocity being measured. Averaging
within a gap-free burst makes that impossible by construction, which is better
than detecting and correcting it.

The 2048 is a *processing* block, not a storage buffer: it is the chunk the buoy
uses for its onboard directional wave spectrum. The evidence is in
`DRIFTER_ANALYSIS.md` §2 — burst duration is rigid (MAD exactly 0) while the
interval varies, which is "collect N, then process," not a clock scheduler.

## 4.3 Never decimate by striding

This is the single most expensive mistake made during this project, and it is
worth understanding precisely because the broken version looked fine.

To reduce a 2 Hz pressure record to 1-minute values there are two obvious
options:

**Striding** — keep every 120th sample:
$$y_k = x_{120k}$$

**Block averaging** — average each group of 120:
$$y_k = \frac{1}{120}\sum_{j=0}^{119} x_{120k+j}$$

Striding is a pure resampling with **no anti-alias filter**. By the sampling
theorem, everything above the new Nyquist frequency $1/120\ \mathrm{s^{-1}} =
8.3$ mHz folds back into the retained band. Swell at 0.1 Hz has enormous
variance, and it lands somewhere in the tidal band determined by the exact
ratio of periods.

Block averaging is a boxcar low-pass followed by resampling. Its transfer
function is

$$H(f) = \frac{\sin(\pi f N \Delta t)}{N\sin(\pi f \Delta t)}$$

which has zeros at every multiple of $1/(N\Delta t)$ and rolls off as $1/f$. It
is not a sharp filter, but it suppresses the wave band by orders of magnitude,
and at tidal frequencies $H \approx 1$ so the signal passes untouched.

### What it cost

| | strided | block-averaged |
|---|---|---|
| variance explained by the tidal fit | 90–97 % | **98.9–99.3 %** |
| post-fit residual | 7.5–14.4 cm | **3.8–4.5 cm** |
| M2 phase spread across the array | **149°** | **1.7°** |

The 149° phase spread is the diagnostic that mattered. Twelve gauges spread
over 11 km in a coherent tidal regime **cannot** have M2 phases differing by
149° — that would be most of a tidal cycle across a few kilometers. The number
was physically impossible, which is what forced the mistake into the open.

**Lesson**: a 90 % variance explained looked perfectly respectable in
isolation. It was the *cross-gauge consistency* that exposed the error. Demand
self-consistency across all data bearing on a situation; any anomaly is
diagnostic.

## 4.4 The other phase trap

Related, and found the same way: phases were initially referenced to **each
record's own start time**. Since the gauges were deployed on different days,
that alone produces a meaningless phase spread — it encodes deployment
schedules, not ocean physics.

All phases must be referenced to a **single common epoch**. The code uses a
module-level `PHASE_EPOCH` so it cannot be forgotten per-record.

Both bugs produced a large phase scatter, and it would have been easy to "fix"
the second and declare victory while the first was still present. They were
separated by checking that the corrected spread was not merely smaller but
*physically plausible* (1.7° over 11 km implies a wave speed, which §5.5 then
turns into an independent test).

## 4.5 Choosing an averaging length

For the drifters, the requirement is:

$$T_{\text{wave}} \ll T_{\text{average}} \ll T_{\text{eddy}}$$

With $T_{\text{wave}} \approx 4.4$ s and $T_{\text{eddy}} = 2.9$ h = 10440 s,
anything from ~30 s to ~1000 s satisfies it. The default of 256 s sits in the
middle on a log scale.

The test that this is not a tuned parameter is in §3.6: the answer moves 0.6 %
across a 16× range of averaging length. When a filter choice barely moves the
answer, the scale separation is doing the work, not the tuning.

## 4.6 Gaps and interpolation

`assemble` puts all four drifters on a common time grid by linear
interpolation, then marks a grid point invalid for a drifter if it is further
than `max_gap_s` (default 1800 s) from any real observation of that drifter.

On this data set the largest actual gap is **421 s**, so no meaningful
interpolation happens — all 355 epochs are genuine four-drifter quadrilaterals.
But note the threshold semantics: since `near` is the distance to the *nearest*
observation, a gap of duration $G$ has its midpoint $G/2$ from each side, so
the default tolerates gaps up to **3600 s** before flagging anything. That is
a third of an eddy rotation. If this pipeline is ever pointed at a record with
real dropouts, lower it.

## 4.7 Aliasing in the pressure records themselves

Even correctly block-averaged, the 1-minute pressure product cannot see
anything above 8.3 mHz. That is deliberate — it is a ~180× size reduction that
loses nothing above the wave band, and every tidal, residual, co-tidal and
gradient result comes from it.

Two short full-rate segments are vendored so the lessons that genuinely need
raw data stay runnable:

- `Pe3_2023-05-22_fullrate.npz` — 16 Hz × 6 h, for the $\cosh(kh)$ wave
  attenuation demonstration (§6.3)
- `HBN_2023-05-22_fullrate.npz` — 2 Hz × 24 h, for reproducing the
  stride-versus-average comparison above

## References

- Emery, W. J. and R. E. Thomson (2014). *Data Analysis Methods in Physical
  Oceanography*, 3rd ed. Elsevier. — filtering and aliasing, chapter 5.
- Oppenheim, A. V. and R. W. Schafer (2010). *Discrete-Time Signal Processing*,
  3rd ed. Pearson. — the boxcar transfer function and decimation.
