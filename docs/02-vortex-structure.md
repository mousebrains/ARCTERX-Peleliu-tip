# 2. Vortex structure

From a single area-averaged number to the radial profile of the vortex, and
what the numbers mean dynamically.

Code: `radial_profile`, `fit_oseen`, `eddy_center` in `src/eddy_kinematics.py`.

---

## 2.1 Why an area average is not enough

Every vorticity estimate in §1 is an average over whatever area the drifter
polygon happened to span at that moment:

$$\bar{\zeta}(t) = \frac{1}{A(t)}\iint_{S(t)} \zeta \, dA$$

If the vortex were solid-body, $\zeta$ would be constant and the average would
equal the value. It is not, so part of the apparent time variation in
$\bar{\zeta}(t)$ is **the cluster sampling different radii**, not the vortex
changing. Untangling those requires resolving $\zeta$ as a function of radius.

## 2.2 The radial profile

Given a center $\mathbf{c}(t)$ and the eddy's translation velocity
$\mathbf{U}_{\text{trans}}$, every drifter sample becomes one $(r, v_\theta)$ pair:

$$\mathbf{r}' = \mathbf{r} - \mathbf{c}, \qquad \mathbf{u}' = \mathbf{u} - \mathbf{U}_{\text{trans}}$$
$$v_\theta = -u'\sin\theta + v'\cos\theta, \qquad v_r = u'\cos\theta + v'\sin\theta$$

Pooling ~1400 samples over the record and binning in radius resolves the
structure. Two derived quantities follow:

$$\Gamma(r) = 2\pi r\, v_\theta(r) \qquad \text{(circulation within radius } r\text{)}$$
$$\zeta(r) = \frac{1}{r}\frac{d}{dr}\left(r\, v_\theta\right) \qquad \text{(local vorticity)}$$

The second is the axisymmetric form of $\zeta = \partial_x v - \partial_y u$ in
polar coordinates, and it is where the profile shape enters.

**$v_r$ is a check, not an assumption.** For a coherent, non-dispersing vortex
it should be small compared with $v_\theta$. Observed: $|v_r|/|v_\theta| = 0.04$.
That validates both the fitted center (a wrong center manufactures spurious
radial flow) and the claim that the vortex is not falling apart.

Bin statistics use the **median and a MAD-based error**, not mean and standard
deviation, because a handful of wave-contaminated samples would otherwise
dominate. Bins with fewer than 20 samples are dropped.

## 2.3 Solid body or not?

For solid-body rotation $v_\theta \propto r$ and $\zeta$ is constant. The
observed profile is not that:

Binning $|\zeta_\text{circ}|$ into quartiles of cluster scale $\sqrt{A}$:

| cluster scale (bin median) | $|\zeta|$ |
|---|---|
| 419 m | 1.58 × 10⁻³ s⁻¹ |
| 903 m | 0.98 × 10⁻³ s⁻¹ |

$|\zeta|$ falls with the scale being sampled:

$$r = -0.56,\quad 95\%\ \text{CI}\ [-0.63, -0.48],\quad p = 9\times10^{-30},
\quad n = 343$$

**Two cautions on that number, both found by trying to reproduce it.**

*It depends on what you correlate against.* "Cluster scale" and "distance from
the vortex centre" are different variables, and they give different answers on
this dataset — they are themselves only weakly related ($r = +0.27$):

| $|\zeta|$ correlated against | $r$ | 95 % CI |
|---|---|---|
| cluster scale $\sqrt{A}$ | **−0.56** | [−0.63, −0.48] |
| radius from the fitted centre | **−0.37** | [−0.46, −0.27] |

Both are decisively non-zero, so the *conclusion* does not turn on the choice,
but the value does. Quote the definition with the number, and do not quote a
second decimal: the CI is ±0.07 wide.

*Earlier versions of this document reported −0.59.* That is inside the CI
above, but it is not reproducible from the current pipeline — the closest
configuration returns −0.583 at $Q_{\min} = 0.20$ rather than the 0.10 actually
used. The binned values 1.58/0.98 **are** exactly reproducible. Treat −0.56 as
the number.

### Why this needed a null test

A correlation of −0.56 is suggestive, not conclusive, because the estimator
itself could manufacture it: the quadrature error of §1.5 grows with cluster
size, so a bigger cluster might report a smaller $|\zeta|$ *even for a
solid-body vortex*. If so, the conclusion would be an artifact.

`tests/test_synthetic_recovery.py` runs a **true solid-body vortex** through the
identical machinery. Advecting a single cluster gives:

| field | cluster | corr($|\zeta|$, radius) |
|---|---|---|
| solid body | 387 m | **+0.06** |
| solid body | 774 m | **−0.05** |
| Lamb–Oseen | 934 m | −0.56 |
| Lamb–Oseen | 1105 m | −0.74 |

**That table does less work than it appears to, for two reasons.**

First, a cluster advected through a *solid-body* field is rigid: its size never
varies (CV ≈ 10⁻⁵) and $|\zeta|$ is exactly constant. The two solid-body rows
correlate one constant against another, so they cannot exercise the
size-dependent artifact they are meant to exclude.

Second, those correlations are against **radius from the centre**, while the
headline −0.56 is against **cluster scale** — different footprints. On the
matched variable the synthetic Lamb–Oseen gives $r = +0.38$ and $+0.41$, the
*opposite* sign, because in that synthetic the cluster shrinks as it moves
outward ($r = -0.78$) whereas the real cluster grows slightly ($r = +0.27$).
Comparing the real cluster-scale number against the synthetic radius numbers
compares two different things.

**Test 2b is the one that settles it.** Pool solid-body clusters over a
232–1238 m range of sizes — the comparison the claim actually rests on:

$$\text{corr}(|\zeta|, \sqrt{A}) = +0.012, \qquad
\frac{\max|\zeta| - \min|\zeta|}{\overline{|\zeta|}} = 0.0000\,\%$$

The estimator returns $1.200000\times10^{-3}$ at *every* cluster size from 232
to 1238 m. There is no size-dependent artifact to worry about, so the observed
−0.56 is diagnostic. **The claim survives** — on a stronger test than the one
originally offered.

The rotation shortfall in §1.8 is the same physics seen another way: the
constellation turned 6.86 revolutions where the *core* vorticity predicts 8.83,
because the outer drifters orbit more slowly than the core vorticity implies.
Two independent diagnostics, one conclusion.

### Corroboration from an unrelated vortex

Poulain et al. (2023) hit the same bias in the Cyprus Gyre and state it
plainly: the drifter estimate "is an overestimate because the vorticity, in
absolute value, always decreases with increasing distance from the gyre
center," and the solid-body inversion holds "only ... for a small radius less
than 10–20 km." Same conclusion, a vortex four orders of magnitude away in
Rossby number and two in size — so this is a property of the *estimator*, not
of Palau.

Their wording and ours point in opposite directions, which is worth being
careful about. Both say the same thing. The Stokes average over a disc of
radius $r$,

$$\bar{\zeta}(r) = \frac{\Gamma(r)}{\pi r^2} = 2\omega(r)$$

lies between $\zeta(0)$ and $\zeta(r)$ whenever $|\zeta|$ decreases outward.
So it **over**estimates the local vorticity at the orbit (Poulain's framing)
and **under**estimates the core (ours). One inequality, read from two ends.

The cost is quantified here. Inverting the *observed* constellation rotation
under a solid-body assumption:

$$P_{\text{obs}} = \frac{25.4\ \mathrm{h}}{6.86\ \mathrm{rev}} = 3.70\ \mathrm{h}
\quad\Longrightarrow\quad
\zeta = \frac{4\pi}{P_{\text{obs}}} = 9.4\times10^{-4}\ \mathrm{s^{-1}}$$

against $1.19\times10^{-3}\ \mathrm{s^{-1}}$ from four-drifter circulation —
**21 % low**, matching the 22 % shortfall in §1.8. Neither number is wrong;
they average over different radii, because the orbital loop is larger than the
instantaneous cluster polygon. Both sit on the profile above: $1.19\times10^{-3}$
near the ~640 m median cluster scale, $9.4\times10^{-4}$ out past ~1 km.

**The practical consequence.** Any method that recovers $\zeta$ from a *single*
trajectory — wavelet ridge, rotary Fourier, complex demodulation — can only
reach it through $\zeta = 4\pi/P$, and so inherits this 21 %. That is the
reason the four-drifter circulation estimate is primary and no single-drifter
spectral method was adopted; see entry 7 of `papers/README.md`.

## 2.4 Finding the center

In a frame moving with the eddy at velocity $\mathbf{c}$, the center is where
the velocity vanishes. From the affine fit,

$$\mathbf{A}(\mathbf{x} - \bar{\mathbf{x}}) = \mathbf{c} - \mathbf{U}_0 \quad \Longrightarrow \quad \mathbf{x}_{\text{center}} = \bar{\mathbf{x}} + \mathbf{A}^{-1}(\mathbf{c} - \mathbf{U}_0)$$

This is an **elliptic critical point** only when $\mathbf{A}$ has complex
eigenvalues — equivalently when $OW < 0$ (§1.7). Where $OW \geq 0$ the critical
point is a saddle and the code returns NaN rather than a meaningless "center."

### The degeneracy that cannot be avoided

At a single instant, a uniform background flow and a displacement of the vortex
center are **exactly** degenerate. For a vortex of strength $\omega$ centered at
$\mathbf{c}$,

$$\mathbf{u} = \omega\,\hat{z}\times(\mathbf{x} - \mathbf{c}) = \underbrace{\omega\,\hat{z}\times\mathbf{x}}_{\text{vortex at origin}} \underbrace{-\ \omega\,\hat{z}\times\mathbf{c}}_{\text{a constant}}$$

The second term is a uniform flow. No instantaneous velocity field can
distinguish "vortex here plus background current" from "vortex over there."
This is not a limitation of the estimator; it is a property of the field.

The degeneracy is broken by **time**: averaging the drifter-mean velocity over
several orbital periods (default 6 h, ~2 periods) cancels the rotational part
and leaves the translation. Every center estimate is therefore conditional on
that low-pass choice, and if a center track looks wrong this is the first thing
to suspect.

### Two gates

$\mathbf{A}^{-1}$ scales as $1/|\zeta|$, so where rotation is weak or poorly
determined the inversion amplifies a small velocity mismatch into an enormous
displacement. Hence:

1. **Signal-to-noise**: require $|\zeta| > 3\sigma_\zeta$ before inverting.
2. **Displacement**: reject a center more than 3 cluster-widths away, since
   drifters all sitting on one side of a distant point do not constrain it.

Together these pass 89 % of windows on the real data. On the synthetic test the
recovered center lands within ~1.1 km of truth (median), which is comparable to
the core radius — a reminder that the center is the **weakest** link in the
chain, and why the radial profile is also computed with the center supplied
independently as a check.

## 2.5 The Lamb–Oseen vortex

The Lamb–Oseen vortex is the exact solution of the Navier–Stokes equations for
a line vortex decaying under viscosity. Its azimuthal velocity is

$$v_\theta(r) = \frac{\Gamma}{2\pi r}\left(1 - e^{-r^2/R^2}\right)$$

with vorticity

$$\zeta(r) = \frac{\Gamma}{\pi R^2}\, e^{-r^2/R^2}$$

Two limits explain why it is the natural model:

- $r \ll R$: expanding, $v_\theta \to \Gamma r/(2\pi R^2)$, linear in $r$ —
  **solid-body core**.
- $r \gg R$: $v_\theta \to \Gamma/(2\pi r)$ — **irrotational free vortex**.

So it interpolates between the two idealizations with a single shape parameter
$R$, the core radius, defined as where the transition happens. It is a
two-parameter family $(\Gamma, R)$, which is about as much as ~12 binned points
can support.

### How the fit is done

Brute force over $R$, with $\Gamma$ solved **linearly** at each $R$. Because
$v_\theta$ is linear in $\Gamma$ at fixed $R$, the inner problem has a closed
form:

$$\Gamma(R) = \frac{\sum_i w_i b_i(R)\, v_i}{\sum_i w_i b_i(R)^2}, \qquad b_i(R) = \frac{1 - e^{-r_i^2/R^2}}{2\pi r_i}$$

Scanning $R$ over a 400-point grid and taking the minimum residual means there
is **no starting guess and no local-minimum risk** — a real advantage over
throwing both parameters at a generic optimizer. Weights are $1/\sigma^2$ from
the bin MADs.

The code checks that the returned $R$ is interior to the search range. On the
real profile it returns 1206 m from a 95–3818 m search, comfortably interior.
An earlier configuration *did* return the lower bound exactly, which is why
this check exists: a boundary hit means the fit is unconstrained and the value
is meaningless.

### The fit is the fragile step

Two correct implementations differing only in percentile convention and
moving-average edge handling returned $R = 1040$ m and $R = 1206$ m. Repeating
over a grid of bin counts and radius cutoffs spans 1161–1286 m. The Python and
MATLAB ports differ: 1206 m vs 1116 m.

**Quote the interval, never the point estimate.** The medians in §1 are robust;
this single global fit to ~12 binned points is not. The MATLAB port adds a
1000-replicate bootstrap for exactly this reason, giving 95 % CI [1005, 1605] m
for $R$ and [−9051, −4591] m² s⁻¹ for $\Gamma$. Both point estimates fall inside
the other's interval.

## 2.6 Dynamical interpretation

### Rossby number

$$Ro = \frac{\zeta}{f}, \qquad f = 2\Omega\sin\phi$$

At 6.99 °N, $f = 2 \times 7.292\times10^{-5} \times \sin(6.99°) = 1.776\times10^{-5}\ \mathrm{s^{-1}}$
— small, because Palau is near the equator. With $\zeta = -1.19\times10^{-3}$,

$$Ro = -67$$

$|Ro| \gg 1$ means rotation of the *vortex* overwhelms rotation of the *planet*.
The flow is **cyclostrophic**, not geostrophic: centrifugal force balances the
pressure gradient and the Coriolis term is a 1.5 % correction. This is a
submesoscale island wake, not a mesoscale eddy, and geostrophic intuition does
not apply to it.

**Independently corroborated at this headland.** Johnston et al. (2019)
measured the same eddies at the Peleliu tip three years earlier, from moored
ADCPs and shipboard survey rather than drifters, and report $Ro$ "reaching 80
and 65" for ~1–2 km diameter anticyclonic wake eddies at the separation point.
They use $Ro = \zeta/f$, the same convention as here. Our $-67$ falls inside
that range — a different instrument class, the same headland, the same number.

One caveat on the comparison: their quoted "$Ro \sim 30$" scaling divides a
1 m s⁻¹ velocity *difference* by a 2 km diameter, which is a one-sided shear
rather than $\zeta$, so their scaling and mooring numbers differ by about a
factor of two among themselves. The circulation/Stokes estimate here is the
better-defined quantity.

**Convention warning.** A second definition of $Ro$ is in wide use — the
*angular velocity* over $f$, i.e. $\omega/f = \zeta/2f$, **half** of the one
above. Poulain et al. (2023) use it and say so in their Appendix ("inertial
motions have an angular velocity of $f$, a vorticity of $2f$ ... and a Rossby
number ... equal to 1"). Their Cyprus Gyre core, quoted as "$Ro \sim 0.25$,
mean vorticity $0.5f$", is $Ro = 0.5$ here. Halve ours or double theirs before
any comparison; the factor of two is silent and neither paper is wrong.

### Is it inertially stable?

A reasonable worry at $Ro = -67$: an anticyclone with $|\zeta| \gg f$ has
$f + \zeta < 0$, which is the textbook signature of inertial instability. It is
not unstable, and the reason is worth recording because the naive test gives
the wrong answer.

The correct criterion for a circular vortex keeps both factors — the
generalised Rayleigh discriminant, i.e. Rayleigh's circulation theorem extended
to a rotating frame (Kloosterziel and van Heijst 1991, §4):

$$\Phi = (f + \zeta)\left(f + \frac{2v_\theta}{r}\right) > 0
\quad \text{for stability}$$

For an anticyclone $\zeta$ and $v_\theta/r$ are both negative, so the two
factors flip sign **together** and the product stays positive. With
$\omega = v_\theta/r = -4.72\times10^{-4}\ \mathrm{s^{-1}}$ from the observed
3.70 h orbit,

$$\Phi = (-1.17\times10^{-3})(-9.26\times10^{-4}) = +1.09\times10^{-6} > 0$$

Stable. The associated frequency $\sqrt{\Phi} = 1.04\times10^{-3}\ \mathrm{s^{-1}}$
corresponds to a period of **1.68 h** — set almost entirely by the vortex, and
59× faster than the 98.3 h local inertial period.

This is also why the near-inertial machinery in the drifter literature does not
transfer here. The linearised effective frequency
$f_e = \sqrt{f^2 + f\zeta} \approx f + \zeta/2$ (Kunze 1985; Poulain et al.
2023 Eq. 1) assumes $|\zeta|/f \ll 1$. At $|\zeta|/f = 67$ the radicand is
$-2.08\times10^{-8}$ and the expression is imaginary — a **linearisation
failure, not a physical instability**, as $\Phi > 0$ shows. Combined with an
inertial period 3.9× longer than the entire 25.4 h record, "near-inertial" is
not a usable category for this dataset.

**What this does not establish.** $\Phi > 0$ rules out *inertial/centrifugal*
instability only. Barotropic (shear) instability is a separate criterion, and
Kloosterziel and van Heijst's central result is that the two senses behave
differently — the **anticyclone is the fragile case**, showing "rather explosive
instability behaviour" and splitting into dipoles, where the cyclone decays
gradually into a tripole. Nothing here has tested our vortex against a
shear-instability criterion. Its 25.4 h coherence ($OW < 0$ in 100 % of windows)
is an **observation, not a prediction**, and the record ends while the vortex is
still coherent — so we do not know how it died.

### Rotation period

Treating the core as solid-body with $\zeta = 2\omega$,

$$T = \frac{2\pi}{\omega} = \frac{4\pi}{|\zeta|} = 2.9\ \text{h}$$

Short enough that the 25.4 h record contains ~8 rotations, which is what makes
the constellation-rotation cross-check meaningful.

### Gradient-wind balance and the surface depression

For a circular vortex the radial momentum balance is

$$\frac{v_\theta^2}{r} + f v_\theta = g\frac{\partial \eta}{\partial r}$$

The first term is centrifugal, the second Coriolis. At $Ro = -67$ the ratio
$v_\theta/(fr) $ makes centrifugal beat Coriolis about **24:1**, so this is
essentially cyclostrophic balance. Integrating inward gives a surface
depression at the center of about **4.1 cm** — which is why §6 concludes the
eddy is invisible in bottom pressure: 4.1 cm at the center, ~1 cm at the
nearest gauge's 1.68 km approach, against a 3.3 cm irreducible oceanographic
residual.

### Circulation budget

$\Gamma \approx -5400$ to $-6100\ \mathrm{m^2 s^{-1}}$ requires roughly 2–3×
more forcing than one M2 half-cycle at the measured 0.64 m s⁻¹ tip current
supplies. Either the flow past the tip exceeds the array-averaged value, or
circulation accumulates over several tidal cycles. This is open work; the
Thompson ADCP spot measurements could discriminate.

## References

- Lamb, H. (1932). *Hydrodynamics*, 6th ed. Cambridge University Press, §334a.
- Oseen, C. W. (1912). Über die Wirbelbewegung in einer reibenden Flüssigkeit.
  *Ark. Mat. Astro. Fys.* **7**, 14–21.
- Saffman, P. G. (1992). *Vortex Dynamics*. Cambridge University Press.
- Holton, J. R. and G. J. Hakim (2013). *An Introduction to Dynamic
  Meteorology*, 5th ed. — gradient-wind balance, §3.2.5.
- Kunze, E. (1985). Near-inertial wave propagation in geostrophic shear.
  *J. Phys. Oceanogr.* **15**, 544–565 — the linearised $f_e$ that §2.6 shows
  does not apply here.
- Kloosterziel, R. C. and G. J. F. van Heijst (1991). An experimental study of
  unstable barotropic vortices in a rotating fluid. *J. Fluid Mech.* **223**,
  1–24 — the modified Rayleigh criterion used in §2.6, and the caveat that
  anticyclones are the fragile case under *shear* instability.
- Johnston, T. M. S. *et al.* (2019). Energy and momentum lost to wake eddies
  and lee waves ... at Peleliu, Palau. *Oceanography* **32**(4), 110–125 —
  same headland, $Ro$ 65–80 from moorings (§2.6).
- Poulain, P.-M. *et al.* (2023). Drifter observations of surface currents in
  the Cyprus Gyre. *Front. Mar. Sci.* **10**, 1266040 — independent statement
  of the Stokes area-average bias (§2.3), and the source of the $\omega/f$
  Rossby convention (§2.6).

Full citations with DOIs in `papers/README.md`.
