# Methods documentation

The mathematics behind every number in `DRIFTER_ANALYSIS.md` and
`PRESSURE_ANALYSIS.md`, derived rather than cited, at the level of an
undergraduate who has had vector calculus and one fluids course.

The two analysis documents tell you **what was found**. These tell you **why
the method is the right one and where it breaks**. Read them when a result
looks surprising, when you want to change a parameter and need to know what it
controls, or when you are asked to defend a number.

## Reading order

| | | |
|---|---|---|
| 1 | [Velocity-gradient kinematics](01-velocity-gradient-kinematics.md) | the tensor, vorticity, divergence, strain, Okubo–Weiss, circulation and Stokes' theorem |
| 2 | [Vortex structure](02-vortex-structure.md) | Lamb–Oseen, radial profiles, Rossby number, gradient-wind balance, why "not solid body" is a real claim |
| 3 | [Uncertainty and bias](03-uncertainty-and-bias.md) | formal vs resampling vs systematic error, what each one can and cannot see, and the null tests |
| 4 | [Sampling and filtering](04-sampling-and-filtering.md) | waves, block averaging, aliasing, and why striding destroys a tidal fit |
| 5 | [Tidal harmonic analysis](05-tidal-harmonic-analysis.md) | least-squares constituents, the Rayleigh criterion, conditioning, and degeneracy |
| 6 | [Pressure and sea level](06-pressure-and-sea-level.md) | hydrostatic inversion, depth attenuation, noise floors, and why the eddy is invisible |
| 7 | [Data formats](07-data-formats.md) | the UBX record layout and the CF conventions the netCDF files satisfy |

## Notation

Consistent throughout, and matching the code:

| symbol | meaning | units |
|---|---|---|
| $u, v$ | eastward, northward velocity | m s⁻¹ |
| $x, y$ | local east, north displacement | m |
| $\zeta$ | relative vorticity, $\partial_x v - \partial_y u$ | s⁻¹ |
| $\delta$ | horizontal divergence, $\partial_x u + \partial_y v$ | s⁻¹ |
| $\sigma_n, \sigma_s$ | normal and shear strain rate | s⁻¹ |
| $\Gamma$ | circulation, $\oint \mathbf{u}\cdot d\boldsymbol{\ell}$ | m² s⁻¹ |
| $f$ | Coriolis parameter, $2\Omega\sin\phi$ | s⁻¹ |
| $Ro$ | Rossby number, $\zeta/f$ | – |
| $R$ | vortex core radius | m |
| $A$ | area enclosed by the drifter polygon | m² |

Sign convention: $\zeta > 0$ is counterclockwise (cyclonic in the northern
hemisphere). The Peleliu vortex is **anticyclonic**, so every vorticity in
this project is negative. Keep the sign; it is physics, not bookkeeping.

**$Ro$ is defined two ways in the literature** and this project uses only the
first: $\zeta/f$ here, versus angular velocity over $f$ — $\omega/f = \zeta/2f$
— used by, among others, Poulain et al. (2023). Theirs is half of ours. Convert
before comparing any published Rossby number with one from this repository.
See §2.6.

## A warning about reading the code

Several results in this project were wrong at some point, and the wrong
versions looked *better* than the right ones — lower residuals, tighter error
bars, cleaner plots. The two analysis documents open with those traps for that
reason. The most instructive:

- A tidal fit that **gained** a corrupted constituent had a **lower** residual
  (3.08 cm vs 4.06 cm) than the correct one. See §5.
- Stride-decimating pressure records raised the variance explained from 90 % to
  99 % *after* the mistake was fixed — the aliased version looked fine in
  isolation. See §4.
- The formal least-squares uncertainty on vorticity is about 5× too small, and
  the leave-one-out uncertainty is blind to the largest error term. See §3.

A number that arrives with a small error bar and no null test is a hypothesis.
