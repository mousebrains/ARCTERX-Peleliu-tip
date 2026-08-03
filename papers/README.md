# References for the Peleliu tip-vortex drifter analysis

Source literature behind `eddy_kinematics.py`, `eddy_analysis.py`,
`eddy_kinematics_drifters.m` and the raw-data readers `mwb_dat.py` / `mwb_nc.py`.

Each entry records the full citation, a DOI, the local PDF, and — the part that
matters when you come back to this in six months — **which function actually
implements it**.

## Status

| # | Short cite | PDF | Citation checked against |
|---|-----------|-----|--------------------------|
| 1 | Okubo & Ebbesmeyer 1976 | ✓ `OkuboEbbesmeyer1976_DSR_drogue_vorticity.pdf` | title page |
| 2 | Molinari & Kirwan 1975 | ✓ `MolinariKirwan1975_JPO_differential_kinematics.pdf` | title page |
| 3 | Okubo 1970 | ✓ `Okubo1970_DSR_velocity_singularities.pdf` | title page |
| 4 | Weiss 1991 | ✓ `Weiss1991_PhysicaD_enstrophy_transfer.pdf` | title page |
| 5 | Efron & Gong 1983 | ✓ `EfronGong1983_AmStat_bootstrap_jackknife.pdf` | title page (JSTOR cover) |
| 6 | Saffman 1992 | ⏳ expected | **not yet verified** |
| 7 | Poulain et al. 2023 | ✓ `Poulain2023_FrontMarSci_cyprus_gyre_wavelet_ridge.pdf` | title page |
| 8 | Johnston et al. 2019 (Peleliu) | ✓ `Johnston2019_Oceanography_peleliu_wake_eddies_lee_waves.pdf` | title page |
| 9 | Rudnick et al. 2019 | ✓ `Rudnick2019_Oceanography_vorticity_flow_past_island.pdf` | title page |
| 10 | St. Laurent et al. 2019 | ✓ `StLaurent2019_Oceanography_palau_wake_turbulence_vorticity.pdf` | title page |
| 11 | Siegelman et al. 2019 | ✓ `Siegelman2019_Oceanography_palau_near_inertial_surface.pdf` | title page |
| 12 | Johnston et al. 2019 (FLEAT) | ✓ `Johnston2019_Oceanography_FLEAT_program_overview.pdf` | title page |
| 13 | Essink et al. 2022 | ✓ `Essink2022_JTECH_drifter_cluster_kinematics.pdf` | title page |
| 14 | Kloosterziel & van Heijst 1991 | ✓ `KloosterzielVanHeijst1991_JFM_unstable_barotropic_vortices.pdf` | title page |

Every volume/issue/page figure below was read off the article's own title page
— not from a publisher web record. Saffman (1992) is a book and is still
outstanding; its details remain unverified until the copy arrives.

Entries **1–5** are implemented. Entry **7** is consulted but deliberately
**not** implemented — see the entry for the arithmetic behind that decision.
Entries **8–12** are the site and phenomenon literature: Palau, and in one case
the Peleliu tip itself. Entries **13–14** are method and dynamics support.

**Entries 8–12 are open access** (*Oceanography*, The Oceanography Society) and
entry 14 is a copy posted by its author, so unlike 1–5 these were downloadable
directly. See "A note on redistribution" for what that does and does not
change.

Several important papers are **still missing behind paywalls** — see
[Wanted](#wanted--behind-a-paywall) at the end. One of them, Zeiden et al.
(2022), is close enough to this analysis that it should be read before any
write-up.

`./fetch_papers.sh --list` reports what is missing at any time.

---

## The references in detail

### 1. Okubo & Ebbesmeyer (1976) — the four-drogue cluster method

> Okubo, A., and C. C. Ebbesmeyer, 1976: Determination of vorticity, divergence,
> and deformation rates from analysis of drogue observations. *Deep-Sea
> Research* **23**, 349–352. Pergamon Press.
> doi:[10.1016/0011-7471(76)90875-5](https://doi.org/10.1016/0011-7471(76)90875-5)
> — received 6 Feb 1975, accepted 28 Mar 1975; Contribution No. 133, Marine
> Sciences Research Center, SUNY Stony Brook.

The foundational paper here. Its abstract states the case exactly: with four or
more drogues followed simultaneously, "not only the mean flow, dispersion and
eddy diffusivities but also the field of mean vorticity, divergence, and
deformation rates can be determined as functions of time."

**Its Eqs. (1)–(2) are precisely the model implemented in `fit_gradient()`** —
each drogue's velocity expanded in a Taylor series about the cluster centroid,

```
u_i(k) = u_bar(k) + du/dx [x_i(k) - x_bar(k)] + du/dy [y_i(k) - y_bar(k)] + u_i''(k)
v_i(k) = v_bar(k) + dv/dx [x_i(k) - x_bar(k)] + dv/dy [y_i(k) - y_bar(k)] + v_i''(k)
```

solved by linear regression. Confirming that the code matches the canonical
formulation, term for term, is worth more than any amount of re-derivation.

**Historical note.** Footnote § on p. 349 records that *after* submission the
authors became aware of Molinari & Kirwan (1975), "who independently developed
the idea and analysis similar to ours." Entries 1 and 2 are therefore
independent co-discoveries, not derivative works — cite both.

**Implements:** `eddy_kinematics.fit_gradient()` / `fitGradient()`.

### 2. Molinari & Kirwan (1975) — least-squares differential kinematics

> Molinari, R., and A. D. Kirwan, Jr., 1975: Calculations of differential
> kinematic properties from Lagrangian observations in the western Caribbean
> Sea. *Journal of Physical Oceanography* **5**, 483–491.
> doi:[10.1175/1520-0485(1975)005<0483:CODKPF>2.0.CO;2](https://doi.org/10.1175/1520-0485(1975)005%3C0483:CODKPF%3E2.0.CO;2)
> — manuscript received 12 Dec 1974, revised 4 Mar 1975.

Computes divergence, vorticity, shear and normal deformation from drifter
clusters by **two independent analyses** and compares them — the same
methodological stance taken here.

Their abstract contains a warning that landed on us directly: the time series
are "fairly smooth when the drifters were moving in the Yucatan Current.
Otherwise, the time series are ragged with frequent changes in sign," which they
attribute to "small values of the shear rates relative to random observational
errors or small-scale turbulent processes." That is exactly the failure mode
behind our LSQ-vs-circulation scatter, and it is why the leave-one-out spread
(17 %) rather than the formal error (5 %) is the number to quote.

**Implements:** `fit_gradient()` / `fitGradient()`, and motivates the
estimator-comparison diagnostics.

### 3. Okubo (1970) — the Okubo half of Okubo–Weiss

> Okubo, A., 1970: Horizontal dispersion of floatable particles in the vicinity
> of velocity singularities such as convergences. *Deep-Sea Research* **17**,
> 445–454. Pergamon Press.
> doi:[10.1016/0011-7471(70)90059-8](https://doi.org/10.1016/0011-7471(70)90059-8)
> — received 5 Sep 1969; Contribution No. 143, Chesapeake Bay Institute, Johns
> Hopkins University.

Approximates the velocity field near a singularity as **linear**, then
decomposes it into "the rates of stretching deformation and shearing
deformation, vorticity, and divergence (or negative convergence)" — the same
four quantities the gradient tensor is split into here. It also shows that a
ring of particles around a singularity deforms into an ellipse whose area grows
or shrinks with the sign of the divergence, which is the conceptual basis for
the `d(ln A)/dt` divergence estimator.

**Implements:** the strain/rotation decomposition and the `OW` term in
`fit_gradient()`; conceptually underpins `divergence_from_area()`.

### 4. Weiss (1991) — the Weiss half

> Weiss, J., 1991: The dynamics of enstrophy transfer in two-dimensional
> hydrodynamics. *Physica D* **48**, 273–294. North-Holland.
> doi:[10.1016/0167-2789(91)90088-Q](https://doi.org/10.1016/0167-2789(91)90088-Q)
> — received 21 Sep 1990, accepted 15 Oct 1990, communicated by U. Frisch.

Shows that vorticity gradients grow exponentially where squared strain exceeds
squared vorticity, and decay in the opposite case. With Okubo (1970) this gives

```
OW = sigma_n^2 + sigma_s^2 - zeta^2      OW < 0  =>  coherent vortex
```

On these data `OW < 0` in **100 %** of windows across the full 25.4 h — the
primary evidence that the four drifters stayed inside one coherent vortex rather
than drifting through unrelated strain.

**Bibliographic note.** The editor's footnote on p. 273 states the manuscript
"was written in 1981 as a La Jolla Institute preprint and never published in the
open literature" before this printing. That is why some papers cite "Weiss
(1981)" for the same result. Both refer to this work; **cite the 1991 Physica D
version**, which is the version of record.

**Implements:** the `OW` term in `fit_gradient()` / `fitGradient()`, and the
`OW < 0` gate in `eddy_center()` / `eddyCenter()`.

### 5. Efron & Gong (1983) — the jackknife error bar

> Efron, B., and G. Gong, 1983: A leisurely look at the bootstrap, the
> jackknife, and cross-validation. *The American Statistician* **37**(1),
> 36–48. Taylor & Francis for the American Statistical Association.
> doi:[10.1080/00031305.1983.10483087](https://doi.org/10.1080/00031305.1983.10483087)
> · JSTOR [2685844](https://www.jstor.org/stable/2685844)

An expository treatment of resampling error estimates. We use the leave-one-out
form: four drifters give four distinct triangles, and for a strictly linear
velocity field all four would return identical vorticity, so their spread
measures noise *plus* curvature of the flow without assuming a noise model.

This matters concretely — the jackknife spread came out **17 %** of |ζ| against
a formal least-squares error of **5 %**. The formal error is roughly three times
too optimistic because it cannot know the flow is curved across the cluster.
**Quote the jackknife number.**

**Implements:** `circulation_jackknife()` / `circulationJackknife()`.

### 6. Saffman (1992) — the Lamb–Oseen vortex ⏳ *awaiting copy*

> Saffman, P. G., 1992: *Vortex Dynamics*. Cambridge University Press.
> ISBN 978-0-521-42058-7

Standard text; the Lamb–Oseen (viscously spreading line) vortex is the profile
fitted to the pooled radial data:

```
v_theta(r) = Gamma/(2 pi r) * (1 - exp(-r^2/R^2))
zeta(r)    = Gamma/(pi R^2) * exp(-r^2/R^2)
```

Fitted values here: R = 1206 m, Γ = −6114 m² s⁻¹. Caveat from the analysis: the
observed profile has not clearly turned over by the largest radius the drifters
sampled (~1900 m), so **R is far better constrained than Γ**, which extrapolates
past the data.

**Citation details are unverified** — page count and edition were quoted from
memory. Check them against the copy when it arrives, and pin the specific
section for the Lamb–Oseen profile so the reference points somewhere useful.

**Implements:** `fit_oseen()` / `fitOseen()`.

### 7. Poulain et al. (2023) — wavelet ridge analysis, *evaluated and not adopted*

> Poulain, P.-M., M. Menna, E. Mauri, A. Pirro, D. R. Hayes, and H. Gildor,
> 2023: Drifter observations of surface currents in the Cyprus Gyre.
> *Frontiers in Marine Science* **10**, 1266040.
> doi:[10.3389/fmars.2023.1266040](https://doi.org/10.3389/fmars.2023.1266040)
> — received 24 Jul 2023, accepted 11 Oct 2023, published 27 Oct 2023.
> Open access, CC BY. Author names are given as initials because that is how
> the title page prints them; do not expand them without a source.

Wavelet ridge analysis (Lilly & Gascard 2006; jLab v1.7.1) applied to SVP
drifters trapped in the Cyprus Gyre, separating the gyre rotation from
near-inertial oscillations. Structurally the closest published analysis to
ours — anticyclonic vortex, drifters as the only in-situ instrument,
vorticity from Stokes' theorem — which is why it was read carefully and why
the negative verdict is recorded here rather than left implicit.

**The method is not adopted.** Wavelet ridges buy one thing over Fourier:
time-resolved instantaneous frequency and amplitude. Poulain needs that — 2–6
month records, gyre period wandering 4→11 days, NIO amplitude intermittent
between 0 and 39 cm/s. Our 25.4 h record cannot use it:

| | Poulain | here |
|---|---|---|
| record | 168 d (drifter 3469260) | 25.4 h |
| cycles at the vortex frequency | 15–42 | **6.86** |
| usable after edge trim (γ=3, β=4 Morse, footprint ≈ P/π = 1.1 cycles per end) | ~13–40 | **~4.6** |
| second band to separate | NIO at ~1 d | none resolvable |

At ~4.6 usable cycles any frequency drift a ridge reports falls inside the
wavelet's own bandwidth — that Morse wavelet has ~70 % relative frequency
resolution (their Appendix: 0.7 cpd at 1 cpd). Rotary Fourier or complex
demodulation at the known 3.7 h frequency is cheaper and better conditioned.
Note also that the feature which would fix the ±90° sign degeneracy in
`explore_drifter_paths.m` is the **rotary** decomposition, not the wavelet;
rotary Fourier spectra supply it, and Poulain use them in their Fig. 8.

**The near-inertial half is inapplicable, not merely marginal.** At 6.99 °N,
f = 1.78 × 10⁻⁵ s⁻¹ and the inertial period is **98.3 h** — the 25.4 h record
spans 0.26 of one cycle. Their Eq. (1) is also undefined on our numbers:

```
fe = sqrt(f^2 + f*zeta)     requires zeta >= -f
zeta = -1.19e-3 = -67 f     f^2 + f*zeta = -2.08e-8  ->  imaginary
```

This is a linearisation failure, **not** an instability. Eq. (1) assumes
|ζ|/f ≪ 1; ours is 67. The unlinearised criterion keeps both factors
negative, so the vortex is inertially stable:

```
(f + zeta)(f + 2V/r) = (-1.17e-3)(-9.26e-4) = +1.09e-6 > 0
effective frequency = 1.04e-3 s^-1  ->  period 1.68 h
```

That effective period is set by the vortex, not by f — 59× faster than the
local inertial period. "Near-inertial" is not a meaningful category here, so
§3.3, §3.4, Table 2 and the damped slab model have no counterpart in this
analysis.

**What it is good for anyway — three things:**

1. **Independent confirmation of our Stokes-theorem trap.** They write that
   the drifter estimate "is an overestimate because the vorticity, in absolute
   value, always decreases with increasing distance from the gyre center," and
   that solid-body inversion holds "only ... for a small radius." That is
   exactly `DRIFTER_ANALYSIS.md` §5 *"each ζ is an area average"*, arrived at
   independently. Our numbers: the observed −6.86 rev in 25.4 h gives
   P = 3.70 h, so the solid-body inversion ζ = 4π/P returns
   9.4 × 10⁻⁴ s⁻¹ against 1.19 × 10⁻³ s⁻¹ from four-drifter circulation —
   **21 % low**, consistent with the 23 % shortfall in §4.
2. **The orbital-speed-versus-radius presentation** (their Figs. 7 and 13,
   with constant-angular-velocity reference lines at f, f/4, f/10, f/20) is
   the same diagnostic as our radial profile and Lamb–Oseen fit, and is a
   clearer way to show departure from solid body than a fitted curve alone.
3. **A convention trap, if their numbers are ever quoted alongside ours.**
   Poulain define the Rossby number as *angular velocity* over f (Appendix:
   "inertial motions have an angular velocity of f, a vorticity of 2f ...
   and a Rossby number ... equal to 1"), i.e. ζ/2f. This project uses ζ/f.
   Their gyre core at "Ro ~ 0.25, mean vorticity 0.5 f" is Ro = 0.5 in our
   convention. **Factor of two — halve ours, or double theirs, before
   comparing.**

For scale: their gyre is Ro ≈ 0.5 (our convention) at ~15–30 km radius;
ours is Ro = −67 at ~1.2 km. Same qualitative object, four orders of
magnitude apart in Rossby number and two in size.

**Implements:** nothing. Consulted, evaluated, method declined; cited for the
Stokes-overestimate confirmation and the convention warning.

### 8. Johnston et al. (2019) — wake eddies at the Peleliu tip ★ *our site*

> Johnston, T. M. S., J. A. MacKinnon, P. L. Colin, P. J. Haley Jr.,
> P. F. J. Lermusiaux, A. J. Lucas, M. A. Merrifield, S. T. Merrifield,
> C. Mirabito, J. D. Nash, C. Y. Ou, M. Siegelman, E. J. Terrill, and
> A. F. Waterhouse, 2019: Energy and momentum lost to wake eddies and lee waves
> generated by the North Equatorial Current and tidal flows at Peleliu, Palau.
> *Oceanography* **32**(4), 110–125.
> doi:[10.5670/oceanog.2019.417](https://doi.org/10.5670/oceanog.2019.417)
> — open access.

**The most directly relevant paper in this directory.** Same island, same
headland, same phenomenon, three field seasons earlier. It reports the eddies
we measured, from moorings, shipboard survey and a model rather than drifters.

Their numbers against ours:

| | Johnston et al. 2019 | this analysis |
|---|---|---|
| eddy diameter at the separation point | ~1 km ("small and intense"), 2 km in surveys B/C | core radius 1.1–1.2 km |
| Rossby number | "Ro ~ 30"; moored ADCP time series "reaching 80 and 65" | **−67** |
| sense | anticyclonic expected and observed with westward flow | anticyclonic |
| shedding | every tidal cycle; intrinsic timescale ~6 h | one event, 25.4 h coherent |

**They use Ro = ζ/f, the same convention as this project** — stated explicitly
as "Ro = ζ/f = (∂ₓv − ∂ᵧu)/f". No factor-of-two conversion, unlike entry 7.
Our −67 sits inside their moored range of 65–80. Treat that as independent
corroboration of the headline number from a different instrument class, with
one caveat: their "Ro ~ 30" scaling divides a 1 m s⁻¹ velocity *difference* by a
2 km diameter, which is a one-sided shear rather than the full ζ, so their
scaling-based and mooring-based numbers differ by about a factor of two among
themselves. Ours is a closed-contour circulation estimate and is the
better-defined quantity.

**It bears directly on the circulation budget in `DRIFTER_ANALYSIS.md` §7.**
That open question asks whether Γ ≈ −5400 to −6100 m² s⁻¹ needs more forcing
than one M2 half-cycle at the array-averaged 0.64 m s⁻¹ tip current supplies —
2–3× more — and offers two horns: the tip flow exceeds the array average, or
circulation accumulates over several cycles. This paper supports the **first**:
it reports the flow "intensifies" at the south point as it is constrained
around topography, and confirms "the velocity difference across the eddy
exceeds 1 m s⁻¹". A local ~1 m s⁻¹ against our array-averaged 0.64 m s⁻¹ is a
factor ~1.6 on Γ ~ UL. Combined with a 6 h shedding timescale against a 6.2 h
M2 half-cycle, that plausibly closes the gap — **hypothesis, not settled**, and
the discriminating measurement is still the Thompson ADCP.

Two dimensionless numbers worth carrying over, both new to this project:

```
island wake parameter (Wolanski et al. 1984)   Ref = H/(Cd L) ~ 100
   -> firmly in the vortex-street shedding regime
Strouhal number                                St = L/(TU) ~ 0.2
   -> with L = 2 km, U = 0.5 m/s, shedding period T ~ 6 h
separation criterion (Garrett 1995)            radius of curvature < H/Cd = 10-100 km
   -> flow separates at essentially every headland here, Peleliu included
```

The 6 h intrinsic shedding timescale is close to the M2 half-cycle, which they
note "may effectively generate eddies" — relevant to the shedding-direction
question in §7 and to `PRESSURE_ANALYSIS.md` §4.

**Implements:** nothing, but it is the physical context for the whole drifter
analysis and should be cited in any write-up of it.

### 9. Rudnick et al. (2019) — vorticity in flow past an island

> Rudnick, D. L., K. L. Zeiden, C. Y. Ou, T. M. S. Johnston, J. A. MacKinnon,
> M. H. Alford, and G. Voet, 2019: Understanding vorticity caused by flow
> passing an island. *Oceanography* **32**(4), 66–73.
> doi:[10.5670/oceanog.2019.412](https://doi.org/10.5670/oceanog.2019.412)
> — open access.

A deliberately pedagogical treatment of exactly the quantity this project
measures: relative vorticity, planetary vorticity, Ro = ζ/f, and why a
solid-body patch has ζ equal to twice its rotation rate. Spans satellite,
glider, ship and mooring observations at Palau, reporting "Ro as large as 30
from eddies of diameter 1 km" and noting vorticity is broad-banded in time,
reaching semidiurnal frequencies under strong westward flow.

Useful here for two reasons: it is the clearest short statement of the
conventions this project uses, and its broad-band point is the same
observation behind our decision not to treat the vortex as a single stationary
oscillation (entry 7).

### 10. St. Laurent et al. (2019) — turbulence and vorticity in the Palau wake

> St. Laurent, L., T. Ijichi, S. T. Merrifield, J. Shapiro, and H. L. Simmons,
> 2019: Turbulence and vorticity in the wake of Palau. *Oceanography* **32**(4),
> 102–109.
> doi:[10.5670/oceanog.2019.416](https://doi.org/10.5670/oceanog.2019.416)
> — open access.

Glider microstructure in the wake. The finding that matters here: direct
wind-driven mixing accounts for only ~10 % of observed turbulence, the rest
coming from shear within the vorticity field of the wake. Context for what the
vortex we measured is doing energetically, and a reminder that the wake is
where dissipation lives.

### 11. Siegelman et al. (2019) — near-inertial surface currents at Palau

> Siegelman, M., M. A. Merrifield, E. Firing, J. A. MacKinnon, M. H. Alford,
> G. Voet, H. W. Wijesekera, T. A. Schramek, K. L. Zeiden, and E. J. Terrill,
> 2019: Observations of near-inertial surface currents at Palau.
> *Oceanography* **32**(4), 74–83.
> doi:[10.5670/oceanog.2019.413](https://doi.org/10.5670/oceanog.2019.413)
> — open access. Pages and DOI read off the article's own citation block.

Held specifically to check the claim in entry 7 that near-inertial analysis is
inapplicable at this latitude. It is the local authority on what near-inertial
motion at Palau looks like and on what timescale, against the 98.3 h inertial
period at 6.99 °N. Consult it before any statement about near-inertial energy
here; do not infer our record can resolve it.

### 12. Johnston et al. (2019) — the FLEAT programme overview

> Johnston, T. M. S., M. C. Schönau, T. Paluszkiewicz, J. A. MacKinnon,
> B. K. Arbic, P. L. Colin, M. H. Alford, M. Andres, L. Centurioni,
> H. C. Graber, K. R. Helfrich, V. Hormann, P. F. J. Lermusiaux,
> R. C. Musgrave, B. S. Powell, B. Qiu, D. L. Rudnick, H. L. Simmons,
> L. St. Laurent, E. J. Terrill, D. S. Trossman, G. Voet, H. W. Wijesekera,
> and K. L. Zeiden, 2019: Flow Encountering Abrupt Topography (FLEAT): A
> multiscale observational and modeling program to understand how topography
> affects flows in the western North Pacific. *Oceanography* **32**(4), 10–21.
> doi:[10.5670/oceanog.2019.407](https://doi.org/10.5670/oceanog.2019.407)
> — open access.

The programme overview. ARCTERX follows FLEAT at the same island, so this is
the map of what was already measured and by whom — read it to find the right
prior dataset rather than for a specific number.

### 13. Essink et al. (2022) — how wrong is a drifter-cluster gradient?

> Essink, S., V. Hormann, L. R. Centurioni, and A. Mahadevan, 2022: On
> characterizing ocean kinematics from surface drifters. *Journal of
> Atmospheric and Oceanic Technology* **39**(8), 1183–1198.
> doi:[10.1175/JTECH-D-21-0068.1](https://doi.org/10.1175/JTECH-D-21-0068.1)
> — obtained from the NOAA Institutional Repository.

The modern successor to entries 1–2, and the paper this project should have had
from the start. It evaluates three cluster methods — Saucier, Kawai, and the
least-squares fit of Molinari & Kirwan (1975), our entry 2 — against synthetic
drifters in a numerical model where truth is known, at exactly our
O(1–10) km scales.

It decomposes the error into (i) the method itself, (ii) **aliasing of
unresolved scales**, and (iii) GPS error, and derives ideal cluster parameters
in the number of drifters *N*, length scale, and aspect ratio. Two things land
squarely on this project:

- **Aliasing of unresolved scales is a named, quantified error term.** That is
  our "the flow is curved across the cluster" — the reason the leave-one-out
  spread is 17 % where the formal error is 5 % (`DRIFTER_ANALYSIS.md` §5). This
  gives it a literature name and an independent magnitude.
- **Aspect ratio drives the error**, and GPS error becomes comparable to the
  other two terms when the aspect ratio is small. That is the justification for
  our isoperimetric gate, currently set by hand.

**Action, not yet taken:** compare our gate threshold (4πA/P²) and our N = 4
against their recommended parameters. If they disagree, theirs is the
better-grounded number.

### 14. Kloosterziel & van Heijst (1991) — the modified Rayleigh criterion

> Kloosterziel, R. C., and G. J. F. van Heijst, 1991: An experimental study of
> unstable barotropic vortices in a rotating fluid. *Journal of Fluid
> Mechanics* **223**, 1–24. — received 13 April 1990.
> Copy posted by R. C. Kloosterziel (Univ. of Hawai‘i).

Source for the stability criterion used in `docs/02-vortex-structure.md` §2.6,
which previously appeared there without one. Laboratory barotropic vortices in
a rotating tank; §4 derives instability criteria for cyclones and anticyclones
"from a modified version of Rayleigh's (circulation) theorem" — the generalised
Rayleigh discriminant

```
Phi = (f + zeta)(f + 2 V/r) > 0   for stability
```

which is what shows our Ro = −67 anticyclone is inertially stable
(Φ = +1.09 × 10⁻⁶) despite f + ζ < 0.

**Read the caveat with it.** Their central result is that cyclones and
anticyclones behave *differently*, and that the anticyclone is the fragile
case — "rather explosive instability behaviour", splitting into dipoles. That
is **barotropic (shear) instability, a different criterion from the inertial
one**. Φ > 0 rules out inertial/centrifugal instability only. Nothing in this
project has tested our vortex against a shear-instability criterion, and its
25.4 h coherence is an observation, not a prediction.

**Implements:** nothing; supports the §2.6 stability argument.

---

## Non-journal sources

Not papers, but the analysis depends on them, so they belong in the index.

| Source | Used for | Link |
|--------|----------|------|
| CF Conventions v1.13 (Dec 2025) | netCDF metadata written by `mwb_nc.py` | <https://cfconventions.org/Data/cf-conventions/cf-conventions-1.13/cf-conventions.html> |
| CF Standard Name Table v94 (2026-06-09) | `standard_name` validation; confirmed no CF name exists for GPS DOP / satellite count / accuracy estimates | <https://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml> |
| ACDD 1.3 | discovery metadata in `mwb_nc.py` | <https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3> |
| u-blox UBX-NAV-PVT message definition | field order/offsets of the 85-byte `.dat` record | u-blox receiver interface description — **generation not identified**; see note |
| CORDC / Scripps | instrument and delivered netCDF provenance | <http://cordc.ucsd.edu/> |

**Note on the u-blox reference.** The 85-byte raw record was decoded
*empirically* and then confirmed bit-exact against the CORDC-delivered netCDF —
it does not rest on the u-blox document. The document is corroborating context
for the fields the netCDF does not carry (`hAcc`, `vAcc`, `sAcc`, `headAcc`,
`tAcc`, `nano`, `iTOW`, `fixType`, `flags`). The exact receiver generation is not
recoverable from the data, so pin down the module before citing a specific
interface-description revision. See the provenance section of `mwb_dat.py`,
which separates what is verified from what is inferred.

---

## Files here

```
README.md                                            this index
references.bib                                       BibTeX, ready for a manuscript
fetch_papers.sh                                      re-attempt any missing download
OkuboEbbesmeyer1976_DSR_drogue_vorticity.pdf         Deep-Sea Res 23, 349-352
MolinariKirwan1975_JPO_differential_kinematics.pdf   J Phys Oceanogr 5, 483-491
Okubo1970_DSR_velocity_singularities.pdf             Deep-Sea Res 17, 445-454
Weiss1991_PhysicaD_enstrophy_transfer.pdf            Physica D 48, 273-294
EfronGong1983_AmStat_bootstrap_jackknife.pdf         Am Statistician 37(1), 36-48
Poulain2023_FrontMarSci_cyprus_gyre_wavelet_ridge.pdf  Front Mar Sci 10, 1266040
Johnston2019_Oceanography_peleliu_wake_eddies_lee_waves.pdf     Oceanography 32(4), 110-125
Rudnick2019_Oceanography_vorticity_flow_past_island.pdf         Oceanography 32(4), 66-73
StLaurent2019_Oceanography_palau_wake_turbulence_vorticity.pdf  Oceanography 32(4), 102-109
Siegelman2019_Oceanography_palau_near_inertial_surface.pdf      Oceanography 32(4), 74-83
Johnston2019_Oceanography_FLEAT_program_overview.pdf            Oceanography 32(4), 10-21
Essink2022_JTECH_drifter_cluster_kinematics.pdf                 J Atmos Ocean Tech 39(8), 1183-1198
KloosterzielVanHeijst1991_JFM_unstable_barotropic_vortices.pdf  J Fluid Mech 223, 1-24
```

Naming convention: `Author####_Journal_short-title.pdf`.

## Provenance

- Entries **1–5**: volume, issue and page numbers were read directly off each
  article's title page in this directory. Efron & Gong was previously flagged
  as unverified and is now confirmed (*Am. Statistician* **37**(1), Feb 1983,
  pp. 36–48).
- Entry **6** (Saffman) is **not yet verified** — the copy is outstanding.
- Entry **7** (Poulain et al.): volume, article number, DOI and the
  received/accepted/published dates were read off the article's own first
  page. Author names are recorded as the initials printed there.
- Entries **8–12** (*Oceanography* 32(4), the FLEAT special issue): volume,
  issue, pages and DOI read off each article's own first page. Siegelman's were
  taken from the article's self-citation block and cross-checked against the
  printed page footers, after an initial guess of 92–101 / `.415` proved wrong
  — the correct values are **74–83** / `.413`.
- Entry **13** (Essink) is the NOAA Institutional Repository copy.
- Entry **14** (Kloosterziel & van Heijst): volume and pages read off the title
  page. The author-posted copy carries no DOI, so the DOI alone was confirmed
  against Crossref rather than the article.
- The **Wanted** table is Crossref metadata only and is **unverified** — no
  title page has been seen for any of it.
- The Weiss PDF is a copy posted publicly by a UCSD course; cite the Physica D
  version of record.
- No paywall was circumvented in assembling this directory. The publisher sites
  for the Wanted entries returned HTTP 403 to an ordinary `curl` and were left
  alone; no user agent was spoofed and no proxy was used.

## Wanted — behind a paywall

Identified as relevant and **not obtained**. All five publisher sites returned
HTTP 403 to an ordinary `curl`, which — as `fetch_papers.sh` says — usually
means the request was blocked as automated, not that access is lacking. Opening
these in a browser on a machine with institutional access should work. Metadata
below is from **Crossref, not from a title page**, so treat volume/issue/pages
as unverified until the PDF is in hand.

| Priority | Citation | DOI |
|---|---|---|
| **1** | Zeiden, K. L., D. L. Rudnick, J. A. MacKinnon, V. Hormann, and L. Centurioni, 2022: Vorticity in the wake of Palau from Lagrangian surface drifters. *J. Phys. Oceanogr.* **52**(9), 2237–2255. | [10.1175/JPO-D-21-0252.1](https://doi.org/10.1175/JPO-D-21-0252.1) |
| **2** | Huntley, H. S., M. Berta, G. Esposito, A. Griffa, B. Mourre, and L. Centurioni, 2022: Conditions for reliable divergence estimates from drifter triplets. *J. Atmos. Oceanic Technol.* **39**(10), 1499–1523. | [10.1175/JTECH-D-21-0161.1](https://doi.org/10.1175/JTECH-D-21-0161.1) |
| 3 | Ohlmann, J. C., M. J. Molemaker, B. Baschek, B. Holt, G. Marmorino, and G. Smith, 2017: Drifter observations of submesoscale flow kinematics in the coastal ocean. *Geophys. Res. Lett.* **44**(1), 330–337. | [10.1002/2016GL071537](https://doi.org/10.1002/2016GL071537) |
| 4 | MacKinnon, J. A., M. H. Alford, G. Voet, K. L. Zeiden, T. M. S. Johnston, M. Siegelman, S. Merrifield, and others, 2019: Eddy wake generation from broadband currents near Palau. *J. Geophys. Res. Oceans* **124**(7), 4891–4903. | [10.1029/2019JC014945](https://doi.org/10.1029/2019JC014945) |
| 5 | Zeiden, K. L., D. L. Rudnick, and J. A. MacKinnon, 2019: Glider observations of a mesoscale oceanic island wake. *J. Phys. Oceanogr.* **49**(9), 2217–2235. | [10.1175/JPO-D-18-0233.1](https://doi.org/10.1175/JPO-D-18-0233.1) |
| 6 | Siegelman, M. N., E. Firing, M. A. Merrifield, J. M. Becker, and R. C. Musgrave, 2023: Near-inertial surface currents around islands. *J. Phys. Oceanogr.* **53**(2), 433–455. | [10.1175/JPO-D-21-0310.1](https://doi.org/10.1175/JPO-D-21-0310.1) |

**Why #1 matters most.** Zeiden et al. (2022) is Palau, SVP drifter *clusters*,
relative vorticity — and it "compares estimates of vorticity from both velocity
spatial gradients (least squares fitting) and velocity time series (wavelet
analysis)". That is the same estimator comparison this project ran, on the same
island, and it is the published answer to the question entry 7 was evaluated
against. Two things to check on arrival:

- Which estimator they preferred and why, against our conclusion that the
  circulation/Stokes contour estimate is primary.
- Their clusters were **~5 km** and reported vorticity up to **6f** in
  submesoscale eddies. Ours is 67f on a ~0.64 km cluster. If the area-average
  argument in `docs/02` §2.3 is right, a 5 km cluster around a ~1.2 km core
  *should* report a far smaller |ζ| — so their 6f and our 67f may be the same
  vortex population seen at different footprints. **That is a quantitative,
  falsifiable prediction and it is worth testing explicitly**; if the numbers do
  not reconcile, the area-average explanation is incomplete.

Also worth chasing but not yet located: Wolanski, E., J. Imberger, and
M. L. Heron, 1984: Island wakes in shallow coastal waters, *J. Geophys. Res.*
**89**(C6), 10553–10569 ([10.1029/JC089iC06p10553](https://doi.org/10.1029/JC089iC06p10553))
— the source of the island wake parameter quoted in entry 8; and Dong, C.,
J. C. McWilliams, and A. F. Shchepetkin, 2007: Island wakes in deep water,
*J. Phys. Oceanogr.* **37**(4).

## A note on redistribution

`.gitignore` excludes `*.pdf` and `LICENSE-DATA` states that journal articles
here are not redistributed. That blanket claim is now **stricter than the
licences require**, and the gap has grown:

| Entries | Licence | Redistributable? |
|---|---|---|
| 1–5 | publisher copyright | **No** |
| 7 (Poulain) | CC BY 4.0, Frontiers | Yes, with attribution |
| 8–12 (*Oceanography*) | **CC BY 4.0** — "This is an open access article made available under the terms of the Creative Commons Attribution 4.0 International License" | Yes, with attribution |
| 13 (Essink) | AMS; obtained via the NOAA Institutional Repository | Check before redistributing |
| 14 (Kloosterziel) | author-posted copy, JFM/CUP copyright | **No** |

So **six** of the local PDFs could legally be committed, not one. They are all
still excluded, because a blanket rule is easier to keep correct than a
per-file exception, and because the repository's value is the analysis, not a
PDF mirror.

If you do want the CC BY ones committed, change `.gitignore` and the
`papers/*.pdf` clause in `LICENSE-DATA` **together** — that clause is an
outward-facing licence claim and must not go stale. Entries 1–5 and 14 must
stay excluded either way.
