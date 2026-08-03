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
| 15 | Zeiden et al. 2022 ★★ | ✓ `Zeiden2022_JPO_palau_wake_vorticity_drifters.pdf` | title page |
| 16 | Huntley et al. 2022 ★ | ✓ `Huntley2022_JTECH_drifter_triplet_divergence.pdf` | title page |
| 17 | MacKinnon et al. 2019 ★ | ✓ `MacKinnon2019_JGR_palau_eddy_wake_broadband.pdf` | title page |
| 18 | Ohlmann et al. 2017 | ✓ `Ohlmann2017_GRL_submesoscale_drifter_kinematics.pdf` | title page |
| 19 | Zeiden et al. 2019 | ✓ `Zeiden2019_JPO_glider_island_wake_palau.pdf` | title page |
| 20 | Siegelman et al. 2023 | ✓ `Siegelman2023_JPO_near_inertial_around_islands.pdf` | title page |
| 21 | Lilly & Gascard 2006 | ✓ `LillyGascard2006_NPG_wavelet_ridge_elliptical.pdf` | title page |
| 22 | Lilly & Pérez-Brunius 2021 | ✓ `LillyPerezBrunius2021_NPG_wavelet_ridge_eddy_detection.pdf` | title page |
| 23 | Spydell et al. 2019 ★ | ✓ `Spydell2019_JTECH_drifter_gps_error_vorticity.pdf` | title page |
| 24 | Rudnick et al. 2015 | ✓ `Rudnick2015_JPO_gulf_mexico_cyclonic_eddies_gliders.pdf` | title page |
| 25 | Pattiaratchi et al. 1987 | ✓ `Pattiaratchi1987_JGR_island_wakes_headland_eddies.pdf` | Crossref (see note) |
| 26 | Dong et al. 2007 | ✓ `Dong2007_JPO_island_wakes_deep_water.pdf` | title page |
| 27 | Shcherbina et al. 2013 | ✓ `Shcherbina2013_GRL_submesoscale_statistics.pdf` | title page |
| 28 | LaCasce 2008 | ✓ `LaCasce2008_ProgOcean_lagrangian_statistics.pdf` | title page |

Every volume/issue/page figure below was read off the article's own title page
— not from a publisher web record. Saffman (1992) is a book and is still
outstanding; its details remain unverified until the copy arrives.

How to read the numbering:

| Entries | What they are |
|---|---|
| **1–5** | **Implemented** — the methods the code actually runs |
| 6 | Lamb–Oseen reference, still outstanding |
| 7 | Consulted, method **declined**; see the entry for the arithmetic |
| **8–12, 17, 19** | **Site and phenomenon** — Palau, and in entry 8's case the Peleliu tip itself |
| 25, 26 | Island-wake theory and modelling: headland vs island, shallow vs deep |
| 27, 28 | Submesoscale reference statistics; Lagrangian statistics (unused) |
| **13, 15, 16, 18, 23, 24** | **Method** — drifter-cluster estimator accuracy and error |
| 14, 20 | Dynamics support: vortex stability, near-inertial theory |
| 21–22 | The wavelet-ridge primary sources behind the entry 7 decision |

The ones marked ★ are to read first if you read nothing else. **23**
(Spydell et al. 2019) argues our cluster-shape gate uses the wrong variable and
supplies the right one, and is the strongest independent support for quoting the
jackknife rather than the formal error.
**15** (Zeiden et al. 2022 ★★) is the closest published analysis to this one and
ran the same estimator comparison; **16** (Huntley et al. 2022) calibrates the
cluster-shape gate we set by hand; **17** (MacKinnon et al. 2019) supplies a
mechanism for the circulation shortfall in `DRIFTER_ANALYSIS.md` §7.

**Everything previously listed as missing has now been obtained.** The next
round of targets is in [Wanted](#wanted--next-round) at the end — headed by
Spydell et al. (2019) and Rudnick et al. (2015), both of which would change
code rather than prose.

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

### 15. Zeiden et al. (2022) — Palau wake vorticity from drifter clusters ★★

> Zeiden, K. L., D. L. Rudnick, J. A. MacKinnon, V. Hormann, and L. Centurioni,
> 2022: Vorticity in the wake of Palau from Lagrangian surface drifters.
> *Journal of Physical Oceanography* **52**(9), 2237–2255.
> doi:[10.1175/JPO-D-21-0252.1](https://doi.org/10.1175/JPO-D-21-0252.1)
> — received 20 Nov 2021, final form 9 May 2022.

**The closest published analysis to this one.** Same island, SVP drifter
clusters, relative vorticity, Okubo–Weiss, and the same two founding references
(Okubo & Ebbesmeyer 1976; Molinari & Kirwan 1975 — our entries 1 and 2). 19
clusters of five drifters, ~5 km scale, released over 2 years; 15 entrained in
the wake.

**It answers the estimator question of entry 7 empirically.** They ran the same
comparison and found a **factor of 2** between methods: for cluster C7 the
wavelet estimate from a *single* drifter's velocity time series gives Ro ~ 6,
while the least-squares fit over *all* drifters gives Ro ~ 3.

That is the same effect documented in `docs/02-vortex-structure.md` §2.3, and
the signs are consistent once you ask which estimator averages over the larger
area:

| | larger-area estimator | reads |
|---|---|---|
| Zeiden C7 | 5 km cluster ≫ eddy core → LS fit | **low** (Ro 3 vs 6) |
| here | orbital loop > 0.64 km cluster → constellation rotation | **low** (9.4 vs 11.9 × 10⁻⁴) |

Different geometry, same rule: **whichever estimator integrates over more area
reads lower, because |ζ| falls with radius.** Our 21 % and their factor of 2 are
one phenomenon. This is the strongest external support the area-average
argument has.

**Three methods here are better than ours and should be considered:**

1. **A solid-body-constrained least-squares fit** (their Eq. 2, after Rudnick
   et al. 2015) alongside the plane fit. It solves for `[u0, v0, ζ/2]` in one
   system, and returns **the eddy center as a fit parameter**. Ours instead
   inverts `A⁻¹(c − U₀)`, which `docs/02` §2.4 calls the weakest link and which
   needs two gates. Theirs is the more robust construction.
2. **The solid-body fit is immune to cluster shape.** They give the error as
   σ_ζ² ≈ (σ_U²/N)(1/σ_L²) — **no aspect-ratio term**, "because in this case the
   velocity gradient does not depend on direction", against
   σ_ζ² ≈ (σ_U²/N)(1 + 1/a²)(1/b²) for the plane fit (Spydell et al. 2019;
   confirmed by Essink et al. 2022, our entry 13). That directly addresses the
   near-collinear blow-up that our isoperimetric gate exists to suppress: they
   use an estimator that does not care, where we discard data.
3. **The solid-body constraint acts as a tide/inertial filter.** For C7 after
   15 days the semidiurnal band held 48 % of the plane-fit vorticity variance
   but only 18 % of the solid-body-fit variance.

**Corroborations of things we found independently:**

- **Outward spiral.** They observe "drifters are moving outwards from within
  their eddy cores", with U_r ~ 0.03–0.05 m s⁻¹, "an order of magnitude weaker
  than the azimuthal velocity". We report radial velocity at 4 % of azimuthal
  and three of four drifters roughly doubling their orbital radius. Same sign,
  same order.
- **Not solid-body.** They fit both a Rankine vortex and a "line vortex ...
  initially irrotational point vortex allowed to decay in time due to friction"
  — that second one is Lamb–Oseen, our `fit_oseen()` — and conclude the Rankine
  form "overpredicts the velocity and vorticity in the vicinity of the eddy
  radius", with "vorticity decreas[ing] with increasing radial distance even
  within the core".
- **Okubo–Weiss as the coherence gate**, defined identically and used for the
  same purpose.

**Two differences that matter, and are not contradictions:**

- **Their wake eddies are cyclonic; ours is anticyclonic.** Theirs were released
  at the **north** end of the archipelago, where westward flow generates
  positive shear; we are at the **south** tip of Peleliu. Johnston et al.
  (entry 8) explicitly expect "a sequence of cyclonic/anticyclonic eddies ...
  formed at the north/south" points. Consistent, not conflicting — but it means
  **their rotary trick does not transfer to us.** They separate wake vorticity
  from internal waves by sense of rotation, since NH internal waves are
  anticyclonic; our vortex is anticyclonic too, so it sits in the same rotary
  half-plane as the wave field. We separate by frequency instead, which is easy
  here (2.9–3.7 h against a 98.3 h inertial period).
- **Scale.** Their eddies are ~40 km diameter with Ro ≳ 1 decaying as t⁻¹; ours
  is ~2.4 km with Ro = −67. The bridge is Johnston et al. at our own headland
  (1–2 km, Ro 65–80), not this paper.

**One quantitative comparison to make, carefully.** They report an enstrophy
spectrum ζ² ∝ k^1.9 ≈ k², i.e. |ζ| ∝ 1/L. Our profile gives |ζ| ∝ L^−0.48
(1.58 × 10⁻³ at ~350 m to 0.98 × 10⁻³ at ~950 m), i.e. ζ² ∝ k^0.96. These are
**not commensurable as they stand** — theirs is an ensemble over many clusters,
times and eddies across 5–50 km in a turbulent wake; ours is one vortex sampled
at different cluster scales inside a single core. Do not report this as a
discrepancy without matching the footprints first.

### 16. Huntley et al. (2022) — shape thresholds for cluster estimates ★

> Huntley, H. S., M. Berta, G. Esposito, A. Griffa, B. Mourre, and
> L. Centurioni, 2022: Conditions for reliable divergence estimates from drifter
> triplets. *Journal of Atmospheric and Oceanic Technology* **39**(10),
> 1499–1523.
> doi:[10.1175/JTECH-D-21-0161.1](https://doi.org/10.1175/JTECH-D-21-0161.1)

Derives, against a high-resolution ROMS simulation where truth is known, the
threshold at which a drifter triangle is too degenerate to trust. Two metrics
are found "equally effective", especially **at scales of 5 km and below** —
which is our regime:

```
scaled aspect ratio  Lambda = 0.20     (Lambda = 1 equilateral, 0 collinear)
largest interior angle  theta = 0.86 pi
```

**This is directly actionable and our gate is looser than theirs.** We gate on
the isoperimetric quotient 4πA/P² with `min_quality=0.10`
(`src/eddy_kinematics.py:336`). Normalising both to the equilateral triangle,
where 4πA/P² = π√3/9 = 0.6046:

| | our gate | Huntley |
|---|---|---|
| as Λ (equilateral = 1) | **0.165** | **0.20** |
| as 4πA/P² | 0.10 | 0.121 |

So our hand-set threshold is ~20 % more permissive than the model-calibrated
recommendation — and they note their own criteria are already "less stringent
than some of the ad hoc criteria previously used".

A second issue this exposes: we apply **one absolute threshold to both the
four-drifter polygon and the leave-one-out triangles**, but the quotient
maxes at 0.6046 for a triangle and 0.7854 for a square. The same 0.10 is
therefore 16.5 % of maximum for triangles and 12.7 % for quadrilaterals — the
gate is relatively *more* permissive exactly where we have fewer constraints.

**Action, not taken:** raising `min_quality` changes published numbers, so it
is left alone pending a decision. The §5 spike to ζ = 0.27 s⁻¹ that motivated
the gate suggests the current value works in practice; the point is that it is
uncalibrated, and now need not be.

They also warn that shape-based discarding "necessarily biases the distribution
of divergence estimates slightly toward positive values" — relevant to our
"divergence < 2 % of |ζ|" claim, though small and opposite to the Lagrangian
convergence-sampling bias.

### 17. MacKinnon et al. (2019) — tidal boost to wake vorticity flux ★

> MacKinnon, J. A., M. H. Alford, G. Voet, K. L. Zeiden, T. M. S. Johnston,
> M. Siegelman, S. Merrifield, and M. Merrifield, 2019: Eddy wake generation
> from broadband currents near Palau. *Journal of Geophysical Research: Oceans*
> **124**(7), 4891–4903.
> doi:[10.1029/2019JC014945](https://doi.org/10.1029/2019JC014945)

Measurements at the flow-separation point at the **north** end of Palau, where
"energetic tides and vertically sheared low-frequency flows are both present" —
the same broadband forcing regime as Peleliu. Small-scale (~1 km) wake eddies of
both vorticity signs on either side of the separation point, evolving over
several tidal periods.

**Its central result speaks straight to our circulation budget.** Most wake
work treats either steady flow or purely tidal flow; this paper treats both
together, and finds that including high-frequency oscillatory currents "may
boost the net flux of vorticity into the ocean interior by a depth dependent
factor of **2 to 25**", concluding that models omitting them "may not
accurately infer the net momentum or energy losses".

`DRIFTER_ANALYSIS.md` §7 records that Γ needs **2–3×** more forcing than one M2
half-cycle at the array-averaged 0.64 m s⁻¹ supplies. A tidal-rectification
boost of 2–25× is a *mechanism* for exactly that shortfall, and it is the
second of the two horns in §7 ("circulation accumulates over several cycles")
rather than the first. Together with entry 8's local flow intensification,
**both horns now have supporting evidence and they are not exclusive.** Neither
is confirmed for our event.

### 18. Ohlmann et al. (2017) — submesoscale kinematics from drifter clusters

> Ohlmann, J. C., M. J. Molemaker, B. Baschek, B. Holt, G. Marmorino, and
> G. Smith, 2017: Drifter observations of submesoscale flow kinematics in the
> coastal ocean. *Geophysical Research Letters* **44**(1), 330–337.
> doi:[10.1002/2016GL071537](https://doi.org/10.1002/2016GL071537)

Nine drifters in a 3 × 3 grid at 1 km spacing, deployed onto features spotted
from the air, with kinematics computed from "all possible four-drifter
clusters" — the same leave-one-out-style construction we use. Reports mean
divergence and vorticity that "can exceed 5f", noted as the largest observed in
the field at the time, and an explicit departure from geostrophy.

Two things to take:

- **Precedent for our regime.** It establishes that O(1 km) drifter clusters
  measuring |Ro| ≫ 1 is a real, published result, not an artifact — useful when
  defending Ro = −67, though our value is an order of magnitude beyond theirs.
- **Their aspect-ratio definition**, α = L_minor/L_major over the cluster, is
  a third shape metric alongside ours and Huntley's. If the gate is ever
  recalibrated, pick one of the three and state which.

### 19. Zeiden et al. (2019) — the mesoscale wake at Palau, from gliders

> Zeiden, K. L., D. L. Rudnick, and J. A. MacKinnon, 2019: Glider observations
> of a mesoscale oceanic island wake. *Journal of Physical Oceanography*
> **49**(9), 2217–2235.
> doi:[10.1175/JPO-D-18-0233.1](https://doi.org/10.1175/JPO-D-18-0233.1)
> — received 6 Nov 2018, final form 31 May 2019.

Two years of glider velocity profiles to 1000 m, east and west of Palau. The
incident NEC accelerates around the island from 0.1 to 0.2 m s⁻¹ at the
surface; the lee shows elevated variability and return flow indicating boundary
layer separation. Mean wake vorticity reaches **0.3f** near the surface, with
instantaneous values exceeding f during sustained strong westward flow, so
"ageostrophic effects become important to first order".

This is the **background flow** our vortex was shed into — the island-scale
context, an order of magnitude larger and weaker than the tip vortex. Also
reports that eastward flow produces an asymmetric wake, which is the
independent-variable side of the shedding-direction question in
`DRIFTER_ANALYSIS.md` §7 and `PRESSURE_ANALYSIS.md` §4.

### 20. Siegelman et al. (2023) — near-inertial currents around islands

> Siegelman, M. N., E. Firing, M. A. Merrifield, J. M. Becker, and
> R. C. Musgrave, 2023: Near-inertial surface currents around islands.
> *Journal of Physical Oceanography* **53**(2), 433–455.
> doi:[10.1175/JPO-D-21-0310.1](https://doi.org/10.1175/JPO-D-21-0310.1)
> — received 17 Dec 2021, final form 2 Sep 2022.

Theory for how an island modifies wind-generated NIOs, motivated by the Palau
observations in entry 11. Small islands enhance near-inertial currents via
island-trapped waves; large islands suppress them through interference between
incident and reflected Poincaré waves.

**Held as the check on a negative claim.** Entry 7 and `docs/02` §2.6 assert
that near-inertial analysis is inapplicable to this dataset. That rests on
record length (98.3 h inertial period against 25.4 h) and on Ro = −67 breaking
the linearised f_e — **not** on any claim that NIOs are absent at Palau. This
paper shows they are present and island-modified. The two statements are
compatible; keep them distinct, and cite this before saying anything about
near-inertial energy here.

### 21–22. Lilly & Gascard (2006); Lilly & Pérez-Brunius (2021) — wavelet ridges

> Lilly, J. M., and J.-C. Gascard, 2006: Wavelet ridge diagnosis of
> time-varying elliptical signals with application to an oceanic eddy.
> *Nonlinear Processes in Geophysics* **13**, 467–483.
> doi:[10.5194/npg-13-467-2006](https://doi.org/10.5194/npg-13-467-2006)
>
> Lilly, J. M., and P. Pérez-Brunius, 2021: Extracting statistically significant
> eddy signals from large Lagrangian datasets using wavelet ridge analysis, with
> application to the Gulf of Mexico. *Nonlinear Processes in Geophysics*
> **28**, 181–212.
> doi:[10.5194/npg-28-181-2021](https://doi.org/10.5194/npg-28-181-2021)
> — both open access (CC BY).

The method behind entry 7, and behind the wavelet half of entry 15, held here
so the decision not to adopt it rests on primary sources rather than on someone
else's description of them. The 2021 paper is the more useful: it is about
extracting *statistically significant* eddy signals from large Lagrangian
datasets, and its significance testing is the part worth reading if the
single-drifter approach is ever revisited.

Neither changes the conclusion in entry 7. The binding constraint here is 6.9
orbital cycles of record, not the quality of the method.

### 23. Spydell et al. (2019) — what actually sets the vorticity error ★

> Spydell, M. S., F. Feddersen, and J. MacMahan, 2019: The effect of drifter GPS
> errors on estimates of submesoscale vorticity. *Journal of Atmospheric and
> Oceanic Technology* **36**(11), 2101–2119.
> doi:[10.1175/JTECH-D-19-0108.1](https://doi.org/10.1175/JTECH-D-19-0108.1)
> — received 25 Jun 2019, final form 22 Aug 2019.

Derives the **a priori** vorticity error for a drifter cluster, and validates it
against two stationary GPS experiments where the true vorticity is exactly zero
— a null test of the same kind this project relies on. Their Eq. (16):

```
sigma_zeta^2 = (1/N) (sigma_u^2 / la^2) (1 + la^2/lb^2) (1 - rho_u1u2)
```

with `la <= lb` the minor and major axes of the cluster (from the position
covariance matrix), `N` the drifter count, `sigma_u` the velocity error and
`rho` the cross-drifter error correlation. Error falls only as `N^-1/2`, so
"large numbers of drifters are required to reduce vorticity error
substantially" — with N = 4 we are near the floor of what helps.

**Their headline conclusion contradicts how we gate.** In their words:
"Previously, cluster area or ellipticity were used as criteria to distinguish
error. We show that the drifter cluster **minor axis** (narrowness) is a key
time-dependent factor affecting vorticity error." They find vorticity error
exceeding 5f once the minor axis drops below 50 m, even at velocity errors under
0.004 m s⁻¹. We gate on the isoperimetric quotient — an area/ellipticity
measure, exactly the class they argue against.

**Applied to our cluster** (`data/eddy_kinematics.npz`, 355 windows):

| | |
|---|---|
| minor axis λ_a, median | **236 m** (p5 = 64 m, min 12 m) |
| major axis λ_b, median | 741 m |
| windows with λ_a < 50 m | 12 (3.4 %) |
| — caught by our existing gate | 9 |
| — **leaked through** | **3** |
| correlation, our quality vs λ_a | **r = 0.82** |

So the isoperimetric gate is largely doing the right thing *by proxy* — it
correlates strongly with narrowness — but it is not the right variable and
three narrow windows survive it. Those three carry median |ζ| = 3.2 × 10⁻³ s⁻¹,
**2.7× the record median**, which is the signature Spydell predicts.

Adding `λ_a >= 50 m` as a second gate:

```
current gate        n=343   median zeta = -1.190e-3   Ro = -67.0
+ minor-axis gate   n=340   median zeta = -1.187e-3   Ro = -66.9   (0.23 % change)
p95 |zeta|          2.42e-3  ->  2.31e-3   (-4.5 %)
```

**The headline is untouched and the tail tightens.** This is a cheap,
principled improvement; it is not applied, because changing gates changes
published numbers and that is a decision, not a detail. The 0.23 % result is
also a useful null test in its own right: it confirms the medians are robust to
the gate, which is what `DRIFTER_ANALYSIS.md` §5 claims.

**It also independently vindicates quoting the jackknife.** Eq. (16) propagates
*instrument* error only — the authors say so explicitly. Evaluated on our
geometry it predicts:

| σ_u | median σ_ζ | as % of \|ζ\| |
|---|---|---|
| 0.010 m s⁻¹ (SVP GPS, per Essink) | 2.3 × 10⁻⁵ | 1.9 % |
| 0.004 m s⁻¹ | 9.0 × 10⁻⁶ | 0.8 % |
| 0.002 m s⁻¹ | 4.5 × 10⁻⁶ | 0.4 % |

against our **5 %** formal error and **17 %** jackknife spread. Instrument error
is roughly a tenth of the jackknife, so the 17 % really is dominated by flow
curvature and unresolved scales — Essink's "aliasing of unresolved scales" —
not by GPS noise. Two caveats, both conservative: ρ was set to 0, and a positive
error correlation would make σ_ζ *smaller* still; and our buoys report Doppler
velocity block-averaged over 2048 samples, so σ_u is likely at the low end.

**Action:** if the error budget is ever revisited, this is the paper to
implement — report λ_a alongside the quotient, and gate on it.

### 24. Rudnick et al. (2015) — the solid-body fit

> Rudnick, D. L., G. Gopalakrishnan, and B. D. Cornuelle, 2015: Cyclonic eddies
> in the Gulf of Mexico: observations by underwater gliders and simulations by
> numerical model. *Journal of Physical Oceanography* **45**(1), 313–326.
> doi:[10.1175/JPO-D-14-0138.1](https://doi.org/10.1175/JPO-D-14-0138.1)
> — received 16 Jul 2014, final form 20 Oct 2014.

The source Zeiden et al. (entry 15) cite for the solid-body-constrained fit.
Their §3 fits "a model of solid body rotation in which the eastward u and
northward v velocity components are given by" a centre plus a rigid rotation,
with "z ... relative vertical vorticity (twice the rotation rate)", and Table 1
tabulates the results. They also test the assumption rather than assume it,
noting the model describes "most, but not all" of the observed structure and
relaxing it to check.

Held for the method, not the region. Combined with entries 15 and 23 it is the
recipe for the estimator this project does not yet have: a fit that returns the
eddy centre directly and whose error does not blow up as the polygon flattens.

### 25. Pattiaratchi et al. (1987) — island wakes vs headland eddies

> Pattiaratchi, C., A. James, and M. Collins, 1987: Island wakes and headland
> eddies: a comparison between remotely sensed data and laboratory experiments.
> *Journal of Geophysical Research: Oceans* **92**(C1), 783–794.
> doi:[10.1029/JC092iC01p00783](https://doi.org/10.1029/JC092iC01p00783)

**Citation caveat.** The PDF's own running header reads "JANUARY 15, **1986**",
which is wrong — volume 92 is 1987, and Crossref gives an issue date of
15 January 1987. This is a typo in the printed journal, not a mis-download.
**Cite 1987.** Recorded here because a header is a claim, not a fact, and this
one would otherwise propagate.

Identifies island wakes and headland eddies in visible-band satellite and
airborne imagery of the Bristol and English Channels — high tidal currents,
high turbidity, suspended sediment as a passive tracer — and compares them
against laboratory experiments. Directly relevant as the **headland** case: our
vortex is shed from a headland (the Peleliu tip), not from an isolated island,
and this is the paper that treats the two side by side.

### 26. Dong et al. (2007) — island wakes in deep water

> Dong, C., J. C. McWilliams, and A. F. Shchepetkin, 2007: Island wakes in deep
> water. *Journal of Physical Oceanography* **37**(4), 962–981.
> doi:[10.1175/JPO3047.1](https://doi.org/10.1175/JPO3047.1)
> — received 21 Oct 2005, final form 12 Sep 2006.

The deep-water counterpart to Wolanski's shallow-water regime, and the source of
the Strouhal/shedding expectations invoked in entry 8. Idealised ROMS
simulations of flow past an island: wake instability, coherent vortex formation,
and the mesoscale-to-submesoscale eddy field that follows. This is the modelling
reference against which the Peleliu observations should be read, and the natural
place to check whether our single 25.4 h event is a typical shed vortex or an
unusual one.

### 27. Shcherbina et al. (2013) — the submesoscale reference distribution

> Shcherbina, A. Y., E. A. D'Asaro, C. M. Lee, J. M. Klymak, M. J. Molemaker,
> and J. C. McWilliams, 2013: Statistics of vertical vorticity, divergence, and
> strain in a developed submesoscale turbulence field. *Geophysical Research
> Letters* **40**(17), 4706–4711.
> doi:[10.1002/grl.50919](https://doi.org/10.1002/grl.50919)
> — author copy, staff.washington.edu.

The first consistent sampling of the **full horizontal velocity gradient
tensor** at O(1 km) in the open ocean, from two vessels running parallel tracks.
It supplies what a single case study cannot: a *population* of ζ, δ and strain
at our scale, against which Ro = −67 can be placed rather than merely asserted.
Also the source of the k⁻² near-surface KE spectrum Zeiden et al. use to argue
that unresolved sub-cluster currents are O(0.01) m s⁻¹ — the number that sets
σ_u in entry 23.

### 28. LaCasce (2008) — Lagrangian statistics review

> LaCasce, J. H., 2008: Statistics from Lagrangian observations. *Progress in
> Oceanography* **77**(1), 1–29.
> doi:[10.1016/j.pocean.2008.02.002](https://doi.org/10.1016/j.pocean.2008.02.002)
> — course copy, pordlabs.ucsd.edu.

The standard review of single- and multi-particle Lagrangian statistics.
**Not currently used**: this project measures kinematics, not dispersion, and
25.4 h of four drifters is far too short for meaningful dispersion statistics.
Held because Zeiden et al. lean on it to interpret cluster-scale growth
(distinguishing turbulent cascade from lateral shear), which is the natural
next question if the drifter record is ever extended.


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
Zeiden2022_JPO_palau_wake_vorticity_drifters.pdf                J Phys Oceanogr 52(9), 2237-2255
Huntley2022_JTECH_drifter_triplet_divergence.pdf                J Atmos Ocean Tech 39(10), 1499-1523
MacKinnon2019_JGR_palau_eddy_wake_broadband.pdf                 J Geophys Res Oceans 124(7), 4891-4903
Ohlmann2017_GRL_submesoscale_drifter_kinematics.pdf             Geophys Res Lett 44(1), 330-337
Zeiden2019_JPO_glider_island_wake_palau.pdf                     J Phys Oceanogr 49(9), 2217-2235
Siegelman2023_JPO_near_inertial_around_islands.pdf              J Phys Oceanogr 53(2), 433-455
LillyGascard2006_NPG_wavelet_ridge_elliptical.pdf               Nonlin Processes Geophys 13, 467-483
LillyPerezBrunius2021_NPG_wavelet_ridge_eddy_detection.pdf      Nonlin Processes Geophys 28, 181-212
Spydell2019_JTECH_drifter_gps_error_vorticity.pdf               J Atmos Ocean Tech 36(11), 2101-2119
Rudnick2015_JPO_gulf_mexico_cyclonic_eddies_gliders.pdf         J Phys Oceanogr 45(1), 313-326
Pattiaratchi1987_JGR_island_wakes_headland_eddies.pdf           J Geophys Res 92(C1), 783-794
Dong2007_JPO_island_wakes_deep_water.pdf                        J Phys Oceanogr 37(4), 962-981
Shcherbina2013_GRL_submesoscale_statistics.pdf                  Geophys Res Lett 40(17), 4706-4711
LaCasce2008_ProgOcean_lagrangian_statistics.pdf                 Prog Oceanogr 77(1), 1-29
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
- Entries **15–20** were retrieved by the maintainer under institutional access
  after the publisher sites refused automated requests. Volume, issue, pages and
  the received/accepted dates were read off each article's own title page. The
  files arrived under publisher or ad hoc names and were renamed to the
  convention here; two carried no useful embedded title, so each was identified
  from its title page rather than its filename.
- Entries **21–22** are the Copernicus open-access PDFs.
- Entries **23–26** were retrieved by the maintainer under institutional access;
  **27–28** are author and course copies. Volume, issue, pages and dates were
  read off each title page, except entry 25 — see the next line.
- Entry **25** (Pattiaratchi) is the one case where the article's own header is
  **wrong**: it prints "JANUARY 15, 1986" on a volume-92 paper, which is 1987.
  Date and DOI taken from Crossref; volume, issue and pages agree with the
  printed header. **Cite 1987.**
- Entry **14** (Kloosterziel & van Heijst): volume and pages read off the title
  page. The author-posted copy carries no DOI, so the DOI alone was confirmed
  against Crossref rather than the article.
- The **Wanted** table is Crossref metadata and citing-paper reference lists,
  and is **unverified** — no title page has been seen for any of it.
- The Weiss PDF is a copy posted publicly by a UCSD course; cite the Physica D
  version of record.
- No paywall was circumvented in assembling this directory. The publisher sites
  for the Wanted entries returned HTTP 403 to an ordinary `curl` and were left
  alone; no user agent was spoofed and no proxy was used.

## Wanted — two left

| Citation | DOI / ISBN | Why, and where to look |
|---|---|---|
| Wolanski, E., J. Imberger, and M. L. Heron, 1984: Island wakes in shallow coastal waters. *J. Geophys. Res.* **89**(C6), 10553–10569. | [10.1029/JC089iC06p10553](https://doi.org/10.1029/JC089iC06p10553) | Origin of the island wake parameter Ref = H/(Cd L) quoted in entry 8 and proposed for `DRIFTER_ANALYSIS.md` §7; currently cited at second hand. AGU/Wiley only — Crossref lists no alternative host, and pre-1997 JGR is not in an open archive. Wiley 403s automated requests but serves it in a browser. Wolanski was at AIMS; their repository is worth a look if Wiley fails. |
| Saffman, P. G., 1992: *Vortex Dynamics*. Cambridge University Press. | ISBN 978-0-521-42058-7 | Outstanding since the beginning — entry 6. The Lamb–Oseen citation stays unverified until a copy arrives, and the page count in `references.bib` was quoted from memory. A library copy settles it; nothing online will. |

Everything else identified in this project's literature search has been
obtained. **The four that were hard to find were all findable** — the DOIs were
correct; the barrier was hosting, not identification:

| | Route that worked |
|---|---|
| LaCasce 2008 | course copy, `pordlabs.ucsd.edu` |
| Shcherbina 2013 | author copy, `staff.washington.edu/shcher` |
| Dong 2007 | retrieved by the maintainer under institutional access |
| Pattiaratchi 1987 | retrieved by the maintainer under institutional access |

General lesson for the next one of these: when a publisher 403s, try the
**author's institutional page** and **course reading lists** before assuming the
paper is out of reach. Both hits above came from those two routes, not from the
publisher and not from a repository aggregator.

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
| 15–20 | AMS / AGU, publisher copyright; retrieved by the maintainer under institutional access | **No** |
| 21–22 (*Nonlin. Processes Geophys.*) | **CC BY** (Copernicus) | Yes, with attribution |

So **eight** of the local PDFs could legally be committed, not one. They are all
still excluded, because a blanket rule is easier to keep correct than a
per-file exception, and because the repository's value is the analysis, not a
PDF mirror.

If you do want the CC BY ones committed, change `.gitignore` and the
`papers/*.pdf` clause in `LICENSE-DATA` **together** — that clause is an
outward-facing licence claim and must not go stale. Entries 1–5 and 14 must
stay excluded either way.
