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

All six journal articles are held locally and every volume/issue/page figure
below was read off the article's own title page — not from a publisher web
record. Saffman (1992) is a book and is still outstanding; its details remain
unverified until the copy arrives.

Entries **1–5** are implemented. Entry **7** is consulted but deliberately
**not** implemented — see the entry for the arithmetic behind that decision.

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
- The Weiss PDF is a copy posted publicly by a UCSD course; cite the Physica D
  version of record.
- No paywall was circumvented in assembling this directory.

## A note on redistribution

`.gitignore` excludes `*.pdf` and `LICENSE-DATA` states that journal articles
here are not redistributed. **Entry 7 is the one exception in principle** —
Poulain et al. (2023) is Frontiers open access under CC BY, so it *may* be
redistributed with attribution. It is still excluded, because the blanket rule
is simpler to keep correct than a per-file exception. If you want it committed,
change `.gitignore` and the `papers/*.pdf` clause in `LICENSE-DATA` together —
that clause is an outward-facing licence claim and must not go stale.
