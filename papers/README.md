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

All five journal articles are held locally and every volume/issue/page figure
below was read off the article's own title page — not from a publisher web
record. Saffman (1992) is a book and is still outstanding; its details remain
unverified until the copy arrives.

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
```

Naming convention: `Author####_Journal_short-title.pdf`.

## Provenance

- Entries **1–5**: volume, issue and page numbers were read directly off each
  article's title page in this directory. Efron & Gong was previously flagged
  as unverified and is now confirmed (*Am. Statistician* **37**(1), Feb 1983,
  pp. 36–48).
- Entry **6** (Saffman) is **not yet verified** — the copy is outstanding.
- The Weiss PDF is a copy posted publicly by a UCSD course; cite the Physica D
  version of record.
- No paywall was circumvented in assembling this directory.
