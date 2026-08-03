# 5. Tidal harmonic analysis

Fitting known astronomical frequencies to a sea-level record, and the several
ways that goes wrong without announcing itself.

Code: `wamos_tpw.tides` (the `wamos_tpw` repository), used here through
`src/pressure_array.py`.

---

## 5.1 The model

The tide is a sum of sinusoids at frequencies fixed by astronomy. Only the
amplitude and phase of each are unknown, and both enter linearly if you write
each constituent in quadrature:

$$\eta(t) = a_0 + \sum_{k} \left[ A_k\cos(\omega_k t) + B_k\sin(\omega_k t) \right]$$

Then amplitude and phase come out as

$$H_k = \sqrt{A_k^2 + B_k^2}, \qquad g_k = \arctan\!\left(\frac{-B_k}{A_k}\right)$$

Because $\omega_k$ is known and only $A_k, B_k$ are free, this is **ordinary
linear least squares** — no iteration, no starting guess, closed form. That is
the entire reason harmonic analysis is the standard tool.

## 5.2 The constituents

| name | period (h) | origin |
|---|---|---|
| M2 | 12.4206 | principal lunar semidiurnal |
| S2 | 12.0000 | principal solar semidiurnal |
| N2 | 12.6583 | larger lunar elliptic |
| K2 | 11.9672 | lunisolar semidiurnal |
| K1 | 23.9345 | lunisolar diurnal |
| O1 | 25.8193 | principal lunar diurnal |
| P1 | 24.0659 | principal solar diurnal |
| Q1 | 26.8684 | larger lunar elliptic diurnal |
| M4, MS4 | 6.21, 6.10 | shallow-water overtides |

At Palau the tide is mixed, mainly semidiurnal: M2 dominates at
**0.502 ± 0.004 m** across all twelve gauges.

## 5.3 The Rayleigh criterion

Two constituents at $\omega_1$ and $\omega_2$ can only be separated if the
record is long enough for their relative phase to complete a full cycle:

$$T_{\text{record}} \geq \frac{2\pi}{|\omega_1 - \omega_2|} = \frac{1}{|f_1 - f_2|}$$

This is the tidal-analysis form of the Rayleigh resolution criterion, and it is
unforgiving. Examples that matter here:

| pair | separation | record needed |
|---|---|---|
| M2 / S2 | 1.016 cyc/day | 14.8 days |
| K1 / O1 | 0.082 cyc/day | 13.7 days |
| **K1 / P1** | 0.0027 cyc/day | **182.6 days** |
| **Mf / MSf** | — | **182.6 days** |

A 29-day record **cannot** separate K1 from P1 by the Rayleigh criterion. The
standard remedy is **inference**: fix the P1/K1 amplitude ratio and phase lag
at their equilibrium-tide values, or at values from a nearby long record, and
solve for one instead of two.

`harmonic_fit` takes an `infer=` argument for this, refuses unresolvable pairs
by default (`allow_unresolvable=False`), and reports `rayleigh_days` and
`unresolvable_pairs` so the caller cannot silently ask for the impossible.

**That guard earned its place on its first real run**, catching Mf and MSf
sitting in the author's own default constituent list — they need 183 days and
the records are 29. MSf was dropped.

## 5.4 Conditioning, and the trap that cost the most

Rayleigh is necessary but **not sufficient**. What actually breaks a fit is the
conditioning of the design matrix, and it can be catastrophic for reasons
Rayleigh does not capture.

### The S1 disaster

Adding S1 (period exactly 24.000 h) to a solve that already contains
K1 (23.934 h) and P1 (24.066 h) creates a **three-way near-degeneracy**. The
three columns are nearly linearly dependent over a 29-day record, and the
condition number goes from $2.9\times10^{3}$ to $7.7\times10^{4}$.

What makes this genuinely dangerous:

> **The corrupted fit has the *lower* residual — 3.08 cm against 4.06 cm.**

Adding a nearly-degenerate column always reduces the residual, because it adds
a degree of freedom. The fit "improves" while the individual constituent
amplitudes become meaningless, trading enormous cancelling values between three
columns that are nearly the same vector. Any model-selection rule based on
residual alone selects the broken model.

This is why `TidalFit` carries `condition` as a first-class field and warns
above `CONDITION_WARN = 1.0e4`. **Watch the conditioning, not the residual.**

### Getting the diagnosis right took three tries

Worth recording because the wrong explanations were plausible:

1. First claim: "K1/P1 are unresolvable in 29 days." *False* — with white noise
   down to 8 cm they separate perfectly well; Rayleigh is conservative.
2. Second claim: "red noise defeats them." *Also false* — tested and rejected.
3. Actual cause: the three-way degeneracy created by **S1**, which was not even
   the constituent being blamed.

The lesson generalizes: when a fit misbehaves, suspect the *design matrix* as a
whole, not the constituent you happen to be looking at.

## 5.5 The co-tidal chart: turning phase into physics

With phases referenced to a common epoch (§4.4), the M2 phase across the array
gives a propagation direction and speed. Fitting a plane to phase versus
position:

$$g(x,y) = g_0 + \nabla g \cdot \mathbf{x}$$

gives a gradient of **0.230 ± 0.071 ° km⁻¹** toward **297°** — 3.2σ from zero,
so resolved, but only just. Converting to a phase speed,

$$c = \frac{360°/T}{|\nabla g|} = 35\ \mathrm{m\,s^{-1}}$$

Now the test. For a shallow-water wave $c = \sqrt{gh}$, so

$$h_{\text{eff}} = \frac{c^2}{g} = \frac{35^2}{9.81} \approx 125\ \mathrm{m}$$

The bank tops are ~19 m ($c = 14\ \mathrm{m\,s^{-1}}$) and the channel 1500 m
($c = 121\ \mathrm{m\,s^{-1}}$). The observed 35 m s⁻¹ sits between them and
matches neither. **The honest reading is that this array cannot separate the two
regimes**; it is 12 km across and the M2 phase varies by 2.2° over it.

> **Retraction.** This section previously reported **0.761 ± 0.102 ° km⁻¹ at
> 314°, 7.5σ, 11 m s⁻¹**, and concluded that $h_{\text{eff}} \approx 11$ m
> identified the wave as bank-top controlled. That is withdrawn — the gradient
> was a **short-baseline artifact**. The Angaur sub-array alone, 2.7 km across,
> returns 0.758 ° km⁻¹; every subset spanning ≥5 km returns 0.18–0.23. A real
> gradient does not depend on array size. `PRESSURE_ANALYSIS.md` §8.1 has the
> full diagnosis and `src/pressure_analysis.py` recomputes it.
>
> The methodological lesson survives intact, and is worth more than the result
> was: **a phase gradient divided by a short baseline is noise amplification**,
> and the check is to refit on sub-arrays and confirm the answer does not move.
> That test is now part of the driver.

This is the kind of result that only exists because the phases were referenced
correctly. With the 149° spread of §4.3 it would have been noise.

## 5.6 What harmonic analysis cannot do here

The array was used to estimate **tidal currents from pressure gradients**, via
the linearized momentum balance. That result was **withdrawn**. The history is
instructive:

- Initially reported 63.6 / 45.9 / 20.3 cm s⁻¹ for the three gauge triangles.
- Validation against the C05 ADCP gave $r = -0.18$ — no relationship.
- The first comparison was *also* methodologically wrong (tide-only prediction
  against total measured current). Redoing it constituent-by-constituent
  confirmed the failure: 2–7× too large, wrong phases.

Root cause: the Hydrographer Bank pressure-gradient **signal-to-noise ratio is
1.0**, and only 17 % / 10 % of the measured current is phase-locked to the
tide. The method fails by **bias, not noise** — averaging longer does not help,
because the thing being measured is mostly not there.

`PRESSURE_ANALYSIS.md` §4 keeps the withdrawal rather than deleting it, because
the failure mode is reusable: a pressure-gradient current estimate needs the
barotropic fraction to be large, and here 83 % of the flow is not barotropic.

## References

- Doodson, A. T. (1921). The harmonic development of the tide-generating
  potential. *Proc. R. Soc. Lond. A* **100**, 305–329.
- Foreman, M. G. G. (1977). *Manual for Tidal Heights Analysis and Prediction*.
  Pacific Marine Science Report 77-10.
- Pawlowicz, R., B. Beardsley and S. Lentz (2002). Classical tidal harmonic
  analysis including error estimates in MATLAB using T_TIDE.
  *Comput. Geosci.* **28**, 929–937.
- Codiga, D. L. (2011). *Unified Tidal Analysis and Prediction Using the UTide
  MATLAB Functions*. GSO Technical Report 2011-01.

Full citations with DOIs in `papers/README.md`.
