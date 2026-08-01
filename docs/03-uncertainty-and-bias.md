# 3. Uncertainty and bias

Three different things get called "the error bar," they measure different
things, and the one usually quoted is the least useful.

Code: `circulation_jackknife` in `src/eddy_kinematics.py`;
`tests/test_synthetic_recovery.py`.

---

## 3.1 The three errors

| | what it measures | what it is blind to |
|---|---|---|
| **formal** (least squares) | how much the fitted slope moves under the assumed noise model | model error, correlated noise, anything systematic |
| **resampling** (jackknife) | how much the answer moves when you drop data | anything common to every subsample |
| **systematic** (bias) | the offset between the estimator and truth | nothing — but you need truth to see it |

The project quotes the second. This section shows why the first is wrong by 5×,
why the second is also optimistic, and measures the third directly.

## 3.2 Formal error, and why it is too small

The least-squares covariance is

$$\mathrm{Cov}(\hat{\mathbf{c}}) = s^2 (\mathbf{M}^T\mathbf{M})^{-1}, \qquad s^2 = \frac{\|\mathbf{r}\|^2}{n - p}$$

Since $\zeta = \partial_x v - \partial_y u$ draws $\partial_x v$ from the $v$
system and $\partial_y u$ from the $u$ system, and those are fit independently,
the code adds variances:

$$\sigma_\zeta = \sqrt{\mathrm{Var}(\partial_x v) + \mathrm{Var}(\partial_y u)}$$

This gives a median $\sigma_\zeta = 6.0\times10^{-5}$, i.e. **5 % of $|\zeta|$**.

It is too small for three reasons, in increasing order of importance:

1. **It assumes the $u$ and $v$ errors are uncorrelated**, so the cross-term is
   dropped. Directional surface waves correlate them: a wave from a given
   direction perturbs $u$ and $v$ together. The omitted term is not zero.

2. **It assumes independent samples.** Within a 30-minute window the residuals
   are serially correlated, so the effective degrees of freedom are fewer than
   $n - p$.

3. **It cannot know the model is wrong.** The formal error describes scatter
   *about the fitted affine plane*. It has no way to express "the true field is
   curved and no plane fits it." That is the dominant term, and it is invisible
   to the covariance by construction.

Molinari and Kirwan hit the same wall in 1975: their vorticity series were
"ragged with frequent changes in sign" wherever shear was small relative to
observational error, which is what an underestimated error bar looks like from
the outside.

**Measured against synthetic truth**, the formal error is **5.2× too small**
(`tests/test_synthetic_recovery.py`, and the standalone run in §3.5).

## 3.3 The jackknife, and what it cannot see

Four drifters give four distinct triangles. For a strictly affine field every
triangle returns the same vorticity, so the spread across them measures the
combined effect of measurement noise and of curvature across the cluster —
assuming neither linearity nor a noise model.

The code uses **median and MAD**, not mean and standard deviation:

$$\hat{\sigma} = 1.4826 \times \mathrm{median}_k\left|\zeta_k - \mathrm{median}_j\,\zeta_j\right|$$

With only four values a single near-degenerate triangle would otherwise
dominate both moments. The 1.4826 converts MAD to a Gaussian-equivalent
standard deviation.

Result: **17 % of $|\zeta|$**, about three times the formal error. That is the
number `DRIFTER_ANALYSIS.md` quotes, and it is the right one to quote of the
two.

### But it is still optimistic

**A jackknife cannot see an error that is common to every subsample.** All four
triangles share the same quadrature approximation, the same center, the same
preprocessing. Whatever those get wrong, they get wrong together, and the
spread across subsamples says nothing about it.

Measured on synthetic data where truth is known, with velocity noise added:

| velocity noise | LOO spread | actual \|error\| | ratio |
|---|---|---|---|
| 0.00 m/s | 4.41 × 10⁻⁶ | 1.23 × 10⁻⁴ | **0.04×** |
| 0.01 | 1.22 × 10⁻⁵ | 1.30 × 10⁻⁴ | 0.09× |
| 0.02 | 2.25 × 10⁻⁵ | 1.33 × 10⁻⁴ | 0.17× |
| 0.05 | 4.60 × 10⁻⁵ | 1.22 × 10⁻⁴ | 0.38× |
| 0.10 | 1.06 × 10⁻⁴ | 1.72 × 10⁻⁴ | 0.62× |

With no noise the jackknife reports essentially zero uncertainty while being
wrong by $1.2\times10^{-4}$ — because the error is entirely systematic and
entirely common-mode. As noise grows the jackknife starts to see it, but even
at 0.10 m s⁻¹ (well above the observed 0.064 m s⁻¹ fit residual) it captures
only ~62 % of the true error.

**Interpretation for this project**: the quoted ±17 % is a *precision*
estimate. Total accuracy is roughly a factor of two worse, i.e. of order
±30 %, once the common-mode terms of §3.4 are included.

## 3.4 The systematic bias, measured

The circulation estimator evaluates the closed line integral with the trapezoid
rule on four vertex velocities. For a curved flow that is approximate, and the
error is a *bias*, not noise — it does not average out over 355 epochs.

Measuring it requires knowing the true field, so the test assumes the fitted
Lamb–Oseen ($\Gamma = -6114$, $R = 1206$ m) is truth, then for **each real
polygon** compares:

- the trapezoid rule on the four vertices (what the code does), against
- 200-point Gauss–Legendre quadrature along the same four straight edges (the
  exact line integral).

Both integrate over the identical path, so the difference isolates the
quadrature error and nothing else.

**Result: $|\zeta|$ is underestimated by 3.8 %** (median; IQR −7.5 % to −1.2 %,
n = 354). Correcting:

$$\zeta = -1.190\times10^{-3} \;\longrightarrow\; -1.236\times10^{-3}, \qquad Ro = -67.0 \;\longrightarrow\; -69.6$$

Two things are worth noting:

- **The direction is predictable.** Along a chord through a vortex the
  tangential velocity is larger at the midpoint than at the endpoints, so a
  two-point trapezoid rule underestimates the integral and hence $|\Gamma|$.
- **It is driven by proximity to the core, not polygon shape.** The cluster
  sits a median 534 m from the center, well inside $R = 1206$ m, where $\zeta$
  varies strongly across the polygon. Splitting by aspect ratio gives −3.75 %
  for round polygons and −3.64 % for elongated ones — essentially identical,
  which disproves the natural guess that elongation drives it.

3.8 % is comfortably inside the ±17 % precision, so **the headline number
stands**. But it is a bias, so it belongs in the accuracy budget rather than
being absorbed into the noise.

## 3.5 The null tests

`tests/test_synthetic_recovery.py` runs four checks on flows whose answer is
known by construction. It is the only place in the project where "correct" is
defined independently of the estimator.

The design matters: trajectories are **integrated** (RK4) through the synthetic
field starting from the real initial positions, so positions and velocities are
mutually consistent. An earlier version imposed a straight-line center track on
the *real* drifter positions; the drifters then never orbited that center
(median radius 2461 m against a 1200 m core) and every downstream number
measured the mistake rather than the code. Consistency is not a detail here.

| test | checks | result |
|---|---|---|
| 1. solid-body exactness | an affine field must be recovered exactly | **0.0000 % error** |
| 2. solid-body null | the "not solid body" diagnostic must not fire on a solid body | **+0.06 / −0.05** vs −0.56/−0.74 for Lamb–Oseen |
| 3. Lamb–Oseen recovery | $\Gamma$ and $R$ come back from a known vortex | $\Gamma$ +2.7 %, $R$ −1.7 % |
| 4. quadrature bias | trapezoid vs exact on real geometry | **−3.8 %** in $|\zeta|$ |

Test 1 is the sharpest: for a solid-body field the trapezoid rule is exact and
Stokes' theorem is exact, so **any** deviation is an implementation bug rather
than a modeling limitation. Returning 0.0000 % is a strong statement that the
circulation machinery is correctly coded.

Test 2 is the one that protects a published conclusion. Without it, −0.59 is
just a number that agrees with a hypothesis.

## 3.6 Other free parameters

An answer that depends on a choice nobody justified is not a measurement. The
two main choices, both varied over their plausible range:

**Block-averaging length** (waves must be removed; how hard?):

| blocks/burst | averaging | median $\zeta$ | change |
|---|---|---|---|
| 2 | 511.8 s | −1.1890 × 10⁻³ | +0.09 % |
| 4 | 255.9 s | −1.1901 × 10⁻³ | — |
| 8 | 127.9 s | −1.1969 × 10⁻³ | −0.57 % |
| 16 | 64.0 s | −1.1922 × 10⁻³ | −0.17 % |

**0.6 % across a 16× range.** The wave band and the eddy band are cleanly
separated, which is the whole reason this works.

**Polygon quality gate**: 1.5 % across $Q_{\min} \in [0, 0.3]$ — see §1.6.

Neither is doing any work in setting the answer.

## 3.7 What to write in a paper

> Relative vorticity $\zeta = -1.19\times10^{-3}\ \mathrm{s^{-1}}$
> (Rossby $-67$), with a leave-one-out precision of 17 % and an estimated
> systematic underestimate of 4 % from the finite-difference quadrature over a
> four-vertex polygon comparable in size to the vortex core. Total accuracy is
> of order 30 %.

Not "$-1.19 \pm 0.06 \times 10^{-3}$", which is the formal error and is wrong
by a factor of five.

## References

- Efron, B. and R. J. Tibshirani (1993). *An Introduction to the Bootstrap*.
  Chapman & Hall. — jackknife and bootstrap, chapters 11 and 12.
- Molinari, R. and A. D. Kirwan (1975). *J. Phys. Oceanogr.* **5**, 483–491.

Full citations with DOIs in `papers/README.md`.
