# 1. Velocity-gradient kinematics

How four drifters become a vorticity measurement.

Code: `src/eddy_kinematics.py`. MATLAB port: `matlab/eddy_kinematics_drifters.m`.

---

## 1.1 The problem

Vorticity is a *derivative* of the velocity field:

$$\zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}$$

You cannot measure a derivative at a point. You can only measure differences
between points, and we have exactly four of them, drifting. Everything below
is about extracting a derivative from four moving samples without fooling
yourself.

## 1.2 The local expansion

Take any smooth velocity field and expand it about the cluster centroid
$(\bar{x}, \bar{y})$ at time $t_c$:

$$u(x,y,t) = U_0 + \frac{\partial u}{\partial t}(t - t_c) + \frac{\partial u}{\partial x}(x - \bar{x}) + \frac{\partial u}{\partial y}(y - \bar{y}) + O(\Delta^2)$$

and likewise for $v$. Truncating after the linear terms is the **affine
approximation**. It is exact for solid-body rotation and for any linear shear;
it is wrong at second order in the cluster size $\Delta$ divided by the scale
over which the flow curves.

That error term is not academic here. The cluster spans ~640 m and the vortex
core is ~1200 m, so $\Delta/R \approx 0.5$ and the neglected term is not small.
§3 measures exactly what it costs.

## 1.3 The velocity-gradient tensor

Collect the four spatial derivatives into

$$\mathbf{A} = \begin{pmatrix} \partial_x u & \partial_y u \\ \partial_x v & \partial_y v \end{pmatrix}$$

Any 2×2 matrix decomposes uniquely into a trace, an antisymmetric part, and a
traceless symmetric part. Those three pieces are the entire local kinematics:

$$\mathbf{A} = \underbrace{\frac{\delta}{2}\begin{pmatrix}1&0\\0&1\end{pmatrix}}_{\text{divergence}} + \underbrace{\frac{\zeta}{2}\begin{pmatrix}0&-1\\1&0\end{pmatrix}}_{\text{rotation}} + \underbrace{\frac{1}{2}\begin{pmatrix}\sigma_n&\sigma_s\\\sigma_s&-\sigma_n\end{pmatrix}}_{\text{strain}}$$

with

| quantity | definition | what it does to a fluid blob |
|---|---|---|
| divergence $\delta$ | $\partial_x u + \partial_y v$ | changes its area |
| vorticity $\zeta$ | $\partial_x v - \partial_y u$ | spins it |
| normal strain $\sigma_n$ | $\partial_x u - \partial_y v$ | stretches along $x$, squeezes along $y$ |
| shear strain $\sigma_s$ | $\partial_x v + \partial_y u$ | stretches along the diagonals |

This is worth internalizing: **four numbers, four independent things a flow can
do to a small patch**. Nothing else is possible at first order.

## 1.4 Solving for the tensor: least squares

Each drifter sample in a time window contributes one row to a linear system.
With $n$ samples the design matrix is

$$\mathbf{M} = \begin{pmatrix} 1 & \Delta x_1 & \Delta y_1 & \Delta t_1 \\ \vdots & \vdots & \vdots & \vdots \\ 1 & \Delta x_n & \Delta y_n & \Delta t_n \end{pmatrix}$$

and we solve $\mathbf{M}\mathbf{c}_u = \mathbf{u}$ and $\mathbf{M}\mathbf{c}_v = \mathbf{v}$.
The two systems share $\mathbf{M}$ but not their coefficients, so they are
independent least-squares problems solved together only for convenience.

### Why scale the columns

The columns carry different units: $1$, meters, meters, seconds. A condition
number computed on the raw matrix is dominated by that unit mismatch and
reports ~2000, which tells you nothing. Normalizing each column to unit norm
before solving leaves the fit unchanged but makes $\mathrm{cond}(\mathbf{M})$ a
statement about **how well the drifters span the plane** — which is what you
actually want to monitor. After scaling it is ~2.2, meaning the cluster is a
healthy triangle-ish shape, not a sliver.

This matters because $\mathbf{A}$ is what gets inverted to find the eddy center
(§2.4), and an ill-conditioned $\mathbf{A}$ throws the center to infinity.

### The time column

Including $\Delta t$ lets the window absorb a uniform temporal acceleration of
the whole cluster, so a steadily strengthening background flow is not
misattributed to a spatial gradient. It costs one degree of freedom.

## 1.5 The independent route: circulation and Stokes' theorem

Stokes' theorem in the plane says that for any closed curve $C$ bounding a
region $S$,

$$\oint_C \mathbf{u}\cdot d\boldsymbol{\ell} = \iint_S \zeta \, dA$$

So the **area-averaged vorticity** is exactly the circulation divided by the
area:

$$\bar{\zeta} = \frac{1}{A}\oint_C \mathbf{u}\cdot d\boldsymbol{\ell}$$

This is remarkable and it is why it is the primary estimator here: it makes
**no assumption whatsoever** about the interior flow. The field can be as
complicated as it likes; if you can walk the boundary and integrate the
tangential velocity, you get the exact area average. For a cluster comparable
in size to the vortex — which is our situation — that is worth a great deal
more than an affine fit.

The divergence theorem gives the companion result with the outward normal:

$$\bar{\delta} = \frac{1}{A}\oint_C \mathbf{u}\cdot \hat{\mathbf{n}}\, d\ell$$

### What the code actually computes

We do not have the boundary; we have four points on it. The code takes the
quadrilateral joining them and applies the trapezoid rule to each straight
edge:

$$\Gamma \approx \sum_{i} \tfrac{1}{2}(\mathbf{u}_i + \mathbf{u}_{i+1})\cdot(\mathbf{r}_{i+1} - \mathbf{r}_i)$$

and the shoelace formula for the area:

$$A = \tfrac{1}{2}\left|\sum_i (x_i y_{i+1} - x_{i+1} y_i)\right|$$

**This introduces a quadrature error, and it is the largest systematic error in
the project.** For a strictly affine field the trapezoid rule is exact, so the
circulation and least-squares estimates agree to machine precision — verified
in `tests/test_synthetic_recovery.py`, which returns 0.0000 % error on a
solid-body field. For a curved field they differ, and the difference is a
measurement of how badly the affine assumption fails. §3.4 quantifies it:
about **−3.8 %** in $|\zeta|$ for the real geometry.

### Vertex ordering

The vertices must be walked in consistent rotational order or the polygon
self-intersects and the area is meaningless. The code sorts by
$\mathrm{atan2}$ of each vertex about the centroid. This is correct for convex
and mildly non-convex quadrilaterals; it would fail for a strongly re-entrant
polygon, which four drifters in a coherent vortex do not produce.

## 1.6 Gating on polygon shape

Vorticity from circulation is $\Gamma/A$. As the four drifters approach
collinearity, $A \to 0$ while $\Gamma$ stays finite and noisy, so the estimate
diverges. Before gating, leave-one-out triangles produced spikes to
$\zeta = 0.27\ \mathrm{s^{-1}}$ — a Rossby number of 15,000, which is nonsense.

The natural quantity to gate on is the **isoperimetric quotient**

$$Q = \frac{4\pi A}{P^2}$$

where $P$ is the perimeter. It is 1 for a circle, $\pi/4 \approx 0.785$ for a
square, and $\to 0$ as the vertices become collinear. Crucially it is
**scale-free**, so it rejects degenerate *shapes* without rejecting a merely
small cluster — which a bare area threshold would do.

The threshold is 0.10. It is a free parameter, so it must be shown not to
matter:

| $Q_{\min}$ | median $\zeta$ | change | epochs kept |
|---|---|---|---|
| 0.00 | −1.1885 × 10⁻³ | +0.14 % | 355 |
| 0.05 | −1.1893 × 10⁻³ | +0.07 % | 352 |
| **0.10** | **−1.1901 × 10⁻³** | — | 343 |
| 0.20 | −1.1843 × 10⁻³ | +0.49 % | 315 |
| 0.30 | −1.1721 × 10⁻³ | +1.52 % | 252 |

The answer moves by 1.5 % while a third of the data is discarded. The gate is
protecting against outliers, not manufacturing the result.

## 1.7 Okubo–Weiss: is it a vortex or a strain field?

Define

$$OW = \sigma_n^2 + \sigma_s^2 - \zeta^2$$

The sign of $OW$ determines the character of the local flow, and the reason is
worth seeing. The eigenvalues of the traceless part of $\mathbf{A}$ are

$$\lambda = \pm\tfrac{1}{2}\sqrt{\sigma_n^2 + \sigma_s^2 - \zeta^2} = \pm\tfrac{1}{2}\sqrt{OW}$$

- $OW > 0$: real eigenvalues, so there are two directions along which fluid
  moves in straight lines — a **saddle**. Strain dominates, and particles
  separate exponentially.
- $OW < 0$: imaginary eigenvalues, so trajectories are **closed ellipses**.
  Rotation dominates and material is trapped.

So $OW < 0$ is the mathematical statement of "this is a coherent vortex, not a
strained filament." This is Okubo (1970) and Weiss (1991).

For the Peleliu vortex, $OW < 0$ in **100 % of windows**. Read that with care:
the strain-to-vorticity ratio is 0.36, so

$$OW = \zeta^2(0.36^2 - 1) = -0.87\,\zeta^2$$

which is negative by a wide margin at every epoch. The 100 % is therefore not a
knife-edge result that happened to fall the right way — it is a robust
statement that rotation exceeds strain by roughly a factor of three throughout.

## 1.8 The third estimator: constellation rotation

For solid-body rotation at angular rate $\omega$, the vorticity is $\zeta = 2\omega$.
The factor of two catches people out. It follows directly: for
$\mathbf{u} = \omega\hat{z}\times\mathbf{r}$, so $u = -\omega y$ and $v = \omega x$,

$$\zeta = \partial_x v - \partial_y u = \omega - (-\omega) = 2\omega$$

So tracking the mean angle of the drifters about their own centroid, unwrapped,
gives a cumulative rotation that can be compared with $\int \zeta/2 \, dt$.
**This uses positions only** — no velocity data at all — so it is genuinely
independent of the other two estimators, which share the same $u, v$.

Over the record the constellation turned **−6.86 revolutions** where $\zeta/2$
integrated predicts **−8.95**. That 23 % shortfall is not an error; it is the
signature of a non-solid-body vortex, and §2.3 shows why.

### A bug worth knowing about

`np.unwrap` propagates a NaN through everything after it. The original
implementation unwrapped each drifter's raw angle column and summed across
drifters, so a single missing sample in one drifter destroyed the whole
series — one injected NaN moved the total from −6.86 to −1.24 turns and left
only 14 % of the series finite, **silently**. It did not bite on this data set
because there are no gaps larger than 421 s and none produce NaNs, but the
number is published, so the function now unwraps each drifter over its own
finite samples and references it to its own start.

## 1.9 Why three estimators

They fail differently:

| estimator | assumes | information used |
|---|---|---|
| circulation / Stokes | only that the boundary is walkable | $u, v$ on the boundary |
| least-squares affine | the field is affine across the cluster | $u, v$ everywhere in the window |
| constellation rotation | solid-body, for the comparison | positions only |

Agreement between the first two tests the affine assumption. Agreement with the
third tests whether the velocities and positions tell the same story, and the
*disagreement* (6.86 vs 8.95 turns) is itself the measurement in §2.3.

Divergence gets the same treatment, with a third estimator
$\delta = d(\ln A)/dt$ from the cluster area alone — again positions only.

## References

- Okubo, A. (1970). Horizontal dispersion of floatable particles in the
  vicinity of velocity singularities such as convergences. *Deep-Sea Res.*
  **17**, 445–454.
- Weiss, J. (1991). The dynamics of enstrophy transfer in two-dimensional
  hydrodynamics. *Physica D* **48**, 273–294.
- Okubo, A. and C. C. Ebbesmeyer (1976). Determination of vorticity, divergence
  and deformation rates from analysis of drogue observations. *Deep-Sea Res.*
  **23**, 349–352. — Eqs. (1)–(2) are exactly the least-squares model
  implemented here.
- Molinari, R. and A. D. Kirwan (1975). Calculations of differential kinematic
  properties from Lagrangian observations in the western Caribbean Sea.
  *J. Phys. Oceanogr.* **5**, 483–491.

Full citations with DOIs in `papers/README.md`.
